#pragma once
#include "global.hpp"
#include "batch.hpp"
#include "container.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <zstd.h>

// Merge N .lhb batch files → .lhg (global, v5 inverted) + .lhgi (index).
//
// v5 inverts the per-accession .lhb blocks into position-centric records.
// The .lhb format stays v4 (unchanged); inversion happens here at merge time.
//
//   Pass 1: scan each .lhb sequentially (one open at a time) → collect
//           BatchBlockRef (offsets only, no payload read). Sort by (hog_id, acc_id).
//
//   Between passes: build the global accession registry — collect all unique
//           acc_ids, sort alphabetically → acc_id_to_idx (stable integer IDs).
//
//   Pass 2: iterate HOG groups. For each HOG:
//             For each block (acc_id, batch_file_idx, offset):
//               pread ShardBlockHeader + compressed payload → decompress.
//               Parse local contig dict; deserialize_varnt_block().
//               For each record r, for each obs:
//                 hog_pos   = r.sstart + obs.hog_offset
//                 codon_idx = obs.packed_codon >> 2
//                 unitig_id = local_contigs[r.contig_idx]
//                 inverted[hog_pos].push_back({acc_idx, codon_idx, unitig_id})
//             Sort each position's obs by acc_idx.
//             Build per-HOG local unitig dict; serialize inverted block.
//             HOG-level zstd compress (level 19). Write entry.
//
//   Pass 3: write HOG index + accession registry → .lhgi; backfill header.
//
// FDs: at most 2 open simultaneously (one source .lhb + the output .lhg).

namespace lhi {

inline void merge_batches(const std::vector<std::string>& batch_paths,
                          const std::string& out_lhg,
                          const std::string& out_lhgi,
                          int zstd_level = 19) {

    // ── Pass 1: scan all .lhb, collect refs ───────────────────────────────────
    std::vector<BatchBlockRef> refs;
    refs.reserve(batch_paths.size() * 512);

    for (size_t fi = 0; fi < batch_paths.size(); ++fi) {
        bool ok = scan_batch_file(batch_paths[fi], fi,
            [](const std::string&) {},
            [&](BatchBlockRef br) { refs.push_back(std::move(br)); }
        );
        if (!ok) throw std::runtime_error("cannot open batch: " + batch_paths[fi]);
        if (batch_paths.size() > 100 && (fi + 1) % 1000 == 0)
            std::cerr << "  scanned " << (fi + 1) << "/" << batch_paths.size() << " batches\r";
    }
    if (batch_paths.size() > 100) std::cerr << "\n";
    std::cerr << "merge: " << refs.size() << " blocks from "
              << batch_paths.size() << " batches\n";

    // Sort by (hog_id, acc_id) so HOG groups are contiguous.
    std::stable_sort(refs.begin(), refs.end(),
        [](const BatchBlockRef& a, const BatchBlockRef& b) {
            if (a.hog_id != b.hog_id) return a.hog_id < b.hog_id;
            return a.acc_id < b.acc_id;
        });

    // ── Build global accession registry ───────────────────────────────────────
    std::vector<std::string> accessions;
    accessions.reserve(refs.size());
    for (const auto& r : refs) accessions.push_back(r.acc_id);
    std::sort(accessions.begin(), accessions.end());
    accessions.erase(std::unique(accessions.begin(), accessions.end()), accessions.end());

    std::unordered_map<std::string, uint32_t> acc_id_to_idx;
    acc_id_to_idx.reserve(accessions.size() * 2);
    for (uint32_t i = 0; i < accessions.size(); ++i) acc_id_to_idx[accessions[i]] = i;

    // ── Open output file (sequential writes via write_all) ────────────────────
    int out_fd = open(out_lhg.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (out_fd < 0) throw std::runtime_error("cannot create: " + out_lhg);

    auto write_all = [&](const void* p, size_t n) {
        const char* buf = (const char*)p;
        size_t done = 0;
        while (done < n) {
            ssize_t r = ::write(out_fd, buf + done, n - done);
            if (r <= 0) throw std::runtime_error("write failed: " + out_lhg);
            done += r;
        }
    };

    // File header placeholder (index_offset filled in later via pwrite).
    uint8_t file_hdr[LHG_HEADER_SZ] = {};
    memcpy(file_hdr, LHG_FILE_MAGIC, 4);
    file_hdr[4] = LHG_VERSION;
    write_all(file_hdr, LHG_HEADER_SZ);
    uint64_t write_pos = LHG_HEADER_SZ;

    // ── Pass 2: HOG-streaming inversion + compression ─────────────────────────
    std::vector<HogIndexEntry> index_entries;
    index_entries.reserve(4096);

    std::vector<uint8_t> cbuf_in;     // temp buffer for reading compressed block
    std::vector<uint8_t> raw_block;   // decompressed single block
    std::vector<uint8_t> inv_raw;     // serialized inverted block
    std::vector<uint8_t> hog_cbuf;    // output HOG compressed buffer

    size_t group_start = 0;
    while (group_start < refs.size()) {
        size_t group_end = group_start + 1;
        while (group_end < refs.size() && refs[group_end].hog_id == refs[group_start].hog_id)
            ++group_end;

        const std::string& hog_id = refs[group_start].hog_id;
        size_t n_blocks = group_end - group_start;

        // hog_pos → observation list (built by inverting all blocks).
        std::unordered_map<uint32_t, std::vector<InvObs>> inverted;
        // per-HOG local unitig dict: store contig numbers, dedup by full string.
        std::vector<uint32_t> unitigs;
        std::unordered_map<std::string, uint32_t> unitig_to_idx;
        auto unitig_idx_of = [&](const std::string& uid) -> uint32_t {
            auto it = unitig_to_idx.find(uid);
            if (it != unitig_to_idx.end()) return it->second;
            uint32_t idx = uint32_t(unitigs.size());
            auto pos = uid.rfind('_');
            uint32_t cnum = uint32_t(std::stoul(uid.substr(pos + 1)));
            unitigs.push_back(cnum);
            unitig_to_idx.emplace(uid, idx);
            return idx;
        };

        {
            size_t cur_file_idx = SIZE_MAX;
            int    src_fd       = -1;
            struct FdGuard {
                int& fd;
                ~FdGuard() { if (fd >= 0) { close(fd); fd = -1; } }
            } guard{src_fd};

            for (size_t bi = 0; bi < n_blocks; ++bi) {
                const BatchBlockRef& ref = refs[group_start + bi];

                if (ref.batch_file_idx != cur_file_idx) {
                    if (src_fd >= 0) { close(src_fd); src_fd = -1; }
                    src_fd = open(batch_paths[ref.batch_file_idx].c_str(), O_RDONLY);
                    if (src_fd < 0)
                        throw std::runtime_error("cannot reopen: " + batch_paths[ref.batch_file_idx]);
                    cur_file_idx = ref.batch_file_idx;
                }

                ShardBlockHeader shard_hdr;
                ssize_t nr = ::pread(src_fd, &shard_hdr, sizeof(shard_hdr), ref.shard_hdr_offset);
                if (nr != ssize_t(sizeof(shard_hdr)))
                    throw std::runtime_error("pread header failed: " + batch_paths[ref.batch_file_idx]);

                cbuf_in.resize(shard_hdr.compressed_sz);
                nr = ::pread(src_fd, cbuf_in.data(), shard_hdr.compressed_sz,
                             ref.shard_hdr_offset + sizeof(shard_hdr));
                if (nr != ssize_t(shard_hdr.compressed_sz))
                    throw std::runtime_error("pread payload failed: " + batch_paths[ref.batch_file_idx]);

                raw_block.resize(shard_hdr.raw_sz);
                size_t rz = ZSTD_decompress(raw_block.data(), shard_hdr.raw_sz,
                                            cbuf_in.data(), shard_hdr.compressed_sz);
                if (ZSTD_isError(rz))
                    throw std::runtime_error(std::string("zstd decompress: ") + ZSTD_getErrorName(rz));

                // Parse local contig dict (the unitig_id strings for this block).
                const uint8_t* p   = raw_block.data();
                const uint8_t* end = p + rz;
                uint32_t n_contigs = 0;
                int n = read_varint(p, end, &n_contigs);
                if (!n) throw std::runtime_error("corrupt contig dict in block for HOG " + hog_id);
                p += n;
                std::vector<std::string> local_contigs(n_contigs);
                for (uint32_t j = 0; j < n_contigs; ++j) {
                    uint32_t len = 0;
                    n = read_varint(p, end, &len);
                    if (!n || p + n + len > end) throw std::runtime_error("truncated contig dict");
                    p += n;
                    local_contigs[j].assign(reinterpret_cast<const char*>(p), len);
                    p += len;
                }

                std::vector<VarNTRecord> recs;
                if (!deserialize_varnt_block(p, end, recs))
                    throw std::runtime_error("corrupt varnt block for HOG " + hog_id);

                uint32_t acc_idx = acc_id_to_idx.at(ref.acc_id);
                for (const auto& r : recs) {
                    const std::string& unitig_id = local_contigs[r.contig_idx];
                    uint32_t u_idx = unitig_idx_of(unitig_id);
                    for (const auto& o : r.vars) {
                        uint32_t hog_pos   = r.sstart + o.hog_offset;
                        uint8_t  codon_idx = uint8_t(o.packed_codon >> 2);
                        inverted[hog_pos].push_back({acc_idx, codon_idx, u_idx});
                    }
                }
            }
        } // FdGuard closes src_fd

        // Sort positions ascending; sort each position's obs by acc_idx.
        std::vector<InvPosition> positions;
        positions.reserve(inverted.size());
        for (auto& kv : inverted) {
            std::sort(kv.second.begin(), kv.second.end(),
                [](const InvObs& a, const InvObs& b) { return a.acc_idx < b.acc_idx; });
            positions.push_back({kv.first, std::move(kv.second)});
        }
        std::sort(positions.begin(), positions.end(),
            [](const InvPosition& a, const InvPosition& b) { return a.hog_pos < b.hog_pos; });

        // Serialize inverted block.
        inv_raw.clear();
        serialize_inverted_block(inv_raw, unitigs, positions);

        // HOG-level zstd compress.
        size_t bound = ZSTD_compressBound(inv_raw.size());
        hog_cbuf.resize(bound);
        size_t csz = ZSTD_compress(hog_cbuf.data(), bound,
                                   inv_raw.data(), inv_raw.size(), zstd_level);
        if (ZSTD_isError(csz))
            throw std::runtime_error(std::string("zstd HOG compress: ") + ZSTD_getErrorName(csz));

        // ── Write HOG entry: "LHHE" + stored_sz + payload ─────────────────────
        // Raw fallback: if zstd doesn't shrink the block, store raw (high bit set).
        bool use_raw = (csz >= inv_raw.size());
        uint32_t stored_sz = use_raw ? (uint32_t(inv_raw.size()) | 0x80000000u)
                                     : uint32_t(csz);
        const uint8_t* payload = use_raw ? inv_raw.data() : hog_cbuf.data();
        size_t payload_sz = use_raw ? inv_raw.size() : csz;

        uint64_t hog_data_offset = write_pos;
        write_all(LHG_HOG_ENTRY_MAGIC, 4); write_pos += 4;
        uint8_t csz4[4]; for (int i = 0; i < 4; ++i) csz4[i] = uint8_t(stored_sz >> (8 * i));
        write_all(csz4, 4); write_pos += 4;
        write_all(payload, payload_sz); write_pos += payload_sz;
        uint64_t hog_data_length = write_pos - hog_data_offset;

        // Distinct accessions contributing to this HOG.
        uint32_t n_accs = 0;
        std::string prev_acc;
        for (size_t bi = 0; bi < n_blocks; ++bi) {
            if (refs[group_start + bi].acc_id != prev_acc) {
                prev_acc = refs[group_start + bi].acc_id; ++n_accs;
            }
        }

        index_entries.push_back({hog_id, hog_data_offset, hog_data_length, n_accs});
        group_start = group_end;
    }

    uint64_t index_offset = write_pos;

    // ── Pass 3: write index + accession registry; fill in file header ─────────
    auto idx_bytes = build_index_bytes(index_entries, accessions);
    write_all(idx_bytes.data(), idx_bytes.size());

    uint8_t final_hdr[LHG_HEADER_SZ] = {};
    memcpy(final_hdr, LHG_FILE_MAGIC, 4);
    final_hdr[4] = LHG_VERSION;
    for (int i = 0; i < 8; ++i) final_hdr[8+i] = uint8_t(index_offset >> (8*i));
    ssize_t pw = ::pwrite(out_fd, final_hdr, LHG_HEADER_SZ, 0);
    if (pw != ssize_t(LHG_HEADER_SZ)) throw std::runtime_error("pwrite header failed");
    close(out_fd);

    int idx_fd = open(out_lhgi.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (idx_fd < 0) throw std::runtime_error("cannot create: " + out_lhgi);
    if (::write(idx_fd, idx_bytes.data(), idx_bytes.size()) != ssize_t(idx_bytes.size()))
        throw std::runtime_error("write .lhgi failed");
    close(idx_fd);

    std::cerr << "merge done: " << index_entries.size() << " HOGs, "
              << accessions.size() << " accessions → " << out_lhg << "\n"
              << "index:      " << out_lhgi
              << " (" << idx_bytes.size() / 1024 << " KB)\n";
}

} // namespace lhi
