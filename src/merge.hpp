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
#include <charconv>
#include <optional>
#include <set>
#include <fstream>
#include <zstd.h>

// Merge N inputs (.lhb per-accession batches and/or already-inverted .lhg
// shards) into one position-centric .lhg plus its .lhgi index. For .lhb inputs
// inversion happens here; .lhg inputs are re-inverted into the global acc space.

namespace lhi {

// Unified per-HOG source reference covering both input kinds.
struct MergeRef {
    enum class Kind { LHB, LHG };
    Kind     kind;
    std::string hog_id;
    std::vector<std::string> acc_ids;   // all accs contributing to this HOG in this source
    size_t   source_file_idx = 0;
    // LHB: offset of the ShardBlockHeader ("LHSB") in the .lhb file.
    uint64_t lhb_shard_hdr_offset = 0;
    // LHG: byte range of the HOG entry ("LHHE"…) in the .lhg file.
    uint64_t lhg_data_offset = 0;
    uint64_t lhg_data_length = 0;
};

// Detect input kind by first-4-byte file magic.
//   LHB_FILE_MAGIC "LHGB" → .lhb ;  LHG_FILE_MAGIC "LHGG" → .lhg
inline MergeRef::Kind detect_input_kind(const std::string& path) {
    UniqueFd fd(open(path.c_str(), O_RDONLY));
    if (fd < 0) throw std::runtime_error("cannot open input: " + path);
    uint8_t m[4];
    if (!fd_read_exact(fd, m, 4)) throw std::runtime_error("truncated input: " + path);
    if (memcmp(m, LHB_FILE_MAGIC, 4) == 0) return MergeRef::Kind::LHB;
    if (memcmp(m, LHG_FILE_MAGIC, 4) == 0) return MergeRef::Kind::LHG;
    throw std::runtime_error("unrecognized magic (not .lhb/.lhg): " + path);
}

// Merge multiple InvBlocks (one per source .lhg) into a single InvBlock.
// acc_remap[i] maps source i's local acc_idx → global acc_idx.
// Unitig dicts are local per-block; output gets a merged local dict whose
// numeric contig values are deduped per (global_acc_idx, contig_num).
inline InvBlock merge_inv_blocks(
        std::vector<std::pair<InvBlock, std::vector<uint32_t>>> sources) {

    // hog_pos → merged obs (acc_idx already global; unitig_idx into merged dict).
    std::unordered_map<uint32_t, std::vector<InvObs>> merged;
    std::vector<uint32_t> out_unitigs;
    std::unordered_map<uint64_t, uint32_t> uni_key_to_idx;  // (acc<<32|cnum) → idx
    auto unitig_idx_of = [&](uint32_t gacc, uint32_t cnum) -> uint32_t {
        uint64_t key = (uint64_t(gacc) << 32) | cnum;
        auto it = uni_key_to_idx.find(key);
        if (it != uni_key_to_idx.end()) return it->second;
        uint32_t idx = uint32_t(out_unitigs.size());
        out_unitigs.push_back(cnum);
        uni_key_to_idx.emplace(key, idx);
        return idx;
    };

    for (auto& [blk, remap] : sources) {
        const uint8_t* p   = blk.pos_ptr;
        const uint8_t* end = blk.end;
        uint32_t n_positions = 0;
        int n = read_varint(p, end, &n_positions);
        if (!n) throw std::runtime_error("merge_inv_blocks: corrupt position count");
        p += n;
        std::vector<InvObs> obs;
        uint32_t prev_pos = 0;
        for (uint32_t pi = 0; pi < n_positions; ++pi) {
            uint32_t hog_pos = 0;
            if (!decode_position(p, end, prev_pos, hog_pos, obs))
                throw std::runtime_error("merge_inv_blocks: corrupt position record");
            auto& dst = merged[hog_pos];
            for (const auto& o : obs) {
                if (o.acc_idx >= remap.size())
                    throw std::runtime_error("merge_inv_blocks: acc_idx out of remap range");
                uint32_t gacc = remap[o.acc_idx];
                if (o.unitig_idx >= blk.unitigs.size())
                    throw std::runtime_error("merge_inv_blocks: unitig_idx OOB");
                uint32_t cnum = blk.unitigs[o.unitig_idx];
                dst.push_back({gacc, o.codon_idx, unitig_idx_of(gacc, cnum)});
            }
        }
    }

    std::vector<InvPosition> positions;
    positions.reserve(merged.size());
    for (auto& kv : merged) {
        std::sort(kv.second.begin(), kv.second.end(),
            [](const InvObs& a, const InvObs& b) { return a.acc_idx < b.acc_idx; });
        positions.push_back({kv.first, std::move(kv.second)});
    }
    std::sort(positions.begin(), positions.end(),
        [](const InvPosition& a, const InvPosition& b) { return a.hog_pos < b.hog_pos; });

    InvBlock out;
    out.unitigs = std::move(out_unitigs);
    serialize_inverted_block(out.raw, out.unitigs, positions);
    out.pos_ptr = nullptr;  // serialized form; callers read out.raw / out.unitigs
    out.end     = out.raw.data() + out.raw.size();
    return out;
}

// accregistry subcommand: collect unique acc_ids from .lhgi LHGA sections,
// sort, write one per line (plain text, no header).
inline void build_acc_registry(const std::vector<std::string>& lhgi_paths,
                               const std::string& out_path) {
    std::set<std::string> accs;
    for (const auto& path : lhgi_paths) {
        GlobalIndex idx;
        if (!idx.load(path) && !idx.load_from_lhg(path))
            throw std::runtime_error("cannot load index: " + path);
        for (const auto& a : idx.accessions) accs.insert(a);
    }
    std::ofstream out(out_path, std::ios::trunc);
    if (!out) throw std::runtime_error("cannot create: " + out_path);
    for (const auto& a : accs) out << a << "\n";
    std::cerr << "accregistry: " << accs.size() << " accessions from "
              << lhgi_paths.size() << " indexes → " << out_path << "\n";
}

// Load a global acc registry (one acc_id per line; position = global acc_idx).
inline std::vector<std::string> load_acc_registry(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open acc registry: " + path);
    std::vector<std::string> accs;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) accs.push_back(line);
    }
    return accs;
}

inline void merge_batches(const std::vector<std::string>& input_paths,
                          const std::string& out_lhg,
                          const std::string& out_lhgi,
                          int zstd_level = 19,
                          const std::string& hog_range_start = "",
                          const std::string& hog_range_end   = "",
                          const std::string& acc_registry_path = "") {

    // Pass 1: scan all inputs, collect per-HOG MergeRefs (offsets only).
    std::vector<MergeRef> refs;
    refs.reserve(input_paths.size() * 512);

    for (size_t fi = 0; fi < input_paths.size(); ++fi) {
        MergeRef::Kind kind = detect_input_kind(input_paths[fi]);
        if (kind == MergeRef::Kind::LHB) {
            std::string cur_acc;
            bool ok = scan_batch_file(input_paths[fi], fi,
                [&](const std::string& a) { cur_acc = a; },
                [&](BatchBlockRef br) {
                    MergeRef r;
                    r.kind = MergeRef::Kind::LHB;
                    r.hog_id = std::move(br.hog_id);
                    r.acc_ids = {br.acc_id};
                    r.source_file_idx = fi;
                    r.lhb_shard_hdr_offset = uint64_t(br.shard_hdr_offset);
                    refs.push_back(std::move(r));
                });
            if (!ok) throw std::runtime_error("cannot open batch: " + input_paths[fi]);
        } else {
            // .lhg: read its LHGI index; one MergeRef per HOG entry, carrying
            // the full local accession registry (per-HOG accs unknown without
            // decoding, so all source accs are recorded for remapping).
            GlobalIndex idx;
            if (!idx.load_from_lhg(input_paths[fi]))
                throw std::runtime_error("cannot load .lhg index: " + input_paths[fi]);
            for (const auto& e : idx.entries) {
                MergeRef r;
                r.kind = MergeRef::Kind::LHG;
                r.hog_id = e.hog_id;
                r.acc_ids = idx.accessions;   // source's full local acc registry
                r.source_file_idx = fi;
                r.lhg_data_offset = e.data_offset;
                r.lhg_data_length = e.data_length;
                refs.push_back(std::move(r));
            }
        }
        if (input_paths.size() > 100 && (fi + 1) % 1000 == 0)
            std::cerr << "  scanned " << (fi + 1) << "/" << input_paths.size() << " inputs\r";
    }
    if (input_paths.size() > 100) std::cerr << "\n";
    std::cerr << "merge: " << refs.size() << " HOG-blocks from "
              << input_paths.size() << " inputs\n";

    // (hog_id, source_file_idx) order makes HOG groups contiguous for streaming.
    std::stable_sort(refs.begin(), refs.end(),
        [](const MergeRef& a, const MergeRef& b) {
            if (a.hog_id != b.hog_id) return a.hog_id < b.hog_id;
            return a.source_file_idx < b.source_file_idx;
        });

    // Global accession registry.
    std::vector<std::string> accessions;
    if (!acc_registry_path.empty()) {
        accessions = load_acc_registry(acc_registry_path);
    } else {
        std::set<std::string> uniq;
        for (const auto& r : refs)
            for (const auto& a : r.acc_ids) uniq.insert(a);
        accessions.assign(uniq.begin(), uniq.end());
    }

    std::unordered_map<std::string, uint32_t> acc_id_to_idx;
    acc_id_to_idx.reserve(accessions.size() * 2);
    for (uint32_t i = 0; i < accessions.size(); ++i) acc_id_to_idx[accessions[i]] = i;

    // Per-source acc_idx → global acc_idx remap tables (for .lhg inversion).
    // Validate every input acc exists in the registry (error if registry given).
    auto global_acc = [&](const std::string& acc_id) -> uint32_t {
        auto it = acc_id_to_idx.find(acc_id);
        if (it == acc_id_to_idx.end())
            throw std::runtime_error("acc_id not in registry: " + acc_id);
        return it->second;
    };

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

    // Pass 2: HOG-streaming inversion + compression.
    std::vector<HogIndexEntry> index_entries;
    index_entries.reserve(4096);

    std::vector<uint8_t> cbuf_in;     // temp buffer for reading compressed block
    std::vector<uint8_t> raw_block;   // decompressed single block
    std::vector<uint8_t> inv_raw;     // serialized inverted block
    std::vector<uint8_t> hog_cbuf;    // output HOG compressed buffer

    // Reused zstd contexts (one alloc for the whole merge, not per HOG/block).
    struct ZstdCtx {
        ZSTD_CCtx* c = ZSTD_createCCtx();
        ZSTD_DCtx* d = ZSTD_createDCtx();
        ~ZstdCtx() { ZSTD_freeCCtx(c); ZSTD_freeDCtx(d); }
    } zc;
    ZSTD_CCtx_setParameter(zc.c, ZSTD_c_compressionLevel, zstd_level);

    size_t group_start = 0;
    while (group_start < refs.size()) {
        size_t group_end = group_start + 1;
        while (group_end < refs.size() && refs[group_end].hog_id == refs[group_start].hog_id)
            ++group_end;

        const std::string& hog_id = refs[group_start].hog_id;
        size_t n_blocks = group_end - group_start;

        // --hog-range: skip HOGs outside [start, end) (lexicographic).
        if ((!hog_range_start.empty() && hog_id < hog_range_start) ||
            (!hog_range_end.empty()   && hog_id >= hog_range_end)) {
            group_start = group_end;
            continue;
        }

        // hog_pos → observation list (built by inverting all blocks).
        std::unordered_map<uint32_t, std::vector<InvObs>> inverted;
        // per-HOG local unitig dict: store contig numbers, dedup by
        // (global_acc_idx, contig_num) — the unitig identity.
        std::vector<uint32_t> unitigs;
        std::unordered_map<uint64_t, uint32_t> unitig_to_idx;
        auto unitig_idx_of = [&](uint32_t gacc, uint32_t cnum) -> uint32_t {
            uint64_t key = (uint64_t(gacc) << 32) | cnum;
            auto it = unitig_to_idx.find(key);
            if (it != unitig_to_idx.end()) return it->second;
            uint32_t idx = uint32_t(unitigs.size());
            unitigs.push_back(cnum);
            unitig_to_idx.emplace(key, idx);
            return idx;
        };
        // Distinct global accessions actually contributing to this HOG.
        std::set<uint32_t> contributing_accs;

        {
            size_t cur_file_idx = SIZE_MAX;
            std::optional<UniqueFd> src;

            for (size_t bi = 0; bi < n_blocks; ++bi) {
                const MergeRef& ref = refs[group_start + bi];

                if (ref.source_file_idx != cur_file_idx) {
                    src.emplace(open(input_paths[ref.source_file_idx].c_str(), O_RDONLY));
                    if (*src < 0)
                        throw std::runtime_error("cannot reopen: " + input_paths[ref.source_file_idx]);
                    cur_file_idx = ref.source_file_idx;
                }
                int src_fd = *src;

                if (ref.kind == MergeRef::Kind::LHB) {
                    uint32_t acc_idx = global_acc(ref.acc_ids.front());

                    uint8_t hdr_bytes[28];
                    ssize_t nr = ::pread(src_fd, hdr_bytes, sizeof(hdr_bytes), off_t(ref.lhb_shard_hdr_offset));
                    if (nr != ssize_t(sizeof(hdr_bytes)))
                        throw std::runtime_error("pread header failed: " + input_paths[ref.source_file_idx]);
                    if (hdr_bytes[4] != SHARD_BLOCK_VERSION)
                        throw std::runtime_error("unsupported shard version in: " + input_paths[ref.source_file_idx]);
                    uint32_t compressed_sz = read_u32_le(hdr_bytes + 8);
                    uint32_t raw_sz        = read_u32_le(hdr_bytes + 12);

                    if (compressed_sz > 256u * 1024 * 1024)
                        throw std::runtime_error("compressed_sz OOB for HOG " + hog_id);
                    if (raw_sz > 512u * 1024 * 1024)
                        throw std::runtime_error("raw_sz OOB for HOG " + hog_id);

                    cbuf_in.resize(compressed_sz);
                    nr = ::pread(src_fd, cbuf_in.data(), compressed_sz,
                                 off_t(ref.lhb_shard_hdr_offset) + sizeof(hdr_bytes));
                    if (nr != ssize_t(compressed_sz))
                        throw std::runtime_error("pread payload failed: " + input_paths[ref.source_file_idx]);

                    raw_block.resize(raw_sz);
                    size_t rz = ZSTD_decompressDCtx(zc.d, raw_block.data(), raw_sz,
                                                    cbuf_in.data(), compressed_sz);
                    if (ZSTD_isError(rz))
                        throw std::runtime_error(std::string("zstd decompress: ") + ZSTD_getErrorName(rz));

                    // Parse local contig dict (the unitig_id strings for this block).
                    const uint8_t* p   = raw_block.data();
                    const uint8_t* end = p + rz;
                    uint32_t n_contigs = 0;
                    int n = read_varint(p, end, &n_contigs);
                    if (!n) throw std::runtime_error("corrupt contig dict in block for HOG " + hog_id);
                    p += n;
                    if (n_contigs > 65536)
                        throw std::runtime_error("n_contigs OOB for HOG " + hog_id);
                    std::vector<uint32_t> contig_cnum(n_contigs);
                    for (uint32_t j = 0; j < n_contigs; ++j) {
                        uint32_t len = 0;
                        n = read_varint(p, end, &len);
                        if (!n || p + n + len > end) throw std::runtime_error("truncated contig dict");
                        p += n;
                        std::string_view uid(reinterpret_cast<const char*>(p), len);
                        p += len;
                        auto us = uid.rfind('_');
                        if (us == std::string_view::npos)
                            throw std::runtime_error("unitig_id has no '_': " + std::string(uid));
                        uint32_t cnum = 0;
                        auto res = std::from_chars(uid.data() + us + 1, uid.data() + uid.size(), cnum);
                        if (res.ec != std::errc{})
                            throw std::runtime_error("unitig_id non-numeric suffix: " + std::string(uid));
                        contig_cnum[j] = cnum;
                    }

                    std::vector<VarNTRecord> recs;
                    if (!deserialize_varnt_block(p, end, recs))
                        throw std::runtime_error("corrupt varnt block for HOG " + hog_id);

                    contributing_accs.insert(acc_idx);
                    for (const auto& r : recs) {
                        if (r.contig_idx >= contig_cnum.size())
                            throw std::runtime_error("contig_idx OOB in HOG " + hog_id);
                        uint32_t u_idx = unitig_idx_of(acc_idx, contig_cnum[r.contig_idx]);
                        for (const auto& o : r.vars) {
                            uint32_t hog_pos   = r.sstart + o.hog_offset;
                            uint8_t  codon_idx = uint8_t(o.packed_codon >> 2);
                            inverted[hog_pos].push_back({acc_idx, codon_idx, u_idx});
                        }
                    }
                } else {
                    // .lhg source: decode the pre-inverted HOG block and re-invert
                    // into the global acc space.
                    HogIndexEntry e;
                    e.hog_id = ref.hog_id;
                    e.data_offset = ref.lhg_data_offset;
                    e.data_length = ref.lhg_data_length;
                    InvBlock blk;
                    if (!read_hog_inverted(input_paths[ref.source_file_idx], e, blk))
                        throw std::runtime_error("failed to read .lhg HOG " + hog_id);

                    // Per-source acc remap: source's local acc registry → global.
                    std::vector<uint32_t> remap(ref.acc_ids.size());
                    for (size_t i = 0; i < ref.acc_ids.size(); ++i)
                        remap[i] = global_acc(ref.acc_ids[i]);

                    const uint8_t* p   = blk.pos_ptr;
                    const uint8_t* end = blk.end;
                    uint32_t n_positions = 0;
                    int n = read_varint(p, end, &n_positions);
                    if (!n) throw std::runtime_error("corrupt position count for HOG " + hog_id);
                    p += n;
                    std::vector<InvObs> obs;
                    uint32_t prev_pos = 0;
                    for (uint32_t pi = 0; pi < n_positions; ++pi) {
                        uint32_t hog_pos = 0;
                        if (!decode_position(p, end, prev_pos, hog_pos, obs))
                            throw std::runtime_error("corrupt position record for HOG " + hog_id);
                        auto& dst = inverted[hog_pos];
                        for (const auto& o : obs) {
                            if (o.acc_idx >= remap.size())
                                throw std::runtime_error("acc_idx out of source registry for HOG " + hog_id);
                            uint32_t gacc = remap[o.acc_idx];
                            if (o.unitig_idx >= blk.unitigs.size())
                                throw std::runtime_error("unitig_idx OOB for HOG " + hog_id);
                            uint32_t cnum = blk.unitigs[o.unitig_idx];
                            contributing_accs.insert(gacc);
                            dst.push_back({gacc, o.codon_idx, unitig_idx_of(gacc, cnum)});
                        }
                    }
                }
            }
        } // src closes here

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

        // HOG-level zstd compress (reused CCtx; level set once before the loop).
        size_t bound = ZSTD_compressBound(inv_raw.size());
        hog_cbuf.resize(bound);
        size_t csz = ZSTD_compress2(zc.c, hog_cbuf.data(), bound,
                                    inv_raw.data(), inv_raw.size());
        if (ZSTD_isError(csz))
            throw std::runtime_error(std::string("zstd HOG compress: ") + ZSTD_getErrorName(csz));

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
        uint32_t n_accs = uint32_t(contributing_accs.size());

        index_entries.push_back({hog_id, hog_data_offset, hog_data_length, n_accs});
        group_start = group_end;
    }

    uint64_t index_offset = write_pos;

    // Pass 3: write index + accession registry; backfill file header.
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
