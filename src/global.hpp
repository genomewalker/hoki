#pragma once
#include "format.hpp"
#include "aa.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <zstd.h>
#include <stdexcept>
#include <cstring>
#include <cstdio>
#include <optional>

// .lhg — LHG Global file.
// .lhgi — LHG Index: companion small file; load once, range-GET .lhg per HOG.
//
// .lhg layout:
//   "LHGG" (4)
//   version (1) = LHG_VERSION
//   flags   (1) = 0
//   pad[2]
//   index_offset (8)  ← byte offset of "LHGI" section from file start
//
//   Data section (HOGs sorted lexicographically):
//     For each HOG entry:
//       "LHHE" (4) magic
//       stored_sz (4)  ← uint32 LE; bit31=raw · bit30=ZSTD
//       payload bytes
//
//   Decompressed HOG payload — LHG_VERSION 8 (current):
//     [local acc dict + pident]
//       varint n_local
//       for each: varint delta_gacc, uint8 pident_u8
//     [coverage dict]
//       varint hog_length
//       n_local × varint covered_aa
//     [position headers]
//       varint n_positions
//       n_positions × (varint pos_delta, varint n_obs)
//     [acc column — all positions concatenated]
//       total_obs × varint delta_local_acc_idx
//     [cnum column — all positions concatenated]
//       total_obs × varint cnum
//     [codon column — per position]
//       for each position:
//         uint8 consensus (or 0xFF = no consensus)
//         if != 0xFF: varint n_var, n_var×varint(delta_obs_ordinal), n_var×uint8(codon)
//         else: n_obs × uint8 codon
//
//   Index section (at index_offset):
//     "LHGI" (4)
//     n_hogs (4)
//     for each HOG (sorted):
//       varint(len) + hog_id bytes
//       data_offset (8)
//       data_length (8)
//       n_accessions (4)  ← distinct accessions contributing to this HOG
//     "LHGA" (4)  ← accession registry
//     n_accs (4)
//     for each acc (sorted, so acc_idx = position):
//       varint(len) + acc_id bytes
//
// .lhgi = standalone copy of the index section (same bytes, starts at "LHGI").
// S3 pattern: download .lhgi once; range-GET .lhg[offset:offset+length] per HOG.

namespace lhi {

constexpr uint8_t LHG_FILE_MAGIC[4]      = {'L','H','G','G'};
constexpr uint8_t LHG_INDEX_MAGIC[4]     = {'L','H','G','I'};
constexpr uint8_t LHG_ACC_MAGIC[4]       = {'L','H','G','A'};
constexpr uint8_t LHG_HOG_ENTRY_MAGIC[4] = {'L','H','H','E'};
constexpr uint8_t LHG_VERSION            = 8;
constexpr size_t  LHG_HEADER_SZ          = 16; // magic(4)+ver(1)+flags(1)+pad(2)+index_offset(8)

struct HogIndexEntry {
    std::string hog_id;
    uint64_t    data_offset   = 0;
    uint64_t    data_length   = 0;
    uint32_t    n_accessions  = 0;  // distinct accessions with ≥1 observation
};

struct GlobalIndex {
    std::vector<HogIndexEntry> entries;
    std::vector<std::string>   accessions;  // acc_idx → acc_id string

    const HogIndexEntry* find(const std::string& hog_id) const {
        auto it = std::lower_bound(entries.begin(), entries.end(), hog_id,
            [](const HogIndexEntry& e, const std::string& k) { return e.hog_id < k; });
        return (it != entries.end() && it->hog_id == hog_id) ? &*it : nullptr;
    }

    bool load(const std::string& path) {
        UniqueFd fd(open(path.c_str(), O_RDONLY));
        if (fd < 0) return false;
        return load_from_fd(fd);
    }

    // Load from .lhg by seeking to index_offset stored in the file header.
    bool load_from_lhg(const std::string& path) {
        UniqueFd fd(open(path.c_str(), O_RDONLY));
        if (fd < 0) return false;
        uint8_t hdr[LHG_HEADER_SZ];
        if (::read(fd, hdr, LHG_HEADER_SZ) != ssize_t(LHG_HEADER_SZ) ||
            memcmp(hdr, LHG_FILE_MAGIC, 4) != 0)
            return false;
        if (hdr[4] != LHG_VERSION) return false;
        uint64_t idx_off = 0;
        for (int i = 0; i < 8; ++i) idx_off |= uint64_t(hdr[8 + i]) << (8 * i);
        if (lseek(fd, off_t(idx_off), SEEK_SET) < 0) return false;
        return load_from_fd(fd);
    }

private:
    bool load_from_fd(int fd) {
        // Accumulate every consumed index byte so we can validate the trailing Adler-32.
        std::vector<uint8_t> consumed;
        auto read_exact = [&](void* buf, size_t n) -> bool {
            char* p = (char*)buf; size_t done = 0;
            while (done < n) { ssize_t r = ::read(fd, p + done, n - done); if (r <= 0) return false; done += r; }
            consumed.insert(consumed.end(), (uint8_t*)buf, (uint8_t*)buf + n);
            return true;
        };
        auto read_varint_fd = [&](uint32_t& out) -> bool {
            out = 0; int sh = 0;
            for (;;) {
                uint8_t x; if (::read(fd, &x, 1) != 1) return false;
                consumed.push_back(x);
                if (sh == 28) {                          // 5th byte: only low 4 bits fit uint32
                    if (x & 0x80) return false;          // 6th continuation byte ⇒ overflow
                    if (x > 0x0F) return false;          // value > UINT32_MAX
                    out |= uint32_t(x) << 28; return true;
                }
                out |= uint32_t(x & 0x7F) << sh; sh += 7;
                if (!(x & 0x80)) break;
            }
            return true;
        };

        uint8_t magic[4]; uint32_t n_hogs = 0;
        if (!read_exact(magic, 4) || memcmp(magic, LHG_INDEX_MAGIC, 4) != 0) return false;
        if (!read_exact(&n_hogs, 4)) return false;
        if (n_hogs > 50000000) return false;
        entries.resize(n_hogs);
        for (auto& e : entries) {
            uint32_t len = 0;
            if (!read_varint_fd(len)) return false;
            e.hog_id.resize(len);
            if (!read_exact(e.hog_id.data(), len)) return false;
            uint64_t lo = 0, hi = 0;
            uint32_t na = 0;
            if (!read_exact(&lo, 8) || !read_exact(&hi, 8) || !read_exact(&na, 4)) return false;
            e.data_offset = lo; e.data_length = hi; e.n_accessions = na;
        }

        // Accession registry section.
        uint8_t amagic[4]; uint32_t n_accs = 0;
        if (!read_exact(amagic, 4) || memcmp(amagic, LHG_ACC_MAGIC, 4) != 0) return false;
        if (!read_exact(&n_accs, 4)) return false;
        if (n_accs > 50000000) return false;
        accessions.resize(n_accs);
        for (auto& a : accessions) {
            uint32_t len = 0;
            if (!read_varint_fd(len)) return false;
            a.resize(len);
            if (!read_exact(a.data(), len)) return false;
        }
        // Validate trailing Adler-32 over all preceding index bytes (raw read, not checksummed).
        uint8_t crc_bytes[4]; size_t got = 0;
        while (got < 4) { ssize_t r = ::read(fd, crc_bytes + got, 4 - got); if (r <= 0) return false; got += size_t(r); }
        if (adler32(consumed.data(), consumed.size()) != read_u32_le(crc_bytes)) return false;
        return true;
    }
};

// ── In-memory model ───────────────────────────────────────────────────────────

struct InvObs {
    uint32_t acc_idx;    // local acc idx in InvBlock::local_accs (v7) or global (write path)
    uint8_t  codon_idx; // packed_codon >> 2 (0–63)
    uint32_t cnum;      // direct contig number (no unitig dict)
};

// One decoded position: its observation list (sorted by acc_idx ascending).
struct InvPosition {
    uint32_t            hog_pos;
    std::vector<InvObs> obs;
};

// Serialize an inverted HOG block (v8 format). scratch_* are reusable buffers (cleared on entry).
inline void serialize_inverted_block(std::vector<uint8_t>& raw,
                                     const std::vector<uint32_t>& local_accs,
                                     const std::vector<uint8_t>& local_acc_pident,
                                     const std::vector<InvPosition>& positions,
                                     uint32_t hog_length,
                                     const std::vector<uint32_t>& covered_aa,
                                     std::vector<uint8_t>& hdr_buf,
                                     std::vector<uint8_t>& acc_buf,
                                     std::vector<uint8_t>& cnum_buf,
                                     std::vector<uint8_t>& codon_buf) {
    uint32_t n_local = uint32_t(local_accs.size());

    // Local acc dict + pident: varint(n_local), then per acc: varint(delta_gacc), uint8(pident)
    write_varint(raw, n_local);
    uint32_t prev_gacc = 0;
    for (uint32_t li = 0; li < n_local; ++li) {
        write_varint(raw, local_accs[li] - prev_gacc);
        prev_gacc = local_accs[li];
        uint8_t pid = (li < local_acc_pident.size()) ? local_acc_pident[li] : 0;
        raw.push_back(pid);
    }

    // Coverage dict: hog_length + per-local-acc covered_aa
    write_varint(raw, hog_length);
    for (uint32_t li = 0; li < n_local; ++li) {
        uint32_t caa = (li < covered_aa.size()) ? covered_aa[li] : 0;
        write_varint(raw, caa);
    }

    uint32_t n_positions = uint32_t(positions.size());

    hdr_buf.clear(); acc_buf.clear(); cnum_buf.clear(); codon_buf.clear();

    // Position headers
    write_varint(hdr_buf, n_positions);
    uint32_t prev_pos = 0;
    for (const auto& pos : positions) {
        write_varint(hdr_buf, pos.hog_pos - prev_pos);
        prev_pos = pos.hog_pos;
        write_varint(hdr_buf, uint32_t(pos.obs.size()));
    }

    // Acc column (all positions concatenated) + cnum column (all positions concatenated).
    for (uint32_t pi = 0; pi < n_positions; ++pi) {
        const auto& pos = positions[pi];
        uint32_t n_obs = uint32_t(pos.obs.size());

        uint32_t lacc_ptr = 0;
        uint32_t prev_li = 0;
        for (uint32_t i = 0; i < n_obs; ++i) {
            while (local_accs[lacc_ptr] < pos.obs[i].acc_idx) ++lacc_ptr;
            uint32_t li = lacc_ptr;
            write_varint(acc_buf, li - prev_li);
            prev_li = li;
            write_varint(cnum_buf, pos.obs[i].cnum);
        }
    }

    // Codon column: per position, consensus or explicit; variants keyed by obs ordinal.
    for (uint32_t pi = 0; pi < n_positions; ++pi) {
        const auto& pos = positions[pi];
        uint32_t n_obs = uint32_t(pos.obs.size());

        uint32_t counts[64] = {};
        for (const auto& o : pos.obs) counts[o.codon_idx & 0x3F]++;
        uint8_t best = 0; uint32_t best_cnt = 0;
        for (int c = 0; c < 64; ++c) if (counts[c] > best_cnt) { best_cnt = counts[c]; best = uint8_t(c); }
        uint32_t n_var = n_obs - best_cnt;

        codon_buf.push_back(best);
        write_varint(codon_buf, n_var);
        uint32_t prev_ord = 0;
        for (uint32_t i = 0; i < n_obs; ++i) {
            if (pos.obs[i].codon_idx != best) {
                write_varint(codon_buf, i - prev_ord);
                prev_ord = i;
            }
        }
        for (const auto& o : pos.obs) if (o.codon_idx != best) codon_buf.push_back(o.codon_idx);
    }

    // Concatenate all sections
    raw.insert(raw.end(), hdr_buf.begin(),   hdr_buf.end());
    raw.insert(raw.end(), acc_buf.begin(),   acc_buf.end());
    raw.insert(raw.end(), cnum_buf.begin(),  cnum_buf.end());
    raw.insert(raw.end(), codon_buf.begin(), codon_buf.end());
}


// Read exactly n bytes from fd; returns false on short read.
inline bool fd_read_exact(int fd, void* buf, size_t n) {
    char* p = (char*)buf; size_t done = 0;
    while (done < n) {
        ssize_t r = ::read(fd, p + done, n - done);
        if (r <= 0) return false;
        done += r;
    }
    return true;
}

// Read + decompress one HOG entry's inverted payload into raw bytes.
// Parses local acc dict (with pident) and coverage dict.
// pos_ptr/end span from n_positions varint to end of decompressed payload.
struct InvBlock {
    std::vector<uint8_t>  raw;              // owns decompressed payload
    std::vector<uint32_t> local_accs;      // HOG-local global acc_idx list
    std::vector<uint8_t>  local_acc_pident; // min pident per local acc
    uint32_t              hog_length = 0;  // protein length in AA (0 = unknown)
    std::vector<uint32_t> covered_aa;      // AA covered per local acc
    const uint8_t*        pos_ptr = nullptr;
    const uint8_t*        end     = nullptr;
};

// Decompress raw HOG payload into out.raw and parse the block header (local acc dict
// with pident, coverage dict). Leaves pos_ptr at the n_positions varint.
// Thread-safe: uses pread (no shared fd seek state).
inline bool read_hog_inverted_fd(int fd, const std::string& lhg_path,
                                 const HogIndexEntry& entry, InvBlock& out) {
    (void)lhg_path;
    uint8_t hoe_buf[8];
    if (::pread(fd, hoe_buf, 8, off_t(entry.data_offset)) != 8) return false;
    if (memcmp(hoe_buf, LHG_HOG_ENTRY_MAGIC, 4) != 0)
        throw std::runtime_error("bad HOG entry magic for " + entry.hog_id);
    uint32_t stored = read_u32_le(hoe_buf + 4);
    bool is_raw  = (stored >> 31) & 1;
    bool is_zstd = !is_raw && ((stored >> 30) & 1);
    // Size comes from the uint64 index, not the inline word: a merged super-HOG block
    // can exceed the inline field's 30-bit size range (~1 GiB), at which point the size
    // bits collide with the raw/zstd flag bits. data_length (= 8-byte header + payload)
    // always carries the true length, so flags come from the inline word, size from here.
    if (entry.data_length < 8) return false;
    size_t payload_sz = size_t(entry.data_length - 8);

    std::vector<uint8_t> cbuf_owned(payload_sz);
    { size_t got = 0;                                   // pread caps at ~2 GiB/call; loop for big blocks
      while (got < payload_sz) {
          ssize_t r = ::pread(fd, cbuf_owned.data() + got, payload_sz - got,
                              off_t(entry.data_offset) + 8 + off_t(got));
          if (r <= 0) return false;
          got += size_t(r);
      } }
    const uint8_t* cbuf = cbuf_owned.data();

    if (is_raw) {
        out.raw.assign(cbuf, cbuf + payload_sz);
    } else if (is_zstd) {
        unsigned long long rsz = ZSTD_getFrameContentSize(cbuf, payload_sz);
        if (rsz == ZSTD_CONTENTSIZE_ERROR || rsz == ZSTD_CONTENTSIZE_UNKNOWN)
            throw std::runtime_error("cannot determine ZSTD raw size for HOG " + entry.hog_id);
        out.raw.resize(size_t(rsz));
        size_t dz = ZSTD_decompress(out.raw.data(), out.raw.size(), cbuf, payload_sz);
        if (ZSTD_isError(dz))
            throw std::runtime_error(std::string("zstd HOG decompress: ") + ZSTD_getErrorName(dz));
        out.raw.resize(dz);
    } else {
        throw std::runtime_error("unsupported compression in HOG entry: " + entry.hog_id);
    }

    const uint8_t* p = out.raw.data();
    out.end = p + out.raw.size();
    int n;

    // Local acc dict + pident
    uint32_t n_local = 0;
    n = read_varint(p, out.end, &n_local);
    if (!n) throw std::runtime_error("corrupt local acc dict for HOG " + entry.hog_id);
    p += n;
    // n_local is bounded only by the global accession count; allow the full uint32 range.
    // ceiling: a HOG cannot reference more than UINT32_MAX accessions (acc_idx is uint32);
    // a varint that decodes past that is genuine corruption, already rejected by read_varint.
    out.local_accs.resize(n_local);
    out.local_acc_pident.resize(n_local);
    {
        uint32_t prev = 0;
        for (uint32_t i = 0; i < n_local; ++i) {
            uint32_t d = 0;
            n = read_varint(p, out.end, &d);
            if (!n) throw std::runtime_error("truncated local acc dict");
            p += n;
            prev += d;
            out.local_accs[i] = prev;
            if (p >= out.end) throw std::runtime_error("truncated pident in acc dict");
            out.local_acc_pident[i] = *p++;
        }
    }

    // Coverage dict
    {
        uint32_t hl = 0;
        n = read_varint(p, out.end, &hl);
        if (!n) throw std::runtime_error("corrupt coverage dict for HOG " + entry.hog_id);
        p += n;
        out.hog_length = hl;
        out.covered_aa.resize(n_local);
        for (uint32_t i = 0; i < n_local; ++i) {
            uint32_t caa = 0;
            n = read_varint(p, out.end, &caa);
            if (!n) throw std::runtime_error("truncated coverage dict for HOG " + entry.hog_id);
            p += n;
            out.covered_aa[i] = caa;
        }
    }

    out.pos_ptr = p;
    return true;
}

// Open-and-read variant for query path (saav/freq). Delegates to read_hog_inverted_fd.
inline bool read_hog_inverted(const std::string& lhg_path,
                              const HogIndexEntry& entry,
                              InvBlock& out) {
    UniqueFd fd(open(lhg_path.c_str(), O_RDONLY));
    if (fd < 0) throw std::runtime_error("cannot open: " + lhg_path);
    return read_hog_inverted_fd(fd, lhg_path, entry, out);
}

// Decode all positions in a v8 InvBlock.
// Returns obs with acc_idx = LOCAL acc idx (caller resolves via blk.local_accs).
inline std::vector<InvPosition> decode_block(const InvBlock& blk) {
    const uint8_t* p   = blk.pos_ptr;
    const uint8_t* end = blk.end;
    int n;

    // Position headers
    uint32_t n_positions = 0;
    n = read_varint(p, end, &n_positions);
    if (!n) throw std::runtime_error("decode_block: corrupt position count");
    p += n;

    std::vector<uint32_t> pos_vals(n_positions), n_obs_vec(n_positions);
    uint32_t total_obs = 0;
    uint32_t prev_pos = 0;
    for (uint32_t pi = 0; pi < n_positions; ++pi) {
        uint32_t delta = 0;
        n = read_varint(p, end, &delta); if (!n) throw std::runtime_error("decode_block: corrupt pos delta");
        p += n;
        pos_vals[pi] = prev_pos + delta;
        prev_pos = pos_vals[pi];
        uint32_t nob = 0;
        n = read_varint(p, end, &nob); if (!n) throw std::runtime_error("decode_block: corrupt n_obs");
        p += n;
        n_obs_vec[pi] = nob;
        total_obs += nob;
    }

    // Acc column: delta-encoded local_acc_idx, delta resets per position.
    std::vector<uint32_t> all_lidx(total_obs);
    {
        uint32_t obs_i = 0;
        for (uint32_t pi = 0; pi < n_positions; ++pi) {
            uint32_t prev_li = 0;
            for (uint32_t k = 0; k < n_obs_vec[pi]; ++k) {
                uint32_t d = 0;
                n = read_varint(p, end, &d); if (!n) throw std::runtime_error("decode_block: corrupt acc column");
                p += n;
                prev_li += d;
                all_lidx[obs_i++] = prev_li;
            }
        }
    }

    // Build InvPosition array: fill acc_idx; cnum filled in cnum column pass below.
    std::vector<InvPosition> result(n_positions);
    {
        uint32_t obs_off = 0;
        for (uint32_t pi = 0; pi < n_positions; ++pi) {
            uint32_t nob = n_obs_vec[pi];
            result[pi].hog_pos = pos_vals[pi];
            result[pi].obs.resize(nob);
            for (uint32_t i = 0; i < nob; ++i)
                result[pi].obs[i].acc_idx = all_lidx[obs_off + i];
            obs_off += nob;
        }
    }

    // Cnum column (per-obs, global).
    for (uint32_t pi = 0; pi < n_positions; ++pi) {
        for (uint32_t i = 0; i < n_obs_vec[pi]; ++i) {
            uint32_t c = 0;
            n = read_varint(p, end, &c);
            if (!n) throw std::runtime_error("decode_block: corrupt cnum column");
            p += n;
            result[pi].obs[i].cnum = c;
        }
    }

    // Codon column: all positions sequentially, obs-ordinal-keyed variants.
    {
        static thread_local std::vector<uint32_t> tl_var_ords;
        for (uint32_t pi = 0; pi < n_positions; ++pi) {
            uint32_t nob = n_obs_vec[pi];
            if (p >= end) throw std::runtime_error("decode_block: truncated codon block");
            uint8_t consensus = *p++;
            for (uint32_t i = 0; i < nob; ++i) result[pi].obs[i].codon_idx = consensus;
            uint32_t n_var = 0;
            n = read_varint(p, end, &n_var); if (!n) throw std::runtime_error("decode_block: corrupt n_var");
            p += n;
            if (n_var > nob) throw std::runtime_error("decode_block: n_var > n_obs");
            tl_var_ords.resize(n_var);
            uint32_t prev_ord = 0;
            for (uint32_t vi = 0; vi < n_var; ++vi) {
                uint32_t d = 0;
                n = read_varint(p, end, &d); if (!n) throw std::runtime_error("decode_block: corrupt var ordinal");
                p += n;
                prev_ord += d;
                tl_var_ords[vi] = prev_ord;
            }
            if (p + n_var > end) throw std::runtime_error("decode_block: truncated var codons");
            for (uint32_t vi = 0; vi < n_var; ++vi) {
                uint8_t vc = *p++;
                uint32_t ord = tl_var_ords[vi];
                if (ord >= nob) throw std::runtime_error("decode_block: var ordinal OOB");
                result[pi].obs[ord].codon_idx = vc;
            }
        }
    }
    return result;
}

// ── Queries ──────────────────────────────────────────────────────────────────

// Query a specific SAAV: HOG + position + optional AA filter.
// Output TSV: acc_id \t unitig_id \t hog_pos \t obs_aa \t codon \t pident
void query_saav(const std::string& lhg_path, const GlobalIndex& idx,
                const std::string& hog_id, uint32_t pos,
                std::optional<uint8_t> aa_filter = std::nullopt,
                uint8_t min_pident = 0);

// Per-position codon frequency table for a HOG.
// Output TSV: hog_pos \t codon \t obs_aa \t n_accessions
void query_freq(const std::string& lhg_path, const GlobalIndex& idx,
                const std::string& hog_id,
                uint8_t min_pident = 0);

inline void query_saav(const std::string& lhg_path, const GlobalIndex& idx,
                       const std::string& hog_id, uint32_t pos,
                       std::optional<uint8_t> aa_filter, uint8_t min_pident) {
    const auto* entry = idx.find(hog_id);
    if (!entry) {
        std::fprintf(stderr, "HOG %s not found in index\n", hog_id.c_str());
        return;
    }
    InvBlock blk;
    if (!read_hog_inverted(lhg_path, *entry, blk)) return;

    std::printf("acc_id\tunitig_id\thog_pos\tobs_aa\tcodon\tpident\n");

    auto positions = decode_block(blk);
    for (const auto& ip : positions) {
        if (ip.hog_pos != pos) continue;
        for (const auto& o : ip.obs) {
            uint32_t li = o.acc_idx;
            if (li >= blk.local_accs.size()) continue;
            uint8_t pid = blk.local_acc_pident[li];
            if (pid < min_pident) continue;
            uint32_t gacc = blk.local_accs[li];
            uint8_t packed = uint8_t(o.codon_idx << 2);
            uint8_t aa = codon_to_aa(packed);
            if (aa_filter && aa != *aa_filter) continue;
            char aac = (aa < 20) ? AA_ALPHA[aa] : 'X';
            char cdn[3]; unpack_codon(packed, cdn);
            const std::string& acc = (gacc < idx.accessions.size())
                                   ? idx.accessions[gacc] : std::string();
            std::string uni = acc + "_" + std::to_string(o.cnum);
            std::printf("%s\t%s\t%u\t%c\t%c%c%c\t%u\n",
                        acc.c_str(), uni.c_str(), ip.hog_pos, aac,
                        cdn[0], cdn[1], cdn[2], unsigned(pid));
        }
        break;  // positions are unique; stop after the match
    }
}

inline void query_freq(const std::string& lhg_path, const GlobalIndex& idx,
                       const std::string& hog_id, uint8_t min_pident) {
    const auto* entry = idx.find(hog_id);
    if (!entry) {
        std::fprintf(stderr, "HOG %s not found in index\n", hog_id.c_str());
        return;
    }
    InvBlock blk;
    if (!read_hog_inverted(lhg_path, *entry, blk)) return;

    std::printf("hog_pos\tcodon\tobs_aa\tn_accessions\n");

    auto positions = decode_block(blk);
    for (const auto& ip : positions) {
        // Count distinct accessions per codon at this position.
        std::array<uint32_t, 64> count{};
        std::array<uint32_t, 64> last_acc{};
        bool seen[64] = {false};
        for (const auto& o : ip.obs) {
            uint32_t li = o.acc_idx;
            if (li >= blk.local_accs.size()) continue;
            if (blk.local_acc_pident[li] < min_pident) continue;
            uint32_t gacc = blk.local_accs[li];
            uint8_t c = o.codon_idx & 0x3F;
            if (!seen[c] || last_acc[c] != gacc) {
                ++count[c];
                last_acc[c] = gacc;
                seen[c] = true;
            }
        }
        for (uint8_t c = 0; c < 64; ++c) {
            if (count[c] == 0) continue;
            uint8_t packed = uint8_t(c << 2);
            uint8_t aa = codon_to_aa(packed);
            char aac = (aa < 20) ? AA_ALPHA[aa] : 'X';
            char cdn[3]; unpack_codon(packed, cdn);
            std::printf("%u\t%c%c%c\t%c\t%u\n",
                        ip.hog_pos, cdn[0], cdn[1], cdn[2], aac, count[c]);
        }
    }
}

// Serialize index entries + accession registry to bytes (both .lhg and .lhgi).
inline std::vector<uint8_t> build_index_bytes(const std::vector<HogIndexEntry>& entries,
                                              const std::vector<std::string>& accessions) {
    std::vector<uint8_t> buf;
    buf.insert(buf.end(), LHG_INDEX_MAGIC, LHG_INDEX_MAGIC + 4);
    uint32_t n = uint32_t(entries.size());
    for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(n >> (8 * i)));
    // GlobalIndex::find binary-searches, so the serialized order MUST be sorted by
    // hog_id. merge-shard pre-sorts; merge_batches emits in completion order — sort
    // here so every .lhgi is searchable regardless of the caller. (Sort pointers to
    // avoid copying entries; already-sorted input costs near-nothing.)
    std::vector<const HogIndexEntry*> order; order.reserve(entries.size());
    for (auto& e : entries) order.push_back(&e);
    std::sort(order.begin(), order.end(),
              [](const HogIndexEntry* a, const HogIndexEntry* b){ return a->hog_id < b->hog_id; });
    for (auto* ep : order) {
        const auto& e = *ep;
        write_varint(buf, uint32_t(e.hog_id.size()));
        buf.insert(buf.end(), e.hog_id.begin(), e.hog_id.end());
        for (int i = 0; i < 8; ++i) buf.push_back(uint8_t(e.data_offset  >> (8 * i)));
        for (int i = 0; i < 8; ++i) buf.push_back(uint8_t(e.data_length  >> (8 * i)));
        for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.n_accessions >> (8 * i)));
    }
    buf.insert(buf.end(), LHG_ACC_MAGIC, LHG_ACC_MAGIC + 4);
    uint32_t na = uint32_t(accessions.size());
    for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(na >> (8 * i)));
    for (auto& a : accessions) {
        write_varint(buf, uint32_t(a.size()));
        buf.insert(buf.end(), a.begin(), a.end());
    }
    // Trailing Adler-32 over all preceding index bytes (validated on load).
    write_u32(buf, adler32(buf.data(), buf.size()));
    return buf;
}

} // namespace lhi
