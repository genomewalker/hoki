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
#include <lz4.h>
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
//       stored_sz (4)  ← uint32 LE; bit31=raw · bit30=ZSTD · else=LZ4
//       payload bytes
//
//   Decompressed HOG payload — LHG_VERSION 4:
//     varint n_unitigs; per unitig: varint(contig_num)
//     varint n_positions
//     per position: varint hog_pos_delta · varint n_obs
//       n_obs × varint(delta global_acc_idx)
//       n_obs × uint8(codon_idx)
//       n_obs × varint(unitig_idx)
//
//   Decompressed HOG payload — LHG_VERSION 5:
//     varint n_unitigs; per unitig: varint(contig_num)
//     varint n_hog_accs; per local acc: varint(delta global_acc_idx)   ← HOG-local dict
//     varint n_positions
//     per position: varint hog_pos_delta · varint n_obs
//       n_obs × varint(delta local_acc_idx)
//       n_obs × varint(unitig_idx)
//       uint8 consensus_codon (0xFF = all explicit)
//       if consensus != 0xFF:
//         varint n_var; n_var × varint(delta local_acc_idx); n_var × uint8(codon_idx)
//       else: n_obs × uint8(codon_idx)
//
//   Decompressed HOG payload — LHG_VERSION 6 (this file):
//     varint n_unitigs; per unitig: varint(contig_num)
//     varint n_hog_accs; per local acc: varint(delta global_acc_idx)
//     varint hog_length                              ← protein length in AA (0 = unknown)
//     n_unitigs × varint(covered_aa)                 ← AA positions covered per unitig
//     varint n_positions
//     per position: varint hog_pos_delta · varint n_obs
//       n_obs × varint(delta local_acc_idx)
//       n_obs × varint(unitig_idx)
//       n_obs × uint8(pident_u8)                     ← round(pident); 0 = unknown
//       uint8 consensus_codon (0xFF = all explicit)
//       if consensus != 0xFF:
//         varint n_var; n_var × varint(delta local_acc_idx); n_var × uint8(codon_idx)
//       else: n_obs × uint8(codon_idx)
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
constexpr uint8_t LHG_VERSION            = 6;
constexpr uint8_t LHG_VERSION_MIN        = 4;  // oldest version we can read
constexpr size_t  LHG_HEADER_SZ          = 16; // magic(4)+ver(1)+flags(1)+pad(2)+index_offset(8)

struct HogIndexEntry {
    std::string hog_id;
    uint64_t    data_offset   = 0;
    uint64_t    data_length   = 0;
    uint32_t    n_accessions  = 0;  // distinct accessions with ≥1 observation
};

struct GlobalIndex {
    uint8_t                    file_version = LHG_VERSION;
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
        if (hdr[4] < LHG_VERSION_MIN || hdr[4] > LHG_VERSION) return false;
        file_version = hdr[4];
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

// ── LHG_VERSION 4 inverted-block in-memory model ─────────────────────────────

struct InvObs {
    uint32_t acc_idx;
    uint8_t  codon_idx;   // packed_codon >> 2 (0–63)
    uint32_t unitig_idx;  // index into per-HOG local unitig dict
    uint8_t  pident_u8;   // round(pident); 0 = unknown (v4/v5 source)
};

// One decoded position: its observation list (sorted by acc_idx ascending).
struct InvPosition {
    uint32_t            hog_pos;
    std::vector<InvObs> obs;
};

// Serialize an inverted HOG block (v6 format).
// local_accs: sorted unique global acc_idx values for this HOG (obs remap global→local).
// positions: sorted ascending by hog_pos; each obs sorted ascending by (global) acc_idx.
// hog_length: protein length in AA (0 = unknown, e.g. when merging v5 shards).
// covered_aa: AA positions covered per local unitig (indexed by unitig_idx; 0 = unknown).
inline void serialize_inverted_block(std::vector<uint8_t>& raw,
                                     const std::vector<uint32_t>& unitigs,
                                     const std::vector<uint32_t>& local_accs,
                                     const std::vector<InvPosition>& positions,
                                     uint32_t hog_length = 0,
                                     const std::vector<uint32_t>& covered_aa = {}) {
    // Unitig dict
    write_varint(raw, uint32_t(unitigs.size()));
    for (auto c : unitigs) write_varint(raw, c);

    // Per-HOG local acc dict (delta global acc_idx, sorted)
    write_varint(raw, uint32_t(local_accs.size()));
    uint32_t prev_gacc = 0;
    for (auto g : local_accs) { write_varint(raw, g - prev_gacc); prev_gacc = g; }

    // v6: coverage dict — protein length + per-unitig covered AA count
    write_varint(raw, hog_length);
    for (size_t i = 0; i < unitigs.size(); ++i) {
        uint32_t caa = (i < covered_aa.size()) ? covered_aa[i] : 0;
        write_varint(raw, caa);
    }

    write_varint(raw, uint32_t(positions.size()));
    uint32_t prev_pos = 0;
    for (const auto& pos : positions) {
        write_varint(raw, pos.hog_pos - prev_pos);
        prev_pos = pos.hog_pos;
        uint32_t n_obs = uint32_t(pos.obs.size());
        write_varint(raw, n_obs);

        // Acc column: obs are sorted by acc_idx → forward scan in local_accs
        // (O(n_obs + n_local) total) vs binary search (O(n_obs * log n_local)).
        static thread_local std::vector<uint32_t> tl_lidxs;
        tl_lidxs.resize(n_obs);
        uint32_t lacc_ptr = 0;
        uint32_t prev_li = 0;
        for (uint32_t i = 0; i < n_obs; ++i) {
            while (local_accs[lacc_ptr] < pos.obs[i].acc_idx) ++lacc_ptr;
            uint32_t li = lacc_ptr;
            tl_lidxs[i] = li;
            write_varint(raw, li - prev_li);
            prev_li = li;
        }

        // Unitig column
        for (const auto& o : pos.obs) write_varint(raw, o.unitig_idx);

        // v6: pident column (uint8 per obs; 0 = unknown for v4/v5-sourced data)
        for (const auto& o : pos.obs) raw.push_back(o.pident_u8);

        // Codon column with consensus optimization.
        uint32_t counts[64] = {};
        for (const auto& o : pos.obs) counts[o.codon_idx & 0x3F]++;
        uint8_t best = 0; uint32_t best_cnt = 0;
        for (int c = 0; c < 64; ++c) if (counts[c] > best_cnt) { best_cnt = counts[c]; best = uint8_t(c); }
        uint32_t n_var = n_obs - best_cnt;

        if (n_var == n_obs) {
            // All different — no consensus worth writing
            raw.push_back(0xFF);
            for (const auto& o : pos.obs) raw.push_back(o.codon_idx);
        } else {
            raw.push_back(best);
            write_varint(raw, n_var);
            // Variant local_acc_idx list (sorted, delta-encoded)
            uint32_t prev_vli = 0;
            for (uint32_t i = 0; i < n_obs; ++i) {
                if (pos.obs[i].codon_idx != best) {
                    write_varint(raw, tl_lidxs[i] - prev_vli);
                    prev_vli = tl_lidxs[i];
                }
            }
            // Variant codon values
            for (const auto& o : pos.obs) if (o.codon_idx != best) raw.push_back(o.codon_idx);
        }
    }
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
// Parses local unitig dict (and local acc dict for v5+, coverage dict for v6+).
// pos_ptr/end span from n_positions varint to end of decompressed payload.
struct InvBlock {
    std::vector<uint8_t>     raw;        // owns decompressed payload
    std::vector<uint32_t>    unitigs;    // local unitig dict (contig numbers)
    std::vector<uint32_t>    local_accs; // v5+: HOG-local global acc_idx list (empty = v4)
    uint32_t                 hog_length = 0;   // v6+: protein length in AA (0 = unknown)
    std::vector<uint32_t>    covered_aa;        // v6+: AA covered per local unitig
    const uint8_t*           pos_ptr = nullptr;
    const uint8_t*           end     = nullptr;
};

// Decompress raw HOG payload into out.raw and parse the header fields (unitig dict
// and, for v5, per-HOG local acc list). Leaves pos_ptr at the n_positions varint.
// file_version: pass the .lhg file's version byte (4 or 5).
// Thread-safe: uses pread (no shared fd seek state).
inline bool read_hog_inverted_fd(int fd, const std::string& lhg_path,
                                 const HogIndexEntry& entry, InvBlock& out,
                                 uint8_t file_version = LHG_VERSION) {
    (void)lhg_path;
    uint8_t hoe_buf[8];
    if (::pread(fd, hoe_buf, 8, off_t(entry.data_offset)) != 8) return false;
    if (memcmp(hoe_buf, LHG_HOG_ENTRY_MAGIC, 4) != 0)
        throw std::runtime_error("bad HOG entry magic for " + entry.hog_id);
    uint32_t stored = read_u32_le(hoe_buf + 4);
    bool is_raw  = (stored >> 31) & 1;
    // v4 blocks: bit31=raw, else ZSTD (no LZ4, no bit30 flag).
    // v5+ blocks: bit31=raw, bit30=ZSTD, else LZ4.
    bool is_zstd = !is_raw && (file_version < 5 || ((stored >> 30) & 1));
    uint32_t payload_sz = stored & (is_raw ? 0x7FFFFFFFu : 0x3FFFFFFFu);
    if (payload_sz > 256u * 1024 * 1024) return false;

    std::vector<uint8_t> cbuf_owned(payload_sz);
    if (::pread(fd, cbuf_owned.data(), payload_sz, off_t(entry.data_offset) + 8) != ssize_t(payload_sz))
        return false;
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
        if (payload_sz < 4) throw std::runtime_error("LZ4 payload truncated: " + entry.hog_id);
        uint32_t raw_sz = read_u32_le(cbuf);
        if (raw_sz > 256u * 1024 * 1024) throw std::runtime_error("LZ4 raw_sz OOB: " + entry.hog_id);
        out.raw.resize(raw_sz);
        int dz = LZ4_decompress_safe(
            reinterpret_cast<const char*>(cbuf + 4),
            reinterpret_cast<char*>(out.raw.data()),
            int(payload_sz - 4), int(raw_sz));
        if (dz < 0)
            throw std::runtime_error("LZ4 HOG decompress failed: " + entry.hog_id);
    }

    const uint8_t* p = out.raw.data();
    out.end = p + out.raw.size();
    int n;
    uint32_t n_unitigs = 0;
    n = read_varint(p, out.end, &n_unitigs);
    if (!n) throw std::runtime_error("corrupt unitig dict for HOG " + entry.hog_id);
    p += n;
    if (n_unitigs > 65536) throw std::runtime_error("n_unitigs OOB for HOG " + entry.hog_id);
    out.unitigs.resize(n_unitigs);
    for (uint32_t i = 0; i < n_unitigs; ++i) {
        uint32_t cnum = 0;
        n = read_varint(p, out.end, &cnum);
        if (!n) throw std::runtime_error("truncated unitig dict");
        p += n;
        out.unitigs[i] = cnum;
    }

    // v5+: parse per-HOG local acc dict (global acc_idx values, delta-encoded).
    out.local_accs.clear();
    if (file_version >= 5) {
        uint32_t n_local = 0;
        n = read_varint(p, out.end, &n_local);
        if (!n) throw std::runtime_error("corrupt local acc dict for HOG " + entry.hog_id);
        p += n;
        out.local_accs.resize(n_local);
        uint32_t prev = 0;
        for (uint32_t i = 0; i < n_local; ++i) {
            uint32_t d = 0;
            n = read_varint(p, out.end, &d);
            if (!n) throw std::runtime_error("truncated local acc dict");
            p += n;
            prev += d;
            out.local_accs[i] = prev;
        }
    }

    // v6+: parse coverage dict — hog_length + per-unitig covered_aa.
    out.hog_length = 0; out.covered_aa.clear();
    if (file_version >= 6) {
        uint32_t hl = 0;
        n = read_varint(p, out.end, &hl);
        if (!n) throw std::runtime_error("corrupt coverage dict for HOG " + entry.hog_id);
        p += n;
        out.hog_length = hl;
        out.covered_aa.resize(n_unitigs);
        for (uint32_t i = 0; i < n_unitigs; ++i) {
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
                              InvBlock& out,
                              uint8_t file_version = LHG_VERSION) {
    UniqueFd fd(open(lhg_path.c_str(), O_RDONLY));
    if (fd < 0) throw std::runtime_error("cannot open: " + lhg_path);
    return read_hog_inverted_fd(fd, lhg_path, entry, out, file_version);
}

// Decode one position's observation columns. Advances `p`.
// local_accs: if non-empty (v5+), maps local_acc_idx → global acc_idx.
//             if empty (v4), local_acc_idx IS the (shard-local) global acc_idx.
// file_version: pass the source .lhg version byte to select the correct column layout.
inline bool decode_position(const uint8_t*& p, const uint8_t* end,
                            uint32_t& prev_pos, uint32_t& hog_pos,
                            std::vector<InvObs>& obs,
                            const std::vector<uint32_t>& local_accs = {},
                            uint8_t file_version = LHG_VERSION) {
    int n;
    uint32_t delta = 0;
    n = read_varint(p, end, &delta); if (!n) return false; p += n;
    hog_pos = prev_pos + delta;
    prev_pos = hog_pos;
    uint32_t n_obs = 0;
    n = read_varint(p, end, &n_obs); if (!n) return false; p += n;
    obs.resize(n_obs);

    if (local_accs.empty()) {
        // v4 fast path: no tl_local_idxs, no remapping overhead.
        uint32_t prev = 0;
        for (uint32_t i = 0; i < n_obs; ++i) {
            uint32_t d = 0;
            n = read_varint(p, end, &d); if (!n) return false; p += n;
            prev += d;
            obs[i].acc_idx = prev;
        }
        if (p + n_obs > end) return false;
        for (uint32_t i = 0; i < n_obs; ++i) obs[i].codon_idx = *p++;
        for (uint32_t i = 0; i < n_obs; ++i) {
            uint32_t u = 0;
            n = read_varint(p, end, &u); if (!n) return false; p += n;
            obs[i].unitig_idx = u;
        }
        for (uint32_t i = 0; i < n_obs; ++i) obs[i].pident_u8 = 0;
    } else {
        // v5/v6: acc remap via local dict, unitig before codon, consensus codon encoding.
        static thread_local std::vector<uint32_t> tl_local_idxs;
        tl_local_idxs.resize(n_obs);
        uint32_t prev_local = 0;
        for (uint32_t i = 0; i < n_obs; ++i) {
            uint32_t d = 0;
            n = read_varint(p, end, &d); if (!n) return false; p += n;
            prev_local += d;
            tl_local_idxs[i] = prev_local;
            obs[i].acc_idx = local_accs[prev_local];
        }
        for (uint32_t i = 0; i < n_obs; ++i) {
            uint32_t u = 0;
            n = read_varint(p, end, &u); if (!n) return false; p += n;
            obs[i].unitig_idx = u;
        }
        // v6: pident column; v5 sources have no pident → fill 0
        if (file_version >= 6) {
            if (p + n_obs > end) return false;
            for (uint32_t i = 0; i < n_obs; ++i) obs[i].pident_u8 = *p++;
        } else {
            for (uint32_t i = 0; i < n_obs; ++i) obs[i].pident_u8 = 0;
        }
        if (p >= end) return false;
        uint8_t consensus = *p++;
        if (consensus == 0xFF) {
            if (p + n_obs > end) return false;
            for (uint32_t i = 0; i < n_obs; ++i) obs[i].codon_idx = *p++;
        } else {
            for (uint32_t i = 0; i < n_obs; ++i) obs[i].codon_idx = consensus;
            uint32_t n_var = 0;
            n = read_varint(p, end, &n_var); if (!n) return false; p += n;
            if (n_var > n_obs) return false;
            static thread_local std::vector<uint32_t> tl_var_locals;
            tl_var_locals.resize(n_var);
            uint32_t prev_li = 0;
            for (uint32_t vi = 0; vi < n_var; ++vi) {
                uint32_t d = 0;
                n = read_varint(p, end, &d); if (!n) return false; p += n;
                prev_li += d;
                tl_var_locals[vi] = prev_li;
            }
            if (p + n_var > end) return false;
            uint32_t j = 0;
            for (uint32_t vi = 0; vi < n_var; ++vi) {
                uint8_t vc = *p++;
                while (j < n_obs && tl_local_idxs[j] < tl_var_locals[vi]) ++j;
                if (j >= n_obs || tl_local_idxs[j] != tl_var_locals[vi]) return false;
                obs[j].codon_idx = vc;
                ++j;
            }
        }
    }
    return true;
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
    if (!read_hog_inverted(lhg_path, *entry, blk, idx.file_version)) return;

    const uint8_t* p   = blk.pos_ptr;
    const uint8_t* end = blk.end;
    uint32_t n_positions = 0;
    int n = read_varint(p, end, &n_positions);
    if (!n) throw std::runtime_error("corrupt position count for HOG " + hog_id);
    p += n;

    std::printf("acc_id\tunitig_id\thog_pos\tobs_aa\tcodon\tpident\n");

    std::vector<InvObs> obs;
    uint32_t prev_pos = 0;
    for (uint32_t pi = 0; pi < n_positions; ++pi) {
        uint32_t hog_pos = 0;
        if (!decode_position(p, end, prev_pos, hog_pos, obs, blk.local_accs, idx.file_version))
            throw std::runtime_error("corrupt position record for HOG " + hog_id);
        if (hog_pos != pos) continue;
        for (const auto& o : obs) {
            if (o.unitig_idx >= blk.unitigs.size()) continue;
            if (o.pident_u8 < min_pident) continue;
            uint8_t packed = uint8_t(o.codon_idx << 2);
            uint8_t aa = codon_to_aa(packed);
            if (aa_filter && aa != *aa_filter) continue;
            char aac = (aa < 20) ? AA_ALPHA[aa] : 'X';
            char cdn[3]; unpack_codon(packed, cdn);
            const std::string& acc = (o.acc_idx < idx.accessions.size())
                                   ? idx.accessions[o.acc_idx] : std::string();
            std::string uni = acc + "_" + std::to_string(blk.unitigs[o.unitig_idx]);
            std::printf("%s\t%s\t%u\t%c\t%c%c%c\t%u\n",
                        acc.c_str(), uni.c_str(), hog_pos, aac,
                        cdn[0], cdn[1], cdn[2], unsigned(o.pident_u8));
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
    if (!read_hog_inverted(lhg_path, *entry, blk, idx.file_version)) return;

    const uint8_t* p   = blk.pos_ptr;
    const uint8_t* end = blk.end;
    uint32_t n_positions = 0;
    int n = read_varint(p, end, &n_positions);
    if (!n) throw std::runtime_error("corrupt position count for HOG " + hog_id);
    p += n;

    std::printf("hog_pos\tcodon\tobs_aa\tn_accessions\n");

    std::vector<InvObs> obs;
    uint32_t prev_pos = 0;
    for (uint32_t pi = 0; pi < n_positions; ++pi) {
        uint32_t hog_pos = 0;
        if (!decode_position(p, end, prev_pos, hog_pos, obs, blk.local_accs, idx.file_version))
            throw std::runtime_error("corrupt position record for HOG " + hog_id);
        // Count distinct accessions per codon at this position.
        // obs is sorted by acc_idx, but multiple unitigs may repeat an acc per
        // codon; count distinct acc_idx per codon_idx (0–63).
        std::array<uint32_t, 64> count{};
        std::array<uint32_t, 64> last_acc{};
        bool seen[64] = {false};
        for (const auto& o : obs) {
            if (o.pident_u8 < min_pident) continue;
            uint8_t c = o.codon_idx & 0x3F;
            if (!seen[c] || last_acc[c] != o.acc_idx) {
                ++count[c];
                last_acc[c] = o.acc_idx;
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
                        hog_pos, cdn[0], cdn[1], cdn[2], aac, count[c]);
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
    for (auto& e : entries) {
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
