#pragma once
#include "format.hpp"
#include "aa.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <fcntl.h>
#include <unistd.h>
#include <zstd.h>
#include <stdexcept>
#include <cstring>
#include <cstdio>
#include <optional>

// .lhg — LHG Global file (LHG_VERSION 4: position-centric inverted format).
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
//       stored_sz (4)  ← uint32 LE; high bit set = raw fallback, low 31 = payload size
//       payload (stored_sz & 0x7FFFFFFF bytes; zstd frame or raw inverted block)
//
//   Decompressed HOG payload (inverted block):
//     varint n_unitigs
//     for each unitig: varint(contig_num)  ← uint32 numeric suffix only
//     varint n_positions
//     for each position (sorted ascending by hog_pos):
//       varint hog_pos (delta from previous; first = absolute)
//       varint n_obs
//       acc_idx   column: n_obs delta-varints (first absolute, then deltas)
//       codon_idx column: n_obs raw bytes (value 0-63 = packed_codon >> 2)
//       unitig_idx column: n_obs varints (index into local unitig dict)
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
constexpr uint8_t LHG_VERSION            = 4;
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
        if (hdr[4] != LHG_VERSION) return false;  // reject unknown version
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
};

// One decoded position: its observation list (sorted by acc_idx ascending).
struct InvPosition {
    uint32_t            hog_pos;
    std::vector<InvObs> obs;
};

// Serialize an inverted HOG block to raw bytes (pre-compression).
// unitigs = per-HOG local dict; positions sorted ascending by hog_pos,
// each position's obs sorted ascending by acc_idx.
inline void serialize_inverted_block(std::vector<uint8_t>& raw,
                                     const std::vector<uint32_t>& unitigs,
                                     const std::vector<InvPosition>& positions) {
    write_varint(raw, uint32_t(unitigs.size()));
    for (const auto& cnum : unitigs) write_varint(raw, cnum);
    write_varint(raw, uint32_t(positions.size()));
    uint32_t prev_pos = 0;
    for (const auto& pos : positions) {
        write_varint(raw, pos.hog_pos - prev_pos);  // delta; first = absolute
        prev_pos = pos.hog_pos;
        write_varint(raw, uint32_t(pos.obs.size()));
        uint32_t prev = 0;
        for (const auto& o : pos.obs) {
            write_varint(raw, o.acc_idx - prev);
            prev = o.acc_idx;
        }
        for (const auto& o : pos.obs) raw.push_back(o.codon_idx);
        for (const auto& o : pos.obs) write_varint(raw, o.unitig_idx);
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
// Parses local unitig dict; leaves `p`/`end` positioned at the n_positions varint.
struct InvBlock {
    std::vector<uint8_t>     raw;       // owns decompressed payload
    std::vector<uint32_t>    unitigs;   // local unitig dict (contig numbers)
    const uint8_t*           pos_ptr;   // → n_positions varint
    const uint8_t*           end;
};

// pread one HOG entry's compressed payload from an already-open .lhg fd.
// Thread-safe: pread carries no seek state. `dctx` may be nullptr → ZSTD_decompress.
inline bool read_hog_inverted_fd(int fd, const std::string& lhg_path,
                                 const HogIndexEntry& entry, InvBlock& out,
                                 ZSTD_DCtx* dctx = nullptr) {
    (void)lhg_path;  // kept for signature symmetry / future diagnostics
    uint8_t hoe[8];
    if (::pread(fd, hoe, 8, off_t(entry.data_offset)) != 8 ||
        memcmp(hoe, LHG_HOG_ENTRY_MAGIC, 4) != 0)
        throw std::runtime_error("bad HOG entry magic for " + entry.hog_id);
    uint32_t stored = read_u32_le(hoe + 4);
    bool is_raw = (stored >> 31) & 1;
    uint32_t payload_sz = stored & 0x7FFFFFFFu;
    if (payload_sz > 256u * 1024 * 1024) return false;
    std::vector<uint8_t> cbuf(payload_sz);
    if (::pread(fd, cbuf.data(), payload_sz, off_t(entry.data_offset) + 8) != ssize_t(payload_sz))
        return false;

    if (is_raw) {
        out.raw.assign(cbuf.begin(), cbuf.end());
    } else {
        unsigned long long rsz = ZSTD_getFrameContentSize(cbuf.data(), payload_sz);
        if (rsz == ZSTD_CONTENTSIZE_ERROR || rsz == ZSTD_CONTENTSIZE_UNKNOWN)
            throw std::runtime_error("cannot determine raw size for HOG " + entry.hog_id);
        out.raw.resize(size_t(rsz));
        size_t dz = dctx
            ? ZSTD_decompressDCtx(dctx, out.raw.data(), out.raw.size(), cbuf.data(), payload_sz)
            : ZSTD_decompress(out.raw.data(), out.raw.size(), cbuf.data(), payload_sz);
        if (ZSTD_isError(dz))
            throw std::runtime_error(std::string("zstd HOG decompress: ") + ZSTD_getErrorName(dz));
        out.raw.resize(dz);
    }

    const uint8_t* p = out.raw.data();
    out.end = p + out.raw.size();
    uint32_t n_unitigs = 0;
    int n = read_varint(p, out.end, &n_unitigs);
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
    out.pos_ptr = p;
    return true;
}

inline bool read_hog_inverted(const std::string& lhg_path,
                              const HogIndexEntry& entry,
                              InvBlock& out) {
    bool is_raw = false;
    uint32_t payload_sz = 0;
    std::vector<uint8_t> cbuf;
    {
        UniqueFd fd(open(lhg_path.c_str(), O_RDONLY));
        if (fd < 0) throw std::runtime_error("cannot open: " + lhg_path);
        uint8_t fhdr[LHG_HEADER_SZ];
        if (::read(fd, fhdr, LHG_HEADER_SZ) != ssize_t(LHG_HEADER_SZ) ||
            memcmp(fhdr, LHG_FILE_MAGIC, 4) != 0 || fhdr[4] != LHG_VERSION)
            throw std::runtime_error("bad or incompatible .lhg: " + lhg_path);
        if (lseek(fd, off_t(entry.data_offset), SEEK_SET) < 0)
            throw std::runtime_error("seek failed in: " + lhg_path);
        uint8_t hoe_magic[4];
        if (!fd_read_exact(fd, hoe_magic, 4) || memcmp(hoe_magic, LHG_HOG_ENTRY_MAGIC, 4) != 0)
            throw std::runtime_error("bad HOG entry magic for " + entry.hog_id);
        uint32_t stored = 0;
        if (!fd_read_exact(fd, &stored, 4)) return false;
        is_raw = (stored >> 31) & 1;
        payload_sz = stored & 0x7FFFFFFFu;
        if (payload_sz > 256u * 1024 * 1024) return false;
        cbuf.resize(payload_sz);
        if (!fd_read_exact(fd, cbuf.data(), payload_sz)) return false;
    }

    if (is_raw) {
        out.raw.assign(cbuf.begin(), cbuf.end());
    } else {
        unsigned long long rsz = ZSTD_getFrameContentSize(cbuf.data(), payload_sz);
        if (rsz == ZSTD_CONTENTSIZE_ERROR || rsz == ZSTD_CONTENTSIZE_UNKNOWN)
            throw std::runtime_error("cannot determine raw size for HOG " + entry.hog_id);
        out.raw.resize(size_t(rsz));
        size_t dz = ZSTD_decompress(out.raw.data(), out.raw.size(), cbuf.data(), payload_sz);
        if (ZSTD_isError(dz))
            throw std::runtime_error(std::string("zstd HOG decompress: ") + ZSTD_getErrorName(dz));
        out.raw.resize(dz);
    }

    const uint8_t* p = out.raw.data();
    out.end = p + out.raw.size();
    uint32_t n_unitigs = 0;
    int n = read_varint(p, out.end, &n_unitigs);
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
    out.pos_ptr = p;
    return true;
}

// Decode one position's observation columns. Advances `p`.
inline bool decode_position(const uint8_t*& p, const uint8_t* end,
                            uint32_t& prev_pos, uint32_t& hog_pos,
                            std::vector<InvObs>& obs) {
    int n;
    uint32_t delta = 0;
    n = read_varint(p, end, &delta); if (!n) return false; p += n;
    hog_pos = prev_pos + delta;
    prev_pos = hog_pos;
    uint32_t n_obs = 0;
    n = read_varint(p, end, &n_obs); if (!n) return false; p += n;
    if (size_t(n_obs) > size_t(end - p)) return false;  // each obs ≥1 byte; bound before resize
    obs.resize(n_obs);
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
    return true;
}

// ── Queries ──────────────────────────────────────────────────────────────────

// Query a specific SAAV: HOG + position + optional AA filter.
// Output TSV: acc_id \t unitig_id \t hog_pos \t obs_aa \t codon
void query_saav(const std::string& lhg_path, const GlobalIndex& idx,
                const std::string& hog_id, uint32_t pos,
                std::optional<uint8_t> aa_filter = std::nullopt);

// Per-position codon frequency table for a HOG.
// Output TSV: hog_pos \t codon \t obs_aa \t n_accessions
void query_freq(const std::string& lhg_path, const GlobalIndex& idx,
                const std::string& hog_id);

inline void query_saav(const std::string& lhg_path, const GlobalIndex& idx,
                       const std::string& hog_id, uint32_t pos,
                       std::optional<uint8_t> aa_filter) {
    const auto* entry = idx.find(hog_id);
    if (!entry) {
        std::fprintf(stderr, "HOG %s not found in index\n", hog_id.c_str());
        return;
    }
    InvBlock blk;
    if (!read_hog_inverted(lhg_path, *entry, blk)) return;

    const uint8_t* p   = blk.pos_ptr;
    const uint8_t* end = blk.end;
    uint32_t n_positions = 0;
    int n = read_varint(p, end, &n_positions);
    if (!n) throw std::runtime_error("corrupt position count for HOG " + hog_id);
    p += n;

    std::printf("acc_id\tunitig_id\thog_pos\tobs_aa\tcodon\n");

    std::vector<InvObs> obs;
    uint32_t prev_pos = 0;
    for (uint32_t pi = 0; pi < n_positions; ++pi) {
        uint32_t hog_pos = 0;
        if (!decode_position(p, end, prev_pos, hog_pos, obs))
            throw std::runtime_error("corrupt position record for HOG " + hog_id);
        if (hog_pos != pos) continue;
        for (const auto& o : obs) {
            if (o.unitig_idx >= blk.unitigs.size()) continue;  // OOB unitig ref ⇒ skip
            uint8_t packed = uint8_t(o.codon_idx << 2);
            uint8_t aa = codon_to_aa(packed);
            if (aa_filter && aa != *aa_filter) continue;
            char aac = (aa < 20) ? AA_ALPHA[aa] : 'X';
            char cdn[3]; unpack_codon(packed, cdn);
            const std::string& acc = (o.acc_idx < idx.accessions.size())
                                   ? idx.accessions[o.acc_idx] : std::string();
            std::string uni = acc + "_" + std::to_string(blk.unitigs[o.unitig_idx]);
            std::printf("%s\t%s\t%u\t%c\t%c%c%c\n",
                        acc.c_str(), uni.c_str(), hog_pos, aac,
                        cdn[0], cdn[1], cdn[2]);
        }
        break;  // positions are unique; stop after the match
    }
}

inline void query_freq(const std::string& lhg_path, const GlobalIndex& idx,
                       const std::string& hog_id) {
    const auto* entry = idx.find(hog_id);
    if (!entry) {
        std::fprintf(stderr, "HOG %s not found in index\n", hog_id.c_str());
        return;
    }
    InvBlock blk;
    if (!read_hog_inverted(lhg_path, *entry, blk)) return;

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
        if (!decode_position(p, end, prev_pos, hog_pos, obs))
            throw std::runtime_error("corrupt position record for HOG " + hog_id);
        // Count distinct accessions per codon at this position.
        // obs is sorted by acc_idx, but multiple unitigs may repeat an acc per
        // codon; count distinct acc_idx per codon_idx (0–63).
        std::array<uint32_t, 64> count{};
        std::array<uint32_t, 64> last_acc{};
        bool seen[64] = {false};
        for (const auto& o : obs) {
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
