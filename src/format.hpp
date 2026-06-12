#pragma once
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <span>

namespace lhi {

// ── File header (64 bytes, fixed) ──────────────────────────────────────────

constexpr uint8_t MAGIC[8] = {'L','H','I','v','1','\0','\0','\0'};
constexpr uint16_t FORMAT_VERSION = 1;

enum class BlockType : uint8_t {
    HOGDict    = 0,  // string dictionary: HOG IDs (zstd-compressed newline-delimited)
    SampleDict = 1,  // string dictionary: SRA accessions
    Alignments = 2,  // variable-length alignment records (zstd-compressed)
    Pileup     = 3,  // per-(hog_idx, position) AA frequency tables
};

#pragma pack(push, 1)

struct FileHeader {
    uint8_t  magic[8];
    uint16_t version;
    uint16_t flags;           // bit0=sorted, bit1=has_cluster_sidecar
    uint32_t n_blocks;
    uint8_t  reserved[48];   // = 64 bytes total
};
static_assert(sizeof(FileHeader) == 64);

struct BlockHeader {
    uint8_t  block_type;     // BlockType
    uint8_t  flags;          // bit0=last_block
    uint32_t compressed_sz;
    uint32_t raw_sz;
    uint32_t n_records;
    uint32_t min_hog_idx;
    uint32_t max_hog_idx;
    uint32_t min_sstart;     // HOG coordinate
    uint32_t max_send;
    uint8_t  min_pident_q;   // quantized pident (×2.55)
    uint8_t  max_pident_q;
    int8_t   min_evalue_q;   // quantized log10(evalue), int8 clamped [-128,0]
    int8_t   max_evalue_q;   // best (most negative) evalue in block
    uint32_t crc32c;
    uint8_t  reserved[26];  // = 64 bytes total
};
static_assert(sizeof(BlockHeader) == 64);

#pragma pack(pop)

// ── Alignment record (variable, in-memory only) ────────────────────────────
// On-disk layout inside a compressed block (all little-endian):
//
//  varint  hog_idx     (dict index)
//  varint  sample_idx  (dict index)
//  uint24  sstart      (1-based HOG position)
//  uint24  send
//  uint8   pident_q
//  int8    evalue_q
//  int8    qframe      (±1..±3, 0=unknown; sign encodes strand)
//  uint8   n_aas_lo    }  n_aas = (send-sstart+1), stored as
//                      }  implicit — no explicit length field;
//                      }  decoder computes from sstart/send
//  [n_aas bytes]       AA array (one byte per HOG position)

struct AlignmentRecord {
    uint32_t hog_idx;
    uint32_t sample_idx;
    uint32_t sstart;   // 1-based
    uint32_t send;
    uint8_t  pident_q;
    int8_t   evalue_q;
    int8_t   qframe;
    std::vector<uint8_t> aas;  // length = send - sstart + 1

    uint32_t aa_span() const { return send - sstart + 1; }
};

// ── Varint codec (unsigned LEB128) ─────────────────────────────────────────

inline void write_varint(std::vector<uint8_t>& buf, uint32_t v) {
    while (v >= 0x80) {
        buf.push_back(static_cast<uint8_t>((v & 0x7F) | 0x80));
        v >>= 7;
    }
    buf.push_back(static_cast<uint8_t>(v));
}

// returns bytes consumed, writes value to *out
inline int read_varint(const uint8_t* p, const uint8_t* end, uint32_t* out) {
    uint32_t v = 0; int shift = 0;
    const uint8_t* start = p;
    while (p < end) {
        uint8_t b = *p++;
        v |= static_cast<uint32_t>(b & 0x7F) << shift;
        shift += 7;
        if (!(b & 0x80)) break;
    }
    *out = v;
    return static_cast<int>(p - start);
}

// ── uint24 codec ──────────────────────────────────────────────────────────

inline void write_u24(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back(static_cast<uint8_t>(v));
    buf.push_back(static_cast<uint8_t>(v >> 8));
    buf.push_back(static_cast<uint8_t>(v >> 16));
}

inline uint32_t read_u24(const uint8_t* p) {
    return static_cast<uint32_t>(p[0])
         | (static_cast<uint32_t>(p[1]) << 8)
         | (static_cast<uint32_t>(p[2]) << 16);
}

// ── Record serialization ──────────────────────────────────────────────────

inline void serialize_record(std::vector<uint8_t>& buf, const AlignmentRecord& r) {
    write_varint(buf, r.hog_idx);
    write_varint(buf, r.sample_idx);
    write_u24(buf, r.sstart);
    write_u24(buf, r.send);
    buf.push_back(r.pident_q);
    buf.push_back(static_cast<uint8_t>(r.evalue_q));
    buf.push_back(static_cast<uint8_t>(r.qframe));
    buf.insert(buf.end(), r.aas.begin(), r.aas.end());
}

// returns bytes consumed (0 on error/end)
inline int deserialize_record(const uint8_t* p, const uint8_t* end, AlignmentRecord& r) {
    const uint8_t* start = p;
    if (p + 12 > end) return 0;  // minimum viable size check

    int n;
    n = read_varint(p, end, &r.hog_idx);    p += n;
    n = read_varint(p, end, &r.sample_idx); p += n;
    if (p + 9 > end) return 0;
    r.sstart   = read_u24(p); p += 3;
    r.send     = read_u24(p); p += 3;
    r.pident_q = *p++;
    r.evalue_q = static_cast<int8_t>(*p++);
    r.qframe   = static_cast<int8_t>(*p++);

    uint32_t span = r.send - r.sstart + 1;
    if (p + span > end) return 0;
    r.aas.assign(p, p + span);
    p += span;
    return static_cast<int>(p - start);
}

} // namespace lhi
