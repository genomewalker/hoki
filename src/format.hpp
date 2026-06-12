#pragma once
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// LHI v2 — Logan HOG Index, lossless bio format
//
// Dropped from original diamond blastx TSV: full_qseq (raw nt contig, ~80% of size)
// Kept exact: qseqid, qstart, qend, qlen, qstrand, sseqid, sstart, send, slen,
//             pident (float32), evalue (float64), qframe, AA array from qseq_translated
//
// Diamond blastx outfmt 6 columns stored:
//   qseqid qstart qend qlen qstrand sseqid sstart send slen pident evalue cigar qseq_translated full_qseq
//   [0]    [1]    [2]  [3]  [4]     [5]    [6]    [7]  [8]  [9]    [10]   [11]  [12]            [13:DROPPED]

namespace lhi {

constexpr uint8_t  MAGIC[8]        = {'L','H','I','v','2','\0','\0','\0'};
constexpr uint16_t FORMAT_VERSION  = 2;

enum class BlockType : uint8_t {
    HOGDict    = 0,  // HOG ID strings, zstd-compressed newline-delimited
    ContigDict = 1,  // full qseqid strings (contig IDs from Logan assembly)
    Alignments = 2,  // packed alignment records, zstd-compressed
    Pileup     = 3,  // materialized per-(hog,pos) AA freq table
};

#pragma pack(push, 1)

struct FileHeader {
    uint8_t  magic[8];
    uint16_t version;
    uint16_t flags;       // bit0=blocks_sorted_by_hog
    uint32_t n_blocks;
    uint8_t  reserved[48];
};
static_assert(sizeof(FileHeader) == 64);

// Block header: 64 bytes. Min/max evalue stored as float32 log10 for skip logic.
// Exact per-record evalue is float64 inside the compressed payload.
struct BlockHeader {
    uint8_t  block_type;
    uint8_t  flags;
    uint32_t compressed_sz;
    uint32_t raw_sz;
    uint32_t n_records;
    uint32_t min_hog_idx;
    uint32_t max_hog_idx;
    uint32_t min_sstart;
    uint32_t max_send;
    float    min_pident;       // exact float32
    float    max_pident;
    float    min_evalue_log;   // log10(evalue), worst (largest) in block
    float    max_evalue_log;   // log10(evalue), best (smallest) in block
    uint32_t crc32c;
    uint8_t  reserved[14];
};
static_assert(sizeof(BlockHeader) == 64);

#pragma pack(pop)

// ── Alignment record (in-memory) ─────────────────────────────────────────────
// On-disk layout inside compressed block (little-endian):
//
//  varint   hog_idx       dict index → HOG ID string
//  varint   contig_idx    dict index → full qseqid
//  uint24   sstart        1-based HOG protein position
//  uint24   send
//  uint24   slen          HOG protein total length
//  uint32   qstart        1-based nt position in contig
//  uint32   qend
//  uint32   qlen          contig total nt length
//  float32  pident        exact percent identity
//  float64  evalue        exact e-value (IEEE 754 double)
//  int8     qframe        reading frame ±1..±3 (sign = strand), 0=unknown
//  [send-sstart+1 bytes]  AA array: one byte per HOG position (encode_aa output)

struct AlignmentRecord {
    uint32_t hog_idx;
    uint32_t contig_idx;
    uint32_t sstart;     // 1-based
    uint32_t send;
    uint32_t slen;
    uint32_t qstart;     // 1-based nt
    uint32_t qend;
    uint32_t qlen;       // contig length in nt
    float    pident;
    double   evalue;
    int8_t   qframe;
    std::vector<uint8_t> aas;  // length = send - sstart + 1

    uint32_t aa_span() const { return send - sstart + 1; }
};

// ── Varint (unsigned LEB128) ──────────────────────────────────────────────────

inline void write_varint(std::vector<uint8_t>& buf, uint32_t v) {
    while (v >= 0x80) {
        buf.push_back(static_cast<uint8_t>((v & 0x7F) | 0x80));
        v >>= 7;
    }
    buf.push_back(static_cast<uint8_t>(v));
}

inline int read_varint(const uint8_t* p, const uint8_t* end, uint32_t* out) {
    uint32_t v = 0; int shift = 0;
    const uint8_t* s = p;
    while (p < end) {
        uint8_t b = *p++;
        v |= static_cast<uint32_t>(b & 0x7F) << shift;
        shift += 7;
        if (!(b & 0x80)) break;
    }
    *out = v;
    return static_cast<int>(p - s);
}

// ── uint24 / uint32 little-endian ────────────────────────────────────────────

inline void write_u24(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back(static_cast<uint8_t>(v));
    buf.push_back(static_cast<uint8_t>(v >> 8));
    buf.push_back(static_cast<uint8_t>(v >> 16));
}
inline uint32_t read_u24(const uint8_t* p) {
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16);
}

inline void write_u32(std::vector<uint8_t>& buf, uint32_t v) {
    buf.push_back(static_cast<uint8_t>(v));
    buf.push_back(static_cast<uint8_t>(v >> 8));
    buf.push_back(static_cast<uint8_t>(v >> 16));
    buf.push_back(static_cast<uint8_t>(v >> 24));
}
inline uint32_t read_u32(const uint8_t* p) {
    return uint32_t(p[0]) | (uint32_t(p[1])<<8) | (uint32_t(p[2])<<16) | (uint32_t(p[3])<<24);
}

// ── IEEE 754 float/double LE ──────────────────────────────────────────────────

inline void write_f32(std::vector<uint8_t>& buf, float v) {
    uint32_t bits; std::memcpy(&bits, &v, 4); write_u32(buf, bits);
}
inline float read_f32(const uint8_t* p) {
    uint32_t bits = read_u32(p); float v; std::memcpy(&v, &bits, 4); return v;
}

inline void write_f64(std::vector<uint8_t>& buf, double v) {
    uint64_t bits; std::memcpy(&bits, &v, 8);
    for (int i = 0; i < 8; ++i) buf.push_back(static_cast<uint8_t>(bits >> (8*i)));
}
inline double read_f64(const uint8_t* p) {
    uint64_t bits = 0;
    for (int i = 0; i < 8; ++i) bits |= uint64_t(p[i]) << (8*i);
    double v; std::memcpy(&v, &bits, 8); return v;
}

// ── Record serialization ──────────────────────────────────────────────────────

inline void serialize_record(std::vector<uint8_t>& buf, const AlignmentRecord& r) {
    write_varint(buf, r.hog_idx);
    write_varint(buf, r.contig_idx);
    write_u24(buf, r.sstart);
    write_u24(buf, r.send);
    write_u24(buf, r.slen);
    write_u32(buf, r.qstart);
    write_u32(buf, r.qend);
    write_u32(buf, r.qlen);
    write_f32(buf, r.pident);
    write_f64(buf, r.evalue);
    buf.push_back(static_cast<uint8_t>(r.qframe));
    buf.insert(buf.end(), r.aas.begin(), r.aas.end());
}

// returns bytes consumed, 0 on underflow
inline int deserialize_record(const uint8_t* p, const uint8_t* end, AlignmentRecord& r) {
    const uint8_t* s = p;
    int n;
    n = read_varint(p, end, &r.hog_idx);    p += n;
    n = read_varint(p, end, &r.contig_idx); p += n;
    // fixed: 3+3+3+4+4+4+4+8+1 = 34 bytes
    if (p + 34 > end) return 0;
    r.sstart  = read_u24(p); p += 3;
    r.send    = read_u24(p); p += 3;
    r.slen    = read_u24(p); p += 3;
    r.qstart  = read_u32(p); p += 4;
    r.qend    = read_u32(p); p += 4;
    r.qlen    = read_u32(p); p += 4;
    r.pident  = read_f32(p); p += 4;
    r.evalue  = read_f64(p); p += 8;
    r.qframe  = static_cast<int8_t>(*p++);
    uint32_t span = r.send - r.sstart + 1;
    if (p + span > end) return 0;
    r.aas.assign(p, p + span);
    p += span;
    return static_cast<int>(p - s);
}

} // namespace lhi
