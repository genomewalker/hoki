#pragma once
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// LHI v3 — Logan HOG Index, fully lossless alignment format
//
// Dropped from diamond blastx TSV: full_qseq (raw nt contig, ~80-85% of size)
//   Pointer-lossless: reconstruct via S3 using accession in qseqid + nt_hash verification
//   Local-lossless: optional NT block type (2-bit+zstd, stored on --store-nt)
//
// CIGAR is fully reconstructible from: aas[] (M/D positions) + inserts[] (I events)
// qseq_translated reconstructible from: aas[], inserts[] (re-encode AA→char, insert at I positions)
//
// Diamond blastx outfmt 6 columns:
//   qseqid qstart qend qlen qstrand sseqid sstart send slen pident evalue cigar qseq_translated full_qseq
//     0      1     2    3     4       5      6      7    8     9     10     11         12            13

namespace lhi {

constexpr uint8_t  MAGIC[8]       = {'L','H','I','v','3','\0','\0','\0'};
constexpr uint16_t FORMAT_VERSION = 3;

// Contract flags (FileHeader.flags)
constexpr uint16_t FLAG_VALUE_LOSSLESS = 0x0001;  // default: typed fields exact, TSV lexemes may differ
constexpr uint16_t FLAG_BYTE_LOSSLESS  = 0x0002;  // exception stream present for TSV-exact round-trip
constexpr uint16_t FLAG_HAS_NT_BLOCKS  = 0x0004;  // NT blocks present (local-lossless for full_qseq)

enum class BlockType : uint8_t {
    HOGDict    = 0,
    ContigDict = 1,
    Alignments = 2,
    Pileup     = 3,
    NT         = 4,
    Variants   = 5,  // sparse: only M-position (hog_pos, obs_aa) pairs per record
    VarNT      = 6,  // variant+codon+tiling block (--varnt mode)
};

#pragma pack(push, 1)

struct FileHeader {
    uint8_t  magic[8];
    uint16_t version;
    uint16_t flags;
    uint32_t n_blocks;
    uint8_t  reserved[48];
};
static_assert(sizeof(FileHeader) == 64);

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
    float    min_pident;
    float    max_pident;
    float    min_evalue_log;   // log10, worst (largest) in block
    float    max_evalue_log;   // log10, best (smallest) in block
    uint32_t crc32c;
    uint8_t  reserved[14];
};
static_assert(sizeof(BlockHeader) == 64);

#pragma pack(pop)

// ── Alignment record (in-memory) ─────────────────────────────────────────────
// On-disk layout inside compressed Alignments block:
//
//  varint   hog_idx
//  varint   contig_idx
//  uint24   sstart          1-based HOG position
//  uint24   send
//  uint32   qstart          1-based nt position in contig
//  uint32   qend
//  uint32   qlen            contig total nt length
//  float32  pident
//  float64  evalue
//  int8     qframe          ±1..±3 (sign=strand, value=frame), 0=unknown
//  uint8    nt_hash[16]     xxh3-128 of full_qseq (pointer-lossless verification)
//  uint8    n_aas = send-sstart+1  (implicit, not stored; decoder computes from sstart/send)
//  [n_aas]  aas[]           AA at each HOG position: M→encode_aa, D→AA_GAP
//  varint   n_inserts
//  per insertion:
//    uint24   before_hog_pos
//    uint8    n_ins_aas
//    [n_ins_aas] ins_aas[]

struct InsertionRecord {
    uint32_t before_hog_pos;
    std::vector<uint8_t> aas;
};

struct AlignmentRecord {
    uint32_t hog_idx;
    uint32_t contig_idx;
    uint32_t sstart, send;
    uint32_t qstart, qend, qlen;
    float    pident;
    double   evalue;
    int8_t   qframe;
    uint8_t  nt_hash[16];           // xxh3-128(full_qseq)
    std::vector<uint8_t>       aas; // length = send-sstart+1
    std::vector<InsertionRecord> inserts;

    uint32_t aa_span() const { return send - sstart + 1; }
};

// ── Varint ────────────────────────────────────────────────────────────────────

inline void write_varint(std::vector<uint8_t>& b, uint32_t v) {
    while (v >= 0x80) { b.push_back(uint8_t((v&0x7F)|0x80)); v>>=7; }
    b.push_back(uint8_t(v));
}
inline int read_varint(const uint8_t* p, const uint8_t* end, uint32_t* out) {
    uint32_t v=0; int sh=0; const uint8_t* s=p;
    while (p<end) { uint8_t x=*p++; v|=uint32_t(x&0x7F)<<sh; sh+=7; if(!(x&0x80))break; }
    *out=v; return int(p-s);
}

// ── Fixed-width LE helpers ────────────────────────────────────────────────────

inline void write_u24(std::vector<uint8_t>& b, uint32_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v>>8)); b.push_back(uint8_t(v>>16));
}
inline uint32_t read_u24(const uint8_t* p) {
    return uint32_t(p[0])|(uint32_t(p[1])<<8)|(uint32_t(p[2])<<16);
}
inline void write_u32(std::vector<uint8_t>& b, uint32_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v>>8));
    b.push_back(uint8_t(v>>16)); b.push_back(uint8_t(v>>24));
}
inline uint32_t read_u32(const uint8_t* p) {
    return uint32_t(p[0])|(uint32_t(p[1])<<8)|(uint32_t(p[2])<<16)|(uint32_t(p[3])<<24);
}
inline void write_f32(std::vector<uint8_t>& b, float v) {
    uint32_t x; std::memcpy(&x,&v,4); write_u32(b,x);
}
inline float read_f32(const uint8_t* p) {
    uint32_t x=read_u32(p); float v; std::memcpy(&v,&x,4); return v;
}
inline void write_f64(std::vector<uint8_t>& b, double v) {
    uint64_t x; std::memcpy(&x,&v,8);
    for(int i=0;i<8;++i) b.push_back(uint8_t(x>>(8*i)));
}
inline double read_f64(const uint8_t* p) {
    uint64_t x=0; for(int i=0;i<8;++i) x|=uint64_t(p[i])<<(8*i);
    double v; std::memcpy(&v,&x,8); return v;
}

// ── Serialization ─────────────────────────────────────────────────────────────

inline void serialize_record(std::vector<uint8_t>& b, const AlignmentRecord& r) {
    write_varint(b, r.hog_idx);
    write_varint(b, r.contig_idx);
    write_u24(b, r.sstart);
    write_u24(b, r.send);
    write_u32(b, r.qstart);
    write_u32(b, r.qend);
    write_u32(b, r.qlen);
    write_f32(b, r.pident);
    write_f64(b, r.evalue);
    b.push_back(uint8_t(r.qframe));
    b.insert(b.end(), r.nt_hash, r.nt_hash + 16);
    b.insert(b.end(), r.aas.begin(), r.aas.end());
    write_varint(b, uint32_t(r.inserts.size()));
    for (auto& ins : r.inserts) {
        write_u24(b, ins.before_hog_pos);
        b.push_back(uint8_t(ins.aas.size()));
        b.insert(b.end(), ins.aas.begin(), ins.aas.end());
    }
}

inline int deserialize_record(const uint8_t* p, const uint8_t* end, AlignmentRecord& r) {
    const uint8_t* s = p;
    int n;
    n = read_varint(p, end, &r.hog_idx);    p += n;
    n = read_varint(p, end, &r.contig_idx); p += n;
    // fixed: 3+3+4+4+4+4+8+1+16 = 47 bytes
    if (p + 47 > end) return 0;
    r.sstart  = read_u24(p); p+=3;
    r.send    = read_u24(p); p+=3;
    r.qstart  = read_u32(p); p+=4;
    r.qend    = read_u32(p); p+=4;
    r.qlen    = read_u32(p); p+=4;
    r.pident  = read_f32(p); p+=4;
    r.evalue  = read_f64(p); p+=8;
    r.qframe  = int8_t(*p++);
    std::memcpy(r.nt_hash, p, 16); p+=16;
    uint32_t span = r.send - r.sstart + 1;
    if (p + span > end) return 0;
    r.aas.assign(p, p+span); p+=span;
    uint32_t n_ins = 0;
    n = read_varint(p, end, &n_ins); p += n;
    r.inserts.resize(n_ins);
    for (auto& ins : r.inserts) {
        if (p + 4 > end) return 0;
        ins.before_hog_pos = read_u24(p); p+=3;
        uint8_t ni = *p++;
        if (p + ni > end) return 0;
        ins.aas.assign(p, p+ni); p+=ni;
    }
    return int(p - s);
}

// ── Variant record (sparse observation, --saav mode) ─────────────────────────
// On-disk layout inside compressed Variants block:
//
//  varint   contig_idx
//  varint   hog_idx
//  uint24   sstart          first HOG position covered by this alignment
//  uint24   send            last HOG position covered
//  float32  pident
//  float64  evalue
//  varint   n_obs           M positions with a valid (non-UNK) AA
//  [n_obs × 4 bytes]:
//    uint24   hog_pos
//    uint8    obs_aa         encoded (0-19)
//
// sstart/send let queries distinguish "not covered" from "covers but matches ref".

struct VariantObs {
    uint32_t hog_pos;
    uint8_t  obs_aa;
};

struct VariantRecord {
    uint32_t contig_idx;
    uint32_t hog_idx;
    uint32_t sstart, send;
    float    pident;
    double   evalue;
    std::vector<VariantObs> obs;
};

inline void serialize_variant(std::vector<uint8_t>& b, const VariantRecord& r) {
    write_varint(b, r.contig_idx);
    write_varint(b, r.hog_idx);
    write_u24(b, r.sstart);
    write_u24(b, r.send);
    write_f32(b, r.pident);
    write_f64(b, r.evalue);
    write_varint(b, uint32_t(r.obs.size()));
    for (auto& o : r.obs) { write_u24(b, o.hog_pos); b.push_back(o.obs_aa); }
}

inline int deserialize_variant(const uint8_t* p, const uint8_t* end, VariantRecord& r) {
    const uint8_t* s = p;
    int n;
    n = read_varint(p, end, &r.contig_idx); p += n;
    n = read_varint(p, end, &r.hog_idx);    p += n;
    if (p + 18 > end) return 0;  // 3+3+4+8
    r.sstart = read_u24(p); p+=3;
    r.send   = read_u24(p); p+=3;
    r.pident = read_f32(p); p+=4;
    r.evalue = read_f64(p); p+=8;
    uint32_t n_obs = 0;
    n = read_varint(p, end, &n_obs); p += n;
    r.obs.resize(n_obs);
    for (auto& o : r.obs) {
        if (p + 4 > end) return 0;
        o.hog_pos = read_u24(p); p+=3;
        o.obs_aa  = *p++;
    }
    return int(p - s);
}


// ── VarNT record (all M-position AA + codon, --varnt mode) ──────────────────
// On-disk layout inside compressed VarNT block:
//
//  varint   contig_idx
//  varint   hog_idx
//  uint24   sstart          HOG coverage extent (tiling gap detection)
//  uint24   send
//  uint32   qstart          nt coords in contig (tiling + codon offset)
//  uint32   qend
//  uint32   qlen
//  int8     qframe
//  float32  pident
//  float64  evalue
//  varint   n_obs
//  per obs (~3 bytes typical):
//    varint  delta_offset    hog_pos delta from previous obs (first = hog_pos - sstart)
//    uint8   obs_aa          encoded (0-19)
//    uint8   packed_codon    2-bit per nt: bits[7:6]=nt0, [5:4]=nt1, [3:2]=nt2  A=0,C=1,G=2,T=3

// 2-bit codon pack/unpack (ACGTN → ACGT, N rounds to T)
inline uint8_t pack_codon(const uint8_t* c) {
    auto e = [](uint8_t n) -> uint8_t {
        switch(n) { case 'A':return 0; case 'C':return 1; case 'G':return 2; default:return 3; }
    };
    return uint8_t((e(c[0])<<6)|(e(c[1])<<4)|(e(c[2])<<2));
}
inline void unpack_codon(uint8_t pk, char* out) {
    constexpr char B[] = {'A','C','G','T'};
    out[0]=B[(pk>>6)&3]; out[1]=B[(pk>>4)&3]; out[2]=B[(pk>>2)&3];
}

struct VarNTObs {
    uint32_t hog_offset;  // in-memory: absolute offset from sstart
    uint8_t  obs_aa;
    uint8_t  packed_codon;
};

struct VarNTRecord {
    uint32_t contig_idx;
    uint32_t hog_idx;
    uint32_t sstart, send;
    uint32_t qstart, qend, qlen;
    int8_t   qframe;
    float    pident;
    double   evalue;
    std::vector<VarNTObs> vars;
};

inline void serialize_varnt(std::vector<uint8_t>& b, const VarNTRecord& r) {
    write_varint(b, r.contig_idx);
    write_varint(b, r.hog_idx);
    write_u24(b, r.sstart);
    write_u24(b, r.send);
    write_u32(b, r.qstart);
    write_u32(b, r.qend);
    write_u32(b, r.qlen);
    b.push_back(uint8_t(r.qframe));
    write_f32(b, r.pident);
    write_f64(b, r.evalue);
    write_varint(b, uint32_t(r.vars.size()));
    uint32_t prev = 0;
    for (auto& v : r.vars) {
        write_varint(b, v.hog_offset - prev);
        prev = v.hog_offset;
        b.push_back(v.obs_aa);
        b.push_back(v.packed_codon);
    }
}

inline int deserialize_varnt(const uint8_t* p, const uint8_t* end, VarNTRecord& r) {
    const uint8_t* s = p;
    int n;
    n = read_varint(p, end, &r.contig_idx); p += n;
    n = read_varint(p, end, &r.hog_idx);    p += n;
    if (p + 31 > end) return 0;  // 3+3+4+4+4+1+4+8
    r.sstart = read_u24(p); p+=3;
    r.send   = read_u24(p); p+=3;
    r.qstart = read_u32(p); p+=4;
    r.qend   = read_u32(p); p+=4;
    r.qlen   = read_u32(p); p+=4;
    r.qframe = int8_t(*p++);
    r.pident = read_f32(p); p+=4;
    r.evalue = read_f64(p); p+=8;
    uint32_t n_vars = 0;
    n = read_varint(p, end, &n_vars); p += n;
    r.vars.resize(n_vars);
    uint32_t cur = 0;
    for (auto& v : r.vars) {
        uint32_t delta = 0;
        n = read_varint(p, end, &delta); p += n;
        cur += delta;
        v.hog_offset  = cur;
        if (p + 2 > end) return 0;
        v.obs_aa       = *p++;
        v.packed_codon = *p++;
    }
    return int(p - s);
}

} // namespace lhi
