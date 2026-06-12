#pragma once
#include "aa.hpp"
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <array>

// LHI v4 — Logan HOG Index, fully lossless alignment format
// v4 vs v3: VarNT blocks use columnar encoding: within each compressed block,
// record headers, presence bitmaps, and codon bytes are stored as three contiguous
// streams (better zstd entropy on homogeneous data).  obs_aa is dropped from disk
// and derived at query time from packed_codon via CODON_TO_AA[].

namespace lhi {

constexpr uint8_t  MAGIC[8]       = {'L','H','I','v','4','\0','\0','\0'};
constexpr uint16_t FORMAT_VERSION = 4;

// Contract flags (FileHeader.flags)
constexpr uint16_t FLAG_VALUE_LOSSLESS = 0x0001;
constexpr uint16_t FLAG_BYTE_LOSSLESS  = 0x0002;
constexpr uint16_t FLAG_HAS_NT_BLOCKS  = 0x0004;

enum class BlockType : uint8_t {
    HOGDict    = 0,
    ContigDict = 1,
    Alignments = 2,
    Pileup     = 3,
    NT         = 4,
    Variants   = 5,  // sparse: only M-position (hog_pos, obs_aa) pairs per record
    VarNT      = 6,  // variant+codon+tiling: columnar block (v4)
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
    float    min_evalue_log;
    float    max_evalue_log;
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
    uint8_t  nt_hash[16];
    std::vector<uint8_t>        aas;
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
//  uint24   sstart
//  uint24   send
//  float32  pident
//  float64  evalue
//  varint   n_obs
//  [n_obs × 4 bytes]:
//    uint24   hog_pos
//    uint8    obs_aa

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
    if (p + 18 > end) return 0;
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

// ── Standard genetic code ─────────────────────────────────────────────────────
// codon_idx = packed_codon >> 2 = (nt0<<4)|(nt1<<2)|nt2  (A=0, C=1, G=2, T=3)
// Returns encoded AA (0-19) matching AA_ALPHA="ACDEFGHIKLMNPQRSTVWY",
// or 0xFE (AA_UNK) for stop codons (TAA=48, TAG=50, TGA=56).

constexpr std::array<uint8_t, 64> CODON_TO_AA = []() constexpr {
    std::array<uint8_t, 64> t{};
    for (auto& x : t) x = uint8_t(0xFE);  // AA_UNK default
    // A** (nt0=A, idx 0-15)
    t[ 0]=8;  t[ 1]=11; t[ 2]=8;  t[ 3]=11; // AAA=K AAC=N AAG=K AAT=N
    t[ 4]=16; t[ 5]=16; t[ 6]=16; t[ 7]=16; // ACA=T ACC=T ACG=T ACT=T
    t[ 8]=14; t[ 9]=15; t[10]=14; t[11]=15; // AGA=R AGC=S AGG=R AGT=S
    t[12]=7;  t[13]=7;  t[14]=10; t[15]=7;  // ATA=I ATC=I ATG=M ATT=I
    // C** (nt0=C, idx 16-31)
    t[16]=13; t[17]=6;  t[18]=13; t[19]=6;  // CAA=Q CAC=H CAG=Q CAT=H
    t[20]=12; t[21]=12; t[22]=12; t[23]=12; // CCA=P CCC=P CCG=P CCT=P
    t[24]=14; t[25]=14; t[26]=14; t[27]=14; // CGA=R CGC=R CGG=R CGT=R
    t[28]=9;  t[29]=9;  t[30]=9;  t[31]=9;  // CTA=L CTC=L CTG=L CTT=L
    // G** (nt0=G, idx 32-47)
    t[32]=3;  t[33]=2;  t[34]=3;  t[35]=2;  // GAA=E GAC=D GAG=E GAT=D
    t[36]=0;  t[37]=0;  t[38]=0;  t[39]=0;  // GCA=A GCC=A GCG=A GCT=A
    t[40]=5;  t[41]=5;  t[42]=5;  t[43]=5;  // GGA=G GGC=G GGG=G GGT=G
    t[44]=17; t[45]=17; t[46]=17; t[47]=17; // GTA=V GTC=V GTG=V GTT=V
    // T** (nt0=T, idx 48-63): TAA=48, TAG=50, TGA=56 are stops (left 0xFE)
    t[49]=19; t[51]=19;                      // TAC=Y TAT=Y
    t[52]=15; t[53]=15; t[54]=15; t[55]=15; // TCA=S TCC=S TCG=S TCT=S
    t[57]=1;  t[58]=18; t[59]=1;             // TGC=C TGG=W TGT=C
    t[60]=9;  t[61]=4;  t[62]=9;  t[63]=4;  // TTA=L TTC=F TTG=L TTT=F
    return t;
}();

inline uint8_t codon_to_aa(uint8_t packed_codon) {
    return CODON_TO_AA[packed_codon >> 2];
}

// ── VarNT record (all M-position codon + tiling, --varnt mode) ───────────────
// Columnar block layout (v4):
//
//   varint   n_records
//   varint   hdr_section_bytes
//   varint   bmp_section_bytes
//   varint   cdn_section_bytes
//   --- hdr section ---
//   per record:
//     varint  contig_idx
//     varint  hog_idx
//     uint24  sstart          1-based HOG coverage start
//     uint24  send            1-based HOG coverage end
//     uint32  qstart, qend, qlen    nt coords in query contig
//     int8    qframe
//     float32 pident
//     float64 evalue
//     varint  n_M             M-position count (= popcount of bitmap)
//   --- bmp section ---
//   per record: ceil((send-sstart+1)/8) bytes, LSB-first
//     bit i = 1 → HOG position (sstart+i) is an M position
//   --- cdn section ---
//   per record: n_M bytes
//     packed_codon per M position in ascending offset order
//     (A=0,C=1,G=2,T=3; bits[7:6]=nt0, [5:4]=nt1, [3:2]=nt2; bits[1:0] unused)
//
// obs_aa is NOT stored on disk; derived at query time as codon_to_aa(packed_codon).

// 2-bit codon pack/unpack
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
    uint32_t hog_offset;   // offset from sstart; derived from bitmap bit index
    uint8_t  obs_aa;       // derived in memory: codon_to_aa(packed_codon) — not on disk
    uint8_t  packed_codon; // bits[7:2]: nt0/nt1/nt2 packed 2-bit; bits[1:0] unused
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

// Serialize all records in a block into three columnar streams then concatenate.
inline void serialize_varnt_block(std::vector<uint8_t>& raw,
                                   const std::vector<VarNTRecord>& recs) {
    std::vector<uint8_t> hdr_buf, bmp_buf, cdn_buf;

    for (auto& r : recs) {
        write_varint(hdr_buf, r.contig_idx);
        write_varint(hdr_buf, r.hog_idx);
        write_u24(hdr_buf, r.sstart);
        write_u24(hdr_buf, r.send);
        write_u32(hdr_buf, r.qstart);
        write_u32(hdr_buf, r.qend);
        write_u32(hdr_buf, r.qlen);
        hdr_buf.push_back(uint8_t(r.qframe));
        write_f32(hdr_buf, r.pident);
        write_f64(hdr_buf, r.evalue);
        write_varint(hdr_buf, uint32_t(r.vars.size()));

        uint32_t span   = r.send - r.sstart + 1;
        uint32_t bmp_sz = (span + 7) / 8;
        size_t   bmp_off = bmp_buf.size();
        bmp_buf.resize(bmp_off + bmp_sz, 0);
        for (auto& v : r.vars)
            bmp_buf[bmp_off + v.hog_offset/8] |= uint8_t(1u << (v.hog_offset & 7));

        for (auto& v : r.vars)
            cdn_buf.push_back(v.packed_codon);
    }

    write_varint(raw, uint32_t(recs.size()));
    write_varint(raw, uint32_t(hdr_buf.size()));
    write_varint(raw, uint32_t(bmp_buf.size()));
    write_varint(raw, uint32_t(cdn_buf.size()));
    raw.insert(raw.end(), hdr_buf.begin(), hdr_buf.end());
    raw.insert(raw.end(), bmp_buf.begin(), bmp_buf.end());
    raw.insert(raw.end(), cdn_buf.begin(), cdn_buf.end());
}

// Deserialize a full VarNT block into out[].  Returns false on truncation.
inline bool deserialize_varnt_block(const uint8_t* p, const uint8_t* end,
                                     std::vector<VarNTRecord>& out) {
    int n;
    uint32_t n_recs=0, hdr_bytes=0, bmp_bytes=0, cdn_bytes=0;
    n = read_varint(p, end, &n_recs);   if (!n) return false; p += n;
    n = read_varint(p, end, &hdr_bytes); if (!n) return false; p += n;
    n = read_varint(p, end, &bmp_bytes); if (!n) return false; p += n;
    n = read_varint(p, end, &cdn_bytes); if (!n) return false; p += n;

    const uint8_t* hdr_ptr = p;
    const uint8_t* hdr_end = p + hdr_bytes;
    const uint8_t* bmp_ptr = hdr_end;
    const uint8_t* bmp_end = bmp_ptr + bmp_bytes;
    const uint8_t* cdn_ptr = bmp_end;
    const uint8_t* cdn_end = cdn_ptr + cdn_bytes;
    if (cdn_end > end) return false;

    out.resize(n_recs);
    std::vector<uint32_t> n_Ms(n_recs);

    for (uint32_t i = 0; i < n_recs; ++i) {
        auto& r = out[i];
        n = read_varint(hdr_ptr, hdr_end, &r.contig_idx); if (!n) return false; hdr_ptr += n;
        n = read_varint(hdr_ptr, hdr_end, &r.hog_idx);    if (!n) return false; hdr_ptr += n;
        if (hdr_ptr + 31 > hdr_end) return false;  // 3+3+4+4+4+1+4+8
        r.sstart = read_u24(hdr_ptr); hdr_ptr += 3;
        r.send   = read_u24(hdr_ptr); hdr_ptr += 3;
        r.qstart = read_u32(hdr_ptr); hdr_ptr += 4;
        r.qend   = read_u32(hdr_ptr); hdr_ptr += 4;
        r.qlen   = read_u32(hdr_ptr); hdr_ptr += 4;
        r.qframe = int8_t(*hdr_ptr++);
        r.pident = read_f32(hdr_ptr); hdr_ptr += 4;
        r.evalue = read_f64(hdr_ptr); hdr_ptr += 8;
        n = read_varint(hdr_ptr, hdr_end, &n_Ms[i]); if (!n) return false; hdr_ptr += n;
    }

    for (uint32_t i = 0; i < n_recs; ++i) {
        auto& r = out[i];
        uint32_t span   = r.send - r.sstart + 1;
        uint32_t bmp_sz = (span + 7) / 8;
        uint32_t n_M    = n_Ms[i];
        if (bmp_ptr + bmp_sz > bmp_end) return false;
        if (cdn_ptr + n_M   > cdn_end)  return false;

        r.vars.clear();
        r.vars.reserve(n_M);
        for (uint32_t b = 0; b < span; ++b) {
            if (bmp_ptr[b/8] & (1u << (b & 7))) {
                VarNTObs obs;
                obs.hog_offset   = b;
                obs.packed_codon = *cdn_ptr++;
                obs.obs_aa       = codon_to_aa(obs.packed_codon);
                r.vars.push_back(obs);
            }
        }
        bmp_ptr += bmp_sz;
    }
    return true;
}

} // namespace lhi
