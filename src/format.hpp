#pragma once
#include "aa.hpp"
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <array>

// hoki — HOG codon Index
// Wire format helpers shared by container.hpp and cigar.hpp.

namespace lhi {

// ── Varint (LEB128) ──────────────────────────────────────────────────────────

inline void write_varint(std::vector<uint8_t>& b, uint32_t v) {
    while (v >= 0x80) { b.push_back(uint8_t((v & 0x7F) | 0x80)); v >>= 7; }
    b.push_back(uint8_t(v));
}
inline int read_varint(const uint8_t* p, const uint8_t* end, uint32_t* out) {
    uint32_t v = 0; int sh = 0; const uint8_t* s = p;
    while (p < end) { uint8_t x = *p++; v |= uint32_t(x & 0x7F) << sh; sh += 7; if (!(x & 0x80)) break; }
    *out = v; return int(p - s);
}

// ── Fixed-width LE helpers ────────────────────────────────────────────────────

inline void write_u24(std::vector<uint8_t>& b, uint32_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v >> 8)); b.push_back(uint8_t(v >> 16));
}
inline uint32_t read_u24(const uint8_t* p) {
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16);
}
inline void write_u32(std::vector<uint8_t>& b, uint32_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v >> 8));
    b.push_back(uint8_t(v >> 16)); b.push_back(uint8_t(v >> 24));
}
inline uint32_t read_u32(const uint8_t* p) {
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}
inline void write_f32(std::vector<uint8_t>& b, float v) {
    uint32_t x; std::memcpy(&x, &v, 4); write_u32(b, x);
}
inline float read_f32(const uint8_t* p) {
    uint32_t x = read_u32(p); float v; std::memcpy(&v, &x, 4); return v;
}
inline void write_f64(std::vector<uint8_t>& b, double v) {
    uint64_t x; std::memcpy(&x, &v, 8);
    for (int i = 0; i < 8; ++i) b.push_back(uint8_t(x >> (8 * i)));
}
inline double read_f64(const uint8_t* p) {
    uint64_t x = 0; for (int i = 0; i < 8; ++i) x |= uint64_t(p[i]) << (8 * i);
    double v; std::memcpy(&v, &x, 8); return v;
}

// ── Standard genetic code ─────────────────────────────────────────────────────
// codon_idx = packed_codon >> 2 = (nt0<<4)|(nt1<<2)|nt2  (A=0, C=1, G=2, T=3)
// Returns encoded AA (0–19) matching AA_ALPHA = "ACDEFGHIKLMNPQRSTVWY",
// or 0xFE (AA_UNK) for stop codons (TAA=48, TAG=50, TGA=56).

constexpr std::array<uint8_t, 64> CODON_TO_AA = []() constexpr {
    std::array<uint8_t, 64> t{};
    for (auto& x : t) x = uint8_t(0xFE);
    // A**
    t[ 0]=8;  t[ 1]=11; t[ 2]=8;  t[ 3]=11; // AAA=K AAC=N AAG=K AAT=N
    t[ 4]=16; t[ 5]=16; t[ 6]=16; t[ 7]=16; // ACA=T ACC=T ACG=T ACT=T
    t[ 8]=14; t[ 9]=15; t[10]=14; t[11]=15; // AGA=R AGC=S AGG=R AGT=S
    t[12]=7;  t[13]=7;  t[14]=10; t[15]=7;  // ATA=I ATC=I ATG=M ATT=I
    // C**
    t[16]=13; t[17]=6;  t[18]=13; t[19]=6;  // CAA=Q CAC=H CAG=Q CAT=H
    t[20]=12; t[21]=12; t[22]=12; t[23]=12; // CCA=P CCC=P CCG=P CCT=P
    t[24]=14; t[25]=14; t[26]=14; t[27]=14; // CGA=R CGC=R CGG=R CGT=R
    t[28]=9;  t[29]=9;  t[30]=9;  t[31]=9;  // CTA=L CTC=L CTG=L CTT=L
    // G**
    t[32]=3;  t[33]=2;  t[34]=3;  t[35]=2;  // GAA=E GAC=D GAG=E GAT=D
    t[36]=0;  t[37]=0;  t[38]=0;  t[39]=0;  // GCA=A GCC=A GCG=A GCT=A
    t[40]=5;  t[41]=5;  t[42]=5;  t[43]=5;  // GGA=G GGC=G GGG=G GGT=G
    t[44]=17; t[45]=17; t[46]=17; t[47]=17; // GTA=V GTC=V GTG=V GTT=V
    // T** (TAA=48, TAG=50, TGA=56 stay 0xFE = stop)
    t[49]=19; t[51]=19;                      // TAC=Y TAT=Y
    t[52]=15; t[53]=15; t[54]=15; t[55]=15; // TCA=S TCC=S TCG=S TCT=S
    t[57]=1;  t[58]=18; t[59]=1;             // TGC=C TGG=W TGT=C
    t[60]=9;  t[61]=4;  t[62]=9;  t[63]=4;  // TTA=L TTC=F TTG=L TTT=F
    return t;
}();

inline uint8_t codon_to_aa(uint8_t packed_codon) {
    return CODON_TO_AA[packed_codon >> 2];
}

// ── Codon pack / unpack ───────────────────────────────────────────────────────
// bits[7:6]=nt0, [5:4]=nt1, [3:2]=nt2, [1:0]=unused  (A=0,C=1,G=2,T=3)

inline uint8_t pack_codon(const uint8_t* c) {
    auto e = [](uint8_t n) -> uint8_t {
        switch (n) { case 'A': return 0; case 'C': return 1; case 'G': return 2; default: return 3; }
    };
    return uint8_t((e(c[0]) << 6) | (e(c[1]) << 4) | (e(c[2]) << 2));
}
inline void unpack_codon(uint8_t pk, char* out) {
    constexpr char B[] = {'A','C','G','T'};
    out[0] = B[(pk >> 6) & 3]; out[1] = B[(pk >> 4) & 3]; out[2] = B[(pk >> 2) & 3];
}

// ── VarNT record (codon-resolved alignment) ───────────────────────────────────
//
// On-disk inside a shard block payload (see container.hpp):
//
//   VarNT columnar block:
//     varint  n_records
//     varint  hdr_section_bytes
//     varint  bmp_section_bytes
//     varint  cdn_section_bytes
//     [hdr section] per record:
//       varint  contig_idx     — index into block-local contig list
//       varint  hog_idx        — always 0 for single-HOG shards
//       uint24  sstart         — 1-based HOG MSA position, coverage start
//       uint24  send           — coverage end
//       uint32  qstart/qend/qlen  — nt coords in query contig
//       int8    qframe         — ±1..±3 (sign=strand, magnitude=reading frame)
//       float32 pident
//       float64 evalue
//       varint  n_M            — number of codon observations in this record
//     [bmp section] per record:
//       ceil((send-sstart+1)/8) bytes, LSB-first
//       bit i = 1 → HOG position (sstart+i) has an observed codon
//     [cdn section]:
//       one packed_codon byte per set bit, ascending position order
//
// obs_aa is NOT stored; derived at query time: codon_to_aa(packed_codon).

struct VarNTObs {
    uint32_t hog_offset;    // bit index into presence bitmap = (hog_pos - sstart)
    uint8_t  obs_aa;        // derived in memory, not on disk
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

inline void serialize_varnt_block(std::vector<uint8_t>& raw,
                                   const std::vector<VarNTRecord>& recs) {
    std::vector<uint8_t> hdr_buf, bmp_buf, cdn_buf;
    for (auto& r : recs) {
        write_varint(hdr_buf, r.contig_idx);
        write_varint(hdr_buf, r.hog_idx);
        write_u24(hdr_buf, r.sstart); write_u24(hdr_buf, r.send);
        write_u32(hdr_buf, r.qstart); write_u32(hdr_buf, r.qend); write_u32(hdr_buf, r.qlen);
        hdr_buf.push_back(uint8_t(r.qframe));
        write_f32(hdr_buf, r.pident); write_f64(hdr_buf, r.evalue);
        write_varint(hdr_buf, uint32_t(r.vars.size()));

        uint32_t span   = r.send - r.sstart + 1;
        uint32_t bmp_sz = (span + 7) / 8;
        size_t   bmp_off = bmp_buf.size();
        bmp_buf.resize(bmp_off + bmp_sz, 0);
        for (auto& v : r.vars)
            bmp_buf[bmp_off + v.hog_offset / 8] |= uint8_t(1u << (v.hog_offset & 7));
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

inline bool deserialize_varnt_block(const uint8_t* p, const uint8_t* end,
                                     std::vector<VarNTRecord>& out) {
    int n;
    uint32_t n_recs = 0, hdr_bytes = 0, bmp_bytes = 0, cdn_bytes = 0;
    n = read_varint(p, end, &n_recs);    if (!n) return false; p += n;
    n = read_varint(p, end, &hdr_bytes); if (!n) return false; p += n;
    n = read_varint(p, end, &bmp_bytes); if (!n) return false; p += n;
    n = read_varint(p, end, &cdn_bytes); if (!n) return false; p += n;

    const uint8_t* hdr_ptr = p,            *hdr_end = p + hdr_bytes;
    const uint8_t* bmp_ptr = hdr_end,      *bmp_end = hdr_end + bmp_bytes;
    const uint8_t* cdn_ptr = bmp_end,      *cdn_end = bmp_end + cdn_bytes;
    if (cdn_end > end) return false;

    out.resize(n_recs);
    std::vector<uint32_t> n_Ms(n_recs);

    for (uint32_t i = 0; i < n_recs; ++i) {
        auto& r = out[i];
        n = read_varint(hdr_ptr, hdr_end, &r.contig_idx); if (!n) return false; hdr_ptr += n;
        n = read_varint(hdr_ptr, hdr_end, &r.hog_idx);    if (!n) return false; hdr_ptr += n;
        if (hdr_ptr + 31 > hdr_end) return false;
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
        if (bmp_ptr + bmp_sz > bmp_end || cdn_ptr + n_M > cdn_end) return false;
        r.vars.clear(); r.vars.reserve(n_M);
        for (uint32_t b = 0; b < span; ++b) {
            if (bmp_ptr[b / 8] & (1u << (b & 7))) {
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
