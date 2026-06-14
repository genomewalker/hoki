#pragma once
#include "aa.hpp"
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

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
// uint16 LE
inline void write_u16(std::vector<uint8_t>& b, uint16_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v >> 8));
}
inline uint16_t read_u16(const uint8_t* p) {
    return uint16_t(p[0]) | uint16_t(uint16_t(p[1]) << 8);
}
// int16 LE
inline void write_i16(std::vector<uint8_t>& b, int16_t v) {
    write_u16(b, uint16_t(v));
}
inline int16_t read_i16(const uint8_t* p) {
    return int16_t(read_u16(p));
}
// Encode pident (0..100 float) → uint16 (×10, range 0–1000)
inline uint16_t encode_pident(float p)  { return uint16_t(std::lroundf(p * 10.0f)); }
inline float    decode_pident(uint16_t q) { return q / 10.0f; }
// Encode evalue → int16 (floor(log10(e)×100), clamped to [−32767, 0]).
// Reconstruct via pow(10, q/100.0).  Precision: ~3 sig-figs of the exponent.
inline int16_t encode_evalue(double e) {
    if (e <= 0.0) return int16_t(-32767);
    double x = std::round(std::log10(e) * 100.0);
    if (x < -32767.0) x = -32767.0;
    if (x >     0.0)  x =     0.0;
    return int16_t(x);
}
inline double decode_evalue(int16_t q) { return std::pow(10.0, q / 100.0); }

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

// ── 6-bit codon packing ───────────────────────────────────────────────────────
// codon_idx = packed_codon >> 2 (6 bits, 64 values).  4 indices fit in 3 bytes.
// Bit layout for a group of 4 values a,b,c,d packed into bytes B0,B1,B2:
//   B0 = (a<<2)|(b>>4)   B1 = (b<<4)|(c>>2)   B2 = (c<<6)|d
// n values need ceil(n*6/8) bytes.

inline size_t six_bit_packed_size(size_t n) { return (n * 6 + 7) / 8; }

inline void pack_6bit(uint8_t* dst, const uint8_t* src, size_t n) {
    size_t o = 0;
    for (size_t i = 0; i < n; i += 4) {
        size_t rem = n - i;
        uint8_t a = src[i];
        uint8_t b = rem > 1 ? src[i+1] : 0;
        uint8_t c = rem > 2 ? src[i+2] : 0;
        uint8_t d = rem > 3 ? src[i+3] : 0;
        dst[o++] = uint8_t((a << 2) | (b >> 4));
        if (rem > 1) dst[o++] = uint8_t((b << 4) | (c >> 2));
        if (rem > 2) dst[o++] = uint8_t((c << 6) |  d);
    }
}

inline void unpack_6bit(const uint8_t* src, uint8_t* dst, size_t n) {
    size_t in = 0;
    for (size_t i = 0; i < n; i += 4) {
        size_t rem = n - i;
        uint8_t b0 = src[in++];
        uint8_t b1 = rem > 1 ? src[in++] : 0;
        uint8_t b2 = rem > 2 ? src[in++] : 0;
        dst[i]   = (b0 >> 2) & 0x3F;
        if (rem > 1) dst[i+1] = uint8_t(((b0 & 3) << 4) | (b1 >> 4));
        if (rem > 2) dst[i+2] = uint8_t(((b1 & 0xF) << 2) | (b2 >> 6));
        if (rem > 3) dst[i+3] = b2 & 0x3F;
    }
}

// ── VarNT record (codon-resolved alignment) ───────────────────────────────────
//
// On-disk inside a shard block payload (see container.hpp):
//
// Format v4 — SoA HDR, sstart-sorted, delta sstart + span, raw-byte CDN:
//   varint  n_records
//   varint  contig_col_bytes  — size of contig_idx varint column
//   varint  sstart_col_bytes  — size of sstart delta-varint column
//   varint  span_col_bytes    — size of span varint column
//   varint  bmp_section_bytes
//   varint  n_codons          — total codon count across all records (= cdn bytes)
//   [contig_idx column: n_records varints, contig_col_bytes total]
//   [sstart column:  n_records delta-varints (first = absolute), sstart_col_bytes]
//   [span column:    n_records varints (span = send - sstart + 1), span_col_bytes]
//   [qframe column:  n_records × 1 byte (int8)]
//   [pident column:  n_records × 2 bytes (uint16 ×10)]
//   [evalue column:  n_records × 2 bytes (int16, floor(log10(e)×100))]
//   [bmp section:    per-record bitmap, ceil(span/8) bytes each]
//   [cdn section:    raw codon_idx bytes (one per codon), n_codons total]
//
// Records are sorted by sstart ascending before serialization, so the sstart
// delta column is monotonically non-decreasing (typically 1 varint byte each).
// n_M (observations per record) is NOT stored; derived at query time as the
// popcount of the record's bitmap slice.
// obs_aa is NOT stored; derived at query time: codon_to_aa(packed_codon).

struct VarNTObs {
    uint32_t hog_offset;    // bit index into presence bitmap = (hog_pos - sstart)
    uint8_t  obs_aa;        // derived in memory, not on disk
    uint8_t  packed_codon;
};

struct VarNTRecord {
    uint32_t contig_idx;
    uint32_t sstart, send;
    int8_t   qframe;   // ±1..±3: sign=strand, magnitude=reading frame
    float    pident;
    double   evalue;
    std::vector<VarNTObs> vars;
};

inline void serialize_varnt_block(std::vector<uint8_t>& raw,
                                   std::vector<VarNTRecord> recs) {
    const size_t N = recs.size();

    // Sort by sstart ascending → monotone sstart column compresses far better
    // and enables delta encoding.  recs is taken by value, caller untouched.
    std::stable_sort(recs.begin(), recs.end(),
        [](const VarNTRecord& a, const VarNTRecord& b) { return a.sstart < b.sstart; });

    // SoA: build each column separately so zstd sees homogeneous data per field.
    std::vector<uint8_t> contig_col, sstart_col, span_col, qframe_col,
                         pident_col, evalue_col, bmp_buf;
    std::vector<uint8_t> cdn_idx;  // raw codon indices (one byte each, value 0–63)
    contig_col.reserve(N * 2);
    sstart_col.reserve(N * 2); span_col.reserve(N * 2);
    qframe_col.reserve(N);
    pident_col.reserve(N * 2); evalue_col.reserve(N * 2);
    bmp_buf.reserve(N * 16);
    cdn_idx.reserve(N * 8);

    uint32_t prev_sstart = 0;
    for (const auto& r : recs) {
        write_varint(contig_col, r.contig_idx);
        write_varint(sstart_col, r.sstart - prev_sstart);  // delta (first = absolute)
        prev_sstart = r.sstart;
        uint32_t span = r.send - r.sstart + 1;
        write_varint(span_col, span);
        qframe_col.push_back(uint8_t(r.qframe));
        write_u16(pident_col, encode_pident(r.pident));
        write_i16(evalue_col, encode_evalue(r.evalue));

        uint32_t bmp_sz = (span + 7) / 8;
        size_t   bmp_off = bmp_buf.size();
        bmp_buf.resize(bmp_off + bmp_sz, 0);
        for (const auto& v : r.vars) {
            bmp_buf[bmp_off + v.hog_offset / 8] |= uint8_t(1u << (v.hog_offset & 7));
            cdn_idx.push_back(uint8_t(v.packed_codon >> 2));  // raw 6-bit value in a byte
        }
    }

    // Header: sizes of variable-width columns; fixed columns derivable from N.
    write_varint(raw, uint32_t(N));
    write_varint(raw, uint32_t(contig_col.size()));
    write_varint(raw, uint32_t(sstart_col.size()));
    write_varint(raw, uint32_t(span_col.size()));
    write_varint(raw, uint32_t(bmp_buf.size()));
    write_varint(raw, uint32_t(cdn_idx.size()));  // n_codons (= cdn section bytes)

    // Columns in order: contig, sstart, span (all varint), fixed (1+2+2 ×N), bmp, cdn.
    raw.insert(raw.end(), contig_col.begin(), contig_col.end());
    raw.insert(raw.end(), sstart_col.begin(), sstart_col.end());
    raw.insert(raw.end(), span_col.begin(),   span_col.end());
    raw.insert(raw.end(), qframe_col.begin(), qframe_col.end());
    raw.insert(raw.end(), pident_col.begin(), pident_col.end());
    raw.insert(raw.end(), evalue_col.begin(), evalue_col.end());
    raw.insert(raw.end(), bmp_buf.begin(),    bmp_buf.end());
    raw.insert(raw.end(), cdn_idx.begin(),    cdn_idx.end());
}

inline bool deserialize_varnt_block(const uint8_t* p, const uint8_t* end,
                                     std::vector<VarNTRecord>& out) {
    int n;
    uint32_t n_recs = 0, contig_bytes = 0, sstart_bytes = 0, span_bytes = 0,
             bmp_bytes = 0, n_codons = 0;
    n = read_varint(p, end, &n_recs);       if (!n) return false; p += n;
    n = read_varint(p, end, &contig_bytes); if (!n) return false; p += n;
    n = read_varint(p, end, &sstart_bytes); if (!n) return false; p += n;
    n = read_varint(p, end, &span_bytes);   if (!n) return false; p += n;
    n = read_varint(p, end, &bmp_bytes);    if (!n) return false; p += n;
    n = read_varint(p, end, &n_codons);     if (!n) return false; p += n;

    // Fixed-column sizes are derivable from n_recs (qframe+pident+evalue).
    size_t fixed_sz = size_t(n_recs) * (1 + 2 + 2);
    if (p + contig_bytes + sstart_bytes + span_bytes + fixed_sz + bmp_bytes + n_codons > end)
        return false;

    const uint8_t* contig_ptr = p;                   p += contig_bytes;
    const uint8_t* contig_end = p;
    const uint8_t* sstart_ptr = p;                   p += sstart_bytes;
    const uint8_t* sstart_end = p;
    const uint8_t* span_ptr   = p;                   p += span_bytes;
    const uint8_t* span_end   = p;
    const uint8_t* qframe_ptr = p;                   p += 1 * n_recs;
    const uint8_t* pident_ptr = p;                   p += 2 * n_recs;
    const uint8_t* evalue_ptr = p;                   p += 2 * n_recs;
    const uint8_t* bmp_ptr    = p;                   p += bmp_bytes;
    const uint8_t* bmp_end_p  = p;
    const uint8_t* cdn_ptr    = p;                   p += n_codons;  // raw byte per codon
    size_t cdn_off = 0;

    out.resize(n_recs);
    uint32_t prev_sstart = 0;
    for (uint32_t i = 0; i < n_recs; ++i) {
        auto& r = out[i];
        uint32_t ci = 0;
        n = read_varint(contig_ptr, contig_end, &ci); if (!n) return false;
        contig_ptr += n;
        r.contig_idx = ci;

        uint32_t dss = 0;
        n = read_varint(sstart_ptr, sstart_end, &dss); if (!n) return false;
        sstart_ptr += n;
        r.sstart = prev_sstart + dss;  // cumulative sum of deltas
        prev_sstart = r.sstart;

        uint32_t span = 0;
        n = read_varint(span_ptr, span_end, &span); if (!n) return false; span_ptr += n;
        if (span == 0) return false;          // corrupt: span ≥ 1
        r.send = r.sstart + span - 1;

        r.qframe = int8_t(*qframe_ptr++);
        r.pident = decode_pident(read_u16(pident_ptr)); pident_ptr += 2;
        r.evalue = decode_evalue(read_i16(evalue_ptr)); evalue_ptr += 2;

        uint32_t bmp_sz = (span + 7) / 8;
        if (bmp_ptr + bmp_sz > bmp_end_p) return false;  // bitmap overrun guard
        // n_M derived as popcount of this record's bitmap slice.
        uint32_t n_M = 0;
        for (uint32_t b = 0; b < bmp_sz; ++b) n_M += uint32_t(__builtin_popcount(bmp_ptr[b]));
        r.vars.clear(); r.vars.reserve(n_M);
        for (uint32_t b = 0; b < span; ++b) {
            if (bmp_ptr[b / 8] & (1u << (b & 7))) {
                if (cdn_off >= n_codons) return false;
                uint8_t ci6 = cdn_ptr[cdn_off++];
                VarNTObs obs;
                obs.hog_offset   = b;
                obs.packed_codon = uint8_t(ci6 << 2);  // 6→8 bit (lower 2 bits always 0)
                obs.obs_aa       = codon_to_aa(obs.packed_codon);
                r.vars.push_back(obs);
            }
        }
        bmp_ptr += bmp_sz;
    }
    return true;
}

} // namespace lhi
