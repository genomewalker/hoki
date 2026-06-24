#pragma once
#include "aa.hpp"
#include <cstdint>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <unistd.h>

namespace lhi {

// RAII wrapper for a POSIX fd: closes on scope exit, move-only.
struct UniqueFd {
    int fd = -1;
    explicit UniqueFd(int f) : fd(f) {}
    ~UniqueFd() { if (fd >= 0) { ::close(fd); fd = -1; } }
    UniqueFd(const UniqueFd&) = delete;
    UniqueFd& operator=(const UniqueFd&) = delete;
    UniqueFd(UniqueFd&& o) noexcept : fd(o.fd) { o.fd = -1; }
    UniqueFd& operator=(UniqueFd&& o) noexcept {
        if (this != &o) { if (fd >= 0) ::close(fd); fd = o.fd; o.fd = -1; }
        return *this;
    }
    operator int() const { return fd; }
};

inline void write_varint(std::vector<uint8_t>& b, uint32_t v) {
    while (v >= 0x80) { b.push_back(uint8_t((v & 0x7F) | 0x80)); v >>= 7; }
    b.push_back(uint8_t(v));
}
inline int read_varint(const uint8_t* p, const uint8_t* end, uint32_t* out) {
    uint32_t v = 0; int sh = 0; const uint8_t* s = p;
    while (p < end && sh < 35) {
        uint8_t x = *p++;
        if (sh == 28) {                          // 5th byte: only low 4 bits fit uint32
            if (x & 0x80) return 0;              // 6th continuation byte ⇒ overflow
            if (x > 0x0F) return 0;              // value > UINT32_MAX
            v |= uint32_t(x) << 28;
            *out = v; return int(p - s);
        }
        v |= uint32_t(x & 0x7F) << sh;
        sh += 7;
        if (!(x & 0x80)) { *out = v; return int(p - s); }
    }
    return 0;  // truncated or overflow
}

inline void write_u32(std::vector<uint8_t>& b, uint32_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v >> 8));
    b.push_back(uint8_t(v >> 16)); b.push_back(uint8_t(v >> 24));
}
inline uint32_t read_u32_le(const uint8_t* p) {
    return uint32_t(p[0]) | (uint32_t(p[1]) << 8) |
           (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}

// Adler-32 over a byte range (self-contained; no extra dependency).
// Used to checksum the raw (non-zstd) index section of a .lhg/.lhgi file.
inline uint32_t adler32(const uint8_t* p, size_t n) {
    uint32_t a = 1, b = 0;
    constexpr uint32_t MOD = 65521;
    for (size_t i = 0; i < n; ++i) { a = (a + p[i]) % MOD; b = (b + a) % MOD; }
    return (b << 16) | a;
}

inline void write_u16(std::vector<uint8_t>& b, uint16_t v) {
    b.push_back(uint8_t(v)); b.push_back(uint8_t(v >> 8));
}
inline uint16_t read_u16(const uint8_t* p) {
    return uint16_t(p[0]) | uint16_t(uint16_t(p[1]) << 8);
}
inline void write_i16(std::vector<uint8_t>& b, int16_t v) {
    write_u16(b, uint16_t(v));
}
inline int16_t read_i16(const uint8_t* p) {
    return int16_t(read_u16(p));
}
inline uint16_t encode_pident(float p)  { return uint16_t(std::lroundf(p * 10.0f)); }
inline float    decode_pident(uint16_t q) { return q / 10.0f; }
inline int16_t encode_evalue(double e) {
    if (e <= 0.0) return int16_t(-32767);
    double x = std::round(std::log10(e) * 100.0);
    if (x < -32767.0) x = -32767.0;
    if (x >     0.0)  x =     0.0;
    return int16_t(x);
}
inline double decode_evalue(int16_t q) { return std::pow(10.0, q / 100.0); }

// codon_idx = packed_codon >> 2 = (nt0<<4)|(nt1<<2)|nt2 (A=0,C=1,G=2,T=3).
// Maps to AA 0–19 (AA_ALPHA order) or AA_UNK (0xFE) for stop codons.
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

// packed_codon bits: [7:6]=nt0 [5:4]=nt1 [3:2]=nt2 [1:0]=unused (A=0,C=1,G=2,T=3)
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

// VarNT block (format v4): SoA columns, sstart-sorted, written by
// serialize_varnt_block below. n_M is derived at read time, not stored.
struct VarNTObs {
    uint32_t hog_offset;    // bit index into presence bitmap = (hog_pos - sstart)
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
                                   std::vector<VarNTRecord>&& recs) {
    const size_t N = recs.size();

    // (contig_idx,sstart)-sorted → contig_col collapses to near-zero deltas; sstart_col
    // stays monotone within each contig run (cross-boundary reset handled in loop below).
    std::stable_sort(recs.begin(), recs.end(), [](const VarNTRecord& a, const VarNTRecord& b) {
        return a.contig_idx < b.contig_idx ||
               (a.contig_idx == b.contig_idx && a.sstart < b.sstart);
    });

    std::vector<uint8_t> contig_col, sstart_col, span_col, qframe_col,
                         pident_col, evalue_col, bmp_buf;
    std::vector<uint8_t> cdn_idx;
    contig_col.reserve(N * 2);
    sstart_col.reserve(N * 2); span_col.reserve(N * 2);
    qframe_col.reserve(N);
    pident_col.reserve(N * 2); evalue_col.reserve(N * 2);
    bmp_buf.reserve(N * 16);
    cdn_idx.reserve(N * 8);

    uint32_t prev_sstart = 0, prev_contig = 0;
    for (const auto& r : recs) {
        uint32_t dc = r.contig_idx - prev_contig;
        write_varint(contig_col, dc);
        if (dc) { prev_sstart = 0; prev_contig = r.contig_idx; }
        write_varint(sstart_col, r.sstart - prev_sstart);
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
            cdn_idx.push_back(uint8_t(v.packed_codon >> 2));
        }
    }

    write_varint(raw, uint32_t(N));
    write_varint(raw, uint32_t(contig_col.size()));
    write_varint(raw, uint32_t(sstart_col.size()));
    write_varint(raw, uint32_t(span_col.size()));
    write_varint(raw, uint32_t(bmp_buf.size()));
    write_varint(raw, uint32_t(cdn_idx.size()));

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

    size_t fixed_sz = size_t(n_recs) * (1 + 2 + 2);  // qframe + pident + evalue
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
    const uint8_t* cdn_ptr    = p;                   p += n_codons;
    size_t cdn_off = 0;

    out.resize(n_recs);
    uint32_t prev_sstart = 0, prev_contig = 0;
    for (uint32_t i = 0; i < n_recs; ++i) {
        auto& r = out[i];
        uint32_t dc = 0;
        n = read_varint(contig_ptr, contig_end, &dc); if (!n) return false;
        contig_ptr += n;
        prev_contig += dc;
        r.contig_idx = prev_contig;
        if (dc) prev_sstart = 0;

        uint32_t dss = 0;
        n = read_varint(sstart_ptr, sstart_end, &dss); if (!n) return false;
        sstart_ptr += n;
        r.sstart = prev_sstart + dss;
        prev_sstart = r.sstart;

        uint32_t span = 0;
        n = read_varint(span_ptr, span_end, &span); if (!n) return false; span_ptr += n;
        if (span == 0) return false;          // span ≥ 1 invariant; 0 ⇒ corrupt
        r.send = r.sstart + span - 1;

        r.qframe = int8_t(*qframe_ptr++);
        r.pident = decode_pident(read_u16(pident_ptr)); pident_ptr += 2;
        r.evalue = decode_evalue(read_i16(evalue_ptr)); evalue_ptr += 2;

        uint32_t bmp_sz = (span + 7) / 8;
        if (bmp_ptr + bmp_sz > bmp_end_p) return false;
        uint32_t n_M = 0;  // observations = popcount of this record's bitmap slice
        for (uint32_t b = 0; b < bmp_sz; ++b) n_M += uint32_t(__builtin_popcount(bmp_ptr[b]));
        r.vars.clear(); r.vars.reserve(n_M);
        for (uint32_t b = 0; b < span; ++b) {
            if (bmp_ptr[b / 8] & (1u << (b & 7))) {
                if (cdn_off >= n_codons) return false;
                uint8_t ci6 = cdn_ptr[cdn_off++];
                VarNTObs obs;
                obs.hog_offset   = b;
                obs.packed_codon = uint8_t(ci6 << 2);  // 6→8 bit; low 2 bits always 0
                r.vars.push_back(obs);
            }
        }
        bmp_ptr += bmp_sz;
    }
    return true;
}

// Column-major variant: bitmap and codon bytes come from separate frame sections (v4).
// p/other_end cover [varnt_header | contig | sstart | span | qframe | pident | evalue].
// bmp_ptr / cdn_ptr point into the frame's bitmap / codon sections for this block.
inline bool deserialize_varnt_block_split(
        const uint8_t* p, const uint8_t* other_end,
        const uint8_t* bmp_ptr, const uint8_t* cdn_ptr,
        std::vector<VarNTRecord>& out) {
    int n;
    uint32_t n_recs = 0, contig_bytes = 0, sstart_bytes = 0, span_bytes = 0,
             bmp_bytes = 0, n_codons = 0;
    n = read_varint(p, other_end, &n_recs);       if (!n) return false; p += n;
    n = read_varint(p, other_end, &contig_bytes); if (!n) return false; p += n;
    n = read_varint(p, other_end, &sstart_bytes); if (!n) return false; p += n;
    n = read_varint(p, other_end, &span_bytes);   if (!n) return false; p += n;
    n = read_varint(p, other_end, &bmp_bytes);    if (!n) return false; p += n;
    n = read_varint(p, other_end, &n_codons);     if (!n) return false; p += n;

    size_t fixed_sz = size_t(n_recs) * (1 + 2 + 2);
    if (p + contig_bytes + sstart_bytes + span_bytes + fixed_sz > other_end)
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
    const uint8_t* bmp_end_p  = bmp_ptr + bmp_bytes;
    size_t cdn_off = 0;

    out.resize(n_recs);
    uint32_t prev_sstart = 0, prev_contig = 0;
    for (uint32_t i = 0; i < n_recs; ++i) {
        auto& r = out[i];
        uint32_t dc = 0;
        n = read_varint(contig_ptr, contig_end, &dc); if (!n) return false;
        contig_ptr += n;
        prev_contig += dc;
        r.contig_idx = prev_contig;
        if (dc) prev_sstart = 0;

        uint32_t dss = 0;
        n = read_varint(sstart_ptr, sstart_end, &dss); if (!n) return false;
        sstart_ptr += n;
        r.sstart = prev_sstart + dss;
        prev_sstart = r.sstart;

        uint32_t span = 0;
        n = read_varint(span_ptr, span_end, &span); if (!n) return false; span_ptr += n;
        if (span == 0) return false;
        r.send = r.sstart + span - 1;

        r.qframe = int8_t(*qframe_ptr++);
        r.pident = decode_pident(read_u16(pident_ptr)); pident_ptr += 2;
        r.evalue = decode_evalue(read_i16(evalue_ptr)); evalue_ptr += 2;

        uint32_t bmp_sz = (span + 7) / 8;
        if (bmp_ptr + bmp_sz > bmp_end_p) return false;
        uint32_t n_M = 0;
        for (uint32_t b = 0; b < bmp_sz; ++b) n_M += uint32_t(__builtin_popcount(bmp_ptr[b]));
        r.vars.clear(); r.vars.reserve(n_M);
        for (uint32_t b = 0; b < span; ++b) {
            if (bmp_ptr[b / 8] & (1u << (b & 7))) {
                if (cdn_off >= n_codons) return false;
                uint8_t ci6 = cdn_ptr[cdn_off++];
                VarNTObs obs;
                obs.hog_offset   = b;
                obs.packed_codon = uint8_t(ci6 << 2);
                r.vars.push_back(obs);
            }
        }
        bmp_ptr += bmp_sz;
    }
    return (cdn_off == n_codons);
}

} // namespace lhi
