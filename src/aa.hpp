#pragma once
#include <cstdint>
#include <string_view>
#include <array>
#include <cmath>

// 20-AA index + gap + unknown sentinels
namespace lhi {

constexpr uint8_t AA_GAP  = 0xFF;  // D in CIGAR (query gap)
constexpr uint8_t AA_UNK  = 0xFE;  // X or unparseable

// ACDEFGHIKLMNPQRSTVWY — alphabetical, standard 20
constexpr std::string_view AA_ALPHA = "ACDEFGHIKLMNPQRSTVWY";

constexpr std::array<uint8_t, 128> build_aa_table() {
    std::array<uint8_t, 128> t{};
    t.fill(AA_UNK);
    for (uint8_t i = 0; i < 20; ++i)
        t[static_cast<uint8_t>(AA_ALPHA[i])] = i;
    // common synonyms
    t['B'] = 11; // N or D → N
    t['Z'] = 13; // Q or E → Q
    t['J'] = 9;  // I or L → L
    t['U'] = 1;  // selenocysteine → C
    t['O'] = 8;  // pyrrolysine → K
    t['X'] = AA_UNK;
    t['*'] = AA_UNK;
    t['-'] = AA_GAP;
    return t;
}

constexpr auto AA_TABLE = build_aa_table();

inline uint8_t encode_aa(char c) {
    auto uc = static_cast<unsigned char>(c);
    return (uc < 128) ? AA_TABLE[uc] : AA_UNK;
}

inline char decode_aa(uint8_t v) {
    if (v < 20) return AA_ALPHA[v];
    if (v == AA_GAP) return '-';
    return 'X';
}

// evalue → int8: round(log10(evalue)), clamped [-128, 0]
inline int8_t quantize_evalue(double ev) {
    if (ev <= 0.0) return -128;
    double l = std::log10(ev);
    if (l > 0.0) l = 0.0;
    if (l < -128.0) l = -128.0;
    return static_cast<int8_t>(std::round(l));
}

// pident → uint8: round(pident * 2.55), recovers at ~0.4% resolution
inline uint8_t quantize_pident(float p) {
    if (p < 0.0f) return 0;
    if (p > 100.0f) return 255;
    return static_cast<uint8_t>(std::round(p * 2.55f));
}
inline float recover_pident(uint8_t q) { return q / 2.55f; }

} // namespace lhi
