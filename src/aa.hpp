#pragma once
#include <cstdint>
#include <string_view>
#include <array>

namespace lhi {

constexpr uint8_t AA_GAP = 0xFF;  // D in CIGAR (query gap at HOG position)
constexpr uint8_t AA_UNK = 0xFE;  // X or unparseable

// ACDEFGHIKLMNPQRSTVWY — alphabetical standard 20
constexpr std::string_view AA_ALPHA = "ACDEFGHIKLMNPQRSTVWY";

constexpr std::array<uint8_t, 128> build_aa_table() {
    std::array<uint8_t, 128> t{};
    t.fill(AA_UNK);
    for (uint8_t i = 0; i < 20; ++i)
        t[static_cast<uint8_t>(AA_ALPHA[i])] = i;
    t['B'] = 11; t['Z'] = 13; t['J'] = 9;
    t['U'] = 1;  t['O'] = 8;
    t['X'] = AA_UNK; t['*'] = AA_UNK; t['-'] = AA_GAP;
    return t;
}
constexpr auto AA_TABLE = build_aa_table();

inline uint8_t encode_aa(char c) {
    auto u = static_cast<unsigned char>(c);
    return (u < 128) ? AA_TABLE[u] : AA_UNK;
}

} // namespace lhi
