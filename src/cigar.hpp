#pragma once
#include "aa.hpp"
#include <string_view>
#include <vector>
#include <stdexcept>
#include <cstdint>

namespace lhi {

// Parse CIGAR + qseq → per-HOG-position AA bytes.
// Returns exactly (send - sstart + 1) bytes (0-indexed: index i = HOG pos sstart+i).
// Gaps in query (D ops) → AA_GAP. Insertions in query (I ops) consume qseq but emit nothing.
inline std::vector<uint8_t> cigar_to_aas(
    std::string_view cigar,
    std::string_view qseq,
    uint32_t sstart,   // 1-based, inclusive
    uint32_t send)     // 1-based, inclusive
{
    const uint32_t span = send - sstart + 1;
    std::vector<uint8_t> aas;
    aas.reserve(span);

    size_t q_off = 0;  // offset into qseq (only M and I consume query)

    const char* p = cigar.data();
    const char* end = p + cigar.size();

    while (p < end && aas.size() < span) {
        // parse length
        uint32_t len = 0;
        while (p < end && *p >= '0' && *p <= '9')
            len = len * 10 + static_cast<uint32_t>(*p++ - '0');
        if (p >= end) break;
        char op = *p++;

        if (op == 'M') {
            for (uint32_t i = 0; i < len && aas.size() < span; ++i, ++q_off) {
                uint8_t aa = (q_off < qseq.size()) ? encode_aa(qseq[q_off]) : AA_UNK;
                aas.push_back(aa);
            }
            // If we consumed more M than needed (shouldn't happen with correct data), skip rest
            if (q_off > qseq.size()) q_off = qseq.size();
        } else if (op == 'I') {
            q_off += len;  // skip query AAs that inserted relative to subject
        } else if (op == 'D') {
            for (uint32_t i = 0; i < len && aas.size() < span; ++i)
                aas.push_back(AA_GAP);
        }
        // S, H, N etc: ignore (shouldn't appear in DIAMOND output)
    }

    // Pad to span if CIGAR was short (malformed data)
    while (aas.size() < span)
        aas.push_back(AA_UNK);

    return aas;
}

} // namespace lhi
