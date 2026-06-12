#pragma once
#include "format.hpp"  // InsertionRecord lives here
#include "aa.hpp"
#include <string_view>
#include <vector>

namespace lhi {

struct AlignedResult {
    std::vector<uint8_t>         aas;     // M/D positions; length = send-sstart+1
    std::vector<InsertionRecord> inserts; // I ops, sparse (fully lossless)
};

// Parse CIGAR + qseq_translated → AA array + insertion events
inline AlignedResult cigar_parse(
    std::string_view cigar,
    std::string_view qseq,
    uint32_t sstart, uint32_t send)
{
    const uint32_t span = send - sstart + 1;
    AlignedResult r;
    r.aas.reserve(span);

    size_t   q_off   = 0;
    uint32_t hog_pos = sstart;

    const char* p = cigar.data(), *end = p + cigar.size();
    while (p < end && r.aas.size() < span) {
        uint32_t len = 0;
        while (p < end && *p >= '0' && *p <= '9') len = len*10 + uint32_t(*p++ - '0');
        if (p >= end || !len) break;
        char op = *p++;

        if (op == 'M') {
            for (uint32_t i = 0; i < len && r.aas.size() < span; ++i, ++q_off, ++hog_pos)
                r.aas.push_back(q_off < qseq.size() ? encode_aa(qseq[q_off]) : AA_UNK);
        } else if (op == 'I') {
            InsertionRecord ev;
            ev.before_hog_pos = hog_pos;
            ev.aas.reserve(len);
            for (uint32_t i = 0; i < len; ++i, ++q_off)
                ev.aas.push_back(q_off < qseq.size() ? encode_aa(qseq[q_off]) : AA_UNK);
            r.inserts.push_back(std::move(ev));
        } else if (op == 'D') {
            for (uint32_t i = 0; i < len && r.aas.size() < span; ++i, ++hog_pos)
                r.aas.push_back(AA_GAP);
        }
    }
    while (r.aas.size() < span) r.aas.push_back(AA_UNK);
    return r;
}

} // namespace lhi
