#pragma once
#include "format.hpp"
#include "aa.hpp"
#include <string_view>
#include <vector>

namespace lhi {

struct AlignedResult {
    std::vector<uint8_t>  aas;           // M/D positions; length = send-sstart+1
    std::vector<uint32_t> qseq_offsets;  // parallel to aas[]; UINT32_MAX for D positions
};

// Parse CIGAR + qseq_translated → AA array + per-position query offsets.
inline AlignedResult cigar_parse(
    std::string_view cigar,
    std::string_view qseq,
    uint32_t sstart, uint32_t send)
{
    const uint32_t span = send - sstart + 1;
    AlignedResult r;
    r.aas.reserve(span);
    r.qseq_offsets.reserve(span);

    size_t   q_off   = 0;
    uint32_t hog_pos = sstart;

    const char* p = cigar.data(), *end = p + cigar.size();
    while (p < end && r.aas.size() < span) {
        uint32_t len = 0;
        while (p < end && *p >= '0' && *p <= '9') len = len * 10 + uint32_t(*p++ - '0');
        if (p >= end || !len) break;
        char op = *p++;

        if (op == 'M') {
            for (uint32_t i = 0; i < len && r.aas.size() < span; ++i, ++hog_pos) {
                r.qseq_offsets.push_back(uint32_t(q_off));
                r.aas.push_back(q_off < qseq.size() ? encode_aa(qseq[q_off]) : AA_UNK);
                ++q_off;
            }
        } else if (op == 'I') {
            q_off += len;  // skip inserted residues (not in subject space)
        } else if (op == 'D') {
            for (uint32_t i = 0; i < len && r.aas.size() < span; ++i, ++hog_pos) {
                r.qseq_offsets.push_back(UINT32_MAX);
                r.aas.push_back(AA_GAP);
            }
        }
    }
    while (r.aas.size() < span) {
        r.qseq_offsets.push_back(UINT32_MAX);
        r.aas.push_back(AA_UNK);
    }
    return r;
}

// In-place variant: clears and fills a pre-allocated AlignedResult, avoiding
// per-record vector allocation (save ~2 allocs per cigar_parse call).
inline void cigar_parse_inplace(
    std::string_view cigar, std::string_view qseq,
    uint32_t sstart, uint32_t send, AlignedResult& r)
{
    const uint32_t span = send - sstart + 1;
    r.aas.clear();
    r.qseq_offsets.clear();
    r.aas.reserve(span);
    r.qseq_offsets.reserve(span);
    size_t q_off = 0; uint32_t hog_pos = sstart;
    const char* p = cigar.data(), *end = p + cigar.size();
    while (p < end && r.aas.size() < span) {
        uint32_t len = 0;
        while (p < end && *p >= '0' && *p <= '9') len = len * 10 + uint32_t(*p++ - '0');
        if (p >= end || !len) break;
        char op = *p++;
        if (op == 'M') {
            for (uint32_t i = 0; i < len && r.aas.size() < span; ++i, ++hog_pos) {
                r.qseq_offsets.push_back(uint32_t(q_off));
                r.aas.push_back(q_off < qseq.size() ? encode_aa(qseq[q_off]) : AA_UNK);
                ++q_off;
            }
        } else if (op == 'I') {
            q_off += len;
        } else if (op == 'D') {
            for (uint32_t i = 0; i < len && r.aas.size() < span; ++i, ++hog_pos) {
                r.qseq_offsets.push_back(UINT32_MAX);
                r.aas.push_back(AA_GAP);
            }
        }
    }
    while (r.aas.size() < span) {
        r.qseq_offsets.push_back(UINT32_MAX);
        r.aas.push_back(AA_UNK);
    }
}

} // namespace lhi
