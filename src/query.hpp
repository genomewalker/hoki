#pragma once
#include "format.hpp"
#include "aa.hpp"
#include "container.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <cstdio>

namespace lhi {

struct QueryOptions {
    std::string hog_id;
    uint32_t    pos        = 0;
    uint8_t     variant_aa = AA_UNK;
    float       min_pident = 0.0f;
    double      max_evalue = 1.0;
};

// Output TSV: contig_id \t hog_pos \t obs_aa \t codon \t pident \t evalue
inline void query_varnt(const std::string& container_dir, const QueryOptions& opts) {
    auto shard = std::filesystem::path(container_dir) / "shards"
               / hog_to_filename(opts.hog_id);

    std::printf("contig_id\thog_pos\tobs_aa\tcodon\tpident\tevalue\n");

    bool found = read_shard_file(shard.string(),
        [&](const std::vector<std::string>& contigs, const std::vector<VarNTRecord>& recs) {
            for (auto& vr : recs) {
                if (vr.pident < opts.min_pident) continue;
                if (vr.evalue > opts.max_evalue)  continue;
                for (auto& o : vr.vars) {
                    uint32_t hog_pos = vr.sstart + o.hog_offset;
                    if (opts.pos && hog_pos != opts.pos) continue;
                    if (opts.variant_aa != AA_UNK && o.obs_aa != opts.variant_aa) continue;
                    char aac = (o.obs_aa < 20) ? AA_ALPHA[o.obs_aa] : 'X';
                    char cdn[3]; unpack_codon(o.packed_codon, cdn);
                    std::printf("%s\t%u\t%c\t%c%c%c\t%.2f\t%g\n",
                        contigs[vr.contig_idx].c_str(),
                        hog_pos, aac, cdn[0], cdn[1], cdn[2],
                        vr.pident, vr.evalue);
                }
            }
        });

    if (!found)
        std::cerr << "HOG " << opts.hog_id << " not in container " << container_dir << "\n";
}

// Find tiling gaps across accessions for a HOG.
// Output TSV: accession \t gap_start \t gap_end \t gap_aa_size \t n_exon_records
inline void query_tile(const std::string& container_dir, const std::string& hog_id,
                        float min_pident, double max_evalue, uint32_t min_acc) {
    auto shard = std::filesystem::path(container_dir) / "shards"
               / hog_to_filename(hog_id);

    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>> acc_spans;

    bool found = read_shard_file(shard.string(),
        [&](const std::vector<std::string>& contigs, const std::vector<VarNTRecord>& recs) {
            for (auto& vr : recs) {
                if (vr.pident < min_pident) continue;
                if (vr.evalue > max_evalue)  continue;
                const std::string& contig = contigs[vr.contig_idx];
                auto dot = contig.rfind('.');
                std::string acc = (dot != std::string::npos) ? contig.substr(0, dot) : contig;
                acc_spans[acc].emplace_back(vr.sstart, vr.send);
            }
        });

    if (!found) { std::cerr << "HOG " << hog_id << " not in container " << container_dir << "\n"; return; }

    std::unordered_map<uint64_t, uint32_t> gap_acc_count;

    std::printf("accession\tgap_start\tgap_end\tgap_aa_size\tn_exon_records\n");
    for (auto& [acc, spans] : acc_spans) {
        std::sort(spans.begin(), spans.end());
        for (size_t i = 1; i < spans.size(); ++i) {
            if (spans[i].first > spans[i-1].second + 1) {
                uint32_t gs = spans[i-1].second, ge = spans[i].first;
                uint64_t key = (uint64_t(gs) << 32) | uint64_t(ge);
                gap_acc_count[key]++;
                if (gap_acc_count[key] >= min_acc)
                    std::printf("%s\t%u\t%u\t%u\t%zu\n", acc.c_str(), gs, ge, ge-gs-1, spans.size());
            }
        }
    }

    std::printf("\n# Gap summary (position -> n_accessions)\n");
    std::printf("gap_start\tgap_end\tgap_aa_size\tn_accessions\n");
    std::vector<std::pair<uint64_t,uint32_t>> sorted_gaps(gap_acc_count.begin(), gap_acc_count.end());
    std::sort(sorted_gaps.begin(), sorted_gaps.end(), [](auto& a, auto& b){ return a.first < b.first; });
    for (auto& [key, cnt] : sorted_gaps) {
        if (cnt < min_acc) continue;
        uint32_t gs = uint32_t(key >> 32), ge = uint32_t(key & 0xFFFFFFFF);
        std::printf("%u\t%u\t%u\t%u\n", gs, ge, ge-gs-1, cnt);
    }
}

} // namespace lhi
