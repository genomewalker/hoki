#pragma once
#include "format.hpp"
#include "cigar.hpp"
#include "aa.hpp"
#include "container.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <limits>

namespace lhi {

struct ConvertOptions {
    int    zstd_level = 3;
    size_t block_recs = 50000;  // records per shard flush batch
    float  min_pident = 0.0f;
    double max_evalue = 1.0;
    bool   verbose    = false;
};

// diamond blastx outfmt 6 columns (0-based):
//   qseqid qstart qend qlen qstrand sseqid sstart send slen pident evalue cigar qseq_translated full_qseq
namespace col {
    constexpr int qseqid=0, qstart=1, qend=2, qlen=3, qstrand=4;
    constexpr int sseqid=5, sstart=6, send=7,  slen=8;
    constexpr int pident=9, evalue=10, cigar=11, qseq_aa=12, full_qseq=13;
}

inline std::string extract_hog(std::string_view sv) {
    auto p = sv.rfind('|');
    return std::string(p != std::string_view::npos ? sv.substr(p + 1) : sv);
}

inline int8_t make_qframe(std::string_view qstrand, uint32_t qstart) {
    int frame = int((qstart > 0 ? qstart - 1 : 0) % 3) + 1;
    return int8_t(qstrand.empty() || qstrand[0] != '-' ? frame : -frame);
}

inline double parse_double(std::string_view sv) {
    char buf[64];
    size_t n = std::min(sv.size(), sizeof(buf) - 1);
    std::memcpy(buf, sv.data(), n); buf[n] = '\0';
    char* end; return std::strtod(buf, &end);
}

// Reverse-complement 3 nt bytes in place.
inline void revcomp_codon(uint8_t c[3]) {
    auto rc = [](uint8_t b) -> uint8_t {
        switch (b) { case 'A':return 'T'; case 'T':return 'A'; case 'C':return 'G'; case 'G':return 'C'; default:return 'N'; }
    };
    uint8_t t[3] = {rc(c[2]), rc(c[1]), rc(c[0])};
    c[0]=t[0]; c[1]=t[1]; c[2]=t[2];
}

// Convert diamond TSV → HOG-sharded container directory.
// Only VarNT (codon-resolved, 100% pident skipped) records are written.
// Multiple concurrent writers are safe: each HOG shard uses flock(LOCK_EX).
inline void convert(const std::string& tsv_path, const std::string& container_dir,
                    ConvertOptions opts) {

    std::filesystem::create_directories(
        std::filesystem::path(container_dir) / "shards");

    std::unordered_map<std::string, uint32_t> hog_dict, contig_dict;
    std::vector<std::string> hog_strings, contig_strings;

    auto intern = [](std::unordered_map<std::string,uint32_t>& d,
                     std::vector<std::string>& v, const std::string& s) -> uint32_t {
        auto [it, ins] = d.emplace(s, uint32_t(v.size()));
        if (ins) v.push_back(s);
        return it->second;
    };

    // Accumulate VarNT records per HOG index.
    std::unordered_map<uint32_t, std::vector<VarNTRecord>> batches;

    auto flush_hog = [&](uint32_t hog_idx) {
        auto it = batches.find(hog_idx);
        if (it == batches.end() || it->second.empty()) return;
        auto shard = std::filesystem::path(container_dir) / "shards"
                   / hog_to_filename(hog_strings[hog_idx]);
        flush_hog_shard(shard, contig_strings, it->second, opts.zstd_level);
        it->second.clear();
    };

    std::istream* in = &std::cin;
    std::ifstream fin;
    if (tsv_path != "-") {
        fin.open(tsv_path);
        if (!fin) throw std::runtime_error("cannot open: " + tsv_path);
        in = &fin;
    }

    std::string line;
    uint64_t lineno = 0, n_written = 0, n_skipped = 0;

    while (std::getline(*in, line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;

        std::array<std::string_view, 14> f{};
        size_t fi = 0;
        const char* p = line.data(), *ep = p + line.size(), *fp = p;
        while (p <= ep && fi < 14) {
            if (p == ep || *p == '\t') { f[fi++] = {fp, size_t(p-fp)}; fp = p+1; }
            ++p;
        }
        if (fi < 13) { if (opts.verbose) std::cerr << "skip L" << lineno << ": " << fi << " fields\n"; ++n_skipped; continue; }
        bool has_nt = (fi == 14 && !f[col::full_qseq].empty());

        float  pident = float(parse_double(f[col::pident]));
        double ev     = parse_double(f[col::evalue]);

        if (pident < opts.min_pident) { ++n_skipped; continue; }
        if (ev     > opts.max_evalue) { ++n_skipped; continue; }
        if (pident == 100.0f)         { ++n_skipped; continue; }  // no variant info at exact identity

        uint32_t sstart, send, qstart, qend, qlen;
        try {
            sstart = uint32_t(std::stoul(std::string(f[col::sstart])));
            send   = uint32_t(std::stoul(std::string(f[col::send])));
            qstart = uint32_t(std::stoul(std::string(f[col::qstart])));
            qend   = uint32_t(std::stoul(std::string(f[col::qend])));
            qlen   = uint32_t(std::stoul(std::string(f[col::qlen])));
        } catch (...) { if (opts.verbose) std::cerr << "parse error L" << lineno << "\n"; ++n_skipped; continue; }
        if (sstart > send) { ++n_skipped; continue; }

        std::string hog_id     = extract_hog(f[col::sseqid]);
        uint32_t    hog_idx    = intern(hog_dict, hog_strings, hog_id);
        uint32_t    contig_idx = intern(contig_dict, contig_strings, std::string(f[col::qseqid]));
        int8_t      qframe     = make_qframe(f[col::qstrand], qstart);

        auto ar = cigar_parse(f[col::cigar], f[col::qseq_aa], sstart, send);

        VarNTRecord vr;
        vr.contig_idx = contig_idx;
        vr.hog_idx    = hog_idx;
        vr.sstart     = sstart; vr.send = send;
        vr.qstart     = qstart; vr.qend = qend; vr.qlen = qlen;
        vr.qframe     = qframe;
        vr.pident     = pident; vr.evalue = ev;

        std::string_view full_nt = has_nt ? f[col::full_qseq] : std::string_view{};
        uint32_t span = send - sstart + 1;

        for (uint32_t i = 0; i < span; ++i) {
            uint8_t obs_aa = ar.aas[i];
            if (obs_aa == AA_GAP || obs_aa == AA_UNK) continue;
            uint32_t q_off = ar.qseq_offsets[i];
            if (q_off == UINT32_MAX || full_nt.empty()) continue;

            uint8_t c0, c1, c2;
            if (qframe > 0) {
                size_t cs = size_t(qstart-1) + size_t(q_off)*3;
                if (cs+2 >= full_nt.size()) continue;
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
            } else {
                if (size_t(q_off)*3+2 > size_t(qend-1)) continue;
                size_t cs = size_t(qend-1) - size_t(q_off)*3 - 2;
                if (cs+2 >= full_nt.size()) continue;
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
                uint8_t tmp[3]={c0,c1,c2}; revcomp_codon(tmp); c0=tmp[0]; c1=tmp[1]; c2=tmp[2];
            }

            uint8_t raw3[3]={c0,c1,c2};
            vr.vars.push_back({i, obs_aa, pack_codon(raw3)});
        }

        if (!vr.vars.empty()) {
            batches[hog_idx].push_back(std::move(vr));
            if (batches[hog_idx].size() >= opts.block_recs) flush_hog(hog_idx);
        }
        ++n_written;
    }

    // Final flush: all remaining HOGs.
    for (auto& [hog_idx, _] : batches) flush_hog(hog_idx);

    std::cerr << "wrote " << n_written << " records, " << n_skipped
              << " skipped → " << container_dir << "\n";
}

} // namespace lhi
