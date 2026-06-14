#pragma once
#include "format.hpp"
#include "cigar.hpp"
#include "aa.hpp"
#include "batch.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <charconv>

namespace lhi {

struct ConvertOptions {
    std::string acc_id;          // SRA/ENA accession (required for .lhb output)
    int    zstd_level = 9;
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

inline int8_t make_qframe(std::string_view qstrand, uint32_t qstart,
                          uint32_t qend, uint32_t qlen) {
    bool minus = (!qstrand.empty() && qstrand[0] == '-');
    uint32_t base = minus ? (qlen > qend ? qlen - qend : 0)
                          : (qstart > 0  ? qstart - 1  : 0);
    int frame = int(base % 3) + 1;
    return int8_t(minus ? -frame : frame);
}

inline void revcomp_codon(uint8_t c[3]) {
    auto rc = [](uint8_t b) -> uint8_t {
        switch (b) { case 'A':return 'T'; case 'T':return 'A'; case 'C':return 'G'; case 'G':return 'C'; default:return 'N'; }
    };
    uint8_t t[3] = {rc(c[2]), rc(c[1]), rc(c[0])};
    c[0]=t[0]; c[1]=t[1]; c[2]=t[2];
}

// Convert a diamond blastx TSV (outfmt 6 + qseq_translated + full_qseq columns)
// to a .lhb batch file.  All HOGs go into a single output file; no per-HOG shards.
inline void convert(const std::string& tsv_path, const std::string& lhb_path,
                    ConvertOptions opts) {
    if (opts.acc_id.empty())
        throw std::runtime_error("acc_id is required for .lhb output (use -a ACC)");

    BatchWriter batch(lhb_path, opts.acc_id);

    std::unordered_map<std::string, uint32_t> hog_dict, contig_dict;
    std::vector<std::string> hog_strings, contig_strings;

    auto intern = [](std::unordered_map<std::string,uint32_t>& d,
                     std::vector<std::string>& v, const std::string& s) -> uint32_t {
        auto [it, ins] = d.emplace(s, uint32_t(v.size()));
        if (ins) v.push_back(s);
        return it->second;
    };

    std::unordered_map<uint32_t, std::vector<VarNTRecord>> batches;

    auto flush_hog = [&](uint32_t hog_idx) {
        auto it = batches.find(hog_idx);
        if (it == batches.end() || it->second.empty()) return;
        batch.write_block(hog_strings[hog_idx], contig_strings, it->second, opts.zstd_level);
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

        float  pident = 0.0f;
        double ev     = 0.0;
        std::from_chars(f[col::pident].data(), f[col::pident].data() + f[col::pident].size(), pident);
        std::from_chars(f[col::evalue].data(),  f[col::evalue].data()  + f[col::evalue].size(),  ev);

        if (pident < opts.min_pident) { ++n_skipped; continue; }
        if (ev     > opts.max_evalue) { ++n_skipped; continue; }
        if (pident == 100.0f)         { ++n_skipped; continue; }

        uint32_t sstart = 0, send = 0, qstart = 0, qend = 0, qlen = 0;
        {
            auto fc = [](std::string_view sv, uint32_t& out) -> bool {
                auto [p, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), out);
                return ec == std::errc{};
            };
            if (!fc(f[col::sstart], sstart) || !fc(f[col::send], send) ||
                !fc(f[col::qstart], qstart) || !fc(f[col::qend], qend) ||
                !fc(f[col::qlen],   qlen)) {
                if (opts.verbose) std::cerr << "parse error L" << lineno << "\n";
                ++n_skipped; continue;
            }
        }
        if (sstart > send) { ++n_skipped; continue; }

        std::string hog_id     = extract_hog(f[col::sseqid]);
        uint32_t    hog_idx    = intern(hog_dict, hog_strings, hog_id);
        uint32_t    contig_idx = intern(contig_dict, contig_strings, std::string(f[col::qseqid]));
        int8_t      qframe     = make_qframe(f[col::qstrand], qstart, qend, qlen);

        auto ar = cigar_parse(f[col::cigar], f[col::qseq_aa], sstart, send);

        VarNTRecord vr;
        vr.contig_idx = contig_idx;
        vr.sstart     = sstart; vr.send = send;
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
                // Minus strand: diamond reports qstart > qend; codon i sits at
                // [(qstart-1) - q_off*3 - 2 .. (qstart-1) - q_off*3], revcomp'd.
                if (size_t(q_off)*3 + 3 > size_t(qstart)) continue;
                size_t cs = size_t(qstart-1) - size_t(q_off)*3 - 2;
                if (cs+2 >= full_nt.size()) continue;
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
                uint8_t tmp[3]={c0,c1,c2}; revcomp_codon(tmp); c0=tmp[0]; c1=tmp[1]; c2=tmp[2];
            }

            // Reject ambiguous NT (N, lowercase, non-ACGT) — pack_codon would coerce to T.
            auto is_acgt = [](uint8_t b) {
                return b=='A'||b=='C'||b=='G'||b=='T';
            };
            if (!is_acgt(c0) || !is_acgt(c1) || !is_acgt(c2)) { ++n_skipped; continue; }

            uint8_t raw3[3]={c0,c1,c2};
            uint8_t packed = pack_codon(raw3);

            // Round-trip: codon→AA must agree with diamond's translated AA.
            // Disagreement indicates a strand/frame extraction error; discard the observation.
            if (codon_to_aa(packed) != obs_aa) {
                if (opts.verbose) std::cerr << "codon/AA mismatch L" << lineno
                    << " (diamond=" << char(AA_ALPHA[obs_aa])
                    << " nt=" << char(c0) << char(c1) << char(c2) << ")\n";
                ++n_skipped; continue;
            }

            vr.vars.push_back({i, obs_aa, packed});
        }

        if (!vr.vars.empty())
            batches[hog_idx].push_back(std::move(vr));
        ++n_written;
    }

    // Flush in HOG-ID order so .lhb blocks are lexicographically sorted, which
    // lets hoki merge do a streaming k-way merge.
    std::vector<uint32_t> sorted_idxs;
    sorted_idxs.reserve(batches.size());
    for (auto& [hog_idx, recs] : batches)
        if (!recs.empty()) sorted_idxs.push_back(hog_idx);
    std::sort(sorted_idxs.begin(), sorted_idxs.end(),
        [&](uint32_t a, uint32_t b) { return hog_strings[a] < hog_strings[b]; });
    for (auto hog_idx : sorted_idxs) flush_hog(hog_idx);
    batch.finalize();

    std::cerr << "convert: " << n_written << " records, " << n_skipped
              << " skipped → " << lhb_path << "\n";
}

} // namespace lhi
