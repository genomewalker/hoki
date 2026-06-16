#pragma once
#include "convert.hpp"
#include "container.hpp"
#include "partition.hpp"

namespace lhi {

// TSV → partition dir directly. No intermediate .lhb files.
// Supports -a auto (acc from qseqid prefix) or -a ACC (single accession).
// Writes out_dir/t0.lhp + out_dir/partition.idx + out_dir/acc.registry.
inline void ingest(const std::string& tsv_path, const std::string& out_dir,
                   ConvertOptions opts) {
    std::filesystem::create_directories(out_dir);

    struct AccState {
        uint32_t acc_idx;
        SvDict contig_dict;
        std::vector<std::string> contig_strings;
        std::unordered_map<uint32_t, std::vector<VarNTRecord>> batches;
    };

    bool auto_acc = (opts.acc_id == "auto");
    std::unordered_map<std::string, AccState> acc_states;
    SvDict hog_dict;
    std::vector<std::string> hog_strings;

    auto intern = [](SvDict& d, std::vector<std::string>& v, std::string_view s) -> uint32_t {
        auto it = d.find(s);
        if (it != d.end()) return it->second;
        uint32_t idx = uint32_t(v.size());
        v.emplace_back(s);
        d.emplace(v.back(), idx);
        return idx;
    };

    TsvReader reader(tsv_path);
    std::string line;
    uint64_t lineno = 0, n_written = 0, n_skipped = 0, n_obs_dropped = 0;

    while (reader.getline(line)) {
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
                !fc(f[col::qstart], qstart) || !fc(f[col::qend],  qend) ||
                !fc(f[col::qlen],   qlen)) {
                if (opts.verbose) std::cerr << "parse error L" << lineno << "\n";
                ++n_skipped; continue;
            }
        }
        if (sstart > send) { ++n_skipped; continue; }

        std::string_view qseqid = f[col::qseqid];
        std::string acc;
        if (auto_acc) {
            auto us = qseqid.find('_');
            acc = std::string(us != std::string_view::npos ? qseqid.substr(0, us) : qseqid);
        } else {
            acc = opts.acc_id;
        }

        AccState& st = acc_states[acc];
        if (st.contig_strings.empty() && st.contig_dict.empty())
            st.acc_idx = uint32_t(acc_states.size() - 1); // assigned below in sorted pass

        uint32_t hog_idx    = intern(hog_dict, hog_strings, extract_hog(f[col::sseqid]));
        uint32_t contig_idx = intern(st.contig_dict, st.contig_strings, qseqid);
        int8_t   qframe     = make_qframe(f[col::qstrand], qstart, qend, qlen);
        auto     ar         = cigar_parse(f[col::cigar], f[col::qseq_aa], sstart, send);

        VarNTRecord vr;
        vr.contig_idx = contig_idx;
        vr.sstart = sstart; vr.send = send;
        vr.qframe = qframe;
        vr.pident = pident; vr.evalue = ev;

        std::string_view full_nt = has_nt ? f[col::full_qseq] : std::string_view{};
        uint32_t span = send - sstart + 1;

        for (uint32_t i = 0; i < span; ++i) {
            uint8_t obs_aa = ar.aas[i];
            if (obs_aa == AA_GAP || obs_aa == AA_UNK) { ++n_obs_dropped; continue; }
            uint32_t q_off = ar.qseq_offsets[i];
            if (q_off == UINT32_MAX || full_nt.empty()) { ++n_obs_dropped; continue; }

            uint8_t c0, c1, c2;
            if (qframe > 0) {
                size_t cs = size_t(qstart-1) + size_t(q_off)*3;
                if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
            } else {
                if (size_t(q_off)*3 + 3 > size_t(qstart)) { ++n_obs_dropped; continue; }
                size_t cs = size_t(qstart-1) - size_t(q_off)*3 - 2;
                if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
                uint8_t tmp[3]={c0,c1,c2}; revcomp_codon(tmp); c0=tmp[0]; c1=tmp[1]; c2=tmp[2];
            }
            auto is_acgt = [](uint8_t b) { return b=='A'||b=='C'||b=='G'||b=='T'; };
            if (!is_acgt(c0) || !is_acgt(c1) || !is_acgt(c2)) { ++n_obs_dropped; continue; }

            uint8_t raw3[3]={c0,c1,c2};
            uint8_t packed = pack_codon(raw3);
            if (codon_to_aa(packed) != obs_aa) {
                if (opts.verbose) std::cerr << "codon/AA mismatch L" << lineno << "\n";
                ++n_obs_dropped; continue;
            }
            vr.vars.push_back({i, packed});
        }
        if (!vr.vars.empty())
            st.batches[hog_idx].push_back(std::move(vr));
        ++n_written;
    }

    // Assign sorted acc_idx (alphabetical, for deterministic merge across jobs).
    std::vector<std::string> accessions;
    accessions.reserve(acc_states.size());
    for (auto& [acc, _] : acc_states) accessions.push_back(acc);
    std::sort(accessions.begin(), accessions.end());
    for (uint32_t i = 0; i < uint32_t(accessions.size()); ++i)
        acc_states[accessions[i]].acc_idx = i;

    write_acc_registry(accessions, out_dir + "/acc.registry");

    // Write all HOG blocks for all accessions to t0.lhp in HOG-ID order.
    PartitionWriter pw(out_dir + "/t0.lhp");

    // Collect and sort HOG indices
    std::vector<uint32_t> sorted_hogs;
    for (auto& [acc, st] : acc_states)
        for (auto& [hog_idx, _] : st.batches)
            sorted_hogs.push_back(hog_idx);
    std::sort(sorted_hogs.begin(), sorted_hogs.end(),
        [&](uint32_t a, uint32_t b) { return hog_strings[a] < hog_strings[b]; });
    sorted_hogs.erase(std::unique(sorted_hogs.begin(), sorted_hogs.end()), sorted_hogs.end());

    for (uint32_t hog_idx : sorted_hogs) {
        const std::string& hog_id = hog_strings[hog_idx];
        for (auto& [acc, st] : acc_states) {
            auto it = st.batches.find(hog_idx);
            if (it == st.batches.end() || it->second.empty()) continue;
            auto sp = compress_shard_payload(st.contig_strings, it->second, opts.zstd_level);
            PartitionEntry ent;
            ent.acc_idx       = st.acc_idx;
            ent.compressed_sz = uint32_t(sp.cbuf.size());
            ent.raw_sz        = sp.raw_sz;
            ent.n_records     = sp.n_records;
            pw.append(hog_id, ent, sp.cbuf.data());
        }
    }
    pw.flush_all();

    // Build and write partition index (single thread → thread_idx always 0)
    std::map<std::string, std::vector<PartitionIndexExtent>> global_idx;
    for (const auto& [hog, extents] : pw.index())
        for (const auto& [off, len] : extents)
            global_idx[hog].push_back({0u, off, len});
    write_partition_index(global_idx, 1u, out_dir + "/partition.idx");

    std::cerr << "ingest: " << n_written << " records, " << n_skipped << " skipped, "
              << n_obs_dropped << " obs dropped, "
              << accessions.size() << " accessions, "
              << global_idx.size() << " HOGs → " << out_dir << "\n";
}

} // namespace lhi
