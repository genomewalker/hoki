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

    // Create writer before the read loop so mid-stream flushes can write to it.
    PartitionWriter pw(out_dir + "/t0.lhp");

    const size_t flush_threshold = [&]() -> size_t {
        if (opts.flush_bytes != 0) return opts.flush_bytes;
        // cgroup v2: walk hierarchy from /proc/self/cgroup, find constrained memory.max
        auto read_mem_limit = []() -> unsigned long long {
            FILE* fg = fopen("/proc/self/cgroup", "r");
            if (!fg) return 0;
            char cgpath[512] = {};
            char line[512];
            while (fgets(line, sizeof(line), fg)) {
                // format: 0::/user.slice/...
                if (line[0] == '0' && line[1] == ':' && line[2] == ':') {
                    size_t n = strlen(line);
                    while (n > 0 && (line[n-1]=='\n'||line[n-1]=='\r')) line[--n]='\0';
                    snprintf(cgpath, sizeof(cgpath), "%s", line+3);
                    break;
                }
            }
            fclose(fg);
            if (!cgpath[0]) return 0;
            // Walk up from /sys/fs/cgroup/<cgpath> to root
            char base[640];
            snprintf(base, sizeof(base), "/sys/fs/cgroup%s", cgpath);
            while (strlen(base) > strlen("/sys/fs/cgroup")) {
                char memmax[700];
                snprintf(memmax, sizeof(memmax), "%s/memory.max", base);
                FILE* fm = fopen(memmax, "r");
                if (fm) {
                    char val[64] = {};
                    fgets(val, sizeof(val), fm);
                    fclose(fm);
                    if (val[0] != 'm') { // not "max"
                        unsigned long long v = 0;
                        sscanf(val, "%llu", &v);
                        if (v > 0) return v;
                    }
                }
                char* slash = strrchr(base, '/');
                if (!slash || slash == base) break;
                *slash = '\0';
            }
            return 0;
        };
        // cgroup v1: walk up from the task's memory cgroup until we find a real limit
        auto read_v1 = []() -> unsigned long long {
            FILE* fg = fopen("/proc/self/cgroup", "r");
            if (!fg) return 0;
            char cgrel[512] = {};
            char line[512];
            while (fgets(line, sizeof(line), fg)) {
                char* p1 = strchr(line, ':');
                if (!p1) continue;
                char* p2 = strchr(p1+1, ':');
                if (!p2) continue;
                *p2 = '\0';
                if (strstr(p1+1, "memory")) {
                    size_t n = strlen(p2+1);
                    while (n > 0 && ((p2+1)[n-1]=='\n'||(p2+1)[n-1]=='\r')) (p2+1)[--n]='\0';
                    snprintf(cgrel, sizeof(cgrel), "%s", p2+1);
                    break;
                }
            }
            fclose(fg);
            if (!cgrel[0]) return 0;
            // Walk up: /slurm/uid_X/job_Y/step_Z/task_0 → job_Y → ... → root
            char base[640];
            snprintf(base, sizeof(base), "/sys/fs/cgroup/memory%s", cgrel);
            while (strlen(base) > strlen("/sys/fs/cgroup/memory")) {
                char fpath[700];
                snprintf(fpath, sizeof(fpath), "%s/memory.limit_in_bytes", base);
                FILE* f = fopen(fpath, "r");
                if (f) {
                    unsigned long long v = 0; fscanf(f, "%llu", &v); fclose(f);
                    if (v > 0 && v < (1ull<<62)) return v;
                }
                char* slash = strrchr(base, '/');
                if (!slash || slash == base) break;
                *slash = '\0';
            }
            return 0;
        };
        unsigned long long limit = read_mem_limit();
        if (!limit) limit = read_v1();
        if (limit) {
            size_t flush = size_t(limit * 7 / 10);
            std::cerr << "ingest: cgroup mem " << (limit>>20) << " MiB → flush at "
                      << (flush>>20) << " MiB\n";
            return flush;
        }
        return 2048ull * 1024 * 1024; // 2 GiB fallback
    }();
    size_t mem_bytes = 0;
    uint64_t n_flushes = 0;

    auto flush_batches = [&]() {
        for (auto& [acc, st] : acc_states) {
            for (auto& [hog_idx, batch] : st.batches) {
                if (batch.empty()) continue;
                auto sp = compress_shard_payload(st.contig_strings, batch, opts.zstd_level);
                PartitionEntry ent{st.acc_idx, uint32_t(sp.cbuf.size()), sp.raw_sz, sp.n_records};
                pw.append(hog_strings[hog_idx], ent, sp.cbuf.data());
            }
            st.batches.clear();
        }
        mem_bytes = 0;
        ++n_flushes;
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

        bool inserted = (acc_states.count(acc) == 0);
        AccState& st = acc_states[acc];
        if (inserted) st.acc_idx = uint32_t(acc_states.size() - 1); // discovery order

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
        if (!vr.vars.empty()) {
            mem_bytes += sizeof(VarNTRecord) + vr.vars.size() * sizeof(vr.vars[0]);
            st.batches[hog_idx].push_back(std::move(vr));
        }
        ++n_written;

        if (mem_bytes > flush_threshold)
            flush_batches();
    }

    if (mem_bytes > 0) flush_batches();
    pw.flush_all();

    std::vector<std::string> accessions(acc_states.size());
    for (auto& [acc, st] : acc_states) accessions[st.acc_idx] = acc;
    write_acc_registry(accessions, out_dir + "/acc.registry");

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
