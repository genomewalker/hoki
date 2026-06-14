#include "convert.hpp"
#include "merge.hpp"
#include "global.hpp"
#include "aa.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <optional>
#include <sys/stat.h>

static void usage(const char* prog) {
    std::cerr
        << "hoki — HOG codon Index\n"
        << "  " << prog << " convert -a ACC [-z LVL] [-p MINPID] [-e MAXEV] [-v] <in.tsv> <out.lhb>\n"
        << "  " << prog << " merge   [-z LVL=3] [--acc-registry file.acc] [--hog-range START END] <out.lhg> <out.lhgi> <input1> [input2 ...]\n"
        << "                  inputs may be .lhb or .lhg, mixed\n"
        << "  " << prog << " accregistry <out.acc> <shard.lhgi> [shard.lhgi ...]\n"
        << "  " << prog << " saav    <global.lhg> <global.lhgi> <HOG_ID> <POS> [AA]\n"
        << "  " << prog << " freq    <global.lhg> <global.lhgi> <HOG_ID>\n"
        << "  " << prog << " stat    <file.lhb | global.lhg [global.lhgi]>\n";
}

int main(int argc, char* argv[]) {
  try {
    if (argc < 2) { usage(argv[0]); return 1; }
    std::string mode = argv[1];

    if (mode == "convert") {
        lhi::ConvertOptions opts;
        std::string in_path, lhb_path;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-a" && i+1 < argc) opts.acc_id     = argv[++i];
            else if (a == "-z" && i+1 < argc) opts.zstd_level = std::stoi(argv[++i]);
            else if (a == "-p" && i+1 < argc) opts.min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue = std::stod(argv[++i]);
            else if (a == "-v")               opts.verbose    = true;
            else if (in_path.empty())         in_path         = a;
            else if (lhb_path.empty())        lhb_path        = a;
        }
        if (lhb_path.empty() || opts.acc_id.empty()) { usage(argv[0]); return 1; }
        lhi::convert(in_path, lhb_path, opts);
        return 0;
    }

    if (mode == "merge") {
        std::string acc_registry, hog_start, hog_end;
        int zstd_level = 3;
        std::vector<std::string> pos;  // out_lhg, out_lhgi, inputs...
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "--acc-registry" && i+1 < argc) acc_registry = argv[++i];
            else if (a == "--hog-range" && i+2 < argc) { hog_start = argv[++i]; hog_end = argv[++i]; }
            else if (a == "-z" && i+1 < argc)          zstd_level = std::stoi(argv[++i]);
            else pos.emplace_back(a);
        }
        if (pos.size() < 3) { usage(argv[0]); return 1; }
        std::string out_lhg  = pos[0];
        std::string out_lhgi = pos[1];
        std::vector<std::string> inputs(pos.begin() + 2, pos.end());
        lhi::merge_batches(inputs, out_lhg, out_lhgi, zstd_level, hog_start, hog_end, acc_registry);
        return 0;
    }

    if (mode == "accregistry") {
        if (argc < 4) { usage(argv[0]); return 1; }
        std::string out_acc = argv[2];
        std::vector<std::string> lhgis;
        for (int i = 3; i < argc; ++i) lhgis.emplace_back(argv[i]);
        lhi::build_acc_registry(lhgis, out_acc);
        return 0;
    }

    if (mode == "saav") {
        if (argc < 6) { usage(argv[0]); return 1; }
        std::string lhg_path  = argv[2];
        std::string lhgi_path = argv[3];
        std::string hog_id    = argv[4];
        uint32_t    pos       = uint32_t(std::stoul(argv[5]));
        std::optional<uint8_t> aa;
        if (argc >= 7) aa = lhi::encode_aa(argv[6][0]);

        lhi::GlobalIndex idx;
        if (!idx.load(lhgi_path) && !idx.load_from_lhg(lhg_path)) {
            std::cerr << "cannot load index from " << lhgi_path << "\n";
            return 1;
        }
        lhi::query_saav(lhg_path, idx, hog_id, pos, aa);
        return 0;
    }

    if (mode == "freq") {
        if (argc < 5) { usage(argv[0]); return 1; }
        std::string lhg_path  = argv[2];
        std::string lhgi_path = argv[3];
        std::string hog_id    = argv[4];

        lhi::GlobalIndex idx;
        if (!idx.load(lhgi_path) && !idx.load_from_lhg(lhg_path)) {
            std::cerr << "cannot load index from " << lhgi_path << "\n";
            return 1;
        }
        lhi::query_freq(lhg_path, idx, hog_id);
        return 0;
    }

    if (mode == "stat") {
        if (argc < 3) { usage(argv[0]); return 1; }
        std::string path = argv[2];

        // If .lhgi provided or path ends with .lhg → global stats.
        if (argc >= 4 || (path.size() > 4 && path.substr(path.size()-4) == ".lhg")) {
            std::string lhgi_path = (argc >= 4) ? argv[3]
                                  : path.substr(0, path.size()-4) + ".lhgi";
            lhi::GlobalIndex idx;
            bool ok = idx.load(lhgi_path);
            if (!ok) ok = idx.load_from_lhg(path);
            if (!ok) { std::cerr << "cannot load index from " << lhgi_path << "\n"; return 1; }

            struct stat st;
            if (::stat(path.c_str(), &st) == 0)
                std::fprintf(stdout, "file_size_mb\t%.1f\n", double(st.st_size) / (1024*1024));
            std::fprintf(stdout, "n_hogs\t%zu\n", idx.entries.size());
            std::fprintf(stdout, "n_accessions\t%zu\n", idx.accessions.size());

            uint64_t total_pairs = 0;
            for (auto& e : idx.entries) total_pairs += e.n_accessions;
            std::fprintf(stdout, "total_hog_accession_pairs\t%llu\n",
                         (unsigned long long)total_pairs);

            std::printf("hog_id\tn_accessions\tdata_length_kb\n");
            for (auto& e : idx.entries)
                std::printf("%s\t%u\t%.1f\n", e.hog_id.c_str(),
                            e.n_accessions,
                            double(e.data_length) / 1024.0);
        } else {
            // .lhb stats
            size_t n_blocks = 0;
            std::string acc;
            std::unordered_map<std::string, size_t> hog_blocks;
            lhi::scan_batch_file(path, 0,
                [&](const std::string& a) { acc = a; },
                [&](lhi::BatchBlockRef ref) {
                    ++n_blocks; hog_blocks[ref.hog_id]++;
                });
            std::printf("acc_id\t%s\n", acc.c_str());
            std::printf("n_hogs\t%zu\n", hog_blocks.size());
            std::printf("n_blocks\t%zu\n", n_blocks);
            std::printf("hog_id\tn_blocks\n");
            std::vector<std::pair<std::string,size_t>> hv(hog_blocks.begin(), hog_blocks.end());
            std::sort(hv.begin(), hv.end());
            for (auto& [h, n] : hv) std::printf("%s\t%zu\n", h.c_str(), n);
        }
        return 0;
    }

    usage(argv[0]);
    return 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
