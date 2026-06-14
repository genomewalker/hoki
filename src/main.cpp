#include "convert.hpp"
#include "merge.hpp"
#include "global.hpp"
#include "aa.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sys/stat.h>

static void usage(const char* prog) {
    std::cerr
        << "hoki — HOG codon Index (v5 position-centric inverted format)\n"
        << "Usage:\n"
        << "  " << prog << " convert -a ACC [-z LVL] [-b RECS] [-p MINPID] [-e MAXEV] [-v] <in.tsv> <out.lhb>\n"
        << "  " << prog << " merge   <out.lhg> <out.lhgi> <batch1.lhb> [batch2.lhb ...]\n"
        << "  " << prog << " saav    <global.lhg> <global.lhgi> <HOG_ID> <POS> [AA]\n"
        << "  " << prog << " freq    <global.lhg> <global.lhgi> <HOG_ID>\n"
        << "  " << prog << " stat    <file.lhb>   — show HOG/block counts for a batch file\n"
        << "  " << prog << " stat    <global.lhg> [<global.lhgi>] — list all HOGs with offsets\n"
        << "\n"
        << "convert: diamond blastx TSV → .lhb batch file (one per accession, v4 format)\n"
        << "  Column layout: qseqid qstart qend qlen qstrand sseqid sstart send slen\n"
        << "                 pident evalue cigar qseq_translated full_qseq\n"
        << "  Stores codon-resolved VarNT records. 100%% pident alignments skipped.\n"
        << "  -a ACC   accession ID (SRA/ENA run ID, required)\n"
        << "  -z LVL   zstd compression level (default 9)\n"
        << "  -b RECS  records per HOG flush (default 50000)\n"
        << "  -p PCT   minimum pident %% (default 0)\n"
        << "  -e EV    max evalue (default 1.0)\n"
        << "  -v       verbose parse errors\n"
        << "\n"
        << "merge: merge N .lhb batch files → one v5 inverted .lhg + .lhgi index\n"
        << "  Inverts per-accession blocks into position-centric records.\n"
        << "  HOGs are sorted lexicographically; .lhgi carries the accession registry.\n"
        << "\n"
        << "saav: query all accessions observed at HOG position POS\n"
        << "  POS   HOG position (matches r.sstart + obs offset)\n"
        << "  AA    optional amino-acid letter filter\n"
        << "  Output TSV: acc_id\\tunitig_id\\thog_pos\\tobs_aa\\tcodon\n"
        << "\n"
        << "freq: per-position codon frequency table for a HOG\n"
        << "  Output TSV: hog_pos\\tcodon\\tobs_aa\\tn_accessions\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) { usage(argv[0]); return 1; }
    std::string mode = argv[1];

    if (mode == "convert") {
        lhi::ConvertOptions opts;
        std::string in_path, lhb_path;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-a" && i+1 < argc) opts.acc_id     = argv[++i];
            else if (a == "-z" && i+1 < argc) opts.zstd_level = std::stoi(argv[++i]);
            else if (a == "-b" && i+1 < argc) opts.block_recs = std::stoul(argv[++i]);
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
        if (argc < 5) { usage(argv[0]); return 1; }
        std::string out_lhg  = argv[2];
        std::string out_lhgi = argv[3];
        std::vector<std::string> batches;
        for (int i = 4; i < argc; ++i) batches.emplace_back(argv[i]);
        if (batches.empty()) { usage(argv[0]); return 1; }
        lhi::merge_batches(batches, out_lhg, out_lhgi);
        return 0;
    }

    if (mode == "saav") {
        // hoki saav <global.lhg> <global.lhgi> <HOG_ID> <POS> [AA]
        if (argc < 6) { usage(argv[0]); return 1; }
        std::string lhg_path  = argv[2];
        std::string lhgi_path = argv[3];
        std::string hog_id    = argv[4];
        uint32_t    pos       = uint32_t(std::stoul(argv[5]));
        uint8_t     aa        = lhi::AA_UNK;
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
        // hoki freq <global.lhg> <global.lhgi> <HOG_ID>
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
        // hoki stat <file.lhb|global.lhg> [<global.lhgi>]
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
}
