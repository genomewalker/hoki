#include "convert.hpp"
#include "query.hpp"
#include <iostream>
#include <string>

static void usage(const char* prog) {
    std::cerr
        << "hoki — HOG codon Index\n"
        << "Usage:\n"
        << "  " << prog << " convert [-z LEVEL] [-b RECS] [-p MINPID] [-e MAXEV] [-v] <in.tsv> <container.lhc>\n"
        << "  " << prog << " varnt   [-p MINPID] [-e MAXEV] <container.lhc> <HOG_ID> [POS [OBS_AA]]\n"
        << "  " << prog << " tile    [-p MINPID] [-e MAXEV] [-n MIN_ACC] <container.lhc> <HOG_ID>\n"
        << "\n"
        << "convert: diamond blastx TSV → HOG-sharded container directory\n"
        << "  Column layout: qseqid qstart qend qlen qstrand sseqid sstart send slen\n"
        << "                 pident evalue cigar qseq_translated full_qseq\n"
        << "  Writes VarNT records (codon-resolved). 100%% pident alignments are skipped.\n"
        << "  Safe for concurrent writers: each HOG shard uses flock(LOCK_EX).\n"
        << "  -z LEVEL   zstd compression level (default 3)\n"
        << "  -b RECS    records per shard flush (default 50000)\n"
        << "  -p MINPID  minimum pident %% (default 0)\n"
        << "  -e MAXEV   max evalue (default 1.0)\n"
        << "  -v         verbose: print parse errors\n"
        << "\n"
        << "varnt: query codon-level variants for a HOG\n"
        << "  POS        1-based HOG position (omit for all positions)\n"
        << "  OBS_AA     filter to a specific amino acid letter\n"
        << "  Output TSV: contig_id\\thog_pos\\tobs_aa\\tcodon\\tpident\\tevalue\n"
        << "\n"
        << "tile: find tiling gaps across accessions for a HOG\n"
        << "  -n MIN_ACC   min accessions showing gap to report (default 2)\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) { usage(argv[0]); return 1; }
    std::string mode = argv[1];

    if (mode == "convert") {
        lhi::ConvertOptions opts;
        std::string in_path, container_dir;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-z" && i+1 < argc) opts.zstd_level = std::stoi(argv[++i]);
            else if (a == "-b" && i+1 < argc) opts.block_recs = std::stoul(argv[++i]);
            else if (a == "-p" && i+1 < argc) opts.min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue = std::stod(argv[++i]);
            else if (a == "-v")               opts.verbose    = true;
            else if (in_path.empty())         in_path       = a;
            else if (container_dir.empty())   container_dir = a;
        }
        if (container_dir.empty()) { usage(argv[0]); return 1; }
        lhi::convert(in_path, container_dir, opts);
        return 0;
    }

    if (mode == "varnt") {
        lhi::QueryOptions opts;
        std::string container_dir;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-p" && i+1 < argc) opts.min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue = std::stod(argv[++i]);
            else if (container_dir.empty())   container_dir = a;
            else if (opts.hog_id.empty())     opts.hog_id   = a;
            else if (!opts.pos)               opts.pos      = uint32_t(std::stoul(a));
            else                              opts.variant_aa = lhi::encode_aa(a[0]);
        }
        if (container_dir.empty() || opts.hog_id.empty()) { usage(argv[0]); return 1; }
        lhi::query_varnt(container_dir, opts);
        return 0;
    }

    if (mode == "tile") {
        float    min_pident = 0.0f;
        double   max_evalue = 1.0;
        uint32_t min_acc    = 2;
        std::string container_dir, hog_id;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-p" && i+1 < argc) min_pident    = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) max_evalue    = std::stod(argv[++i]);
            else if (a == "-n" && i+1 < argc) min_acc       = uint32_t(std::stoul(argv[++i]));
            else if (container_dir.empty())   container_dir = a;
            else if (hog_id.empty())          hog_id        = a;
        }
        if (container_dir.empty() || hog_id.empty()) { usage(argv[0]); return 1; }
        lhi::query_tile(container_dir, hog_id, min_pident, max_evalue, min_acc);
        return 0;
    }

    usage(argv[0]);
    return 1;
}
