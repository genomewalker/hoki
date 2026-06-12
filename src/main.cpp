#include "convert.hpp"
#include "query.hpp"
#include <iostream>
#include <string>

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " convert [-z LEVEL] [-b RECS] [-p MINPID] [-e MAXEV] <in.tsv> <out.lhi>\n"
        << "  " << prog << " query   [-p MINPID] [-e MAXEV] [-C] <in.lhi> <HOG_ID> [POS [VARIANT_AA]]\n"
        << "  " << prog << " stats   <in.lhi>\n"
        << "\n"
        << "convert: diamond blastx TSV → LHI v3 (stdin if in.tsv is '-')\n"
        << "  Column layout: qseqid qstart qend qlen qstrand sseqid sstart send slen\n"
        << "                 pident evalue cigar qseq_translated full_qseq\n"
        << "  full_qseq hashed (xxh3-128) for pointer-lossless verification; not stored.\n"
        << "  -z LEVEL   zstd compression level (default 3, max 19)\n"
        << "  -b RECS    records per block (default 50000)\n"
        << "  -p MINPID  minimum pident % (default 0)\n"
        << "  -e MAXEV   max evalue (default 1.0, e.g. 1e-5)\n"
        << "  -v         verbose: print parse errors\n"
        << "  --saav     write sparse Variant blocks only (compact, for DB import)\n"
        << "  --varnt    write VarNT blocks (all aligned AA + codon per position)\n"
        << "\n"
        << "varnt: query VarNT blocks for a HOG position\n"
        << "  lhi varnt [-p MINPID] [-e MAXEV] <in.lhi> <HOG_ID> [POS [OBS_AA]]\n"
        << "  Outputs TSV: contig_id\\thog_pos\\tobs_aa\\tcodon\\tpident\\tevalue\n"
        << "\n"
        << "tile: find tiling gaps across accessions for a HOG\n"
        << "  lhi tile [-p MINPID] [-e MAXEV] [-n MIN_ACC] <in.lhi> <HOG_ID>\n"
        << "  -n MIN_ACC   min accessions showing gap to report (default 2)\n"
        << "\n"
        << "query: AA frequency table at a HOG position\n"
        << "  POS        1-based HOG position (omit for all positions)\n"
        << "  VARIANT_AA AA letter to collect contig IDs for (requires POS and -C)\n"
        << "  -C         list contig IDs for variant matches\n"
        << "  -p/-e      same filters as convert\n"
        << "\n"
        << "stats: block summary (no decompression)\n"
        << "\n"
        << "saav: query Variant blocks for a specific (HOG, position, AA)\n"
        << "  lhi saav <in.lhi> <HOG_ID> <POS> [OBS_AA]\n"
        << "  Outputs TSV: contig_id\\thog_pos\\tobs_aa\\tpident\\tevalue\n"
        << "  POS        1-based HOG position\n"
        << "  OBS_AA     filter to a specific amino acid letter (optional)\n"
        << "  -p/-e      same filters as convert\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) { usage(argv[0]); return 1; }
    std::string mode = argv[1];

    if (mode == "convert") {
        lhi::ConvertOptions opts;
        std::string in_path, out_path;

        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-z" && i+1 < argc) opts.zstd_level = std::stoi(argv[++i]);
            else if (a == "-b" && i+1 < argc) opts.block_recs = std::stoul(argv[++i]);
            else if (a == "-p" && i+1 < argc) opts.min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue = std::stod(argv[++i]);
            else if (a == "-v")               opts.verbose = true;
            else if (a == "--saav")           opts.saav_only = true;
            else if (a == "--varnt")          opts.varnt = true;
            else if (in_path.empty())         in_path  = a;
            else if (out_path.empty())        out_path = a;
        }
        if (out_path.empty()) { usage(argv[0]); return 1; }
        lhi::convert(in_path, out_path, opts);
        return 0;
    }

    if (mode == "query") {
        lhi::QueryOptions opts;
        std::string lhi_path;

        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-p" && i+1 < argc) opts.min_pident   = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue   = std::stod(argv[++i]);
            else if (a == "-C")               opts.list_contigs = true;
            else if (lhi_path.empty())        lhi_path   = a;
            else if (opts.hog_id.empty())     opts.hog_id = a;
            else if (!opts.pos)               opts.pos    = uint32_t(std::stoul(a));
            else                              opts.variant_aa = lhi::encode_aa(a[0]);
        }
        if (lhi_path.empty() || opts.hog_id.empty()) { usage(argv[0]); return 1; }
        lhi::query_position(lhi_path, opts);
        return 0;
    }

    if (mode == "saav") {
        if (argc < 5) { usage(argv[0]); return 1; }
        lhi::QueryOptions opts;
        std::string lhi_path;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-p" && i+1 < argc) opts.min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue = std::stod(argv[++i]);
            else if (lhi_path.empty())        lhi_path    = a;
            else if (opts.hog_id.empty())     opts.hog_id = a;
            else if (!opts.pos)               opts.pos    = uint32_t(std::stoul(a));
            else                              opts.variant_aa = lhi::encode_aa(a[0]);
        }
        if (lhi_path.empty() || opts.hog_id.empty() || !opts.pos) { usage(argv[0]); return 1; }
        lhi::query_variants(lhi_path, opts);
        return 0;
    }

    if (mode == "varnt") {
        lhi::QueryOptions opts;
        std::string lhi_path;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-p" && i+1 < argc) opts.min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue = std::stod(argv[++i]);
            else if (lhi_path.empty())        lhi_path    = a;
            else if (opts.hog_id.empty())     opts.hog_id = a;
            else if (!opts.pos)               opts.pos    = uint32_t(std::stoul(a));
            else                              opts.variant_aa = lhi::encode_aa(a[0]);
        }
        if (lhi_path.empty() || opts.hog_id.empty()) { usage(argv[0]); return 1; }
        lhi::query_varnt(lhi_path, opts);
        return 0;
    }

    if (mode == "tile") {
        float    min_pident = 0.0f;
        double   max_evalue = 1.0;
        uint32_t min_acc    = 2;
        std::string lhi_path, hog_id;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-p" && i+1 < argc) min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) max_evalue = std::stod(argv[++i]);
            else if (a == "-n" && i+1 < argc) min_acc    = uint32_t(std::stoul(argv[++i]));
            else if (lhi_path.empty())        lhi_path   = a;
            else if (hog_id.empty())          hog_id     = a;
        }
        if (lhi_path.empty() || hog_id.empty()) { usage(argv[0]); return 1; }
        lhi::query_tile(lhi_path, hog_id, min_pident, max_evalue, min_acc);
        return 0;
    }

    if (mode == "stats") {
        if (argc < 3) { usage(argv[0]); return 1; }
        lhi::LHIReader rdr(argv[2]);
        rdr.load_dicts();
        std::cout << "HOGs:    " << rdr.hog_strings.size()   << "\n";
        std::cout << "Contigs: " << rdr.contig_strings.size() << "\n";
        std::cout << "Blocks:  " << rdr.n_blocks << "\n";
        lhi::BlockHeader bh; std::vector<uint8_t> raw;
        uint64_t total_recs = 0; int aln_blocks = 0;
        while (rdr.read_block(bh, raw)) {
            if (lhi::BlockType(bh.block_type) == lhi::BlockType::Alignments) {
                total_recs += bh.n_records; ++aln_blocks;
            }
        }
        std::cout << "Aln blocks: " << aln_blocks  << "\n";
        std::cout << "Total recs: " << total_recs << "\n";
        return 0;
    }

    usage(argv[0]);
    return 1;
}
