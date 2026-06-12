#include "convert.hpp"
#include "query.hpp"
#include <iostream>
#include <string>
#include <cstring>

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " convert [-z LEVEL] [-b RECS] [-p MINPID] [-e MAXEQ] [-s SAMPLE] <in.tsv> <out.lhi>\n"
        << "  " << prog << " query   [-p MINPID] [-e MAXEQ] [-S] <in.lhi> <HOG_ID> [POS [VARIANT_AA]]\n"
        << "  " << prog << " stats   <in.lhi>\n"
        << "\n"
        << "convert: parse Logan DIAMOND TSV → LHI binary (reads stdin if in.tsv is '-')\n"
        << "  -z LEVEL   zstd level (default 3)\n"
        << "  -b RECS    records per block (default 50000)\n"
        << "  -p MINPID  minimum pident (default 0)\n"
        << "  -e MAXEQ   max log10(evalue) integer to keep (default 0 = all)\n"
        << "  -s SAMPLE  override sample accession (default: parse from qseqid)\n"
        << "\n"
        << "query: print per-position AA frequency table for one HOG\n"
        << "  POS        1-based HOG position (omit for all positions)\n"
        << "  VARIANT_AA single AA letter to list samples (requires POS)\n"
        << "  -S         list sample IDs for variant matches\n"
        << "\n"
        << "stats: print block metadata summary\n";
}

int main(int argc, char* argv[]) {
    if (argc < 2) { usage(argv[0]); return 1; }
    std::string mode = argv[1];

    if (mode == "convert") {
        lhi::ConvertOptions opts;
        std::string in_path, out_path, sample_override;

        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-z" && i+1 < argc) { opts.zstd_level = std::stoi(argv[++i]); }
            else if (a == "-b" && i+1 < argc) { opts.block_recs = std::stoul(argv[++i]); }
            else if (a == "-p" && i+1 < argc) { opts.min_pident = std::stof(argv[++i]); }
            else if (a == "-e" && i+1 < argc) { opts.max_evalue_q = std::stoi(argv[++i]); }
            else if (a == "-s" && i+1 < argc) { sample_override = argv[++i]; }
            else if (a == "-v") { opts.verbose = true; }
            else if (in_path.empty())  { in_path = a; }
            else if (out_path.empty()) { out_path = a; }
        }
        if (out_path.empty()) { usage(argv[0]); return 1; }
        lhi::convert(in_path, out_path, sample_override, opts);
        return 0;
    }

    if (mode == "query") {
        lhi::QueryOptions opts;
        std::string lhi_path;

        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-p" && i+1 < argc) { opts.min_pident = std::stof(argv[++i]); }
            else if (a == "-e" && i+1 < argc) { opts.max_evalue_q = std::stoi(argv[++i]); }
            else if (a == "-S") { opts.list_samples = true; }
            else if (lhi_path.empty())     { lhi_path = a; }
            else if (opts.hog_id.empty())  { opts.hog_id = a; }
            else if (!opts.pos)            { opts.pos = static_cast<uint32_t>(std::stoul(a)); }
            else { opts.variant_aa = lhi::encode_aa(a[0]); }
        }
        if (lhi_path.empty() || opts.hog_id.empty()) { usage(argv[0]); return 1; }
        lhi::query_position(lhi_path, opts);
        return 0;
    }

    if (mode == "stats") {
        if (argc < 3) { usage(argv[0]); return 1; }
        lhi::LHIReader rdr(argv[2]);
        std::cout << "HOG dict: loading...\n";
        rdr.load_dicts();
        std::cout << "HOGs:    " << rdr.hog_strings.size() << "\n";
        std::cout << "Samples: " << rdr.sample_strings.size() << "\n";
        std::cout << "Blocks:  " << rdr.n_blocks << "\n";
        lhi::BlockHeader bh;
        std::vector<uint8_t> raw;
        uint64_t total_recs = 0;
        int aln_blocks = 0;
        while (rdr.read_block(bh, raw)) {
            if (static_cast<lhi::BlockType>(bh.block_type) == lhi::BlockType::Alignments) {
                total_recs += bh.n_records;
                ++aln_blocks;
            }
        }
        std::cout << "Aln blocks: " << aln_blocks << "\n";
        std::cout << "Total recs: " << total_recs << "\n";
        return 0;
    }

    usage(argv[0]);
    return 1;
}
