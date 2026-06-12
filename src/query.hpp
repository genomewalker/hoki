#pragma once
#include "format.hpp"
#include "aa.hpp"
#include <zstd.h>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cstdio>

namespace lhi {

struct PosFreq {
    uint32_t hog_idx;
    uint32_t pos;      // 1-based HOG position
    uint32_t counts[22]{};  // AA_TABLE index 0-19, GAP=20, UNK=21
    uint32_t total = 0;

    void add(uint8_t aa) {
        uint8_t idx = (aa == AA_GAP) ? 20 : (aa == AA_UNK ? 21 : aa);
        if (idx < 22) counts[idx]++;
        ++total;
    }
};

struct QueryOptions {
    std::string hog_id;     // required: e.g. "N0.HOG0047149"
    uint32_t    pos;        // 1-based HOG position (0 = all positions)
    uint8_t     variant_aa; // AA_UNK = any non-reference; specific = filter to this AA
    float       min_pident = 0.0f;
    int8_t      max_evalue_q = 0;
    bool        list_samples = false;  // print sample IDs for matching hits
};

// Load dicts from open file positioned at start (after FileHeader)
struct LHIReader {
    std::vector<std::string> hog_strings;
    std::vector<std::string> sample_strings;
    std::FILE* f;
    uint32_t n_blocks;
    int blocks_read = 0;

    explicit LHIReader(const std::string& path) {
        f = std::fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open " + path);

        FileHeader fh;
        if (std::fread(&fh, sizeof(fh), 1, f) != 1)
            throw std::runtime_error("short file header");
        if (std::memcmp(fh.magic, MAGIC, 8) != 0)
            throw std::runtime_error("bad magic");
        if (fh.version != FORMAT_VERSION)
            throw std::runtime_error("unsupported version");
        n_blocks = fh.n_blocks;
    }

    ~LHIReader() { if (f) std::fclose(f); }

    // Read and decompress one block header + payload. Returns false at EOF.
    bool read_block(BlockHeader& bh, std::vector<uint8_t>& raw) {
        if (std::fread(&bh, sizeof(bh), 1, f) != 1) return false;
        std::vector<uint8_t> cbuf(bh.compressed_sz);
        if (std::fread(cbuf.data(), 1, bh.compressed_sz, f) != bh.compressed_sz)
            throw std::runtime_error("truncated block");
        raw.resize(bh.raw_sz);
        size_t r = ZSTD_decompress(raw.data(), bh.raw_sz, cbuf.data(), bh.compressed_sz);
        if (ZSTD_isError(r))
            throw std::runtime_error(std::string("zstd decomp: ") + ZSTD_getErrorName(r));
        ++blocks_read;
        return true;
    }

    void load_dicts() {
        BlockHeader bh;
        std::vector<uint8_t> raw;
        // HOG dict
        if (!read_block(bh, raw))
            throw std::runtime_error("missing HOG dict block");
        parse_dict(raw, hog_strings);
        // Sample dict
        if (!read_block(bh, raw))
            throw std::runtime_error("missing sample dict block");
        parse_dict(raw, sample_strings);
    }

    uint32_t find_hog(const std::string& id) const {
        for (uint32_t i = 0; i < hog_strings.size(); ++i)
            if (hog_strings[i] == id) return i;
        return UINT32_MAX;
    }

private:
    static void parse_dict(const std::vector<uint8_t>& raw, std::vector<std::string>& out) {
        out.clear();
        const char* p = reinterpret_cast<const char*>(raw.data());
        const char* end = p + raw.size();
        while (p < end) {
            const char* nl = static_cast<const char*>(std::memchr(p, '\n', end - p));
            if (!nl) nl = end;
            out.emplace_back(p, nl - p);
            p = nl + 1;
        }
    }
};

// Query LHI file for AA distribution at a specific HOG position
inline void query_position(const std::string& lhi_path, const QueryOptions& opts) {
    LHIReader rdr(lhi_path);
    rdr.load_dicts();

    uint32_t hog_idx = rdr.find_hog(opts.hog_id);
    if (hog_idx == UINT32_MAX) {
        std::cerr << "HOG " << opts.hog_id << " not found in " << lhi_path << "\n";
        return;
    }

    std::unordered_map<uint32_t, PosFreq> freqs;  // pos → freq
    std::vector<std::string> hit_samples;

    BlockHeader bh;
    std::vector<uint8_t> raw;
    while (rdr.read_block(bh, raw)) {
        if (static_cast<BlockType>(bh.block_type) != BlockType::Alignments) continue;

        // Block-level skip: HOG range
        if (bh.min_hog_idx > hog_idx || bh.max_hog_idx < hog_idx) continue;
        // Block-level skip: position range
        if (opts.pos && (bh.min_sstart > opts.pos || bh.max_send < opts.pos)) continue;
        // Block-level skip: evalue
        if (bh.max_evalue_q > opts.max_evalue_q) continue;  // best evalue in block still too weak

        const uint8_t* p = raw.data();
        const uint8_t* end = p + raw.size();
        AlignmentRecord r;
        while (p < end) {
            int consumed = deserialize_record(p, end, r);
            if (consumed <= 0) break;
            p += consumed;

            if (r.hog_idx != hog_idx) continue;
            if (r.pident_q < quantize_pident(opts.min_pident)) continue;
            if (r.evalue_q > opts.max_evalue_q) continue;

            if (opts.pos) {
                if (opts.pos < r.sstart || opts.pos > r.send) continue;
                uint8_t aa = r.aas[opts.pos - r.sstart];
                if (aa == AA_GAP || aa == AA_UNK) continue;
                freqs[opts.pos].hog_idx = hog_idx;
                freqs[opts.pos].pos = opts.pos;
                freqs[opts.pos].add(aa);
                if (opts.list_samples && aa == opts.variant_aa)
                    hit_samples.push_back(rdr.sample_strings[r.sample_idx]);
            } else {
                // all covered positions
                for (uint32_t i = 0; i < r.aa_span(); ++i) {
                    uint8_t aa = r.aas[i];
                    uint32_t p2 = r.sstart + i;
                    freqs[p2].hog_idx = hog_idx;
                    freqs[p2].pos = p2;
                    freqs[p2].add(aa);
                }
            }
        }
    }

    // output TSV
    std::cout << "hog_id\tpos\taa\tcount\ttotal\tfreq\n";
    for (auto& [pos, pf] : freqs) {
        for (uint8_t a = 0; a < 22; ++a) {
            if (pf.counts[a] == 0) continue;
            char aac = (a < 20) ? AA_ALPHA[a] : (a == 20 ? '-' : 'X');
            double freq = pf.total ? static_cast<double>(pf.counts[a]) / pf.total : 0.0;
            std::printf("%s\t%u\t%c\t%u\t%u\t%.4f\n",
                        opts.hog_id.c_str(), pos, aac, pf.counts[a], pf.total, freq);
        }
    }
    if (opts.list_samples && !hit_samples.empty()) {
        std::cerr << "samples with variant:\n";
        for (auto& s : hit_samples) std::cerr << "  " << s << "\n";
    }
}

} // namespace lhi
