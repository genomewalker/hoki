#pragma once
#include "format.hpp"
#include "aa.hpp"
#include <zstd.h>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <cstdio>
#include <cstring>

namespace lhi {

struct PosFreq {
    uint32_t hog_idx;
    uint32_t pos;
    uint32_t counts[22]{};  // 0-19=AA, 20=gap, 21=unk
    uint32_t total = 0;

    void add(uint8_t aa) {
        uint8_t i = (aa == AA_GAP) ? 20 : (aa == AA_UNK ? 21 : aa);
        if (i < 22) counts[i]++;
        ++total;
    }
};

struct QueryOptions {
    std::string hog_id;
    uint32_t    pos          = 0;
    uint8_t     variant_aa   = AA_UNK;
    float       min_pident   = 0.0f;
    double      max_evalue   = 1.0;
    bool        list_contigs = false;
};

struct LHIReader {
    std::vector<std::string> hog_strings;
    std::vector<std::string> contig_strings;
    std::FILE* f;
    uint32_t n_blocks;

    explicit LHIReader(const std::string& path) {
        f = std::fopen(path.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot open: " + path);
        FileHeader fh;
        if (std::fread(&fh, sizeof(fh), 1, f) != 1)
            throw std::runtime_error("short header");
        if (std::memcmp(fh.magic, MAGIC, 8) != 0)
            throw std::runtime_error("bad magic — not an LHI v3 file");
        if (fh.version != FORMAT_VERSION)
            throw std::runtime_error("unsupported version " + std::to_string(fh.version));
        n_blocks = fh.n_blocks;
    }
    ~LHIReader() { if (f) std::fclose(f); }

    bool read_block(BlockHeader& bh, std::vector<uint8_t>& raw) {
        if (std::fread(&bh, sizeof(bh), 1, f) != 1) return false;
        std::vector<uint8_t> cbuf(bh.compressed_sz);
        if (std::fread(cbuf.data(), 1, bh.compressed_sz, f) != bh.compressed_sz)
            throw std::runtime_error("truncated block");
        raw.resize(bh.raw_sz);
        size_t r = ZSTD_decompress(raw.data(), bh.raw_sz, cbuf.data(), bh.compressed_sz);
        if (ZSTD_isError(r))
            throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(r));
        return true;
    }

    void load_dicts() {
        BlockHeader bh; std::vector<uint8_t> raw;
        if (!read_block(bh, raw)) throw std::runtime_error("missing HOG dict");
        parse_dict(raw, hog_strings);
        if (!read_block(bh, raw)) throw std::runtime_error("missing contig dict");
        parse_dict(raw, contig_strings);
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
        const char* e = p + raw.size();
        while (p < e) {
            const char* nl = static_cast<const char*>(std::memchr(p, '\n', e - p));
            if (!nl) nl = e;
            out.emplace_back(p, nl - p);
            p = nl + 1;
        }
    }
};

inline void query_position(const std::string& lhi_path, const QueryOptions& opts) {
    LHIReader rdr(lhi_path);
    rdr.load_dicts();

    uint32_t hog_idx = rdr.find_hog(opts.hog_id);
    if (hog_idx == UINT32_MAX) {
        std::cerr << "HOG " << opts.hog_id << " not in " << lhi_path << "\n";
        return;
    }

    std::unordered_map<uint32_t, PosFreq> freqs;
    std::vector<std::string> hit_contigs;

    BlockHeader bh; std::vector<uint8_t> raw;
    while (rdr.read_block(bh, raw)) {
        if (BlockType(bh.block_type) != BlockType::Alignments) continue;
        if (bh.min_hog_idx > hog_idx || bh.max_hog_idx < hog_idx) continue;
        if (opts.pos && (bh.min_sstart > opts.pos || bh.max_send < opts.pos)) continue;
        if (opts.max_evalue < 1.0) {
            float log_thresh = evalue_log(opts.max_evalue);
            if (bh.max_evalue_log > log_thresh) continue;
        }

        const uint8_t* p = raw.data(), *end = p + raw.size();
        AlignmentRecord r;
        while (p < end) {
            int n = deserialize_record(p, end, r);
            if (n <= 0) break;
            p += n;
            if (r.hog_idx != hog_idx) continue;
            if (r.pident < opts.min_pident) continue;
            if (r.evalue > opts.max_evalue) continue;

            if (opts.pos) {
                if (opts.pos < r.sstart || opts.pos > r.send) continue;
                uint8_t aa = r.aas[opts.pos - r.sstart];
                if (aa == AA_GAP || aa == AA_UNK) continue;
                auto& pf = freqs[opts.pos];
                pf.hog_idx = hog_idx; pf.pos = opts.pos;
                pf.add(aa);
                if (opts.list_contigs && aa == opts.variant_aa)
                    hit_contigs.push_back(rdr.contig_strings[r.contig_idx]);
            } else {
                for (uint32_t i = 0; i < r.aa_span(); ++i) {
                    uint8_t aa = r.aas[i];
                    if (aa == AA_UNK) continue;
                    uint32_t p2 = r.sstart + i;
                    auto& pf = freqs[p2];
                    pf.hog_idx = hog_idx; pf.pos = p2;
                    pf.add(aa);
                }
            }
        }
    }

    std::printf("hog_id\tpos\taa\tcount\ttotal\tfreq\n");
    for (auto& [pos, pf] : freqs) {
        for (uint8_t a = 0; a < 22; ++a) {
            if (!pf.counts[a]) continue;
            char aac = (a < 20) ? AA_ALPHA[a] : (a == 20 ? '-' : 'X');
            std::printf("%s\t%u\t%c\t%u\t%u\t%.6f\n",
                opts.hog_id.c_str(), pos, aac,
                pf.counts[a], pf.total,
                pf.total ? double(pf.counts[a]) / pf.total : 0.0);
        }
    }
    if (opts.list_contigs && !hit_contigs.empty()) {
        std::fprintf(stderr, "contigs with variant:\n");
        for (auto& c : hit_contigs) std::fprintf(stderr, "  %s\n", c.c_str());
    }
}

} // namespace lhi
