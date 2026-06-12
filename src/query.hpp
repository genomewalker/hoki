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

inline void query_variants(const std::string& lhi_path, const QueryOptions& opts) {
    LHIReader rdr(lhi_path);
    rdr.load_dicts();

    uint32_t hog_idx = rdr.find_hog(opts.hog_id);
    if (hog_idx == UINT32_MAX) {
        std::cerr << "HOG " << opts.hog_id << " not in " << lhi_path << "\n";
        return;
    }

    std::printf("contig_id\thog_pos\tobs_aa\tpident\tevalue\n");

    BlockHeader bh; std::vector<uint8_t> raw;
    while (rdr.read_block(bh, raw)) {
        if (BlockType(bh.block_type) != BlockType::Variants) continue;
        if (bh.min_hog_idx > hog_idx || bh.max_hog_idx < hog_idx) continue;
        if (bh.min_sstart > opts.pos  || bh.max_send   < opts.pos)  continue;

        const uint8_t* p = raw.data(), *end = p + raw.size();
        VariantRecord vr;
        while (p < end) {
            int n = deserialize_variant(p, end, vr);
            if (n <= 0) break;
            p += n;
            if (vr.hog_idx != hog_idx) continue;
            if (vr.pident < opts.min_pident) continue;
            if (vr.evalue > opts.max_evalue) continue;
            if (opts.pos < vr.sstart || opts.pos > vr.send) continue;
            for (auto& o : vr.obs) {
                if (o.hog_pos != opts.pos) continue;
                if (opts.variant_aa != AA_UNK && o.obs_aa != opts.variant_aa) continue;
                char aac = (o.obs_aa < 20) ? AA_ALPHA[o.obs_aa] : 'X';
                std::printf("%s\t%u\t%c\t%.2f\t%g\n",
                    rdr.contig_strings[vr.contig_idx].c_str(),
                    o.hog_pos, aac, vr.pident, vr.evalue);
            }
        }
    }
}


// query_varnt: scan VarNT blocks for a HOG position + optional AA filter
// Output TSV: contig_id \t hog_pos \t obs_aa \t codon \t pident \t evalue
inline void query_varnt(const std::string& lhi_path, const QueryOptions& opts) {
    LHIReader rdr(lhi_path);
    rdr.load_dicts();

    uint32_t hog_idx = rdr.find_hog(opts.hog_id);
    if (hog_idx == UINT32_MAX) {
        std::cerr << "HOG " << opts.hog_id << " not in " << lhi_path << "\n";
        return;
    }

    std::printf("contig_id\thog_pos\tobs_aa\tcodon\tpident\tevalue\n");

    BlockHeader bh; std::vector<uint8_t> raw;
    while (rdr.read_block(bh, raw)) {
        if (BlockType(bh.block_type) != BlockType::VarNT) continue;
        if (bh.min_hog_idx > hog_idx || bh.max_hog_idx < hog_idx) continue;
        if (opts.pos && (bh.min_sstart > opts.pos || bh.max_send < opts.pos)) continue;

        const uint8_t* p = raw.data(), *end = p + raw.size();
        VarNTRecord vr;
        while (p < end) {
            int n = deserialize_varnt(p, end, vr);
            if (n <= 0) break;
            p += n;
            if (vr.hog_idx != hog_idx) continue;
            if (vr.pident < opts.min_pident) continue;
            if (vr.evalue > opts.max_evalue) continue;
            for (auto& o : vr.vars) {
                uint32_t hog_pos = vr.sstart + o.hog_offset;
                if (opts.pos && hog_pos != opts.pos) continue;
                if (opts.variant_aa != AA_UNK && o.obs_aa != opts.variant_aa) continue;
                char aac = (o.obs_aa < 20) ? AA_ALPHA[o.obs_aa] : 'X';
                std::printf("%s\t%u\t%c\t%c%c%c\t%.2f\t%g\n",
                    rdr.contig_strings[vr.contig_idx].c_str(),
                    hog_pos, aac,
                    char(o.codon[0]), char(o.codon[1]), char(o.codon[2]),
                    vr.pident, vr.evalue);
            }
        }
    }
}

// query_tile: find tiling gaps across accessions for a HOG
// Groups VarNT records by accession prefix (strip after last '.' in contig_id)
// Output TSV: accession \t gap_start \t gap_end \t gap_aa_size \t n_exon_records
inline void query_tile(const std::string& lhi_path, const std::string& hog_id,
                        float min_pident, double max_evalue, uint32_t min_acc) {
    LHIReader rdr(lhi_path);
    rdr.load_dicts();

    uint32_t hog_idx = rdr.find_hog(hog_id);
    if (hog_idx == UINT32_MAX) {
        std::cerr << "HOG " << hog_id << " not in " << lhi_path << "\n";
        return;
    }

    // accession → sorted list of (sstart, send)
    std::unordered_map<std::string, std::vector<std::pair<uint32_t,uint32_t>>> acc_spans;

    BlockHeader bh; std::vector<uint8_t> raw;
    while (rdr.read_block(bh, raw)) {
        if (BlockType(bh.block_type) != BlockType::VarNT) continue;
        if (bh.min_hog_idx > hog_idx || bh.max_hog_idx < hog_idx) continue;

        const uint8_t* p = raw.data(), *end = p + raw.size();
        VarNTRecord vr;
        while (p < end) {
            int n = deserialize_varnt(p, end, vr);
            if (n <= 0) break;
            p += n;
            if (vr.hog_idx != hog_idx) continue;
            if (vr.pident < min_pident) continue;
            if (vr.evalue > max_evalue) continue;
            const std::string& contig = rdr.contig_strings[vr.contig_idx];
            auto dot = contig.rfind('.');
            std::string acc = (dot != std::string::npos) ? contig.substr(0, dot) : contig;
            acc_spans[acc].emplace_back(vr.sstart, vr.send);
        }
    }

    // gap positions → count of accessions showing a gap there
    std::unordered_map<uint64_t, uint32_t> gap_acc_count;  // key = (gap_start<<32)|gap_end

    std::printf("accession\tgap_start\tgap_end\tgap_aa_size\tn_exon_records\n");
    for (auto& [acc, spans] : acc_spans) {
        std::sort(spans.begin(), spans.end());
        for (size_t i = 1; i < spans.size(); ++i) {
            if (spans[i].first > spans[i-1].second + 1) {
                uint32_t gs = spans[i-1].second;
                uint32_t ge = spans[i].first;
                uint32_t gap_sz = ge - gs - 1;
                uint64_t key = (uint64_t(gs) << 32) | uint64_t(ge);
                gap_acc_count[key]++;
                if (gap_acc_count[key] >= min_acc) {
                    std::printf("%s\t%u\t%u\t%u\t%zu\n",
                        acc.c_str(), gs, ge, gap_sz, spans.size());
                }
            }
        }
    }

    // Summary: HOG position ranges with gaps seen in >= min_acc accessions
    std::printf("\n# Gap summary (position -> n_accessions)\n");
    std::printf("gap_start\tgap_end\tgap_aa_size\tn_accessions\n");
    std::vector<std::pair<uint64_t,uint32_t>> sorted_gaps(gap_acc_count.begin(), gap_acc_count.end());
    std::sort(sorted_gaps.begin(), sorted_gaps.end(), [](auto& a, auto& b){ return a.first < b.first; });
    for (auto& [key, cnt] : sorted_gaps) {
        if (cnt < min_acc) continue;
        uint32_t gs = uint32_t(key >> 32);
        uint32_t ge = uint32_t(key & 0xFFFFFFFF);
        std::printf("%u\t%u\t%u\t%u\n", gs, ge, ge - gs - 1, cnt);
    }
}

} // namespace lhi
