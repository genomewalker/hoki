#pragma once
#include "format.hpp"
#include "cigar.hpp"
#include "aa.hpp"
#include <xxhash.h>
#include <zstd.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <limits>

namespace lhi {

struct ConvertOptions {
    int    zstd_level  = 3;
    size_t block_recs  = 50000;
    float  min_pident  = 0.0f;
    double max_evalue  = 1.0;
    bool   verbose     = false;
    bool   saav_only   = false;  // write Variant blocks only (no full AA arrays)
};

// diamond blastx outfmt 6 column indices (0-based):
//   qseqid qstart qend qlen qstrand sseqid sstart send slen pident evalue cigar qseq_translated full_qseq
//     0      1     2    3     4       5      6      7    8     9     10     11         12            13
namespace col {
    constexpr int qseqid   = 0, qstart = 1, qend = 2, qlen = 3, qstrand = 4;
    constexpr int sseqid   = 5, sstart = 6, send  = 7, slen = 8;
    constexpr int pident   = 9, evalue = 10, cigar = 11, qseq_aa = 12;
    constexpr int full_qseq = 13;  // nt contig — hashed (xxh3-128), not stored
}

inline std::string extract_hog(std::string_view sv) {
    auto p = sv.rfind('|');
    return std::string(p != std::string_view::npos ? sv.substr(p + 1) : sv);
}

inline int8_t make_qframe(std::string_view qstrand, uint32_t qstart) {
    uint32_t base = (qstart > 0 ? qstart - 1 : 0);
    int frame = static_cast<int>(base % 3) + 1;
    int sign  = (qstrand.empty() || qstrand[0] != '-') ? 1 : -1;
    return static_cast<int8_t>(sign * frame);
}

inline double parse_double(std::string_view sv) {
    char buf[64];
    size_t n = std::min(sv.size(), sizeof(buf) - 1);
    std::memcpy(buf, sv.data(), n);
    buf[n] = '\0';
    char* end;
    return std::strtod(buf, &end);
}

class BlockWriter {
    std::vector<uint8_t> raw_buf_;
    std::vector<AlignmentRecord> records_;
    std::unordered_map<std::string, uint32_t>& hog_dict_;
    std::unordered_map<std::string, uint32_t>& contig_dict_;
    std::vector<std::string>& hog_strings_;
    std::vector<std::string>& contig_strings_;
    std::FILE* out_;
    ConvertOptions opts_;
    size_t n_blocks_ = 0;

public:
    uint64_t n_written = 0, n_skipped = 0;

    BlockWriter(std::FILE* out,
                std::unordered_map<std::string, uint32_t>& hd,
                std::unordered_map<std::string, uint32_t>& cd,
                std::vector<std::string>& hs,
                std::vector<std::string>& cs,
                ConvertOptions opts)
        : hog_dict_(hd), contig_dict_(cd), hog_strings_(hs), contig_strings_(cs),
          out_(out), opts_(opts) {}

    uint32_t intern_hog(const std::string& s) {
        auto [it, ins] = hog_dict_.emplace(s, uint32_t(hog_strings_.size()));
        if (ins) hog_strings_.push_back(s);
        return it->second;
    }
    uint32_t intern_contig(const std::string& s) {
        auto [it, ins] = contig_dict_.emplace(s, uint32_t(contig_strings_.size()));
        if (ins) contig_strings_.push_back(s);
        return it->second;
    }

    void add(AlignmentRecord&& r) {
        records_.push_back(std::move(r));
        if (records_.size() >= opts_.block_recs) flush();
    }

    void flush() {
        if (records_.empty()) return;
        raw_buf_.clear();
        for (auto& r : records_) serialize_record(raw_buf_, r);

        BlockHeader bh{};
        bh.block_type     = uint8_t(BlockType::Alignments);
        bh.n_records      = uint32_t(records_.size());
        bh.min_hog_idx    = UINT32_MAX;
        bh.max_hog_idx    = 0;
        bh.min_sstart     = UINT32_MAX;
        bh.max_send       = 0;
        bh.min_pident     = 100.0f;
        bh.max_pident     = 0.0f;
        bh.min_evalue_log = 0.0f;
        bh.max_evalue_log = -300.0f;

        for (auto& r : records_) {
            bh.min_hog_idx    = std::min(bh.min_hog_idx,  r.hog_idx);
            bh.max_hog_idx    = std::max(bh.max_hog_idx,  r.hog_idx);
            bh.min_sstart     = std::min(bh.min_sstart,   r.sstart);
            bh.max_send       = std::max(bh.max_send,     r.send);
            bh.min_pident     = std::min(bh.min_pident,   r.pident);
            bh.max_pident     = std::max(bh.max_pident,   r.pident);
            float el = evalue_log(r.evalue);
            bh.min_evalue_log = std::max(bh.min_evalue_log, el);
            bh.max_evalue_log = std::min(bh.max_evalue_log, el);
        }
        records_.clear();

        size_t bound = ZSTD_compressBound(raw_buf_.size());
        std::vector<uint8_t> cbuf(bound);
        size_t csz = ZSTD_compress(cbuf.data(), bound, raw_buf_.data(), raw_buf_.size(), opts_.zstd_level);
        if (ZSTD_isError(csz))
            throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(csz));

        bh.raw_sz        = uint32_t(raw_buf_.size());
        bh.compressed_sz = uint32_t(csz);
        std::fwrite(&bh, sizeof(bh), 1, out_);
        std::fwrite(cbuf.data(), 1, csz, out_);
        ++n_blocks_;
    }

    size_t n_blocks() const { return n_blocks_; }
};

class VariantBlockWriter {
    std::vector<uint8_t>   raw_buf_;
    std::vector<VariantRecord> records_;
    std::FILE* out_;
    ConvertOptions opts_;
    size_t n_blocks_ = 0;

public:
    VariantBlockWriter(std::FILE* out, ConvertOptions opts) : out_(out), opts_(opts) {}

    void add(VariantRecord&& r) {
        records_.push_back(std::move(r));
        if (records_.size() >= opts_.block_recs) flush();
    }

    void flush() {
        if (records_.empty()) return;
        raw_buf_.clear();
        for (auto& r : records_) serialize_variant(raw_buf_, r);

        BlockHeader bh{};
        bh.block_type  = uint8_t(BlockType::Variants);
        bh.n_records   = uint32_t(records_.size());
        bh.min_hog_idx = UINT32_MAX; bh.max_hog_idx = 0;
        bh.min_sstart  = UINT32_MAX; bh.max_send    = 0;
        bh.min_pident  = 100.0f;     bh.max_pident  = 0.0f;
        bh.min_evalue_log = 0.0f;    bh.max_evalue_log = -300.0f;
        for (auto& r : records_) {
            bh.min_hog_idx    = std::min(bh.min_hog_idx, r.hog_idx);
            bh.max_hog_idx    = std::max(bh.max_hog_idx, r.hog_idx);
            bh.min_sstart     = std::min(bh.min_sstart,  r.sstart);
            bh.max_send       = std::max(bh.max_send,    r.send);
            bh.min_pident     = std::min(bh.min_pident,  r.pident);
            bh.max_pident     = std::max(bh.max_pident,  r.pident);
            float el = evalue_log(r.evalue);
            bh.min_evalue_log = std::max(bh.min_evalue_log, el);
            bh.max_evalue_log = std::min(bh.max_evalue_log, el);
        }
        records_.clear();

        size_t bound = ZSTD_compressBound(raw_buf_.size());
        std::vector<uint8_t> cbuf(bound);
        size_t csz = ZSTD_compress(cbuf.data(), bound, raw_buf_.data(), raw_buf_.size(), opts_.zstd_level);
        if (ZSTD_isError(csz))
            throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(csz));
        bh.raw_sz = uint32_t(raw_buf_.size()); bh.compressed_sz = uint32_t(csz);
        std::fwrite(&bh, sizeof(bh), 1, out_);
        std::fwrite(cbuf.data(), 1, csz, out_);
        ++n_blocks_;
    }

    size_t n_blocks() const { return n_blocks_; }
};

inline void write_dict_block(std::FILE* out, BlockType type,
                              const std::vector<std::string>& strings, int zstd_level) {
    std::string raw;
    raw.reserve(strings.size() * 32);
    for (auto& s : strings) { raw += s; raw += '\n'; }

    size_t bound = ZSTD_compressBound(raw.size());
    std::vector<uint8_t> cbuf(bound);
    size_t csz = ZSTD_compress(cbuf.data(), bound, raw.data(), raw.size(), zstd_level);
    if (ZSTD_isError(csz))
        throw std::runtime_error(std::string("zstd dict: ") + ZSTD_getErrorName(csz));

    BlockHeader bh{};
    bh.block_type    = uint8_t(type);
    bh.n_records     = uint32_t(strings.size());
    bh.raw_sz        = uint32_t(raw.size());
    bh.compressed_sz = uint32_t(csz);
    std::fwrite(&bh, sizeof(bh), 1, out);
    std::fwrite(cbuf.data(), 1, csz, out);
}

inline void convert(const std::string& tsv_path, const std::string& out_path,
                    ConvertOptions opts) {

    std::unordered_map<std::string, uint32_t> hog_dict, contig_dict;
    std::vector<std::string> hog_strings, contig_strings;

    std::string tmp_path  = out_path + ".blocks.tmp";
    std::string vtmp_path = out_path + ".vblocks.tmp";
    std::FILE* tmp  = std::fopen(tmp_path.c_str(), "wb");
    if (!tmp) throw std::runtime_error("cannot open tmp: " + tmp_path);
    std::FILE* vtmp = opts.saav_only ? std::fopen(vtmp_path.c_str(), "wb") : nullptr;
    if (opts.saav_only && !vtmp) throw std::runtime_error("cannot open vtmp: " + vtmp_path);

    BlockWriter        bw(tmp,  hog_dict, contig_dict, hog_strings, contig_strings, opts);
    VariantBlockWriter vw(vtmp ? vtmp : tmp, opts);

    std::istream* in = &std::cin;
    std::ifstream fin;
    if (tsv_path != "-") {
        fin.open(tsv_path);
        if (!fin) throw std::runtime_error("cannot open: " + tsv_path);
        in = &fin;
    }

    std::string line;
    uint64_t lineno = 0;

    while (std::getline(*in, line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;

        // Parse all 14 tab-separated fields; full_qseq is field[13]
        std::array<std::string_view, 14> f;
        size_t fi = 0;
        const char* p = line.data(), *ep = p + line.size(), *fp = p;
        while (p <= ep && fi < 14) {
            if (p == ep || *p == '\t') {
                f[fi++] = {fp, size_t(p - fp)};
                fp = p + 1;
            }
            ++p;
        }
        if (fi < 13) {
            if (opts.verbose)
                std::cerr << "skip line " << lineno << ": " << fi << " fields\n";
            ++bw.n_skipped; continue;
        }
        // field[13] may be absent (older diamond without full_qseq) — hash becomes all-zero
        bool has_nt = (fi == 14 && !f[col::full_qseq].empty());

        float  pident = static_cast<float>(parse_double(f[col::pident]));
        double ev     = parse_double(f[col::evalue]);

        if (pident < opts.min_pident) { ++bw.n_skipped; continue; }
        if (ev     > opts.max_evalue) { ++bw.n_skipped; continue; }

        uint32_t sstart, send, qstart, qend, qlen;
        try {
            sstart = uint32_t(std::stoul(std::string(f[col::sstart])));
            send   = uint32_t(std::stoul(std::string(f[col::send])));
            qstart = uint32_t(std::stoul(std::string(f[col::qstart])));
            qend   = uint32_t(std::stoul(std::string(f[col::qend])));
            qlen   = uint32_t(std::stoul(std::string(f[col::qlen])));
        } catch (const std::exception& e) {
            if (opts.verbose)
                std::cerr << "parse error line " << lineno << ": " << e.what() << "\n";
            ++bw.n_skipped; continue;
        }
        if (sstart > send) { ++bw.n_skipped; continue; }

        AlignmentRecord r;
        r.hog_idx    = bw.intern_hog(extract_hog(f[col::sseqid]));
        r.contig_idx = bw.intern_contig(std::string(f[col::qseqid]));
        r.sstart     = sstart;
        r.send       = send;
        r.qstart     = qstart;
        r.qend       = qend;
        r.qlen       = qlen;
        r.pident     = pident;
        r.evalue     = ev;
        r.qframe     = make_qframe(f[col::qstrand], qstart);

        // xxh3-128 of full_qseq — pointer-lossless fingerprint
        if (has_nt) {
            XXH128_hash_t h = XXH3_128bits(f[col::full_qseq].data(), f[col::full_qseq].size());
            // store as LE bytes: low 8 bytes first
            uint64_t lo = h.low64, hi = h.high64;
            for (int i = 0; i < 8; ++i) { r.nt_hash[i]   = uint8_t(lo >> (8*i)); }
            for (int i = 0; i < 8; ++i) { r.nt_hash[8+i] = uint8_t(hi >> (8*i)); }
        } else {
            std::memset(r.nt_hash, 0, 16);
        }

        // CIGAR parse — lossless: aas[] + inserts[]
        auto ar   = cigar_parse(f[col::cigar], f[col::qseq_aa], sstart, send);
        r.aas     = std::move(ar.aas);
        r.inserts = std::move(ar.inserts);

        if (opts.saav_only) {
            // build sparse variant record from M positions only
            VariantRecord vr;
            vr.contig_idx = r.contig_idx;
            vr.hog_idx    = r.hog_idx;
            vr.sstart     = r.sstart;
            vr.send       = r.send;
            vr.pident     = r.pident;
            vr.evalue     = r.evalue;
            vr.obs.reserve(r.aas.size());
            for (uint32_t i = 0; i < uint32_t(r.aas.size()); ++i) {
                uint8_t aa = r.aas[i];
                if (aa != AA_GAP && aa != AA_UNK)
                    vr.obs.push_back({r.sstart + i, aa});
            }
            vw.add(std::move(vr));
        } else {
            bw.add(std::move(r));
        }
        ++bw.n_written;
    }
    bw.flush();
    if (opts.saav_only) vw.flush();

    // Assemble: file header → HOG dict → contig dict → data blocks
    std::FILE* fout = std::fopen(out_path.c_str(), "wb");
    if (!fout) {
        std::fclose(tmp);
        if (vtmp) std::fclose(vtmp);
        throw std::runtime_error("cannot open: " + out_path);
    }

    size_t data_blocks = opts.saav_only ? vw.n_blocks() : bw.n_blocks();

    FileHeader fh{};
    std::memcpy(fh.magic, MAGIC, 8);
    fh.version  = FORMAT_VERSION;
    fh.flags    = FLAG_VALUE_LOSSLESS;
    fh.n_blocks = uint32_t(2 + data_blocks);
    std::fwrite(&fh, sizeof(fh), 1, fout);

    write_dict_block(fout, BlockType::HOGDict,    hog_strings,    opts.zstd_level);
    write_dict_block(fout, BlockType::ContigDict, contig_strings, opts.zstd_level);

    auto copy_tmp = [](std::FILE* src, std::FILE* dst, const std::string& p) {
        std::fclose(src);
        std::FILE* t = std::fopen(p.c_str(), "rb");
        if (!t) return;
        std::vector<uint8_t> buf(1 << 20);
        size_t n;
        while ((n = std::fread(buf.data(), 1, buf.size(), t)) > 0)
            std::fwrite(buf.data(), 1, n, dst);
        std::fclose(t);
    };

    if (opts.saav_only) {
        std::fclose(tmp); std::filesystem::remove(tmp_path);
        copy_tmp(vtmp, fout, vtmp_path);
        std::filesystem::remove(vtmp_path);
    } else {
        copy_tmp(tmp, fout, tmp_path);
        std::filesystem::remove(tmp_path);
    }
    std::fclose(fout);

    std::cerr << "wrote " << bw.n_written << " records, "
              << bw.n_skipped << " skipped, "
              << data_blocks << " blocks → " << out_path << "\n";
}

} // namespace lhi
