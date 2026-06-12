#pragma once
#include "format.hpp"
#include "cigar.hpp"
#include "aa.hpp"
#include <zstd.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <cstdio>

namespace lhi {

struct ConvertOptions {
    int zstd_level   = 3;
    size_t block_recs = 50000;  // records per block before flush
    float min_pident  = 0.0f;
    int   max_evalue_q = 0;     // max log10(evalue) to keep (0 = keep all)
    bool  verbose     = false;
};

// Extract SRA accession from qseqid: "DRR000001_1_ka:..." → "DRR000001"
inline std::string extract_sample(std::string_view qseqid) {
    auto pos = qseqid.find('_');
    return std::string(pos != std::string_view::npos ? qseqid.substr(0, pos) : qseqid);
}

// Extract HOG id from sseqid: "panbarley.bpgv2|N0.HOG0047149" → "N0.HOG0047149"
inline std::string extract_hog(std::string_view sseqid) {
    auto pos = sseqid.rfind('|');
    return std::string(pos != std::string_view::npos ? sseqid.substr(pos + 1) : sseqid);
}

// Extract qframe from Logan qseqid: "DRR000001_1_ka:f:17.750____L:-:550:+"
// Format appears to be: ...L:<strand>:<start>:<frame>
// Returns signed frame: +1/+2/+3 or -1/-2/-3, 0 if unknown
inline int8_t extract_qframe(std::string_view qseqid) {
    // find "L:" or "l:" segment
    auto lpos = qseqid.rfind("L:");
    if (lpos == std::string_view::npos) return 0;
    auto seg = qseqid.substr(lpos + 2);
    // seg = "-:550:+" or "+:550:1" etc.
    if (seg.empty()) return 0;
    int sign = (seg[0] == '-') ? -1 : 1;
    // find last ':' for frame
    auto fpos = seg.rfind(':');
    if (fpos == std::string_view::npos) return 0;
    auto fstr = seg.substr(fpos + 1);
    if (fstr.empty()) return 0;
    int frame = fstr[0] - '0';
    if (frame < 1 || frame > 3) return 0;
    return static_cast<int8_t>(sign * frame);
}

class BlockWriter {
    std::vector<uint8_t> raw_buf_;
    std::vector<AlignmentRecord> records_;
    std::unordered_map<std::string, uint32_t>& hog_dict_;
    std::unordered_map<std::string, uint32_t>& sample_dict_;
    std::vector<std::string>& hog_strings_;
    std::vector<std::string>& sample_strings_;
    std::FILE* out_;
    ConvertOptions opts_;
    size_t n_blocks_ = 0;

public:
    uint64_t n_written = 0, n_skipped = 0;

    BlockWriter(std::FILE* out,
                std::unordered_map<std::string, uint32_t>& hd,
                std::unordered_map<std::string, uint32_t>& sd,
                std::vector<std::string>& hs,
                std::vector<std::string>& ss,
                ConvertOptions opts)
        : hog_dict_(hd), sample_dict_(sd), hog_strings_(hs), sample_strings_(ss),
          out_(out), opts_(opts) {}

    uint32_t intern_hog(const std::string& s) {
        auto [it, inserted] = hog_dict_.emplace(s, static_cast<uint32_t>(hog_strings_.size()));
        if (inserted) hog_strings_.push_back(s);
        return it->second;
    }

    uint32_t intern_sample(const std::string& s) {
        auto [it, inserted] = sample_dict_.emplace(s, static_cast<uint32_t>(sample_strings_.size()));
        if (inserted) sample_strings_.push_back(s);
        return it->second;
    }

    void add(AlignmentRecord&& r) {
        records_.push_back(std::move(r));
        if (records_.size() >= opts_.block_recs)
            flush();
    }

    void flush() {
        if (records_.empty()) return;
        raw_buf_.clear();
        for (auto& r : records_)
            serialize_record(raw_buf_, r);

        // build block header stats
        BlockHeader bh{};
        bh.block_type    = static_cast<uint8_t>(BlockType::Alignments);
        bh.n_records     = static_cast<uint32_t>(records_.size());
        bh.min_hog_idx   = UINT32_MAX;
        bh.max_hog_idx   = 0;
        bh.min_sstart    = UINT32_MAX;
        bh.max_send      = 0;
        bh.min_pident_q  = 255;
        bh.max_pident_q  = 0;
        bh.min_evalue_q  = 0;     // worst (closest to 0)
        bh.max_evalue_q  = -128;  // best (most negative)
        for (auto& r : records_) {
            bh.min_hog_idx  = std::min(bh.min_hog_idx,  r.hog_idx);
            bh.max_hog_idx  = std::max(bh.max_hog_idx,  r.hog_idx);
            bh.min_sstart   = std::min(bh.min_sstart,   r.sstart);
            bh.max_send     = std::max(bh.max_send,     r.send);
            bh.min_pident_q = std::min(bh.min_pident_q, r.pident_q);
            bh.max_pident_q = std::max(bh.max_pident_q, r.pident_q);
            bh.min_evalue_q = std::max(bh.min_evalue_q, r.evalue_q);  // worst = max
            bh.max_evalue_q = std::min(bh.max_evalue_q, r.evalue_q);  // best = min
        }
        records_.clear();

        // compress
        size_t bound = ZSTD_compressBound(raw_buf_.size());
        std::vector<uint8_t> cbuf(bound);
        size_t csz = ZSTD_compress(cbuf.data(), bound, raw_buf_.data(), raw_buf_.size(), opts_.zstd_level);
        if (ZSTD_isError(csz))
            throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(csz));

        bh.raw_sz        = static_cast<uint32_t>(raw_buf_.size());
        bh.compressed_sz = static_cast<uint32_t>(csz);

        std::fwrite(&bh, sizeof(bh), 1, out_);
        std::fwrite(cbuf.data(), 1, csz, out_);
        ++n_blocks_;
    }

    size_t n_blocks() const { return n_blocks_; }
};

// Write dict block (HOG or sample strings, newline-delimited, zstd-compressed)
inline void write_dict_block(std::FILE* out, BlockType type,
                              const std::vector<std::string>& strings, int zstd_level) {
    std::string raw;
    for (auto& s : strings) { raw += s; raw += '\n'; }
    size_t bound = ZSTD_compressBound(raw.size());
    std::vector<uint8_t> cbuf(bound);
    size_t csz = ZSTD_compress(cbuf.data(), bound, raw.data(), raw.size(), zstd_level);
    if (ZSTD_isError(csz))
        throw std::runtime_error(std::string("zstd dict: ") + ZSTD_getErrorName(csz));

    BlockHeader bh{};
    bh.block_type    = static_cast<uint8_t>(type);
    bh.n_records     = static_cast<uint32_t>(strings.size());
    bh.raw_sz        = static_cast<uint32_t>(raw.size());
    bh.compressed_sz = static_cast<uint32_t>(csz);
    std::fwrite(&bh, sizeof(bh), 1, out);
    std::fwrite(cbuf.data(), 1, csz, out);
}

inline void convert(const std::string& tsv_path, const std::string& out_path,
                    const std::string& sample_acc_override, ConvertOptions opts) {

    std::unordered_map<std::string, uint32_t> hog_dict, sample_dict;
    std::vector<std::string> hog_strings, sample_strings;

    // temp file for blocks (dicts go first in final file, blocks follow)
    std::string tmp_path = out_path + ".blocks.tmp";
    std::FILE* tmp = std::fopen(tmp_path.c_str(), "wb");
    if (!tmp) throw std::runtime_error("cannot open " + tmp_path);

    BlockWriter bw(tmp, hog_dict, sample_dict, hog_strings, sample_strings, opts);

    std::istream* in_ptr = &std::cin;
    std::ifstream fin;
    if (tsv_path != "-") {
        fin.open(tsv_path);
        if (!fin) throw std::runtime_error("cannot open " + tsv_path);
        in_ptr = &fin;
    }

    std::string line;
    uint64_t lineno = 0;
    while (std::getline(*in_ptr, line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;

        // split on tabs (14 fields expected)
        std::array<std::string_view, 14> f;
        size_t fi = 0;
        const char* p = line.data();
        const char* ep = p + line.size();
        const char* fp = p;
        while (p <= ep && fi < 14) {
            if (p == ep || *p == '\t') {
                f[fi++] = std::string_view(fp, p - fp);
                fp = p + 1;
            }
            ++p;
        }
        if (fi < 13) {
            if (opts.verbose)
                std::cerr << "skip line " << lineno << ": only " << fi << " fields\n";
            ++bw.n_skipped;
            continue;
        }

        // parse fields
        auto qseqid = f[0];
        auto sseqid = f[5];
        uint32_t sstart, send;
        float pident; double evalue;
        try {
            sstart = static_cast<uint32_t>(std::stoul(std::string(f[6])));
            send   = static_cast<uint32_t>(std::stoul(std::string(f[7])));
            pident = std::stof(std::string(f[9]));
        } catch (const std::exception& e) {
            if (opts.verbose)
                std::cerr << "skip line " << lineno << ": " << e.what() << "\n";
            ++bw.n_skipped; continue;
        }
        // stod throws out_of_range for subnormal doubles (e.g. 8e-316 < DBL_MIN)
        // treat as 0.0 → quantizes to -128 (best evalue bucket)
        try { evalue = std::stod(std::string(f[10])); }
        catch (const std::out_of_range&) { evalue = 0.0; }
        auto     cigar   = f[11];
        auto     qseq    = f[12];

        // filters
        if (pident < opts.min_pident) { ++bw.n_skipped; continue; }
        int8_t eq = quantize_evalue(evalue);
        if (eq > opts.max_evalue_q) { ++bw.n_skipped; continue; }

        if (sstart > send) { ++bw.n_skipped; continue; }

        std::string hog    = extract_hog(sseqid);
        std::string sample = sample_acc_override.empty()
                           ? extract_sample(qseqid)
                           : sample_acc_override;

        AlignmentRecord r;
        r.hog_idx    = bw.intern_hog(hog);
        r.sample_idx = bw.intern_sample(sample);
        r.sstart     = sstart;
        r.send       = send;
        r.pident_q   = quantize_pident(pident);
        r.evalue_q   = eq;
        r.qframe     = extract_qframe(qseqid);
        r.aas        = cigar_to_aas(cigar, qseq, sstart, send);

        ++bw.n_written;
        bw.add(std::move(r));
    }
    bw.flush();

    // write dicts + blocks to final file
    std::FILE* fout = std::fopen(out_path.c_str(), "wb");
    if (!fout) { std::fclose(tmp); throw std::runtime_error("cannot open " + out_path); }

    FileHeader fh{};
    std::memcpy(fh.magic, MAGIC, 8);
    fh.version  = FORMAT_VERSION;
    fh.n_blocks = static_cast<uint32_t>(2 + bw.n_blocks());  // 2 dict blocks + N aln blocks
    std::fwrite(&fh, sizeof(fh), 1, fout);

    write_dict_block(fout, BlockType::HOGDict,    hog_strings,    opts.zstd_level);
    write_dict_block(fout, BlockType::SampleDict, sample_strings, opts.zstd_level);

    // append block data from tmp
    std::fclose(tmp);
    std::FILE* ttmp = std::fopen(tmp_path.c_str(), "rb");
    if (ttmp) {
        std::vector<uint8_t> buf(1 << 20);
        size_t n;
        while ((n = std::fread(buf.data(), 1, buf.size(), ttmp)) > 0)
            std::fwrite(buf.data(), 1, n, fout);
        std::fclose(ttmp);
    }
    std::fclose(fout);
    std::filesystem::remove(tmp_path);

    std::cerr << "wrote " << bw.n_written << " records, "
              << bw.n_skipped << " skipped, "
              << bw.n_blocks() << " blocks → " << out_path << "\n";
}

} // namespace lhi
