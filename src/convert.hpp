#pragma once
#include "format.hpp"
#include "cigar.hpp"
#include "aa.hpp"
#include "batch.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <charconv>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <map>
#include <sys/stat.h>
#include <zlib.h>

namespace lhi {

struct ConvertOptions {
    std::string acc_id;          // SRA/ENA accession; "auto" = extract from qseqid prefix
    int    zstd_level = 3;
    float  min_pident = 0.0f;
    double max_evalue = 1.0;
    bool   verbose    = false;
    size_t flush_bytes = 0; // 0 = auto-detect from cgroup; else explicit threshold in bytes
    size_t buckets    = 0;  // spill bucket count B (0 = auto from input size + flush budget)
};

// diamond blastx outfmt 6 columns (0-based):
//   qseqid qstart qend qlen qstrand sseqid sstart send slen pident evalue cigar qseq_translated full_qseq
namespace col {
    constexpr int qseqid=0, qstart=1, qend=2, qlen=3, qstrand=4;
    constexpr int sseqid=5, sstart=6, send=7,  slen=8;
    constexpr int pident=9, evalue=10, cigar=11, qseq_aa=12, full_qseq=13;
}

// Returns a view into the caller's line buffer (alive for the call duration).
inline std::string_view extract_hog(std::string_view sv) {
    auto p = sv.rfind('|');
    return p != std::string_view::npos ? sv.substr(p + 1) : sv;
}

// Contig unitig-id → numeric cnum: the integer field after the first '_' (e.g.
// "chr2H_18749_..." → 18749). Returns `fallback` when no parseable field is present.
// Shared by the merge-shard spill encoder (spill_extent_into) and the fused ingest
// flush so the two paths produce identical cnum for parseable contigs.
inline uint32_t parse_cnum(std::string_view uid, uint32_t fallback) {
    auto fs = uid.find('_');
    if (fs != std::string_view::npos && fs + 1 < uid.size()) {
        auto fe = uid.find('_', fs + 1);
        auto part = (fe != std::string_view::npos) ? uid.substr(fs + 1, fe - fs - 1)
                                                   : uid.substr(fs + 1);
        uint32_t v = 0;
        auto cr = std::from_chars(part.data(), part.data() + part.size(), v);
        if (cr.ec == std::errc{}) return v;
    }
    return fallback;
}

// Transparent hash/equality so unordered_map<std::string,...> can be looked up
// by std::string_view without constructing a temporary std::string.
struct SvHash {
    using is_transparent = void;
    size_t operator()(std::string_view sv) const { return std::hash<std::string_view>{}(sv); }
    size_t operator()(const std::string& s)  const { return std::hash<std::string_view>{}(s); }
};

inline int8_t make_qframe(std::string_view qstrand, uint32_t qstart,
                          uint32_t qend, uint32_t qlen) {
    bool minus = (!qstrand.empty() && qstrand[0] == '-');
    uint32_t base = minus ? (qlen > qend ? qlen - qend : 0)
                          : (qstart > 0  ? qstart - 1  : 0);
    int frame = int(base % 3) + 1;
    return int8_t(minus ? -frame : frame);
}

inline void revcomp_codon(uint8_t c[3]) {
    auto rc = [](uint8_t b) -> uint8_t {
        switch (b) { case 'A':return 'T'; case 'T':return 'A'; case 'C':return 'G'; case 'G':return 'C'; default:return 'N'; }
    };
    uint8_t t[3] = {rc(c[2]), rc(c[1]), rc(c[0])};
    c[0]=t[0]; c[1]=t[1]; c[2]=t[2];
}

using SvDict = std::unordered_map<std::string, uint32_t, SvHash, std::equal_to<>>;

// Parallel multi-frame zstd decoder. A zstd stream made of >1 independent frames can be
// decoded concurrently (one frame per thread). One reader thread streams the compressed
// file once, locates frame boundaries by walking block headers (no decompression), and
// hands each frame's compressed bytes to a decoder pool; decoded frames are delivered to
// the consumer IN ORDER, so a line spanning a frame boundary is reassembled by the caller's
// getline. Memory is bounded by compressed- and decompressed-inflight byte caps.
// Single-frame .zst can't be split → TsvReader uses its serial path instead.
struct ParZstd {
    struct Job { size_t idx; std::vector<char> comp; };
    int fd;
    size_t nthreads, cap_comp, cap_dec;
    std::mutex m;
    std::condition_variable cv_job, cv_ready, cv_room;
    std::deque<Job> jobq;
    std::map<size_t, std::vector<char>> ready;   // idx -> decoded
    size_t comp_inflight = 0, dec_inflight = 0, next_consume = 0, decoding = 0;
    bool reader_done = false, stop_ = false, err = false;
    std::string emsg;
    std::thread reader;
    std::vector<std::thread> decs;

    ParZstd(int fd_, size_t nt, size_t cap_comp_, size_t cap_dec_)
        : fd(fd_), nthreads(std::max<size_t>(1, nt)),
          cap_comp(std::max<size_t>(cap_comp_, size_t(64) << 20)),
          cap_dec(std::max<size_t>(cap_dec_, size_t(64) << 20)) {
        reader = std::thread([this] { read_loop(); });
        for (size_t i = 0; i < nthreads; ++i) decs.emplace_back([this] { decode_loop(); });
    }
    ~ParZstd() {
        { std::lock_guard<std::mutex> l(m); stop_ = true; }
        cv_job.notify_all(); cv_ready.notify_all(); cv_room.notify_all();
        if (reader.joinable()) reader.join();
        for (auto& t : decs) if (t.joinable()) t.join();
    }
    ParZstd(const ParZstd&) = delete; ParZstd& operator=(const ParZstd&) = delete;

    void set_err(const std::string& e) {
        std::lock_guard<std::mutex> l(m);
        if (!err) { err = true; emsg = e; }
        cv_job.notify_all(); cv_ready.notify_all(); cv_room.notify_all();
    }

    static ssize_t pread_all(int fd, void* buf, size_t n, uint64_t off) {
        char* p = static_cast<char*>(buf); size_t done = 0;
        while (done < n) { ssize_t r = ::pread(fd, p + done, n - done, off_t(off + done));
                           if (r < 0) return r; if (r == 0) break; done += size_t(r); }
        return ssize_t(done);
    }

    static constexpr uint64_t CAP_EXCEEDED = ~uint64_t(0);
    // Length of the zstd frame at `off` (compressed bytes), via block-header walk. 0 = no frame.
    // If cap>0 and the frame exceeds `cap` bytes, returns CAP_EXCEEDED early (used by the peek).
    static uint64_t frame_len(int fd, uint64_t off, uint64_t fsize, uint64_t cap = 0) {
        uint8_t h[18];
        if (pread_all(fd, h, 4, off) != 4) return 0;
        // Skippable frame: magic 0x184D2A50..5F → size in next 4 bytes.
        if (h[3] == 0x18 && h[2] == 0x4d && (h[1] & 0xf0) == 0x20 && (h[0] & 0xf0) == 0x50) {
            uint8_t s[4]; if (pread_all(fd, s, 4, off + 4) != 4) return 0;
            return 8 + (uint64_t(s[0]) | s[1] << 8 | s[2] << 16 | uint64_t(s[3]) << 24);
        }
        if (!(h[0] == 0x28 && h[1] == 0xb5 && h[2] == 0x2f && h[3] == 0xfd))
            throw std::runtime_error("zstd: bad frame magic");
        ssize_t hn = pread_all(fd, h, 18, off); if (hn < 6) throw std::runtime_error("zstd: short frame header");
        uint8_t fhd = h[4];
        int fcs = fhd >> 6, ss = (fhd >> 5) & 1, cc = (fhd >> 2) & 1, did = fhd & 3;
        uint64_t p = off + 5;
        if (!ss) p += 1;
        p += (did == 0 ? 0 : did == 1 ? 1 : did == 2 ? 2 : 4);
        p += (fcs == 0 ? (ss ? 1 : 0) : fcs == 1 ? 2 : fcs == 2 ? 4 : 8);
        for (;;) {
            uint8_t b[3]; if (pread_all(fd, b, 3, p) != 3) throw std::runtime_error("zstd: truncated block header");
            uint32_t v = uint32_t(b[0]) | b[1] << 8 | b[2] << 16;
            uint32_t last = v & 1, bt = (v >> 1) & 3, bsz = v >> 3;
            p += 3 + (bt == 1 ? 1 : bsz);
            if (p > fsize) throw std::runtime_error("zstd: block overruns file");
            if (cap && p - off > cap) return CAP_EXCEEDED;
            if (last) break;
        }
        if (cc) p += 4;
        return p - off;
    }

    void read_loop() {
        try {
            struct stat st{}; if (fstat(fd, &st) != 0) throw std::runtime_error("zstd: fstat");
            uint64_t fsize = uint64_t(st.st_size), off = 0; size_t idx = 0;
            while (off < fsize) {
                uint64_t flen = frame_len(fd, off, fsize);
                if (flen == 0) break;
                std::vector<char> comp(static_cast<size_t>(flen));
                if (pread_all(fd, comp.data(), size_t(flen), off) != ssize_t(flen))
                    throw std::runtime_error("zstd: short frame read");
                off += flen;
                std::unique_lock<std::mutex> l(m);
                cv_room.wait(l, [&] { return stop_ || err || comp_inflight + comp.size() <= cap_comp || jobq.empty(); });
                if (stop_ || err) return;
                comp_inflight += comp.size();
                jobq.push_back({idx++, std::move(comp)});
                cv_job.notify_one();
            }
            std::lock_guard<std::mutex> l(m); reader_done = true; cv_job.notify_all(); cv_ready.notify_all();
        } catch (const std::exception& e) { set_err(e.what()); }
    }

    void decode_loop() {
        for (;;) {
            Job job;
            {
                std::unique_lock<std::mutex> l(m);
                for (;;) {                       // re-check everything after each wake
                    if (stop_ || err) return;
                    if (!jobq.empty()) {
                        // admit if decompressed inflight is under budget, OR this is the exact
                        // frame the consumer is blocked on (deadlock guard; jobq is FIFO by idx).
                        if (dec_inflight <= cap_dec || jobq.front().idx == next_consume) {
                            job = std::move(jobq.front()); jobq.pop_front();
                            comp_inflight -= job.comp.size();
                            ++decoding;
                            cv_room.notify_all();
                            break;
                        }
                        cv_room.wait(l);         // inflight full, front not urgent → wait for room
                        continue;
                    }
                    if (reader_done) return;
                    cv_job.wait(l);
                }
            }
            try {
                unsigned long long usz = ZSTD_getFrameContentSize(job.comp.data(), job.comp.size());
                std::vector<char> out;
                if (usz != ZSTD_CONTENTSIZE_UNKNOWN && usz != ZSTD_CONTENTSIZE_ERROR) {
                    out.resize(size_t(usz));
                    size_t r = ZSTD_decompress(out.data(), out.size(), job.comp.data(), job.comp.size());
                    if (ZSTD_isError(r) || r != usz) throw std::runtime_error("zstd: frame decompress");
                } else {
                    // unknown size: stream-grow
                    ZSTD_DStream* ds = ZSTD_createDStream(); ZSTD_initDStream(ds);
                    ZSTD_inBuffer in{job.comp.data(), job.comp.size(), 0};
                    std::vector<char> buf(size_t(4) << 20);
                    for (;;) {
                        ZSTD_outBuffer ob{buf.data(), buf.size(), 0};
                        size_t rc = ZSTD_decompressStream(ds, &ob, &in);
                        if (ZSTD_isError(rc)) { ZSTD_freeDStream(ds); throw std::runtime_error("zstd: stream frame"); }
                        out.insert(out.end(), buf.data(), buf.data() + ob.pos);
                        if (in.pos >= in.size) break;
                    }
                    ZSTD_freeDStream(ds);
                }
                std::lock_guard<std::mutex> l(m);
                dec_inflight += out.size();
                --decoding;
                ready.emplace(job.idx, std::move(out));
                cv_ready.notify_all();
            } catch (const std::exception& e) { set_err(e.what()); return; }
        }
    }

    // Pull the next decoded frame in order. Returns false at clean end.
    bool next(std::vector<char>& out) {
        std::unique_lock<std::mutex> l(m);
        for (;;) {
            if (err) throw std::runtime_error("parallel zstd: " + emsg);
            auto it = ready.find(next_consume);
            if (it != ready.end()) {
                out = std::move(it->second); ready.erase(it);
                dec_inflight -= out.size(); ++next_consume;
                cv_room.notify_all();
                return true;
            }
            // Clean end: reader finished, no queued jobs, none being decoded → frame can't appear.
            if (reader_done && jobq.empty() && decoding == 0) return false;
            cv_ready.wait(l);
        }
    }
};

// Transparent line reader: plain file, .gz file, .zst file, or stdin ("-").
// Provides getline(std::string&) → bool.
struct TsvReader {
    enum class Kind { Plain, Gz, Zst, Stdin };
    Kind kind;
    std::ifstream fin;
    gzFile gz = nullptr;
    std::vector<char> gz_buf;

    // zstd streaming state
    int          zst_fd  = -1;
    ZSTD_DStream* zst_dctx = nullptr;
    std::vector<char> zst_in;   // compressed input buffer
    std::vector<char> zst_out;  // decompressed output buffer (also holds a parallel frame)
    size_t zst_in_pos  = 0;
    size_t zst_in_end  = 0;
    size_t zst_out_pos = 0;
    size_t zst_out_end = 0;
    bool   zst_eof     = false; // no more compressed bytes
    std::unique_ptr<ParZstd> par;   // set iff multi-frame .zst → parallel decode
    bool   zst_par = false;

    explicit TsvReader(const std::string& path, size_t threads = 1) {
        if (path == "-") {
            kind = Kind::Stdin;
        } else if (path.size() > 4 && path.compare(path.size()-4, 4, ".zst") == 0) {
            kind = Kind::Zst;
            zst_fd = ::open(path.c_str(), O_RDONLY);
            if (zst_fd < 0) throw std::runtime_error("cannot open: " + path);
            posix_fadvise(zst_fd, 0, 0, POSIX_FADV_SEQUENTIAL);
            // Peek frame 0 (cap the walk so a huge single frame bails fast). If it ends well
            // before EOF, the stream is multi-frame → decode frames in parallel. Otherwise
            // (single/huge frame) → serial streaming (can't split one frame).
            struct stat st{}; uint64_t fsize = (fstat(zst_fd, &st) == 0) ? uint64_t(st.st_size) : 0;
            uint64_t f0 = 0;
            try { f0 = ParZstd::frame_len(zst_fd, 0, fsize, uint64_t(256) << 20); } catch (...) { f0 = 0; }
            if (threads > 1 && f0 != 0 && f0 != ParZstd::CAP_EXCEEDED && f0 < fsize) {
                zst_par = true;
                size_t cap_comp = std::max<size_t>(size_t(256) << 20, threads * (size_t(32) << 20));
                size_t cap_dec  = std::max<size_t>(size_t(512) << 20, threads * (size_t(64) << 20));
                par = std::make_unique<ParZstd>(zst_fd, threads, cap_comp, cap_dec);
            } else {
                zst_dctx = ZSTD_createDStream();
                if (!zst_dctx) throw std::runtime_error("ZSTD_createDStream failed");
                zst_in.resize(ZSTD_DStreamInSize());
                zst_out.resize(ZSTD_DStreamOutSize() * 4); // 4× recommended for fewer refills
            }
        } else if (path.size() > 3 && path.compare(path.size()-3, 3, ".gz") == 0) {
            kind = Kind::Gz;
            gz = gzopen(path.c_str(), "r");
            if (!gz) throw std::runtime_error("cannot open: " + path);
            gzbuffer(gz, 1u << 20);
            gz_buf.resize(1u << 16);
        } else {
            kind = Kind::Plain;
            fin.open(path);
            if (!fin) throw std::runtime_error("cannot open: " + path);
        }
    }
    ~TsvReader() {
        par.reset();   // join decode threads (they pread zst_fd) BEFORE closing it
        if (gz) gzclose(gz);
        if (zst_dctx) ZSTD_freeDStream(zst_dctx);
        if (zst_fd >= 0) ::close(zst_fd);
    }
    TsvReader(const TsvReader&) = delete;
    TsvReader& operator=(const TsvReader&) = delete;

    // Decompress more bytes into zst_out[0..zst_out_end). Returns false only at real EOF.
    // Loops over frame boundaries so a multi-frame stream is read in full (a frame ending with
    // no new output and input still available must NOT be mistaken for EOF — that was a bug that
    // silently truncated multi-frame .zst at the first frame).
    bool zst_refill() {
        for (;;) {
            if (zst_eof && zst_in_pos >= zst_in_end) return false;
            if (zst_in_pos >= zst_in_end && !zst_eof) {
                ssize_t r = ::read(zst_fd, zst_in.data(), zst_in.size());
                if (r <= 0) { zst_eof = true; return false; }
                zst_in_pos = 0; zst_in_end = size_t(r);
            }
            ZSTD_inBuffer  in  = {zst_in.data() + zst_in_pos, zst_in_end - zst_in_pos, 0};
            ZSTD_outBuffer out = {zst_out.data(), zst_out.size(), 0};
            size_t ret = ZSTD_decompressStream(zst_dctx, &out, &in);
            if (ZSTD_isError(ret))
                throw std::runtime_error(std::string("zstd decompress: ") + ZSTD_getErrorName(ret));
            zst_in_pos  += in.pos;
            zst_out_pos  = 0;
            zst_out_end  = out.pos;
            if (zst_out_end > 0) return true;
            // No output this round (frame boundary / needs more input) — loop and continue.
        }
    }

    // Get the next decompressed window into zst_out[0..zst_out_end). Parallel path pulls the next
    // whole frame (in order); serial path streams. Returns false at end.
    bool zst_next_window() {
        if (zst_par) {
            std::vector<char> b;
            if (!par->next(b)) return false;
            zst_out = std::move(b);
            zst_out_pos = 0; zst_out_end = zst_out.size();
            return true;
        }
        return zst_refill();
    }

    bool getline(std::string& line) {
        if (kind == Kind::Zst) {
            line.clear();
            while (true) {
                // Scan current output window for '\n'
                if (zst_out_pos < zst_out_end) {
                    const char* base = zst_out.data() + zst_out_pos;
                    size_t avail = zst_out_end - zst_out_pos;
                    const char* nl = static_cast<const char*>(memchr(base, '\n', avail));
                    if (nl) {
                        size_t len = size_t(nl - base);
                        line.append(base, len);
                        if (!line.empty() && line.back() == '\r') line.pop_back();
                        zst_out_pos += len + 1;
                        return true;
                    }
                    line.append(base, avail);
                    zst_out_pos = zst_out_end = 0;
                }
                // Need more decompressed data
                if (!zst_next_window()) return !line.empty();
            }
        }
        if (kind == Kind::Gz) {
            line.clear();
            while (true) {
                if (!gzgets(gz, gz_buf.data(), int(gz_buf.size()))) return !line.empty();
                line += gz_buf.data();
                if (!line.empty() && line.back() == '\n') {
                    line.pop_back();
                    if (!line.empty() && line.back() == '\r') line.pop_back();
                    return true;
                }
            }
        }
        std::istream& s = (kind == Kind::Stdin) ? std::cin : static_cast<std::istream&>(fin);
        return bool(std::getline(s, line));
    }
};

// Convert a diamond blastx TSV (outfmt 6 + qseq_translated + full_qseq columns)
// to a .lhb batch file.  All HOGs go into a single output file; no per-HOG shards.
inline void convert(const std::string& tsv_path, const std::string& lhb_path,
                    ConvertOptions opts) {
    if (opts.acc_id.empty())
        throw std::runtime_error("acc_id is required for .lhb output (use -a ACC)");

    BatchWriter batch(lhb_path, opts.acc_id);

    SvDict hog_dict, contig_dict;
    std::vector<std::string> hog_strings, contig_strings;

    // Lookup by string_view (no allocation on hit); allocate only on first sight.
    auto intern = [](SvDict& d, std::vector<std::string>& v, std::string_view s) -> uint32_t {
        auto it = d.find(s);
        if (it != d.end()) return it->second;
        uint32_t idx = uint32_t(v.size());
        v.emplace_back(s);
        d.emplace(v.back(), idx);
        return idx;
    };

    std::unordered_map<uint32_t, std::vector<VarNTRecord>> batches;

    auto flush_hog = [&](uint32_t hog_idx) {
        auto it = batches.find(hog_idx);
        if (it == batches.end() || it->second.empty()) return;
        batch.write_block(hog_strings[hog_idx], contig_strings, it->second, opts.zstd_level);
        it->second.clear();
    };

    TsvReader reader(tsv_path);
    std::string line;
    uint64_t lineno = 0, n_written = 0, n_skipped = 0, n_obs_dropped = 0;

    while (reader.getline(line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;

        std::array<std::string_view, 14> f{};
        size_t fi = 0;
        const char* p = line.data(), *ep = p + line.size(), *fp = p;
        while (p <= ep && fi < 14) {
            if (p == ep || *p == '\t') { f[fi++] = {fp, size_t(p-fp)}; fp = p+1; }
            ++p;
        }
        if (fi < 13) { if (opts.verbose) std::cerr << "skip L" << lineno << ": " << fi << " fields\n"; ++n_skipped; continue; }
        bool has_nt = (fi == 14 && !f[col::full_qseq].empty());

        float  pident = 0.0f;
        double ev     = 0.0;
        std::from_chars(f[col::pident].data(), f[col::pident].data() + f[col::pident].size(), pident);
        std::from_chars(f[col::evalue].data(),  f[col::evalue].data()  + f[col::evalue].size(),  ev);

        if (pident < opts.min_pident) { ++n_skipped; continue; }
        if (ev     > opts.max_evalue) { ++n_skipped; continue; }
        if (pident == 100.0f)         { ++n_skipped; continue; }

        uint32_t sstart = 0, send = 0, qstart = 0, qend = 0, qlen = 0;
        {
            auto fc = [](std::string_view sv, uint32_t& out) -> bool {
                auto [p, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), out);
                return ec == std::errc{};
            };
            if (!fc(f[col::sstart], sstart) || !fc(f[col::send], send) ||
                !fc(f[col::qstart], qstart) || !fc(f[col::qend], qend) ||
                !fc(f[col::qlen],   qlen)) {
                if (opts.verbose) std::cerr << "parse error L" << lineno << "\n";
                ++n_skipped; continue;
            }
        }
        if (sstart > send) { ++n_skipped; continue; }

        uint32_t    hog_idx    = intern(hog_dict, hog_strings, extract_hog(f[col::sseqid]));
        uint32_t    contig_idx = intern(contig_dict, contig_strings, f[col::qseqid]);
        int8_t      qframe     = make_qframe(f[col::qstrand], qstart, qend, qlen);

        auto ar = cigar_parse(f[col::cigar], f[col::qseq_aa], sstart, send);

        VarNTRecord vr;
        vr.contig_idx = contig_idx;
        vr.sstart     = sstart; vr.send = send;
        vr.qframe     = qframe;
        vr.pident     = pident; vr.evalue = ev;

        std::string_view full_nt = has_nt ? f[col::full_qseq] : std::string_view{};
        uint32_t span = send - sstart + 1;

        for (uint32_t i = 0; i < span; ++i) {
            uint8_t obs_aa = ar.aas[i];
            if (obs_aa == AA_GAP || obs_aa == AA_UNK) { ++n_obs_dropped; continue; }
            uint32_t q_off = ar.qseq_offsets[i];
            if (q_off == UINT32_MAX || full_nt.empty()) { ++n_obs_dropped; continue; }

            uint8_t c0, c1, c2;
            if (qframe > 0) {
                size_t cs = size_t(qstart-1) + size_t(q_off)*3;
                if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
            } else {
                // Minus strand: diamond reports qstart > qend; codon i sits at
                // [(qstart-1) - q_off*3 - 2 .. (qstart-1) - q_off*3], revcomp'd.
                if (size_t(q_off)*3 + 3 > size_t(qstart)) { ++n_obs_dropped; continue; }
                size_t cs = size_t(qstart-1) - size_t(q_off)*3 - 2;
                if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
                uint8_t tmp[3]={c0,c1,c2}; revcomp_codon(tmp); c0=tmp[0]; c1=tmp[1]; c2=tmp[2];
            }

            // Reject ambiguous NT (N, lowercase, non-ACGT) — pack_codon would coerce to T.
            auto is_acgt = [](uint8_t b) {
                return b=='A'||b=='C'||b=='G'||b=='T';
            };
            if (!is_acgt(c0) || !is_acgt(c1) || !is_acgt(c2)) { ++n_obs_dropped; continue; }

            uint8_t raw3[3]={c0,c1,c2};
            uint8_t packed = pack_codon(raw3);

            // Round-trip: codon→AA must agree with diamond's translated AA.
            // Disagreement indicates a strand/frame extraction error; discard the observation.
            if (codon_to_aa(packed) != obs_aa) {
                if (opts.verbose) std::cerr << "codon/AA mismatch L" << lineno
                    << " (diamond=" << char(AA_ALPHA[obs_aa])
                    << " nt=" << char(c0) << char(c1) << char(c2) << ")\n";
                ++n_obs_dropped; continue;
            }

            vr.vars.push_back({i, packed});
        }

        if (!vr.vars.empty())
            batches[hog_idx].push_back(std::move(vr));
        ++n_written;
    }

    // Flush in HOG-ID order so .lhb blocks are lexicographically sorted, which
    // lets hoki merge do a streaming k-way merge.
    std::vector<uint32_t> sorted_idxs;
    sorted_idxs.reserve(batches.size());
    for (auto& [hog_idx, recs] : batches)
        if (!recs.empty()) sorted_idxs.push_back(hog_idx);
    std::sort(sorted_idxs.begin(), sorted_idxs.end(),
        [&](uint32_t a, uint32_t b) { return hog_strings[a] < hog_strings[b]; });
    for (auto hog_idx : sorted_idxs) flush_hog(hog_idx);
    batch.finalize();

    std::cerr << "convert: " << n_written << " records, " << n_skipped
              << " records skipped, " << n_obs_dropped << " observations dropped → "
              << lhb_path << "\n";
}

// Multi-accession convert: acc_id detected from qseqid prefix (before first '_').
// out_dir is created if it doesn't exist; writes out_dir/ACC.lhb per accession.
inline void convert_multi(const std::string& tsv_path, const std::string& out_dir,
                          ConvertOptions opts) {
    if (mkdir(out_dir.c_str(), 0755) != 0 && errno != EEXIST)
        throw std::runtime_error("cannot create output dir: " + out_dir);

    struct AccState {
        std::unique_ptr<BatchWriter> writer;
        SvDict contig_dict;
        std::vector<std::string> contig_strings;
        std::unordered_map<uint32_t, std::vector<VarNTRecord>> batches;
    };
    std::unordered_map<std::string, AccState> acc_states;

    SvDict hog_dict;
    std::vector<std::string> hog_strings;

    TsvReader reader(tsv_path);
    std::string line;
    uint64_t lineno = 0, n_written = 0, n_skipped = 0, n_obs_dropped = 0;

    auto intern = [](SvDict& d, std::vector<std::string>& v, std::string_view s) -> uint32_t {
        auto it = d.find(s);
        if (it != d.end()) return it->second;
        uint32_t idx = uint32_t(v.size());
        v.emplace_back(s);
        d.emplace(v.back(), idx);
        return idx;
    };

    while (reader.getline(line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;

        std::array<std::string_view, 14> f{};
        size_t fi = 0;
        const char* p = line.data(), *ep = p + line.size(), *fp = p;
        while (p <= ep && fi < 14) {
            if (p == ep || *p == '\t') { f[fi++] = {fp, size_t(p-fp)}; fp = p+1; }
            ++p;
        }
        if (fi < 13) { if (opts.verbose) std::cerr << "skip L" << lineno << ": " << fi << " fields\n"; ++n_skipped; continue; }
        bool has_nt = (fi == 14 && !f[col::full_qseq].empty());

        float  pident = 0.0f;
        double ev     = 0.0;
        std::from_chars(f[col::pident].data(), f[col::pident].data() + f[col::pident].size(), pident);
        std::from_chars(f[col::evalue].data(),  f[col::evalue].data()  + f[col::evalue].size(),  ev);

        if (pident < opts.min_pident) { ++n_skipped; continue; }
        if (ev     > opts.max_evalue) { ++n_skipped; continue; }
        if (pident == 100.0f)         { ++n_skipped; continue; }

        uint32_t sstart = 0, send = 0, qstart = 0, qend = 0, qlen = 0;
        {
            auto fc = [](std::string_view sv, uint32_t& out) -> bool {
                auto [p, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), out);
                return ec == std::errc{};
            };
            if (!fc(f[col::sstart], sstart) || !fc(f[col::send], send) ||
                !fc(f[col::qstart], qstart) || !fc(f[col::qend], qend) ||
                !fc(f[col::qlen],   qlen)) {
                if (opts.verbose) std::cerr << "parse error L" << lineno << "\n";
                ++n_skipped; continue;
            }
        }
        if (sstart > send) { ++n_skipped; continue; }

        // Extract accession from qseqid prefix before first '_'
        std::string_view qseqid = f[col::qseqid];
        auto us = qseqid.find('_');
        std::string acc(us != std::string_view::npos ? qseqid.substr(0, us) : qseqid);

        AccState& st = acc_states[acc];
        if (!st.writer) {
            std::string path = out_dir + "/" + acc + ".lhb";
            st.writer = std::make_unique<BatchWriter>(path, acc);
        }

        uint32_t hog_idx    = intern(hog_dict, hog_strings, extract_hog(f[col::sseqid]));
        uint32_t contig_idx = intern(st.contig_dict, st.contig_strings, qseqid);
        int8_t   qframe     = make_qframe(f[col::qstrand], qstart, qend, qlen);
        auto     ar         = cigar_parse(f[col::cigar], f[col::qseq_aa], sstart, send);

        VarNTRecord vr;
        vr.contig_idx = contig_idx;
        vr.sstart     = sstart; vr.send = send;
        vr.qframe     = qframe;
        vr.pident     = pident; vr.evalue = ev;

        std::string_view full_nt = has_nt ? f[col::full_qseq] : std::string_view{};
        uint32_t span = send - sstart + 1;

        for (uint32_t i = 0; i < span; ++i) {
            uint8_t obs_aa = ar.aas[i];
            if (obs_aa == AA_GAP || obs_aa == AA_UNK) { ++n_obs_dropped; continue; }
            uint32_t q_off = ar.qseq_offsets[i];
            if (q_off == UINT32_MAX || full_nt.empty()) { ++n_obs_dropped; continue; }

            uint8_t c0, c1, c2;
            if (qframe > 0) {
                size_t cs = size_t(qstart-1) + size_t(q_off)*3;
                if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
            } else {
                if (size_t(q_off)*3 + 3 > size_t(qstart)) { ++n_obs_dropped; continue; }
                size_t cs = size_t(qstart-1) - size_t(q_off)*3 - 2;
                if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
                uint8_t tmp[3]={c0,c1,c2}; revcomp_codon(tmp); c0=tmp[0]; c1=tmp[1]; c2=tmp[2];
            }

            auto is_acgt = [](uint8_t b) { return b=='A'||b=='C'||b=='G'||b=='T'; };
            if (!is_acgt(c0) || !is_acgt(c1) || !is_acgt(c2)) { ++n_obs_dropped; continue; }

            uint8_t raw3[3]={c0,c1,c2};
            uint8_t packed = pack_codon(raw3);
            if (codon_to_aa(packed) != obs_aa) {
                if (opts.verbose) std::cerr << "codon/AA mismatch L" << lineno << "\n";
                ++n_obs_dropped; continue;
            }
            vr.vars.push_back({i, packed});
        }

        if (!vr.vars.empty())
            st.batches[hog_idx].push_back(std::move(vr));
        ++n_written;
    }

    // Flush each accession's batches in HOG-ID order
    for (auto& [acc, st] : acc_states) {
        std::vector<uint32_t> sorted_idxs;
        sorted_idxs.reserve(st.batches.size());
        for (auto& [hog_idx, recs] : st.batches)
            if (!recs.empty()) sorted_idxs.push_back(hog_idx);
        std::sort(sorted_idxs.begin(), sorted_idxs.end(),
            [&](uint32_t a, uint32_t b) { return hog_strings[a] < hog_strings[b]; });
        for (auto hog_idx : sorted_idxs) {
            auto it = st.batches.find(hog_idx);
            if (it != st.batches.end() && !it->second.empty())
                st.writer->write_block(hog_strings[hog_idx], st.contig_strings, it->second, opts.zstd_level);
        }
        st.writer->finalize();
    }

    std::cerr << "convert: " << n_written << " records, " << n_skipped
              << " records skipped, " << n_obs_dropped << " observations dropped, "
              << acc_states.size() << " accessions → " << out_dir << "/\n";
}

// Dispatcher: routes to convert_multi when acc_id == "auto", else single-acc convert.
inline void convert_dispatch(const std::string& tsv_path, const std::string& out_path,
                             ConvertOptions opts) {
    if (opts.acc_id == "auto")
        convert_multi(tsv_path, out_path, opts);
    else
        convert(tsv_path, out_path, opts);
}

} // namespace lhi
