#pragma once
#include "global.hpp"
#include "batch.hpp"
#include "container.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <cstring>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <stdexcept>
#include <thread>
#include <atomic>
#include <mutex>
#include <filesystem>

// Phase 2 — partition: scatter every ShardBlock from N .lhb inputs into per-thread
// sequential files tN.lhp + a partition.idx index. Each .lhb is read exactly once.
// Memory is O(n_threads × 8 MB write buffer + in-memory sidecar index).
//
// Design: each worker holds exactly ONE open fd for its tN.lhp file.  All HOGs
// interleave into one sequential stream keyed by (hog_id_len u16 | hog_id | PartitionEntry | payload).
// A sidecar index records (entry_offset, entry_len) per HOG so merge-shard can
// pread individual HOG extents without any per-HOG file opens.

namespace lhi {

constexpr uint8_t LHP_FILE_MAGIC[4]   = {'L','H','G','P'};
constexpr uint8_t LHP_VERSION         = 4;  // v4: column-major frames [other|bitmap|codon]
constexpr size_t  LHP_HEADER_SZ       = 8;  // magic(4)+ver(1)+flags(1)+pad(2)

constexpr uint8_t LHP_INDEX_MAGIC[4]  = {'L','H','P','I'};
constexpr uint8_t LHP_INDEX_VERSION   = 4;  // v4 matches LHP_VERSION 4

// v4: column-major frame layout. Each decompressed frame is:
//   [other_sec_len(u32) | bmp_sec_len(u32) | other_section | bitmap_section | codon_section]
// Extent fields point into these sections (offsets relative to section start).
struct PartitionIndexExtent {
    uint32_t thread_idx;
    uint64_t frame_off;    // byte offset of compressed frame in tN.lhp (after 4-byte csz prefix)
    uint32_t frame_csz;    // compressed frame size
    uint32_t other_off;    // offset in decompressed frame's other-section
    uint32_t other_len;    // length in other-section (contig_dict + varnt_hdr + cols except bmp/cdn)
    uint32_t bmp_off;      // offset in frame's bitmap-section
    uint32_t cdn_off;      // offset in frame's codon-section
    uint32_t acc_idx;
    uint32_t n_records;
};  // 4+8+4+4+4+4+4+4+4 = 40 bytes

// Parse a raw VarNT block and return byte offsets splitting it into
// [other (contig_dict+hdr+cols)] | [bitmaps] | [codons].
inline bool split_raw_block(const uint8_t* raw, uint32_t raw_len,
                            uint32_t& other_len, uint32_t& bmp_len, uint32_t& cdn_len) {
    const uint8_t* p   = raw;
    const uint8_t* end = raw + raw_len;
    uint32_t v = 0; int n;
    // skip contig dict
    n = read_varint(p, end, &v); if (!n) return false; p += n;
    for (uint32_t i = 0; i < v; ++i) {
        uint32_t nlen = 0;
        n = read_varint(p, end, &nlen); if (!n) return false; p += n + nlen;
    }
    // read varnt header: n_recs, contig_b, sstart_b, span_b, bmp_b, cdn_b
    uint32_t n_recs = 0, contig_b = 0, sstart_b = 0, span_b = 0, bmp_b = 0, cdn_b = 0;
    n = read_varint(p, end, &n_recs);   if (!n) return false; p += n;
    n = read_varint(p, end, &contig_b); if (!n) return false; p += n;
    n = read_varint(p, end, &sstart_b); if (!n) return false; p += n;
    n = read_varint(p, end, &span_b);   if (!n) return false; p += n;
    n = read_varint(p, end, &bmp_b);    if (!n) return false; p += n;
    n = read_varint(p, end, &cdn_b);    if (!n) return false; p += n;
    // p now points at start of contig column; other = everything before bitmap
    uint32_t fixed_sz = n_recs * 5;  // qframe(1) + pident(2) + evalue(2)
    other_len = uint32_t(p - raw) + contig_b + sstart_b + span_b + fixed_sz;
    bmp_len   = bmp_b;
    cdn_len   = cdn_b;
    return (other_len + bmp_len + cdn_len == raw_len);
}

// Per-thread sequential writer. Batches raw VarNT payloads into column-major groups,
// then compresses each group as a single zstd frame.
// Frame layout: [other_sec_len(u32) | bmp_sec_len(u32) | other_sec | bmp_sec | cdn_sec]
class PartitionWriter {
    static constexpr size_t BUF_CAP          = 8u * 1024 * 1024;
    static constexpr size_t FRAME_RAW_TARGET = 2u * 1024 * 1024;

    std::string path_;
    int         fd_    = -1;
    uint64_t    fpos_  = 0;
    std::vector<uint8_t> buf_;

    std::vector<uint8_t> other_raw_, bmp_raw_, cdn_raw_;
    struct PendingHog {
        std::string hog_id;
        uint32_t acc_idx, n_records;
        uint32_t other_off, other_len, bmp_off, cdn_off;
    };
    std::vector<PendingHog> pending_;
    std::vector<uint8_t>    cbuf_, frame_assemble_;

    ZSTD_CCtx* cctx_  = nullptr;
    int        level_  = 3;

    std::unordered_map<std::string, std::vector<PartitionIndexExtent>> idx_;

    void do_write(const uint8_t* p, size_t n) {
        size_t done = 0;
        while (done < n) {
            ssize_t r = ::write(fd_, p + done, n - done);
            if (r <= 0) throw std::runtime_error("write failed: " + path_);
            done += size_t(r);
        }
    }
    void maybe_flush() {
        if (buf_.size() >= BUF_CAP) {
            do_write(buf_.data(), buf_.size());
            fpos_ += buf_.size();
            buf_.clear();
        }
    }
    void flush_frame() {
        if (pending_.empty()) return;
        uint32_t other_sec_len = uint32_t(other_raw_.size());
        uint32_t bmp_sec_len   = uint32_t(bmp_raw_.size());
        // assemble column-major frame: [other_sec_len | bmp_sec_len | other | bmp | cdn]
        frame_assemble_.clear();
        frame_assemble_.reserve(8 + other_raw_.size() + bmp_raw_.size() + cdn_raw_.size());
        for (int i = 0; i < 4; ++i) frame_assemble_.push_back(uint8_t(other_sec_len >> (8*i)));
        for (int i = 0; i < 4; ++i) frame_assemble_.push_back(uint8_t(bmp_sec_len   >> (8*i)));
        frame_assemble_.insert(frame_assemble_.end(), other_raw_.begin(), other_raw_.end());
        frame_assemble_.insert(frame_assemble_.end(), bmp_raw_.begin(),   bmp_raw_.end());
        frame_assemble_.insert(frame_assemble_.end(), cdn_raw_.begin(),   cdn_raw_.end());

        size_t bound = ZSTD_compressBound(frame_assemble_.size());
        cbuf_.resize(bound);
        ZSTD_CCtx_setParameter(cctx_, ZSTD_c_compressionLevel, level_);
        ZSTD_CCtx_setParameter(cctx_, ZSTD_c_checksumFlag, 1);
        size_t csz = ZSTD_compress2(cctx_, cbuf_.data(), bound,
                                     frame_assemble_.data(), frame_assemble_.size());
        if (ZSTD_isError(csz))
            throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(csz));
        cbuf_.resize(csz);

        uint64_t frame_off = fpos_ + uint64_t(buf_.size()) + 4;
        uint32_t u32csz = uint32_t(csz);
        const uint8_t* lp = reinterpret_cast<const uint8_t*>(&u32csz);
        buf_.insert(buf_.end(), lp, lp + 4);
        buf_.insert(buf_.end(), cbuf_.data(), cbuf_.data() + csz);
        maybe_flush();

        for (auto& ph : pending_)
            idx_[ph.hog_id].push_back({0u, frame_off, u32csz,
                                        ph.other_off, ph.other_len,
                                        ph.bmp_off, ph.cdn_off,
                                        ph.acc_idx, ph.n_records});
        other_raw_.clear(); bmp_raw_.clear(); cdn_raw_.clear();
        pending_.clear();
    }

public:
    explicit PartitionWriter(const std::string& path, int level = 3)
        : path_(path), level_(level)
    {
        fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd_ < 0)
            throw std::runtime_error("cannot create: " + path + ": " + strerror(errno));
        cctx_ = ZSTD_createCCtx();
        buf_.reserve(BUF_CAP + 256 * 1024);
        other_raw_.reserve(FRAME_RAW_TARGET);
        bmp_raw_.reserve(FRAME_RAW_TARGET / 8);
        cdn_raw_.reserve(FRAME_RAW_TARGET);
    }
    ~PartitionWriter() {
        if (fd_ >= 0) {
            try { flush_all(); } catch (...) {}
            ::close(fd_);
            fd_ = -1;
        }
        if (cctx_) { ZSTD_freeCCtx(cctx_); cctx_ = nullptr; }
    }
    PartitionWriter(const PartitionWriter&) = delete;
    PartitionWriter& operator=(const PartitionWriter&) = delete;

    void append_raw(const std::string& hog_id, uint32_t acc_idx, uint32_t n_records,
                    const uint8_t* raw, uint32_t raw_len) {
        uint32_t other_len = 0, bmp_len = 0, cdn_len = 0;
        if (!split_raw_block(raw, raw_len, other_len, bmp_len, cdn_len))
            throw std::runtime_error("split_raw_block failed for HOG " + hog_id);
        pending_.push_back({hog_id, acc_idx, n_records,
                            uint32_t(other_raw_.size()), other_len,
                            uint32_t(bmp_raw_.size()),
                            uint32_t(cdn_raw_.size())});
        other_raw_.insert(other_raw_.end(), raw,                          raw + other_len);
        bmp_raw_.insert  (bmp_raw_.end(),   raw + other_len,              raw + other_len + bmp_len);
        cdn_raw_.insert  (cdn_raw_.end(),   raw + other_len + bmp_len,    raw + raw_len);
        size_t total = other_raw_.size() + bmp_raw_.size() + cdn_raw_.size();
        if (total >= FRAME_RAW_TARGET) flush_frame();
    }

    void flush_all() {
        flush_frame();
        if (!buf_.empty()) {
            do_write(buf_.data(), buf_.size());
            fpos_ += buf_.size();
            buf_.clear();
        }
    }

    const std::unordered_map<std::string, std::vector<PartitionIndexExtent>>&
    index() const { return idx_; }
};

// Write the partition index: sorted hog_id → extents across N thread files.
// Binary format v4: LHP_INDEX_MAGIC(4) + version(1) + n_threads(u32) + n_hogs(u32)
//   per HOG: hog_id_len(u16) + hog_id + n_extents(u32)
//     per extent: thread_idx(u32)+frame_off(u64)+frame_csz(u32)+other_off(u32)+other_len(u32)
//                 +bmp_off(u32)+cdn_off(u32)+acc_idx(u32)+n_records(u32) [40 bytes]
inline void write_partition_index(
        const std::map<std::string, std::vector<PartitionIndexExtent>>& idx,
        uint32_t n_threads,
        const std::string& out_path) {
    std::vector<uint8_t> buf;
    buf.reserve(idx.size() * 52 + 16);

    buf.insert(buf.end(), LHP_INDEX_MAGIC, LHP_INDEX_MAGIC + 4);
    buf.push_back(LHP_INDEX_VERSION);
    for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(n_threads >> (8 * i)));
    uint32_t n_hogs = uint32_t(idx.size());
    for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(n_hogs >> (8 * i)));

    for (const auto& [hog, exts] : idx) {
        uint16_t hlen = uint16_t(hog.size());
        const uint8_t* hp = reinterpret_cast<const uint8_t*>(&hlen);
        buf.insert(buf.end(), hp, hp + 2);
        buf.insert(buf.end(), hog.begin(), hog.end());
        uint32_t ne = uint32_t(exts.size());
        for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(ne >> (8*i)));
        for (const auto& e : exts) {
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.thread_idx >> (8*i)));
            for (int i = 0; i < 8; ++i) buf.push_back(uint8_t(e.frame_off  >> (8*i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.frame_csz  >> (8*i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.other_off  >> (8*i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.other_len  >> (8*i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.bmp_off    >> (8*i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.cdn_off    >> (8*i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.acc_idx    >> (8*i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.n_records  >> (8*i)));
        }
    }

    UniqueFd fd(::open(out_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644));
    if (fd < 0) throw std::runtime_error("cannot create partition index: " + out_path);
    const uint8_t* p = buf.data(); size_t rem = buf.size(), done = 0;
    while (done < rem) {
        ssize_t r = ::write(fd, p + done, rem - done);
        if (r <= 0) throw std::runtime_error("write failed: " + out_path);
        done += size_t(r);
    }
}

// Load partition index. n_threads_out receives the number of tN.lhp files.
inline std::map<std::string, std::vector<PartitionIndexExtent>>
load_partition_index(const std::string& path, uint32_t& n_threads_out) {
    UniqueFd fd(::open(path.c_str(), O_RDONLY));
    if (fd < 0) throw std::runtime_error("cannot open partition index: " + path);
    struct stat st;
    if (fstat(fd, &st) != 0) throw std::runtime_error("fstat failed: " + path);
    std::vector<uint8_t> buf(size_t(st.st_size));
    if (st.st_size > 0 && !fd_read_exact(fd, buf.data(), buf.size()))
        throw std::runtime_error("truncated partition index: " + path);
    if (buf.size() < 13 || memcmp(buf.data(), LHP_INDEX_MAGIC, 4) != 0)
        throw std::runtime_error("bad partition index magic: " + path);
    if (buf[4] == 1 || buf[4] == 2 || buf[4] == 3)
        throw std::runtime_error("partition index v1/v2/v3 not supported — re-run ingest: " + path);
    if (buf[4] != LHP_INDEX_VERSION)
        throw std::runtime_error("unsupported partition index version " + std::to_string(buf[4]) + ": " + path);

    const uint8_t* cur = buf.data() + 5;
    const uint8_t* eof = buf.data() + buf.size();
    n_threads_out = read_u32_le(cur); cur += 4;
    uint32_t n_hogs = read_u32_le(cur); cur += 4;

    std::map<std::string, std::vector<PartitionIndexExtent>> idx;
    for (uint32_t h = 0; h < n_hogs; ++h) {
        if (cur + 2 > eof) throw std::runtime_error("truncated partition index (hog_id_len)");
        uint16_t hlen = uint16_t(cur[0]) | (uint16_t(cur[1]) << 8); cur += 2;
        if (cur + hlen > eof) throw std::runtime_error("truncated partition index (hog_id)");
        std::string hog_id(reinterpret_cast<const char*>(cur), hlen); cur += hlen;
        if (cur + 4 > eof) throw std::runtime_error("truncated partition index (n_extents)");
        uint32_t ne = read_u32_le(cur); cur += 4;
        auto& exts = idx[hog_id];
        exts.reserve(ne);
        for (uint32_t ei = 0; ei < ne; ++ei) {
            if (cur + 40 > eof) throw std::runtime_error("truncated partition index (extent)");
            PartitionIndexExtent ext;
            ext.thread_idx = read_u32_le(cur); cur += 4;
            ext.frame_off  = 0;
            for (int i = 0; i < 8; ++i) ext.frame_off |= (uint64_t(cur[i]) << (8*i));
            cur += 8;
            ext.frame_csz  = read_u32_le(cur); cur += 4;
            ext.other_off  = read_u32_le(cur); cur += 4;
            ext.other_len  = read_u32_le(cur); cur += 4;
            ext.bmp_off    = read_u32_le(cur); cur += 4;
            ext.cdn_off    = read_u32_le(cur); cur += 4;
            ext.acc_idx    = read_u32_le(cur); cur += 4;
            ext.n_records  = read_u32_le(cur); cur += 4;
            exts.push_back(ext);
        }
    }
    return idx;
}

// Write the standalone LHGA accession registry (varint format, Adler-32 trailer).
inline void write_acc_registry(const std::vector<std::string>& accessions,
                               const std::string& out_path) {
    std::vector<uint8_t> buf;
    buf.insert(buf.end(), LHG_ACC_MAGIC, LHG_ACC_MAGIC + 4);
    uint32_t na = uint32_t(accessions.size());
    for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(na >> (8 * i)));
    for (auto& a : accessions) {
        write_varint(buf, uint32_t(a.size()));
        buf.insert(buf.end(), a.begin(), a.end());
    }
    write_u32(buf, adler32(buf.data(), buf.size()));
    UniqueFd fd(open(out_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644));
    if (fd < 0) throw std::runtime_error("cannot create: " + out_path);
    const uint8_t* p = buf.data(); size_t rem = buf.size(), done = 0;
    while (done < rem) {
        ssize_t r = ::write(fd, p + done, rem - done);
        if (r <= 0) throw std::runtime_error("write failed: " + out_path);
        done += size_t(r);
    }
}

// Load the standalone LHGA accession registry (acc_idx = position).
inline std::vector<std::string> load_partition_acc_registry(const std::string& path) {
    UniqueFd fd(open(path.c_str(), O_RDONLY));
    if (fd < 0) throw std::runtime_error("cannot open acc registry: " + path);
    struct stat st;
    if (fstat(fd, &st) != 0) throw std::runtime_error("fstat failed: " + path);
    std::vector<uint8_t> buf(size_t(st.st_size));
    if (!fd_read_exact(fd, buf.data(), buf.size()))
        throw std::runtime_error("truncated acc registry: " + path);
    if (buf.size() < 8 + 4 || memcmp(buf.data(), LHG_ACC_MAGIC, 4) != 0)
        throw std::runtime_error("bad acc registry magic: " + path);
    if (adler32(buf.data(), buf.size() - 4) != read_u32_le(buf.data() + buf.size() - 4))
        throw std::runtime_error("acc registry checksum mismatch: " + path);
    const uint8_t* p   = buf.data() + 4;
    const uint8_t* end = buf.data() + buf.size() - 4;
    uint32_t n_accs = read_u32_le(p); p += 4;
    std::vector<std::string> accs(n_accs);
    for (auto& a : accs) {
        uint32_t len = 0;
        int n = read_varint(p, end, &len);
        if (!n || p + n + len > end) throw std::runtime_error("truncated acc registry: " + path);
        p += n;
        a.assign(reinterpret_cast<const char*>(p), len);
        p += len;
    }
    return accs;
}

// Scatter all .lhb inputs into out_dir/t{N}.lhp + out_dir/partition.idx + out_dir/acc.registry.
// No per-HOG files; one FD open per thread throughout; no concat step.
inline void partition_lhbs(const std::vector<std::string>& input_paths,
                           const std::string& out_dir,
                           int n_threads_override = 0) {
    std::filesystem::create_directories(out_dir);

    size_t n_threads = std::max<size_t>(1, std::min<size_t>(
        n_threads_override > 0 ? size_t(n_threads_override) : std::thread::hardware_concurrency(),
        std::max<size_t>(1, input_paths.size())));

    std::vector<std::string> accessions;
    {
        std::set<std::string> uniq;
        std::mutex uniq_mtx;
        std::atomic<size_t> hn{0};
        auto hw = [&]() {
            for (;;) {
                size_t fi = hn.fetch_add(1, std::memory_order_relaxed);
                if (fi >= input_paths.size()) break;
                std::string acc;
                if (!scan_batch_file(input_paths[fi], fi,
                        [&](const std::string& a) { acc = a; },
                        [&](BatchBlockRef) {}))
                    throw std::runtime_error("cannot open batch: " + input_paths[fi]);
                std::lock_guard<std::mutex> lk(uniq_mtx);
                uniq.insert(acc);
            }
        };
        std::vector<std::thread> hts;
        hts.reserve(n_threads - 1);
        for (size_t i = 1; i < n_threads; ++i) hts.emplace_back(hw);
        hw();
        for (auto& t : hts) t.join();
        accessions.assign(uniq.begin(), uniq.end());
    }
    std::unordered_map<std::string, uint32_t> acc_to_idx;
    acc_to_idx.reserve(accessions.size() * 2);
    for (uint32_t i = 0; i < uint32_t(accessions.size()); ++i) acc_to_idx[accessions[i]] = i;

    write_acc_registry(accessions, out_dir + "/acc.registry");

    // Scatter pass: each worker writes to out_dir/t{tid}.lhp sequentially.
    std::vector<std::unique_ptr<PartitionWriter>> writers(n_threads);
    std::atomic<size_t> next{0};
    std::atomic<bool>   failed{false};
    std::string         err;
    std::mutex          err_mtx;

    auto worker = [&](size_t tid) {
        writers[tid] = std::make_unique<PartitionWriter>(
            out_dir + "/t" + std::to_string(tid) + ".lhp");
        auto& w = *writers[tid];
        std::vector<uint8_t> cbuf, raw;
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        for (;;) {
            size_t fi = next.fetch_add(1, std::memory_order_relaxed);
            if (fi >= input_paths.size() || failed.load(std::memory_order_relaxed)) break;
            try {
                UniqueFd sfd(open(input_paths[fi].c_str(), O_RDONLY));
                if (sfd < 0) throw std::runtime_error("cannot open: " + input_paths[fi]);
                posix_fadvise(sfd, 0, 0, POSIX_FADV_SEQUENTIAL);
                uint32_t cur_idx = 0;
                bool ok = scan_batch_file(input_paths[fi], fi,
                    [&](const std::string& a) {
                        auto it = acc_to_idx.find(a);
                        if (it == acc_to_idx.end())
                            throw std::runtime_error("acc_id not in registry: " + a);
                        cur_idx = it->second;
                    },
                    [&](BatchBlockRef br) {
                        cbuf.resize(br.compressed_sz);
                        off_t pay_off = br.shard_hdr_offset + 28;
                        if (::pread(sfd, cbuf.data(), br.compressed_sz, pay_off)
                                != ssize_t(br.compressed_sz))
                            throw std::runtime_error("pread payload failed: " + input_paths[fi]);
                        raw.resize(br.raw_sz);
                        size_t rz = ZSTD_decompressDCtx(dctx, raw.data(), br.raw_sz,
                                                         cbuf.data(), br.compressed_sz);
                        if (ZSTD_isError(rz))
                            throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(rz));
                        w.append_raw(br.hog_id, cur_idx, 0, raw.data(), uint32_t(rz));
                    });
                if (!ok) throw std::runtime_error("cannot open batch: " + input_paths[fi]);
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (!failed.exchange(true)) err = e.what();
                break;
            }
            if (input_paths.size() > 100 && (fi + 1) % 1000 == 0)
                std::fprintf(stderr, "  partitioned %zu/%zu inputs\r", fi + 1, input_paths.size());
        }
        w.flush_all();
        ZSTD_freeDCtx(dctx);
    };

    {
        std::vector<std::thread> pool;
        pool.reserve(n_threads - 1);
        for (size_t t = 1; t < n_threads; ++t) pool.emplace_back(worker, t);
        worker(0);
        for (auto& t : pool) t.join();
    }
    if (failed.load()) throw std::runtime_error("partition failed: " + err);
    if (input_paths.size() > 100) std::cerr << "\n";

    // Merge per-thread sidecars into one sorted global index
    std::map<std::string, std::vector<PartitionIndexExtent>> global_idx;
    for (size_t t = 0; t < n_threads; ++t) {
        for (const auto& [hog, extents] : writers[t]->index())
            for (auto e : extents) { e.thread_idx = uint32_t(t); global_idx[hog].push_back(e); }
    }
    write_partition_index(global_idx, uint32_t(n_threads), out_dir + "/partition.idx");

    std::cerr << "partition done: " << global_idx.size() << " HOGs, "
              << accessions.size() << " accessions from " << input_paths.size()
              << " inputs → " << out_dir << "\n";
}

} // namespace lhi
