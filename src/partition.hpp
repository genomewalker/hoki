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
constexpr uint8_t LHP_VERSION         = 2;
constexpr size_t  LHP_HEADER_SZ       = 8;  // magic(4)+ver(1)+flags(1)+pad(2)

constexpr uint8_t LHP_INDEX_MAGIC[4]  = {'L','H','P','I'};
constexpr uint8_t LHP_INDEX_VERSION   = 1;

#pragma pack(push, 1)
struct PartitionEntry {
    uint32_t acc_idx;
    uint32_t compressed_sz;
    uint32_t raw_sz;
    uint32_t n_records;
};
static_assert(sizeof(PartitionEntry) == 16);
#pragma pack(pop)

// Location of a single HOG's PartitionEntry+payload within a thread file.
// entry_offset: byte offset inside tT.lhp of the PartitionEntry struct.
// entry_len:    sizeof(PartitionEntry) + compressed_sz (pread length).
struct PartitionIndexExtent {
    uint32_t thread_idx;
    uint64_t entry_offset;
    uint32_t entry_len;
};

// Per-thread sequential writer. Keeps one fd open for the whole scatter pass;
// HOGs interleave into one sequential 8 MB-buffered stream.  A sidecar index
// records extents so merge-shard can pread individual HOGs without opening 30K files.
class PartitionWriter {
    static constexpr size_t BUF_CAP = 8u * 1024u * 1024u;

    std::string path_;
    int         fd_      = -1;
    uint64_t    fpos_    = 0;   // bytes committed to disk
    std::vector<uint8_t> buf_;

    // hog_id → [(entry_offset, entry_len)]
    std::unordered_map<std::string, std::vector<std::pair<uint64_t, uint32_t>>> idx_;

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

public:
    explicit PartitionWriter(const std::string& path) : path_(path) {
        fd_ = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd_ < 0)
            throw std::runtime_error("cannot create: " + path + ": " + strerror(errno));
        buf_.reserve(BUF_CAP + 256 * 1024);
    }
    ~PartitionWriter() {
        if (fd_ >= 0) {
            try { flush_all(); } catch (...) {}
            ::close(fd_);
            fd_ = -1;
        }
    }
    PartitionWriter(const PartitionWriter&) = delete;
    PartitionWriter& operator=(const PartitionWriter&) = delete;

    void append(const std::string& hog_id, const PartitionEntry& ent,
                const uint8_t* payload) {
        uint16_t hlen = uint16_t(hog_id.size());
        // entry_offset: where PartitionEntry lands in the file (after hog_id prefix)
        uint64_t entry_offset = fpos_ + uint64_t(buf_.size()) + 2 + hlen;
        uint32_t entry_len    = uint32_t(sizeof(PartitionEntry)) + ent.compressed_sz;
        idx_[hog_id].emplace_back(entry_offset, entry_len);

        const uint8_t* hp = reinterpret_cast<const uint8_t*>(&hlen);
        buf_.insert(buf_.end(), hp, hp + 2);
        buf_.insert(buf_.end(),
                    reinterpret_cast<const uint8_t*>(hog_id.data()),
                    reinterpret_cast<const uint8_t*>(hog_id.data()) + hlen);
        const uint8_t* ep = reinterpret_cast<const uint8_t*>(&ent);
        buf_.insert(buf_.end(), ep, ep + sizeof(PartitionEntry));
        buf_.insert(buf_.end(), payload, payload + ent.compressed_sz);
        maybe_flush();
    }

    void flush_all() {
        if (!buf_.empty()) {
            do_write(buf_.data(), buf_.size());
            fpos_ += buf_.size();
            buf_.clear();
        }
    }

    const std::unordered_map<std::string, std::vector<std::pair<uint64_t, uint32_t>>>&
    index() const { return idx_; }
};

// Write the partition index: sorted hog_id → extents across N thread files.
// Binary format: LHP_INDEX_MAGIC(4) + version(1) + n_threads(u32) + n_hogs(u32)
//   per HOG: hog_id_len(u16) + hog_id + n_extents(u16)
//     per extent: thread_idx(u32) + entry_offset(u64) + entry_len(u32)  [16 bytes]
inline void write_partition_index(
        const std::map<std::string, std::vector<PartitionIndexExtent>>& idx,
        uint32_t n_threads,
        const std::string& out_path) {
    std::vector<uint8_t> buf;
    buf.reserve(idx.size() * 36 + 16);

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
        uint16_t ne = uint16_t(exts.size());
        const uint8_t* np = reinterpret_cast<const uint8_t*>(&ne);
        buf.insert(buf.end(), np, np + 2);
        for (const auto& e : exts) {
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.thread_idx   >> (8 * i)));
            for (int i = 0; i < 8; ++i) buf.push_back(uint8_t(e.entry_offset >> (8 * i)));
            for (int i = 0; i < 4; ++i) buf.push_back(uint8_t(e.entry_len    >> (8 * i)));
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
    if (buf[4] != LHP_INDEX_VERSION)
        throw std::runtime_error("unsupported partition index version: " + path);

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
        if (cur + 2 > eof) throw std::runtime_error("truncated partition index (n_extents)");
        uint16_t ne = uint16_t(cur[0]) | (uint16_t(cur[1]) << 8); cur += 2;
        auto& exts = idx[hog_id];
        exts.reserve(ne);
        for (uint16_t ei = 0; ei < ne; ++ei) {
            if (cur + 16 > eof) throw std::runtime_error("truncated partition index (extent)");
            PartitionIndexExtent ext;
            ext.thread_idx   = read_u32_le(cur); cur += 4;
            ext.entry_offset = 0;
            for (int i = 0; i < 8; ++i) ext.entry_offset |= (uint64_t(cur[i]) << (8 * i));
            cur += 8;
            ext.entry_len    = read_u32_le(cur); cur += 4;
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
        std::vector<uint8_t> payload;
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
                        PartitionEntry ent;
                        ent.acc_idx       = cur_idx;
                        ent.compressed_sz = br.compressed_sz;
                        ent.raw_sz        = br.raw_sz;
                        ent.n_records     = 0;
                        payload.resize(br.compressed_sz);
                        off_t pay_off = br.shard_hdr_offset + 28;
                        if (::pread(sfd, payload.data(), br.compressed_sz, pay_off)
                                != ssize_t(br.compressed_sz))
                            throw std::runtime_error("pread payload failed: " + input_paths[fi]);
                        w.append(br.hog_id, ent, payload.data());
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
        for (const auto& [hog, extents] : writers[t]->index()) {
            auto& ge = global_idx[hog];
            for (const auto& [off, len] : extents)
                ge.push_back({uint32_t(t), off, len});
        }
    }
    write_partition_index(global_idx, uint32_t(n_threads), out_dir + "/partition.idx");

    std::cerr << "partition done: " << global_idx.size() << " HOGs, "
              << accessions.size() << " accessions from " << input_paths.size()
              << " inputs → " << out_dir << "\n";
}

} // namespace lhi
