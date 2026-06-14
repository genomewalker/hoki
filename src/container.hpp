#pragma once
#include "format.hpp"
#include <filesystem>
#include <fcntl.h>
#include <unistd.h>
#include <zstd.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>
#include <cstring>

// LHI Container: HOG-sharded append-only directory.
//
// Layout:
//   container.lhc/
//     shards/<HOG_ID>.lhs   — one file per HOG, sequence of shard blocks
//
// Each .lhs file is a sequence of ShardBlock records written independently
// by concurrent workers.  Writers hold a POSIX write lock (fcntl F_SETLKW)
// for the duration of a single append — NFS-safe (unlike flock).
// No global index is needed: the shard path is derived deterministically
// from the HOG ID string.
//
// Shard block wire format:
//   ShardBlockHeader (28 bytes)
//   zstd-compressed payload (compressed_sz bytes)
//
// Payload (decompressed raw_sz bytes):
//   varint  n_contigs
//   for each contig: varint len + len bytes
//   VarNT columnar block (format.hpp: serialize_varnt_block)
//     contig_idx references the local contig list above

namespace lhi {

constexpr uint8_t SHARD_BLOCK_MAGIC[4] = {'L','H','S','B'};
constexpr uint8_t SHARD_BLOCK_VERSION  = 4;

#pragma pack(push, 1)
struct ShardBlockHeader {
    uint8_t  magic[4];       // "LHSB"
    uint8_t  version;        // SHARD_BLOCK_VERSION
    uint8_t  flags;          // reserved, must be 0
    uint8_t  reserved[2];    // reserved, must be 0
    uint32_t compressed_sz;
    uint32_t raw_sz;
    uint32_t n_records;
    uint32_t min_sstart;
    uint32_t max_send;
};
static_assert(sizeof(ShardBlockHeader) == 28);
#pragma pack(pop)

// POSIX write lock (fcntl F_SETLKW) — NFS-safe, unlike flock().
inline void posix_wlock(int fd) {
    struct flock fl{};
    fl.l_type = F_WRLCK; fl.l_whence = SEEK_SET; fl.l_start = 0; fl.l_len = 0;
    if (fcntl(fd, F_SETLKW, &fl) != 0)
        throw std::runtime_error("fcntl F_WRLCK failed");
}
inline void posix_unlock(int fd) {
    struct flock fl{};
    fl.l_type = F_UNLCK; fl.l_whence = SEEK_SET; fl.l_start = 0; fl.l_len = 0;
    fcntl(fd, F_SETLK, &fl);
}

// Convert a HOG ID string to a safe filename (replaces '/' ':' ' ' with '_').
inline std::string hog_to_filename(const std::string& hog_id) {
    std::string s = hog_id;
    for (char& c : s) if (c == '/' || c == ':' || c == ' ') c = '_';
    return s + ".lhs";
}

// Serialize local contig dict + VarNT block, compress, write to open fd.
// Caller must hold flock(LOCK_EX) on fd before calling.
inline void write_shard_block(int fd,
                               const std::vector<std::string>& global_contigs,
                               std::vector<VarNTRecord>& recs,
                               int zstd_level) {
    if (recs.empty()) return;

    // Build local contig dict covering only this batch.
    std::unordered_map<uint32_t, uint32_t> g2l;
    std::vector<std::string> local_contigs;
    for (auto& r : recs) {
        if (!g2l.count(r.contig_idx)) {
            g2l[r.contig_idx] = uint32_t(local_contigs.size());
            local_contigs.push_back(global_contigs[r.contig_idx]);
        }
    }

    // Remap contig_idx to local space in-place (caller clears batch after this call).
    for (auto& r : recs) r.contig_idx = g2l.at(r.contig_idx);

    // Payload = local contig dict + VarNT columnar block.
    size_t raw_est = 5;
    for (auto& s : local_contigs) raw_est += 5 + s.size();
    raw_est += recs.size() * 64;
    std::vector<uint8_t> raw;
    raw.reserve(raw_est);
    write_varint(raw, uint32_t(local_contigs.size()));
    for (auto& s : local_contigs) {
        write_varint(raw, uint32_t(s.size()));
        raw.insert(raw.end(), s.begin(), s.end());
    }
    serialize_varnt_block(raw, recs);

    // Compress payload.
    size_t bound = ZSTD_compressBound(raw.size());
    std::vector<uint8_t> cbuf(bound);
    size_t csz = ZSTD_compress(cbuf.data(), bound, raw.data(), raw.size(), zstd_level);
    if (ZSTD_isError(csz))
        throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(csz));

    ShardBlockHeader hdr{};
    std::memcpy(hdr.magic, SHARD_BLOCK_MAGIC, 4);
    hdr.version       = SHARD_BLOCK_VERSION;
    hdr.compressed_sz = uint32_t(csz);
    hdr.raw_sz        = uint32_t(raw.size());
    hdr.n_records     = uint32_t(recs.size());
    hdr.min_sstart    = UINT32_MAX; hdr.max_send = 0;
    for (const auto& r : recs) {
        hdr.min_sstart = std::min(hdr.min_sstart, r.sstart);
        hdr.max_send   = std::max(hdr.max_send,   r.send);
    }

    if (::write(fd, &hdr,        sizeof(hdr)) != ssize_t(sizeof(hdr)) ||
        ::write(fd, cbuf.data(), csz)         != ssize_t(csz))
        throw std::runtime_error("shard write error");
}

// Open (or create) shard_path, lock, append one block, unlock, close.
inline void flush_hog_shard(const std::filesystem::path& shard_path,
                              const std::vector<std::string>& global_contigs,
                              std::vector<VarNTRecord>& recs,
                              int zstd_level) {
    if (recs.empty()) return;
    int fd = open(shard_path.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);
    if (fd < 0) throw std::runtime_error("cannot open shard: " + shard_path.string());
    posix_wlock(fd);
    try { write_shard_block(fd, global_contigs, recs, zstd_level); }
    catch (...) { posix_unlock(fd); close(fd); throw; }
    posix_unlock(fd);
    close(fd);
}

// Read all shard blocks from a .lhs file.
// For each block calls: cb(local_contigs, records)
// Returns false if the file does not exist (HOG absent from container).
template<typename Cb>
inline bool read_shard_file(const std::string& path, Cb cb) {
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) return false;

    ShardBlockHeader hdr;
    while (true) {
        ssize_t nr = ::read(fd, &hdr, sizeof(hdr));
        if (nr == 0) break;
        if (nr != ssize_t(sizeof(hdr)))
            throw std::runtime_error("truncated shard header: " + path);
        if (std::memcmp(hdr.magic, SHARD_BLOCK_MAGIC, 4) != 0)
            throw std::runtime_error("bad shard magic in: " + path);
        if (hdr.version != SHARD_BLOCK_VERSION)
            throw std::runtime_error("unsupported shard version " +
                std::to_string(hdr.version) + " in: " + path);

        std::vector<uint8_t> cbuf(hdr.compressed_sz);
        if (::read(fd, cbuf.data(), hdr.compressed_sz) != ssize_t(hdr.compressed_sz))
            throw std::runtime_error("truncated shard data: " + path);

        std::vector<uint8_t> raw(hdr.raw_sz);
        size_t rz = ZSTD_decompress(raw.data(), hdr.raw_sz, cbuf.data(), hdr.compressed_sz);
        if (ZSTD_isError(rz))
            throw std::runtime_error(std::string("shard zstd: ") + ZSTD_getErrorName(rz));

        const uint8_t* p = raw.data(), *end = p + raw.size();

        uint32_t n_contigs = 0;
        int n = read_varint(p, end, &n_contigs);
        if (!n) throw std::runtime_error("corrupt contig dict in: " + path);
        p += n;
        std::vector<std::string> local_contigs(n_contigs);
        for (uint32_t i = 0; i < n_contigs; ++i) {
            uint32_t len = 0;
            n = read_varint(p, end, &len);
            if (!n || p + n + len > end) throw std::runtime_error("truncated contig dict");
            p += n;
            local_contigs[i].assign(reinterpret_cast<const char*>(p), len);
            p += len;
        }

        std::vector<VarNTRecord> recs;
        if (!deserialize_varnt_block(p, end, recs))
            throw std::runtime_error("corrupt varnt block in: " + path);

        cb(local_contigs, recs);
    }
    close(fd);
    return true;
}

} // namespace lhi
