#pragma once
#include "container.hpp"
#include "format.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <vector>
#include <cstring>
#include <stdexcept>

// .lhb — LHG Batch file: single-file per-accession output from hoki convert.
//
// Format:
//   "LHGB" (4)  file magic
//   version (1) = 1
//   flags   (1) = 0
//   pad[2]
//   acc_id: varint(len) + bytes
//
//   Repeated until EOF:
//     "LHBT" (4) block tag
//     hog_id: varint(len) + bytes
//     ShardBlockHeader (28 bytes, magic "LHSB")
//     compressed payload (header.compressed_sz bytes)
//
//   "LHBE" (4) EOF sentinel

namespace lhi {

constexpr uint8_t LHB_FILE_MAGIC[4]  = {'L','H','G','B'};
constexpr uint8_t LHB_BLOCK_MAGIC[4] = {'L','H','B','T'};
constexpr uint8_t LHB_EOF_MAGIC[4]   = {'L','H','B','E'};
constexpr uint8_t LHB_VERSION        = 1;

struct BatchWriter {
    int         fd = -1;
    std::string acc_id;

    BatchWriter(const std::string& path, const std::string& acc) : acc_id(acc) {
        fd = open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) throw std::runtime_error("cannot create batch: " + path);
        std::vector<uint8_t> buf;
        buf.insert(buf.end(), LHB_FILE_MAGIC, LHB_FILE_MAGIC + 4);
        buf.push_back(LHB_VERSION);
        buf.push_back(0); buf.push_back(0); buf.push_back(0);
        write_varint(buf, uint32_t(acc.size()));
        buf.insert(buf.end(), acc.begin(), acc.end());
        write_all(buf.data(), buf.size());
    }
    ~BatchWriter() { if (fd >= 0) close(fd); }
    BatchWriter(const BatchWriter&) = delete;
    BatchWriter& operator=(const BatchWriter&) = delete;

    void write_block(const std::string& hog_id,
                     const std::vector<std::string>& global_contigs,
                     std::vector<VarNTRecord>& recs,
                     int zstd_level) {
        if (recs.empty()) return;
        std::vector<uint8_t> tag;
        tag.insert(tag.end(), LHB_BLOCK_MAGIC, LHB_BLOCK_MAGIC + 4);
        write_varint(tag, uint32_t(hog_id.size()));
        tag.insert(tag.end(), hog_id.begin(), hog_id.end());
        write_all(tag.data(), tag.size());
        write_shard_block(fd, global_contigs, recs, zstd_level);
    }

    void finalize() {
        write_all(LHB_EOF_MAGIC, 4);
    }

private:
    void write_all(const void* p, size_t n) {
        if (::write(fd, p, n) != ssize_t(n))
            throw std::runtime_error("batch write failed");
    }
};

// Per-block metadata collected during a scan pass (no data read, just offsets).
struct BatchBlockRef {
    std::string hog_id;
    std::string acc_id;
    size_t      batch_file_idx;
    off_t       shard_hdr_offset;  // offset of ShardBlockHeader ("LHSB") in the .lhb file
    uint32_t    compressed_sz;
    uint32_t    raw_sz;
    uint32_t    min_sstart;
    uint32_t    max_send;
};

// Scan a .lhb file without reading compressed payloads.
// header_cb(acc_id) called once.
// block_cb(BatchBlockRef) called per block.
// Returns false if file does not exist.
template<typename HdrCb, typename BlkCb>
inline bool scan_batch_file(const std::string& path, size_t file_idx,
                             HdrCb header_cb, BlkCb block_cb) {
    UniqueFd fd(open(path.c_str(), O_RDONLY));
    if (fd < 0) return false;

    auto read_exact = [&](void* buf, size_t n) -> bool {
        char* p = (char*)buf;
        size_t done = 0;
        while (done < n) {
            ssize_t r = ::read(fd, p + done, n - done);
            if (r <= 0) return false;
            done += r;
        }
        return true;
    };
    auto read_varint_fd = [&](uint32_t& out) -> bool {
        out = 0; int sh = 0;
        for (;;) {
            uint8_t x;
            if (::read(fd, &x, 1) != 1) return false;
            if (sh == 28) {                          // 5th byte: only low 4 bits fit uint32
                if (x & 0x80) return false;          // 6th continuation byte ⇒ overflow
                if (x > 0x0F) return false;          // value > UINT32_MAX
                out |= uint32_t(x) << 28; return true;
            }
            out |= uint32_t(x & 0x7F) << sh; sh += 7;
            if (!(x & 0x80)) break;
        }
        return true;
    };
    auto read_str = [&](std::string& s) -> bool {
        uint32_t len = 0;
        if (!read_varint_fd(len)) return false;
        s.resize(len);
        return read_exact(s.data(), len);
    };

    uint8_t magic[4], ver, flags, pad[2];
    if (!read_exact(magic, 4) || memcmp(magic, LHB_FILE_MAGIC, 4) != 0)
        throw std::runtime_error("bad .lhb magic: " + path);
    if (!read_exact(&ver, 1) || !read_exact(&flags, 1) || !read_exact(pad, 2))
        throw std::runtime_error("truncated .lhb header: " + path);
    if (ver != LHB_VERSION)
        throw std::runtime_error("unsupported .lhb version " + std::to_string(ver) + " in: " + path);
    std::string acc_id;
    if (!read_str(acc_id)) throw std::runtime_error("truncated acc_id: " + path);
    header_cb(acc_id);

    for (;;) {
        uint8_t blk_magic[4];
        ssize_t r = ::read(fd, blk_magic, 4);
        if (r == 0) break;
        if (r != 4) throw std::runtime_error("truncated block tag in: " + path);
        if (memcmp(blk_magic, LHB_EOF_MAGIC, 4) == 0) break;
        if (memcmp(blk_magic, LHB_BLOCK_MAGIC, 4) != 0)
            throw std::runtime_error("bad block magic in: " + path);

        std::string hog_id;
        if (!read_str(hog_id)) throw std::runtime_error("truncated hog_id in: " + path);

        off_t shard_hdr_off = lseek(fd, 0, SEEK_CUR);
        uint8_t hdr[28];
        if (!read_exact(hdr, sizeof(hdr))) throw std::runtime_error("truncated shard hdr in: " + path);
        if (memcmp(hdr, SHARD_BLOCK_MAGIC, 4) != 0)
            throw std::runtime_error("bad shard magic in: " + path);
        if (hdr[4] != SHARD_BLOCK_VERSION)
            throw std::runtime_error("unsupported shard version in: " + path);

        uint32_t compressed_sz = read_u32_le(hdr + 8);
        block_cb(BatchBlockRef{hog_id, acc_id, file_idx, shard_hdr_off,
                               compressed_sz, read_u32_le(hdr + 12),
                               read_u32_le(hdr + 20), read_u32_le(hdr + 24)});

        if (lseek(fd, compressed_sz, SEEK_CUR) < 0)
            throw std::runtime_error("lseek failed in: " + path);
    }

    return true;
}

} // namespace lhi
