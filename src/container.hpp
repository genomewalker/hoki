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

// Shard block: the per-HOG record unit written into .lhb batch files.
//
// Wire format:
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


// Serialize local contig dict + VarNT block, compress, write to open fd.
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

    // Header stats computed before recs is moved into serialize_varnt_block.
    uint32_t n_records = uint32_t(recs.size());
    uint32_t min_sstart = UINT32_MAX, max_send = 0;
    for (const auto& r : recs) {
        min_sstart = std::min(min_sstart, r.sstart);
        max_send   = std::max(max_send,   r.send);
    }

    serialize_varnt_block(raw, std::move(recs));

    // Compress payload (content checksum enabled so the reader validates on decompress).
    size_t bound = ZSTD_compressBound(raw.size());
    std::vector<uint8_t> cbuf(bound);
    ZSTD_CCtx* cctx = ZSTD_createCCtx();
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, zstd_level);
    ZSTD_CCtx_setParameter(cctx, ZSTD_c_checksumFlag, 1);
    size_t csz = ZSTD_compress2(cctx, cbuf.data(), bound, raw.data(), raw.size());
    ZSTD_freeCCtx(cctx);
    if (ZSTD_isError(csz))
        throw std::runtime_error(std::string("zstd: ") + ZSTD_getErrorName(csz));

    // Serialize the 28-byte header explicitly little-endian (host-endian-independent).
    std::vector<uint8_t> hdr_bytes;
    hdr_bytes.reserve(28);
    hdr_bytes.insert(hdr_bytes.end(), SHARD_BLOCK_MAGIC, SHARD_BLOCK_MAGIC + 4);
    hdr_bytes.push_back(SHARD_BLOCK_VERSION);
    hdr_bytes.push_back(0);                 // flags
    hdr_bytes.push_back(0); hdr_bytes.push_back(0);  // reserved[2]
    write_u32(hdr_bytes, uint32_t(csz));
    write_u32(hdr_bytes, uint32_t(raw.size()));
    write_u32(hdr_bytes, n_records);
    write_u32(hdr_bytes, min_sstart);
    write_u32(hdr_bytes, max_send);

    if (::write(fd, hdr_bytes.data(), hdr_bytes.size()) != ssize_t(hdr_bytes.size()) ||
        ::write(fd, cbuf.data(), csz)                   != ssize_t(csz))
        throw std::runtime_error("shard write error");
}



} // namespace lhi
