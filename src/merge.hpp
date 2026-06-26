#pragma once
#include <ctime>
#include "global.hpp"
#include "batch.hpp"
#include "container.hpp"
#include "partition.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <fcntl.h>
#include <unistd.h>
#include <malloc.h>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <charconv>
#include <optional>
#include <set>
#include <fstream>
#include <queue>
#include <numeric>
#include <deque>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <cstdint>
#include <cstdio>
#include <zstd.h>

static inline uint64_t clock_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return uint64_t(ts.tv_sec) * 1000000000ull + uint64_t(ts.tv_nsec);
}

namespace lhi {

// Fast resident-set read via /proc/self/statm (single line: field 2 = resident pages).
// Cheap enough to poll in the build-admission predicate (vs parsing /proc/self/status).
inline size_t fast_rss_bytes() {
    int fd = ::open("/proc/self/statm", O_RDONLY);
    if (fd < 0) return 0;
    char b[64]; ssize_t n = ::read(fd, b, sizeof(b) - 1); ::close(fd);
    if (n <= 0) return 0;
    b[n] = 0;
    const char* p = b; while (*p && *p != ' ') ++p;   // skip total size, reach resident
    size_t pages = strtoull(p, nullptr, 10);
    return pages * size_t(sysconf(_SC_PAGESIZE));
}

// 8 MB write buffer: reduces per-HOG write() syscalls from 3 to ~0.003
// (one flush per ~2730 HOGs instead of 3 syscalls per HOG).
struct WriteBuffer {
    static constexpr size_t CAP = 8u * 1024u * 1024u;
    std::vector<uint8_t> buf;
    WriteBuffer() { buf.reserve(CAP); }
    void flush_to(int fd, const std::string& path) {
        if (buf.empty()) return;
        const char* p = reinterpret_cast<const char*>(buf.data());
        size_t rem = buf.size(), done = 0;
        while (done < rem) {
            ssize_t r = ::write(fd, p + done, rem - done);
            if (r <= 0) throw std::runtime_error("write failed: " + path);
            done += size_t(r);
        }
        buf.clear();
    }
    void append(int fd, const std::string& path, const void* data, size_t n) {
        if (buf.size() + n > CAP) flush_to(fd, path);
        if (n > CAP) {
            const char* p = reinterpret_cast<const char*>(data);
            size_t rem = n, done = 0;
            while (done < rem) {
                ssize_t r = ::write(fd, p + done, rem - done);
                if (r <= 0) throw std::runtime_error("write failed: " + path);
                done += size_t(r);
            }
        } else {
            buf.insert(buf.end(),
                       reinterpret_cast<const uint8_t*>(data),
                       reinterpret_cast<const uint8_t*>(data) + n);
        }
    }
};

} // namespace lhi (WriteBuffer)

// Merge N inputs (.lhb per-accession batches and/or already-inverted .lhg
// shards) into one position-centric .lhg plus its .lhgi index. For .lhb inputs
// inversion happens here; .lhg inputs are re-inverted into the global acc space.

namespace lhi {

// Unified per-HOG source reference covering both input kinds.
struct MergeRef {
    enum class Kind { LHB, LHG };
    Kind     kind;
    std::string hog_id;
    // LHB only: the single accession owning this block (lookup avoids the
    // per-HOG acc_ids vector; .lhg uses source_accs[source_file_idx] instead).
    std::string lhb_acc_id;
    size_t   source_file_idx = 0;
    // LHB: offset of the ShardBlockHeader ("LHSB") in the .lhb file.
    uint64_t lhb_shard_hdr_offset = 0;
    // LHG: byte range of the HOG entry ("LHHE"…) in the .lhg file.
    uint64_t lhg_data_offset = 0;
    uint64_t lhg_data_length = 0;
};

// Detect input kind by first-4-byte file magic.
//   LHB_FILE_MAGIC "LHGB" → .lhb ;  LHG_FILE_MAGIC "LHGG" → .lhg
inline MergeRef::Kind detect_input_kind(const std::string& path) {
    UniqueFd fd(open(path.c_str(), O_RDONLY));
    if (fd < 0) throw std::runtime_error("cannot open input: " + path);
    uint8_t m[4];
    if (!fd_read_exact(fd, m, 4)) throw std::runtime_error("truncated input: " + path);
    if (memcmp(m, LHB_FILE_MAGIC, 4) == 0) return MergeRef::Kind::LHB;
    if (memcmp(m, LHG_FILE_MAGIC, 4) == 0) return MergeRef::Kind::LHG;
    throw std::runtime_error("unrecognized magic (not .lhb/.lhg): " + path);
}

// Result produced by parallel compute workers; written serially by the main thread.
struct ShardResult {
    std::string hog_id;
    uint32_t stored_sz = 0;  // with flag bits (0x40000000=zstd, 0x80000000=raw)
    uint32_t n_accs    = 0;
    bool     is_v9     = true;   // finalize emits the v9 (uid+trailer) layout
    std::vector<uint8_t> payload;
};

// merge_shard_compute* deleted: read the old per-HOG .lhp format which partition no
// longer writes. Active merge-shard path: run_merge_shard_spill (below).


// Per-thread scratch for the spill build (build_inverted_from_scratch). Allocate once;
// pass by reference so the 152 MB seen_epoch vector and dctx are reused across HOGs.
struct ShardScratch {
    ZSTD_DCtx* dctx = nullptr;
    // Epoch-based seen-accession tracker: seen_epoch[i]==epoch means acc i seen this HOG.
    // Bumping epoch instead of zeroing the vector gives O(1) reset per HOG.
    std::vector<uint32_t> seen_epoch;
    uint32_t epoch = 0;
    std::vector<uint32_t> seen_accs;
    std::vector<std::pair<uint32_t, InvObs>> flat_inv;
    std::unordered_map<uint32_t, std::vector<std::pair<uint32_t,uint32_t>>> acc_intervals_map;
    std::unordered_map<uint32_t, uint8_t> acc_pident_map;
    std::vector<uint8_t> raw_block, extent_buf;
    std::vector<uint32_t> contig_cnum;
    std::vector<VarNTRecord> recs;
    std::vector<uint8_t> inv_raw, hdr_buf, acc_b, cnum_b, codon_b, hog_cbuf;
    // Frame-decode reuse: extents sorted by (thread_idx, frame_off); raw_block holds the
    // last decoded frame so adjacent same-frame extents skip the pread+zstd. Persists
    // across HOG calls on the same worker (only ever a hit, never wrong).
    std::vector<uint32_t> ext_order;
    uint64_t cached_frame_key = ~uint64_t(0);

    explicit ShardScratch(size_t n_acc) : seen_epoch(n_acc, 0) {
        dctx = ZSTD_createDCtx();
    }
    ~ShardScratch() { if (dctx) ZSTD_freeDCtx(dctx); }
    ShardScratch(const ShardScratch&) = delete;
    ShardScratch& operator=(const ShardScratch&) = delete;
};



// ════════════════════════════════════════════════════════════════════════════
// Memory-bounded merge-shard: decode-once + hash-partition spill.
//
// HOGs and frames are densely many-to-many (a HOG spans hundreds of the shard's
// distinct frames), so per-HOG decoding re-decompresses each shared frame thousands
// of times, while decoding every frame into RAM holds the whole working set. Instead:
//   Pass 1: decode each distinct frame ONCE (parallel), deserialize its extents'
//           records, append them to one of B on-disk bucket files keyed by hog_idx%B.
//           Peak RAM = per-thread write buffers.
//   Pass 2: process one bucket at a time — read it (<= budget), group by hog, run the
//           invert+coverage+compress build, emit each HOG. Peak RAM = one bucket.
// Each frame is decoded exactly once; the working set spills to scratch once.
// Spill record (little-endian), appended to bucket[hog_idx % B]:
//   hog_idx u32 | acc_idx u32 | cnum u32 | sstart u32 | send u32 | pu8 u8 | n_obs u32
//   then n_obs * (hog_offset u32, codon u8)
// ════════════════════════════════════════════════════════════════════════════

// Parse one extent from an already-decoded frame; append its spill records to `out`.
inline void spill_extent_into(const uint8_t* frame_ptr, size_t frame_len,
                              const PartitionIndexExtent& ext, uint32_t hog_idx,
                              const std::vector<const std::vector<uint32_t>*>& fd_acc_remap,
                              std::vector<uint32_t>& contig_cnum,
                              std::vector<VarNTRecord>& recs,
                              std::vector<uint8_t>& out) {
    uint32_t acc_idx = ext.acc_idx;
    if (!fd_acc_remap.empty() && ext.thread_idx < fd_acc_remap.size()
            && fd_acc_remap[ext.thread_idx]) {
        const auto& rm = *fd_acc_remap[ext.thread_idx];
        if (acc_idx < rm.size()) acc_idx = rm[acc_idx];
    }
    if (frame_len < 8) throw std::runtime_error("frame too small in spill");
    uint32_t other_sec_len = read_u32_le(frame_ptr);
    uint32_t bmp_sec_len   = read_u32_le(frame_ptr + 4);
    const uint8_t* other_sec = frame_ptr + 8;
    const uint8_t* bmp_sec   = other_sec + other_sec_len;
    const uint8_t* cdn_sec   = bmp_sec   + bmp_sec_len;
    const uint8_t* p   = other_sec + ext.other_off;
    const uint8_t* end = p + ext.other_len;
    uint32_t n_contigs = 0;
    int n = read_varint(p, end, &n_contigs);
    if (!n) throw std::runtime_error("corrupt contig dict in spill");
    p += n;
    contig_cnum.clear(); contig_cnum.resize(n_contigs);
    for (uint32_t j = 0; j < n_contigs; ++j) {
        uint32_t len = 0;
        n = read_varint(p, end, &len);
        if (!n || p + n + len > end) throw std::runtime_error("truncated contig dict in spill");
        p += n;
        std::string_view uid(reinterpret_cast<const char*>(p), len);
        p += len;
        uint32_t cnum = j;
        auto fs = uid.find('_');
        if (fs != std::string_view::npos && fs + 1 < uid.size()) {
            auto fe = uid.find('_', fs + 1);
            auto part = (fe != std::string_view::npos) ? uid.substr(fs+1, fe-fs-1) : uid.substr(fs+1);
            uint32_t v = 0;
            auto cr = std::from_chars(part.data(), part.data() + part.size(), v);
            if (cr.ec == std::errc{}) cnum = v;
        }
        contig_cnum[j] = cnum;
    }
    recs.clear();
    if (!deserialize_varnt_block_split(p, end, bmp_sec + ext.bmp_off, cdn_sec + ext.cdn_off, recs))
        throw std::runtime_error("corrupt varnt block in spill");
    for (const auto& r : recs) {
        if (r.contig_idx >= contig_cnum.size())
            throw std::runtime_error("contig_idx OOB in spill");
        uint32_t cnum = contig_cnum[r.contig_idx];
        uint8_t pu8 = uint8_t(std::min(100.0f, r.pident + 0.5f));
        write_u32(out, hog_idx); write_u32(out, acc_idx); write_u32(out, cnum);
        write_u32(out, r.sstart); write_u32(out, r.send); out.push_back(pu8);
        write_u32(out, uint32_t(r.vars.size()));
        for (const auto& o : r.vars) {
            write_u32(out, o.hog_offset);
            out.push_back(uint8_t(o.packed_codon >> 2));
        }
    }
}

// Build inverted+compressed payload from state already accumulated in sc
// (flat_inv / acc_intervals_map / acc_pident_map / seen_accs).
// Amino-acid coverage of a set of [start,end] intervals after merging overlaps. Shared by
// the merge-shard spill build and the merge_batches reduce build.
inline uint32_t merged_interval_cov(std::vector<std::pair<uint32_t,uint32_t>>& ivals) {
    if (ivals.empty()) return 0;
    std::sort(ivals.begin(), ivals.end());
    uint32_t cov = 0, cs = ivals[0].first, ce = ivals[0].second;
    for (size_t k = 1; k < ivals.size(); ++k) {
        if (ivals[k].first <= ce) { ce = std::max(ce, ivals[k].second); }
        else { cov += ce - cs + 1; cs = ivals[k].first; ce = ivals[k].second; }
    }
    return cov + ce - cs + 1;
}

// Serialize a HOG's inverted block and zstd-compress it; sets stored_sz (with the raw/zstd
// flag bits) and a view of the payload (inv_raw if stored raw, else hog_cbuf). cctx!=null
// uses ZSTD_compress2 (persistent ctx); else one-shot ZSTD_compress at `level`. If
// ser_done_ns!=null it is set right after serialize (before compress) so callers can split
// build vs compress timing. Shared by the spill build and the merge_batches reduce build.
// Streaming serialize+compress: emits the SAME byte layout as serialize_inverted_block
// (acc-dict, coverage, position headers, then acc/cnum/codon columns) but feeds it through
// ZSTD_compressStream2 in <=CHUNK pieces, so neither the whole serialized block nor a
// compressBound output buffer is ever materialized. This is what lets a super-HOG (e.g.
// chr2H18749, ~10 GiB serialized) build within --flush instead of needing ~2x that resident.
// `chunk` is reused scratch (bounded to ~CHUNK); compressed bytes are appended to `out`.
inline void serialize_compress_streamed(
        const std::vector<uint32_t>& local_accs, const std::vector<uint8_t>& local_acc_pident,
        const std::vector<InvPosition>& positions, uint32_t hog_length,
        const std::vector<uint32_t>& covered_aa, ZSTD_CCtx* cctx, int level,
        std::vector<uint8_t>& chunk, std::vector<uint8_t>& out) {
    static constexpr size_t CHUNK = 8u << 20;
    (void)level;   // the persistent cctx already carries the configured level; session-only
    out.clear(); chunk.clear();
    ZSTD_CCtx_reset(cctx, ZSTD_reset_session_only);   // reset keeps compressionLevel param
    std::vector<uint8_t> ob(ZSTD_CStreamOutSize());
    auto push = [&](bool last) {
        ZSTD_inBuffer in{ chunk.data(), chunk.size(), 0 };
        for (;;) {
            ZSTD_outBuffer o{ ob.data(), ob.size(), 0 };
            size_t r = ZSTD_compressStream2(cctx, &o, &in, last ? ZSTD_e_end : ZSTD_e_continue);
            if (ZSTD_isError(r)) throw std::runtime_error(std::string("zstd stream: ") + ZSTD_getErrorName(r));
            out.insert(out.end(), ob.data(), ob.data() + o.pos);
            if (last) { if (r == 0) break; }
            else if (in.pos == in.size) break;
        }
        chunk.clear();
    };
    auto flush_if_big = [&]{ if (chunk.size() >= CHUNK) push(false); };

    uint32_t n_local = uint32_t(local_accs.size());
    write_varint(chunk, n_local);
    uint32_t prev_gacc = 0;
    for (uint32_t li = 0; li < n_local; ++li) {
        write_varint(chunk, local_accs[li] - prev_gacc); prev_gacc = local_accs[li];
        chunk.push_back(li < local_acc_pident.size() ? local_acc_pident[li] : 0);
        flush_if_big();
    }
    write_varint(chunk, hog_length);
    for (uint32_t li = 0; li < n_local; ++li) {
        write_varint(chunk, li < covered_aa.size() ? covered_aa[li] : 0);
        flush_if_big();
    }
    uint32_t n_pos = uint32_t(positions.size());
    write_varint(chunk, n_pos);
    uint32_t prev_pos = 0;
    for (const auto& pos : positions) {
        write_varint(chunk, pos.hog_pos - prev_pos); prev_pos = pos.hog_pos;
        write_varint(chunk, uint32_t(pos.obs.size()));
        flush_if_big();
    }
    for (const auto& pos : positions) {                       // acc column (li deltas, reset per pos)
        uint32_t lacc = 0, prev_li = 0;
        for (const auto& o : pos.obs) { while (local_accs[lacc] < o.acc_idx) ++lacc;
            write_varint(chunk, lacc - prev_li); prev_li = lacc; }
        flush_if_big();
    }
    for (const auto& pos : positions) {                       // cnum column
        for (const auto& o : pos.obs) write_varint(chunk, o.cnum);
        flush_if_big();
    }
    for (const auto& pos : positions) {                       // codon column
        uint32_t n_obs = uint32_t(pos.obs.size());
        uint32_t counts[64] = {};
        for (const auto& o : pos.obs) counts[o.codon_idx & 0x3F]++;
        uint8_t best = 0; uint32_t best_cnt = 0;
        for (int c = 0; c < 64; ++c) if (counts[c] > best_cnt) { best_cnt = counts[c]; best = uint8_t(c); }
        chunk.push_back(best); write_varint(chunk, n_obs - best_cnt);
        uint32_t prev_ord = 0;
        for (uint32_t i = 0; i < n_obs; ++i) if (pos.obs[i].codon_idx != best) { write_varint(chunk, i - prev_ord); prev_ord = i; }
        for (const auto& o : pos.obs) if (o.codon_idx != best) chunk.push_back(o.codon_idx);
        flush_if_big();
    }
    push(true);
}

// ── Memory-bounded streaming merge for big all-.lhg HOGs (v9 format) ──────────
// A position-major cursor walks one decompressed InvBlock.raw WITHOUT decode_block,
// pre-scanning to locate the acc/cnum/codon columns, then yielding a position's obs
// (merged-local acc, cnum, codon) on demand. v9 replaces the acc+cnum columns with a
// unitig trailer + a Δ-uid column (see build below).
struct LhgPosCursor {
    std::vector<uint32_t> pos_vals, nobs;
    size_t n_pos = 0, cur = 0;
    uint32_t base = 0;                                   // merged-local base for this block
    bool v9 = false;                                     // input block layout
    const std::vector<uint32_t> *tr_acc = nullptr, *tr_cnum = nullptr;  // v9 trailer
    const uint8_t *end = nullptr, *acc0 = nullptr, *cnum0 = nullptr, *codon0 = nullptr;
    const uint8_t *pa = nullptr, *pc = nullptr, *pk = nullptr;

    void init(const InvBlock& blk, uint32_t base_) {
        base = base_; v9 = blk.is_v9; tr_acc = &blk.trailer_acc; tr_cnum = &blk.trailer_cnum;
        const uint8_t* q = blk.pos_ptr; end = blk.end;
        uint32_t np = 0; int n = read_varint(q, end, &np); if (!n) throw std::runtime_error("cursor: pos count"); q += n;
        n_pos = np; pos_vals.resize(np); nobs.resize(np);
        uint32_t prev = 0; uint64_t tot = 0;
        for (uint32_t i = 0; i < np; ++i) {
            uint32_t d = 0; n = read_varint(q, end, &d); if (!n) throw std::runtime_error("cursor: pos delta"); q += n;
            prev += d; pos_vals[i] = prev;
            uint32_t no = 0; n = read_varint(q, end, &no); if (!n) throw std::runtime_error("cursor: n_obs"); q += n;
            nobs[i] = no; tot += no;
        }
        acc0 = q; const uint8_t* z = q;
        // v8: acc col + cnum col + codon col. v9: uid col + codon col (acc/cnum via trailer).
        for (uint64_t j = 0; j < tot; ++j) { uint32_t d; int m = read_varint(z, end, &d); if (!m) throw std::runtime_error("cursor: col skip"); z += m; }
        cnum0 = z;
        if (!v9) {
            for (uint64_t j = 0; j < tot; ++j) { uint32_t d; int m = read_varint(z, end, &d); if (!m) throw std::runtime_error("cursor: cnum skip"); z += m; }
        }
        codon0 = z;
        reset();
    }
    void reset() { pa = acc0; pc = cnum0; pk = codon0; cur = 0; }
    bool done() const { return cur >= n_pos; }
    uint32_t pos() const { return pos_vals[cur]; }
    // Append the current position's obs (merged-local acc, cnum, codon) to parallel buckets; advance.
    void emit(std::vector<uint32_t>& oa, std::vector<uint32_t>& oc, std::vector<uint8_t>& ok) {
        uint32_t no = nobs[cur];
        if (v9) {
            uint32_t prev_uid = 0;
            for (uint32_t i = 0; i < no; ++i) {
                uint32_t d = 0; int m = read_varint(pa, end, &d); pa += m; prev_uid += d;
                if (prev_uid >= tr_acc->size()) throw std::runtime_error("cursor: uid OOB");
                oa.push_back(base + (*tr_acc)[prev_uid]); oc.push_back((*tr_cnum)[prev_uid]);
            }
        } else {
            uint32_t prev_li = 0;
            for (uint32_t i = 0; i < no; ++i) { uint32_t d = 0; int m = read_varint(pa, end, &d); pa += m; prev_li += d; oa.push_back(base + prev_li); }
            for (uint32_t i = 0; i < no; ++i) { uint32_t c = 0; int m = read_varint(pc, end, &c); pc += m; oc.push_back(c); }
        }
        if (pk >= end) throw std::runtime_error("cursor: codon trunc");
        uint8_t cons = *pk++; size_t b = ok.size(); ok.resize(b + no, cons);
        uint32_t nvar = 0; int m = read_varint(pk, end, &nvar); pk += m;
        if (nvar > no) throw std::runtime_error("cursor: nvar > nobs");
        static thread_local std::vector<uint32_t> ords; ords.resize(nvar);
        uint32_t prev_ord = 0;
        for (uint32_t v = 0; v < nvar; ++v) { uint32_t d = 0; int mm = read_varint(pk, end, &d); pk += mm; prev_ord += d; ords[v] = prev_ord; }
        for (uint32_t v = 0; v < nvar; ++v) { if (ords[v] >= no) throw std::runtime_error("cursor: codon ord OOB"); uint8_t vc = *pk++; ok[b + ords[v]] = vc; }
        ++cur;
    }
};

// v9 build: union N position-sorted .lhg blocks into one payload. A unitig (acc,cnum) is
// constant across its ~150 positions, so we fold acc+cnum into a per-HOG trailer of distinct
// (acc,cnum) pairs (sorted) and store, per obs, a Δ-encoded uid (index into the trailer, obs
// sorted by (acc,cnum) within position). The Δ-uid pattern repeats across a unitig's positions
// so it compresses like the acc column — ~60% smaller blocks, lossless. Memory stays
// O(Σ raw blocks + unitig dict + one-position window).
inline void build_inverted_streamed_lhg(
        std::vector<InvBlock>& blocks,
        const std::vector<const std::vector<uint32_t>*>& remaps,
        const std::vector<uint32_t>& local_accs, const std::vector<uint8_t>& local_acc_pident,
        const std::vector<uint32_t>& covered_aa, uint32_t hog_length,
        ZSTD_CCtx* cctx, std::vector<uint8_t>& out) {
    const size_t N = blocks.size();
    uint32_t n_local = uint32_t(local_accs.size());

    std::vector<uint8_t> head;
    // Header: acc-dict + coverage (sections 1-2, unchanged from v8).
    write_varint(head, n_local);
    { uint32_t prev = 0; for (uint32_t li = 0; li < n_local; ++li) {
        write_varint(head, local_accs[li] - prev); prev = local_accs[li];
        head.push_back(li < local_acc_pident.size() ? local_acc_pident[li] : 0); } }
    write_varint(head, hog_length);
    for (uint32_t li = 0; li < n_local; ++li) write_varint(head, li < covered_aa.size() ? covered_aa[li] : 0);

    std::vector<uint32_t> base(N, 0);
    for (size_t bi = 0; bi < N; ++bi) {
        const auto& rm = *remaps[bi];
        uint32_t g0 = rm.empty() ? 0 : rm[0];
        base[bi] = uint32_t(std::lower_bound(local_accs.begin(), local_accs.end(), g0) - local_accs.begin());
    }

    std::vector<LhgPosCursor> cur(N);
    for (size_t bi = 0; bi < N; ++bi) cur[bi].init(blocks[bi], base[bi]);

    struct HItem { uint32_t pos; size_t src; };
    auto cmp = [](const HItem& a, const HItem& b) { return a.pos > b.pos; };
    auto run_pass = [&](auto per_pos) {
        std::priority_queue<HItem, std::vector<HItem>, decltype(cmp)> heap(cmp);
        for (size_t bi = 0; bi < N; ++bi) if (!cur[bi].done()) heap.push({cur[bi].pos(), bi});
        std::vector<std::vector<uint32_t>> ba(N), bc(N);
        std::vector<std::vector<uint8_t>> bk(N);
        while (!heap.empty()) {
            uint32_t cp = heap.top().pos;
            for (size_t bi = 0; bi < N; ++bi) { ba[bi].clear(); bc[bi].clear(); bk[bi].clear(); }
            while (!heap.empty() && heap.top().pos == cp) {
                size_t bi = heap.top().src; heap.pop();
                cur[bi].emit(ba[bi], bc[bi], bk[bi]);
                if (!cur[bi].done()) heap.push({cur[bi].pos(), bi});
            }
            per_pos(cp, ba, bc, bk);
        }
    };
    auto pack = [](uint32_t a, uint32_t c) -> uint64_t { return (uint64_t(a) << 32) | c; };

    // Pass 1: collect distinct unitigs (merged-local acc, cnum).
    std::unordered_set<uint64_t> uni;
    run_pass([&](uint32_t, std::vector<std::vector<uint32_t>>& ba, std::vector<std::vector<uint32_t>>& bc,
                 std::vector<std::vector<uint8_t>>&) {
        for (size_t bi = 0; bi < N; ++bi)
            for (size_t i = 0; i < ba[bi].size(); ++i) uni.insert(pack(ba[bi][i], bc[bi][i]));
    });

    // Sort unitigs -> dense uid; build trailer (Δacc, cnum) + a packed->uid lookup.
    std::vector<uint64_t> us(uni.begin(), uni.end());
    std::sort(us.begin(), us.end());
    std::unordered_set<uint64_t>().swap(uni);
    std::unordered_map<uint64_t, uint32_t> uidmap; uidmap.reserve(us.size() * 2);
    std::vector<uint8_t> trailer; write_varint(trailer, uint32_t(us.size()));
    { uint32_t prev_acc = 0;
      for (uint32_t u = 0; u < us.size(); ++u) {
          uint32_t a = uint32_t(us[u] >> 32), c = uint32_t(us[u] & 0xFFFFFFFF);
          write_varint(trailer, a - prev_acc); prev_acc = a; write_varint(trailer, c);
          uidmap.emplace(us[u], u);
      } }
    std::vector<uint64_t>().swap(us);

    // Pass 2: emit Δ-uid column + codon column (obs sorted by (acc,cnum)=uid within position).
    for (size_t bi = 0; bi < N; ++bi) cur[bi].reset();
    std::vector<uint8_t> pos_hdr, uid_col, codon_col;
    uint32_t n_pos = 0, prev_pos = 0;
    struct Obs { uint32_t acc, cnum; uint8_t codon; };
    std::vector<Obs> win; win.reserve(1u << 20);
    run_pass([&](uint32_t cp, std::vector<std::vector<uint32_t>>& ba, std::vector<std::vector<uint32_t>>& bc,
                 std::vector<std::vector<uint8_t>>& bk) {
        win.clear();
        for (size_t bi = 0; bi < N; ++bi)
            for (size_t i = 0; i < ba[bi].size(); ++i) win.push_back({ba[bi][i], bc[bi][i], bk[bi][i]});
        std::sort(win.begin(), win.end(), [](const Obs& x, const Obs& y) {
            return x.acc != y.acc ? x.acc < y.acc : x.cnum < y.cnum; });
        uint32_t no = uint32_t(win.size());
        write_varint(pos_hdr, cp - prev_pos); prev_pos = cp; write_varint(pos_hdr, no);
        uint32_t prev_uid = 0;
        for (uint32_t i = 0; i < no; ++i) { uint32_t uid = uidmap[pack(win[i].acc, win[i].cnum)]; write_varint(uid_col, uid - prev_uid); prev_uid = uid; }
        uint32_t counts[64] = {};
        for (uint32_t i = 0; i < no; ++i) counts[win[i].codon & 0x3F]++;
        uint8_t best = 0; uint32_t best_cnt = 0;
        for (int c = 0; c < 64; ++c) if (counts[c] > best_cnt) { best_cnt = counts[c]; best = uint8_t(c); }
        codon_col.push_back(best); write_varint(codon_col, no - best_cnt);
        uint32_t prev_ord = 0;
        for (uint32_t i = 0; i < no; ++i) if (win[i].codon != best) { write_varint(codon_col, i - prev_ord); prev_ord = i; }
        for (uint32_t i = 0; i < no; ++i) if (win[i].codon != best) codon_col.push_back(win[i].codon);
        ++n_pos;
    });

    for (auto& b : blocks) { std::vector<uint8_t>().swap(b.raw); b.pos_ptr = b.end = nullptr; }

    std::vector<uint8_t> ph; write_varint(ph, n_pos);
    ph.insert(ph.end(), pos_hdr.begin(), pos_hdr.end());

    // Layout: [acc-dict+coverage][trailer][pos headers][Δ-uid column][codon column].
    out.clear();
    ZSTD_CCtx_reset(cctx, ZSTD_reset_session_only);
    std::vector<uint8_t> ob(ZSTD_CStreamOutSize());
    auto feed = [&](const std::vector<uint8_t>& buf, bool last) {
        ZSTD_inBuffer in{ buf.data(), buf.size(), 0 };
        for (;;) {
            ZSTD_outBuffer o{ ob.data(), ob.size(), 0 };
            size_t r = ZSTD_compressStream2(cctx, &o, &in, last ? ZSTD_e_end : ZSTD_e_continue);
            if (ZSTD_isError(r)) throw std::runtime_error(std::string("zstd merge stream: ") + ZSTD_getErrorName(r));
            out.insert(out.end(), ob.data(), ob.data() + o.pos);
            if (last) { if (r == 0) break; } else if (in.pos == in.size) break;
        }
    };
    feed(head, false); feed(trailer, false); feed(ph, false); feed(uid_col, false); feed(codon_col, true);
}

// v9 serialize from a materialized `positions` vector (obs hold GLOBAL acc_idx). Builds the
// same uid-trailer + Δ-uid layout as build_inverted_streamed_lhg, mapping global acc -> local
// via local_accs. Used by finalize for the small-HOG (merge) and merge-shard build paths.
inline void serialize_inverted_block_v9(
        std::vector<uint8_t>& raw,
        const std::vector<uint32_t>& local_accs, const std::vector<uint8_t>& local_acc_pident,
        const std::vector<InvPosition>& positions, uint32_t hog_length,
        const std::vector<uint32_t>& covered_aa) {
    raw.clear();
    uint32_t n_local = uint32_t(local_accs.size());
    write_varint(raw, n_local);
    { uint32_t prev = 0; for (uint32_t li = 0; li < n_local; ++li) {
        write_varint(raw, local_accs[li] - prev); prev = local_accs[li];
        raw.push_back(li < local_acc_pident.size() ? local_acc_pident[li] : 0); } }
    write_varint(raw, hog_length);
    for (uint32_t li = 0; li < n_local; ++li) write_varint(raw, li < covered_aa.size() ? covered_aa[li] : 0);

    auto to_local = [&](uint32_t g) -> uint32_t {
        return uint32_t(std::lower_bound(local_accs.begin(), local_accs.end(), g) - local_accs.begin()); };
    auto pack = [](uint32_t a, uint32_t c) -> uint64_t { return (uint64_t(a) << 32) | c; };

    std::unordered_set<uint64_t> uni;
    for (const auto& pos : positions)
        for (const auto& o : pos.obs) uni.insert(pack(to_local(o.acc_idx), o.cnum));
    std::vector<uint64_t> us(uni.begin(), uni.end());
    std::sort(us.begin(), us.end());
    std::unordered_set<uint64_t>().swap(uni);
    std::unordered_map<uint64_t, uint32_t> uidmap; uidmap.reserve(us.size() * 2);
    write_varint(raw, uint32_t(us.size()));
    { uint32_t prev_acc = 0;
      for (uint32_t u = 0; u < us.size(); ++u) {
          uint32_t a = uint32_t(us[u] >> 32), c = uint32_t(us[u] & 0xFFFFFFFF);
          write_varint(raw, a - prev_acc); prev_acc = a; write_varint(raw, c);
          uidmap.emplace(us[u], u);
      } }
    std::vector<uint64_t>().swap(us);

    std::vector<uint8_t> ph, uc, kc;
    write_varint(ph, uint32_t(positions.size()));
    uint32_t prev_pos = 0;
    std::vector<std::pair<uint64_t, uint8_t>> obs;
    for (const auto& pos : positions) {
        obs.clear(); obs.reserve(pos.obs.size());
        for (const auto& o : pos.obs) obs.push_back({pack(to_local(o.acc_idx), o.cnum), o.codon_idx});
        std::sort(obs.begin(), obs.end(), [](const std::pair<uint64_t,uint8_t>& x, const std::pair<uint64_t,uint8_t>& y){ return x.first < y.first; });
        write_varint(ph, pos.hog_pos - prev_pos); prev_pos = pos.hog_pos;
        write_varint(ph, uint32_t(obs.size()));
        uint32_t prev_uid = 0;
        for (auto& x : obs) { uint32_t uid = uidmap[x.first]; write_varint(uc, uid - prev_uid); prev_uid = uid; }
        uint32_t no = uint32_t(obs.size());
        uint32_t counts[64] = {};
        for (auto& x : obs) counts[x.second & 0x3F]++;
        uint8_t best = 0; uint32_t best_cnt = 0;
        for (int c = 0; c < 64; ++c) if (counts[c] > best_cnt) { best_cnt = counts[c]; best = uint8_t(c); }
        kc.push_back(best); write_varint(kc, no - best_cnt);
        uint32_t prev_ord = 0;
        for (uint32_t i = 0; i < no; ++i) if (obs[i].second != best) { write_varint(kc, i - prev_ord); prev_ord = i; }
        for (uint32_t i = 0; i < no; ++i) if (obs[i].second != best) kc.push_back(obs[i].second);
    }
    raw.insert(raw.end(), ph.begin(), ph.end());
    raw.insert(raw.end(), uc.begin(), uc.end());
    raw.insert(raw.end(), kc.begin(), kc.end());
}

// Always emits a v9 block (uid trailer + Δ-uid column). The caller writes the LHH2 magic.
inline void finalize_inverted_payload(
        const std::vector<uint32_t>& local_accs, const std::vector<uint8_t>& local_acc_pident,
        const std::vector<InvPosition>& positions, uint32_t hog_length,
        const std::vector<uint32_t>& covered_aa_v,
        std::vector<uint8_t>& inv_raw, std::vector<uint8_t>& hog_cbuf,
        std::vector<uint8_t>& sb_hdr, std::vector<uint8_t>& sb_acc,
        std::vector<uint8_t>& sb_cnum, std::vector<uint8_t>& sb_codon,
        ZSTD_CCtx* cctx, int level, uint64_t* ser_done_ns,
        uint32_t& stored_sz, const uint8_t*& payload, size_t& payload_sz) {
    (void)sb_hdr; (void)sb_acc; (void)sb_cnum; (void)sb_codon;
    inv_raw.clear();
    serialize_inverted_block_v9(inv_raw, local_accs, local_acc_pident, positions, hog_length, covered_aa_v);
    if (ser_done_ns) *ser_done_ns = clock_ns();
    size_t raw_sz = inv_raw.size();
    size_t bound = ZSTD_compressBound(raw_sz);
    hog_cbuf.resize(bound);
    size_t csz = cctx ? ZSTD_compress2(cctx, hog_cbuf.data(), bound, inv_raw.data(), raw_sz)
                      : ZSTD_compress(hog_cbuf.data(), bound, inv_raw.data(), raw_sz, level);
    bool use_raw = ZSTD_isError(csz) || csz >= raw_sz;
    // Inline word carries flags only; the true (uint64) length lives in the index, so
    // blocks may exceed the old 30-bit inline size limit (merged super-HOGs hit ~2 GiB).
    if (use_raw) { stored_sz = 0x80000000u; payload = inv_raw.data();  payload_sz = raw_sz; }
    else         { stored_sz = 0x40000000u; payload = hog_cbuf.data(); payload_sz = csz; }
}

inline ShardResult build_inverted_from_scratch(const std::string& hog_id, ShardScratch& sc,
                                               uint32_t n_accs, uint32_t hog_length,
                                               int out_zstd_level) {
    ShardResult res; res.hog_id = hog_id;
    auto& flat_inv = sc.flat_inv;
    uint32_t max_pos = hog_length > 0 ? hog_length - 1 : 0;
    for (const auto& e : flat_inv) if (e.first > max_pos) max_pos = e.first;
    std::vector<uint32_t> cs_count(max_pos + 1, 0), cs_pidx(max_pos + 1, UINT32_MAX);
    for (const auto& e : flat_inv) ++cs_count[e.first];
    std::vector<InvPosition> positions;
    for (uint32_t pp = 0; pp <= max_pos; ++pp) {
        if (!cs_count[pp]) continue;
        cs_pidx[pp] = uint32_t(positions.size());
        InvPosition ip; ip.hog_pos = pp; ip.obs.reserve(cs_count[pp]);
        positions.push_back(std::move(ip));
    }
    for (auto& e : flat_inv) positions[cs_pidx[e.first]].obs.push_back(std::move(e.second));
    for (auto& pos : positions)
        std::sort(pos.obs.begin(), pos.obs.end(),
                  [](const InvObs& a, const InvObs& b) { return a.acc_idx < b.acc_idx; });
    std::vector<uint32_t> local_accs(sc.seen_accs);
    std::sort(local_accs.begin(), local_accs.end());
    std::vector<uint8_t>  local_acc_pident(local_accs.size(), 0);
    std::vector<uint32_t> covered_aa_v(local_accs.size(), 0);
    for (size_t li = 0; li < local_accs.size(); ++li) {
        uint32_t gacc = local_accs[li];
        auto it = sc.acc_pident_map.find(gacc);
        if (it != sc.acc_pident_map.end()) local_acc_pident[li] = it->second;
        auto it3 = sc.acc_intervals_map.find(gacc);
        if (it3 != sc.acc_intervals_map.end()) covered_aa_v[li] += merged_interval_cov(it3->second);
    }
    uint32_t stored_sz; const uint8_t* payload; size_t payload_sz;
    finalize_inverted_payload(local_accs, local_acc_pident, positions, hog_length, covered_aa_v,
                              sc.inv_raw, sc.hog_cbuf, sc.hdr_buf, sc.acc_b, sc.cnum_b, sc.codon_b,
                              nullptr, out_zstd_level, nullptr, stored_sz, payload, payload_sz);
    res.stored_sz = stored_sz; res.n_accs = n_accs;
    res.payload.assign(payload, payload + payload_sz);
    return res;
}

using ShardHogList =
    std::vector<const std::pair<const std::string, std::vector<PartitionIndexExtent>>*>;

// Memory-bounded merge of one shard set. emit() is called once per built HOG and
// must be thread-safe. Peak RAM is bounded by `budget` regardless of shard size.
inline void run_merge_shard_spill(
        const ShardHogList& hog_list,
        const std::vector<int>& tfd_ints,
        const std::vector<const std::vector<uint32_t>*>& fd_acc_remap,
        size_t n_accessions, size_t budget, size_t nt, int out_zstd_level,
        const std::string& spill_dir, bool do_profile,
        const std::function<void(ShardResult&&)>& emit) {

    // ── group extents by frame ───────────────────────────────────────────────
    struct FExt { uint32_t hog_idx; const PartitionIndexExtent* ext; };
    struct FrameJob { int tfd; uint64_t off; uint32_t csz; std::vector<FExt> exts; };
    std::vector<FrameJob> frames;
    std::unordered_map<uint64_t, uint32_t> fmap;
    uint64_t sum_csz = 0;
    for (uint32_t hi = 0; hi < uint32_t(hog_list.size()); ++hi)
        for (const auto& e : hog_list[hi]->second) {
            uint64_t key = (uint64_t(e.thread_idx) << 48)
                         | (uint64_t(e.frame_off) & 0xFFFFFFFFFFFFull);
            auto it = fmap.find(key);
            uint32_t fi;
            if (it == fmap.end()) {
                fi = uint32_t(frames.size()); fmap.emplace(key, fi);
                frames.push_back({tfd_ints[e.thread_idx], e.frame_off, e.frame_csz, {}});
                sum_csz += e.frame_csz;
            } else fi = it->second;
            frames[fi].exts.push_back({hi, &e});
        }

    // Pass-2 peak RAM = the bucket read whole (`data`) + the `groups` offset-map (~0.3-0.5x
    // bucket) + nt per-HOG build scratch — all resident together. --flush must bound the SUM,
    // not just the bucket, and the spill-size estimate is imprecise. So size buckets to ~1/3 of
    // --flush, reserving the rest for the groups map + build scratch + estimate error.
    size_t est_spill = sum_csz * 10;   // spill record format expands ~8x
    size_t bucket_budget = std::max<size_t>(budget / 3, size_t(64) << 20);
    size_t B = std::max<size_t>(1, (est_spill + bucket_budget - 1) / bucket_budget);
    std::vector<UniqueFd> bfd; bfd.reserve(B);
    std::vector<std::unique_ptr<std::mutex>> bmtx;
    for (size_t b = 0; b < B; ++b) {
        std::string bp = spill_dir + "/bucket." + std::to_string(b);
        int fd = open(bp.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) throw std::runtime_error("cannot create spill bucket: " + bp);
        bfd.emplace_back(fd);
        bmtx.push_back(std::make_unique<std::mutex>());
    }
    std::cerr << "merge-shard: spill " << frames.size() << " frames -> " << B
              << " buckets (budget " << (budget >> 20) << " MiB, -t " << nt << ")\n";

    std::atomic<bool> failed{false}; std::string err; std::mutex err_mtx;
    uint64_t _t0 = clock_ns();

    // ── pass 1: decode each frame once, scatter records to bucket files ───────
    {
        std::atomic<size_t> fnext{0};
        auto p1 = [&]() {
            ZSTD_DCtx* dctx = ZSTD_createDCtx();
            std::vector<uint8_t> cbuf, frame;
            std::vector<uint32_t> ccnum; std::vector<VarNTRecord> recs;
            std::vector<std::vector<uint8_t>> lbuf(B);
            const size_t CHUNK = 8u * 1024 * 1024;
            auto flush_b = [&](size_t b) {
                if (lbuf[b].empty()) return;
                std::lock_guard<std::mutex> lk(*bmtx[b]);
                const uint8_t* d = lbuf[b].data(); size_t rem = lbuf[b].size(), done = 0;
                while (done < rem) {
                    ssize_t w = ::write(bfd[b], d + done, rem - done);
                    if (w <= 0) throw std::runtime_error("bucket write failed");
                    done += size_t(w);
                }
                lbuf[b].clear();
            };
            try {
                for (;;) {
                    size_t i = fnext.fetch_add(1, std::memory_order_relaxed);
                    if (i >= frames.size() || failed.load(std::memory_order_relaxed)) break;
                    auto& fj = frames[i];
                    cbuf.resize(fj.csz);
                    if (::pread(fj.tfd, cbuf.data(), fj.csz, off_t(fj.off)) != ssize_t(fj.csz))
                        throw std::runtime_error("pread frame failed");
                    uint64_t rsz = ZSTD_getFrameContentSize(cbuf.data(), fj.csz);
                    if (rsz == ZSTD_CONTENTSIZE_UNKNOWN || rsz == ZSTD_CONTENTSIZE_ERROR)
                        throw std::runtime_error("bad frame content size");
                    frame.resize(size_t(rsz));
                    size_t rz = ZSTD_decompressDCtx(dctx, frame.data(), size_t(rsz), cbuf.data(), fj.csz);
                    if (ZSTD_isError(rz)) throw std::runtime_error(ZSTD_getErrorName(rz));
                    for (auto& fe : fj.exts) {
                        size_t b = fe.hog_idx % B;
                        spill_extent_into(frame.data(), frame.size(), *fe.ext, fe.hog_idx,
                                          fd_acc_remap, ccnum, recs, lbuf[b]);
                        if (lbuf[b].size() >= CHUNK) flush_b(b);
                    }
                }
                for (size_t b = 0; b < B; ++b) flush_b(b);
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (!failed.exchange(true)) err = e.what();
            }
            ZSTD_freeDCtx(dctx);
        };
        std::vector<std::thread> pool;
        for (size_t t = 1; t < nt; ++t) pool.emplace_back(p1);
        p1();
        for (auto& t : pool) t.join();
        if (failed.load()) throw std::runtime_error("spill pass 1: " + err);
    }
    uint64_t _t_p1 = clock_ns() - _t0;
    uint64_t _t_read = 0, _t_build = 0;

    // ── pass 2: build one bucket at a time (parallel HOGs within a bucket) ────
    for (size_t b = 0; b < B; ++b) {
        off_t sz = lseek(bfd[b], 0, SEEK_END);
        if (sz <= 0) continue;
        std::vector<uint8_t> data; data.resize(size_t(sz));
        uint64_t _r0 = clock_ns();
        { size_t got = 0;
          while (got < size_t(sz)) {
              ssize_t r = ::pread(bfd[b], data.data() + got, size_t(sz) - got, off_t(got));
              if (r <= 0) throw std::runtime_error("bucket read failed");
              got += size_t(r);
          } }
        _t_read += clock_ns() - _r0;
        // group record spans by hog_idx
        std::unordered_map<uint32_t, std::vector<uint64_t>> groups;  // hog_idx -> record start offsets (64-bit: buckets exceed 4GB at scale)
        std::unordered_map<uint32_t, uint64_t> group_bytes;          // hog_idx -> spill bytes (build-scratch estimate)
        size_t off = 0;
        while (off + 25 <= data.size()) {
            uint64_t rstart = off;
            uint32_t hog_idx = read_u32_le(&data[off]);
            off += 4 + 4 + 4 + 4 + 4 + 1;  // hog,acc,cnum,sstart,send,pu8
            uint32_t nobs = read_u32_le(&data[off]); off += 4;
            off += size_t(nobs) * 5;
            groups[hog_idx].push_back(rstart);
            group_bytes[hog_idx] += off - rstart;
        }
        std::vector<uint32_t> hkeys; hkeys.reserve(groups.size());
        for (auto& kv : groups) hkeys.push_back(kv.first);
        // LPT scheduling: biggest HOGs first so a giant one can't strand a thread
        // at the tail while the rest idle (the pass2-build Amdahl bottleneck).
        std::sort(hkeys.begin(), hkeys.end(),
                  [&](uint32_t a, uint32_t b){ return groups[a].size() > groups[b].size(); });

        // Bound concurrent per-HOG build scratch to a share of --flush so pass2 RSS (bucket +
        // groups + builds) stays within budget. est ≈ hog spill bytes × FACTOR (flat_inv +
        // inv_raw + cbuf vs the ~5 B/obs spill record); a lone worker always proceeds.
        const size_t build_budget = std::max<size_t>(budget / 2, size_t(64) << 20);
        constexpr size_t MS_BUILD_FACTOR = 6;
        size_t build_in_use = 0;
        std::mutex build_mtx;
        std::condition_variable build_cv;

        std::atomic<size_t> gnext{0};
        auto p2 = [&]() {
            ShardScratch sc(n_accessions);
            try {
                for (;;) {
                    size_t gi = gnext.fetch_add(1, std::memory_order_relaxed);
                    if (gi >= hkeys.size() || failed.load(std::memory_order_relaxed)) break;
                    uint32_t hidx = hkeys[gi];
                    if (++sc.epoch == 0) { std::fill(sc.seen_epoch.begin(), sc.seen_epoch.end(), 0); sc.epoch = 1; }
                    uint32_t epoch = sc.epoch, n_accs = 0, hog_length = 0;
                    sc.seen_accs.clear(); sc.flat_inv.clear();
                    sc.acc_intervals_map.clear(); sc.acc_pident_map.clear();
                    size_t est = std::max<size_t>(group_bytes[hidx] * MS_BUILD_FACTOR, size_t(1) << 20);
                    {
                        std::unique_lock<std::mutex> lk(build_mtx);
                        build_cv.wait(lk, [&]{ return build_in_use == 0 ||
                                                      build_in_use + est <= build_budget ||
                                                      failed.load(std::memory_order_relaxed); });
                        build_in_use += est;
                    }
                    for (uint64_t ro : groups[hidx]) {
                        const uint8_t* q = &data[ro];
                        q += 4;                              // skip hog_idx
                        uint32_t acc  = read_u32_le(q); q += 4;
                        uint32_t cnum = read_u32_le(q); q += 4;
                        uint32_t ss   = read_u32_le(q); q += 4;
                        uint32_t se   = read_u32_le(q); q += 4;
                        uint8_t  pu8  = *q++;
                        uint32_t nobs = read_u32_le(q); q += 4;
                        if (sc.seen_epoch[acc] != epoch) { sc.seen_epoch[acc] = epoch; ++n_accs; sc.seen_accs.push_back(acc); }
                        hog_length = std::max(hog_length, se + 1);
                        sc.acc_intervals_map[acc].emplace_back(ss, se);
                        { auto [it2, ins] = sc.acc_pident_map.emplace(acc, pu8); if (!ins && pu8 < it2->second) it2->second = pu8; }
                        for (uint32_t k = 0; k < nobs; ++k) {
                            uint32_t hoff = read_u32_le(q); q += 4;
                            uint8_t cod = *q++;
                            sc.flat_inv.emplace_back(ss + hoff, InvObs{acc, cod, cnum});
                        }
                    }
                    ShardResult r = build_inverted_from_scratch(hog_list[hidx]->first, sc,
                                                                n_accs, hog_length, out_zstd_level);
                    emit(std::move(r));
                    // Release grown per-HOG scratch so nt workers don't each retain a big HOG's
                    // capacity (the slab-retention pattern); small HOGs keep buffers for reuse.
                    auto rel = [](auto& v){ if (v.capacity() > (size_t(32) << 20)) { std::decay_t<decltype(v)>().swap(v); } };
                    rel(sc.flat_inv); rel(sc.inv_raw); rel(sc.hog_cbuf);
                    { std::lock_guard<std::mutex> lk(build_mtx); build_in_use -= est; }
                    build_cv.notify_all();
                }
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (!failed.exchange(true)) err = e.what();
            }
        };
        size_t bt = std::max<size_t>(1, std::min(nt, hkeys.size()));
        uint64_t _b0 = clock_ns();
        std::vector<std::thread> pool;
        for (size_t t = 1; t < bt; ++t) pool.emplace_back(p2);
        p2();
        for (auto& t : pool) t.join();
        _t_build += clock_ns() - _b0;
        if (failed.load()) throw std::runtime_error("spill pass 2: " + err);
        // free this bucket's backing store
        std::vector<uint8_t>().swap(data);
        ftruncate(bfd[b], 0);
    }
    if (do_profile)
        std::fprintf(stderr, "spill-prof: pass1(decode+spill) %.1fs | pass2-read %.1fs | pass2-build+compress %.1fs\n",
                     _t_p1/1e9, _t_read/1e9, _t_build/1e9);
}

inline void merge_batches(const std::vector<std::string>& input_paths,
                          const std::string& out_lhg,
                          const std::string& out_lhgi,
                          const std::string& hog_range_start = "",
                          const std::string& hog_range_end   = "",
                          int n_buckets = 1,
                          int n_threads_override = 0,
                          bool do_profile = false,
                          int hot_threshold = 100,
                          int out_zstd_level = 6,
                          size_t mem_budget = 16ull << 30) {
    // out_zstd_level: ZSTD compression level (default 3)

    // Pin a fixed mmap threshold so the per-HOG build buffers (inv_raw/hog_cbuf, multi-GB
    // for super-HOGs) are always mmap-backed and returned to the OS on free. Without this,
    // glibc's *dynamic* threshold ratchets up after freeing large mmap blocks, routes
    // subsequent big allocations through the sbrk arena, and never returns them — so RSS
    // climbs unbounded even though live memory stays within the build budget.
    mallopt(M_MMAP_THRESHOLD, 4 * 1024 * 1024);
    mallopt(M_TRIM_THRESHOLD, 16 * 1024 * 1024);

    // Pass 1: parallel scan of all inputs into per-file ref lists; avoids single-threaded 12s serial cost.
    size_t scan_threads = std::max<size_t>(1, std::min<size_t>(
        n_threads_override > 0 ? size_t(n_threads_override) : std::thread::hardware_concurrency(),
        input_paths.size()));

    std::vector<std::vector<MergeRef>> per_file_refs(input_paths.size());
    std::vector<std::vector<std::string>> source_accs(input_paths.size());
    std::atomic<size_t> scan_next{0};
    std::atomic<bool>   scan_failed{false};
    std::string         scan_error;
    std::mutex          scan_err_mtx;

    auto scan_worker = [&]() {
        for (;;) {
            size_t fi = scan_next.fetch_add(1, std::memory_order_relaxed);
            if (fi >= input_paths.size() || scan_failed.load(std::memory_order_relaxed)) break;
            try {
                MergeRef::Kind kind = detect_input_kind(input_paths[fi]);
                if (kind == MergeRef::Kind::LHB) {
                    std::string cur_acc;
                    bool ok = scan_batch_file(input_paths[fi], fi,
                        [&](const std::string& a) { cur_acc = a; },
                        [&](BatchBlockRef br) {
                            MergeRef r;
                            r.kind = MergeRef::Kind::LHB;
                            r.hog_id = std::move(br.hog_id);
                            r.lhb_acc_id = br.acc_id;
                            r.source_file_idx = fi;
                            r.lhb_shard_hdr_offset = uint64_t(br.shard_hdr_offset);
                            per_file_refs[fi].push_back(std::move(r));
                        });
                    if (!ok) throw std::runtime_error("cannot open batch: " + input_paths[fi]);
                } else {
                    GlobalIndex idx;
                    if (!idx.load_from_lhg(input_paths[fi]))
                        throw std::runtime_error("cannot load .lhg index: " + input_paths[fi]);
                    source_accs[fi] = std::move(idx.accessions);
                    for (const auto& e : idx.entries) {
                        MergeRef r;
                        r.kind = MergeRef::Kind::LHG;
                        r.hog_id = e.hog_id;
                        r.source_file_idx = fi;
                        r.lhg_data_offset = e.data_offset;
                        r.lhg_data_length = e.data_length;
                        per_file_refs[fi].push_back(std::move(r));
                    }
                }
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lk(scan_err_mtx);
                if (!scan_failed.exchange(true)) scan_error = e.what();
            }
            if (input_paths.size() > 100 && (fi + 1) % 1000 == 0)
                std::fprintf(stderr, "  scanned %zu/%zu inputs\r", fi + 1, input_paths.size());
        }
    };

    {
        std::vector<std::thread> st;
        st.reserve(scan_threads - 1);
        for (size_t i = 1; i < scan_threads; ++i) st.emplace_back(scan_worker);
        scan_worker();
        for (auto& t : st) t.join();
    }
    if (scan_failed.load()) throw std::runtime_error("scan failed: " + scan_error);
    if (input_paths.size() > 100) std::cerr << "\n";

    // K-way merge: LHBs are pre-sorted by hog_id string (convert.hpp:215), so each
    // per_file_refs[fi] is already sorted. O(N log K) replaces O(N log N) sort.
    size_t total_refs = 0;
    for (auto& v : per_file_refs) total_refs += v.size();
    std::vector<MergeRef> refs;
    refs.reserve(total_refs);
    {
        struct KM { std::string hid; size_t fi, idx; };
        auto cmp = [](const KM& a, const KM& b) {
            int c = a.hid.compare(b.hid);
            if (c != 0) return c > 0;
            return a.fi > b.fi;
        };
        std::priority_queue<KM, std::vector<KM>, decltype(cmp)> pq(cmp);
        for (size_t fi = 0; fi < per_file_refs.size(); ++fi)
            if (!per_file_refs[fi].empty())
                pq.push({per_file_refs[fi][0].hog_id, fi, 0});
        while (!pq.empty()) {
            size_t fi  = pq.top().fi;
            size_t idx = pq.top().idx;
            pq.pop();
            refs.push_back(std::move(per_file_refs[fi][idx]));
            size_t nxt = idx + 1;
            if (nxt < per_file_refs[fi].size())
                pq.push({per_file_refs[fi][nxt].hog_id, fi, nxt});
        }
        for (auto& v : per_file_refs) { v.clear(); v.shrink_to_fit(); }
    }
    std::cerr << "merge: " << refs.size() << " HOG-blocks from "
              << input_paths.size() << " inputs\n";

    // Global accession registry.
    std::set<std::string> uniq;
    for (const auto& r : refs)
        if (r.kind == MergeRef::Kind::LHB) uniq.insert(r.lhb_acc_id);
    for (const auto& accs : source_accs)
        for (const auto& a : accs) uniq.insert(a);
    std::vector<std::string> accessions(uniq.begin(), uniq.end());

    std::unordered_map<std::string, uint32_t> acc_id_to_idx;
    acc_id_to_idx.reserve(accessions.size() * 2);
    for (uint32_t i = 0; i < accessions.size(); ++i) acc_id_to_idx[accessions[i]] = i;

    auto global_acc = [&](const std::string& acc_id) -> uint32_t {
        auto it = acc_id_to_idx.find(acc_id);
        if (it == acc_id_to_idx.end())
            throw std::runtime_error("acc_id not in registry: " + acc_id);
        return it->second;
    };

    // Remap tables: one per source, computed once, shared read-only across threads.
    std::vector<std::vector<uint32_t>> src_remap(input_paths.size());
    for (size_t fi = 0; fi < input_paths.size(); ++fi) {
        const auto& accs = source_accs[fi];
        src_remap[fi].resize(accs.size());
        for (size_t i = 0; i < accs.size(); ++i) src_remap[fi][i] = global_acc(accs[i]);
    }

    // One fd per source, opened lazily, kept alive for all of pass-2.
    struct SrcFile {
        UniqueFd fd{-1};
        std::vector<uint8_t> buf;  // entire file in RAM for small files (LHBs); avoids pread per HOG
        SrcFile() = default;
        SrcFile(const SrcFile&) = delete;
        SrcFile& operator=(const SrcFile&) = delete;
        SrcFile(SrcFile&&) = delete;
    };
    std::vector<SrcFile> src_files(input_paths.size());
    std::vector<std::atomic<bool>> fd_opened(input_paths.size());
    std::vector<std::mutex> fd_open_mtx(input_paths.size());
    auto get_src_file = [&](size_t fi) -> SrcFile& {
        if (!fd_opened[fi].load(std::memory_order_acquire)) {
            std::lock_guard<std::mutex> lk(fd_open_mtx[fi]);
            if (!fd_opened[fi].load(std::memory_order_relaxed)) {
                int fd = open(input_paths[fi].c_str(), O_RDONLY);
                if (fd < 0) throw std::runtime_error("cannot open: " + input_paths[fi]);
                posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
                src_files[fi].fd = UniqueFd(fd);
                // Buffer small files (LHBs) entirely in RAM — eliminates 60M pread syscalls.
                struct stat st; off_t fsz = 0;
                if (fstat(fd, &st) == 0) fsz = st.st_size;
                if (fsz > 0 && fsz <= off_t(64 * 1024 * 1024)) {
                    src_files[fi].buf.resize(size_t(fsz));
                    ssize_t nr = ::read(fd, src_files[fi].buf.data(), size_t(fsz));
                    if (nr != ssize_t(fsz))
                        throw std::runtime_error("buffered read failed: " + input_paths[fi]);
                }
                fd_opened[fi].store(true, std::memory_order_release);
            }
        }
        return src_files[fi];
    };

    // Pre-scan contiguous HOG groups (refs are stable_sorted by hog_id).
    std::vector<std::pair<size_t, size_t>> groups;
    {
        size_t gs = 0;
        while (gs < refs.size()) {
            size_t ge = gs + 1;
            while (ge < refs.size() && refs[ge].hog_id == refs[gs].hog_id) ++ge;
            const std::string& hid = refs[gs].hog_id;
            bool skip = (!hog_range_start.empty() && hid < hog_range_start) ||
                        (!hog_range_end.empty()   && hid >= hog_range_end);
            if (!skip) groups.emplace_back(gs, ge);
            gs = ge;
        }
    }

    // Bucket assignment: deterministic hash of hog_id.
    if (n_buckets < 1) n_buckets = 1;
    size_t n_buckets_sz = size_t(n_buckets);
    size_t n_groups = groups.size();
    auto bucket_of = [&](const std::string& hog_id) -> size_t {
        uint64_t h = 1469598103934665603ull;  // FNV-1a 64-bit
        for (unsigned char c : hog_id) { h ^= c; h *= 1099511628211ull; }
        return size_t(h % n_buckets_sz);
    };
    std::vector<size_t> group_bucket(n_groups);
    for (size_t g = 0; g < n_groups; ++g)
        group_bucket[g] = bucket_of(refs[groups[g].first].hog_id);

    // Largest-first dispatch.
    std::vector<uint64_t> group_sz(n_groups, 0);
    for (size_t g = 0; g < n_groups; ++g)
        for (size_t i = groups[g].first; i < groups[g].second; ++i)
            group_sz[g] += refs[i].lhg_data_length;
    std::vector<size_t> dispatch_order(n_groups);
    std::iota(dispatch_order.begin(), dispatch_order.end(), 0);
    std::sort(dispatch_order.begin(), dispatch_order.end(),
              [&](size_t a, size_t b) { return group_sz[a] > group_sz[b]; });

    struct ProfCounters {
        std::atomic<uint64_t> ns_decode{0}, ns_build{0}, ns_compress{0}, n_groups{0};
    } prof;

    // Hot-HOG parallel decode state.
    struct HotDecodeTask { size_t g; size_t bi; };
    struct HotHogState {
        std::vector<InvBlock> blocks;
        std::atomic<int>      pending{0};
    };
    std::vector<std::unique_ptr<HotHogState>> hot_states(groups.size());

    std::deque<HotDecodeTask> dq_decode;
    std::mutex                dq_mtx;

    struct GroupResult {
        std::string          hog_id;
        std::vector<uint8_t> payload;
        uint32_t             stored_sz = 0;
        uint32_t             n_accs = 0;
        bool                 is_v9 = false;   // streaming build emits the v9 (uid+trailer) layout
    };
    std::vector<GroupResult> results(groups.size());

    struct BucketOut {
        UniqueFd               lhg_fd{-1};
        UniqueFd               lhgi_fd{-1};
        std::string            lhg_path;
        std::string            lhgi_path;
        uint64_t               write_pos = LHG_HEADER_SZ;
        std::vector<HogIndexEntry> index_entries;
        WriteBuffer            wbuf;
    };
    auto fmt3 = [](size_t b) {
        char s[16]; std::snprintf(s, sizeof(s), "%03zu", b); return std::string(s);
    };
    std::vector<BucketOut> buckets(n_buckets_sz);
    for (size_t b = 0; b < n_buckets_sz; ++b) {
        BucketOut& bk = buckets[b];
        bk.lhg_path  = (n_buckets_sz == 1) ? out_lhg  : out_lhg  + ".b" + fmt3(b);
        bk.lhgi_path = (n_buckets_sz == 1) ? out_lhgi : out_lhgi + ".b" + fmt3(b);
        bk.lhg_fd  = UniqueFd(open(bk.lhg_path.c_str(),  O_WRONLY | O_CREAT | O_TRUNC, 0644));
        bk.lhgi_fd = UniqueFd(open(bk.lhgi_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644));
        if (bk.lhg_fd  < 0) throw std::runtime_error("cannot create: " + bk.lhg_path);
        if (bk.lhgi_fd < 0) throw std::runtime_error("cannot create: " + bk.lhgi_path);
        uint8_t fhdr[LHG_HEADER_SZ] = {};
        memcpy(fhdr, LHG_FILE_MAGIC, 4); fhdr[4] = LHG_VERSION;
        if (::write(bk.lhg_fd, fhdr, LHG_HEADER_SZ) != ssize_t(LHG_HEADER_SZ))
            throw std::runtime_error("write header failed: " + bk.lhg_path);
    }

    std::atomic<bool> failed{false};
    std::string first_error;
    std::mutex err_mtx;

    // Completion-order, memory-bounded output: workers push finished group ids; the writer
    // drains them in COMPLETION order (not index order) and frees each payload immediately.
    // A byte budget throttles workers when too much output is pending unwritten, so peak RAM
    // is bounded regardless of the final .lhg size (vs the old all-results-resident model).
    std::deque<size_t>      done_q;
    std::mutex              done_mtx;
    std::condition_variable done_cv;     // writer wakeup
    std::condition_variable space_cv;    // worker backpressure release
    size_t                  unwritten_bytes = 0;

    // --flush is the TOTAL RSS budget (like ingest/merge-shard). Two consumers share it:
    // unwritten OUTPUT payloads (out_budget) and per-HOG BUILD scratch — decoded input
    // blocks + serialized inv_raw + compress buffer — which dominates (it was ~560 GB,
    // 32 workers × multi-GB, while output stayed ~24 GB). Bound concurrent builds so their
    // summed estimated cost stays under build_budget. est = compressed block bytes ×
    // MEM_FACTOR (decoded ~4x + inv_raw ~4x + cbuf). A lone worker always proceeds, so a
    // single HOG larger than the budget runs (the effective --flush floor is the biggest
    // HOG's scratch). ceiling: a HOG whose scratch exceeds RAM still OOMs; upgrade: stream
    // the per-HOG serialize/compress instead of materializing inv_raw whole.
    // est = compressed bytes × MEM_FACTOR is only a herd-limiter (it prevents 32 workers
    // admitting at once); the real ceiling is a hard ACTUAL-RSS guard, because per-HOG
    // scratch is dominated by the *decompressed* inv_raw (5-8× the compressed input) and
    // no static factor predicts it reliably. A lone worker always proceeds.
    const size_t out_budget   = std::max<size_t>(mem_budget / 4, size_t(64) << 20);
    const size_t build_budget = mem_budget > out_budget ? mem_budget - out_budget : mem_budget;
    const size_t rss_cap      = mem_budget;
    // Streaming build peak ≈ Σ decompressed inputs + columns ≈ compressed bytes × ~8.
    static constexpr size_t MEM_FACTOR = 8;
    size_t                  build_in_use = 0;
    std::mutex              build_mtx;
    std::condition_variable build_cv;

    // Per-worker reusable containers — cleared (not reallocated) per group.
    struct WorkerScratch {
        // Distinct-acc counter: seen_epoch[gacc]==epoch ⇒ already counted this group.
        std::vector<uint32_t>   seen_epoch;
        uint32_t                epoch = 0;
        uint32_t                n_accs = 0;
        std::vector<uint32_t>   seen_accs;   // distinct gaccs this group; cleared, not realloc'd
        std::vector<InvPosition> positions;
        // all_lhg path
        std::vector<InvBlock>   blocks;
        std::vector<std::vector<InvPosition>> src_pos;
        std::vector<InvObs>     obs;
        std::vector<std::vector<InvObs>> per_block_obs;
        // mixed/lhb path — flat (pos, obs) vector, sorted in one pass after all blocks
        std::vector<std::pair<uint32_t, InvObs>> flat_inv;
        std::vector<VarNTRecord> recs;
        std::vector<uint32_t>   contig_cnum;
        // serialize scratch: cleared per group, reused across groups to avoid allocs
        std::vector<uint8_t> ser_hdr_buf, ser_acc_buf, ser_cnum_buf, ser_codon_buf;
        // counting sort scratch for grouping flat_inv by hog_pos (O(N + max_pos) vs O(N log N))
        std::vector<uint32_t> cs_count;   // cs_count[pos] = obs count at that position
        std::vector<uint32_t> cs_pidx;    // cs_pidx[pos]  = positions[] index for that position
        // v8: per-HOG coverage/pident accumulators keyed by global acc_idx
        uint32_t hog_length = 0;
        // .lhb sources: store raw HSP intervals; merged at serialize time to avoid double-counting overlaps
        std::unordered_map<uint32_t, std::vector<std::pair<uint32_t,uint32_t>>> acc_intervals_map;
        // .lhg sources: pre-merged covered_aa summed directly (shards are HOG-disjoint, no overlap)
        std::unordered_map<uint32_t, uint32_t> acc_covered_aa_map;
        std::unordered_map<uint32_t, uint8_t>  acc_pident_map;      // gacc → min pident_u8
    };

    // decode_one_block: pread + decompress one source block into hot_states[g]->blocks[bi].
    auto decode_one_block = [&](size_t g, size_t bi) {
        size_t group_start = groups[g].first;
        const std::string& hog_id = refs[group_start].hog_id;
        const MergeRef& ref = refs[group_start + bi];
        SrcFile& sf = get_src_file(ref.source_file_idx);
        HotHogState& hs = *hot_states[g];

        HogIndexEntry e;
        e.hog_id      = ref.hog_id;
        e.data_offset = ref.lhg_data_offset;
        e.data_length = ref.lhg_data_length;
        if (!read_hog_inverted_fd(sf.fd, input_paths[ref.source_file_idx], e, hs.blocks[bi]))
            throw std::runtime_error("failed to read .lhg HOG " + hog_id + " block " + std::to_string(bi));
        hs.pending.fetch_sub(1, std::memory_order_release);
    };

    auto process_group = [&](size_t g, ZSTD_CCtx* cctx, ZSTD_DCtx* dctx,
                             std::vector<uint8_t>& cbuf_in,
                             std::vector<uint8_t>& raw_block,
                             std::vector<uint8_t>& inv_raw,
                             std::vector<uint8_t>& hog_cbuf,
                             WorkerScratch& sc,
                             uint64_t& tl_decode, uint64_t& tl_build, uint64_t& tl_compress) {
        size_t group_start = groups[g].first;
        size_t group_end   = groups[g].second;
        const std::string& hog_id = refs[group_start].hog_id;
        size_t n_blocks = group_end - group_start;

        // Distinct global accessions contributing to this HOG.
        if (sc.seen_epoch.size() < accessions.size()) sc.seen_epoch.assign(accessions.size(), 0);
        ++sc.epoch;
        sc.n_accs = 0;
        sc.seen_accs.clear();
        auto mark_acc = [&sc](uint32_t gacc) {
            if (sc.seen_epoch[gacc] != sc.epoch) {
                sc.seen_epoch[gacc] = sc.epoch;
                ++sc.n_accs;
                sc.seen_accs.push_back(gacc);
            }
        };

        bool all_lhg = true;
        for (size_t bi = 0; bi < n_blocks && all_lhg; ++bi)
            if (refs[group_start + bi].kind != MergeRef::Kind::LHG) all_lhg = false;

        auto& positions = sc.positions; positions.clear();

        uint64_t t0 = clock_ns();
        uint64_t t_dec_end = 0, t_build_end = 0;

        // Big all-.lhg HOGs take the memory-bounded streaming merge (no materialized obs);
        // stream_blocks then points at the decoded blocks the cursors walk at finalize.
        bool streamed_big = false;
        std::vector<InvBlock>* stream_blocks = nullptr;

        if (all_lhg) {
            // Shared tail for both decode strategies: `blocks` holds the n_blocks decoded
            // InvBlocks; accumulate coverage/pident, remap to global acc space, then K-way
            // heap-merge the per-block position streams into `positions`.
            auto merge_lhg_blocks = [&](std::vector<InvBlock>& blocks) {
                auto& src_pos = sc.src_pos; src_pos.clear(); src_pos.resize(n_blocks);
                sc.hog_length = 0; sc.acc_intervals_map.clear(); sc.acc_covered_aa_map.clear(); sc.acc_pident_map.clear();

                for (size_t bi = 0; bi < n_blocks; ++bi) {
                    const auto& remap = src_remap[refs[group_start + bi].source_file_idx];
                    const InvBlock& blk = blocks[bi];

                    sc.hog_length = std::max(sc.hog_length, blk.hog_length);
                    for (size_t li = 0; li < blk.local_accs.size(); ++li) {
                        if (li >= remap.size()) continue;
                        uint32_t gacc = remap[li];
                        if (li < blk.covered_aa.size())
                            sc.acc_covered_aa_map[gacc] += blk.covered_aa[li];
                        if (li < blk.local_acc_pident.size()) {
                            uint8_t pid = blk.local_acc_pident[li];
                            auto it = sc.acc_pident_map.find(gacc);
                            if (it == sc.acc_pident_map.end() || pid < it->second) sc.acc_pident_map[gacc] = pid;
                        }
                    }

                    auto blk_positions = decode_block(blk);
                    src_pos[bi].reserve(blk_positions.size());
                    for (auto& ip : blk_positions) {
                        InvPosition out_ip;
                        out_ip.hog_pos = ip.hog_pos;
                        out_ip.obs.reserve(ip.obs.size());
                        for (const auto& o : ip.obs) {
                            uint32_t li = o.acc_idx;
                            if (li >= remap.size())
                                throw std::runtime_error("acc_idx out of source registry for HOG " + hog_id);
                            uint32_t gacc = remap[li];
                            mark_acc(gacc);
                            out_ip.obs.push_back({gacc, o.codon_idx, o.cnum});
                        }
                        src_pos[bi].push_back(std::move(out_ip));
                    }
                }

                // K-way min-heap merge of the per-block position streams.
                struct HeapItem { uint32_t hog_pos; size_t src; size_t cur; };
                auto cmp = [](const HeapItem& a, const HeapItem& b) { return a.hog_pos > b.hog_pos; };
                std::priority_queue<HeapItem, std::vector<HeapItem>, decltype(cmp)> heap(cmp);
                for (size_t bi = 0; bi < n_blocks; ++bi)
                    if (!src_pos[bi].empty())
                        heap.push({src_pos[bi][0].hog_pos, bi, 0});
                sc.per_block_obs.resize(n_blocks);
                while (!heap.empty()) {
                    uint32_t cur_pos = heap.top().hog_pos;
                    InvPosition out_pos;
                    out_pos.hog_pos = cur_pos;
                    while (!heap.empty() && heap.top().hog_pos == cur_pos) {
                        HeapItem it = heap.top(); heap.pop();
                        auto& obs = src_pos[it.src][it.cur].obs;
                        sc.per_block_obs[it.src].insert(sc.per_block_obs[it.src].end(),
                                                        std::make_move_iterator(obs.begin()),
                                                        std::make_move_iterator(obs.end()));
                        size_t next = it.cur + 1;
                        if (next < src_pos[it.src].size())
                            heap.push({src_pos[it.src][next].hog_pos, it.src, next});
                    }
                    { size_t tot = 0; for (size_t bi = 0; bi < n_blocks; ++bi) tot += sc.per_block_obs[bi].size(); out_pos.obs.reserve(tot); }
                    for (size_t bi = 0; bi < n_blocks; ++bi) {
                        out_pos.obs.insert(out_pos.obs.end(),
                                           std::make_move_iterator(sc.per_block_obs[bi].begin()),
                                           std::make_move_iterator(sc.per_block_obs[bi].end()));
                        sc.per_block_obs[bi].clear();
                    }
                    positions.push_back(std::move(out_pos));
                }
            };

            // Big-HOG path: build only the acc-union + coverage/pident from headers (no obs
            // decode), so build_inverted_streamed_lhg can window-merge the blocks at finalize.
            auto build_acc_maps_only = [&](std::vector<InvBlock>& blocks) {
                sc.hog_length = 0; sc.acc_intervals_map.clear();
                sc.acc_covered_aa_map.clear(); sc.acc_pident_map.clear();
                for (size_t bi = 0; bi < n_blocks; ++bi) {
                    const auto& remap = src_remap[refs[group_start + bi].source_file_idx];
                    const InvBlock& blk = blocks[bi];
                    sc.hog_length = std::max(sc.hog_length, blk.hog_length);
                    for (size_t li = 0; li < blk.local_accs.size(); ++li) {
                        if (li >= remap.size()) continue;
                        uint32_t gacc = remap[li];
                        if (li < blk.covered_aa.size()) sc.acc_covered_aa_map[gacc] += blk.covered_aa[li];
                        if (li < blk.local_acc_pident.size()) {
                            uint8_t pid = blk.local_acc_pident[li];
                            auto it = sc.acc_pident_map.find(gacc);
                            if (it == sc.acc_pident_map.end() || pid < it->second) sc.acc_pident_map[gacc] = pid;
                        }
                        mark_acc(gacc);
                    }
                }
            };
            // Dispatch: >64 MiB decompressed input → streaming merge; else the in-RAM path.
            auto process_blocks = [&](std::vector<InvBlock>& blocks) {
                size_t raw_tot = 0; for (auto& b : blocks) raw_tot += b.raw.size();
                if (raw_tot > (64u << 20)) { build_acc_maps_only(blocks); streamed_big = true; stream_blocks = &blocks; }
                else                       { merge_lhg_blocks(blocks); }
            };

            bool is_hot = (int(n_blocks) > hot_threshold);
            if (is_hot) {
                // Hot path: fan the per-block pread+decompress out to the worker pool.
                hot_states[g] = std::make_unique<HotHogState>();
                hot_states[g]->blocks.resize(n_blocks);
                hot_states[g]->pending.store(int(n_blocks), std::memory_order_relaxed);
                {
                    std::lock_guard<std::mutex> lk(dq_mtx);
                    for (size_t bi = 0; bi < n_blocks; ++bi)
                        dq_decode.push_back({g, bi});
                }
                while (hot_states[g]->pending.load(std::memory_order_acquire) > 0) {
                    HotDecodeTask dt{g, 0};
                    bool got = false;
                    {
                        std::lock_guard<std::mutex> lk(dq_mtx);
                        if (!dq_decode.empty()) { dt = dq_decode.front(); dq_decode.pop_front(); got = true; }
                    }
                    if (got) decode_one_block(dt.g, dt.bi);
                    else     std::this_thread::yield();
                }
                t_dec_end = clock_ns();
                process_blocks(hot_states[g]->blocks);
            } else {
                // Cold path: decode all blocks inline on this thread.
                auto& blocks = sc.blocks; blocks.clear(); blocks.resize(n_blocks);
                for (size_t bi = 0; bi < n_blocks; ++bi) {
                    const MergeRef& ref = refs[group_start + bi];
                    SrcFile& sf = get_src_file(ref.source_file_idx);
                    HogIndexEntry e;
                    e.hog_id = ref.hog_id;
                    e.data_offset = ref.lhg_data_offset;
                    e.data_length = ref.lhg_data_length;
                    if (!read_hog_inverted_fd(sf.fd, input_paths[ref.source_file_idx], e, blocks[bi]))
                        throw std::runtime_error("failed to read .lhg HOG " + hog_id);
                }
                t_dec_end = clock_ns();
                process_blocks(blocks);
            }
        } else {
            // Mixed/LHB groups: accumulate into a flat (pos, obs) vector, sort once.
            auto& flat_inv = sc.flat_inv; flat_inv.clear();
            sc.hog_length = 0; sc.acc_intervals_map.clear(); sc.acc_covered_aa_map.clear(); sc.acc_pident_map.clear();

            for (size_t bi = 0; bi < n_blocks; ++bi) {
                const MergeRef& ref = refs[group_start + bi];
                SrcFile& sf = get_src_file(ref.source_file_idx);
                int src_fd = sf.fd;

                if (ref.kind == MergeRef::Kind::LHB) {
                    uint32_t acc_idx = global_acc(ref.lhb_acc_id);
                    const auto& sbuf = sf.buf;
                    const size_t hdr_off = ref.lhb_shard_hdr_offset;

                    uint8_t hdr_bytes[28];
                    if (!sbuf.empty()) {
                        if (hdr_off + sizeof(hdr_bytes) > sbuf.size())
                            throw std::runtime_error("buf hdr OOB: " + input_paths[ref.source_file_idx]);
                        memcpy(hdr_bytes, sbuf.data() + hdr_off, sizeof(hdr_bytes));
                    } else {
                        ssize_t nr = ::pread(src_fd, hdr_bytes, sizeof(hdr_bytes), off_t(hdr_off));
                        if (nr != ssize_t(sizeof(hdr_bytes)))
                            throw std::runtime_error("pread header failed: " + input_paths[ref.source_file_idx]);
                    }
                    if (hdr_bytes[4] != SHARD_BLOCK_VERSION)
                        throw std::runtime_error("unsupported shard version in: " + input_paths[ref.source_file_idx]);
                    uint32_t compressed_sz = read_u32_le(hdr_bytes + 8);
                    uint32_t raw_sz        = read_u32_le(hdr_bytes + 12);

                    if (compressed_sz > 1024u * 1024 * 1024)
                        throw std::runtime_error("compressed_sz OOB for HOG " + hog_id);
                    if (raw_sz > 4096u * 1024 * 1024ull)
                        throw std::runtime_error("raw_sz OOB for HOG " + hog_id);

                    const size_t pay_off = hdr_off + sizeof(hdr_bytes);
                    const uint8_t* cdata;
                    if (!sbuf.empty()) {
                        if (pay_off + compressed_sz > sbuf.size())
                            throw std::runtime_error("buf payload OOB: " + input_paths[ref.source_file_idx]);
                        cdata = sbuf.data() + pay_off;
                    } else {
                        cbuf_in.resize(compressed_sz);
                        ssize_t nr = ::pread(src_fd, cbuf_in.data(), compressed_sz, off_t(pay_off));
                        if (nr != ssize_t(compressed_sz))
                            throw std::runtime_error("pread payload failed: " + input_paths[ref.source_file_idx]);
                        cdata = cbuf_in.data();
                    }

                    raw_block.resize(raw_sz);
                    size_t rz = ZSTD_decompressDCtx(dctx, raw_block.data(), raw_sz,
                                                    cdata, compressed_sz);
                    if (ZSTD_isError(rz))
                        throw std::runtime_error(std::string("zstd decompress: ") + ZSTD_getErrorName(rz));

                    // Parse local contig dict
                    const uint8_t* p   = raw_block.data();
                    const uint8_t* end = p + rz;
                    uint32_t n_contigs = 0;
                    int n = read_varint(p, end, &n_contigs);
                    if (!n) throw std::runtime_error("corrupt contig dict in block for HOG " + hog_id);
                    p += n;
                    if (n_contigs > 65536)
                        throw std::runtime_error("n_contigs OOB for HOG " + hog_id);
                    auto& contig_cnum = sc.contig_cnum; contig_cnum.clear(); contig_cnum.resize(n_contigs);
                    for (uint32_t j = 0; j < n_contigs; ++j) {
                        uint32_t len = 0;
                        n = read_varint(p, end, &len);
                        if (!n || p + n + len > end) throw std::runtime_error("truncated contig dict");
                        p += n;
                        std::string_view uid(reinterpret_cast<const char*>(p), len);
                        p += len;
                        auto us = uid.rfind('_');
                        if (us == std::string_view::npos)
                            throw std::runtime_error("unitig_id has no '_': " + std::string(uid));
                        uint32_t cnum = 0;
                        auto res = std::from_chars(uid.data() + us + 1, uid.data() + uid.size(), cnum);
                        if (res.ec != std::errc{})
                            throw std::runtime_error("unitig_id non-numeric suffix: " + std::string(uid));
                        contig_cnum[j] = cnum;
                    }

                    auto& recs = sc.recs; recs.clear();
                    if (!deserialize_varnt_block(p, end, recs))
                        throw std::runtime_error("corrupt varnt block for HOG " + hog_id);

                    mark_acc(acc_idx);
                    for (const auto& r : recs) {
                        if (r.contig_idx >= contig_cnum.size())
                            throw std::runtime_error("contig_idx OOB in HOG " + hog_id);
                        uint32_t cnum = contig_cnum[r.contig_idx];
                        uint8_t pu8 = uint8_t(std::min(100.0f, r.pident + 0.5f));
                        sc.hog_length = std::max(sc.hog_length, r.send + 1);
                        sc.acc_intervals_map[acc_idx].emplace_back(r.sstart, r.send);
                        {
                            auto [it2, ins] = sc.acc_pident_map.emplace(acc_idx, pu8);
                            if (!ins && pu8 < it2->second) it2->second = pu8;
                        }
                        for (const auto& o : r.vars) {
                            uint32_t hog_pos   = r.sstart + o.hog_offset;
                            uint8_t  codon_idx = uint8_t(o.packed_codon >> 2);
                            flat_inv.emplace_back(hog_pos, InvObs{acc_idx, codon_idx, cnum});
                        }
                    }
                } else {
                    // .lhg source: decode and re-invert into global acc space.
                    HogIndexEntry e;
                    e.hog_id = ref.hog_id;
                    e.data_offset = ref.lhg_data_offset;
                    e.data_length = ref.lhg_data_length;
                    InvBlock blk;
                    if (!read_hog_inverted_fd(src_fd, input_paths[ref.source_file_idx], e, blk))
                        throw std::runtime_error("failed to read .lhg HOG " + hog_id);

                    const auto& remap = src_remap[ref.source_file_idx];
                    sc.hog_length = std::max(sc.hog_length, blk.hog_length);
                    for (size_t li = 0; li < blk.local_accs.size(); ++li) {
                        if (li >= remap.size()) continue;
                        uint32_t gacc = remap[li];
                        if (li < blk.covered_aa.size())
                            sc.acc_covered_aa_map[gacc] += blk.covered_aa[li];
                        if (li < blk.local_acc_pident.size()) {
                            uint8_t pid = blk.local_acc_pident[li];
                            auto it = sc.acc_pident_map.find(gacc);
                            if (it == sc.acc_pident_map.end() || pid < it->second) sc.acc_pident_map[gacc] = pid;
                        }
                    }

                    auto blk_positions = decode_block(blk);
                    for (const auto& ip : blk_positions) {
                        for (const auto& o : ip.obs) {
                            uint32_t li = o.acc_idx;
                            if (li >= remap.size())
                                throw std::runtime_error("acc_idx out of source registry for HOG " + hog_id);
                            uint32_t gacc = remap[li];
                            mark_acc(gacc);
                            flat_inv.emplace_back(ip.hog_pos, InvObs{gacc, o.codon_idx, o.cnum});
                        }
                    }
                }
            }

            t_dec_end = clock_ns();

            // Counting sort by hog_pos (O(N + hog_length)) — better than O(N log N) flat sort
            // since hog_length ≤ 2000 typically and N = all obs across all blocks.
            uint32_t max_pos = sc.hog_length > 0 ? sc.hog_length - 1 : 0;
            for (const auto& e : flat_inv) if (e.first > max_pos) max_pos = e.first;
            auto& cs_count = sc.cs_count; cs_count.assign(max_pos + 1, 0);
            auto& cs_pidx  = sc.cs_pidx;  cs_pidx.assign(max_pos + 1, UINT32_MAX);
            for (const auto& e : flat_inv) ++cs_count[e.first];
            // Build InvPosition vector (sorted by pos; reserve obs counts)
            for (uint32_t p = 0; p <= max_pos; ++p) {
                if (!cs_count[p]) continue;
                cs_pidx[p] = uint32_t(positions.size());
                InvPosition ip; ip.hog_pos = p; ip.obs.reserve(cs_count[p]);
                positions.push_back(std::move(ip));
            }
            // Fill obs into positions
            for (auto& e : flat_inv)
                positions[cs_pidx[e.first]].obs.push_back(std::move(e.second));
            // Sort each position's obs by acc_idx
            for (auto& pos : positions)
                std::sort(pos.obs.begin(), pos.obs.end(),
                          [](const InvObs& a, const InvObs& b) { return a.acc_idx < b.acc_idx; });
        }

        // Build per-HOG local acc dict from seen_accs (populated by mark_acc above).
        std::vector<uint32_t> local_accs = std::move(sc.seen_accs);
        std::sort(local_accs.begin(), local_accs.end());

        // Build per-local-acc pident and covered_aa vectors from maps.
        // covered_aa = merged .lhb intervals + direct .lhg sum (shards are HOG-disjoint).
        std::vector<uint8_t>  local_acc_pident(local_accs.size(), 0);
        std::vector<uint32_t> covered_aa_v(local_accs.size(), 0);
        for (size_t li = 0; li < local_accs.size(); ++li) {
            uint32_t gacc = local_accs[li];
            auto it = sc.acc_pident_map.find(gacc);
            if (it != sc.acc_pident_map.end()) local_acc_pident[li] = it->second;
            auto it2 = sc.acc_covered_aa_map.find(gacc);
            if (it2 != sc.acc_covered_aa_map.end()) covered_aa_v[li] += it2->second;
            auto it3 = sc.acc_intervals_map.find(gacc);
            if (it3 != sc.acc_intervals_map.end()) covered_aa_v[li] += merged_interval_cov(it3->second);
        }

        // Serialize + compress the inverted block (shared with the spill build).
        uint32_t stored_sz; const uint8_t* payload; size_t payload_sz;
        if (streamed_big) {
            // Memory-bounded windowed merge: no materialized obs vectors (frees input raw
            // before compress). Output compressed payload lands in hog_cbuf.
            std::vector<const std::vector<uint32_t>*> remaps; remaps.reserve(n_blocks);
            for (size_t bi = 0; bi < n_blocks; ++bi)
                remaps.push_back(&src_remap[refs[group_start + bi].source_file_idx]);
            build_inverted_streamed_lhg(*stream_blocks, remaps, local_accs, local_acc_pident,
                                        covered_aa_v, sc.hog_length, cctx, hog_cbuf);
            t_build_end = clock_ns();
            stored_sz = 0x40000000u; payload = hog_cbuf.data(); payload_sz = hog_cbuf.size();
        } else {
            finalize_inverted_payload(local_accs, local_acc_pident, positions, sc.hog_length, covered_aa_v,
                                      inv_raw, hog_cbuf, sc.ser_hdr_buf, sc.ser_acc_buf, sc.ser_cnum_buf,
                                      sc.ser_codon_buf, cctx, 0, &t_build_end, stored_sz, payload, payload_sz);
        }
        tl_decode   += t_dec_end - t0;
        tl_build    += t_build_end - t_dec_end;
        tl_compress += clock_ns() - t_build_end;

        GroupResult& gr = results[g];
        gr.hog_id = hog_id;
        gr.n_accs = sc.n_accs;
        gr.stored_sz = stored_sz;
        gr.is_v9 = true;                  // both the streaming build and finalize emit v9
        gr.payload.assign(payload, payload + payload_sz);
        size_t psz = gr.payload.size();
        {
            std::unique_lock<std::mutex> lk(done_mtx);
            space_cv.wait(lk, [&]{ return unwritten_bytes <= out_budget ||
                                          failed.load(std::memory_order_relaxed); });
            done_q.push_back(g);
            unwritten_bytes += psz;
        }
        done_cv.notify_one();
    };

    std::atomic<size_t> next_group{0};
    // (failed/first_error/err_mtx hoisted above process_group)
    size_t hw = std::thread::hardware_concurrency();
    size_t default_threads = std::min<size_t>(32u, hw);
    size_t n_threads = (n_threads_override > 0)
                       ? std::min<size_t>(size_t(n_threads_override), groups.size())
                       : std::min<size_t>(default_threads, groups.size());
    if (n_threads == 0) n_threads = 1;

    auto worker = [&]() {
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, out_zstd_level);
        std::vector<uint8_t> cbuf_in, raw_block, inv_raw, hog_cbuf;
        WorkerScratch scratch;
        uint64_t tl_decode = 0, tl_build = 0, tl_compress = 0;
        for (;;) {
            size_t di = next_group.fetch_add(1, std::memory_order_relaxed);
            if (di < dispatch_order.size() && !failed.load(std::memory_order_relaxed)) {
                size_t g = dispatch_order[di];
                // Reserve build-scratch budget before decoding so concurrent builds stay
                // within --flush; a lone worker proceeds regardless (oversized-HOG floor).
                size_t est = std::max<size_t>(group_sz[g] * MEM_FACTOR, size_t(1) << 20);
                {
                    std::unique_lock<std::mutex> lk(build_mtx);
                    build_cv.wait(lk, [&]{
                        if (failed.load(std::memory_order_relaxed)) return true;
                        if (build_in_use == 0) return true;                 // lone runner → progress
                        return build_in_use + est <= build_budget &&        // herd limit
                               fast_rss_bytes() < rss_cap;                  // hard RSS ceiling
                    });
                    build_in_use += est;
                }
                bool group_ok = true;
                try {
                    process_group(g, cctx, dctx, cbuf_in, raw_block, inv_raw, hog_cbuf,
                                  scratch, tl_decode, tl_build, tl_compress);
                } catch (const std::exception& e) {
                    group_ok = false;
                    { std::lock_guard<std::mutex> lk(err_mtx);
                      if (!failed.exchange(true)) first_error = e.what(); }
                    done_cv.notify_all(); space_cv.notify_all();
                }
                { std::lock_guard<std::mutex> lk(build_mtx); build_in_use -= est; }
                build_cv.notify_all();
                // Return grown scratch to the OS so 32 workers don't retain GBs of capacity
                // after a big HOG (the slab-retention pattern). Small HOGs keep their buffers
                // for reuse; only oversized ones (>64 MiB) are released.
                auto shrink_big = [](std::vector<uint8_t>& v){
                    if (v.capacity() > (size_t(64) << 20)) std::vector<uint8_t>().swap(v); };
                shrink_big(inv_raw); shrink_big(hog_cbuf); shrink_big(cbuf_in); shrink_big(raw_block);
                // After a big HOG, reclaim arena pages (sub-mmap-threshold allocations) so RSS
                // tracks live memory. Gated to large groups to avoid the trim cost on small ones.
                if (est > (size_t(512) << 20)) malloc_trim(0);
                if (!group_ok) break;
                continue;
            }
            HotDecodeTask dt{0, 0};
            bool got = false;
            {
                std::lock_guard<std::mutex> lk(dq_mtx);
                if (!dq_decode.empty()) { dt = dq_decode.front(); dq_decode.pop_front(); got = true; }
            }
            if (got) {
                try { decode_one_block(dt.g, dt.bi); }
                catch (const std::exception& e) {
                    { std::lock_guard<std::mutex> lk(err_mtx);
                      if (!failed.exchange(true)) first_error = e.what(); }
                    done_cv.notify_all(); space_cv.notify_all();
                    break;
                }
                continue;
            }
            bool hog_done = (next_group.load(std::memory_order_relaxed) >= dispatch_order.size());
            bool hot_empty;
            { std::lock_guard<std::mutex> lk(dq_mtx); hot_empty = dq_decode.empty(); }
            if (hog_done && hot_empty) break;
            std::this_thread::yield();
        }
        prof.ns_decode.fetch_add(tl_decode, std::memory_order_relaxed);
        prof.ns_build.fetch_add(tl_build,   std::memory_order_relaxed);
        prof.ns_compress.fetch_add(tl_compress, std::memory_order_relaxed);
        prof.n_groups.fetch_add(1, std::memory_order_relaxed);
        ZSTD_freeCCtx(cctx);
        ZSTD_freeDCtx(dctx);
    };

    std::string writer_error;
    auto writer_fn = [&]() {
        try {
            size_t written = 0;
            while (written < n_groups) {
                size_t g;
                {
                    std::unique_lock<std::mutex> lk(done_mtx);
                    done_cv.wait(lk, [&]{ return !done_q.empty() ||
                                                 failed.load(std::memory_order_relaxed); });
                    if (done_q.empty()) return;   // failed with nothing left to drain
                    g = done_q.front(); done_q.pop_front();
                }
                ++written;
                GroupResult& gr = results[g];
                BucketOut&   bk = buckets[group_bucket[g]];
                uint64_t off = bk.write_pos;
                uint8_t hdr8[8];
                memcpy(hdr8, gr.is_v9 ? LHG_HOG_ENTRY_MAGIC_V2 : LHG_HOG_ENTRY_MAGIC, 4);
                for (int i = 0; i < 4; ++i) hdr8[4+i] = uint8_t(gr.stored_sz >> (8*i));
                bk.wbuf.append(bk.lhg_fd, bk.lhg_path, hdr8, 8);
                bk.write_pos += 8;
                bk.wbuf.append(bk.lhg_fd, bk.lhg_path, gr.payload.data(), gr.payload.size());
                bk.write_pos += gr.payload.size();
                bk.index_entries.push_back(
                    {std::move(gr.hog_id), off, bk.write_pos - off, gr.n_accs});
                size_t psz = gr.payload.size();
                std::vector<uint8_t>().swap(gr.payload);
                { std::lock_guard<std::mutex> lk(done_mtx); unwritten_bytes -= psz; }
                space_cv.notify_all();
            }
            for (auto& bk : buckets) bk.wbuf.flush_to(bk.lhg_fd, bk.lhg_path);
        } catch (const std::exception& e) {
            writer_error = e.what();
            failed.store(true, std::memory_order_relaxed);
            space_cv.notify_all();
        }
    };

    std::thread writer_thread(writer_fn);
    if (groups.size() == 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        pool.reserve(n_threads);
        for (size_t t = 0; t < n_threads; ++t) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
    }
    writer_thread.join();
    if (failed.load()) {
        throw std::runtime_error(writer_error.empty() ? first_error : writer_error);
    }

    if (do_profile) {
        uint64_t d = prof.ns_decode.load();
        uint64_t b = prof.ns_build.load();
        uint64_t c = prof.ns_compress.load();
        uint64_t total = d + b + c;
        if (total == 0) total = 1;
        auto pct = [&](uint64_t x) { return int(x * 100 / total); };
        auto secs = [](uint64_t x) { return double(x) / 1e9; };
        std::fprintf(stderr, "profile[hogs]     n=%zu threads=%zu\n", n_groups, n_threads);
        std::fprintf(stderr, "profile[decode]   %.1fs (%d%%)\n", secs(d), pct(d));
        std::fprintf(stderr, "profile[build]    %.1fs (%d%%)\n", secs(b), pct(b));
        std::fprintf(stderr, "profile[compress] %.1fs (%d%%)\n", secs(c), pct(c));
    }

    // Finalize each bucket: LHGI section into the .lhg, backfill header, standalone .lhgi.
    size_t total_hogs = 0;
    for (size_t b = 0; b < n_buckets_sz; ++b) {
        BucketOut& bk = buckets[b];
        uint64_t index_offset = bk.write_pos;
        auto idx_bytes = build_index_bytes(bk.index_entries, accessions);
        { WriteBuffer wb; wb.append(bk.lhg_fd, bk.lhg_path, idx_bytes.data(), idx_bytes.size()); wb.flush_to(bk.lhg_fd, bk.lhg_path); }

        uint8_t final_hdr[LHG_HEADER_SZ] = {};
        memcpy(final_hdr, LHG_FILE_MAGIC, 4);
        final_hdr[4] = LHG_VERSION;
        for (int i = 0; i < 8; ++i) final_hdr[8+i] = uint8_t(index_offset >> (8*i));
        if (::pwrite(bk.lhg_fd, final_hdr, LHG_HEADER_SZ, 0) != ssize_t(LHG_HEADER_SZ))
            throw std::runtime_error("pwrite header failed: " + bk.lhg_path);

        { WriteBuffer wb; wb.append(bk.lhgi_fd, bk.lhgi_path, idx_bytes.data(), idx_bytes.size()); wb.flush_to(bk.lhgi_fd, bk.lhgi_path); }

        total_hogs += bk.index_entries.size();
        std::cerr << "merge done: " << bk.index_entries.size() << " HOGs, "
                  << accessions.size() << " accessions → " << bk.lhg_path << "\n"
                  << "index:      " << bk.lhgi_path
                  << " (" << idx_bytes.size() / 1024 << " KB)\n";
    }
    if (n_buckets_sz > 1)
        std::cerr << "buckets: " << n_buckets_sz << " files, " << total_hogs << " HOGs total\n";
}

} // namespace lhi
