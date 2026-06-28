#pragma once
#include "convert.hpp"
#include "container.hpp"
#include "partition.hpp"
#include <fstream>
#include <iomanip>
#include <numeric>
#include <malloc.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <deque>
#include <functional>
#include <future>
#include <chrono>
#include <immintrin.h>
#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

namespace lhi {

namespace { constexpr auto k_acgt_lut = []() noexcept {
    std::array<uint8_t,256> t{};
    t['A']=t['C']=t['G']=t['T']=1; return t; }(); }

inline size_t read_rss_bytes() {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            size_t kb = 0;
            sscanf(line.c_str() + 6, "%zu", &kb);
            return kb * 1024;
        }
    }
    return 0;
}

// Compute flush threshold from cgroup / SLURM / MemTotal.
// Returns threshold in bytes; emits one diagnostic line to stderr.
inline size_t compute_flush_threshold(size_t explicit_bytes) {
    if (explicit_bytes != 0) return explicit_bytes;

    auto chomp = [](char* s) {
        size_t n = strlen(s);
        while (n && (s[n-1] == '\n' || s[n-1] == '\r')) s[--n] = '\0';
    };

    auto read_v1 = [&]() -> unsigned long long {
        FILE* fg = fopen("/proc/self/cgroup", "r");
        if (!fg) return 0;
        char cgrel[512] = {}, line[512];
        while (fgets(line, sizeof(line), fg)) {
            char* p1 = strchr(line, ':');  if (!p1) continue;
            char* p2 = strchr(p1+1, ':'); if (!p2) continue;
            *p2 = '\0';
            if (strstr(p1+1, "memory")) { chomp(p2+1); snprintf(cgrel, sizeof(cgrel), "%s", p2+1); break; }
        }
        fclose(fg);
        if (!cgrel[0]) return 0;
        char stat[700];
        snprintf(stat, sizeof(stat), "/sys/fs/cgroup/memory%s/memory.stat", cgrel);
        FILE* fs = fopen(stat, "r");
        if (!fs) return 0;
        unsigned long long v = 0;
        while (fgets(line, sizeof(line), fs))
            if (strncmp(line, "hierarchical_memory_limit ", 26) == 0) { sscanf(line+26, "%llu", &v); break; }
        fclose(fs);
        return (v > 0 && v < (1ull<<62)) ? v : 0;
    };

    auto read_v2 = [&]() -> unsigned long long {
        FILE* fg = fopen("/proc/self/cgroup", "r");
        if (!fg) return 0;
        char cgpath[512] = {}, line[512];
        while (fgets(line, sizeof(line), fg))
            if (line[0]=='0' && line[1]==':' && line[2]==':') { chomp(line); snprintf(cgpath, sizeof(cgpath), "%s", line+3); break; }
        fclose(fg);
        if (!cgpath[0]) return 0;
        char base[640];
        snprintf(base, sizeof(base), "/sys/fs/cgroup%s", cgpath);
        unsigned long long best = 0;
        while (strlen(base) > strlen("/sys/fs/cgroup")) {
            char fpath[700];
            snprintf(fpath, sizeof(fpath), "%s/memory.max", base);
            FILE* f = fopen(fpath, "r");
            if (f) {
                char val[64] = {}; fgets(val, sizeof(val), f); fclose(f);
                if (val[0] != 'm') { unsigned long long v = 0; sscanf(val, "%llu", &v); if (v > 0 && (!best || v < best)) best = v; }
            }
            char* slash = strrchr(base, '/');
            if (!slash || slash == base) break;
            *slash = '\0';
        }
        return best;
    };

    // MemAvailable reflects real system pressure better than MemTotal.
    auto read_available = []() -> unsigned long long {
        FILE* f = fopen("/proc/meminfo", "r");
        if (!f) return 0;
        char line[128];
        while (fgets(line, sizeof(line), f))
            if (strncmp(line, "MemAvailable:", 13) == 0) {
                unsigned long long kb = 0; sscanf(line+13, "%llu", &kb); fclose(f); return kb * 1024;
            }
        fclose(f); return 0;
    };

    unsigned long long limit = read_v1();
    if (!limit) limit = read_v2();
    if (!limit) limit = read_available();
    size_t flush = limit ? size_t(limit * 7 / 10) : 2048ull*1024*1024;
    std::cerr << "ingest: mem avail " << (limit>>20) << " MiB → flush at " << (flush>>20) << " MiB\n";
    return flush;
}

// ── Multithreaded ingest ─────────────────────────────────────────────────────
//
// Reader thread dispatches lines to N workers by round-robin (acc_idx % N).
// Dispatcher owns the global acc registry (no lock needed — sole accessor).
// Each worker writes to its own tN.lhp; partition.idx is merged post-join.
//
// Queue transport: arena-batch — reader copies raw line bytes into a per-batch
// vector<char> arena; WorkItems store (acc_idx, byte-offset, len). Workers
// parse string_views directly from the arena. This reduces cross-thread frees
// from 845M individual string allocs to ~825K arena allocs (1 per batch).

inline void ingest_mt(const std::string& tsv_path, const std::string& out_dir,
                      ConvertOptions opts, size_t N, size_t flush_threshold) {
    std::filesystem::create_directories(out_dir);

    // Limit glibc per-thread arenas to N+1 to cap RSS inflation.
    mallopt(M_ARENA_MAX, int(N) + 1);
    // Pin the mmap threshold (8 MiB). By default glibc RAISES it whenever a large mmap'd
    // block is freed (assuming reuse), which then routes the big per-flush snapshot buffers
    // into the arenas, where freeing them only returns the memory to glibc — not the OS.
    // Setting it explicitly also disables that dynamic growth, so buffers ≥8 MiB stay
    // mmap-backed and munmap straight back to the OS on free. The buffers are per-flush
    // (thousands), not per-record, so the extra mmap/munmap syscalls are negligible.
    mallopt(M_MMAP_THRESHOLD, 8 << 20);

    // ── Global acc registry (dispatcher-owned, no lock needed) ───────────────
    SvDict             acc_map;
    std::vector<std::string> acc_vec;

    // ── Arena-batched bounded queue (one per worker) ─────────────────────────
    static constexpr size_t BATCH_SZ  = 1024;
    static constexpr size_t QUEUE_CAP = 64;

    struct WorkItem { uint32_t acc_idx, offset, len; };
    struct Batch {
        std::vector<char>     arena;
        std::vector<WorkItem> items;
        bool full() const { return items.size() >= BATCH_SZ; }
        void push(uint32_t acc_idx, const char* ptr, size_t n) {
            uint32_t off = uint32_t(arena.size());
            arena.insert(arena.end(), ptr, ptr + n);
            items.push_back({acc_idx, off, uint32_t(n)});
        }
    };


    struct WorkQueue {
        std::deque<Batch>       q;
        std::mutex              mtx;
        std::condition_variable cv_push, cv_pop;
        bool                    done = false;

        void push(Batch b) {
            std::unique_lock lk(mtx);
            cv_push.wait(lk, [&]{ return q.size() < QUEUE_CAP; });
            q.push_back(std::move(b));
            cv_pop.notify_one();
        }
        bool pop(Batch& out) {
            std::unique_lock lk(mtx);
            cv_pop.wait(lk, [&]{ return !q.empty() || done; });
            if (q.empty()) return false;
            out = std::move(q.front());
            q.pop_front();
            cv_push.notify_one();
            return true;
        }
        void finish() {
            std::lock_guard lk(mtx);
            done = true;
            cv_pop.notify_all();
        }
    };

    std::vector<std::unique_ptr<WorkQueue>> queues(N);
    for (auto& q : queues) q = std::make_unique<WorkQueue>();

    // ── Fused spill-bucket count B (shared by all workers; baked into the on-disk
    // bucket each record lands in). Any B≥1 is correct — it only trades pass2 bucket
    // size (memory) against bucket count (read overhead). Size buckets to ~1/3 of the
    // flush budget, estimating total spill from the compressed input.
    const size_t full_budget = (opts.flush_bytes != 0 ? opts.flush_bytes : flush_threshold);
    size_t B;
    if (opts.buckets) B = opts.buckets;
    else {
        size_t in_sz = 0;
        std::error_code ec; auto fs = std::filesystem::file_size(tsv_path, ec);
        if (!ec) in_sz = size_t(fs);
        size_t bucket_budget = std::max<size_t>(full_budget / 3, size_t(64) << 20);
        // Measured: spill ≈ 1.75× a .zst input (0.5× decompressed). Use ×4 to err toward MORE
        // buckets (smaller pass2 `data` → memory-safe); cost of over-bucketing is only extra
        // bucket reads. ceiling: assumes input compressed; a plain .tsv input over-buckets ~8×
        // (still safe). Pass --buckets to size exactly.
        size_t est_spill = std::max<size_t>(in_sz * 4, bucket_budget);
        B = std::min<size_t>(4096, std::max<size_t>(1, (est_spill + bucket_budget - 1) / bucket_budget));
    }

    // ── Per-worker aggregate counters ────────────────────────────────────────
    std::vector<uint64_t> w_written(N,0), w_skipped(N,0), w_obs_dropped(N,0);
    std::vector<uint64_t> w_stall_ns(N,0), w_queue_ns(N,0);

    // ── Worker lambda ─────────────────────────────────────────────────────────
    auto do_worker = [&](size_t tid) {
        // flat slab: (hog_offset:10 | codon6:6) → ceiling: i < 1024
        struct PackedObs { uint16_t v; };
        struct RecHdr {
            uint32_t contig_idx, sstart, send, obs_off;
            float    pident;
            uint16_t obs_n;
            int8_t   qframe;
            uint8_t  _pad;
            double   evalue;  // last: 8-byte aligned at offset 24
        };  // sizeof = 32

        struct LocalAcc {
            SvDict                   contig_dict;
            std::vector<std::string> contig_strings;
            std::vector<uint32_t>   rec_hog;   // SoA key (counting-sort pass only)
            std::vector<RecHdr>     rec_hdr;   // SoA metadata
            std::vector<PackedObs>  obs_pool;  // all vars, 2B each, contiguous
            size_t tracked_obs_cap = 0, tracked_hdr_cap = 0, tracked_hog_cap = 0;
            uint32_t max_hog_idx = 0;
        };
        std::vector<LocalAcc> local_accs;  // indexed by acc_idx / N (dense per tid)
        SvDict hog_dict;
        std::vector<std::string> hog_strings;

        size_t mem_bytes    = 0;
        uint32_t rss_tick   = 0;
        // --flush is the total RSS budget. The flush is async double-buffered (one batch
        // is compressed+written in the background while the next accumulates), so peak
        // RSS ~= 2x the per-worker batch. Trigger at HALF the per-worker share so the
        // in-flight snapshot + new accumulation together stay within budget; a frequent
        // hard RSS guard (below) catches estimate undercount (dicts, arenas, fragmentation).
        // --flush is the hard RSS CAP. The steady accumulation is then compressed at the
        // final flush, a transient that roughly doubles it — so target accumulation at a
        // QUARTER of the cap per worker (cap/2 total, halved again for the async snapshot).
        // The peak (steady + final-flush spike) then lands at/under --flush.
        const size_t flush_budget = (opts.flush_bytes != 0 ? opts.flush_bytes : flush_threshold);
        // The DATA-flush trigger is sized off the FULL budget — DON'T tighten it for headroom.
        // Lowering thr means smaller batches → more flushes → more frames → a BIGGER partition
        // index, which is non-evictable and grows RAM: headroom that defeats itself, plus flush
        // overhead. Keep thr at full budget.
        const size_t thr     = std::max<size_t>(1, flush_budget / (N * 4));
        // Headroom lives ONLY in the hard RSS guard: fire at 95% of --flush so the index + the
        // ~256-record overshoot between checks land strictly under budget, without changing how
        // often data is flushed. (Measured: peak was 24.1 GB at 100% on a 24 GB budget.)
        const size_t rss_cap = flush_budget * 95 / 100;

        // This worker owns B spill-bucket files (tN.bucket.b). The async flush expands records
        // into spill-record bytes, bucketed by hash(hog_name)%B (consistent across workers →
        // every worker's records for a HOG co-locate in the same bucket index), and drains each
        // bucket buffer once it reaches CHUNK.
        std::vector<UniqueFd> bfd;
        std::vector<int>      bfd_ints;
        const size_t CHUNK = std::max<size_t>(size_t(64) << 10,
                                              (flush_budget / 8) / std::max<size_t>(1, N * B));
        bfd.reserve(B); bfd_ints.reserve(B);
        for (size_t b = 0; b < B; ++b) {
            std::string bp = out_dir + "/t" + std::to_string(tid) + ".bucket." + std::to_string(b);
            int fd = ::open(bp.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (fd < 0) throw std::runtime_error("cannot create spill bucket: " + bp);
            bfd.emplace_back(fd); bfd_ints.push_back(fd);
        }

        auto intern = [](SvDict& d, std::vector<std::string>& v, std::string_view s) -> uint32_t {
            auto it = d.find(s);
            if (it != d.end()) return it->second;
            uint32_t idx = uint32_t(v.size());
            v.emplace_back(s);
            d.emplace(v.back(), idx);
            return idx;
        };

        std::future<void>        flush_future;
        uint64_t& n_stall_ns = w_stall_ns[tid];
        uint64_t& n_queue_ns = w_queue_ns[tid];
        AlignedResult ar_buf;

        auto flush = [&]() {
            // Serialise: wait for the previous async flush before evacuating new data.
            if (flush_future.valid()) { auto _t0=std::chrono::steady_clock::now(); flush_future.get(); n_stall_ns+=uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-_t0).count()); }

            struct AccSnap {
                uint32_t                 acc_idx, max_hog_idx;
                std::vector<uint32_t>    rec_hog;
                std::vector<RecHdr>      rec_hdr;
                std::vector<PackedObs>   obs_pool;
                std::vector<std::string> contig_strings;
            };
            std::vector<AccSnap> snaps;

            for (size_t ai = 0; ai < local_accs.size(); ai++) {
                LocalAcc& la = local_accs[ai];
                if (la.rec_hog.empty()) continue;

                AccSnap s;
                s.acc_idx     = uint32_t(ai * N + tid);
                s.max_hog_idx = la.max_hog_idx;

                // Move the slabs out (leaves la's vectors empty, capacity ~0) and RELEASE
                // that capacity instead of re-reserving it. The old reserve() retained each
                // accession's peak slab capacity for the whole run; with many accessions
                // (contiguous in the input) that accumulated unbounded RSS across all
                // ~50k LocalAccs — the real OOM driver, not --flush. Re-growth on reuse is
                // cheap and only matters for interleaved accessions (rare here).
                s.rec_hog        = std::move(la.rec_hog);
                s.rec_hdr        = std::move(la.rec_hdr);
                s.obs_pool       = std::move(la.obs_pool);
                s.contig_strings = std::move(la.contig_strings);

                la.obs_pool.shrink_to_fit();  la.tracked_obs_cap = 0;
                la.rec_hdr.shrink_to_fit();   la.tracked_hdr_cap = 0;
                la.rec_hog.shrink_to_fit();   la.tracked_hog_cap = 0;
                la.contig_dict.clear();
                la.max_hog_idx = 0;
                snaps.push_back(std::move(s));
            }
            mem_bytes = 0;
            if (snaps.empty()) return;

            // Snapshot hog_strings: worker may extend concurrently after we return.
            // All hog indices in snaps are bounded by max_hog_idx < hog_strings.size()
            // at this point, so the snapshot covers every name the async thread needs.
            auto hog_snap = hog_strings;

            // Expand each record into a spill record (the merge-shard pass2 contract) and
            // append to this worker's bucket files, bucketed by hash(hog_name)%B (consistent
            // across workers → a HOG's records co-locate in one bucket index). The record
            // carries the LOCAL hog idx; merge-shard remaps it to global via this worker's
            // hog registry. No sort/columnar/compress — that all moves to the single pass2.
            flush_future = std::async(std::launch::async,
                [snaps = std::move(snaps), hog = std::move(hog_snap),
                 B, CHUNK, bfds = bfd_ints]() mutable {
                    std::vector<std::vector<uint8_t>> lbuf(B);
                    std::vector<uint8_t> cbuf;
                    // Each drain is one self-describing frame: [u32 csz][u32 usz][csz zstd bytes].
                    // Compress (level 3) on the decode-bound worker cores; merge-shard decompresses
                    // on read.
                    auto drain = [&](size_t b) {
                        size_t usz = lbuf[b].size();
                        size_t bound = ZSTD_compressBound(usz);
                        if (cbuf.size() < bound + 8) cbuf.resize(bound + 8);
                        size_t csz = ZSTD_compress(cbuf.data() + 8, bound, lbuf[b].data(), usz, 3);
                        if (ZSTD_isError(csz)) throw std::runtime_error(std::string("spill zstd: ") + ZSTD_getErrorName(csz));
                        for (int i = 0; i < 4; ++i) cbuf[i]     = uint8_t(uint32_t(csz) >> (8*i));
                        for (int i = 0; i < 4; ++i) cbuf[4 + i] = uint8_t(uint32_t(usz) >> (8*i));
                        const uint8_t* d = cbuf.data(); size_t rem = csz + 8, done = 0;
                        while (done < rem) {
                            ssize_t w = ::write(bfds[b], d + done, rem - done);
                            if (w <= 0) throw std::runtime_error("bucket write failed");
                            done += size_t(w);
                        }
                        lbuf[b].clear();
                    };
                    for (auto& s : snaps) {
                        for (uint32_t r = 0; r < uint32_t(s.rec_hog.size()); ++r) {
                            uint32_t h = s.rec_hog[r];
                            size_t b = std::hash<std::string_view>{}(std::string_view(hog[h])) % B;
                            const auto& hdr = s.rec_hdr[r];
                            uint32_t cnum = parse_cnum(s.contig_strings[hdr.contig_idx], hdr.contig_idx);
                            auto& buf = lbuf[b];
                            write_u32(buf, h);
                            write_u32(buf, s.acc_idx);
                            write_u32(buf, cnum);
                            write_u32(buf, hdr.sstart);
                            write_u32(buf, hdr.send);
                            buf.push_back(uint8_t(std::min(100.0f, hdr.pident + 0.5f)));
                            write_u32(buf, uint32_t(hdr.obs_n));
                            for (uint16_t k = 0; k < hdr.obs_n; ++k) {
                                uint16_t v = s.obs_pool[hdr.obs_off + k].v;
                                write_u32(buf, uint32_t(v >> 6));
                                buf.push_back(uint8_t(v & 0x3F));
                            }
                            if (buf.size() >= CHUNK) drain(b);
                        }
                    }
                    for (size_t b = 0; b < B; ++b) if (!lbuf[b].empty()) drain(b);
                });
        };

        uint64_t& n_written     = w_written[tid];
        uint64_t& n_skipped     = w_skipped[tid];
        uint64_t& n_obs_dropped = w_obs_dropped[tid];

        auto process = [&](uint32_t acc_idx, std::string_view line) {
            std::array<std::string_view, 14> f{};
            size_t fi = 0;
            {
                const char* lp = line.data(); size_t llen = line.size();
                uint32_t tabs[13], n = 0;
#ifdef __AVX2__
                const __m256i vt = _mm256_set1_epi8('\t');
                size_t i = 0;
                for (; i + 32 <= llen && n < 13; i += 32) {
                    uint32_t m = (uint32_t)_mm256_movemask_epi8(
                        _mm256_cmpeq_epi8(_mm256_loadu_si256((const __m256i*)(lp+i)), vt));
                    while (m && n < 13) { tabs[n++] = uint32_t(i + __builtin_ctz(m)); m &= m-1; }
                }
                for (; i < llen && n < 13; ++i) if (lp[i] == '\t') tabs[n++] = uint32_t(i);
#elif defined(__ARM_NEON)
                const uint8x16_t vt = vdupq_n_u8('\t');
                size_t i = 0;
                for (; i + 16 <= llen && n < 13; i += 16) {
                    uint8x16_t chunk = vld1q_u8((const uint8_t*)(lp+i));
                    uint8x16_t eq    = vceqq_u8(chunk, vt);
                    if (vmaxvq_u8(eq))
                        for (size_t j = i; j < i+16 && n < 13; ++j)
                            if ((uint8_t)lp[j] == '\t') tabs[n++] = uint32_t(j);
                }
                for (; i < llen && n < 13; ++i) if ((uint8_t)lp[i] == '\t') tabs[n++] = uint32_t(i);
#else
                const char* fp = lp, *ep = lp + llen;
                while (n < 13) {
                    const char* t = (const char*)memchr(fp, '\t', size_t(ep-fp));
                    if (!t) break;
                    tabs[n++] = uint32_t(t - lp); fp = t + 1;
                }
#endif
                uint32_t prev = 0;
                for (uint32_t t = 0; t < n; ++t) { f[fi++] = {lp+prev, tabs[t]-prev}; prev = tabs[t]+1; }
                if (fi < 14) f[fi++] = {lp+prev, llen-prev};
            }
            if (fi < 13) { ++n_skipped; return; }
            bool has_nt = (fi == 14 && !f[col::full_qseq].empty());

            float pident = 0.0f; double ev = 0.0;
            { auto [ptr,ec] = std::from_chars(f[col::pident].data(), f[col::pident].data()+f[col::pident].size(), pident);
              if (ec != std::errc{}) { ++n_skipped; return; } }
            { auto [ptr,ec] = std::from_chars(f[col::evalue].data(), f[col::evalue].data()+f[col::evalue].size(), ev);
              if (ec != std::errc{}) { ++n_skipped; return; } }
            if (pident < opts.min_pident) { ++n_skipped; return; }
            if (ev     > opts.max_evalue) { ++n_skipped; return; }
            if (pident == 100.0f)         { ++n_skipped; return; }

            uint32_t sstart=0,send=0,qstart=0,qend=0,qlen=0;
            {
                auto fc=[](std::string_view sv,uint32_t& out)->bool{
                    auto [p,ec]=std::from_chars(sv.data(),sv.data()+sv.size(),out);
                    return ec==std::errc{};};
                if(!fc(f[col::sstart],sstart)||!fc(f[col::send],send)||
                   !fc(f[col::qstart],qstart)||!fc(f[col::qend],qend)||
                   !fc(f[col::qlen],qlen)){++n_skipped;return;}
            }
            if (sstart > send) { ++n_skipped; return; }
            if (!has_nt) return;

            size_t la_idx = acc_idx / N;
            if (la_idx >= local_accs.size()) local_accs.resize(la_idx + 1);
            LocalAcc& la = local_accs[la_idx];
            std::string_view qseqid = f[col::qseqid];

            size_t n_hogs_before = hog_strings.size();
            uint32_t hog_idx = intern(hog_dict, hog_strings, extract_hog(f[col::sseqid]));
            if (hog_strings.size() > n_hogs_before) mem_bytes += hog_strings.back().size() + 64;

            size_t n_contigs_before = la.contig_strings.size();
            uint32_t contig_idx = intern(la.contig_dict, la.contig_strings, qseqid);
            if (la.contig_strings.size() > n_contigs_before) mem_bytes += qseqid.size() + 64;

            int8_t qframe = make_qframe(f[col::qstrand], qstart, qend, qlen);
            cigar_parse_inplace(f[col::cigar], f[col::qseq_aa], sstart, send, ar_buf);

            std::string_view full_nt = f[col::full_qseq];
            uint32_t span = send - sstart + 1;
            uint32_t obs0 = uint32_t(la.obs_pool.size());
            for (uint32_t i = 0; i < span; ++i) {
                uint8_t obs_aa = ar_buf.aas[i];
                if (obs_aa == AA_GAP || obs_aa == AA_UNK) { ++n_obs_dropped; continue; }
                uint32_t q_off = ar_buf.qseq_offsets[i];
                if (q_off == UINT32_MAX) { ++n_obs_dropped; continue; }
                uint8_t c0,c1,c2;
                if (qframe > 0) {
                    size_t cs = size_t(qstart-1) + size_t(q_off)*3;
                    if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                    c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
                } else {
                    if (size_t(q_off)*3+3 > size_t(qstart)) { ++n_obs_dropped; continue; }
                    size_t cs = size_t(qstart-1) - size_t(q_off)*3 - 2;
                    if (cs+2 >= full_nt.size()) { ++n_obs_dropped; continue; }
                    c0=uint8_t(full_nt[cs]); c1=uint8_t(full_nt[cs+1]); c2=uint8_t(full_nt[cs+2]);
                    uint8_t tmp[3]={c0,c1,c2}; revcomp_codon(tmp); c0=tmp[0]; c1=tmp[1]; c2=tmp[2];
                }
                if (!(k_acgt_lut[c0] & k_acgt_lut[c1] & k_acgt_lut[c2])) { ++n_obs_dropped; continue; }
                uint8_t raw3[3]={c0,c1,c2};
                uint8_t packed = pack_codon(raw3);
                if (codon_to_aa(packed) != obs_aa) { ++n_obs_dropped; continue; }
                if (i >= 1024) { ++n_obs_dropped; continue; }  // ceiling: 10-bit hog_offset
                la.obs_pool.push_back({uint16_t((i << 6) | (packed >> 2))});
            }
            uint32_t obs_n = uint32_t(la.obs_pool.size()) - obs0;
            if (obs_n > 0) {
                la.rec_hog.push_back(hog_idx);
                la.max_hog_idx = std::max(la.max_hog_idx, hog_idx);
                la.rec_hdr.push_back({contig_idx, sstart, send, obs0, pident, uint16_t(obs_n), qframe, 0, ev});
                ++n_written;
                mem_bytes += obs_n * sizeof(PackedObs) + sizeof(RecHdr) + sizeof(uint32_t);
                if (la.obs_pool.capacity() > la.tracked_obs_cap) {
                    // charge delta × 2: covers old-buffer copy during realloc
                    mem_bytes += (la.obs_pool.capacity() - la.tracked_obs_cap) * sizeof(PackedObs) * 2;
                    la.tracked_obs_cap = la.obs_pool.capacity();
                }
                if (la.rec_hdr.capacity() > la.tracked_hdr_cap) {
                    mem_bytes += (la.rec_hdr.capacity() - la.tracked_hdr_cap) * sizeof(RecHdr);
                    la.tracked_hdr_cap = la.rec_hdr.capacity();
                }
                if (la.rec_hog.capacity() > la.tracked_hog_cap) {
                    mem_bytes += (la.rec_hog.capacity() - la.tracked_hog_cap) * sizeof(uint32_t);
                    la.tracked_hog_cap = la.rec_hog.capacity();
                }
            }

            if (mem_bytes > thr) {
                flush();                                  // normal path: async double-buffer
            } else if ((++rss_tick & 0xFF) == 0 && read_rss_bytes() > rss_cap) {
                // Hard RSS ceiling. The estimate undercounts (contig_strings, dicts, arenas),
                // so this real-RSS check is the true bound. Flush SYNCHRONOUSLY — wait for the
                // write to free the snapshot before accumulating again — otherwise the async
                // snapshot stays resident while new data piles on and RSS blows past budget.
                flush();
                if (flush_future.valid()) flush_future.get();
                // Return the freed slabs to the OS. Without this, glibc keeps them in its
                // per-thread arenas as fragmented free chunks (845M alloc/free churn over 16
                // arenas), so RSS only ever climbs and this guard becomes a no-op — peak RSS
                // settles ~2.7x the working set (65 GB vs a 24 GB budget). Trimming only here
                // (the over-budget path, ~rare) keeps it off the hot flush path.
                malloc_trim(0);
            }
        }; // end process lambda

        Batch batch;
        while (true) {
            auto _qt0 = std::chrono::steady_clock::now();
            if (!queues[tid]->pop(batch)) break;
            n_queue_ns += uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-_qt0).count());
            for (auto& wi : batch.items)
                process(wi.acc_idx, {batch.arena.data() + wi.offset, wi.len});
        }

        if (mem_bytes > 0) flush();
        if (flush_future.valid()) { auto _t0=std::chrono::steady_clock::now(); flush_future.get(); n_stall_ns+=uint64_t(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now()-_t0).count()); }
        // Persist this worker's local hog registry (local idx → name) so merge-shard can
        // remap the local hog ids in this worker's buckets to global ids.
        write_acc_registry(hog_strings, out_dir + "/t" + std::to_string(tid) + ".hog.registry");
    };

    // ── Launch workers ────────────────────────────────────────────────────────
    std::vector<std::thread> pool;
    pool.reserve(N);
    for (size_t t = 0; t < N; ++t) pool.emplace_back(do_worker, t);

    // ── Reader / dispatcher ───────────────────────────────────────────────────
    // Multi-frame .zst → parallel decode + N dispatcher threads (decode no longer the cap, and
    // the per-line acc-extract/route is spread across threads). Else → single reader thread.
    {
        bool auto_acc = (opts.acc_id == "auto");
        if (!auto_acc) { acc_vec.push_back(opts.acc_id); acc_map.emplace(acc_vec.back(), 0u); }
        std::mutex acc_mtx;   // guards the shared accession registry across parallel dispatchers

        // Resolve a line's accession to a global id (per-dispatcher local cache → no lock on hit;
        // shared registry under acc_mtx on miss) and route it to queue[acc%N] via the dispatcher's
        // own pending batches (WorkQueue::push is multi-producer-safe). Shared by both paths.
        auto dispatch_line = [&](std::string_view line, SvDict& lcache, std::vector<Batch>& pend) {
            if (line.empty() || line[0] == '#') return;
            uint32_t acc_idx = 0;
            if (auto_acc) {
                const char* p = line.data();
                const char* tab = static_cast<const char*>(memchr(p, '\t', line.size()));
                const char* qend = tab ? tab : p + line.size();
                const char* us = static_cast<const char*>(memchr(p, '_', size_t(qend - p)));
                std::string_view acc_sv(p, us ? size_t(us - p) : size_t(qend - p));
                auto it = lcache.find(acc_sv);
                if (it != lcache.end()) acc_idx = it->second;
                else {
                    std::lock_guard<std::mutex> lk(acc_mtx);
                    auto git = acc_map.find(acc_sv);
                    if (git != acc_map.end()) acc_idx = git->second;
                    else { acc_idx = uint32_t(acc_vec.size()); acc_vec.emplace_back(acc_sv);
                           acc_map.emplace(acc_vec.back(), acc_idx); }
                    lcache.emplace(std::string(acc_sv), acc_idx);
                }
            }
            size_t tid = size_t(acc_idx) % N;
            pend[tid].push(acc_idx, line.data(), line.size());
            if (pend[tid].full()) queues[tid]->push(std::move(pend[tid]));
        };

        // Detect multi-frame .zst via a capped peek of frame 0.
        bool parallel = false; int zfd = -1;
        if (N > 1 && tsv_path.size() > 4 && tsv_path.compare(tsv_path.size()-4, 4, ".zst") == 0) {
            zfd = ::open(tsv_path.c_str(), O_RDONLY);
            if (zfd >= 0) {
                posix_fadvise(zfd, 0, 0, POSIX_FADV_SEQUENTIAL);
                struct stat st{}; uint64_t fsz = (fstat(zfd, &st) == 0) ? uint64_t(st.st_size) : 0;
                uint64_t f0 = 0;
                try { f0 = ParZstd::frame_len(zfd, 0, fsz, uint64_t(256) << 20); } catch (...) { f0 = 0; }
                parallel = (f0 != 0 && f0 != ParZstd::CAP_EXCEEDED && f0 < fsz);
                if (!parallel) { ::close(zfd); zfd = -1; }
            }
        }

        if (parallel) {
            size_t dthreads = std::max<size_t>(2, N / 2);                 // decode pool
            // Flat, flush-bounded decode footprint: streamed decoded pieces (never whole frames),
            // with the compressed-queue and decoded-inflight caps tied to a fraction of --flush so
            // peak scales with the budget instead of the frame/input size.
            size_t dcap = std::max<size_t>(size_t(128) << 20, flush_threshold / 16);
            ParZstd pz(zfd, dthreads, dcap, dcap);
            std::vector<std::thread> disp;
            bool derr = false; std::string derrmsg; std::mutex demtx;
            for (size_t i = 0; i < N; ++i) disp.emplace_back([&] {
                try {
                    SvDict lcache; std::vector<Batch> pend(N);
                    std::vector<char> chunk, prefix; std::string bl;
                    while (pz.next(chunk, prefix)) {
                        const char* base = chunk.data(); size_t sz = chunk.size();
                        const char* nl = static_cast<const char*>(memchr(base, '\n', sz));
                        if (nl) {
                            bl.assign(prefix.data(), prefix.size());
                            bl.append(base, size_t(nl - base));
                            if (!bl.empty() && bl.back() == '\r') bl.pop_back();
                            dispatch_line(bl, lcache, pend);
                            size_t pos = size_t(nl - base) + 1;
                            while (pos < sz) {
                                const char* nl2 = static_cast<const char*>(memchr(base + pos, '\n', sz - pos));
                                size_t len = nl2 ? size_t(nl2 - (base + pos)) : (sz - pos);
                                std::string_view lv(base + pos, len);
                                if (!lv.empty() && lv.back() == '\r') lv.remove_suffix(1);
                                dispatch_line(lv, lcache, pend);
                                if (!nl2) break;
                                pos += len + 1;
                            }
                        } else {
                            bl.assign(prefix.data(), prefix.size());
                            bl.append(base, sz);
                            if (!bl.empty() && bl.back() == '\r') bl.pop_back();
                            dispatch_line(bl, lcache, pend);
                        }
                    }
                    for (size_t t = 0; t < N; ++t) if (!pend[t].items.empty()) queues[t]->push(std::move(pend[t]));
                } catch (const std::exception& e) {
                    std::lock_guard<std::mutex> lk(demtx); if (!derr) { derr = true; derrmsg = e.what(); }
                }
            });
            for (auto& t : disp) t.join();
            ::close(zfd);
            if (derr) throw std::runtime_error("ingest dispatch: " + derrmsg);
        } else {
            SvDict lcache; std::vector<Batch> pend(N);
            TsvReader reader(tsv_path, 1);   // serial decode (single/huge frame, .gz, plain, stdin)
            std::string_view line;
            while (reader.next_view(line)) dispatch_line(line, lcache, pend);
            for (size_t t = 0; t < N; ++t) if (!pend[t].items.empty()) queues[t]->push(std::move(pend[t]));
        }
        for (auto& q : queues) q->finish();
    }

    for (auto& t : pool) t.join();

    // ── Finalize the fused spill dir: acc registry + spill meta (B, N). The per-worker
    // hog registries and bucket files were already written by each worker on exit. There is
    // no partition index — merge-shard reads buckets directly.
    write_acc_registry(acc_vec, out_dir + "/acc.registry");
    write_spill_meta(out_dir + "/spill.meta", uint32_t(B), uint32_t(N));

    uint64_t tot_written=0, tot_skipped=0, tot_obs=0;
    uint64_t wmax=0, wmin=UINT64_MAX;
    for (size_t t=0;t<N;++t){
        tot_written+=w_written[t]; tot_skipped+=w_skipped[t]; tot_obs+=w_obs_dropped[t];
        wmax=std::max(wmax,w_written[t]); wmin=std::min(wmin,w_written[t]);
    }
    std::cerr << "ingest: " << tot_written << " records, " << tot_skipped << " skipped, "
              << tot_obs << " obs dropped, "
              << acc_vec.size() << " accessions, "
              << B << " buckets → " << out_dir
              << " (" << N << " threads)\n";
    std::cerr << "ingest: worker load [";
    for (size_t t=0;t<N;++t) std::cerr << (t?",":"") << "t" << t << ":" << w_written[t];
    double skew = wmax ? double(wmax)/double(tot_written/N) : 1.0;
    std::cerr << "] skew=" << std::fixed << std::setprecision(2) << skew
              << "× (max/avg, ideal=1.00)\n";
    uint64_t tot_stall=0, tot_queue=0;
    for (size_t t=0;t<N;++t){ tot_stall+=w_stall_ns[t]; tot_queue+=w_queue_ns[t]; }
    std::cerr << "ingest: flush stall " << tot_stall/1000000 << " ms total ("
              << tot_stall/1000000/N << " ms/worker avg)";
    if (tot_queue) std::cerr << "\ningest: input-stall " << tot_queue/1000000
                             << " ms total (" << tot_queue/1000000/N << " ms/worker avg)";
    std::cerr << "\n";
}

// ── Public entry point ────────────────────────────────────────────────────────
inline void ingest(const std::string& tsv_path, const std::string& out_dir,
                   ConvertOptions opts, int n_threads_arg = 0) {
    size_t flush_threshold = compute_flush_threshold(opts.flush_bytes);
    size_t N = (n_threads_arg > 0) ? size_t(n_threads_arg) : 1u;
    ingest_mt(tsv_path, out_dir, opts, N, flush_threshold);
}

} // namespace lhi
