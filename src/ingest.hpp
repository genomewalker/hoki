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

// ── Single-threaded ingest (original code path, used when -t 1 or not set) ──

inline void ingest_st(const std::string& tsv_path, const std::string& out_dir,
                      ConvertOptions opts, size_t flush_threshold) {
    std::filesystem::create_directories(out_dir);

    struct AccState {
        uint32_t acc_idx;
        SvDict contig_dict;
        std::vector<std::string> contig_strings;
        std::unordered_map<uint32_t, std::vector<VarNTRecord>> batches;
    };

    bool auto_acc = (opts.acc_id == "auto");
    std::unordered_map<std::string, AccState> acc_states;
    SvDict hog_dict;
    std::vector<std::string> hog_strings;

    auto intern = [](SvDict& d, std::vector<std::string>& v, std::string_view s) -> uint32_t {
        auto it = d.find(s);
        if (it != d.end()) return it->second;
        uint32_t idx = uint32_t(v.size());
        v.emplace_back(s);
        d.emplace(v.back(), idx);
        return idx;
    };

    PartitionWriter pw(out_dir + "/t0.lhp", opts.zstd_level);
    size_t mem_bytes = 0;
    uint64_t n_flushes = 0;

    auto flush_batches = [&]() {
        for (auto& [acc, st] : acc_states) {
            for (auto& [hog_idx, batch] : st.batches) {
                if (batch.empty()) continue;
                auto raw = serialize_shard_raw(st.contig_strings, batch);
                pw.append_raw(hog_strings[hog_idx], st.acc_idx, uint32_t(batch.size()), raw.data(), uint32_t(raw.size()));
            }
            st.batches.clear();
            st.contig_dict.clear();
            st.contig_strings.clear();
        }
        // HOG dict stays — 166K names × ~15B = ~2.5MB, free to keep across flushes.
        // Only VarNTRecord batches and per-acc contig dicts drive memory growth.
        mem_bytes = 0;
        ++n_flushes;
        malloc_trim(0);
    };

    TsvReader reader(tsv_path);
    std::string line;
    uint64_t lineno = 0, n_written = 0, n_skipped = 0, n_obs_dropped = 0;
    AlignedResult ar_buf;

    while (reader.getline(line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;

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
        if (fi < 13) { if (opts.verbose) std::cerr << "skip L" << lineno << ": " << fi << " fields\n"; ++n_skipped; continue; }
        bool has_nt = (fi == 14 && !f[col::full_qseq].empty());

        float pident = 0.0f; double ev = 0.0;
        { auto [ptr,ec] = std::from_chars(f[col::pident].data(), f[col::pident].data()+f[col::pident].size(), pident);
          if (ec != std::errc{}) { ++n_skipped; continue; } }
        { auto [ptr,ec] = std::from_chars(f[col::evalue].data(),  f[col::evalue].data()+f[col::evalue].size(),  ev);
          if (ec != std::errc{}) { ++n_skipped; continue; } }
        if (pident < opts.min_pident) { ++n_skipped; continue; }
        if (ev     > opts.max_evalue) { ++n_skipped; continue; }
        if (pident == 100.0f)         { ++n_skipped; continue; }

        uint32_t sstart=0, send=0, qstart=0, qend=0, qlen=0;
        {
            auto fc = [](std::string_view sv, uint32_t& out) -> bool {
                auto [p,ec] = std::from_chars(sv.data(), sv.data()+sv.size(), out);
                return ec == std::errc{};
            };
            if (!fc(f[col::sstart],sstart)||!fc(f[col::send],send)||
                !fc(f[col::qstart],qstart)||!fc(f[col::qend],qend)||
                !fc(f[col::qlen],qlen)) { if (opts.verbose) std::cerr << "parse error L" << lineno << "\n"; ++n_skipped; continue; }
        }
        if (sstart > send) { ++n_skipped; continue; }
        if (!has_nt) continue;

        std::string_view qseqid = f[col::qseqid];
        std::string acc;
        if (auto_acc) {
            auto us = qseqid.find('_');
            acc = std::string(us != std::string_view::npos ? qseqid.substr(0,us) : qseqid);
        } else {
            acc = opts.acc_id;
        }

        bool inserted = (acc_states.count(acc) == 0);
        AccState& st = acc_states[acc];
        if (inserted) st.acc_idx = uint32_t(acc_states.size() - 1);

        size_t n_hogs_before = hog_strings.size();
        uint32_t hog_idx    = intern(hog_dict, hog_strings, extract_hog(f[col::sseqid]));
        if (hog_strings.size() > n_hogs_before) mem_bytes += hog_strings.back().size() + 64;

        size_t n_contigs_before = st.contig_strings.size();
        uint32_t contig_idx = intern(st.contig_dict, st.contig_strings, qseqid);
        if (st.contig_strings.size() > n_contigs_before) mem_bytes += qseqid.size() + 64;

        int8_t qframe = make_qframe(f[col::qstrand], qstart, qend, qlen);
        cigar_parse_inplace(f[col::cigar], f[col::qseq_aa], sstart, send, ar_buf);

        VarNTRecord vr;
        vr.contig_idx = contig_idx; vr.sstart = sstart; vr.send = send;
        vr.qframe = qframe; vr.pident = pident; vr.evalue = ev;

        std::string_view full_nt = f[col::full_qseq];
        uint32_t span = send - sstart + 1;
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
            if (codon_to_aa(packed) != obs_aa) {
                if (opts.verbose) std::cerr << "codon/AA mismatch L" << lineno << "\n";
                ++n_obs_dropped; continue;
            }
            vr.vars.push_back({i, packed});
        }
        if (!vr.vars.empty()) {
            auto& bucket = st.batches[hog_idx];
            if (bucket.empty()) mem_bytes += 64;
            mem_bytes += sizeof(VarNTRecord) + vr.vars.capacity() * sizeof(vr.vars[0]) + 16;
            bucket.push_back(std::move(vr));
            ++n_written;
        }

        bool over = mem_bytes > flush_threshold;
        if (!over && (lineno & 0x3FFFu) == 0) {
            size_t rss = read_rss_bytes();
            if (rss > flush_threshold) {
                std::cerr << "ingest: RSS " << (rss>>20) << " MiB > threshold "
                          << (flush_threshold>>20) << " MiB (estimated "
                          << (mem_bytes>>20) << " MiB) — flushing\n";
                over = true;
            }
        }
        if (over) flush_batches();
    }

    if (mem_bytes > 0) flush_batches();
    pw.flush_all();

    std::vector<std::string> accessions(acc_states.size());
    for (auto& [acc, st] : acc_states) accessions[st.acc_idx] = acc;
    write_acc_registry(accessions, out_dir + "/acc.registry");

    std::map<std::string, std::vector<PartitionIndexExtent>> global_idx;
    for (const auto& [hog, extents] : pw.index())
        for (auto e : extents) { e.thread_idx = 0u; global_idx[hog].push_back(e); }
    write_partition_index(global_idx, 1u, out_dir + "/partition.idx");

    std::cerr << "ingest: " << n_written << " records, " << n_skipped << " skipped, "
              << n_obs_dropped << " obs dropped, "
              << accessions.size() << " accessions, "
              << global_idx.size() << " HOGs → " << out_dir << "\n";
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

    // ── Pre-create PartitionWriters so index() is accessible after join ──────
    std::vector<std::unique_ptr<PartitionWriter>> writers(N);
    for (size_t t = 0; t < N; ++t)
        writers[t] = std::make_unique<PartitionWriter>(
            out_dir + "/t" + std::to_string(t) + ".lhp", opts.zstd_level);

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

        PartitionWriter& pw = *writers[tid];
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
        const size_t rss_cap = (opts.flush_bytes != 0 ? opts.flush_bytes : flush_threshold);
        const size_t thr     = std::max<size_t>(1, rss_cap / (N * 4));

        // cctx is now owned by PartitionWriter; no per-worker CCtx needed here.

        auto intern = [](SvDict& d, std::vector<std::string>& v, std::string_view s) -> uint32_t {
            auto it = d.find(s);
            if (it != d.end()) return it->second;
            uint32_t idx = uint32_t(v.size());
            v.emplace_back(s);
            d.emplace(v.back(), idx);
            return idx;
        };

        std::vector<uint32_t>    sort_counts, sort_order, scatter_pos;  // reused across flushes
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

            // Move sort buffers to async thread; worker gets fresh (empty) vectors.
            auto sc  = std::move(sort_counts);
            auto so  = std::move(sort_order);
            auto sp  = std::move(scatter_pos);

            flush_future = std::async(std::launch::async,
                [snaps = std::move(snaps), hog = std::move(hog_snap),
                 sc = std::move(sc), so = std::move(so), sp = std::move(sp),
                 &pw]() mutable {
                    // Reusable buffers — grown to HWM across the entire async batch
                    std::vector<uint32_t> present, sort_scratch;
                    std::vector<uint64_t> sk;                       // Win 2: packed sort keys
                    std::vector<uint32_t> g2l_slot, g2l_gen, local_gids;  // Win 1: flat g2l
                    uint32_t              g2l_cur = 0;
                    std::vector<uint8_t>  raw_buf, cc_col, ss_col, span_col,
                                          qf_col,  pi_col, ev_col, bmp_col, cdn_col;
                    for (auto& s : snaps) {
                        // Grow flat g2l arrays to cover all contig indices this snap
                        if (g2l_slot.size() < s.contig_strings.size()) {
                            g2l_slot.resize(s.contig_strings.size(), 0);
                            g2l_gen.resize(s.contig_strings.size(), 0);
                        }
                        uint32_t cs_bound = s.max_hog_idx + 2;
                        if (sc.size() < cs_bound) sc.resize(cs_bound);
                        if (sp.size() < cs_bound) sp.resize(cs_bound);
                        std::fill(sc.begin(), sc.begin() + cs_bound, 0);
                        for (uint32_t h : s.rec_hog) ++sc[h + 1];
                        present.clear();
                        for (uint32_t i = 1; i < cs_bound; i++) {
                            if (sc[i] > 0) present.push_back(i - 1);
                            sc[i] += sc[i - 1];
                        }
                        std::copy_n(sc.begin(), cs_bound, sp.begin());
                        so.resize(s.rec_hog.size());
                        for (uint32_t r = 0; r < uint32_t(s.rec_hog.size()); r++)
                            so[sp[s.rec_hog[r]]++] = r;

                        for (uint32_t h : present) {
                            uint32_t lo = sc[h], hi = sc[h + 1];
                            uint32_t N  = hi - lo;

                            // Generation-stamped flat g2l — assign local contig slots in
                            // first-seen order. Must precede the sort below (the sort key and
                            // the contig delta-encoding both read g2l_slot).
                            if (++g2l_cur == 0) { std::fill(g2l_gen.begin(), g2l_gen.end(), 0); g2l_cur = 1; }
                            local_gids.clear();
                            for (uint32_t i = lo; i < hi; i++) {
                                uint32_t ci = s.rec_hdr[so[i]].contig_idx;
                                if (g2l_gen[ci] != g2l_cur) { g2l_gen[ci]=g2l_cur; g2l_slot[ci]=uint32_t(local_gids.size()); local_gids.push_back(ci); }
                            }

                            // Sort records by (contig slot, sstart) — the column format's
                            // contract: contig grouped, sstart monotone within each contig run,
                            // so contig deltas are non-negative and sstart resets per contig.
                            sort_scratch.resize(N);
                            for (uint32_t i = 0; i < N; i++) sort_scratch[i] = so[lo + i];
                            std::stable_sort(sort_scratch.begin(), sort_scratch.end(),
                                [&](uint32_t a, uint32_t b) {
                                    uint32_t sa = g2l_slot[s.rec_hdr[a].contig_idx];
                                    uint32_t sb = g2l_slot[s.rec_hdr[b].contig_idx];
                                    if (sa != sb) return sa < sb;
                                    return s.rec_hdr[a].sstart < s.rec_hdr[b].sstart;
                                });

                            // Columnar serialisation directly from SoA. contig is delta-encoded
                            // against the previous slot; sstart resets to 0 at each contig boundary
                            // (mirrors deserialize_varnt_block's accumulation).
                            cc_col.clear(); ss_col.clear(); span_col.clear(); qf_col.clear();
                            pi_col.clear(); ev_col.clear(); bmp_col.clear(); cdn_col.clear();
                            uint32_t prev_ss = 0, prev_slot = 0;
                            for (uint32_t idx : sort_scratch) {
                                const auto& hdr = s.rec_hdr[idx];
                                uint32_t slot = g2l_slot[hdr.contig_idx];
                                uint32_t dc = slot - prev_slot;
                                write_varint(cc_col,   dc);
                                if (dc) { prev_ss = 0; prev_slot = slot; }
                                write_varint(ss_col,   hdr.sstart - prev_ss); prev_ss = hdr.sstart;
                                write_varint(span_col, hdr.send - hdr.sstart + 1);
                                qf_col.push_back(uint8_t(hdr.qframe));
                                write_u16(pi_col, encode_pident(hdr.pident));
                                write_i16(ev_col, encode_evalue(hdr.evalue));
                                uint32_t bmp_sz = (hdr.send - hdr.sstart + 8) / 8;
                                size_t   bmp_off = bmp_col.size();
                                bmp_col.resize(bmp_off + bmp_sz, 0);
                                for (uint16_t k = 0; k < hdr.obs_n; k++) {
                                    uint16_t v = s.obs_pool[hdr.obs_off + k].v;
                                    uint16_t hog_off = v >> 6;
                                    bmp_col[bmp_off + hog_off / 8] |= uint8_t(1u << (hog_off & 7));
                                    cdn_col.push_back(uint8_t(v & 0x3F));
                                }
                            }

                            size_t raw_est = 5;
                            for (uint32_t gid : local_gids) raw_est += 5 + s.contig_strings[gid].size();
                            raw_est += size_t(N) * 48;
                            raw_buf.clear(); raw_buf.reserve(raw_est);
                            write_varint(raw_buf, uint32_t(local_gids.size()));
                            for (uint32_t gid : local_gids) { const auto& cs=s.contig_strings[gid]; write_varint(raw_buf, uint32_t(cs.size())); raw_buf.insert(raw_buf.end(), cs.begin(), cs.end()); }
                            write_varint(raw_buf, N);
                            write_varint(raw_buf, uint32_t(cc_col.size()));
                            write_varint(raw_buf, uint32_t(ss_col.size()));
                            write_varint(raw_buf, uint32_t(span_col.size()));
                            write_varint(raw_buf, uint32_t(bmp_col.size()));
                            write_varint(raw_buf, uint32_t(cdn_col.size()));
                            raw_buf.insert(raw_buf.end(), cc_col.begin(),   cc_col.end());
                            raw_buf.insert(raw_buf.end(), ss_col.begin(),   ss_col.end());
                            raw_buf.insert(raw_buf.end(), span_col.begin(), span_col.end());
                            raw_buf.insert(raw_buf.end(), qf_col.begin(),   qf_col.end());
                            raw_buf.insert(raw_buf.end(), pi_col.begin(),   pi_col.end());
                            raw_buf.insert(raw_buf.end(), ev_col.begin(),   ev_col.end());
                            raw_buf.insert(raw_buf.end(), bmp_col.begin(),  bmp_col.end());
                            raw_buf.insert(raw_buf.end(), cdn_col.begin(),  cdn_col.end());

                            pw.append_raw(hog[h], s.acc_idx, N, raw_buf.data(), uint32_t(raw_buf.size()));
                        }
                    }
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
        pw.flush_all();
    };

    // ── Launch workers ────────────────────────────────────────────────────────
    std::vector<std::thread> pool;
    pool.reserve(N);
    for (size_t t = 0; t < N; ++t) pool.emplace_back(do_worker, t);

    // ── Reader / dispatcher (main thread) ─────────────────────────────────────
    {
        bool auto_acc = (opts.acc_id == "auto");

        if (!auto_acc) {
            acc_vec.push_back(opts.acc_id);
            acc_map.emplace(acc_vec.back(), 0u);
        }

        std::vector<Batch> pending(N);

        auto flush_pending = [&](size_t tid) {
            if (!pending[tid].items.empty()) {
                queues[tid]->push(std::move(pending[tid]));
                // moved-from Batch has empty arena+items; no explicit reset needed
            }
        };

        TsvReader reader(tsv_path);
        std::string line;
        while (reader.getline(line)) {
            if (line.empty() || line[0] == '#') continue;

            uint32_t acc_idx;
            if (auto_acc) {
                const char* p    = line.data();
                const char* tab  = static_cast<const char*>(memchr(p, '\t', line.size()));
                const char* qend = tab ? tab : p + line.size();
                const char* us   = static_cast<const char*>(memchr(p, '_', size_t(qend - p)));
                std::string_view acc_sv(p, us ? size_t(us - p) : size_t(qend - p));

                auto it = acc_map.find(acc_sv);
                if (it != acc_map.end()) {
                    acc_idx = it->second;
                } else {
                    acc_idx = uint32_t(acc_vec.size());
                    acc_vec.emplace_back(acc_sv);
                    acc_map.emplace(acc_vec.back(), acc_idx);
                }
            } else {
                acc_idx = 0;
            }

            size_t tid = size_t(acc_idx) % N;
            pending[tid].push(acc_idx, line.data(), line.size());
            if (pending[tid].full()) flush_pending(tid);
        }

        for (size_t t = 0; t < N; ++t) flush_pending(t);
        for (auto& q : queues) q->finish();
    }

    for (auto& t : pool) t.join();

    // ── Merge per-worker sidecars → partition.idx ─────────────────────────────
    std::map<std::string, std::vector<PartitionIndexExtent>> global_idx;
    for (size_t t = 0; t < N; ++t)
        for (const auto& [hog, extents] : writers[t]->index())
            for (auto e : extents) { e.thread_idx = uint32_t(t); global_idx[hog].push_back(e); }
    write_partition_index(global_idx, uint32_t(N), out_dir + "/partition.idx");
    write_acc_registry(acc_vec, out_dir + "/acc.registry");

    uint64_t tot_written=0, tot_skipped=0, tot_obs=0;
    uint64_t wmax=0, wmin=UINT64_MAX;
    for (size_t t=0;t<N;++t){
        tot_written+=w_written[t]; tot_skipped+=w_skipped[t]; tot_obs+=w_obs_dropped[t];
        wmax=std::max(wmax,w_written[t]); wmin=std::min(wmin,w_written[t]);
    }
    std::cerr << "ingest: " << tot_written << " records, " << tot_skipped << " skipped, "
              << tot_obs << " obs dropped, "
              << acc_vec.size() << " accessions, "
              << global_idx.size() << " HOGs → " << out_dir
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
    size_t N = (n_threads_arg > 0)
               ? size_t(n_threads_arg)
               : 1u; // default single-threaded; caller passes hardware_concurrency() if desired
    if (N <= 1)
        ingest_st(tsv_path, out_dir, opts, flush_threshold);
    else
        ingest_mt(tsv_path, out_dir, opts, N, flush_threshold);
}

} // namespace lhi
