#pragma once
#include <ctime>
#include "global.hpp"
#include "batch.hpp"
#include "container.hpp"
#include <lz4.h>
#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <fcntl.h>
#include <unistd.h>
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
#include <mutex>
#include <cstdint>
#include <cstdio>
#include <zstd.h>

static inline uint64_t clock_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return uint64_t(ts.tv_sec) * 1000000000ull + uint64_t(ts.tv_nsec);
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

// Merge multiple InvBlocks (one per source .lhg) into a single InvBlock.
// acc_remap[i] maps source i's local acc_idx → global acc_idx.
// Unitig dicts are local per-block; output gets a merged local dict whose
// numeric contig values are deduped per (global_acc_idx, contig_num).
inline InvBlock merge_inv_blocks(
        std::vector<std::pair<InvBlock, std::vector<uint32_t>>> sources) {

    // hog_pos → merged obs (acc_idx already global; unitig_idx into merged dict).
    std::unordered_map<uint32_t, std::vector<InvObs>> merged;
    std::vector<uint32_t> out_unitigs;
    std::unordered_map<uint64_t, uint32_t> uni_key_to_idx;  // (acc<<32|cnum) → idx
    auto unitig_idx_of = [&](uint32_t gacc, uint32_t cnum) -> uint32_t {
        uint64_t key = (uint64_t(gacc) << 32) | cnum;
        auto it = uni_key_to_idx.find(key);
        if (it != uni_key_to_idx.end()) return it->second;
        uint32_t idx = uint32_t(out_unitigs.size());
        out_unitigs.push_back(cnum);
        uni_key_to_idx.emplace(key, idx);
        return idx;
    };

    for (auto& [blk, remap] : sources) {
        const uint8_t* p   = blk.pos_ptr;
        const uint8_t* end = blk.end;
        uint32_t n_positions = 0;
        int n = read_varint(p, end, &n_positions);
        if (!n) throw std::runtime_error("merge_inv_blocks: corrupt position count");
        p += n;
        std::vector<InvObs> obs;
        uint32_t prev_pos = 0;
        for (uint32_t pi = 0; pi < n_positions; ++pi) {
            uint32_t hog_pos = 0;
            if (!decode_position(p, end, prev_pos, hog_pos, obs))
                throw std::runtime_error("merge_inv_blocks: corrupt position record");
            auto& dst = merged[hog_pos];
            for (const auto& o : obs) {
                if (o.acc_idx >= remap.size())
                    throw std::runtime_error("merge_inv_blocks: acc_idx out of remap range");
                uint32_t gacc = remap[o.acc_idx];
                if (o.unitig_idx >= blk.unitigs.size())
                    throw std::runtime_error("merge_inv_blocks: unitig_idx OOB");
                uint32_t cnum = blk.unitigs[o.unitig_idx];
                dst.push_back({gacc, o.codon_idx, unitig_idx_of(gacc, cnum)});
            }
        }
    }

    std::vector<InvPosition> positions;
    positions.reserve(merged.size());
    for (auto& kv : merged) {
        std::sort(kv.second.begin(), kv.second.end(),
            [](const InvObs& a, const InvObs& b) { return a.acc_idx < b.acc_idx; });
        positions.push_back({kv.first, std::move(kv.second)});
    }
    std::sort(positions.begin(), positions.end(),
        [](const InvPosition& a, const InvPosition& b) { return a.hog_pos < b.hog_pos; });

    InvBlock out;
    out.unitigs = std::move(out_unitigs);
    serialize_inverted_block(out.raw, out.unitigs, positions);
    out.pos_ptr = nullptr;  // serialized form; callers read out.raw / out.unitigs
    out.end     = out.raw.data() + out.raw.size();
    return out;
}

// accregistry subcommand: collect unique acc_ids from .lhgi LHGA sections,
// sort, write one per line (plain text, no header).
inline void build_acc_registry(const std::vector<std::string>& lhgi_paths,
                               const std::string& out_path) {
    std::set<std::string> accs;
    for (const auto& path : lhgi_paths) {
        GlobalIndex idx;
        if (!idx.load(path) && !idx.load_from_lhg(path))
            throw std::runtime_error("cannot load index: " + path);
        for (const auto& a : idx.accessions) accs.insert(a);
    }
    std::ofstream out(out_path, std::ios::trunc);
    if (!out) throw std::runtime_error("cannot create: " + out_path);
    for (const auto& a : accs) out << a << "\n";
    std::cerr << "accregistry: " << accs.size() << " accessions from "
              << lhgi_paths.size() << " indexes → " << out_path << "\n";
}

// Load a global acc registry (one acc_id per line; position = global acc_idx).
inline std::vector<std::string> load_acc_registry(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open acc registry: " + path);
    std::vector<std::string> accs;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) accs.push_back(line);
    }
    return accs;
}

inline void merge_batches(const std::vector<std::string>& input_paths,
                          const std::string& out_lhg,
                          const std::string& out_lhgi,
                          int zstd_level = 19,
                          const std::string& hog_range_start = "",
                          const std::string& hog_range_end   = "",
                          const std::string& acc_registry_path = "",
                          int n_buckets = 1,
                          int n_threads_override = 0,
                          bool do_profile = false,
                          int hot_threshold = 100,
                          int out_zstd_level = -1) {
    // out_zstd_level: <0 = LZ4 output (fast, default), >=0 = ZSTD level (smaller, for transfer)
    (void)zstd_level;

    // Pass 1: scan all inputs, collect per-HOG MergeRefs (offsets only).
    std::vector<MergeRef> refs;
    refs.reserve(input_paths.size() * 512);

    // Per-source local accession registry (.lhg only). Populated once per
    // source instead of copied into every MergeRef (Fix 2).
    std::vector<std::vector<std::string>> source_accs(input_paths.size());

    for (size_t fi = 0; fi < input_paths.size(); ++fi) {
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
                    refs.push_back(std::move(r));
                });
            if (!ok) throw std::runtime_error("cannot open batch: " + input_paths[fi]);
        } else {
            // .lhg: read its LHGI index; one MergeRef per HOG entry. The source's
            // full local acc registry is stored once in source_accs[fi].
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
                refs.push_back(std::move(r));
            }
        }
        if (input_paths.size() > 100 && (fi + 1) % 1000 == 0)
            std::cerr << "  scanned " << (fi + 1) << "/" << input_paths.size() << " inputs\r";
    }
    if (input_paths.size() > 100) std::cerr << "\n";
    std::cerr << "merge: " << refs.size() << " HOG-blocks from "
              << input_paths.size() << " inputs\n";

    // (hog_id, source_file_idx) order makes HOG groups contiguous for streaming.
    std::stable_sort(refs.begin(), refs.end(),
        [](const MergeRef& a, const MergeRef& b) {
            if (a.hog_id != b.hog_id) return a.hog_id < b.hog_id;
            return a.source_file_idx < b.source_file_idx;
        });

    // Global accession registry.
    std::vector<std::string> accessions;
    if (!acc_registry_path.empty()) {
        accessions = load_acc_registry(acc_registry_path);
    } else {
        std::set<std::string> uniq;
        for (const auto& r : refs)
            if (r.kind == MergeRef::Kind::LHB) uniq.insert(r.lhb_acc_id);
        for (const auto& accs : source_accs)
            for (const auto& a : accs) uniq.insert(a);
        accessions.assign(uniq.begin(), uniq.end());
    }

    std::unordered_map<std::string, uint32_t> acc_id_to_idx;
    acc_id_to_idx.reserve(accessions.size() * 2);
    for (uint32_t i = 0; i < accessions.size(); ++i) acc_id_to_idx[accessions[i]] = i;

    // Per-source acc_idx → global acc_idx remap tables (for .lhg inversion).
    // Validate every input acc exists in the registry (error if registry given).
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
    // pread is positional, so a shared read-only fd is thread-safe on Linux.
    struct SrcFile {
        UniqueFd fd{-1};
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
                fd_opened[fi].store(true, std::memory_order_release);
            }
        }
        return src_files[fi];
    };

    // Fix 3 step 1: pre-scan contiguous HOG groups (refs are stable_sorted by
    // hog_id, then honoring --hog-range filtering up front).
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

    // Bucket assignment: a HOG's bucket is a deterministic hash of its hog_id,
    // NOT its group index. This is the load-bearing invariant — bucket R must
    // hold the identical HOG set in EVERY shard, so a stage-2 job reading
    // "*.lhg.bR from all shards" sees one complete, disjoint HOG partition. An
    // index-relative formula (g*B/n_groups) would split the same HOG into
    // different buckets across shards with differing group counts, double-
    // counting it at stage-2. FNV-1a over the hog_id bytes spreads ~evenly.
    // n_buckets==1 ⇒ every group → bucket 0 (single output, byte-identical).
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

    // Largest-first dispatch: process big HOGs first to minimise tail latency.
    // Precompute per-group total compressed size for largest-first dispatch.
    std::vector<uint64_t> group_sz(n_groups, 0);
    for (size_t g = 0; g < n_groups; ++g)
        for (size_t i = groups[g].first; i < groups[g].second; ++i)
            group_sz[g] += refs[i].lhg_data_length;
    std::vector<size_t> dispatch_order(n_groups);
    std::iota(dispatch_order.begin(), dispatch_order.end(), 0);
    std::sort(dispatch_order.begin(), dispatch_order.end(),
              [&](size_t a, size_t b) { return group_sz[a] > group_sz[b]; });

    // Per-phase profiling counters (summed across all threads).
    struct ProfCounters {
        std::atomic<uint64_t> ns_decode{0}, ns_build{0}, ns_compress{0}, n_groups{0};
    } prof;

    // Hot-HOG parallel decode state.
    // Workers only fill blocks[bi] (pread+ZSTD); src_pos is built by the owner
    // thread in bi order after pending==0, matching the cold path's traversal order
    // so the unitig dict comes out identical and the output is byte-for-byte the same.
    struct HotDecodeTask { size_t g; size_t bi; };
    struct HotHogState {
        std::vector<InvBlock> blocks;
        std::atomic<int>      pending{0};
    };
    std::vector<std::unique_ptr<HotHogState>> hot_states(groups.size());

    std::deque<HotDecodeTask> dq_decode;
    std::mutex                dq_mtx;

    // Fix 3 step 2: per-group result (compressed payload, ready to write).
    struct GroupResult {
        std::string          hog_id;
        std::vector<uint8_t> payload;  // serialized "LHHE" body sans magic/size
        uint32_t             stored_sz = 0;  // size field (high bit = raw fallback)
        uint32_t             n_accs = 0;
    };
    std::vector<GroupResult> results(groups.size());

    // Output files opened here (before process_group + pool) so the async
    // writer thread can start immediately alongside the compute pool.
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

    // Per-result ready flag: worker stores true (release) after results[g] is
    // fully populated; writer loads with acquire before reading.
    std::unique_ptr<std::atomic<bool>[]> result_ready(new std::atomic<bool>[n_groups]);
    for (size_t i = 0; i < n_groups; ++i)
        result_ready[i].store(false, std::memory_order_relaxed);

    // Per-worker reusable containers — cleared (not reallocated) per group.
    struct WorkerScratch {
        std::vector<uint32_t>   unitigs;
        std::unordered_map<uint64_t, uint32_t> unitig_map;
        // Distinct-acc counter: seen_epoch[gacc]==epoch ⇒ already counted this
        // group. Bump epoch per group instead of clearing — O(obs), no per-acc
        // alloc, no full sort (the set<> and the sort+unique both cost more).
        std::vector<uint32_t>   seen_epoch;
        uint32_t                epoch = 0;
        uint32_t                n_accs = 0;
        std::vector<InvPosition> positions;
        // all_lhg path
        std::vector<InvBlock>   blocks;
        std::vector<std::vector<InvPosition>> src_pos;
        std::vector<InvObs>     obs;
        // per-block obs accumulator for the all_lhg and cold heap merge paths:
        // indexed by src/bi; avoids O(n log n) sort since blocks have disjoint acc ranges.
        std::vector<std::vector<InvObs>> per_block_obs;
        // mixed/lhb path
        std::unordered_map<uint32_t, std::vector<InvObs>> inverted;
        std::vector<VarNTRecord> recs;
        std::vector<uint32_t>   contig_cnum;
    };

    // decode_one_block: pread + LZ4 decompress one source block into hot_states[g]->blocks[bi].
    // Safe to call from any worker — each (g,bi) slot written by exactly one thread.
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

    // Per-group inversion + serialize + compress. Self-contained; uses only the
    // caller's read-only state plus worker-owned contexts/buffers.
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

        // per-HOG local unitig dict: store contig numbers, dedup by
        // (global_acc_idx, contig_num) — the unitig identity.
        auto& unitigs    = sc.unitigs;    unitigs.clear();
        auto& unitig_map = sc.unitig_map; unitig_map.clear();
        auto unitig_idx_of = [&](uint32_t gacc, uint32_t cnum) -> uint32_t {
            uint64_t key = (uint64_t(gacc) << 32) | cnum;
            auto it = unitig_map.find(key);
            if (it != unitig_map.end()) return it->second;
            uint32_t idx = uint32_t(unitigs.size());
            unitigs.push_back(cnum);
            unitig_map.emplace(key, idx);
            return idx;
        };
        // Distinct global accessions contributing to this HOG. Epoch-marking
        // counter: count each gacc once, O(1) per obs, no alloc after warmup.
        if (sc.seen_epoch.size() < accessions.size()) sc.seen_epoch.assign(accessions.size(), 0);
        ++sc.epoch;
        sc.n_accs = 0;
        auto mark_acc = [&sc](uint32_t gacc) {
            if (sc.seen_epoch[gacc] != sc.epoch) { sc.seen_epoch[gacc] = sc.epoch; ++sc.n_accs; }
        };

        // Fix 4: heap-merge path when every ref is .lhg. Each source block's
        // positions are already ascending; merge them without a map.
        bool all_lhg = true;
        for (size_t bi = 0; bi < n_blocks && all_lhg; ++bi)
            if (refs[group_start + bi].kind != MergeRef::Kind::LHG) all_lhg = false;

        auto& positions = sc.positions; positions.clear();

        uint64_t t0 = clock_ns();
        uint64_t t_dec_end = 0, t_build_end = 0;

        if (all_lhg) {
            bool is_hot = (int(n_blocks) > hot_threshold);

            if (is_hot) {
                // Hot path: allocate shared state, push decode tasks, help drain queue.
                hot_states[g] = std::make_unique<HotHogState>();
                hot_states[g]->blocks.resize(n_blocks);
                hot_states[g]->pending.store(int(n_blocks), std::memory_order_relaxed);

                {
                    std::lock_guard<std::mutex> lk(dq_mtx);
                    for (size_t bi = 0; bi < n_blocks; ++bi)
                        dq_decode.push_back({g, bi});
                }

                // Owner helps drain until all blocks for this group are done.
                while (hot_states[g]->pending.load(std::memory_order_acquire) > 0) {
                    HotDecodeTask dt{g, 0};
                    bool got = false;
                    {
                        std::lock_guard<std::mutex> lk(dq_mtx);
                        if (!dq_decode.empty()) {
                            dt = dq_decode.front(); dq_decode.pop_front();
                            got = true;
                        }
                    }
                    if (got) decode_one_block(dt.g, dt.bi);
                    else     std::this_thread::yield();
                }

                t_dec_end = clock_ns();

                // Owner thread: parse + remap each pre-decoded block in bi=0..n_blocks-1
                // order — identical to the cold path loop — so unitig_idx_of is called in
                // the same insertion order, producing a byte-identical unitig dict.
                auto& blocks  = hot_states[g]->blocks;
                auto& src_pos = sc.src_pos; src_pos.clear(); src_pos.resize(n_blocks);
                for (size_t bi = 0; bi < n_blocks; ++bi) {
                    const auto& remap = src_remap[refs[group_start + bi].source_file_idx];
                    const InvBlock& blk = blocks[bi];
                    const uint8_t* p   = blk.pos_ptr;
                    const uint8_t* end = blk.end;
                    uint32_t n_positions = 0;
                    int n = read_varint(p, end, &n_positions);
                    if (!n) throw std::runtime_error("corrupt position count for HOG " + hog_id);
                    p += n;
                    src_pos[bi].reserve(n_positions);
                    auto& obs = sc.obs;
                    uint32_t prev_pos = 0;
                    for (uint32_t pi = 0; pi < n_positions; ++pi) {
                        uint32_t hog_pos = 0;
                        if (!decode_position(p, end, prev_pos, hog_pos, obs))
                            throw std::runtime_error("corrupt position record for HOG " + hog_id);
                        InvPosition ip;
                        ip.hog_pos = hog_pos;
                        ip.obs.reserve(obs.size());
                        for (const auto& o : obs) {
                            if (o.acc_idx >= remap.size())
                                throw std::runtime_error("acc_idx out of source registry for HOG " + hog_id);
                            uint32_t gacc = remap[o.acc_idx];
                            if (o.unitig_idx >= blk.unitigs.size())
                                throw std::runtime_error("unitig_idx OOB for HOG " + hog_id);
                            uint32_t cnum = blk.unitigs[o.unitig_idx];
                            mark_acc(gacc);
                            ip.obs.push_back({gacc, o.codon_idx, unitig_idx_of(gacc, cnum)});
                        }
                        src_pos[bi].push_back(std::move(ip));
                    }
                }

                // K-way min-heap merge (identical to cold path).
                // Each source block holds a disjoint acc_idx range (stage-2 invariant),
                // so concatenating per-block in src order produces sorted output without sort().
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
                    for (size_t bi = 0; bi < n_blocks; ++bi) {
                        out_pos.obs.insert(out_pos.obs.end(),
                                           std::make_move_iterator(sc.per_block_obs[bi].begin()),
                                           std::make_move_iterator(sc.per_block_obs[bi].end()));
                        sc.per_block_obs[bi].clear();
                    }
                    positions.push_back(std::move(out_pos));
                }
            } else {
                // Cold path: existing inline decode logic.
                auto& blocks  = sc.blocks;  blocks.clear();  blocks.resize(n_blocks);
                auto& src_pos = sc.src_pos; src_pos.clear(); src_pos.resize(n_blocks);
                for (size_t bi = 0; bi < n_blocks; ++bi) {
                    const MergeRef& ref = refs[group_start + bi];
                    SrcFile& sf = get_src_file(ref.source_file_idx);
                    HogIndexEntry e;
                    e.hog_id = ref.hog_id;
                    e.data_offset = ref.lhg_data_offset;
                    e.data_length = ref.lhg_data_length;
                    if (!read_hog_inverted_fd(sf.fd, input_paths[ref.source_file_idx], e, blocks[bi]))
                        throw std::runtime_error("failed to read .lhg HOG " + hog_id);

                    const auto& remap = src_remap[ref.source_file_idx];

                    const uint8_t* p   = blocks[bi].pos_ptr;
                    const uint8_t* end = blocks[bi].end;
                    uint32_t n_positions = 0;
                    int n = read_varint(p, end, &n_positions);
                    if (!n) throw std::runtime_error("corrupt position count for HOG " + hog_id);
                    p += n;
                    src_pos[bi].reserve(n_positions);
                    auto& obs = sc.obs;
                    uint32_t prev_pos = 0;
                    for (uint32_t pi = 0; pi < n_positions; ++pi) {
                        uint32_t hog_pos = 0;
                        if (!decode_position(p, end, prev_pos, hog_pos, obs))
                            throw std::runtime_error("corrupt position record for HOG " + hog_id);
                        InvPosition ip;
                        ip.hog_pos = hog_pos;
                        ip.obs.reserve(obs.size());
                        for (const auto& o : obs) {
                            if (o.acc_idx >= remap.size())
                                throw std::runtime_error("acc_idx out of source registry for HOG " + hog_id);
                            uint32_t gacc = remap[o.acc_idx];
                            if (o.unitig_idx >= blocks[bi].unitigs.size())
                                throw std::runtime_error("unitig_idx OOB for HOG " + hog_id);
                            uint32_t cnum = blocks[bi].unitigs[o.unitig_idx];
                            mark_acc(gacc);
                            ip.obs.push_back({gacc, o.codon_idx, unitig_idx_of(gacc, cnum)});
                        }
                        src_pos[bi].push_back(std::move(ip));
                    }
                }

                t_dec_end = clock_ns();

                // K-way min-heap over (hog_pos, src_idx, cursor). Concatenate obs at
                // equal hog_pos in src order — blocks have disjoint acc_idx ranges so
                // this is equivalent to sorting by acc_idx without the O(n log n) cost.
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
                    for (size_t bi = 0; bi < n_blocks; ++bi) {
                        out_pos.obs.insert(out_pos.obs.end(),
                                           std::make_move_iterator(sc.per_block_obs[bi].begin()),
                                           std::make_move_iterator(sc.per_block_obs[bi].end()));
                        sc.per_block_obs[bi].clear();
                    }
                    positions.push_back(std::move(out_pos));
                }
            }
        } else {
            // Mixed/LHB groups: accumulate into a map, then sort (existing path).
            auto& inverted = sc.inverted; inverted.clear();
            for (size_t bi = 0; bi < n_blocks; ++bi) {
                const MergeRef& ref = refs[group_start + bi];
                SrcFile& sf = get_src_file(ref.source_file_idx);
                int src_fd = sf.fd;

                if (ref.kind == MergeRef::Kind::LHB) {
                    uint32_t acc_idx = global_acc(ref.lhb_acc_id);

                    uint8_t hdr_bytes[28];
                    ssize_t nr = ::pread(src_fd, hdr_bytes, sizeof(hdr_bytes), off_t(ref.lhb_shard_hdr_offset));
                    if (nr != ssize_t(sizeof(hdr_bytes)))
                        throw std::runtime_error("pread header failed: " + input_paths[ref.source_file_idx]);
                    if (hdr_bytes[4] != SHARD_BLOCK_VERSION)
                        throw std::runtime_error("unsupported shard version in: " + input_paths[ref.source_file_idx]);
                    uint32_t compressed_sz = read_u32_le(hdr_bytes + 8);
                    uint32_t raw_sz        = read_u32_le(hdr_bytes + 12);

                    if (compressed_sz > 256u * 1024 * 1024)
                        throw std::runtime_error("compressed_sz OOB for HOG " + hog_id);
                    if (raw_sz > 512u * 1024 * 1024)
                        throw std::runtime_error("raw_sz OOB for HOG " + hog_id);

                    cbuf_in.resize(compressed_sz);
                    nr = ::pread(src_fd, cbuf_in.data(), compressed_sz,
                                 off_t(ref.lhb_shard_hdr_offset) + sizeof(hdr_bytes));
                    if (nr != ssize_t(compressed_sz))
                        throw std::runtime_error("pread payload failed: " + input_paths[ref.source_file_idx]);

                    raw_block.resize(raw_sz);
                    size_t rz = ZSTD_decompressDCtx(dctx, raw_block.data(), raw_sz,
                                                    cbuf_in.data(), compressed_sz);
                    if (ZSTD_isError(rz))
                        throw std::runtime_error(std::string("zstd decompress: ") + ZSTD_getErrorName(rz));

                    // Parse local contig dict (the unitig_id strings for this block).
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
                        uint32_t u_idx = unitig_idx_of(acc_idx, contig_cnum[r.contig_idx]);
                        for (const auto& o : r.vars) {
                            uint32_t hog_pos   = r.sstart + o.hog_offset;
                            uint8_t  codon_idx = uint8_t(o.packed_codon >> 2);
                            inverted[hog_pos].push_back({acc_idx, codon_idx, u_idx});
                        }
                    }
                } else {
                    // .lhg source: decode the pre-inverted HOG block and re-invert
                    // into the global acc space.
                    HogIndexEntry e;
                    e.hog_id = ref.hog_id;
                    e.data_offset = ref.lhg_data_offset;
                    e.data_length = ref.lhg_data_length;
                    InvBlock blk;
                    if (!read_hog_inverted_fd(src_fd, input_paths[ref.source_file_idx], e, blk))
                        throw std::runtime_error("failed to read .lhg HOG " + hog_id);

                    // Per-source acc remap: source's local acc registry → global.
                    const auto& remap = src_remap[ref.source_file_idx];

                    const uint8_t* p   = blk.pos_ptr;
                    const uint8_t* end = blk.end;
                    uint32_t n_positions = 0;
                    int n = read_varint(p, end, &n_positions);
                    if (!n) throw std::runtime_error("corrupt position count for HOG " + hog_id);
                    p += n;
                    auto& obs = sc.obs;
                    uint32_t prev_pos = 0;
                    for (uint32_t pi = 0; pi < n_positions; ++pi) {
                        uint32_t hog_pos = 0;
                        if (!decode_position(p, end, prev_pos, hog_pos, obs))
                            throw std::runtime_error("corrupt position record for HOG " + hog_id);
                        auto& dst = inverted[hog_pos];
                        for (const auto& o : obs) {
                            if (o.acc_idx >= remap.size())
                                throw std::runtime_error("acc_idx out of source registry for HOG " + hog_id);
                            uint32_t gacc = remap[o.acc_idx];
                            if (o.unitig_idx >= blk.unitigs.size())
                                throw std::runtime_error("unitig_idx OOB for HOG " + hog_id);
                            uint32_t cnum = blk.unitigs[o.unitig_idx];
                            mark_acc(gacc);
                            dst.push_back({gacc, o.codon_idx, unitig_idx_of(gacc, cnum)});
                        }
                    }
                }
            }

            t_dec_end = clock_ns();

            positions.reserve(inverted.size());
            for (auto& kv : inverted) {
                std::sort(kv.second.begin(), kv.second.end(),
                    [](const InvObs& a, const InvObs& b) { return a.acc_idx < b.acc_idx; });
                positions.push_back({kv.first, std::move(kv.second)});
            }
            std::sort(positions.begin(), positions.end(),
                [](const InvPosition& a, const InvPosition& b) { return a.hog_pos < b.hog_pos; });
        }

        // Serialize inverted block (build phase ends here).
        inv_raw.clear();
        serialize_inverted_block(inv_raw, unitigs, positions);
        t_build_end = clock_ns();

        tl_decode += t_dec_end - t0;
        tl_build  += t_build_end - t_dec_end;

        // stored_sz: bit31=raw · bit30=ZSTD · else=LZ4
        // LZ4 payload: [4B raw_sz LE][lz4 bytes]
        size_t raw_sz = inv_raw.size();
        size_t payload_sz = 0;
        uint32_t stored_sz = 0;
        bool use_raw = false;

        if (out_zstd_level >= 0) {
            size_t bound = ZSTD_compressBound(raw_sz);
            hog_cbuf.resize(bound);
            size_t csz = ZSTD_compress2(cctx, hog_cbuf.data(), bound,
                                        inv_raw.data(), raw_sz);
            if (ZSTD_isError(csz))
                throw std::runtime_error(std::string("zstd HOG compress: ") + ZSTD_getErrorName(csz));
            use_raw = (csz >= raw_sz);
            if (!use_raw) { payload_sz = csz; stored_sz = uint32_t(csz) | 0x40000000u; }
        } else {
            size_t lz4_bound = size_t(LZ4_compressBound(int(raw_sz)));
            hog_cbuf.resize(4 + lz4_bound);
            hog_cbuf[0] = uint8_t(raw_sz);       hog_cbuf[1] = uint8_t(raw_sz >> 8);
            hog_cbuf[2] = uint8_t(raw_sz >> 16); hog_cbuf[3] = uint8_t(raw_sz >> 24);
            int lcsz = LZ4_compress_default(reinterpret_cast<const char*>(inv_raw.data()),
                                            reinterpret_cast<char*>(hog_cbuf.data() + 4),
                                            int(raw_sz), int(lz4_bound));
            size_t total = (lcsz > 0) ? size_t(4 + lcsz) : raw_sz;
            use_raw = (lcsz <= 0 || total >= raw_sz);
            if (!use_raw) { payload_sz = total; stored_sz = uint32_t(total); }
        }
        if (use_raw) { payload_sz = raw_sz; stored_sz = uint32_t(raw_sz) | 0x80000000u; }

        tl_compress += clock_ns() - t_build_end;

        GroupResult& gr = results[g];
        gr.hog_id = hog_id;
        gr.n_accs = sc.n_accs;
        gr.stored_sz = stored_sz;
        if (use_raw) gr.payload.assign(inv_raw.begin(), inv_raw.end());
        else         gr.payload.assign(hog_cbuf.begin(), hog_cbuf.begin() + payload_sz);
        // Signal the async writer that this result slot is ready.
        result_ready[g].store(true, std::memory_order_release);
    };

    // Fix 3 step 3: work-stealing across N workers. Each owns its zstd contexts.
    std::atomic<size_t> next_group{0};
    std::atomic<bool> failed{false};
    std::string first_error;
    std::mutex err_mtx;
    // Default: cap at 32 — bandwidth saturation and NUMA effects cause
    // degradation beyond this on typical 2-socket servers (measured: 32T=6.8s,
    // 96T=8.9s on same workload). Override with -t N.
    size_t hw = std::thread::hardware_concurrency();
    size_t default_threads = std::min<size_t>(32u, hw);
    size_t n_threads = (n_threads_override > 0)
                       ? std::min<size_t>(size_t(n_threads_override), groups.size())
                       : std::min<size_t>(default_threads, groups.size());
    if (n_threads == 0) n_threads = 1;

    auto worker = [&]() {
        ZSTD_CCtx* cctx = ZSTD_createCCtx();
        ZSTD_DCtx* dctx = ZSTD_createDCtx();
        if (out_zstd_level >= 0)
            ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, out_zstd_level);
        std::vector<uint8_t> cbuf_in, raw_block, inv_raw, hog_cbuf;
        WorkerScratch scratch;
        uint64_t tl_decode = 0, tl_build = 0, tl_compress = 0;
        for (;;) {
            size_t di = next_group.fetch_add(1, std::memory_order_relaxed);
            if (di < dispatch_order.size() && !failed.load(std::memory_order_relaxed)) {
                size_t g = dispatch_order[di];
                try {
                    process_group(g, cctx, dctx, cbuf_in, raw_block, inv_raw, hog_cbuf,
                                  scratch, tl_decode, tl_build, tl_compress);
                } catch (const std::exception& e) {
                    std::lock_guard<std::mutex> lk(err_mtx);
                    if (!failed.exchange(true)) first_error = e.what();
                    break;
                }
                continue;
            }
            // HOG queue exhausted; try to drain hot-decode tasks.
            HotDecodeTask dt{0, 0};
            bool got = false;
            {
                std::lock_guard<std::mutex> lk(dq_mtx);
                if (!dq_decode.empty()) { dt = dq_decode.front(); dq_decode.pop_front(); got = true; }
            }
            if (got) {
                try { decode_one_block(dt.g, dt.bi); }
                catch (const std::exception& e) {
                    std::lock_guard<std::mutex> lk(err_mtx);
                    if (!failed.exchange(true)) first_error = e.what();
                    break;
                }
                continue;
            }
            // Nothing to do: check if we can exit.
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

    // Async writer: streams results to NFS as workers produce them,
    // overlapping I/O with computation. Single-threaded; no locking needed
    // on buckets/write_pos/index_entries (only this thread touches them).
    std::string writer_error;
    auto writer_fn = [&]() {
        try {
            for (size_t g = 0; g < n_groups; ++g) {
                while (!result_ready[g].load(std::memory_order_acquire)) {
                    if (failed.load(std::memory_order_relaxed)) return;
                    std::this_thread::yield();
                }
                GroupResult& gr = results[g];
                BucketOut&   bk = buckets[group_bucket[g]];
                uint64_t off = bk.write_pos;
                uint8_t hdr8[8];
                memcpy(hdr8, LHG_HOG_ENTRY_MAGIC, 4);
                for (int i = 0; i < 4; ++i) hdr8[4+i] = uint8_t(gr.stored_sz >> (8*i));
                bk.wbuf.append(bk.lhg_fd, bk.lhg_path, hdr8, 8);
                bk.write_pos += 8;
                bk.wbuf.append(bk.lhg_fd, bk.lhg_path, gr.payload.data(), gr.payload.size());
                bk.write_pos += gr.payload.size();
                bk.index_entries.push_back(
                    {std::move(gr.hog_id), off, bk.write_pos - off, gr.n_accs});
                { std::vector<uint8_t>().swap(gr.payload); }  // free RSS immediately
            }
            for (auto& bk : buckets) bk.wbuf.flush_to(bk.lhg_fd, bk.lhg_path);
        } catch (const std::exception& e) {
            writer_error = e.what();
            failed.store(true, std::memory_order_relaxed);
        }
    };

    // Start async writer then compute pool; both run concurrently.
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
