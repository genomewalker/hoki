#include "convert.hpp"
#include "merge.hpp"
#include "partition.hpp"
#include "global.hpp"
#include "aa.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <optional>
#include <thread>
#include <atomic>
#include <mutex>
#include <algorithm>
#include <sys/stat.h>

static void usage(const char* prog) {
    std::cerr
        << "hoki — HOG codon Index\n"
        << "  " << prog << " convert -a ACC [-z LVL] [-p MINPID] [-e MAXEV] [-v] <in.tsv> <out.lhb>\n"
        << "  " << prog << " merge   [-zo LVL=3] [--buckets N=1] [--hog-range START END]\n"
        << "                  [--profile] [--hot-threshold N=100]\n"
        << "                  <out.lhg> <out.lhgi> <input1> [input2 ...]\n"
        << "                  inputs may be .lhb or .lhg, mixed\n"
        << "  " << prog << " partition   [-t N] <out_dir> <input1.lhb> [input2.lhb ...]\n"
        << "  " << prog << " merge-shard [-t N] [--hog-range START END] <partition_dir> <out.lhg> <out.lhgi>\n"
        << "                  partition_dir must contain partition.idx, acc.registry, t0.lhp ...\n"
        << "  " << prog << " saav    <global.lhg> <global.lhgi> <HOG_ID> <POS> [AA] [--min-pident N]\n"
        << "  " << prog << " freq    <global.lhg> <global.lhgi> <HOG_ID> [--min-pident N]\n"
        << "  " << prog << " stat    <file.lhb | global.lhg [global.lhgi]>\n";
}

int main(int argc, char* argv[]) {
  try {
    if (argc < 2) { usage(argv[0]); return 1; }
    std::string mode = argv[1];

    if (mode == "convert") {
        lhi::ConvertOptions opts;
        std::string in_path, lhb_path;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-a" && i+1 < argc) opts.acc_id     = argv[++i];
            else if (a == "-z" && i+1 < argc) opts.zstd_level = std::stoi(argv[++i]);
            else if (a == "-p" && i+1 < argc) opts.min_pident = std::stof(argv[++i]);
            else if (a == "-e" && i+1 < argc) opts.max_evalue = std::stod(argv[++i]);
            else if (a == "-v")               opts.verbose    = true;
            else if (in_path.empty())         in_path         = a;
            else if (lhb_path.empty())        lhb_path        = a;
        }
        if (lhb_path.empty() || opts.acc_id.empty()) { usage(argv[0]); return 1; }
        lhi::convert(in_path, lhb_path, opts);
        return 0;
    }

    if (mode == "merge") {
        std::string hog_start, hog_end;
        int n_buckets = 1;
        int n_threads = 0;
        bool do_profile = false;
        int hot_threshold = 100;
        int out_compress_level = 3;
        std::vector<std::string> pos;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "--hog-range" && i+2 < argc) { hog_start = argv[++i]; hog_end = argv[++i]; }
            else if (a == "-zo" && i+1 < argc)              out_compress_level= std::stoi(argv[++i]);
            else if (a == "--buckets" && i+1 < argc)        n_buckets         = std::stoi(argv[++i]);
            else if (a == "-t"        && i+1 < argc)        n_threads         = std::stoi(argv[++i]);
            else if (a == "--profile")                      do_profile        = true;
            else if (a == "--hot-threshold" && i+1 < argc) hot_threshold      = std::stoi(argv[++i]);
            else pos.emplace_back(a);
        }
        if (pos.size() < 3) { usage(argv[0]); return 1; }
        std::string out_lhg  = pos[0];
        std::string out_lhgi = pos[1];
        std::vector<std::string> inputs(pos.begin() + 2, pos.end());
        lhi::merge_batches(inputs, out_lhg, out_lhgi, hog_start, hog_end,
                           n_buckets, n_threads, do_profile, hot_threshold, out_compress_level);
        return 0;
    }

    if (mode == "partition") {
        int n_threads = 0;
        std::vector<std::string> pos;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-t" && i+1 < argc) n_threads = std::stoi(argv[++i]);
            else pos.emplace_back(a);
        }
        if (pos.size() < 2) { usage(argv[0]); return 1; }
        std::string out_dir = pos[0];
        std::vector<std::string> inputs(pos.begin() + 1, pos.end());
        lhi::partition_lhbs(inputs, out_dir, n_threads);
        return 0;
    }

    if (mode == "merge-shard") {
        int n_threads = 0;
        std::string hog_range_start, hog_range_end;
        std::vector<std::string> pos;
        for (int i = 2; i < argc; ++i) {
            std::string a = argv[i];
            if      (a == "-t" && i+1 < argc)          n_threads       = std::stoi(argv[++i]);
            else if (a == "--hog-range" && i+2 < argc) { hog_range_start = argv[++i]; hog_range_end = argv[++i]; }
            else pos.emplace_back(a);
        }
        if (pos.size() < 3) { usage(argv[0]); return 1; }
        std::string part_dir = pos[0];
        std::string out_lhg  = pos[1];
        std::string out_lhgi = pos[2];
        std::string idx_path = part_dir + "/partition.idx";
        std::string acc_path = part_dir + "/acc.registry";

        uint32_t n_tfiles = 0;
        auto part_idx = lhi::load_partition_index(idx_path, n_tfiles);
        auto accessions = lhi::load_partition_acc_registry(acc_path);

        // Open all thread files read-only; shared by all workers (pread is thread-safe).
        std::vector<lhi::UniqueFd> tfds;
        std::vector<int>           tfd_ints;
        tfds.reserve(n_tfiles); tfd_ints.reserve(n_tfiles);
        for (uint32_t t = 0; t < n_tfiles; ++t) {
            std::string tp = part_dir + "/t" + std::to_string(t) + ".lhp";
            int raw = ::open(tp.c_str(), O_RDONLY);
            if (raw < 0) { std::cerr << "cannot open: " << tp << "\n"; return 1; }
            tfds.emplace_back(raw);
            tfd_ints.push_back(raw);
        }

        // Filter HOGs by range
        using HogEntry = std::pair<const std::string, std::vector<lhi::PartitionIndexExtent>>;
        std::vector<const HogEntry*> hog_list;
        hog_list.reserve(part_idx.size());
        for (const auto& kv : part_idx) {
            if (!hog_range_start.empty() && kv.first < hog_range_start) continue;
            if (!hog_range_end.empty()   && kv.first > hog_range_end)   continue;
            hog_list.push_back(&kv);
        }

        size_t nt = size_t(n_threads > 0 ? n_threads : int(std::thread::hardware_concurrency()));
        nt = std::max<size_t>(1, std::min(nt, hog_list.size()));

        // Parallel compute: workers steal HOGs, store results in pre-allocated slots.
        std::vector<lhi::ShardResult> results(hog_list.size());
        std::atomic<size_t> next_h{0};
        std::atomic<bool>   sfailed{false};
        std::string         serr;
        std::mutex          serr_mtx;

        auto sworker = [&]() {
            for (;;) {
                size_t hi = next_h.fetch_add(1, std::memory_order_relaxed);
                if (hi >= hog_list.size() || sfailed.load(std::memory_order_relaxed)) break;
                try {
                    results[hi] = lhi::merge_shard_compute_extents(
                        hog_list[hi]->first, hog_list[hi]->second, tfd_ints, accessions);
                } catch (const std::exception& e) {
                    std::lock_guard<std::mutex> lk(serr_mtx);
                    if (!sfailed.exchange(true)) serr = e.what();
                    break;
                }
            }
        };

        {
            std::vector<std::thread> pool;
            pool.reserve(nt - 1);
            for (size_t t = 1; t < nt; ++t) pool.emplace_back(sworker);
            sworker();
            for (auto& t : pool) t.join();
        }
        if (sfailed.load()) { std::cerr << "error: " << serr << "\n"; return 1; }

        // Serial write pass (main thread); write in completion order.
        lhi::UniqueFd fd_out(open(out_lhg.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644));
        if (fd_out < 0) { std::cerr << "cannot create: " << out_lhg << "\n"; return 1; }
        uint8_t fhdr[lhi::LHG_HEADER_SZ] = {};
        memcpy(fhdr, lhi::LHG_FILE_MAGIC, 4); fhdr[4] = lhi::LHG_VERSION;
        if (::write(fd_out, fhdr, lhi::LHG_HEADER_SZ) != ssize_t(lhi::LHG_HEADER_SZ)) {
            std::cerr << "write header failed\n"; return 1;
        }

        lhi::GlobalIndex gi;
        gi.accessions = accessions;
        lhi::WriteBuffer wb;
        uint64_t data_pos = lhi::LHG_HEADER_SZ;

        for (auto& r : results) {
            uint8_t hdr8[8];
            memcpy(hdr8, lhi::LHG_HOG_ENTRY_MAGIC, 4);
            for (int i = 0; i < 4; ++i) hdr8[4+i] = uint8_t(r.stored_sz >> (8*i));

            lhi::HogIndexEntry e;
            e.hog_id       = std::move(r.hog_id);
            e.data_offset  = data_pos;
            e.data_length  = uint64_t(8) + r.payload.size();
            e.n_accessions = r.n_accs;
            gi.entries.push_back(std::move(e));

            wb.append(fd_out, out_lhg, hdr8, 8);
            wb.append(fd_out, out_lhg, r.payload.data(), r.payload.size());
            data_pos += uint64_t(8) + r.payload.size();
            { std::vector<uint8_t>().swap(r.payload); }
        }
        wb.flush_to(fd_out, out_lhg);

        // Sort index entries by hog_id (physical order = completion order, may differ)
        std::sort(gi.entries.begin(), gi.entries.end(),
                  [](const lhi::HogIndexEntry& a, const lhi::HogIndexEntry& b) {
                      return a.hog_id < b.hog_id;
                  });

        uint64_t index_offset = data_pos;
        auto idx_bytes = lhi::build_index_bytes(gi.entries, gi.accessions);
        { lhi::WriteBuffer w2; w2.append(fd_out, out_lhg, idx_bytes.data(), idx_bytes.size()); w2.flush_to(fd_out, out_lhg); }

        uint8_t final_hdr[lhi::LHG_HEADER_SZ] = {};
        memcpy(final_hdr, lhi::LHG_FILE_MAGIC, 4); final_hdr[4] = lhi::LHG_VERSION;
        for (int i = 0; i < 8; ++i) final_hdr[8+i] = uint8_t(index_offset >> (8*i));
        if (::pwrite(fd_out, final_hdr, lhi::LHG_HEADER_SZ, 0) != ssize_t(lhi::LHG_HEADER_SZ)) {
            std::cerr << "pwrite header failed\n"; return 1;
        }
        lhi::UniqueFd fi(open(out_lhgi.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644));
        if (fi < 0) { std::cerr << "cannot create: " << out_lhgi << "\n"; return 1; }
        { lhi::WriteBuffer w3; w3.append(fi, out_lhgi, idx_bytes.data(), idx_bytes.size()); w3.flush_to(fi, out_lhgi); }

        std::cerr << "merge-shard done: " << gi.entries.size() << " HOGs, "
                  << gi.accessions.size() << " accessions → " << out_lhg << " (-t " << nt << ")\n";
        return 0;
    }

    if (mode == "saav") {
        if (argc < 6) { usage(argv[0]); return 1; }
        std::string lhg_path  = argv[2];
        std::string lhgi_path = argv[3];
        std::string hog_id    = argv[4];
        uint32_t    pos       = uint32_t(std::stoul(argv[5]));
        std::optional<uint8_t> aa;
        uint8_t min_pident = 0;
        for (int i = 6; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--min-pident" && i+1 < argc) min_pident = uint8_t(std::stoul(argv[++i]));
            else if (aa == std::nullopt && a.size() == 1) aa = lhi::encode_aa(a[0]);
        }

        lhi::GlobalIndex idx;
        if (!idx.load(lhgi_path) && !idx.load_from_lhg(lhg_path)) {
            std::cerr << "cannot load index from " << lhgi_path << "\n";
            return 1;
        }
        lhi::query_saav(lhg_path, idx, hog_id, pos, aa, min_pident);
        return 0;
    }

    if (mode == "freq") {
        if (argc < 5) { usage(argv[0]); return 1; }
        std::string lhg_path  = argv[2];
        std::string lhgi_path = argv[3];
        std::string hog_id    = argv[4];
        uint8_t min_pident = 0;
        for (int i = 5; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "--min-pident" && i+1 < argc) min_pident = uint8_t(std::stoul(argv[++i]));
        }

        lhi::GlobalIndex idx;
        if (!idx.load(lhgi_path) && !idx.load_from_lhg(lhg_path)) {
            std::cerr << "cannot load index from " << lhgi_path << "\n";
            return 1;
        }
        lhi::query_freq(lhg_path, idx, hog_id, min_pident);
        return 0;
    }

    if (mode == "stat") {
        if (argc < 3) { usage(argv[0]); return 1; }
        std::string path = argv[2];

        // If .lhgi provided or path ends with .lhg → global stats.
        if (argc >= 4 || (path.size() > 4 && path.substr(path.size()-4) == ".lhg")) {
            std::string lhgi_path = (argc >= 4) ? argv[3]
                                  : path.substr(0, path.size()-4) + ".lhgi";
            lhi::GlobalIndex idx;
            bool ok = idx.load(lhgi_path);
            if (!ok) ok = idx.load_from_lhg(path);
            if (!ok) { std::cerr << "cannot load index from " << lhgi_path << "\n"; return 1; }

            struct stat st;
            if (::stat(path.c_str(), &st) == 0)
                std::fprintf(stdout, "file_size_mb\t%.1f\n", double(st.st_size) / (1024*1024));
            std::fprintf(stdout, "n_hogs\t%zu\n", idx.entries.size());
            std::fprintf(stdout, "n_accessions\t%zu\n", idx.accessions.size());

            uint64_t total_pairs = 0;
            for (auto& e : idx.entries) total_pairs += e.n_accessions;
            std::fprintf(stdout, "total_hog_accession_pairs\t%llu\n",
                         (unsigned long long)total_pairs);

            std::printf("hog_id\tn_accessions\tdata_length_kb\n");
            for (auto& e : idx.entries)
                std::printf("%s\t%u\t%.1f\n", e.hog_id.c_str(),
                            e.n_accessions,
                            double(e.data_length) / 1024.0);
        } else {
            // .lhb stats
            size_t n_blocks = 0;
            std::string acc;
            std::unordered_map<std::string, size_t> hog_blocks;
            lhi::scan_batch_file(path, 0,
                [&](const std::string& a) { acc = a; },
                [&](lhi::BatchBlockRef ref) {
                    ++n_blocks; hog_blocks[ref.hog_id]++;
                });
            std::printf("acc_id\t%s\n", acc.c_str());
            std::printf("n_hogs\t%zu\n", hog_blocks.size());
            std::printf("n_blocks\t%zu\n", n_blocks);
            std::printf("hog_id\tn_blocks\n");
            std::vector<std::pair<std::string,size_t>> hv(hog_blocks.begin(), hog_blocks.end());
            std::sort(hv.begin(), hv.end());
            for (auto& [h, n] : hv) std::printf("%s\t%zu\n", h.c_str(), n);
        }
        return 0;
    }

    usage(argv[0]);
    return 1;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
