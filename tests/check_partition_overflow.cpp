// Regression test: partition index round-trips correctly with >65535 extents per HOG.
// Previously n_extents was stored as uint16 — this would silently truncate and corrupt.
#include <iostream>
#include "../src/partition.hpp"
#include <cassert>
#include <cstdio>

using namespace lhi;

int main() {
    const char* path = "/tmp/hoki_test_overflow.idx";

    const uint32_t N = 70000; // > UINT16_MAX
    std::map<std::string, std::vector<PartitionIndexExtent>> idx;
    auto& exts = idx["N0.HOG0000001"];
    exts.reserve(N);
    for (uint32_t i = 0; i < N; ++i)
        exts.push_back({0u, uint64_t(i) * 128, 128u});

    write_partition_index(idx, 1u, path);

    uint32_t n_threads = 0;
    auto got = load_partition_index(path, n_threads);

    assert(got.count("N0.HOG0000001"));
    assert(got["N0.HOG0000001"].size() == N);
    for (uint32_t i = 0; i < N; ++i) {
        assert(got["N0.HOG0000001"][i].entry_offset == uint64_t(i) * 128);
        assert(got["N0.HOG0000001"][i].entry_len    == 128u);
    }

    ::remove(path);
    std::puts("check_partition_overflow: PASS");
    return 0;
}
