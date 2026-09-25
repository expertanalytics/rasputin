// Increment 14, T11 (docs/increments/14-adaptive-refinement.md, R7):
// for_each_chunk over std::jthread.
//
// Interface: include/terrain/parallel_util/chunks.hpp.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/parallel_util/chunks.hpp>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

using terrain::parallel_util::for_each_chunk;

namespace {

struct Visit {
    std::vector<std::atomic<int>> hits;
    std::vector<std::pair<std::size_t, std::size_t>> chunks;
    std::mutex m;
    explicit Visit(std::size_t n) : hits(n) {}
};

void run(Visit& v, std::size_t n, unsigned threads) {
    for_each_chunk(n, threads, [&v](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) v.hits[i].fetch_add(1);
        const std::lock_guard lock{v.m};
        v.chunks.emplace_back(begin, end);
    });
}

}  // namespace

TEST_CASE("for_each_chunk visits every index exactly once, in contiguous non-empty chunks",
          "[refinement][chunks]") {
    const std::size_t n = GENERATE(std::size_t{1}, 3, 7, 8, 100);
    const unsigned threads = GENERATE(1u, 2u, 7u, 8u, 0u);
    CAPTURE(n, threads);
    Visit v{n};
    run(v, n, threads);

    for (std::size_t i = 0; i < n; ++i) {
        CAPTURE(i);
        REQUIRE(v.hits[i].load() == 1);
    }
    std::sort(v.chunks.begin(), v.chunks.end());
    REQUIRE(v.chunks.front().first == 0);
    REQUIRE(v.chunks.back().second == n);
    for (std::size_t k = 0; k < v.chunks.size(); ++k) {
        REQUIRE(v.chunks[k].first < v.chunks[k].second);  // never empty
        if (k > 0) REQUIRE(v.chunks[k].first == v.chunks[k - 1].second);
    }
    const unsigned effective =
        threads == 0 ? std::max(1u, std::thread::hardware_concurrency()) : threads;
    REQUIRE(v.chunks.size() <= std::min<std::size_t>(n, effective));
}

TEST_CASE("for_each_chunk with n == 0 never calls fn", "[refinement][chunks]") {
    const unsigned threads = GENERATE(1u, 2u, 0u);
    int calls = 0;
    for_each_chunk(0, threads, [&calls](std::size_t, std::size_t) { ++calls; });
    REQUIRE(calls == 0);
}

TEST_CASE("for_each_chunk with one thread makes a single chunk", "[refinement][chunks]") {
    // threads == 1 is the serial reference T6 compares against.
    Visit v{10};
    run(v, 10, 1);
    REQUIRE(v.chunks.size() == 1);
    REQUIRE(v.chunks.front() == std::pair<std::size_t, std::size_t>{0, 10});
}
