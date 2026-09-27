// Increment 21a, QW1 (docs/increments/21-parallel-refine.md, section 3, and
// "Pinned by the red suite (21a)"): the dynamic block scheduler that replaces
// for_each_chunk for the scan. A new suite; test_refinement_chunks, which pins
// for_each_chunk, is unchanged.
//
// INVARIANT-CRITICAL (section 7): a dropped block leaves a stale scan result
// and a doubled one is a second writer to a result slot, and either breaks the
// tolerance guarantee silently. So every index is visited exactly once for
// every n, thread count, block size and side of the inline threshold, and the
// exception contract is decided here. Mutation-tested in the red step.
//
// Interface pinned here (include/terrain/parallel_util/chunks.hpp):
//
//   struct BlockSchedule {
//       std::size_t block = 0;         // items per block; 0: default_block(n, threads)
//       std::size_t inline_below = ?;  // n < inline_below runs on the calling thread;
//                                      // the default is 21a's, from the sweep
//   };
//   constexpr std::size_t default_block(std::size_t n, unsigned threads) noexcept;
//       // max(1, n / (16 * threads)), threads >= 1
//   template <class Fn>
//   void for_each_block(std::size_t n, unsigned threads, BlockSchedule, Fn&& fn);
//
// threads == 0 means hardware_concurrency, or 1 if that is 0, as for
// for_each_chunk. With b the block size, block k is [k*b, min(n, (k+1)*b)),
// and fn(begin, end) is called exactly once per block: the set of calls is the
// partition, whatever the threads or the timing. With one thread or
// n < inline_below every call is on the calling thread, in ascending order.
// Otherwise at most min(threads, blocks) threads call fn, and two blocks can
// run at once. Every block runs even when some throw; after the join the
// exception of the lowest block index that threw is rethrown.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/parallel_util/chunks.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <set>
#include <thread>
#include <utility>
#include <vector>

using terrain::parallel_util::BlockSchedule;
using terrain::parallel_util::default_block;
using terrain::parallel_util::for_each_block;

namespace {

using Range = std::pair<std::size_t, std::size_t>;

unsigned resolved(unsigned threads) {
    return threads == 0 ? std::max(1u, std::thread::hardware_concurrency()) : threads;
}

// The partition the contract promises, built independently of the scheduler.
std::vector<Range> partition(std::size_t n, std::size_t block) {
    std::vector<Range> blocks;
    for (std::size_t begin = 0; begin < n; begin += block)
        blocks.emplace_back(begin, std::min(n, begin + block));
    return blocks;
}

// Every call, in the order calls were recorded, with the thread that made it.
struct Trace {
    std::vector<std::atomic<int>> hits;
    std::vector<Range> calls;
    std::vector<std::thread::id> by;
    std::mutex m;
    explicit Trace(std::size_t n) : hits(n) {}

    auto sink() {
        return [this](std::size_t begin, std::size_t end) {
            for (std::size_t i = begin; i < end; ++i) hits[i].fetch_add(1);
            const std::lock_guard lock{m};
            calls.emplace_back(begin, end);
            by.push_back(std::this_thread::get_id());
        };
    }
};

struct BlockError {
    std::size_t block;
};

}  // namespace

TEST_CASE("default_block is n / (16 * threads), at least 1", "[refinement][chunks][dynamic]") {
    STATIC_REQUIRE(default_block(0, 1) == 1);
    STATIC_REQUIRE(default_block(1, 1) == 1);
    STATIC_REQUIRE(default_block(15, 1) == 1);
    STATIC_REQUIRE(default_block(16, 1) == 1);
    STATIC_REQUIRE(default_block(47, 1) == 2);
    STATIC_REQUIRE(default_block(1000, 1) == 62);
    STATIC_REQUIRE(default_block(1280, 8) == 10);
    STATIC_REQUIRE(default_block(1279, 8) == 9);
    STATIC_REQUIRE(default_block(10, 8) == 1);
    STATIC_REQUIRE(default_block(100000, 20) == 312);
}

TEST_CASE("for_each_block calls fn once per block of the partition, every index exactly once",
          "[refinement][chunks][dynamic]") {
    // n covers 0, 1, primes, powers of two and their neighbours; block covers
    // the default (0), 1, blocks that do and do not divide n, and one larger
    // than every n. The threshold is put below, at and above n.
    const std::size_t n = GENERATE(std::size_t{0}, 1, 2, 3, 7, 13, 16, 17, 97, 128, 1009);
    const unsigned threads = GENERATE(1u, 2u, 3u, 7u, 8u, 0u);
    const std::size_t block = GENERATE(std::size_t{0}, 1, 2, 3, 7, 64, 5000);
    const int where = GENERATE(-1, 0, 1);  // threshold at n + where, clamped at 0
    const std::size_t inline_below = where < 0 ? (n == 0 ? 0 : n - 1) : n + static_cast<std::size_t>(where);
    CAPTURE(n, threads, block, inline_below);

    Trace t{n};
    const std::thread::id caller = std::this_thread::get_id();
    for_each_block(n, threads, BlockSchedule{.block = block, .inline_below = inline_below}, t.sink());

    for (std::size_t i = 0; i < n; ++i) {
        CAPTURE(i);
        REQUIRE(t.hits[i].load() == 1);
    }
    const std::size_t b = block == 0 ? default_block(n, resolved(threads)) : block;
    const std::vector<Range> expected = partition(n, b);
    std::vector<Range> sorted = t.calls;
    std::sort(sorted.begin(), sorted.end());
    REQUIRE(sorted == expected);

    const bool inline_run = resolved(threads) == 1 || n < inline_below;
    if (inline_run) {
        REQUIRE(t.calls == expected);  // ascending, on the caller
        for (const auto id : t.by) REQUIRE(id == caller);
    }
    const std::set<std::thread::id> used(t.by.begin(), t.by.end());
    REQUIRE(used.size() <= std::max<std::size_t>(1, std::min<std::size_t>(resolved(threads), expected.size())));
}

TEST_CASE("for_each_block with n == 0 never calls fn", "[refinement][chunks][dynamic]") {
    const unsigned threads = GENERATE(1u, 2u, 0u);
    const std::size_t inline_below = GENERATE(std::size_t{0}, 1, 100);
    int calls = 0;
    for_each_block(0, threads, BlockSchedule{.block = 0, .inline_below = inline_below},
                   [&calls](std::size_t, std::size_t) { ++calls; });
    REQUIRE(calls == 0);
}

TEST_CASE("for_each_block at or above the threshold runs blocks concurrently",
          "[refinement][chunks][dynamic]") {
    // Two blocks, two threads, n == inline_below. Each block waits, up to a
    // generous timeout, for the other to have started, which it can only do
    // if they run on different threads at once. A scheduler that ran the
    // round inline (the threshold ignored, or compared the wrong way) makes
    // each wait time out. A correct one never waits more than a thread start.
    const std::size_t inline_below = GENERATE(std::size_t{0}, 2);
    CAPTURE(inline_below);
    std::mutex m;
    std::condition_variable cv;
    int started = 0;
    std::atomic<int> timed_out{0};
    for_each_block(2, 2, BlockSchedule{.block = 1, .inline_below = inline_below},
                   [&](std::size_t, std::size_t) {
                       std::unique_lock lock{m};
                       ++started;
                       cv.notify_all();
                       if (!cv.wait_for(lock, std::chrono::seconds{10}, [&] { return started >= 2; }))
                           timed_out.fetch_add(1);
                   });
    REQUIRE(started == 2);
    REQUIRE(timed_out.load() == 0);
}

TEST_CASE("for_each_block below the threshold runs on the calling thread even with many threads",
          "[refinement][chunks][dynamic]") {
    const std::size_t n = GENERATE(std::size_t{1}, 5, 99);
    Trace t{n};
    const std::thread::id caller = std::this_thread::get_id();
    for_each_block(n, 8, BlockSchedule{.block = 1, .inline_below = n + 1}, t.sink());
    REQUIRE(t.calls == partition(n, 1));
    for (const auto id : t.by) REQUIRE(id == caller);
}

// The two exception-order cases below make arrival order opposite ways, so
// that neither "keep the first exception to arrive" nor "keep the last" can
// pass both. Arrival order is forced with waits on what the other blocks have
// done, not left to sleeps racing each other; every wait has a timeout, so a
// scheduler that serialises the blocks is slow, never hung.

TEST_CASE("for_each_block rethrows the lowest throwing block's exception after running every block",
          "[refinement][chunks][dynamic]") {
    // 64 items in blocks of 4 is 16 blocks. Every block from `first` on throws
    // its own index. On the threaded path block `first` waits until every
    // higher block is about to throw, then a little longer, so it is the LAST
    // to arrive. Kills: keep-the-first-to-arrive, keep-the-highest-index.
    // (Keep-the-last-to-arrive passes this case; the next one kills it.) The
    // inline paths (one thread; n below the threshold) run in ascending order,
    // so there block `first` arrives first; they must rethrow the same block.
    const unsigned threads = GENERATE(1u, 4u, 8u);
    const std::size_t inline_below = GENERATE(std::size_t{0}, 65);
    const std::size_t first = GENERATE(std::size_t{0}, 1, 5, 15);
    CAPTURE(threads, inline_below, first);
    const bool threaded = threads > 1 && inline_below <= 64;
    const int higher = static_cast<int>(15 - first);
    std::atomic<int> calls{0}, short_blocks{0}, about_to_throw{0}, timed_out{0};
    std::size_t caught = 99;
    try {
        for_each_block(64, threads, BlockSchedule{.block = 4, .inline_below = inline_below},
                       [&](std::size_t begin, std::size_t end) {
                           // No Catch2 assertion off the test thread: count instead.
                           if (end - begin != 4) short_blocks.fetch_add(1);
                           calls.fetch_add(1);
                           const std::size_t k = begin / 4;
                           if (k == first && threaded) {
                               const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds{10};
                               while (about_to_throw.load() < higher) {
                                   if (std::chrono::steady_clock::now() > deadline) {
                                       timed_out.fetch_add(1);
                                       break;
                                   }
                                   std::this_thread::yield();
                               }
                               // Let the higher blocks' throws reach the scheduler.
                               std::this_thread::sleep_for(std::chrono::milliseconds{30});
                           }
                           if (k > first) about_to_throw.fetch_add(1);
                           if (k >= first) throw BlockError{k};
                       });
    } catch (const BlockError& e) {
        caught = e.block;
    }
    REQUIRE(timed_out.load() == 0);
    REQUIRE(caught == first);
    REQUIRE(short_blocks.load() == 0);
    REQUIRE(calls.load() == 16);  // no block is skipped after a throw
}

TEST_CASE("for_each_block rethrows the lowest throwing block's exception when it arrives first",
          "[refinement][chunks][dynamic]") {
    // The reverse arrival order: block `low` throws at once, and block `high`
    // waits until `low` is about to throw, then a little longer, so `low` is
    // the FIRST to arrive and `high` the LAST, on every path (the inline paths
    // run `low` first anyway). Kills: keep-the-last-to-arrive, and
    // keep-the-highest-index. `high` is block 15, the last one handed out, or a
    // middle one, so the slow block is not always the scheduler's last.
    const unsigned threads = GENERATE(1u, 2u, 4u, 8u);
    const std::size_t inline_below = GENERATE(std::size_t{0}, 65);
    const std::size_t low = GENERATE(std::size_t{0}, 5);
    const std::size_t high = GENERATE(std::size_t{9}, 15);
    CAPTURE(threads, inline_below, low, high);
    std::atomic<int> calls{0}, timed_out{0};
    std::atomic<bool> low_thrown{false};
    std::size_t caught = 99;
    try {
        for_each_block(64, threads, BlockSchedule{.block = 4, .inline_below = inline_below},
                       [&](std::size_t begin, std::size_t) {
                           calls.fetch_add(1);
                           const std::size_t k = begin / 4;
                           if (k == low) {
                               low_thrown.store(true);
                               throw BlockError{k};
                           }
                           if (k == high) {
                               const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds{10};
                               while (!low_thrown.load()) {
                                   if (std::chrono::steady_clock::now() > deadline) {
                                       timed_out.fetch_add(1);
                                       break;
                                   }
                                   std::this_thread::yield();
                               }
                               // Let block `low`'s throw reach the scheduler first.
                               std::this_thread::sleep_for(std::chrono::milliseconds{30});
                               throw BlockError{k};
                           }
                       });
    } catch (const BlockError& e) {
        caught = e.block;
    }
    REQUIRE(timed_out.load() == 0);
    REQUIRE(caught == low);
    REQUIRE(calls.load() == 16);
}

TEST_CASE("for_each_block's exception is the same whatever the thread count and threshold",
          "[refinement][chunks][dynamic]") {
    // The throwing set is scattered, not a suffix: blocks 3, 6 and 11 of 13
    // (the last one short). Only block 3's exception may surface.
    const unsigned threads = GENERATE(1u, 2u, 3u, 7u, 0u);
    const std::size_t inline_below = GENERATE(std::size_t{0}, 50, 51);
    CAPTURE(threads, inline_below);
    std::atomic<int> calls{0};
    std::size_t caught = 99;
    try {
        for_each_block(50, threads, BlockSchedule{.block = 4, .inline_below = inline_below},
                       [&](std::size_t begin, std::size_t) {
                           calls.fetch_add(1);
                           const std::size_t k = begin / 4;
                           if (k == 3) std::this_thread::sleep_for(std::chrono::milliseconds{20});
                           if (k == 3 || k == 6 || k == 11) throw BlockError{k};
                       });
    } catch (const BlockError& e) {
        caught = e.block;
    }
    REQUIRE(caught == 3);
    REQUIRE(calls.load() == 13);
}
