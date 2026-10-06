#pragma once

// Dynamic-block parallel for over std::jthread, the one the refinement scan
// uses (docs/increments/21-parallel-refine.md, QW1; the static-chunk variant of
// docs/increments/14-adaptive-refinement.md, R7, it replaced is gone).
//
// [0, n) is cut into blocks of `block` items, block k being
// [k*b, min(n, (k+1)*b)), and workers take block indices from one shared atomic
// counter until it runs out, so a slow thread takes fewer blocks. fn(begin,
// end) is called exactly once per block, and every thread is joined before this
// returns. The counter is the only shared write: whether fn is race-free is
// fn's business, and the refinement scan makes it so by writing only its own
// slots; which thread runs a block cannot change what fn writes, so 14 R7's "no
// shared result writes" holds.
//
// No pool. Threads are created per call and joined at its end, so no state
// outlives a call. threads == 0 means hardware_concurrency, or 1 if that is 0.
// With one thread or n < inline_below every block runs on the calling thread in
// ascending order, which spares small rounds the thread starts. Every block
// runs even when some throw, and the exception of the lowest block index that
// threw is rethrown, fixed by the partition, not by timing.

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <thread>
#include <vector>

namespace terrain::parallel_util {

struct BlockSchedule {
    std::size_t block = 0;  // items per block; 0: default_block(n, threads)
    // Rounds of fewer items run inline. From a thread sweep on the 1 m
    // benchmark (docs/increments/21-parallel-refine.md, "21a: inline_below").
    std::size_t inline_below = 256;
};

[[nodiscard]] constexpr std::size_t default_block(std::size_t n, unsigned threads) noexcept {
    return std::max<std::size_t>(1, n / (std::size_t{16} * threads));
}

template <typename Fn>
void for_each_block(std::size_t n, unsigned threads, BlockSchedule schedule, Fn&& fn) {
    if (threads == 0)
        threads = std::max(1u, std::thread::hardware_concurrency());
    const std::size_t b = schedule.block == 0 ? default_block(n, threads) : schedule.block;
    const std::size_t blocks = n / b + (n % b != 0 ? 1 : 0);
    const auto run = [&fn, n, b](std::size_t k) {
        const std::size_t begin = k * b;
        fn(begin, begin + std::min(b, n - begin));
    };
    const std::size_t workers = std::min<std::size_t>(threads, blocks);
    if (workers <= 1 || n < schedule.inline_below) {
        std::exception_ptr first;  // ascending, so the first to throw is the lowest
        for (std::size_t k = 0; k < blocks; ++k) {
            try {
                run(k);
            } catch (...) {
                if (!first)
                    first = std::current_exception();
            }
        }
        if (first)
            std::rethrow_exception(first);
        return;
    }
    // Declared before the workers so they outlive the join; slot k is written
    // only by the thread that took block k.
    std::vector<std::exception_ptr> errors(blocks);
    std::atomic<std::size_t> next{0};
    {
        std::vector<std::jthread> team;
        team.reserve(workers);
        for (std::size_t w = 0; w < workers; ++w)
            team.emplace_back([&] {
                // Relaxed: the counter orders nothing; the join publishes fn's writes.
                for (std::size_t k; (k = next.fetch_add(1, std::memory_order_relaxed)) < blocks;) {
                    try {
                        run(k);
                    } catch (...) {
                        errors[k] = std::current_exception();
                    }
                }
            });
    }  // the jthreads join here
    for (const std::exception_ptr& error : errors)
        if (error)
            std::rethrow_exception(error);
}

}  // namespace terrain::parallel_util
