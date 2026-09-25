#pragma once

// Contiguous-chunk parallel for over std::jthread
// (docs/increments/14-adaptive-refinement.md, R7).
//
// fn(begin, end) is called once per non-empty chunk [begin, end) of [0, n),
// each on its own thread, and every thread is joined before this returns.
// Nothing here is shared between chunks: whether fn is race-free is fn's
// business, and the refinement scan makes it so by writing only its own slots.
//
// No pool. Threads are created per call and joined at its end, so no state
// outlives a call. threads == 0 means hardware_concurrency, or 1 if that is 0.
// With fewer items than threads, fewer threads start: a chunk is never empty.

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

namespace terrain::parallel_util {

template <typename Fn>
void for_each_chunk(std::size_t n, unsigned threads, Fn&& fn) {
    if (threads == 0)
        threads = std::max(1u, std::thread::hardware_concurrency());
    const std::size_t chunks = std::min<std::size_t>(n, threads);
    if (chunks <= 1) {
        if (n > 0)
            fn(std::size_t{0}, n);
        return;
    }
    // The first n % chunks chunks take one extra item.
    const std::size_t base = n / chunks, extra = n % chunks;
    std::vector<std::jthread> workers;
    workers.reserve(chunks);
    std::size_t begin = 0;
    for (std::size_t k = 0; k < chunks; ++k) {
        const std::size_t end = begin + base + (k < extra ? 1 : 0);
        workers.emplace_back([&fn, begin, end] { fn(begin, end); });
        begin = end;
    }
}  // the jthreads join here

}  // namespace terrain::parallel_util
