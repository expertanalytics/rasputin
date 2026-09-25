#pragma once

// Contiguous-chunk parallel for over std::jthread
// (docs/increments/14-adaptive-refinement.md, R7).
//
// fn(begin, end) is called once per non-empty chunk [begin, end) of [0, n),
// each on its own thread, and every thread is joined before this returns.
// Nothing here is shared between chunks: whether fn is race-free is fn's
// business, and the refinement scan makes it so by writing only its own slots.
//
// If fn throws, every chunk still runs to its end and is joined, then the
// exception of the lowest-index chunk that threw is rethrown; the others are
// dropped. Which one surfaces is fixed by the chunking, not by thread timing,
// and the single-chunk path propagates its throw the same way.
//
// No pool. Threads are created per call and joined at its end, so no state
// outlives a call. threads == 0 means hardware_concurrency, or 1 if that is 0.
// With fewer items than threads, fewer threads start: a chunk is never empty.

#include <algorithm>
#include <cstddef>
#include <exception>
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
    // Declared before the workers so it outlives their join, even when
    // starting a later thread throws. Each slot is written by its own chunk.
    std::vector<std::exception_ptr> errors(chunks);
    {
        std::vector<std::jthread> workers;
        workers.reserve(chunks);
        std::size_t begin = 0;
        for (std::size_t k = 0; k < chunks; ++k) {
            const std::size_t end = begin + base + (k < extra ? 1 : 0);
            workers.emplace_back([&fn, &error = errors[k], begin, end] {
                try {
                    fn(begin, end);
                } catch (...) {
                    error = std::current_exception();
                }
            });
            begin = end;
        }
    }  // the jthreads join here
    for (const std::exception_ptr& error : errors)
        if (error)
            std::rethrow_exception(error);
}

}  // namespace terrain::parallel_util
