// SCRATCH micro-benchmark, 21c "A0 and the scaling model" (not production code).
// Measures, for T threads (argv[1], default 8):
//   spawn   -- start T std::jthread with an empty body and join them (what one
//              for_each_chunk / for_each_block call pays in thread management);
//   feb     -- parallel_util::for_each_block over T*256 items, a trivial body
//              (the project's own spawn path, inline_below respected);
//   barrier -- a persistent team of T threads, K phases of std::barrier::arrive_and_wait;
//   spin    -- the same with a sense-reversing spin barrier (atomic counter,
//              spin with a yield every 64 polls).
// Each is repeated in R samples; the median and min/max per event are printed.
// Build: c++ -std=c++20 -O3 -DNDEBUG -I include sync_bench.cpp -o sync_bench
#include <terrain/parallel_util/chunks.hpp>

#include <algorithm>
#include <atomic>
#include <barrier>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>

using clk = std::chrono::steady_clock;
static double us(clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::micro>(b - a).count(); }

static void report(const char* name, std::vector<double> v) {
    std::sort(v.begin(), v.end());
    std::printf("%-8s median %9.3f us  min %9.3f  max %9.3f  (n=%zu)\n", name, v[v.size() / 2], v.front(), v.back(), v.size());
}

struct SpinBarrier {
    explicit SpinBarrier(unsigned n) : n_(n) {}
    void wait(bool& local_sense) {
        local_sense = !local_sense;
        if (count_.fetch_add(1, std::memory_order_acq_rel) == n_ - 1) {
            count_.store(0, std::memory_order_relaxed);
            sense_.store(local_sense, std::memory_order_release);
        } else {
            unsigned k = 0;
            while (sense_.load(std::memory_order_acquire) != local_sense)
                if (++k % 64 == 0)
                    std::this_thread::yield();
        }
    }
    unsigned n_;
    std::atomic<unsigned> count_{0};
    std::atomic<bool> sense_{false};
};

int main(int argc, char** argv) {
    const unsigned T = argc > 1 ? static_cast<unsigned>(std::atoi(argv[1])) : 8;
    constexpr int R = 200, K = 2000, S = 15;
    std::atomic<long> sink{0};
    std::printf("threads %u\n", T);
    {
        std::vector<double> v;
        for (int r = 0; r < R; ++r) {
            const auto a = clk::now();
            {
                std::vector<std::jthread> w;
                w.reserve(T);
                for (unsigned i = 0; i < T; ++i)
                    w.emplace_back([&] { sink.fetch_add(1, std::memory_order_relaxed); });
            }
            v.push_back(us(a, clk::now()));
        }
        report("spawn", v);
    }
    {
        std::vector<double> v;
        std::vector<int> x(T * 256, 1);
        for (int r = 0; r < R; ++r) {
            const auto a = clk::now();
            terrain::parallel_util::for_each_block(x.size(), T, terrain::parallel_util::BlockSchedule{},
                                                   [&](std::size_t b, std::size_t e) {
                                                       long s = 0;
                                                       for (auto i = b; i < e; ++i)
                                                           s += x[i];
                                                       sink.fetch_add(s, std::memory_order_relaxed);
                                                   });
            v.push_back(us(a, clk::now()));
        }
        report("feb", v);
    }
    {
        std::vector<double> v;
        for (int s = 0; s < S; ++s) {
            std::barrier bar(static_cast<std::ptrdiff_t>(T));
            clk::time_point a, b;
            {
                std::vector<std::jthread> w;
                for (unsigned i = 0; i < T; ++i)
                    w.emplace_back([&, i] {
                        bar.arrive_and_wait();
                        if (i == 0)
                            a = clk::now();
                        for (int k = 0; k < K; ++k)
                            bar.arrive_and_wait();
                        if (i == 0)
                            b = clk::now();
                    });
            }
            v.push_back(us(a, b) / K);
        }
        report("barrier", v);
    }
    {
        std::vector<double> v;
        for (int s = 0; s < S; ++s) {
            SpinBarrier bar(T);
            clk::time_point a, b;
            {
                std::vector<std::jthread> w;
                for (unsigned i = 0; i < T; ++i)
                    w.emplace_back([&, i] {
                        bool sense = false;
                        bar.wait(sense);
                        if (i == 0)
                            a = clk::now();
                        for (int k = 0; k < K; ++k)
                            bar.wait(sense);
                        if (i == 0)
                            b = clk::now();
                    });
            }
            v.push_back(us(a, b) / K);
        }
        report("spin", v);
    }
    return sink.load() == -1;
}
