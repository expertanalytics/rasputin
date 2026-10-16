#pragma once

#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <tuple>
#include <vector>

// The Priority-Flood both `upstream` and `accumulate` run, in one place so the
// two cannot drift (docs/increments/29-nve-reference-catchments.md,
// "Accumulation, from the same flood (C++)"; Barnes, Lehman and Mulla 2014a).
namespace terrain::hydrology::detail {

enum : std::uint8_t { kUnreached = 0, kIn = 1, kOut = 2, kNoData = 3 };

// Calls f(j) for each 8-neighbour j of node i, in a fixed order: the order is
// part of the result, since first in, first out breaks ties.
template <typename F>
void each_neighbour(std::size_t i, std::size_t rows, std::size_t cols, F&& f) {
    const std::size_t r = i / cols, c = i % cols;
    const std::size_t r0 = r > 0 ? r - 1 : 0, r1 = r + 1 < rows ? r + 1 : r;
    const std::size_t c0 = c > 0 ? c - 1 : 0, c1 = c + 1 < cols ? c + 1 : c;
    for (std::size_t rr = r0; rr <= r1; ++rr)
        for (std::size_t cc = c0; cc <= c1; ++cc)
            if (rr != r || cc != c)
                f(rr * cols + cc);
}

// Floods z from its outlets (every valid node on the window's edge or beside
// NoData). `state` is resized to the raster, kNoData on NoData; a reached node
// is set to kOut before on_reach(i, j) is called, which may relabel it (never
// to kUnreached). on_reach(j, j) announces outlet j; on_reach(i, j) says popped
// node i reached node j first. Keys are (level, push counter): equal levels
// pop first in, first out, so the result depends on nothing but the input.
// Returns the outlets beside NoData, ascending.
template <raster::RasterSource R, typename F>
std::vector<std::size_t> flood(const R& z, std::vector<std::uint8_t>& state, F&& on_reach) {
    const raster::RasterGeometry& g = z.geometry();
    const std::size_t rows = g.rows(), cols = g.cols(), n = g.size();
    state.assign(n, kUnreached);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c)
            if (z.is_nodata(raster::CellIndex{r, c}))
                state[r * cols + c] = kNoData;

    using Entry = std::tuple<double, std::uint64_t, std::size_t>;
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> queue;
    std::uint64_t counter = 0;
    const auto level_of = [&](std::size_t i) {
        return static_cast<double>(z.value_at(raster::CellIndex{i / cols, i % cols}));
    };

    std::vector<std::size_t> beside_nodata;
    for (std::size_t i = 0; i < n; ++i) {
        if (state[i] == kNoData)
            continue;
        const std::size_t r = i / cols, c = i % cols;
        bool nodata_near = false;
        each_neighbour(i, rows, cols, [&](std::size_t j) { nodata_near |= state[j] == kNoData; });
        if (nodata_near)
            beside_nodata.push_back(i);
        if (nodata_near || r == 0 || c == 0 || r + 1 == rows || c + 1 == cols) {
            state[i] = kOut;
            on_reach(i, i);
            queue.emplace(level_of(i), counter++, i);
        }
    }

    while (!queue.empty()) {
        const auto [level, order, i] = queue.top();
        queue.pop();
        each_neighbour(i, rows, cols, [&](std::size_t j) {
            if (state[j] != kUnreached)
                return;
            state[j] = kOut;
            on_reach(i, j);
            const double zj = level_of(j);
            queue.emplace(zj > level ? zj : level, counter++, j);
        });
    }
    return beside_nodata;
}

} // namespace terrain::hydrology::detail
