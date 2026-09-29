#pragma once

#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <queue>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace terrain::hydrology {

// The catchment of a seed set: every DEM node whose drainage path over the
// filled surface passes through a seed (docs/increments/22-auto-catchment.md,
// "Flow and membership: one flood"; Barnes, Lehman and Mulla 2014a,
// Priority-Flood, labelling one catchment).
struct UpstreamOutcome {
    std::vector<std::uint8_t> mask; // row-major, 1 in, 0 otherwise (out, NoData)
    std::size_t nodes_in{};
    std::size_t row_min{}, row_max{}, col_min{}, col_max{}; // inclusive; nodes_in > 0
    // Where the catchment may continue beyond what the window knows: an
    // in-node that is an edge outlet or an 8-neighbour of one, and likewise
    // for an outlet beside NoData ("The flags").
    bool touches_edge{};
    bool touches_nodata{};
};

namespace detail {

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

} // namespace detail

// `seed` is row-major, one byte per node, non-zero for a seed; a NoData seed
// is ignored. Keys are (level, push counter): equal levels pop first in,
// first out, so the result depends on nothing but the input.
template <raster::RasterSource R>
[[nodiscard]] UpstreamOutcome upstream(const R& z, std::span<const std::uint8_t> seed) {
    using namespace detail;
    const raster::RasterGeometry& g = z.geometry();
    const std::size_t rows = g.rows(), cols = g.cols(), n = g.size();
    if (seed.size() != n)
        throw std::invalid_argument("upstream: the seed mask's size is not the raster's");

    UpstreamOutcome out;
    std::vector<std::uint8_t>& state = out.mask;
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

    // Outlets: every valid node on the window's edge or beside NoData.
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
            state[i] = seed[i] != 0 ? kIn : kOut;
            queue.emplace(level_of(i), counter++, i);
        }
    }

    while (!queue.empty()) {
        const auto [level, order, i] = queue.top();
        queue.pop();
        const std::uint8_t label = state[i];
        each_neighbour(i, rows, cols, [&](std::size_t j) {
            if (state[j] != kUnreached)
                return;
            state[j] = seed[j] != 0 ? static_cast<std::uint8_t>(kIn) : label;
            const double zj = level_of(j);
            queue.emplace(zj > level ? zj : level, counter++, j);
        });
    }

    out.row_min = rows;
    out.col_min = cols;
    for (std::size_t i = 0; i < n; ++i) {
        if (state[i] != kIn)
            continue;
        const std::size_t r = i / cols, c = i % cols;
        ++out.nodes_in;
        out.row_min = r < out.row_min ? r : out.row_min;
        out.row_max = r > out.row_max ? r : out.row_max;
        out.col_min = c < out.col_min ? c : out.col_min;
        out.col_max = c > out.col_max ? c : out.col_max;
    }
    if (out.nodes_in == 0) {
        out.row_min = out.col_min = 0;
    } else {
        out.touches_edge = out.row_min <= 1 || out.col_min <= 1 || out.row_max + 2 >= rows
                           || out.col_max + 2 >= cols;
    }
    for (const std::size_t o : beside_nodata) {
        bool near = state[o] == kIn;
        each_neighbour(o, rows, cols, [&](std::size_t j) { near |= state[j] == kIn; });
        out.touches_nodata = out.touches_nodata || near;
    }
    for (auto& s : state)
        s = s == kIn ? 1 : 0;
    return out;
}

} // namespace terrain::hydrology
