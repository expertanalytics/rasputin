#pragma once

#include <terrain/hydrology/flood.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
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
    const std::vector<std::size_t> beside_nodata = flood(z, state, [&](std::size_t i, std::size_t j) {
        if (seed[j] != 0)
            state[j] = kIn;
        else if (i != j)
            state[j] = state[i];
    });

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
