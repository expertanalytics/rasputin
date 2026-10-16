#pragma once

#include <terrain/hydrology/flood.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace terrain::hydrology {

// How many nodes drain through each node, from the flood `upstream` runs
// (docs/increments/29-nve-reference-catchments.md, "Accumulation, from the
// same flood (C++)"). For every node c with data, count[c] and the two bits of
// reach[c] equal upstream(z, {c})'s nodes_in, touches_edge and touches_nodata.
struct AccumulateOutcome {
    std::vector<std::uint32_t> count; // row-major; 0 on NoData, else the number of
                                      // nodes that drain through the node, itself included
    std::vector<std::uint8_t> reach;  // bit 0: that catchment touches the window's edge,
                                      // bit 1: it touches NoData (upstream's two flags)
    std::vector<std::uint8_t> flow_to; // the neighbour the node drains to (its flooder), as
                                       // 3*(dr+1)+(dc+1) with dr, dc in -1..1 (4 is never
                                       // used); 255 for an outlet and on NoData
};

inline constexpr std::uint8_t kFlowOutlet = 255;

// Refused with std::length_error at 2^32 nodes or more (count and order are
// 32-bit), before any cell is read or any per-node array allocated.
template <raster::RasterSource R>
[[nodiscard]] AccumulateOutcome accumulate(const R& z) {
    using namespace detail;
    const raster::RasterGeometry& g = z.geometry();
    const std::size_t rows = g.rows(), cols = g.cols(), n = g.size();
    if (rows != 0 && cols > std::size_t{0xFFFFFFFF} / rows) // rows * cols >= 2^32, unwrapped
        throw std::length_error("accumulate: the raster has 2^32 nodes or more");

    AccumulateOutcome out;
    out.count.assign(n, 0);
    out.flow_to.assign(n, kFlowOutlet);
    std::vector<std::uint32_t> order; // push order: every node after its flooder
    order.reserve(n);
    // The flood's state lives in `reach` until the bits replace it.
    const std::vector<std::size_t> beside_nodata = flood(z, out.reach, [&](std::size_t i, std::size_t j) {
        out.count[j] = 1;
        order.push_back(static_cast<std::uint32_t>(j));
        if (i == j)
            return;
        const auto dr = static_cast<std::ptrdiff_t>(i / cols) - static_cast<std::ptrdiff_t>(j / cols);
        const auto dc = static_cast<std::ptrdiff_t>(i % cols) - static_cast<std::ptrdiff_t>(j % cols);
        out.flow_to[j] = static_cast<std::uint8_t>(3 * (dr + 1) + (dc + 1));
    });

    // The bits a seed at each node would raise on its own: bit 0 within one
    // node of the edge, bit 1 an outlet beside NoData or an 8-neighbour of one.
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t r = i / cols, c = i % cols;
        const bool near_edge = r <= 1 || c <= 1 || r + 2 >= rows || c + 2 >= cols;
        out.reach[i] = out.count[i] != 0 && near_edge ? 1 : 0;
    }
    for (const std::size_t o : beside_nodata) {
        out.reach[o] |= 2U;
        each_neighbour(o, rows, cols, [&](std::size_t j) {
            if (out.count[j] != 0)
                out.reach[j] |= 2U;
        });
    }

    // Reverse push order visits every node before its flooder.
    for (auto k = order.size(); k-- > 0;) {
        const std::size_t j = order[k];
        const std::uint8_t to = out.flow_to[j];
        if (to == kFlowOutlet)
            continue;
        const std::size_t p = (j / cols + to / 3 - 1) * cols + (j % cols + to % 3 - 1);
        out.count[p] += out.count[j];
        out.reach[p] |= out.reach[j];
    }
    return out;
}

} // namespace terrain::hydrology
