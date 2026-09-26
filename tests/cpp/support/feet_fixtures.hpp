#pragma once

// Test-only fixtures for increment 20b (docs/increments/20b-min-insertion-distance.md):
// a domain whose left side runs at a small angle to a lattice column, so the
// nodes beside it sit a fraction of a cell from a 146 m constraint segment,
// and the DEMs to refine it against. Shared by
// property/prop_refinement_constraint_feet.cpp and by the scratch program that
// recorded F5's increment-20 digests, so the recorded input and the tested
// input are one definition.
//
// Geometry: refinement_fixtures::geometry (dx = 10 m, dy = 5 m, non-square so
// a row/col swap in the world distance cannot pass). cap = 2.5 m and
// floor = 0.05 m there (R3).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include "quality_fixtures.hpp"
#include "refinement_fixtures.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace feet_fixtures {

inline constexpr std::size_t kN = 33;

// The needle node (row, col) and how far, in columns, the left side passes
// to its left at that row: 0.0035 cells, 3.5 cm in world (F2).
inline constexpr std::int64_t kNeedleRow = 16;
inline constexpr std::int64_t kNeedleCol = 8;
inline constexpr double kNeedleOffset = 0.0035;

// Columns gained per row along the left side: 1 / 81.3, about 1.4 degrees in
// world, so nodes on column 8 below the needle sit 0.0158, 0.0281, ... 0.18
// cells from it and none lies on it.
inline constexpr double kSlope = 1.0 / 81.3;

// The quadrilateral A B C D in (col, row), counter-clockwise in world, with
// A-B the tilted left side. Masks 1, 2, 4, 8 on A-B, B-C, C-D, D-A, so a mask
// carried to the wrong piece shows.
inline std::vector<std::array<double, 2>> needle_ring(double slope = kSlope) {
    const double c16 = static_cast<double>(kNeedleCol) - kNeedleOffset;
    const double ra = 1.37, rb = 30.61;
    return {{c16 - slope * (16.0 - ra), ra},
            {c16 + slope * (rb - 16.0), rb},
            {29.3, 30.2},
            {28.7, 1.9}};
}

inline quality_fixtures::Start needle_start(const terrain::raster::RasterGeometry& g,
                                            double slope = kSlope) {
    auto s = quality_fixtures::fan(g, needle_ring(slope));
    // fan() stores edge i as (min, max) of (i, i + 1); give each its side's mask.
    const std::array<std::uint32_t, 4> side{1u, 2u, 4u, 8u};
    for (std::size_t i = 0; i < 4; ++i) s.masks[i] = side[i];
    return s;
}

// A gentle plane plus a bump centred on the needle node: the node is high,
// the side's linear z low, so refinement must put a vertex beside the segment.
// `steep` scales the plane: 0 leaves ε at the cap, 1 makes it slope-bound.
inline terrain::raster::Raster<float> needle_dem(double bump = 6.0, double steep = 0.0) {
    const auto g = refinement_fixtures::geometry(kN, kN);
    std::vector<float> z(kN * kN);
    for (std::size_t r = 0; r < kN; ++r)
        for (std::size_t c = 0; c < kN; ++c) {
            const double dr = static_cast<double>(r) - static_cast<double>(kNeedleRow);
            const double dc = static_cast<double>(c) - static_cast<double>(kNeedleCol);
            const double x = static_cast<double>(c) * g.delta_x(), y = static_cast<double>(r) * g.delta_y();
            z[r * kN + c] = static_cast<float>(0.01 * x + steep * (0.6 * x + 0.35 * y)
                                               + bump * std::exp(-(dr * dr + dc * dc) / 4.5));
        }
    return terrain::raster::Raster<float>{g, std::move(z)};
}

}  // namespace feet_fixtures
