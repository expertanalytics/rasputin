#pragma once

// Test-only helpers for increment 14 (docs/increments/14-adaptive-refinement.md).
//
// Everything here is written independently of include/terrain/refinement/ and
// include/terrain/mesh/: the orientation, the node-set membership and the
// plane are re-derived from the design's R2 text, never borrowed from the
// headers under test. The T3 oracle depends on that.
//
// Frame. A lattice node is (row, col). World is x = x_min + col * dx,
// y = y_max - row * dy, so rows grow DOWNWARD. Orientation here is taken in the
// world's handedness, i.e. on (x, y) = (col, -row), so "counter-clockwise"
// means the same thing for a lattice triangle as for the IndexedMesh2 it came
// from. The design does not say which frame R2/R4 mean (see the handback).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <utility>
#include <vector>

namespace refinement_fixtures {

struct RC {
    std::int64_t row;
    std::int64_t col;
    friend constexpr bool operator==(const RC&, const RC&) = default;
    friend constexpr auto operator<=>(const RC&, const RC&) = default;
};

// Twice the signed area of (p, q, r) in the world's handedness.
constexpr std::int64_t orient(RC p, RC q, RC r) noexcept {
    const std::int64_t ax = q.col - p.col, ay = -(q.row - p.row);
    const std::int64_t bx = r.col - p.col, by = -(r.row - p.row);
    return ax * by - ay * bx;
}

// p in the closed triangle (a, b, c), which must be CCW.
constexpr bool in_closed(RC a, RC b, RC c, RC p) noexcept {
    return orient(a, b, p) >= 0 && orient(b, c, p) >= 0 && orient(c, a, p) >= 0;
}

// p strictly inside the open segment (a, b).
constexpr bool on_open_segment(RC a, RC b, RC p) noexcept {
    if (orient(a, b, p) != 0 || p == a || p == b) return false;
    return std::min(a.row, b.row) <= p.row && p.row <= std::max(a.row, b.row)
        && std::min(a.col, b.col) <= p.col && p.col <= std::max(a.col, b.col);
}

// A UTM-scale geometry with non-square cells, so a row/col swap cannot pass.
inline terrain::raster::RasterGeometry geometry(std::size_t rows, std::size_t cols) {
    return terrain::raster::RasterGeometry{500000.0, 7000000.0, 10.0, 5.0, cols, rows};
}

inline RC to_rc(const terrain::raster::RasterGeometry& g, terrain::Point2 p) {
    return RC{static_cast<std::int64_t>(std::llround((g.y_max() - p.y) / g.delta_y())),
              static_cast<std::int64_t>(std::llround((p.x - g.x_min()) / g.delta_x()))};
}

// A start mesh on the lattice, CCW in world, and its constraint edges.
struct StartMesh {
    terrain::IndexedMesh2 mesh;
    std::vector<std::array<std::uint32_t, 2>> edges;  // vertex indices
    std::vector<std::uint32_t> masks;                 // one per edge
    std::vector<RC> lattice;                          // per vertex
};

// Every `stride`-th node plus the last, two triangles per cell with the
// (top-left, bottom-right) diagonal. The outer ring is constrained; its four
// sides carry masks 1, 2, 4, 8 (top, right, bottom, left) so inheritance is
// visible per side.
inline StartMesh grid_mesh(const terrain::raster::RasterGeometry& g, std::size_t stride) {
    auto axis = [stride](std::size_t n) {
        std::vector<std::int64_t> v;
        for (std::size_t i = 0; i < n; i += stride) v.push_back(static_cast<std::int64_t>(i));
        if (v.back() != static_cast<std::int64_t>(n - 1)) v.push_back(static_cast<std::int64_t>(n - 1));
        return v;
    };
    const auto rows = axis(g.rows());
    const auto cols = axis(g.cols());
    const auto nc = static_cast<std::uint32_t>(cols.size());
    const auto nr = static_cast<std::uint32_t>(rows.size());
    StartMesh s;
    std::vector<terrain::Point2> xy;
    for (auto r : rows)
        for (auto c : cols) {
            xy.push_back(g.node({static_cast<std::size_t>(r), static_cast<std::size_t>(c)}));
            s.lattice.push_back(RC{r, c});
        }
    auto at = [nc](std::uint32_t i, std::uint32_t j) { return i * nc + j; };
    std::vector<terrain::TriangleIndices> tris;
    for (std::uint32_t i = 0; i + 1 < nr; ++i)
        for (std::uint32_t j = 0; j + 1 < nc; ++j) {
            const auto tl = at(i, j), tr = at(i, j + 1), bl = at(i + 1, j), br = at(i + 1, j + 1);
            tris.push_back({tl, bl, br});  // CCW in world: down, then right
            tris.push_back({tl, br, tr});
        }
    std::vector<std::uint8_t> constrained(tris.size(), 0);
    s.mesh = terrain::IndexedMesh2{std::move(xy), std::move(tris), std::move(constrained)};
    auto side = [&](std::uint32_t a, std::uint32_t b, std::uint32_t mask) {
        s.edges.push_back({std::min(a, b), std::max(a, b)});
        s.masks.push_back(mask);
    };
    for (std::uint32_t j = 0; j + 1 < nc; ++j) side(at(0, j), at(0, j + 1), 1);
    for (std::uint32_t i = 0; i + 1 < nr; ++i) side(at(i, nc - 1), at(i + 1, nc - 1), 2);
    for (std::uint32_t j = 0; j + 1 < nc; ++j) side(at(nr - 1, j), at(nr - 1, j + 1), 4);
    for (std::uint32_t i = 0; i + 1 < nr; ++i) side(at(i, 0), at(i + 1, 0), 8);
    return s;
}

// Seeded DEMs. mt19937's raw output is specified by the standard; the
// distributions are not, so values are scaled by hand.
inline std::vector<float> rough_dem(std::size_t rows, std::size_t cols, std::uint32_t seed) {
    std::mt19937 gen{seed};
    std::vector<float> v(rows * cols);
    for (auto& x : v) x = static_cast<float>(gen() % 10000u) / 100.0f;  // [0, 100) in 0.01 steps
    return v;
}

inline std::vector<float> smooth_dem(std::size_t rows, std::size_t cols, std::uint32_t seed) {
    std::mt19937 gen{seed};
    const double p1 = static_cast<double>(gen() % 628u) / 100.0;
    const double p2 = static_cast<double>(gen() % 628u) / 100.0;
    std::vector<float> v(rows * cols);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c)
            v[r * cols + c] = static_cast<float>(
                50.0 + 20.0 * std::sin(0.3 * static_cast<double>(c) + p1)
                     + 15.0 * std::cos(0.25 * static_cast<double>(r) + p2));
    return v;
}

}  // namespace refinement_fixtures
