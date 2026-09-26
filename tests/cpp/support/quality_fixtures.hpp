#pragma once

// Test-only fixtures for increment 20 (docs/increments/20-start-quality.md):
// start meshes whose triangles are bad on purpose, and the DEMs to refine them
// against. Shared by prop_refinement_quality.cpp and by the scratch program that
// recorded Q7's increment-18 digests, so the recorded input and the tested
// input are one definition.
//
// cone() and island() are 14b's T12 DEMs, copied from
// property/prop_refinement_refine.cpp, where they sit in an anonymous
// namespace.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <utility>
#include <vector>

namespace quality_fixtures {

using Edges = std::vector<std::array<std::uint32_t, 2>>;

// Integer x_min, y_max, dx and dy: world coordinates of nodes are exact.
inline terrain::raster::RasterGeometry integral_geometry(std::size_t n, double dx, double dy) {
    return terrain::raster::RasterGeometry{1000.0, 2000.0, dx, dy, n, n};
}

// 14b T12: a cone, 100 m at the centre node, falling 1 m per cell.
inline terrain::raster::Raster<float> cone(std::size_t n, double dx, double dy) {
    const double mid = static_cast<double>(n - 1) / 2.0;
    std::vector<float> z(n * n);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            z[r * n + c] = static_cast<float>(
                100.0 - std::hypot(static_cast<double>(r) - mid, static_cast<double>(c) - mid));
    return terrain::raster::Raster<float>{integral_geometry(n, dx, dy), std::move(z)};
}

// 14b T12: flat sea at 0, land from 3 m on a circular coast rising inland.
inline terrain::raster::Raster<float> island(std::size_t n, double dx, double dy) {
    const double mid = static_cast<double>(n - 1) / 2.0, radius = 0.3 * static_cast<double>(n);
    std::vector<float> z(n * n, 0.0f);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c) {
            const double d = std::hypot(static_cast<double>(r) - mid, static_cast<double>(c) - mid);
            if (d <= radius) z[r * n + c] = static_cast<float>(3.0 + 0.5 * (radius - d));
        }
    return terrain::raster::Raster<float>{integral_geometry(n, dx, dy), std::move(z)};
}

// A start mesh in world coordinates with its constraint edges.
struct Start {
    terrain::IndexedMesh2 mesh;
    Edges edges;
    std::vector<std::uint32_t> masks;
};

inline terrain::Point2 world(const terrain::raster::RasterGeometry& g, double col, double row) {
    return terrain::Point2{g.x_min() + col * g.delta_x(), g.y_max() - row * g.delta_y()};
}

// A convex ring given in (col, row), counter-clockwise in world, fanned from
// vertex 0: the boundary fan increment 20 exists to remove. Every ring edge is
// a constraint with mask 1.
inline Start fan(const terrain::raster::RasterGeometry& g,
                 const std::vector<std::array<double, 2>>& ring) {
    std::vector<terrain::Point2> xy;
    for (const auto& [c, r] : ring) xy.push_back(world(g, c, r));
    const auto n = static_cast<std::uint32_t>(xy.size());
    std::vector<terrain::TriangleIndices> tris;
    for (std::uint32_t i = 1; i + 1 < n; ++i) tris.push_back({0, i, i + 1});
    Start s;
    s.mesh = terrain::IndexedMesh2{std::move(xy), std::move(tris),
                                   std::vector<std::uint8_t>(n - 2, 0)};
    for (std::uint32_t i = 0; i < n; ++i) {
        s.edges.push_back({std::min(i, (i + 1) % n), std::max(i, (i + 1) % n)});
        s.masks.push_back(1);
    }
    return s;
}

// `k` vertices on a circle about (cc, rc) in (col, row), evenly spaced from a
// phase that keeps them off-node, counter-clockwise in world (row grows down,
// so the angle runs the other way in row).
inline std::vector<std::array<double, 2>> circle_ring(double cc, double rc, double radius,
                                                      std::size_t k, double phase = 0.1234) {
    std::vector<std::array<double, 2>> ring;
    for (std::size_t i = 0; i < k; ++i) {
        const double t = phase + 2.0 * std::numbers::pi * static_cast<double>(i) / static_cast<double>(k);
        ring.push_back({cc + radius * std::cos(t), rc - radius * std::sin(t)});
    }
    return ring;
}

// The domain start Q7 pins: a 23-gon, off-node, on a 33 x 33 grid.
inline std::vector<std::array<double, 2>> q7_ring() { return circle_ring(16.21, 15.87, 15.1, 23); }

}  // namespace quality_fixtures
