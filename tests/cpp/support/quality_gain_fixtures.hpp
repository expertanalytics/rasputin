#pragma once

// Test-only scenes for increment 20c, PR 20c-2 (docs/increments/20c-soft-quality.md,
// "Tests @tester writes red first", 20c-2, T-P3): three starts refined with
// the quality start and feet on, shared by tests/cpp/property/prop_quality_gain_refine.cpp
// and by the scratch program that recorded their 20c-1 digests, so the
// recorded input and the tested input are one definition. The options here
// leave min_gain_deg alone: the recording ran on code that had no such field.
//
//   0  a features start: 20c-1's interior line (constraint_foot_fixtures::line_start)
//      on rough ground, 17 x 17 nodes, 1 m cells;
//   1  20's Q7 ring, fanned, on rough ground, 33 x 33 nodes, 10 m x 5 m cells;
//   2  test_quality_gain's line fixture (quality_gain::line_beyond: a long
//      line whose worst triangle's node lies beyond it) on rough ground,
//      33 x 25 nodes, 1 m cells.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>

#include "constraint_foot_fixtures.hpp"
#include "quality_fixtures.hpp"
#include "quality_gain_support.hpp"
#include "refinement_fixtures.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace quality_gain_fixtures {

inline constexpr std::size_t kScenes = 3;

// A quality_gain::Fixture as a start mesh in g's world.
inline quality_fixtures::Start start_of(const quality_gain::Fixture& f, const terrain::raster::RasterGeometry& g) {
    std::vector<terrain::Point2> xy;
    for (const auto& v : f.vertices) xy.push_back(quality_fixtures::world(g, v.col, v.row));
    quality_fixtures::Start s;
    s.mesh = terrain::IndexedMesh2{std::move(xy), f.triangles, std::vector<std::uint8_t>(f.triangles.size(), 0)};
    for (const auto& [e, mask] : f.constraints) {
        s.edges.push_back({e.first, e.second});
        s.masks.push_back(mask);
    }
    return s;
}

struct Scene {
    terrain::raster::Raster<float> dem;
    quality_fixtures::Start start;
    double tolerance;
};

inline Scene scene(std::size_t i) {
    namespace cff = constraint_foot_fixtures;
    namespace rf = refinement_fixtures;
    switch (i) {
    case 0:
        return Scene{terrain::raster::Raster<float>{cff::geometry(), rf::rough_dem(cff::kN, cff::kN, 7)},
                     cff::line_start(4.3), 1.0};
    case 1: {
        terrain::raster::Raster<float> dem{rf::geometry(33, 33), rf::rough_dem(33, 33, 11)};
        auto start = quality_fixtures::fan(dem.geometry(), quality_fixtures::q7_ring());
        return Scene{std::move(dem), std::move(start), 0.5};
    }
    default: {
        const terrain::raster::RasterGeometry g{0.0, 24.0, 1.0, 1.0, 33, 25};
        return Scene{terrain::raster::Raster<float>{g, rf::rough_dem(25, 33, 5)},
                     start_of(quality_gain::line_beyond(), g), 1.0};
    }
    }
}

// The quality start at 25 degrees and feet on, one thread; min_gain_deg as
// the caller leaves it.
inline terrain::refinement::RefineOptions options(const Scene& s, unsigned threads = 1) {
    terrain::refinement::RefineOptions o;
    o.tolerance = s.tolerance;
    o.threads = threads;
    o.min_angle_deg = 25.0;
    o.constraint_feet = true;
    return o;
}

inline terrain::refinement::RefineOutcome run(const Scene& s, const terrain::refinement::RefineOptions& o) {
    return terrain::refinement::refine(s.dem, s.start.mesh,
                                       std::span<const std::array<std::uint32_t, 2>>{s.start.edges},
                                       std::span<const std::uint32_t>{s.start.masks}, o);
}

}  // namespace quality_gain_fixtures
