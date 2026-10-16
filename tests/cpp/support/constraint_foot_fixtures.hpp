#pragma once

// Test-only fixtures for increment 20c, PR 20c-1 (docs/increments/20c-soft-quality.md,
// "Design of PR 20c-1" and "Tests @tester writes red first", 20c-1): the foot
// rule on every insertion path. Shared by the four test_constraint_foot*
// suites and by the scratch program that recorded the off-switch digests from
// c074f900 (master's production code; the branch changed only docs), so the
// recorded input and the tested input are one definition.
//
// The frame is a 17 x 17 node grid with x_min 0, y_max 16 and dx = dy = 1, so
// world (x, y) is (col, 16 - row) and a dyadic (col, row) maps both ways
// exactly. delta = min(dx, dy) / 2 = 0.5 cells on every path (R2.2, R4.4).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/refinement/check_points.hpp>

#include "quality_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace constraint_foot_fixtures {

inline constexpr std::size_t kN = 17;

inline terrain::raster::RasterGeometry geometry() {
    return terrain::raster::RasterGeometry{0.0, 16.0, 1.0, 1.0, kN, kN};
}

inline terrain::Point2 world(double col, double row) { return terrain::Point2{col, 16.0 - row}; }

// A features start: the square [0, 16] x [0, 16] (cols, rows) with one
// interior constraint line A-B along row `row_ab` from col 0.5 to col 15.5,
// in six triangles. Vertices 0..3 are the corners (0,0), (16,0), (16,16),
// (0,16) in (col, row); 4 is A and 5 is B. The square's sides carry mask 1,
// A-B mask 2.
inline quality_fixtures::Start line_start(double row_ab = 4.25) {
    std::vector<terrain::Point2> xy{world(0, 0),  world(16, 0),      world(16, 16),
                                    world(0, 16), world(0.5, row_ab), world(15.5, row_ab)};
    std::vector<terrain::TriangleIndices> tris{{0, 4, 5}, {0, 5, 1}, {1, 5, 2},
                                               {5, 4, 2}, {4, 3, 2}, {0, 3, 4}};
    quality_fixtures::Start s;
    s.mesh = terrain::IndexedMesh2{std::move(xy), std::move(tris), std::vector<std::uint8_t>(6, 0)};
    s.edges = {{0, 1}, {1, 2}, {2, 3}, {0, 3}, {4, 5}};
    s.masks = {1, 1, 1, 1, 2};
    return s;
}

// The start's z and valid flags, bilinear in `dem` at each start vertex.
inline std::pair<std::vector<double>, std::vector<std::uint8_t>> start_z(
    const terrain::raster::Raster<float>& dem, const terrain::IndexedMesh2& mesh) {
    std::pair<std::vector<double>, std::vector<std::uint8_t>> out;
    for (const auto& p : mesh.vertices()) {
        const auto z = terrain::raster::bilinear(dem, p);
        out.first.push_back(z.value_or(0.0));
        out.second.push_back(z ? 1 : 0);
    }
    return out;
}

// A check-point store with exact double positions, for points that a float
// cell offset cannot hold (within r(g) of a line, or a fixed distance from a
// tilted one). refine_points' Store: geometry(), frozen() and for_each_in,
// the cell rule CheckPoints' (floor, the last cell for the far edge), the
// points of a cell in the order given.
class ExactStore {
public:
    explicit ExactStore(terrain::raster::RasterGeometry g) : g_{g} {}

    void add(terrain::mesh::MeshVertex p, float z) { points_.push_back({p, z}); }

    [[nodiscard]] const terrain::raster::RasterGeometry& geometry() const noexcept { return g_; }
    [[nodiscard]] bool frozen() const noexcept { return true; }
    [[nodiscard]] std::size_t size() const noexcept { return points_.size(); }
    [[nodiscard]] const std::vector<std::pair<terrain::mesh::MeshVertex, float>>& points() const noexcept {
        return points_;
    }

    template <class F>
    void for_each_in(std::size_t row, std::size_t c0, std::size_t c1, F&& f) const {
        using terrain::refinement::last_cell;
        for (const auto& [p, z] : points_) {
            const auto r = std::min(static_cast<std::size_t>(std::max(p.row, 0.0)), last_cell(g_.rows()));
            const auto c = std::min(static_cast<std::size_t>(std::max(p.col, 0.0)), last_cell(g_.cols()));
            if (r == row && c >= c0 && c <= c1)
                f(p, z);
        }
    }

private:
    terrain::raster::RasterGeometry g_;
    std::vector<std::pair<terrain::mesh::MeshVertex, float>> points_;
};

// Every node of `dem` with data, as stored points with the node's own z: the
// final check against a source DEM on the target grid itself.
inline ExactStore every_node(const terrain::raster::Raster<float>& dem) {
    ExactStore s{dem.geometry()};
    for (std::size_t r = 0; r < dem.geometry().rows(); ++r)
        for (std::size_t c = 0; c < dem.geometry().cols(); ++c)
            if (!dem.is_nodata({r, c}))
                s.add(terrain::mesh::MeshVertex{static_cast<double>(c), static_cast<double>(r)},
                      dem.value_at({r, c}));
    return s;
}

}  // namespace constraint_foot_fixtures
