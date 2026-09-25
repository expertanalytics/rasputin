#pragma once

// The per-triangle scan: how far a triangle's plane is from the DEM
// (docs/increments/14-adaptive-refinement.md, R2 and R6).
//
// A triangle's node set is every DEM node in the CLOSED triangle except its
// three vertices. Membership is three exact integer orientation tests, so "on
// an edge" is orient == 0, not a tolerance. NoData nodes are skipped. Nodes
// outside the grid cannot occur: every vertex is a node, so the bounding box
// is inside the grid.
//
// For a triangle with three valid vertices the result is the largest
// |z - plane| over the set, the node where it occurs, and where that node
// lies. The plane at p is (o_bc z_a + o_ca z_b + o_ab z_c) / 2A with integer
// orientations. The walk is row-major and only a strictly larger error
// replaces the best, so ties go to the smallest (row, col) whatever order the
// triangles are scanned in.
//
// A triangle with a NoData vertex (a void triangle, R6) has no plane. Its
// result is instead the valid node nearest any NoData vertex, by squared
// lattice distance with the same tie-break, which is where the refinement
// loop carves it; `uncovered` counts the valid nodes in its set.
//
// Pure: reads the DEM and the mesh, writes nothing but its return value, so
// any number of threads may scan one mesh at once.

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/raster.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>

namespace terrain::refinement {

enum class NodeLocation : std::uint8_t { Inside, Edge0, Edge1, Edge2 };

struct ScanResult {
    double max_error = 0.0;                   // 0 when the set is empty, and for a void triangle
    std::optional<mesh::LatticeVertex> node;  // the argmax, or a void triangle's carve point
    NodeLocation where = NodeLocation::Inside;
    bool is_void = false;
    std::size_t uncovered = 0;  // void only: valid nodes in the set
};

template <raster::RasterSource R>
[[nodiscard]] ScanResult scan(const R& dem, const mesh::LatticeMesh& m, std::uint32_t t) {
    using mesh::LatticeVertex;
    const std::array<LatticeVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
    auto cell = [](LatticeVertex p) { return raster::CellIndex{p.row, p.col}; };
    auto z = [&](LatticeVertex p) { return static_cast<double>(dem.value_at(cell(p))); };
    std::array<bool, 3> nodata{};
    for (unsigned k = 0; k < 3; ++k)
        nodata[k] = dem.is_nodata(cell(v[k]));

    ScanResult r;
    r.is_void = nodata[0] || nodata[1] || nodata[2];
    const double two_a = static_cast<double>(mesh::orient(v[0], v[1], v[2]));
    auto dist2 = [&](LatticeVertex p) {
        std::int64_t best = std::numeric_limits<std::int64_t>::max();
        for (unsigned k = 0; k < 3; ++k)
            if (nodata[k]) {
                const auto dr = std::int64_t{p.row} - v[k].row, dc = std::int64_t{p.col} - v[k].col;
                best = std::min(best, dr * dr + dc * dc);
            }
        return best;
    };
    std::int64_t nearest = std::numeric_limits<std::int64_t>::max();

    const auto [r0, r1] = std::minmax({v[0].row, v[1].row, v[2].row});
    const auto [c0, c1] = std::minmax({v[0].col, v[1].col, v[2].col});
    for (std::uint32_t row = r0; row <= r1; ++row)
        for (std::uint32_t col = c0; col <= c1; ++col) {
            const LatticeVertex p{row, col};
            const std::int64_t o0 = mesh::orient(v[0], v[1], p);
            const std::int64_t o1 = mesh::orient(v[1], v[2], p);
            const std::int64_t o2 = mesh::orient(v[2], v[0], p);
            if (o0 < 0 || o1 < 0 || o2 < 0 || p == v[0] || p == v[1] || p == v[2]
                || dem.is_nodata(cell(p)))
                continue;
            const NodeLocation where = o0 == 0   ? NodeLocation::Edge0
                                       : o1 == 0 ? NodeLocation::Edge1
                                       : o2 == 0 ? NodeLocation::Edge2
                                                 : NodeLocation::Inside;
            if (r.is_void) {
                ++r.uncovered;
                if (const auto d = dist2(p); d < nearest) {
                    nearest = d;
                    r.node = p;
                    r.where = where;
                }
                continue;
            }
            const double plane = (static_cast<double>(o1) * z(v[0]) + static_cast<double>(o2) * z(v[1])
                                  + static_cast<double>(o0) * z(v[2]))
                               / two_a;
            if (const double err = std::abs(z(p) - plane); err > r.max_error) {
                r.max_error = err;
                r.node = p;
                r.where = where;
            }
        }
    return r;
}

}  // namespace terrain::refinement
