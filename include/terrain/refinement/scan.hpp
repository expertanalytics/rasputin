#pragma once

// The per-triangle scan: how far a triangle's plane is from the DEM
// (docs/increments/14-adaptive-refinement.md, R2 and R6).
//
// A triangle's node set is every DEM node in the CLOSED triangle except its
// three vertices. Membership is three exact orientation tests, so "on an
// edge" is a zero orientation, not a tolerance: integer arithmetic when all
// three vertices are nodes, as in 14b, and DefaultKernel on (col, -row) when
// one is off-node (docs/increments/16-domain-polygon.md, R2). The box is
// ceil/floor of the fractional extents. NoData nodes are skipped. Every vertex
// is in the node rectangle, so the box is inside the grid.
//
// A vertex's z is value_at for a node and bilinear at its fractional position
// otherwise (R0), refused as raster::bilinear refuses; a vertex without one is
// a NoData vertex.
//
// For a triangle with three valid vertices the result is the largest
// |z - plane| over the set, the node where it occurs, and where that node
// lies. The plane at p is (o_bc z_a + o_ca z_b + o_ab z_c) / 2A, with integer
// orientations among nodes and double ones otherwise; only the error value
// comes out of it. The walk is row-major and only a strictly larger error
// replaces the best, so ties go to the smallest (row, col) whatever order the
// triangles are scanned in.
//
// A triangle with a NoData vertex (a void triangle, R6) has no plane. Its
// result is instead the valid node nearest any NoData vertex, by squared
// distance in the fractional frame with the same tie-break, which is where
// the refinement loop carves it; `uncovered` counts the valid nodes in its set.
//
// Pure: reads the DEM and the mesh, writes nothing but its return value, so
// any number of threads may scan one mesh at once.

#include <terrain/core/point.hpp>
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

// A vertex's height: value_at for a node; otherwise bilinear in the
// fractional frame over the cell holding it (the last cell on the border), and
// nullopt when any of that cell's four corners is NoData, whatever its weight.
template <raster::RasterSource R>
[[nodiscard]] std::optional<double> vertex_z(const R& dem, mesh::MeshVertex v) {
    const raster::RasterGeometry& g = dem.geometry();
    auto at = [&](std::size_t row, std::size_t col) -> std::optional<double> {
        if (dem.is_nodata({row, col}))
            return std::nullopt;
        return static_cast<double>(dem.value_at({row, col}));
    };
    if (v.is_node())
        return at(static_cast<std::size_t>(v.row), static_cast<std::size_t>(v.col));
    if (g.rows() < 2 || g.cols() < 2)
        return std::nullopt;
    const auto r0 = std::min(static_cast<std::size_t>(v.row), g.rows() - 2);
    const auto c0 = std::min(static_cast<std::size_t>(v.col), g.cols() - 2);
    const double ty = v.row - static_cast<double>(r0), tx = v.col - static_cast<double>(c0);
    const auto z00 = at(r0, c0), z01 = at(r0, c0 + 1), z10 = at(r0 + 1, c0),
               z11 = at(r0 + 1, c0 + 1);
    if (!z00 || !z01 || !z10 || !z11)
        return std::nullopt;
    return *z00 * (1.0 - tx) * (1.0 - ty) + *z01 * tx * (1.0 - ty) + *z10 * (1.0 - tx) * ty
         + *z11 * tx * ty;
}

template <raster::RasterSource R>
[[nodiscard]] ScanResult scan(const R& dem, const mesh::LatticeMesh& m, std::uint32_t t) {
    using mesh::LatticeVertex;
    using mesh::MeshVertex;
    const std::array<MeshVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
    const std::array<std::optional<double>, 3> zv{vertex_z(dem, v[0]), vertex_z(dem, v[1]),
                                                  vertex_z(dem, v[2])};
    const bool nodes = v[0].is_node() && v[1].is_node() && v[2].is_node();
    std::array<LatticeVertex, 3> lv{};
    if (nodes)
        lv = {v[0], v[1], v[2]};
    auto cell = [](LatticeVertex p) { return raster::CellIndex{p.row, p.col}; };

    // Twice the signed area of (v[k], v[k+1], p), with its exact sign.
    struct Orient {
        int sign;
        double value;
    };
    auto orient = [&](unsigned k, LatticeVertex p) {
        if (nodes) {
            const std::int64_t o = mesh::orient(lv[k], lv[(k + 1) % 3], p);
            return Orient{(o > 0) - (o < 0), static_cast<double>(o)};
        }
        const Point2 a = v[k].frame(), b = v[(k + 1) % 3].frame(), q = MeshVertex{p}.frame();
        return Orient{mesh::orient_sign(v[k], v[(k + 1) % 3], p),
                      (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x)};
    };

    ScanResult r;
    r.is_void = !zv[0] || !zv[1] || !zv[2];
    const Point2 f0 = v[0].frame(), f1 = v[1].frame(), f2 = v[2].frame();
    const double two_a = nodes ? static_cast<double>(mesh::orient(lv[0], lv[1], lv[2]))
                               : (f1.x - f0.x) * (f2.y - f0.y) - (f1.y - f0.y) * (f2.x - f0.x);
    auto dist2 = [&](LatticeVertex p) {
        double best = std::numeric_limits<double>::infinity();
        for (unsigned k = 0; k < 3; ++k)
            if (!zv[k]) {
                const double dr = p.row - v[k].row, dc = p.col - v[k].col;
                best = std::min(best, dr * dr + dc * dc);
            }
        return best;
    };
    double nearest = std::numeric_limits<double>::infinity();

    const auto [rlo, rhi] = std::minmax({v[0].row, v[1].row, v[2].row});
    const auto [clo, chi] = std::minmax({v[0].col, v[1].col, v[2].col});
    auto up = [](double x) { return static_cast<std::uint32_t>(std::ceil(x)); };
    auto down = [](double x) { return static_cast<std::uint32_t>(std::floor(x)); };
    for (std::uint32_t row = up(rlo); row <= down(rhi); ++row)
        for (std::uint32_t col = up(clo); col <= down(chi); ++col) {
            const LatticeVertex p{row, col};
            const Orient o0 = orient(0, p), o1 = orient(1, p), o2 = orient(2, p);
            if (o0.sign < 0 || o1.sign < 0 || o2.sign < 0 || v[0] == p || v[1] == p || v[2] == p
                || dem.is_nodata(cell(p)))
                continue;
            const NodeLocation where = o0.sign == 0   ? NodeLocation::Edge0
                                       : o1.sign == 0 ? NodeLocation::Edge1
                                       : o2.sign == 0 ? NodeLocation::Edge2
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
            const double plane =
                (o1.value * *zv[0] + o2.value * *zv[1] + o0.value * *zv[2]) / two_a;
            if (const double err = std::abs(static_cast<double>(dem.value_at(cell(p))) - plane);
                err > r.max_error) {
                r.max_error = err;
                r.node = p;
                r.where = where;
            }
        }
    return r;
}

}  // namespace terrain::refinement
