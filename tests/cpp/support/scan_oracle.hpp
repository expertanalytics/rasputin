#pragma once

// THE FROZEN ORACLE for increment 18 (docs/increments/18-row-span-scan.md, R5).
//
// scan_bbox and vertex_z_bbox are include/terrain/refinement/scan.hpp's scan
// and vertex_z COPIED VERBATIM at e0578e5 (increment 17's scan, the design
// commit), before any production change on the branch. The only edits are the
// two names, the namespace, and qualifying the result types, which stay the
// production ones so a result compares field for field. It walks every node of
// the triangle's bounding box and tests three exact orientations per node.
//
// FROZEN: no test may change this file to agree with new code. T1 and T2
// compare the row-span scan against it. bbox_node_set is the same box walk
// with the DEM taken out: the geometric node set, sorted row-major, that T1
// compares for_each_row_span's union of spans against.
//
// Test-only. Nothing under include/ may include it.

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace scan_oracle {

using terrain::Point2;
namespace mesh = terrain::mesh;
namespace raster = terrain::raster;
using terrain::refinement::NodeLocation;
using terrain::refinement::ScanResult;

// A vertex's height: value_at for a node; otherwise bilinear in the
// fractional frame over the cell holding it (the last cell on the border), and
// nullopt when any of that cell's four corners is NoData, whatever its weight.
template <raster::RasterSource R>
[[nodiscard]] std::optional<double> vertex_z_bbox(const R& dem, mesh::MeshVertex v) {
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
[[nodiscard]] ScanResult scan_bbox(const R& dem, const mesh::LatticeMesh& m, std::uint32_t t) {
    using mesh::LatticeVertex;
    using mesh::MeshVertex;
    const std::array<MeshVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
    const std::array<std::optional<double>, 3> zv{vertex_z_bbox(dem, v[0]), vertex_z_bbox(dem, v[1]),
                                                  vertex_z_bbox(dem, v[2])};
    const bool nodes = v[0].is_node() && v[1].is_node() && v[2].is_node();
    std::array<LatticeVertex, 3> lv{};
    if (nodes)
        lv = {*v[0].as_node(), *v[1].as_node(), *v[2].as_node()};
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
            // In an extreme off-node sliver two_a can round to 0 (or below)
            // while its exact sign is positive, and the plane would be NaN,
            // which `err > max_error` skips. p is in the closed triangle, so
            // the plane is a convex combination of the corner heights, and the
            // largest corner difference bounds its error from above.
            const auto z = static_cast<double>(dem.value_at(cell(p)));
            const double err =
                two_a > 0.0
                    ? std::abs(z - (o1.value * *zv[0] + o2.value * *zv[1] + o0.value * *zv[2])
                                       / two_a)
                    : std::max({std::abs(z - *zv[0]), std::abs(z - *zv[1]), std::abs(z - *zv[2])});
            if (err > r.max_error) {
                r.max_error = err;
                r.node = p;
                r.where = where;
            }
        }
    return r;
}

// The box walk's membership test alone: every node p of the bounding box with
// all three orientations >= 0 (exact: integer among nodes, DefaultKernel
// otherwise), except a node equal to a vertex. Row-major, so already sorted.
[[nodiscard]] inline std::vector<mesh::LatticeVertex> bbox_node_set(
    const std::array<mesh::MeshVertex, 3>& v) {
    using mesh::LatticeVertex;
    const bool nodes = v[0].is_node() && v[1].is_node() && v[2].is_node();
    auto sign = [&](unsigned k, LatticeVertex p) {
        if (nodes) {
            const std::int64_t o = mesh::orient(*v[k].as_node(), *v[(k + 1) % 3].as_node(), p);
            return (o > 0) - (o < 0);
        }
        return mesh::orient_sign(v[k], v[(k + 1) % 3], p);
    };
    std::vector<LatticeVertex> out;
    const auto [rlo, rhi] = std::minmax({v[0].row, v[1].row, v[2].row});
    const auto [clo, chi] = std::minmax({v[0].col, v[1].col, v[2].col});
    auto up = [](double x) { return static_cast<std::uint32_t>(std::ceil(x)); };
    auto down = [](double x) { return static_cast<std::uint32_t>(std::floor(x)); };
    for (std::uint32_t row = up(rlo); row <= down(rhi); ++row)
        for (std::uint32_t col = up(clo); col <= down(chi); ++col) {
            const LatticeVertex p{row, col};
            if (sign(0, p) < 0 || sign(1, p) < 0 || sign(2, p) < 0 || v[0] == p || v[1] == p
                || v[2] == p)
                continue;
            out.push_back(p);
        }
    return out;
}

}  // namespace scan_oracle
