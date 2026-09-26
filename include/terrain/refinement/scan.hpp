#pragma once

// The per-triangle scan: how far a triangle's plane is from the DEM
// (docs/increments/14-adaptive-refinement.md, R2 and R6).
//
// A triangle's node set is every DEM node in the CLOSED triangle except its
// three vertices. Membership is three exact orientation tests, so "on an
// edge" is a zero orientation, not a tolerance: integer arithmetic when all
// three vertices are nodes, as in 14b, and DefaultKernel on (col, -row) when
// one is off-node (docs/increments/16-domain-polygon.md, R2). The set is not
// tested node by node: mesh::for_each_row_span yields it exactly as one column
// interval per row, and each interval is walked as contiguous row segments
// (docs/increments/18-row-span-scan.md, R1 and R6). NoData nodes are skipped.
// Every vertex is in the node rectangle, so every span is inside the grid.
//
// A vertex's z is value_at for a node and bilinear at its fractional position
// otherwise (R0), refused as raster::bilinear refuses; a vertex without one is
// a NoData vertex.
//
// For a triangle with three valid vertices the result is the largest
// |z - plane| over the set, the node where it occurs, and where that node
// lies. The plane at p is (o_bc z_a + o_ca z_b + o_ab z_c) / 2A, with integer
// orientations among nodes (advanced exactly along a row) and today's double
// expression per node otherwise, so the error is bit-identical to increment
// 17's box walk (18's C1 (a)); only the error value comes out of it. The walk is row-major and only a strictly larger error
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
#include <terrain/mesh/row_spans.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/row_segments.hpp>

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
    using T = typename R::value_type;
    const std::array<MeshVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
    const std::array<std::optional<double>, 3> zv{vertex_z(dem, v[0]), vertex_z(dem, v[1]),
                                                  vertex_z(dem, v[2])};
    const bool nodes = v[0].is_node() && v[1].is_node() && v[2].is_node();
    std::array<LatticeVertex, 3> lv{};
    if (nodes)
        lv = {*v[0].as_node(), *v[1].as_node(), *v[2].as_node()};

    // The exact sign of (v[k], v[k+1], p), and for an off-node triangle the
    // double value today's plane uses. Among nodes the value is the integer
    // orientation, which advances by exactly step[k] per column.
    auto sign = [&](unsigned k, LatticeVertex p) {
        if (nodes) {
            const std::int64_t o = mesh::orient(lv[k], lv[(k + 1) % 3], p);
            return (o > 0) - (o < 0);
        }
        return mesh::orient_sign(v[k], v[(k + 1) % 3], p);
    };
    auto value = [&](unsigned k, LatticeVertex p) {
        const Point2 a = v[k].frame(), b = v[(k + 1) % 3].frame(), q = MeshVertex{p}.frame();
        return (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x);
    };
    const std::array<std::int64_t, 3> step{std::int64_t{lv[1].row} - lv[0].row,
                                           std::int64_t{lv[2].row} - lv[1].row,
                                           std::int64_t{lv[0].row} - lv[2].row};

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
    const std::optional<T>& nd = dem.nodata();
    auto missing = [&](T z) { return z != z || (nd && z == *nd); };

    // In an extreme off-node sliver two_a can round to 0 (or below) while its
    // exact sign is positive, and the plane would be NaN. p is in the closed
    // triangle, so the plane is a convex combination of the corner heights,
    // and the largest corner difference bounds its error from above.
    auto error = [&](double z, double o0, double o1, double o2) {
        return two_a > 0.0
                   ? std::abs(z - (o1 * *zv[0] + o2 * *zv[1] + o0 * *zv[2]) / two_a)
                   : std::max({std::abs(z - *zv[0]), std::abs(z - *zv[1]), std::abs(z - *zv[2])});
    };
    // NoData needs no branch: a NaN error never beats the best, and a
    // sentinel's error is selected to 0, which never does either.
    auto consider = [&](T z, double err, LatticeVertex p) {
        err = missing(z) ? 0.0 : err;
        if (err > r.max_error) {
            r.max_error = err;
            r.node = p;
        }
    };

    mesh::for_each_row_span(v, [&](mesh::RowSpan span) {
        raster::for_each_row_segment(dem, span.row, span.c0, span.c1, [&](raster::RowSegment<T> s) {
            const std::uint32_t row = span.row, c = s.first_col;
            if (r.is_void) {
                for (std::uint32_t j = 0; j < s.values.size(); ++j) {
                    if (missing(s.values[j]))
                        continue;
                    ++r.uncovered;
                    const LatticeVertex p{row, c + j};
                    if (const auto d = dist2(p); d < nearest) {
                        nearest = d;
                        r.node = p;
                    }
                }
            } else if (nodes) {
                std::array<std::int64_t, 3> o{};
                for (unsigned k = 0; k < 3; ++k)
                    o[k] = mesh::orient(lv[k], lv[(k + 1) % 3], LatticeVertex{row, c});
                for (std::uint32_t j = 0; j < s.values.size(); ++j) {
                    const T z = s.values[j];
                    consider(z,
                             error(static_cast<double>(z), static_cast<double>(o[0]),
                                   static_cast<double>(o[1]), static_cast<double>(o[2])),
                             LatticeVertex{row, c + j});
                    for (unsigned k = 0; k < 3; ++k)
                        o[k] += step[k];
                }
            } else {
                for (std::uint32_t j = 0; j < s.values.size(); ++j) {
                    const LatticeVertex p{row, c + j};
                    const T z = s.values[j];
                    consider(z, error(static_cast<double>(z), value(0, p), value(1, p), value(2, p)), p);
                }
            }
        });
    });
    // Where the recorded node lies, by today's exact zero-tests in today's
    // order; once per triangle rather than per candidate.
    if (r.node)
        r.where = sign(0, *r.node) == 0   ? NodeLocation::Edge0
                  : sign(1, *r.node) == 0 ? NodeLocation::Edge1
                  : sign(2, *r.node) == 0 ? NodeLocation::Edge2
                                          : NodeLocation::Inside;
    return r;
}

}  // namespace terrain::refinement
