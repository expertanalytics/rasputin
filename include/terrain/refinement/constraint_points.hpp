#pragma once

// The edge strip's check points (docs/increments/15f-edge-strip.md, D2, D3):
// where each constraint edge crosses a grid line, and the midpoint between
// each two neighbouring crossings, the edge's ends counting as neighbours.
//
// Per edge, from its lower vertex index P0 to its higher P1, so that the
// points do not depend on the order the pair was given in: the ends go to the
// lattice by refine's own detail::lattice_position; the crossings with the
// integer column and row lines strictly between the ends are computed at the
// parameter t along P0 -> P1, on the line exactly, clamped to the node
// rectangle; a crossing whose rounded position is a node on the edge (exact
// orient_sign) is that node. Sorted by t, then position, a point at the same
// position as the entry before it, or a t not strictly between that entry's
// and 1, is dropped and counted in `duplicates`. Midpoints are taken over the
// list with the ends added, before any NoData drop; one that rounds onto a
// neighbour's position or t is a duplicate too. Every point's z is vertex_z
// there; a point for which it refuses (a NoData corner of its cell, whatever
// its weight) is dropped and counted in `no_data`. The ends are not check
// points. An edge {i, i} is refused.
//
// The store is built only by the generator and immutable afterwards, so any
// number of threads may read it.

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace terrain::refinement {

struct ConstraintPoint {
    mesh::MeshVertex at;  // lattice (col, row), as the generator computed it
    double z;             // vertex_z there
    double s;             // the parameter along its edge, P0 = 0 to P1 = 1
};
static_assert(sizeof(ConstraintPoint) == 32);

class ConstraintCheckPoints;

template <raster::RasterSource R>
[[nodiscard]] ConstraintCheckPoints constraint_check_points(
    const R& dem, std::span<const Point2> vertices,
    std::span<const std::array<std::uint32_t, 2>> edges);

class ConstraintCheckPoints {
public:
    [[nodiscard]] const raster::RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] std::size_t edge_count() const noexcept { return edges_.size(); }
    // (P0, P1), lower index first.
    [[nodiscard]] std::array<std::uint32_t, 2> edge(std::size_t k) const { return edges_.at(k); }
    // The points kept on edge(k), s strictly increasing.
    [[nodiscard]] std::span<const ConstraintPoint> on_edge(std::size_t k) const {
        return std::span{points_}.subspan(offsets_.at(k), offsets_.at(k + 1) - offsets_[k]);
    }
    [[nodiscard]] std::size_t size() const noexcept { return points_.size(); }
    [[nodiscard]] std::size_t no_data() const noexcept { return no_data_; }
    [[nodiscard]] std::size_t duplicates() const noexcept { return duplicates_; }

private:
    explicit ConstraintCheckPoints(const raster::RasterGeometry& g) : geometry_{g} {}

    template <raster::RasterSource R>
    friend ConstraintCheckPoints constraint_check_points(
        const R& dem, std::span<const Point2> vertices,
        std::span<const std::array<std::uint32_t, 2>> edges);

    raster::RasterGeometry geometry_;
    std::vector<std::array<std::uint32_t, 2>> edges_;
    std::vector<std::size_t> offsets_{0};  // edges_.size() + 1, into points_
    std::vector<ConstraintPoint> points_;
    std::size_t no_data_ = 0;
    std::size_t duplicates_ = 0;
};

template <raster::RasterSource R>
[[nodiscard]] ConstraintCheckPoints constraint_check_points(
    const R& dem, std::span<const Point2> vertices,
    std::span<const std::array<std::uint32_t, 2>> edges) {
    const raster::RasterGeometry& g = dem.geometry();
    ConstraintCheckPoints out{g};
    auto end = [&](std::uint32_t i) {
        if (i >= vertices.size())
            throw std::invalid_argument("constraint_check_points: vertex index "
                                        + std::to_string(i) + " is out of range");
        if (!g.cell_of(vertices[i]))
            throw std::invalid_argument("constraint_check_points: vertex " + std::to_string(i)
                                        + " is outside the DEM's node rectangle");
        return detail::lattice_position(g, vertices[i]);
    };
    const double cmax = static_cast<double>(g.cols() - 1);
    const double rmax = static_cast<double>(g.rows() - 1);
    struct Cut {
        double t;
        mesh::MeshVertex at;
    };
    std::vector<Cut> cuts;
    for (const auto& given : edges) {
        const std::array<std::uint32_t, 2> e{std::min(given[0], given[1]),
                                             std::max(given[0], given[1])};
        if (e[0] == e[1])
            throw std::invalid_argument("constraint_check_points: edge (" + std::to_string(e[0])
                                        + ", " + std::to_string(e[1]) + ") is degenerate");
        const mesh::MeshVertex a = end(e[0]), b = end(e[1]);
        cuts.clear();
        // Crossings: col == K exactly for a column line, row == R for a row line.
        auto cross = [&](double t, double col, double row) {
            mesh::MeshVertex at{std::clamp(col, 0.0, cmax), std::clamp(row, 0.0, rmax)};
            const mesh::MeshVertex node{std::round(at.col), std::round(at.row)};
            if (mesh::orient_sign(a, b, node) == 0)
                at = node;
            cuts.push_back({t, at});
        };
        for (double k = std::floor(std::min(a.col, b.col)) + 1; k < std::max(a.col, b.col); ++k) {
            const double t = (k - a.col) / (b.col - a.col);
            cross(t, k, a.row + t * (b.row - a.row));
        }
        for (double r = std::floor(std::min(a.row, b.row)) + 1; r < std::max(a.row, b.row); ++r) {
            const double t = (r - a.row) / (b.row - a.row);
            cross(t, a.col + t * (b.col - a.col), r);
        }
        std::sort(cuts.begin(), cuts.end(), [](const Cut& x, const Cut& y) {
            return std::tie(x.t, x.at.col, x.at.row) < std::tie(y.t, y.at.col, y.at.row);
        });
        // The list from P0 to P1: t strictly inside (t_last, 1), and no two
        // neighbours at one position (D2 step 5, I1-I3).
        std::vector<Cut> list{{0.0, a}};
        for (const Cut& c : cuts) {
            if (c.t <= list.back().t || c.t >= 1.0 || c.at == list.back().at)
                ++out.duplicates_;
            else
                list.push_back(c);
        }
        list.push_back({1.0, b});
        auto keep = [&](mesh::MeshVertex at, double s) {
            if (const auto z = vertex_z(dem, at))
                out.points_.push_back({at, *z, s});
            else
                ++out.no_data_;
        };
        for (std::size_t i = 0; i + 1 < list.size(); ++i) {
            const Cut& l = list[i];
            const Cut& r = list[i + 1];
            if (i > 0)
                keep(l.at, l.t);
            // A midpoint that rounds onto a neighbour (step 6) is a duplicate.
            const mesh::MeshVertex mid{(l.at.col + r.at.col) / 2, (l.at.row + r.at.row) / 2};
            const double s = (l.t + r.t) / 2;
            if (mid == l.at || mid == r.at || s == l.t || s == r.t)
                ++out.duplicates_;
            else
                keep(mid, s);
        }
        out.edges_.push_back(e);
        out.offsets_.push_back(out.points_.size());
    }
    return out;
}

}  // namespace terrain::refinement
