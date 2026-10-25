#pragma once

// The error summary of a mesh as written: the largest and the root-mean-square
// |z - plane| over the DEM's nodes (a raster) or the stored check points (a
// CheckPoints store) inside it, each counted once.
//
// The walk is scan's or scan_points' with an accumulator (scan.hpp), so the
// plane expression is the one refine and refine_points evaluate. A node on an
// edge shared by two triangles is owned by the lower triangle index; a vertex
// is never a point; a NoData node is not a point; a void triangle measures
// nothing, nor does a node it owns.
//
// The raster overload takes the corner heights as scan does (value_at or
// bilinear); the store overload takes z and valid per vertex, as
// refine_points does. A vertex refine_points inserted is written at
// (x_min + col dx, y_max - row dy) from its stored position, and dividing back
// need not give that position bit for bit; an off-node vertex whose stored
// position writes it back exactly is placed there, so the point it was is
// skipped as a corner and the planes are refine_points' own.
//
// Triangles are scanned in parallel into one slot each and reduced serially in
// triangle order, so the result is bit-identical for any thread count. A vertex
// outside the node rectangle or not finite, or a mesh that is not a valid
// counter-clockwise triangulation, is std::invalid_argument; an unfrozen store
// is std::logic_error.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace terrain::refinement {

struct ErrorSummary {
    double max_error = 0.0;
    double sum_squares = 0.0;
    std::size_t count = 0;

    // 0 when nothing was measured.
    [[nodiscard]] double rms() const noexcept {
        return count == 0 ? 0.0 : std::sqrt(sum_squares / static_cast<double>(count));
    }
    void add(double error, bool owned) noexcept {
        if (!owned)
            return;
        ++count;
        sum_squares += error * error;
        max_error = std::max(max_error, error);
    }
};

namespace detail {

// The lattice mesh of `mesh`, without constraints; place(lattice position,
// world point) gives each vertex's position.
template <class Place>
[[nodiscard]] mesh::LatticeMesh summary_lattice(const raster::RasterGeometry& g, const IndexedMesh2& mesh,
                                                Place place) {
    std::vector<mesh::MeshVertex> lattice;
    lattice.reserve(mesh.vertices().size());
    for (std::size_t i = 0; i < mesh.vertices().size(); ++i) {
        const Point2 p = mesh.vertices()[i];
        if (!g.cell_of(p))
            throw std::invalid_argument("error_summary: vertex " + std::to_string(i)
                                        + " is outside the node rectangle");
        lattice.push_back(place(lattice_position(g, p), p));
    }
    const std::size_t nt = mesh.triangle_count();
    auto m = mesh::LatticeMesh::build(std::move(lattice), {mesh.triangles().begin(), mesh.triangles().end()},
                                      std::vector<std::uint8_t>(nt, 0),
                                      std::vector<std::array<std::uint32_t, 3>>(nt, {0, 0, 0}));
    if (!m)
        throw std::invalid_argument("error_summary: a triangle is not counter-clockwise, or the mesh is not a "
                                    "valid triangulation");
    return std::move(*m);
}

// scan_one(t, acc) per triangle in parallel, then the serial reduction.
template <class ScanOne>
[[nodiscard]] ErrorSummary summarise(std::size_t nt, unsigned threads, ScanOne scan_one) {
    std::vector<ErrorSummary> parts(nt);
    parallel_util::for_each_block(nt, threads, parallel_util::BlockSchedule{}, [&](std::size_t b, std::size_t e) {
        for (std::size_t t = b; t < e; ++t)
            scan_one(static_cast<std::uint32_t>(t), parts[t]);
    });
    ErrorSummary total;
    for (const ErrorSummary& p : parts) {
        total.max_error = std::max(total.max_error, p.max_error);
        total.sum_squares += p.sum_squares;
        total.count += p.count;
    }
    return total;
}

}  // namespace detail

template <raster::RasterSource R>
[[nodiscard]] ErrorSummary error_summary(const R& dem, const IndexedMesh2& mesh, unsigned threads) {
    const auto m = detail::summary_lattice(dem.geometry(), mesh, [](mesh::MeshVertex v, Point2) { return v; });
    return detail::summarise(m.triangle_count(), threads, [&](std::uint32_t t, ErrorSummary& acc) {
        (void)scan(dem, m, t, 0.0, static_cast<const NoSlope*>(nullptr), nullptr, &acc);
    });
}

[[nodiscard]] inline ErrorSummary error_summary(const CheckPoints& points, const IndexedMesh2& mesh,
                                                std::span<const double> z, std::span<const std::uint8_t> valid,
                                                unsigned threads) {
    if (!points.frozen())
        throw std::logic_error("error_summary: the check-point store is not frozen");
    const raster::RasterGeometry& g = points.geometry();
    const auto m = detail::summary_lattice(g, mesh, [&](mesh::MeshVertex v, Point2 p) {
        const mesh::MeshVertex s = points.position(p);
        const Point2 back{g.x_min() + s.col * g.delta_x(), g.y_max() - s.row * g.delta_y()};
        return !v.is_node() && back == p ? s : v;
    });
    std::vector<double> zt(m.vertices().size(), std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < zt.size() && i < z.size() && i < valid.size(); ++i)
        if (valid[i] != 0)
            zt[i] = z[i];
    return detail::summarise(m.triangle_count(), threads, [&](std::uint32_t t, ErrorSummary& acc) {
        (void)detail::scan_points(points, m, zt, t, 0.0, static_cast<const NoSlope*>(nullptr), nullptr, &acc);
    });
}

}  // namespace terrain::refinement
