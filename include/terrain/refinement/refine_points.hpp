#pragma once

// The final check (docs/increments/15c-geographic-dem.md, D5): greedy insertion
// against stored check points, starting from phase 1's mesh. refine.hpp's loop
// with a different scanner; refine.hpp itself is not edited (J1).
//
// The store's geometry() is the frame and the buckets, nothing more: phase 2
// reads no raster value. A start vertex keeps the z it is given (NaN where it
// is not valid); an inserted vertex is a check point at its stored position,
// output at (x_min + col dx, y_max - row dy), with the point's own z.
//
// Scan (parallel, read-only, one result per triangle). For each cell row the
// triangle meets, the column range it covers in that band, widened by a cell
// each side, bounds the stored points tested; membership in the CLOSED
// triangle is three exact orient_sign calls. A point equal to a corner is
// skipped. Error is |z - plane| with scan.hpp's off-node double expression
// (the largest corner difference for a sliver whose 2A rounds to <= 0); the
// worst point wins by strictly larger error, so ties go to the first in store
// order. A void triangle (a NaN corner) takes the point nearest a void corner,
// and counts its points in `uncovered`.
//
// Split (serial, triangle-index order), as refine: split_inside, or split_edge
// for a point on an edge, skipped this round when the neighbour across it was
// already touched; legalise_around follows. The start is legalised once first,
// as refine does. Terminates: each insertion is a stored point not yet a
// vertex, and the store is finite.
//
// Coincident points. A stored point equal to a start vertex is a corner of
// every triangle holding it, so it is never scanned nor inserted. After the
// loop each start vertex looks it up in its own cell; matches are counted in
// `coincident`, their largest |z - vertex z| in `coincident_max_error`.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace terrain::refinement {

struct PointRefineOptions {
    double tolerance = 0.0;  // metres, finite and >= 0
    unsigned threads = 0;    // 0: hardware concurrency
};

// RefineOutcome, with max_error over check points, plus the coincident points.
struct PointRefineOutcome : RefineOutcome {
    std::size_t coincident = 0;         // stored points equal to a start vertex
    double coincident_max_error = 0.0;  // their largest |z - vertex z|, valid vertices only
};

namespace detail {

struct PointScan {
    double max_error = 0.0;                  // 0 when no point is in the triangle, and when void
    std::optional<mesh::MeshVertex> point;  // the argmax, or a void triangle's carve point
    float z = 0.0f;                          // its z
    NodeLocation where = NodeLocation::Inside;
    bool is_void = false;
    std::size_t uncovered = 0;  // void only: points in the triangle
};

template <class Store>
[[nodiscard]] PointScan scan_points(const Store& points, const mesh::LatticeMesh& m,
                                    std::span<const double> zt, std::uint32_t t) {
    using mesh::MeshVertex;
    const raster::RasterGeometry& g = points.geometry();
    const auto& tri = m.triangles()[t];
    const std::array<MeshVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
    const std::array<double, 3> zv{zt[tri[0]], zt[tri[1]], zt[tri[2]]};
    PointScan r;
    r.is_void = std::isnan(zv[0]) || std::isnan(zv[1]) || std::isnan(zv[2]);
    const Point2 f0 = v[0].frame(), f1 = v[1].frame(), f2 = v[2].frame();
    const double two_a = (f1.x - f0.x) * (f2.y - f0.y) - (f1.y - f0.y) * (f2.x - f0.x);
    auto value = [&](unsigned k, MeshVertex p) {
        const Point2 a = v[k].frame(), b = v[(k + 1) % 3].frame(), q = p.frame();
        return (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x);
    };
    double nearest = std::numeric_limits<double>::infinity();
    auto visit = [&](MeshVertex p, float zf) {
        if (p == v[0] || p == v[1] || p == v[2])
            return;
        for (unsigned k = 0; k < 3; ++k)
            if (mesh::orient_sign(v[k], v[(k + 1) % 3], p) < 0)
                return;
        const auto z = static_cast<double>(zf);
        if (r.is_void) {
            ++r.uncovered;
            double d = std::numeric_limits<double>::infinity();
            for (unsigned k = 0; k < 3; ++k)
                if (std::isnan(zv[k]))
                    d = std::min(d, (p.row - v[k].row) * (p.row - v[k].row) + (p.col - v[k].col) * (p.col - v[k].col));
            if (d < nearest) {
                nearest = d;
                r.point = p;
                r.z = zf;
            }
            return;
        }
        const double err =
            two_a > 0.0 ? std::abs(z - (value(1, p) * zv[0] + value(2, p) * zv[1] + value(0, p) * zv[2]) / two_a)
                        : std::max({std::abs(z - zv[0]), std::abs(z - zv[1]), std::abs(z - zv[2])});
        if (err > r.max_error) {
            r.max_error = err;
            r.point = p;
            r.z = zf;
        }
    };

    // Per cell row, the triangle's column range in the band [b, b + 1]: its
    // vertices in the band and its edges' crossings of the band's two lines.
    const std::size_t last_row = last_cell(g.rows()), last_col = last_cell(g.cols());
    const auto [rlo, rhi] = std::minmax({v[0].row, v[1].row, v[2].row});
    const auto b_hi = std::min(static_cast<std::size_t>(rhi), last_row);
    for (std::size_t b = std::min(static_cast<std::size_t>(rlo), last_row); b <= b_hi; ++b) {
        const auto top = static_cast<double>(b), bottom = top + 1.0;
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (unsigned k = 0; k < 3; ++k) {
            const MeshVertex a = v[k], c = v[(k + 1) % 3];
            if (a.row >= top && a.row <= bottom) {
                lo = std::min(lo, a.col);
                hi = std::max(hi, a.col);
            }
            for (const double y : {top, bottom})
                if ((a.row - y) * (c.row - y) < 0.0) {
                    const double x = a.col + (c.col - a.col) * (y - a.row) / (c.row - a.row);
                    lo = std::min(lo, x);
                    hi = std::max(hi, x);
                }
        }
        if (!(lo <= hi))
            continue;
        const auto c0 = static_cast<std::size_t>(std::max(std::floor(lo) - 1.0, 0.0));
        const auto c1 = std::min(static_cast<std::size_t>(std::floor(hi) + 1.0), last_col);
        points.for_each_in(b, c0, c1, visit);
    }
    if (r.point)
        r.where = mesh::orient_sign(v[0], v[1], *r.point) == 0   ? NodeLocation::Edge0
                  : mesh::orient_sign(v[1], v[2], *r.point) == 0 ? NodeLocation::Edge1
                  : mesh::orient_sign(v[2], v[0], *r.point) == 0 ? NodeLocation::Edge2
                                                                 : NodeLocation::Inside;
    return r;
}

}  // namespace detail

// Store: CheckPoints, or a test double with geometry(), frozen() and for_each_in.
// An unfrozen store is a programming error (std::logic_error), as add after
// freeze is: RefineStatus has no value for it, and refine.hpp is not edited.
template <class Store>
[[nodiscard]] PointRefineOutcome refine_points(const Store& points, const IndexedMesh2& start,
                                               std::span<const double> z, std::span<const std::uint8_t> valid,
                                               std::span<const std::array<std::uint32_t, 2>> edges,
                                               std::span<const std::uint32_t> masks,
                                               const PointRefineOptions& options) {
    PointRefineOutcome out;
    const auto refuse = [&](RefineOutcome r) {
        static_cast<RefineOutcome&>(out) = std::move(r);
        return out;
    };
    if (!std::isfinite(options.tolerance) || options.tolerance < 0.0)
        return refuse(detail::refusal(RefineStatus::InvalidTolerance,
                                      "refine_points: tolerance must be finite and >= 0"));
    if (!points.frozen())
        throw std::logic_error("refine_points: the check-point store is not frozen");
    const raster::RasterGeometry& g = points.geometry();
    auto built = detail::to_lattice(g, start, edges, masks);
    if (auto* refused = std::get_if<RefineOutcome>(&built))
        return refuse(std::move(*refused));
    auto& m = std::get<mesh::LatticeMesh>(built);

    const std::size_t n0 = start.vertices().size();
    std::vector<double> zt(n0, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < n0 && i < z.size() && i < valid.size(); ++i)
        if (valid[i] != 0)
            zt[i] = z[i];

    const auto frame = mesh::lattice_frame(g.delta_x(), g.delta_y(), g.rows(), g.cols());
    out.flips = mesh::legalise_all<pred::DefaultKernel>(m, frame, [](std::uint32_t) {});
    std::vector<detail::PointScan> results;
    mesh::FlipStack flip_stack;
    std::vector<std::uint32_t> active(m.triangle_count());
    for (std::uint32_t t = 0; t < active.size(); ++t)
        active[t] = t;

    using clock = std::chrono::steady_clock;
    const auto since = [](clock::time_point t0) {
        return std::chrono::duration<double>(clock::now() - t0).count();
    };
    while (true) {
        ++out.rounds;
        results.resize(m.triangle_count());
        auto t0 = clock::now();
        parallel_util::for_each_block(active.size(), options.threads, parallel_util::BlockSchedule{},
                                      [&](std::size_t begin, std::size_t end) {
                                          for (std::size_t i = begin; i < end; ++i)
                                              results[active[i]] = detail::scan_points(points, m, zt, active[i]);
                                      });
        out.scan_seconds += since(t0);
        t0 = clock::now();
        std::vector<char> touched(m.triangle_count(), 0);
        std::vector<std::uint32_t> skipped;
        bool any = false;
        for (const std::uint32_t t : active) {
            const detail::PointScan& r = results[t];
            if (!r.point || !(r.is_void || r.max_error > options.tolerance))
                continue;
            any = true;
            if (touched[t] != 0)
                continue;
            const auto before = static_cast<std::uint32_t>(m.triangle_count());
            std::array<std::uint32_t, 4> seeds{t, before, before + 1, before + 1};
            std::size_t n_seeds = 3;
            std::uint32_t q = 0;
            if (r.where == NodeLocation::Inside) {
                q = m.split_inside(t, *r.point);
            } else {
                const auto e = static_cast<unsigned>(r.where) - 1;
                const std::uint32_t u = m.neighbours(t)[e];
                if (u != mesh::kNoNeighbour && touched[u] != 0) {
                    skipped.push_back(t);
                    continue;
                }
                q = m.split_edge(t, e, *r.point);
                n_seeds = u != mesh::kNoNeighbour ? 4 : 2;
                if (n_seeds == 4)
                    seeds[3] = u;
            }
            zt.push_back(static_cast<double>(r.z));
            touched.resize(m.triangle_count(), 1);
            touched[t] = 1;
            if (n_seeds == 4)
                touched[seeds[3]] = 1;
            out.flips += mesh::legalise_around<pred::DefaultKernel>(
                m, q, std::span<const std::uint32_t>{seeds.data(), n_seeds}, frame, flip_stack,
                [&](std::uint32_t s) { touched[s] = 1; });
            ++out.inserted;
            out.carved += r.is_void ? 1 : 0;
        }
        out.split_seconds += since(t0);
        if (!any)
            break;
        detail::rebuild_active(touched, skipped, active);
    }

    for (const detail::PointScan& r : results) {
        if (r.is_void)
            out.uncovered += r.uncovered;
        else
            out.max_error = std::max(out.max_error, r.max_error);
    }
    for (std::size_t i = 0; i < m.vertices().size(); ++i) {
        const mesh::MeshVertex v = m.vertices()[i];
        if (i < n0) {
            const auto row = std::min(static_cast<std::size_t>(v.row), last_cell(g.rows()));
            const auto col = std::min(static_cast<std::size_t>(v.col), last_cell(g.cols()));
            points.for_each_in(row, col, col, [&](mesh::MeshVertex p, float pz) {
                if (p != v)
                    return;
                ++out.coincident;
                if (!std::isnan(zt[i]))
                    out.coincident_max_error = std::max(out.coincident_max_error, std::abs(static_cast<double>(pz) - zt[i]));
            });
        }
        out.vertices.push_back(i < n0 ? start.vertices()[i]
                                      : Point2{g.x_min() + v.col * g.delta_x(), g.y_max() - v.row * g.delta_y()});
        out.z.push_back(std::isnan(zt[i]) ? 0.0 : zt[i]);
        out.valid.push_back(std::isnan(zt[i]) ? 0 : 1);
    }
    out.triangles.assign(m.triangles().begin(), m.triangles().end());
    std::tie(out.edges, out.masks) = m.constraint_edges();
    return out;
}

}  // namespace terrain::refinement
