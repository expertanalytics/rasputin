#pragma once

// Output oracles for increment 20c, PR 20c-1 (docs/increments/20c-soft-quality.md,
// "Tests @tester writes red first", 20c-1), shared by the test_constraint_foot*
// suites. Each RETURNS its findings and reads only what refine, refine_points
// or refine_strip returned (world vertices, z, valid, triangles, constraint
// edges, masks) and the start it was given: no scan record, no counter. No
// Catch2 include.
//
// The QA rules' section D, on every path:
//   - delaunay_violations: no interior edge that is not a constraint edge has
//     an apex strictly inside the other triangle's circle, by the exact
//     incircle in the producer's frame (col dx, -(row dy)), on the fractional
//     (col, row) recovered from the world point, and (since 20c-1's
//     green-step ruling 2) inside by more than a margin, see below;
//   - node_findings: every valid DEM node in every closed output triangle with
//     three valid vertices is within tolerance of that triangle's plane,
//     recomputed here (orientation and membership exact on (col, -row)).
// And Q3's "the constraint edges as a set of lines are unchanged, bits and
// masks on both halves" (line_findings): every output constraint edge lies on
// exactly one input segment, to kOnLine cells, with its mask, and the pieces
// of each input segment chain from one end to the other.
//
// Scale: kOnLine = 1e-9 cells, for grids up to a few hundred cells a side; a
// foot is a double projection, so it is on its line to rounding, not exactly.

#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace constraint_foot_oracles {

using terrain::Point2;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;

inline constexpr double kOnLine = 1e-9;  // cells

struct Frac {
    double col;
    double row;
};

inline Frac frac(const terrain::raster::RasterGeometry& g, Point2 p) {
    return Frac{(p.x - g.x_min()) / g.delta_x(), (g.y_max() - p.y) / g.delta_y()};
}

// Parameter along a-b and perpendicular distance, in cells.
inline std::pair<double, double> param_dist(Frac a, Frac b, Frac p) {
    const double ux = b.col - a.col, uy = b.row - a.row, len2 = ux * ux + uy * uy;
    return {((p.col - a.col) * ux + (p.row - a.row) * uy) / len2,
            std::abs((p.col - a.col) * uy - (p.row - a.row) * ux) / std::sqrt(len2)};
}

// How far `d` lies inside the circle through a, b, c: R - |d - centre|, in the
// units of the points, centre and radius in double (relative to a).
inline double incircle_depth(Point2 a, Point2 b, Point2 c, Point2 d) {
    const double bx = b.x - a.x, by = b.y - a.y, cx = c.x - a.x, cy = c.y - a.y;
    const double den = 2.0 * (bx * cy - by * cx);
    const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    const double ux = (cy * b2 - by * c2) / den, uy = (bx * c2 - cx * b2) / den;
    return std::hypot(ux, uy) - std::hypot(d.x - a.x - ux, d.y - a.y - uy);
}

// The margin (20c-1's green-step ruling 2): an edge is a violation when the
// exact incircle says Inside AND the apex is inside the circle by more than
// eta = 1e-7 * min(dx, dy) in the producer frame. Scale: 5e-7 m on the 10 m x
// 5 m fixtures, about 500 times the rounding of the (col, row) rebuild at
// y_max = 7e6 m, where it was checked; it holds while half an ulp of the
// largest world coordinate stays under eta / 100 (world coordinates below
// about 4e7 m at 5 m cells). The producer makes near-cocircular quads along a
// footed line (exact incircle +1.77e-8 on terms around 60), which rounding in
// the rebuild or in the producer's contracted arithmetic can tip to Inside; a
// real missed flip is off by a fraction of a cell. Every quad the margin
// excuses has its depth appended to `excused`, if given, for the case to print.
template <typename Outcome>
std::size_t delaunay_violations(const terrain::raster::RasterGeometry& g, const Outcome& out,
                                std::vector<double>* excused = nullptr) {
    const double eta = 1e-7 * std::min(g.delta_x(), g.delta_y());
    auto lf = [&](std::uint32_t i) {
        const Frac f = frac(g, out.vertices[i]);
        return Point2{f.col * g.delta_x(), -(f.row * g.delta_y())};
    };
    std::set<std::pair<std::uint32_t, std::uint32_t>> constrained;
    for (const auto& e : out.edges) constrained.insert(std::minmax(e[0], e[1]));
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> sides;
    for (std::size_t t = 0; t < out.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            sides[std::minmax(out.triangles[t][k], out.triangles[t][(k + 1) % 3])].push_back(t);
    std::size_t bad = 0;
    for (const auto& [e, ts] : sides) {
        if (ts.size() > 2) ++bad;  // an edge in three triangles is no triangulation
        if (ts.size() != 2 || constrained.count(e) != 0) continue;
        for (unsigned s = 0; s < 2; ++s) {
            const auto& tri = out.triangles[ts[s]];
            std::uint32_t apex = 0;
            for (const auto x : out.triangles[ts[1 - s]])
                if (x != e.first && x != e.second) apex = x;
            const Point2 a = lf(tri[0]), b = lf(tri[1]), c = lf(tri[2]), d = lf(apex);
            if (DefaultKernel::incircle(a, b, c, d) != terrain::pred::Incircle::Inside) continue;
            const double depth = incircle_depth(a, b, c, d);
            if (!(depth <= eta))  // a NaN depth (a flat triangle) is not excused
                ++bad;
            else if (excused != nullptr)
                excused->push_back(depth);
        }
    }
    return bad;
}

struct NodeFindings {
    std::size_t over = 0;     // (node, triangle) pairs over tolerance
    std::size_t not_ccw = 0;  // triangles not counter-clockwise in (col, -row)
    double worst = 0.0;
};

template <typename Outcome>
NodeFindings node_findings(const terrain::raster::Raster<float>& dem, const Outcome& out, double tol) {
    const auto& g = dem.geometry();
    std::vector<Point2> fp;
    double zmax = 1.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const Frac f = frac(g, out.vertices[i]);
        fp.push_back(Point2{f.col, -f.row});
        zmax = std::max(zmax, std::abs(out.z[i]));
    }
    auto cross = [](Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); };
    NodeFindings nf;
    for (const auto& tri : out.triangles) {
        const Point2 a = fp[tri[0]], b = fp[tri[1]], c = fp[tri[2]];
        if (DefaultKernel::orient2d(a, b, c) != Orientation::CounterClockwise) ++nf.not_ccw;
        if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
        const double two_a = cross(a, b, c);
        const auto lo_c = static_cast<std::int64_t>(std::ceil(std::min({a.x, b.x, c.x})));
        const auto hi_c = static_cast<std::int64_t>(std::floor(std::max({a.x, b.x, c.x})));
        const auto lo_r = static_cast<std::int64_t>(std::ceil(-std::max({a.y, b.y, c.y})));
        const auto hi_r = static_cast<std::int64_t>(std::floor(-std::min({a.y, b.y, c.y})));
        for (std::int64_t r = std::max<std::int64_t>(lo_r, 0); r <= hi_r && r < static_cast<std::int64_t>(g.rows()); ++r)
            for (std::int64_t col = std::max<std::int64_t>(lo_c, 0); col <= hi_c && col < static_cast<std::int64_t>(g.cols()); ++col) {
                const terrain::raster::CellIndex n{static_cast<std::size_t>(r), static_cast<std::size_t>(col)};
                const Point2 p{static_cast<double>(col), -static_cast<double>(r)};
                if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                    continue;
                if (p == a || p == b || p == c || dem.is_nodata(n)) continue;
                const double plane = (cross(p, b, c) * out.z[tri[0]] + cross(a, p, c) * out.z[tri[1]]
                                      + cross(a, b, p) * out.z[tri[2]]) / two_a;
                const double err = std::abs(static_cast<double>(dem.value_at(n)) - plane);
                nf.worst = std::max(nf.worst, err);
                if (err > tol + 1e-9 * zmax) ++nf.over;
            }
    }
    return nf;
}

struct LineFindings {
    std::size_t unplaced = 0;      // an output constraint edge on no input segment
    std::size_t wrong_mask = 0;    // on one, with another mask
    std::size_t broken_chain = 0;  // an input segment whose pieces do not chain end to end
};

// `start_world` and `edges`, `masks`: the input segments as given.
template <typename Outcome, typename Edges, typename Masks>
LineFindings line_findings(const terrain::raster::RasterGeometry& g, const std::vector<Point2>& start_world,
                           const Edges& edges, const Masks& masks, const Outcome& out) {
    LineFindings lf;
    std::vector<std::vector<std::pair<double, double>>> pieces(edges.size());
    for (std::size_t k = 0; k < out.edges.size(); ++k) {
        const Frac p = frac(g, out.vertices[out.edges[k][0]]), q = frac(g, out.vertices[out.edges[k][1]]);
        std::optional<std::size_t> seg;
        for (std::size_t i = 0; i < edges.size() && !seg; ++i) {
            const Frac a = frac(g, start_world[edges[i][0]]), b = frac(g, start_world[edges[i][1]]);
            const auto [tp, dp] = param_dist(a, b, p);
            const auto [tq, dq] = param_dist(a, b, q);
            if (dp <= kOnLine && dq <= kOnLine && std::min(tp, tq) >= -kOnLine && std::max(tp, tq) <= 1.0 + kOnLine) {
                seg = i;
                pieces[i].push_back(std::minmax(tp, tq));
            }
        }
        if (!seg) {
            ++lf.unplaced;
            continue;
        }
        if (out.masks[k] != masks[*seg]) ++lf.wrong_mask;
    }
    for (auto& ps : pieces) {
        std::sort(ps.begin(), ps.end());
        bool ok = !ps.empty() && std::abs(ps.front().first) <= kOnLine && std::abs(ps.back().second - 1.0) <= kOnLine;
        for (std::size_t j = 1; ok && j < ps.size(); ++j) ok = std::abs(ps[j].first - ps[j - 1].second) <= kOnLine;
        if (!ok) ++lf.broken_chain;
    }
    return lf;
}

// The output vertex at (col, row), to `eps` cells; empty if none.
template <typename Outcome>
std::optional<std::size_t> vertex_at(const terrain::raster::RasterGeometry& g, const Outcome& out, Frac at,
                                     double eps = 1e-9) {
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const Frac f = frac(g, out.vertices[i]);
        if (std::abs(f.col - at.col) <= eps && std::abs(f.row - at.row) <= eps) return i;
    }
    return std::nullopt;
}

}  // namespace constraint_foot_oracles
