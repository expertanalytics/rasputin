#pragma once

// J2's tolerance oracle (docs/increments/15c-geographic-dem.md, J2), shared by
// ES9 (property/prop_refinement_edge_strip.cpp) and RP3 and RP5
// (property/prop_refinement_refine_points.cpp). Harness h17 §4a
// (docs/increments/h17-ci-test-time.md) replaced the two copies those files
// held, each walking every (triangle, check point) pair, with this one, which
// files the check points by DEM cell and walks only the cells a triangle's box
// covers (uniform-grid bucketing; the grid is the DEM's own lattice). The
// verdict is the copies' verdict: the same pairs are counted.
//
// The oracle: every check point that is not a start vertex, in every CLOSED
// output triangle with three valid vertices (exact orientation in the frame
// given), is within tolerance of that triangle's plane, the plane recomputed
// here from the z given (the output's, or a planted copy). It reads no scan
// record and RETURNS its findings, so a caller can plant a defect and show it
// fails.
//
// Frame: (col, -row) in lattice units, as both callers build it from world
// points. A check point outside [0, cols-1] x [0, rows-1] (or NaN) throws
// std::invalid_argument: no point is silently left out. No Catch2 include.

#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace j2_oracle {

using terrain::Point2;

struct Findings {
    std::size_t over = 0;     // (check point, triangle) pairs over tolerance
    std::size_t not_ccw = 0;  // triangles not counter-clockwise, valid or not
};

inline double cross(Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }

namespace detail {

// floor(v) clamped to [0, last] (§4a rule 6). NaN maps to the far end of the
// range on its side. That does not make a triangle with a NaN corner visit
// every bucket: std::min and std::max over the corners return NaN only when the
// NaN comes first among them (min({NaN, 1, 2}) is NaN, min({1, NaN, 2}) is 1),
// so which buckets it visits depends on the corner's position. The verdict is
// the all-pairs copies' all the same, for another reason: a NaN corner makes
// two_a and the plane NaN, and abs(NaN - z) > slack is false, so such a
// triangle counts no pair over tolerance however many points it visits;
// not_ccw is decided before the bucket walk, for every triangle.
inline std::size_t bucket_lo(double v, std::size_t last) {
    const double f = std::floor(v);
    if (!(f >= 0.0)) return 0;
    if (f >= static_cast<double>(last)) return last;
    return static_cast<std::size_t>(f);
}
inline std::size_t bucket_hi(double v, std::size_t last) {
    const double f = std::floor(v);
    if (!(f <= static_cast<double>(last))) return last;
    if (f <= 0.0) return 0;
    return static_cast<std::size_t>(f);
}

}  // namespace detail

// `skip[i]` is 1 when check point i is a start vertex (J2 excludes it); the
// caller decides it once per point (§4a rule 1).
inline Findings violations(std::span<const Point2> vertices, std::span<const double> z,
                           std::span<const std::uint8_t> valid,
                           std::span<const std::array<std::uint32_t, 3>> triangles,
                           std::span<const Point2> points, std::span<const float> point_z,
                           std::span<const std::uint8_t> skip, std::size_t cols, std::size_t rows,
                           double tol) {
    using terrain::pred::DefaultKernel;
    using terrain::pred::Orientation;
    if (z.size() != vertices.size() || valid.size() != vertices.size())
        throw std::invalid_argument("j2_oracle: z and valid must hold one entry per vertex");
    if (point_z.size() != points.size() || skip.size() != points.size())
        throw std::invalid_argument("j2_oracle: point_z and skip must hold one entry per check point");
    if (cols == 0 || rows == 0) throw std::invalid_argument("j2_oracle: the grid has no node");
    for (const auto& t : triangles)
        for (const auto v : t)
            if (v >= vertices.size()) throw std::invalid_argument("j2_oracle: a triangle names no vertex");

    // Rule 2: file the points by (floor(col), floor(row)), row = -y, as a
    // compressed table: start[b] .. start[b + 1] are bucket b's points.
    const std::size_t last_col = cols - 1, last_row = rows - 1;
    std::vector<std::size_t> bucket(points.size());
    std::vector<std::size_t> start(cols * rows + 1, 0);
    for (std::size_t i = 0; i < points.size(); ++i) {
        const double col = points[i].x, row = -points[i].y;
        if (!(col >= 0.0 && col <= static_cast<double>(last_col) && row >= 0.0
              && row <= static_cast<double>(last_row)))
            throw std::invalid_argument("j2_oracle: check point " + std::to_string(i) + " at (col "
                                        + std::to_string(col) + ", row " + std::to_string(row)
                                        + ") is outside the " + std::to_string(cols) + " x "
                                        + std::to_string(rows) + " grid");
        bucket[i] = static_cast<std::size_t>(std::floor(row)) * cols + static_cast<std::size_t>(std::floor(col));
        ++start[bucket[i] + 1];
    }
    for (std::size_t b = 0; b < cols * rows; ++b) start[b + 1] += start[b];
    std::vector<std::size_t> filed(points.size());
    {
        std::vector<std::size_t> next(start.begin(), start.end() - 1);
        for (std::size_t i = 0; i < points.size(); ++i) filed[next[bucket[i]]++] = i;
    }

    // Rule 4: the scale of the slack, over the z given.
    double zmax = 1.0;
    for (const double v : z) zmax = std::max(zmax, std::abs(v));

    Findings f;
    for (const auto& t : triangles) {
        const Point2 a = vertices[t[0]], b = vertices[t[1]], c = vertices[t[2]];
        if (DefaultKernel::orient2d(a, b, c) != Orientation::CounterClockwise) ++f.not_ccw;  // rule 5
        if (!(valid[t[0]] && valid[t[1]] && valid[t[2]])) continue;
        // Rule 3: the buckets from floor of the box's low corner to floor of
        // its high corner, inclusive; rule 6: clamped to the grid.
        const double lo_col = std::min({a.x, b.x, c.x}), hi_col = std::max({a.x, b.x, c.x});
        const double lo_row = -std::max({a.y, b.y, c.y}), hi_row = -std::min({a.y, b.y, c.y});
        const std::size_t c0 = detail::bucket_lo(lo_col, last_col), c1 = detail::bucket_hi(hi_col, last_col);
        const std::size_t r0 = detail::bucket_lo(lo_row, last_row), r1 = detail::bucket_hi(hi_row, last_row);
        const double two_a = cross(a, b, c);
        for (std::size_t r = r0; r <= r1; ++r)
            for (std::size_t col = c0; col <= c1; ++col)
                for (std::size_t k = start[r * cols + col]; k < start[r * cols + col + 1]; ++k) {
                    const std::size_t i = filed[k];
                    if (skip[i]) continue;
                    const Point2 p = points[i];
                    if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                        || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                        || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                        continue;
                    const double plane = (cross(p, b, c) * z[t[0]] + cross(a, p, c) * z[t[1]]
                                          + cross(a, b, p) * z[t[2]]) / two_a;
                    if (std::abs(plane - static_cast<double>(point_z[i])) > tol + 1e-9 * zmax) ++f.over;
                }
    }
    return f;
}

}  // namespace j2_oracle
