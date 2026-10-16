#pragma once

// Test-only oracles for increment 33 (docs/increments/33-feature-tolerance.md,
// sections 3 and 9): the ramp t(d), a point's and a triangle's distance to a
// set of segments by brute force, and the DEM-node guarantee G3 with a
// tolerance that varies from node to node.
//
// Written from section 3's text, never from line_tolerance.hpp: the ramp is
// section 3's formula; a distance is the plain Euclidean one in world metres
// (a point to a segment, clamped projection), computed here with no index and
// no cap. The node oracle is strip_oracle::node_findings' membership and plane
// (exact orientation on (col, -row), barycentric weights) with the uniform
// tolerance replaced by t(d(n)), d(n) the node's distance to the ORIGINAL
// segments (G3). It reads no record of the run and returns its findings, so a
// caller can show it fails.
//
// No Catch2 include.

#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace line_tolerance_oracle {

using terrain::Point2;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using Seg = std::array<double, 4>;  // x0 y0 x1 y1, world metres (the binding's row layout)

// Section 3: N, F, S, E and the margin, in metres.
struct Ramp {
    double near;
    double far;
    double start;
    double end;
    double margin = 0.0;
};

// Section 3's t(max(0, d - margin)). S = E is a step: N up to S, F after.
inline double ramp(const Ramp& r, double d) {
    const double x = std::max(0.0, d - r.margin);
    if (x <= r.start) return r.near;
    if (x >= r.end) return r.far;
    return r.near + (r.far - r.near) * (x - r.start) / (r.end - r.start);
}

// Euclidean distance from p to the closed segment (a, b); a point for a == b.
inline double point_segment(Point2 p, Point2 a, Point2 b) {
    const double ux = b.x - a.x, uy = b.y - a.y, len2 = ux * ux + uy * uy;
    double t = len2 > 0.0 ? ((p.x - a.x) * ux + (p.y - a.y) * uy) / len2 : 0.0;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(p.x - (a.x + t * ux), p.y - (a.y + t * uy));
}

inline double point_lines(Point2 p, std::span<const Seg> segs) {
    double best = std::numeric_limits<double>::infinity();
    for (const Seg& s : segs) best = std::min(best, point_segment(p, Point2{s[0], s[1]}, Point2{s[2], s[3]}));
    return best;
}

inline double orient(Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }

// p in the closed triangle (a, b, c), CCW, by double orientation: the oracle's
// threshold test, not a topology decision (section 4.2 says the same of the
// producer's).
inline bool in_closed(Point2 a, Point2 b, Point2 c, Point2 p) {
    return orient(a, b, p) >= 0.0 && orient(b, c, p) >= 0.0 && orient(c, a, p) >= 0.0;
}

// Proper or touching intersection of closed segments (p, q) and (r, s).
inline bool segments_meet(Point2 p, Point2 q, Point2 r, Point2 s) {
    const double d1 = orient(r, s, p), d2 = orient(r, s, q), d3 = orient(p, q, r), d4 = orient(p, q, s);
    if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) return true;
    auto on = [](Point2 a, Point2 b, Point2 c) {  // c on closed (a, b), given collinear
        return std::min(a.x, b.x) <= c.x && c.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= c.y
            && c.y <= std::max(a.y, b.y);
    };
    return (d1 == 0 && on(r, s, p)) || (d2 == 0 && on(r, s, q)) || (d3 == 0 && on(p, q, r))
        || (d4 == 0 && on(p, q, s));
}

// The closed triangle (a, b, c), CCW, to the closed segment (p, q): 0 when they
// meet, else the smallest of the corner-to-segment and end-to-edge distances.
inline double triangle_segment(Point2 a, Point2 b, Point2 c, Point2 p, Point2 q) {
    if (in_closed(a, b, c, p) || in_closed(a, b, c, q) || segments_meet(a, b, p, q) || segments_meet(b, c, p, q)
        || segments_meet(c, a, p, q))
        return 0.0;
    return std::min({point_segment(a, p, q), point_segment(b, p, q), point_segment(c, p, q),
                     point_segment(p, a, b), point_segment(p, b, c), point_segment(p, c, a),
                     point_segment(q, a, b), point_segment(q, b, c), point_segment(q, c, a)});
}

// What the DEM-node guarantee found over one output.
struct RampFindings {
    std::size_t nodes = 0;    // (node, triangle) pairs checked
    std::size_t over = 0;     // pairs whose error is over t(d(n)) + slack
    std::size_t not_ccw = 0;
    double worst_excess = -std::numeric_limits<double>::infinity();  // max(err - t(d(n)))
};

// G3: every valid DEM node in every CLOSED output triangle with three valid
// vertices is within t(d(n)) of the triangle's plane, d(n) its world distance
// to `segs` (the original segments). Slack 1e-9 max(1, |z|max), as
// node_findings'.
template <class Mesh>
RampFindings ramp_findings(const Raster<float>& dem, const Mesh& out, std::span<const Seg> segs, const Ramp& r) {
    using terrain::pred::DefaultKernel;
    using terrain::pred::Orientation;
    const RasterGeometry& g = dem.geometry();
    std::vector<Point2> fp;
    double zmax = 0.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const strip_oracle::Lat l = strip_oracle::lat(g, out.vertices[i]);
        fp.push_back(Point2{l.col, -l.row});
        if (out.valid[i]) zmax = std::max(zmax, std::abs(out.z[i]));
    }
    RampFindings f;
    for (const auto& tri : out.triangles) {
        const Point2 a = fp[tri[0]], b = fp[tri[1]], c = fp[tri[2]];
        if (DefaultKernel::orient2d(a, b, c) != Orientation::CounterClockwise) {
            ++f.not_ccw;
            continue;
        }
        if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
        const double two_a = orient(a, b, c);
        const auto lo_c = static_cast<std::int64_t>(std::ceil(std::min({a.x, b.x, c.x})));
        const auto hi_c = static_cast<std::int64_t>(std::floor(std::max({a.x, b.x, c.x})));
        const auto lo_r = static_cast<std::int64_t>(std::ceil(-std::max({a.y, b.y, c.y})));
        const auto hi_r = static_cast<std::int64_t>(std::floor(-std::min({a.y, b.y, c.y})));
        for (std::int64_t row = std::max<std::int64_t>(lo_r, 0); row <= hi_r; ++row)
            for (std::int64_t col = std::max<std::int64_t>(lo_c, 0); col <= hi_c; ++col) {
                const Point2 p{static_cast<double>(col), -static_cast<double>(row)};
                if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                    continue;
                const auto ur = static_cast<std::size_t>(row), uc = static_cast<std::size_t>(col);
                if (p == a || p == b || p == c || strip_oracle::nodata(dem, ur, uc)) continue;
                const double plane = (orient(p, b, c) * out.z[tri[0]] + orient(a, p, c) * out.z[tri[1]]
                                      + orient(a, b, p) * out.z[tri[2]]) / two_a;
                const double err = std::abs(strip_oracle::at(dem, ur, uc) - plane);
                const double allowed = ramp(r, point_lines(g.node({ur, uc}), segs));
                ++f.nodes;
                f.worst_excess = std::max(f.worst_excess, err - allowed);
                if (err > allowed + 1e-9 * std::max(1.0, zmax)) ++f.over;
            }
    }
    return f;
}

}  // namespace line_tolerance_oracle
