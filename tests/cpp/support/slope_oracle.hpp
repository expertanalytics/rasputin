#pragma once

// Test-only fixtures and oracles for increment 34
// (docs/increments/34-slope-tolerance.md, sections 3, 7 and 9): a vertical
// tolerance that follows the slope of the DEM at every node.
//
// Independence. Written from section 3's text and the probes' fixtures
// (docs/increments/34-probes/fixture_figures.py, mutant_fixtures.py,
// border_check.py), never from steepness.hpp or slope_tolerance.hpp: Horn's
// slope with section 3's fill of missing neighbours is computed here in long
// double; the class is the smallest half degree at or above it; the ramp is
// section 3's t(s); a cell's class is the largest of its valid corners, 0 when
// none is. The node and point oracles borrow the producer's predicates (exact
// orientation in (col, -row), barycentric weights, CheckPoints' filing rule)
// and never its records. They return findings, so a caller can show they
// fail.
//
// Frame. Every fixture here has x_min = y_max = 0, so a world point is
// (col dx, -row dy) and "metres off the north-west node" is the world point
// itself.
//
// No Catch2 include.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include "line_tolerance_oracle.hpp"
#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace slope_oracle {

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::raster::CellIndex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;

inline constexpr float kNoData = -9999.0f;
inline constexpr std::uint8_t kNoDataClass = 255;  // section 4.5's reserved class
inline constexpr long double kBoundary = 1e-9L;    // degrees: test 2's "within 1e-9 of a class boundary"

// ---------------------------------------------------------------- fixtures

inline RasterGeometry geometry(std::size_t rows, std::size_t cols, double dx, double dy) {
    return RasterGeometry{0.0, 0.0, dx, dy, cols, rows};
}

inline double radians(double deg) { return deg * std::numbers::pi / 180.0; }

// Section 9's valley, as fixture_figures.py's valley(): a flat floor
// (x < 200 m), a 40-degree wall to x = 440 m, a plateau, plus 3 m bumps
// 3 sin(2 pi y / 80) sin(2 pi x / 130), x = col dx and y = row dy; computed
// in double, stored as float. `holes` are NoData nodes (V1n's).
inline Raster<float> valley(std::size_t rows, std::size_t cols, double dx, double dy,
                            std::span<const CellIndex> holes = {}) {
    std::vector<float> z(rows * cols);
    const double wall = std::tan(radians(40.0));
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) {
            const double x = static_cast<double>(c) * dx, y = static_cast<double>(r) * dy;
            const double base = wall * std::clamp(x - 200.0, 0.0, 240.0);
            const double bumps = 3.0 * std::sin(2.0 * std::numbers::pi * y / 80.0)
                               * std::sin(2.0 * std::numbers::pi * x / 130.0);
            z[r * cols + c] = static_cast<float>(base + bumps);
        }
    for (const CellIndex h : holes) z[h.row * cols + h.col] = kNoData;
    return Raster<float>{geometry(rows, cols, dx, dy), std::move(z), kNoData};
}

// V1: 65 x 65 nodes at 10 m. V2: the same valley on micro_tiff's grid, 33 rows
// by 41 columns, 10 m by 5 m. V1n: V1 with one NoData node on the wall, row 32,
// column 30 (x = 300 m).
inline Raster<float> v1() { return valley(65, 65, 10.0, 10.0); }
inline Raster<float> v2() { return valley(33, 41, 10.0, 5.0); }
inline constexpr CellIndex kV1nHole{32, 30};
inline Raster<float> v1n() {
    const std::array<CellIndex, 1> hole{kV1nHole};
    return valley(65, 65, 10.0, 10.0, hole);
}

// border_check.py's plane(): z rising at `deg` towards `azimuth_deg` (0: +x,
// 90: north, -row), x = col dx, y = -row dy; stored as float.
inline Raster<float> plane(std::size_t rows, std::size_t cols, double dx, double dy, double deg,
                           double azimuth_deg, std::span<const CellIndex> holes = {}) {
    std::vector<float> z(rows * cols);
    const double t = std::tan(radians(deg)), a = radians(azimuth_deg);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) {
            const double x = static_cast<double>(c) * dx, y = -static_cast<double>(r) * dy;
            z[r * cols + c] = static_cast<float>(t * (x * std::cos(a) + y * std::sin(a)));
        }
    for (const CellIndex h : holes) z[h.row * cols + h.col] = kNoData;
    return Raster<float>{geometry(rows, cols, dx, dy), std::move(z), kNoData};
}

// A start mesh in world coordinates: a counter-clockwise ring fanned from its
// first vertex, every side a constraint edge with mask 1.
inline strip_oracle::Start fan(const std::vector<Point2>& ring) {
    std::vector<terrain::TriangleIndices> tris;
    for (std::uint32_t i = 1; i + 1 < ring.size(); ++i) tris.push_back({0, i, i + 1});
    strip_oracle::Start s;
    const auto n = static_cast<std::uint32_t>(ring.size());
    for (std::uint32_t i = 0; i < n; ++i) {
        const std::uint32_t j = (i + 1) % n;
        s.edges.push_back({std::min(i, j), std::max(i, j)});
        s.masks.push_back(1);
    }
    const std::size_t nt = tris.size();
    s.mesh = IndexedMesh2{ring, std::move(tris), std::vector<std::uint8_t>(nt, 0)};
    return s;
}

// The node rectangle's outline of a fixture (its four corner nodes).
inline strip_oracle::Start outline(const RasterGeometry& g) {
    const double x1 = g.x_max(), y0 = g.y_min();
    return fan({{0.0, y0}, {x1, y0}, {x1, 0.0}, {0.0, 0.0}});
}

// V1c: V1's outline with its south-east corner cut by the edge from
// (640, -300) to (150, -640) m, across the wall, between nodes along its length
// (mutant_fixtures.py's RING, without the closing vertex).
inline constexpr Point2 kCutA{640.0, -300.0}, kCutB{150.0, -640.0};
inline strip_oracle::Start v1c_start() { return fan({{0.0, -640.0}, kCutB, kCutA, {640.0, 0.0}, {0.0, 0.0}}); }

// V1's line (tests 5 and 7): (100, -50) to (600, -600) m; N = 1, ramp 0 to 200 m.
inline std::vector<line_tolerance_oracle::Seg> v1_line() { return {{100.0, -50.0, 600.0, -600.0}}; }

// ---------------------------------------------------------------- Horn, section 3

// Horn's slope in degrees at every node, missing neighbours (outside the grid,
// or NoData) filled as section 3 says: an edge neighbour by its reflection
// through the node, 2 z(node) - z(opposite), when the opposite is there, else
// z(node); a corner neighbour by its reflection when the opposite corner is
// there, else z(row neighbour) + z(column neighbour) - z(node), from the edge
// neighbours as filled. NaN at a NoData node. `dx_for_both` is mutant M4's
// Horn (dx for both axes), for showing the fixtures can kill it.
inline std::vector<long double> horn(const Raster<float>& dem, bool dx_for_both = false) {
    const RasterGeometry& g = dem.geometry();
    const auto rows = static_cast<std::ptrdiff_t>(g.rows()), cols = static_cast<std::ptrdiff_t>(g.cols());
    const long double dx = g.delta_x(), dy = dx_for_both ? g.delta_x() : g.delta_y();
    std::vector<long double> out(g.size(), std::numeric_limits<long double>::quiet_NaN());
    for (std::ptrdiff_t r = 0; r < rows; ++r)
        for (std::ptrdiff_t c = 0; c < cols; ++c) {
            const auto get = [&](std::ptrdiff_t dr, std::ptrdiff_t dc) -> std::optional<long double> {
                const std::ptrdiff_t rr = r + dr, cc = c + dc;
                if (rr < 0 || cc < 0 || rr >= rows || cc >= cols) return std::nullopt;
                const CellIndex at{static_cast<std::size_t>(rr), static_cast<std::size_t>(cc)};
                if (dem.is_nodata(at)) return std::nullopt;
                return static_cast<long double>(dem.value_at(at));
            };
            const auto z = get(0, 0);
            if (!z) continue;
            const auto edge = [&](std::ptrdiff_t dr, std::ptrdiff_t dc) {
                if (const auto v = get(dr, dc)) return *v;
                if (const auto o = get(-dr, -dc)) return 2.0L * *z - *o;
                return *z;
            };
            const auto corner = [&](std::ptrdiff_t dr, std::ptrdiff_t dc) {
                if (const auto v = get(dr, dc)) return *v;
                if (const auto o = get(-dr, -dc)) return 2.0L * *z - *o;
                return edge(dr, 0) + edge(0, dc) - *z;
            };
            const long double a = corner(-1, -1), b = edge(-1, 0), cc = corner(-1, 1);
            const long double d = edge(0, -1), f = edge(0, 1);
            const long double gg = corner(1, -1), h = edge(1, 0), i = corner(1, 1);
            const long double gx = ((cc + 2 * f + i) - (a + 2 * d + gg)) / (8 * dx);
            const long double gy = ((gg + 2 * h + i) - (a + 2 * b + cc)) / (8 * dy);
            out[static_cast<std::size_t>(r * cols + c)] =
                std::atan(std::hypot(gx, gy)) * 180.0L / std::numbers::pi_v<long double>;
        }
    return out;
}

// The smallest half degree at or above `deg`, as a class (section 3); 255 for NaN.
inline std::uint8_t class_of(long double deg) {
    if (std::isnan(deg)) return kNoDataClass;
    return static_cast<std::uint8_t>(std::ceil(2.0L * deg));
}

// Within kBoundary degrees of a class boundary (a multiple of half a degree).
// A slope of exactly 0 (gx = gy = 0, exact in any arithmetic) is class 0 only.
inline bool near_boundary(long double deg) {
    return !std::isnan(deg) && deg != 0.0L && std::abs(2.0L * deg - std::round(2.0L * deg)) <= 2.0L * kBoundary;
}

inline std::vector<std::uint8_t> classes(const Raster<float>& dem, bool dx_for_both = false) {
    std::vector<std::uint8_t> out;
    for (const long double s : horn(dem, dx_for_both)) out.push_back(class_of(s));
    return out;
}

// Test 2's acceptance: the producer's class is the oracle's, or one above
// where the oracle's slope lies within kBoundary of a class boundary; never
// lower. NoData: 255 exactly.
inline bool class_ok(std::uint8_t got, long double deg) {
    const std::uint8_t want = class_of(deg);
    return got == want || (near_boundary(deg) && got == want + 1);
}

// ---------------------------------------------------------------- the ramp

// Section 3's t(s): N from END up, F up to START (and below END), linear between.
struct Ramp {
    double near, far, start, end;
    [[nodiscard]] double at(double s) const {
        if (s >= end) return near;
        if (s <= start) return far;
        return far + (near - far) * (s - start) / (end - start);
    }
    [[nodiscard]] double of_class(std::uint8_t c) const { return c <= 180 ? at(c / 2.0) : far; }
};

// The largest class of the valid corners of the cell CheckPoints files lattice
// point (col, row) in: (floor(row), floor(col)), clamped to the last cell;
// 0 when no corner is valid (section 4.5).
inline std::uint8_t cell_class(const RasterGeometry& g, std::span<const std::uint8_t> cls, double col, double row) {
    const std::size_t lr = std::max<std::size_t>(g.rows(), 2) - 2, lc = std::max<std::size_t>(g.cols(), 2) - 2;
    const std::size_t r0 = std::min(static_cast<std::size_t>(std::floor(row)), lr);
    const std::size_t c0 = std::min(static_cast<std::size_t>(std::floor(col)), lc);
    std::uint8_t best = 0;
    for (const std::size_t r : {r0, r0 + 1})
        for (const std::size_t c : {c0, c0 + 1})
            if (r < g.rows() && c < g.cols() && cls[r * g.cols() + c] != kNoDataClass)
                best = std::max(best, cls[r * g.cols() + c]);
    return best;
}

// ---------------------------------------------------------------- G3, every DEM node

struct NodeFindings {
    std::size_t nodes = 0;    // (node, triangle) pairs checked
    std::size_t over = 0;     // pairs over their own allowed error + slack
    std::size_t not_ccw = 0;
    double worst_excess = -std::numeric_limits<double>::infinity();  // max(err - allowed)
};

// Every valid DEM node in every CLOSED output triangle with three valid
// vertices is within allowed[node] of the triangle's plane recomputed from
// the output (strip_oracle::node_findings' membership and plane). Slack
// 1e-9 max(1, |z|max), as node_findings'.
template <class Mesh>
NodeFindings node_findings(const Raster<float>& dem, const Mesh& out, std::span<const double> allowed) {
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
    NodeFindings f;
    for (const auto& tri : out.triangles) {
        const Point2 a = fp[tri[0]], b = fp[tri[1]], c = fp[tri[2]];
        if (DefaultKernel::orient2d(a, b, c) != Orientation::CounterClockwise) {
            ++f.not_ccw;
            continue;
        }
        if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
        const double two_a = strip_oracle::cross(a, b, c);
        const auto lo_c = static_cast<std::int64_t>(std::ceil(std::min({a.x, b.x, c.x})));
        const auto hi_c = static_cast<std::int64_t>(std::floor(std::max({a.x, b.x, c.x})));
        const auto lo_r = static_cast<std::int64_t>(std::ceil(-std::max({a.y, b.y, c.y})));
        const auto hi_r = static_cast<std::int64_t>(std::floor(-std::min({a.y, b.y, c.y})));
        for (std::int64_t r = std::max<std::int64_t>(lo_r, 0); r <= hi_r; ++r)
            for (std::int64_t col = std::max<std::int64_t>(lo_c, 0); col <= hi_c; ++col) {
                const Point2 p{static_cast<double>(col), -static_cast<double>(r)};
                if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                    continue;
                const auto ur = static_cast<std::size_t>(r), uc = static_cast<std::size_t>(col);
                if (p == a || p == b || p == c || strip_oracle::nodata(dem, ur, uc)) continue;
                const double plane = (strip_oracle::cross(p, b, c) * out.z[tri[0]]
                                      + strip_oracle::cross(a, p, c) * out.z[tri[1]]
                                      + strip_oracle::cross(a, b, p) * out.z[tri[2]]) / two_a;
                const double err = std::abs(strip_oracle::at(dem, ur, uc) - plane);
                const double t = allowed[ur * g.cols() + uc];
                ++f.nodes;
                f.worst_excess = std::max(f.worst_excess, err - t);
                if (err > t + 1e-9 * std::max(1.0, zmax)) ++f.over;
            }
    }
    return f;
}

// Each node's allowed error by the slope alone: the ramp at the oracle's class.
inline std::vector<double> slope_allowed(std::span<const std::uint8_t> cls, const Ramp& r) {
    std::vector<double> out;
    for (const std::uint8_t c : cls) out.push_back(r.of_class(c));
    return out;
}

// With lines too: the smaller of that and 33's ramp at the node's brute-force
// distance to the original segments (G3).
inline std::vector<double> slope_and_lines_allowed(const RasterGeometry& g, std::span<const std::uint8_t> cls,
                                                   const Ramp& r, std::span<const line_tolerance_oracle::Seg> segs,
                                                   const line_tolerance_oracle::Ramp& lines) {
    std::vector<double> out;
    for (std::size_t i = 0; i < cls.size(); ++i) {
        const Point2 p = g.node({i / g.cols(), i % g.cols()});
        out.push_back(std::min(r.of_class(cls[i]),
                               line_tolerance_oracle::ramp(lines, line_tolerance_oracle::point_lines(p, segs))));
    }
    return out;
}

// ---------------------------------------------------------------- G4, check points

struct Sources {
    std::vector<Point2> xy;
    std::vector<float> z;
};

// Four points per cell at seeded dyadic offsets in (0, 1), z the bilinear
// surface of `surface` there plus up to +-4 m of noise in 1/64 m steps, so it
// is a float exactly: 33's scattered() with four points per cell. (Not
// NumPy's draw: test_core_slope_tolerance.py uses scattered(65, 4, 7) itself.)
// `surface` has no NoData; the points are filed on its geometry.
inline Sources scattered(const Raster<float>& surface, std::uint32_t seed) {
    const RasterGeometry& g = surface.geometry();
    std::mt19937 gen{seed};
    Sources p;
    for (std::size_t r = 0; r + 1 < g.rows(); ++r)
        for (std::size_t c = 0; c + 1 < g.cols(); ++c)
            for (int k = 0; k < 4; ++k) {
                const double col = static_cast<double>(c) + static_cast<double>(1u + gen() % 1023u) / 1024.0;
                const double row = static_cast<double>(r) + static_cast<double>(1u + gen() % 1023u) / 1024.0;
                const double noise = (static_cast<double>(gen() % 513u) - 256.0) / 64.0;
                const auto h = strip_oracle::height(surface, strip_oracle::Lat{col, row});
                p.xy.push_back(strip_oracle::world(g, col, row));
                p.z.push_back(static_cast<float>(std::round((*h + noise) * 64.0) / 64.0));
            }
    return p;
}

struct PointFindings {
    std::size_t checked = 0;
    std::size_t at_vertex = 0;  // at an output vertex (an inserted point): not compared
    std::size_t unheld = 0;  // in no output triangle with three valid vertices
    std::size_t over = 0;    // (point, triangle) pairs over the point's allowed error + slack
    double worst_excess = -std::numeric_limits<double>::infinity();
};

// G4: every point in every CLOSED output triangle with three valid vertices is
// within allowed[i] of the triangle's plane, by exact orientation in
// (col, -row). Points at an output vertex are skipped (a vertex's z is the
// point's own or the start's, not a plane value).
template <class Mesh>
PointFindings point_findings(const RasterGeometry& g, const Sources& src, std::span<const double> allowed,
                             const Mesh& out) {
    using terrain::pred::DefaultKernel;
    using terrain::pred::Orientation;
    std::vector<Point2> fp;
    double zmax = 0.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const strip_oracle::Lat l = strip_oracle::lat(g, out.vertices[i]);
        fp.push_back(Point2{l.col, -l.row});
        if (out.valid[i]) zmax = std::max(zmax, std::abs(out.z[i]));
    }
    PointFindings f;
    for (std::size_t k = 0; k < src.xy.size(); ++k) {
        const strip_oracle::Lat l = strip_oracle::lat(g, src.xy[k]);
        const Point2 p{l.col, -l.row};
        if (std::find(fp.begin(), fp.end(), p) != fp.end()) {
            ++f.at_vertex;
            continue;
        }
        ++f.checked;
        bool held = false;
        for (const auto& tri : out.triangles) {
            const Point2 a = fp[tri[0]], b = fp[tri[1]], c = fp[tri[2]];
            if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
            if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                continue;
            held = true;
            const double plane = (strip_oracle::cross(p, b, c) * out.z[tri[0]]
                                  + strip_oracle::cross(a, p, c) * out.z[tri[1]]
                                  + strip_oracle::cross(a, b, p) * out.z[tri[2]])
                               / strip_oracle::cross(a, b, c);
            const double err = std::abs(static_cast<double>(src.z[k]) - plane);
            f.worst_excess = std::max(f.worst_excess, err - allowed[k]);
            if (err > allowed[k] + 1e-9 * std::max(1.0, zmax)) ++f.over;
        }
        f.unheld += held ? 0 : 1;
    }
    return f;
}

// Each point's allowed error: the ramp at its cell's class (the largest valid corner).
inline std::vector<double> point_allowed(const RasterGeometry& g, std::span<const std::uint8_t> cls,
                                         const Sources& src, const Ramp& r) {
    std::vector<double> out;
    for (const Point2 xy : src.xy) {
        const strip_oracle::Lat l = strip_oracle::lat(g, xy);
        out.push_back(r.of_class(cell_class(g, cls, l.col, l.row)));
    }
    return out;
}

}  // namespace slope_oracle
