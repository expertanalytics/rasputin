// Increment 20b (docs/increments/20b-min-insertion-distance.md): constraint
// feet, a minimum insertion distance from constraint segments scaled by the
// tolerance. The suite the design calls test_refinement_constraint_feet;
// invariant-critical. Registered as test_refinement_constraint_feet from this
// file (property/, since every case goes through refine()).
//
// Interface, as R9 names it:
//
//   RefineOptions::constraint_feet   bool, default false; set by member assignment
//   RefineOutcome::feet              std::size_t, feet inserted, counted in `inserted`
//   RefineOutcome::feet_refused      std::size_t, R2 step 4 refusals (N inserted instead)
//
// CHOSEN HERE, where the design leaves it open:
//   - the output world point of a foot is (x_min + col dx, y_max - row dy) of
//     its stored fractional vertex (R4), so the oracles below map an output
//     point back to (col, row) by the same affine map and allow 1e-9 cells;
//   - its z is compared with bilinear() at the output world point to 1e-9
//     relative, since the stored (col, row) is not visible from outside;
//   - F5's increment-20 references with the quality pass on were recorded from
//     517e0e3's refine.hpp (increment 20, unchanged) with refine_digest's
//     topology digest, before any production change. The pass-off references
//     are Q7's, which are increment 18's = increment 20's with the pass off.
//
// Oracles are independent of scan.hpp and refine.hpp: membership by
// DefaultKernel::orient2d on (col, -row), the plane from barycentric weights,
// every valid DEM node in every triangle with three valid vertices (F1). A
// foot is recognised from the output alone: an off-node vertex that is not a
// start vertex, lying on an input segment, and the orthogonal projection of a
// valid DEM node no further than the cap from it (R2, R3).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/refinement/refine.hpp>

#include "feet_fixtures.hpp"
#include "quality_fixtures.hpp"
#include "refine_digest.hpp"
#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;
using quality_fixtures::Start;

namespace {

struct Frac {
    double col;
    double row;
};

Frac frac(const RasterGeometry& g, Point2 p) {
    return Frac{(p.x - g.x_min()) / g.delta_x(), (g.y_max() - p.y) / g.delta_y()};
}

bool is_node(const RasterGeometry& g, Point2 p) {
    const Frac f = frac(g, p);
    const double c = std::round(f.col), r = std::round(f.row);
    if (!(c >= 0 && r >= 0 && c < static_cast<double>(g.cols()) && r < static_cast<double>(g.rows())))
        return false;
    return g.node({static_cast<std::size_t>(r), static_cast<std::size_t>(c)}) == p;
}

bool nodata(const Raster<float>& dem, std::int64_t r, std::int64_t c) {
    return dem.is_nodata({static_cast<std::size_t>(r), static_cast<std::size_t>(c)});
}

double at(const Raster<float>& dem, std::int64_t r, std::int64_t c) {
    return static_cast<double>(dem.value_at({static_cast<std::size_t>(r), static_cast<std::size_t>(c)}));
}

auto run(const Raster<float>& dem, const Start& s, double tol, bool feet, double min_angle = 0.0,
         unsigned threads = 1) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.min_angle_deg = min_angle;
    o.constraint_feet = feet;
    return refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks}, o);
}

// The increment-20 call, byte for byte: constraint_feet never named.
auto run20(const Raster<float>& dem, const Start& s, double tol, double min_angle) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = 1;
    o.min_angle_deg = min_angle;
    return refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks}, o);
}

// An input constraint segment in (col, row) with its mask.
struct Segment {
    Frac a;
    Frac b;
    std::uint32_t mask;
};

std::vector<Segment> segments(const RasterGeometry& g, const Start& s) {
    std::vector<Segment> out;
    for (std::size_t i = 0; i < s.edges.size(); ++i)
        out.push_back(Segment{frac(g, s.mesh.vertices()[s.edges[i][0]]),
                              frac(g, s.mesh.vertices()[s.edges[i][1]]), s.masks[i]});
    return out;
}

// Parameter along the segment and perpendicular distance, both in cells.
std::pair<double, double> param_dist(const Segment& s, Frac p) {
    const double ux = s.b.col - s.a.col, uy = s.b.row - s.a.row;
    const double len2 = ux * ux + uy * uy;
    const double t = ((p.col - s.a.col) * ux + (p.row - s.a.row) * uy) / len2;
    const double d = std::abs((p.col - s.a.col) * uy - (p.row - s.a.row) * ux) / std::sqrt(len2);
    return {t, d};
}

constexpr double kOnLine = 1e-9;  // cells: "every foot within 1e-9 cells of the input line"

std::optional<std::size_t> segment_of(const std::vector<Segment>& segs, Frac p) {
    for (std::size_t i = 0; i < segs.size(); ++i) {
        const auto [t, d] = param_dist(segs[i], p);
        if (d <= kOnLine && t >= -kOnLine && t <= 1.0 + kOnLine) return i;
    }
    return std::nullopt;
}

// World distance between two (col, row) points.
double world(const RasterGeometry& g, Frac p, Frac q) {
    return std::hypot((p.col - q.col) * g.delta_x(), (p.row - q.row) * g.delta_y());
}

// The orthogonal projection of (col, row) onto the segment's line, in the
// world frame (R2 step 1: distances in (col dx, row dy)).
Frac project(const RasterGeometry& g, const Segment& s, Frac p) {
    const double ux = (s.b.col - s.a.col) * g.delta_x(), uy = (s.b.row - s.a.row) * g.delta_y();
    const double px = (p.col - s.a.col) * g.delta_x(), py = (p.row - s.a.row) * g.delta_y();
    const double t = (px * ux + py * uy) / (ux * ux + uy * uy);
    return Frac{s.a.col + t * (s.b.col - s.a.col), s.a.row + t * (s.b.row - s.a.row)};
}

double cap(const RasterGeometry& g) { return std::min(g.delta_x(), g.delta_y()) / 2.0; }
double floor_eps(const RasterGeometry& g) { return std::min(g.delta_x(), g.delta_y()) / 100.0; }

// What the output says about its feet, recomputed from the output alone.
struct Feet {
    std::vector<std::uint32_t> vertex;           // output index of each foot
    std::vector<std::size_t> segment;            // the input segment it lies on
    std::vector<std::array<std::int64_t, 2>> source;  // the DEM node it is the foot of (row, col)
};

// Every output vertex is a start vertex as given, a DEM node, or a foot (R4):
// on an input segment strictly inside it, and the projection of a valid DEM
// node within the cap of it. z per R0 and R4. Returns the feet.
template <typename Outcome>
Feet classify(const Raster<float>& dem, const Start& start, const Outcome& out) {
    const RasterGeometry& g = dem.geometry();
    const auto segs = segments(g, start);
    const std::size_t n0 = start.mesh.vertices().size();
    REQUIRE(out.vertices.size() == n0 + out.quality_inserted + out.inserted);
    Feet feet;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        CAPTURE(i);
        const Point2 p = out.vertices[i];
        if (i < n0) {
            REQUIRE(p.x == start.mesh.vertices()[i].x);  // F3: input vertices unchanged
            REQUIRE(p.y == start.mesh.vertices()[i].y);
            continue;
        }
        if (is_node(g, p)) {
            const Frac f = frac(g, p);
            const auto r = std::llround(f.row), c = std::llround(f.col);
            REQUIRE(static_cast<bool>(out.valid[i]) == !nodata(dem, r, c));
            if (out.valid[i]) REQUIRE(out.z[i] == at(dem, r, c));
            continue;
        }
        const Frac f = frac(g, p);
        CAPTURE(f.col, f.row);
        const auto s = segment_of(segs, f);
        REQUIRE(s.has_value());  // an off-node inserted vertex is on a constraint
        const auto [t, d] = param_dist(segs[*s], f);
        REQUIRE(t > 0.0);
        REQUIRE(t < 1.0);
        // The node it replaced: within the cap, and f is its projection.
        std::optional<std::array<std::int64_t, 2>> src;
        const auto rc = static_cast<std::int64_t>(std::floor(f.row)), cc = static_cast<std::int64_t>(std::floor(f.col));
        for (std::int64_t r = rc - 1; r <= rc + 2; ++r)
            for (std::int64_t c = cc - 1; c <= cc + 2; ++c) {
                if (r < 0 || c < 0 || r >= static_cast<std::int64_t>(g.rows())
                    || c >= static_cast<std::int64_t>(g.cols()) || nodata(dem, r, c))
                    continue;
                const Frac node{static_cast<double>(c), static_cast<double>(r)};
                const Frac pr = project(g, segs[*s], node);
                if (world(g, pr, f) <= 1e-7 && world(g, node, f) < cap(g) * (1.0 + 1e-12)) src = {r, c};
            }
        REQUIRE(src.has_value());  // R3: ε never exceeds half a cell
        REQUIRE(out.valid[i] == 1);
        const auto z = terrain::raster::bilinear(dem, p);
        REQUIRE(z.has_value());
        REQUIRE(std::abs(out.z[i] - *z) <= 1e-9 * std::max(1.0, std::abs(*z)));  // R4: bilinear, not a node's
        feet.vertex.push_back(static_cast<std::uint32_t>(i));
        feet.segment.push_back(*s);
        feet.source.push_back(*src);
    }
    REQUIRE(feet.vertex.size() == out.feet);
    return feet;
}

// F1: every valid DEM node in every triangle with three valid vertices is
// within tolerance of the output's own plane, recomputed here, not from scan.
template <typename Outcome>
void tolerance_oracle(const Raster<float>& dem, const Outcome& out, double tol) {
    const RasterGeometry& g = dem.geometry();
    std::vector<Point2> fp;
    double zmax = 0.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const Frac f = frac(g, out.vertices[i]);
        fp.push_back(Point2{f.col, -f.row});
        zmax = std::max(zmax, std::abs(out.z[i]));
    }
    auto cross = [](Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); };
    double worst = 0.0;
    for (const auto& tri : out.triangles) {
        const Point2 a = fp[tri[0]], b = fp[tri[1]], c = fp[tri[2]];
        REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::CounterClockwise);
        if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
        const double two_a = cross(a, b, c);
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
                if (p == a || p == b || p == c || nodata(dem, r, col)) continue;
                const double plane = (cross(p, b, c) * out.z[tri[0]] + cross(a, p, c) * out.z[tri[1]]
                                      + cross(a, b, p) * out.z[tri[2]]) / two_a;
                const double err = std::abs(at(dem, r, col) - plane);
                worst = std::max(worst, err);
                CAPTURE(r, col, err, tol);
                REQUIRE(err <= tol + 1e-9 * std::max(1.0, zmax));
            }
    }
    REQUIRE(out.max_error <= tol);
}

// F3: the constrained edges are exactly the input segments cut at their
// vertices: each piece on one segment with that segment's mask, the pieces of
// each segment chaining from one end to the other. F4's spacing: along a
// segment, a foot is at least `floor` from its neighbours (R2 step 2, R6.2).
template <typename Outcome>
void constraints_oracle(const Raster<float>& dem, const Start& start, const Outcome& out,
                        const Feet& feet) {
    const RasterGeometry& g = dem.geometry();
    const auto segs = segments(g, start);
    std::vector<std::vector<std::pair<double, double>>> pieces(segs.size());  // (t0, t1) sorted
    for (std::size_t k = 0; k < out.edges.size(); ++k) {
        const Frac p = frac(g, out.vertices[out.edges[k][0]]), q = frac(g, out.vertices[out.edges[k][1]]);
        CAPTURE(k, p.col, p.row, q.col, q.row);
        std::optional<std::size_t> s;
        for (std::size_t i = 0; i < segs.size() && !s; ++i) {
            const auto [tp, dp] = param_dist(segs[i], p);
            const auto [tq, dq] = param_dist(segs[i], q);
            if (dp <= kOnLine && dq <= kOnLine && std::min(tp, tq) >= -kOnLine && std::max(tp, tq) <= 1 + kOnLine)
                s = i;
        }
        REQUIRE(s.has_value());
        REQUIRE(out.masks[k] == segs[*s].mask);  // bits and masks on both halves
        const double tp = param_dist(segs[*s], p).first, tq = param_dist(segs[*s], q).first;
        pieces[*s].push_back(std::minmax(tp, tq));
    }
    std::set<std::uint32_t> is_foot(feet.vertex.begin(), feet.vertex.end());
    for (std::size_t i = 0; i < segs.size(); ++i) {
        CAPTURE(i);
        auto& ps = pieces[i];
        std::sort(ps.begin(), ps.end());
        REQUIRE(!ps.empty());
        REQUIRE(std::abs(ps.front().first) <= kOnLine);
        REQUIRE(std::abs(ps.back().second - 1.0) <= kOnLine);
        for (std::size_t j = 1; j < ps.size(); ++j) REQUIRE(std::abs(ps[j].first - ps[j - 1].second) <= kOnLine);
    }
    // Feet are at least floor apart from every other vertex on their segment.
    for (std::size_t f = 0; f < feet.vertex.size(); ++f) {
        const Frac p = frac(g, out.vertices[feet.vertex[f]]);
        for (std::size_t v = 0; v < out.vertices.size(); ++v) {
            if (v == feet.vertex[f]) continue;
            const Frac q = frac(g, out.vertices[v]);
            const auto [t, d] = param_dist(segs[feet.segment[f]], q);
            if (d > kOnLine || t < -kOnLine || t > 1 + kOnLine) continue;
            CAPTURE(f, v);
            REQUIRE(world(g, p, q) >= floor_eps(g) * (1.0 - 1e-9));
        }
    }
}

// 14b R10 with feet on: every interior edge that is not a constraint edge has
// neither apex strictly inside the other triangle's circle, in 14b's
// LatticeFrame (col dx, -(row dy)) on the fractional coordinates recovered from
// the output. A foot's split must be legalised like any other (R2 step 3);
// prop_refinement_refine.cpp's T16 oracle, restated here.
template <typename Outcome>
void delaunay_oracle(const RasterGeometry& g, const Outcome& out) {
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
        REQUIRE(ts.size() <= 2);
        if (ts.size() != 2 || constrained.count(e) != 0) continue;
        for (unsigned s = 0; s < 2; ++s) {
            const auto& tri = out.triangles[ts[s]];
            std::uint32_t apex = 0;
            for (const auto x : out.triangles[ts[1 - s]])
                if (x != e.first && x != e.second) apex = x;
            if (DefaultKernel::incircle(lf(tri[0]), lf(tri[1]), lf(tri[2]), lf(apex)) == terrain::pred::Incircle::Inside) {
                UNSCOPED_INFO("edge " << e.first << "-" << e.second << " apex " << apex);
                ++bad;
            }
        }
    }
    REQUIRE(bad == 0);
}

// `delaunay` is false only for the one known pre-20b defect below.
template <typename Outcome>
Feet check(const Raster<float>& dem, const Start& start, const Outcome& out, double tol,
           bool delaunay = true) {
    REQUIRE(out.ok());
    const Feet feet = classify(dem, start, out);
    tolerance_oracle(dem, out, tol);
    constraints_oracle(dem, start, out, feet);
    if (delaunay) delaunay_oracle(dem.geometry(), out);
    // R2 step 5: a node is footed once, so no two feet share a source node.
    REQUIRE(std::set(feet.source.begin(), feet.source.end()).size() == feet.source.size());
    return feet;
}

bool has_vertex(const RasterGeometry& g, const std::vector<Point2>& vs, Frac p) {
    return std::any_of(vs.begin(), vs.end(), [&](Point2 v) { return world(g, frac(g, v), p) <= 1e-7; });
}

double min_angle_deg(const std::vector<Point2>& v, const terrain::TriangleIndices& t) {
    double m = 180.0;
    for (unsigned k = 0; k < 3; ++k) {
        const Point2 p = v[t[k]], q = v[t[(k + 1) % 3]], r = v[t[(k + 2) % 3]];
        const double ux = q.x - p.x, uy = q.y - p.y, wx = r.x - p.x, wy = r.y - p.y;
        m = std::min(m, std::atan2(std::abs(ux * wy - uy * wx), ux * wx + uy * wy) * 180.0 / std::numbers::pi);
    }
    return m;
}

const Frac kNeedle{static_cast<double>(feet_fixtures::kNeedleCol), static_cast<double>(feet_fixtures::kNeedleRow)};

}  // namespace

// ------------------------------------------------------------------------ F1

TEST_CASE("F1: the tolerance holds at every valid node beside a tilted segment on steep ground",
          "[refinement][feet]") {
    const double tol = GENERATE(0.1, 0.5, 2.0);
    const double steep = GENERATE(0.0, 1.0);
    const double slope = GENERATE(feet_fixtures::kSlope, 1.0 / 7.77, std::numbers::sqrt2 / 40.0);
    const double min_angle = GENERATE(0.0, 25.0);
    CAPTURE(tol, steep, slope, min_angle);
    const auto dem = feet_fixtures::needle_dem(6.0, steep);
    const auto start = feet_fixtures::needle_start(dem.geometry(), slope);
    const auto out = run(dem, start, tol, true, min_angle);
    check(dem, start, out, tol);
}

TEST_CASE("F1: rough ground beside the tilted segment stays within tolerance", "[refinement][feet]") {
    // 14b's T3 re-run with the rule on: a DEM with no structure at all.
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    const double tol = GENERATE(0.0, 5.0);
    CAPTURE(seed, tol);
    const Raster<float> dem{refinement_fixtures::geometry(feet_fixtures::kN, feet_fixtures::kN),
                            refinement_fixtures::rough_dem(feet_fixtures::kN, feet_fixtures::kN, seed)};
    const auto start = feet_fixtures::needle_start(dem.geometry());
    check(dem, start, run(dem, start, tol, true), tol);
}

TEST_CASE("F1: the rule fires on the steep fixture", "[refinement][feet]") {
    // So F1 cannot pass with a rule that never runs.
    const auto dem = feet_fixtures::needle_dem(6.0, 1.0);
    const auto start = feet_fixtures::needle_start(dem.geometry());
    const auto out = run(dem, start, 0.5, true);
    REQUIRE(check(dem, start, out, 0.5).vertex.size() >= 1);
}

// ------------------------------------------------------------------------ F2

TEST_CASE("F2: the needle node is resolved by a foot and no needle is left", "[refinement][feet]") {
    const auto dem = feet_fixtures::needle_dem();
    const RasterGeometry& g = dem.geometry();
    const auto start = feet_fixtures::needle_start(g);
    const auto out = run(dem, start, 0.5, true);
    const Feet feet = check(dem, start, out, 0.5);
    REQUIRE(out.feet >= 1);
    REQUIRE(out.feet_refused == 0);  // R2 step 4: measured 0; nothing here is degenerate

    // The foot is the needle node's projection on A-B, and the node is not a vertex.
    const Frac f = project(g, segments(g, start)[0], kNeedle);
    REQUIRE(world(g, f, kNeedle) < 0.05);  // 3.5 cm: the fixture is what it claims
    REQUIRE(has_vertex(g, out.vertices, f));
    REQUIRE_FALSE(has_vertex(g, out.vertices, kNeedle));
    REQUIRE(std::count(feet.source.begin(), feet.source.end(),
                       std::array<std::int64_t, 2>{feet_fixtures::kNeedleRow, feet_fixtures::kNeedleCol}) == 1);

    // No triangle under 0.1 degrees near the needle, where increment 20 left one.
    for (const auto& t : out.triangles) {
        bool near = false;
        for (const auto v : t) near |= world(g, frac(g, out.vertices[v]), kNeedle) <= 2.0 * std::hypot(g.delta_x(), g.delta_y());
        if (!near) continue;
        CAPTURE(t[0], t[1], t[2]);
        REQUIRE(min_angle_deg(out.vertices, t) >= 0.1);
    }
}

TEST_CASE("F2: without feet the fixture does make the needle", "[refinement][feet]") {
    // The converse: F2's fixture is the defect, so F2 is not passing vacuously.
    const auto dem = feet_fixtures::needle_dem();
    const RasterGeometry& g = dem.geometry();
    const auto start = feet_fixtures::needle_start(g);
    const auto out = run(dem, start, 0.5, false);
    REQUIRE(out.ok());
    REQUIRE(out.feet == 0);
    REQUIRE(has_vertex(g, out.vertices, kNeedle));
    double worst = 180.0;
    for (const auto& t : out.triangles) worst = std::min(worst, min_angle_deg(out.vertices, t));
    REQUIRE(worst < 0.1);
}

TEST_CASE("F2: a foot in a NoData cell is refused and the node inserted instead", "[refinement][feet]") {
    // R2 step 4: vertex_z refuses F (a NoData corner in F's cell), N itself is
    // valid, so N goes in as today and the refusal is counted.
    const auto dem0 = feet_fixtures::needle_dem();
    const RasterGeometry g = dem0.geometry();
    std::vector<float> z;
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (const float v : dem0.row(r)) z.push_back(v);
    const float hole = std::numeric_limits<float>::quiet_NaN();
    // Column 7 beside the needle is outside the domain; its nodes are corners of F's cell.
    for (std::int64_t r = feet_fixtures::kNeedleRow - 1; r <= feet_fixtures::kNeedleRow + 1; ++r)
        z[static_cast<std::size_t>(r) * feet_fixtures::kN + 7] = hole;
    const Raster<float> dem{g, std::move(z)};
    const auto start = feet_fixtures::needle_start(g);
    const auto out = run(dem, start, 0.5, true);
    check(dem, start, out, 0.5);
    REQUIRE(out.feet_refused == 1);  // measured: the needle alone is refused
    REQUIRE(has_vertex(g, out.vertices, kNeedle));
}

// ------------------------------------------------------------------------ F3

TEST_CASE("F3: constraints stay covered with their masks on every piece", "[refinement][feet]") {
    const double tol = GENERATE(0.0, 0.5);
    const double slope = GENERATE(feet_fixtures::kSlope, 1.0 / 7.77);
    CAPTURE(tol, slope);
    const auto dem = feet_fixtures::needle_dem(6.0, 0.3);
    const auto start = feet_fixtures::needle_start(dem.geometry(), slope);
    const auto out = run(dem, start, tol, true);
    const Feet feet = check(dem, start, out, tol);
    REQUIRE(!feet.vertex.empty());
    // Every foot splits A-B: the side has more pieces than it started with.
    std::size_t ab = 0;
    for (const auto m : out.masks) ab += m == 1u ? 1 : 0;
    REQUIRE(ab >= 1 + static_cast<std::size_t>(std::count(feet.segment.begin(), feet.segment.end(), 0u)));
}

// ------------------------------------------------------------------------ F4

TEST_CASE("F4: at tolerance 0 the fallback inserts the footed node and the run ends",
          "[refinement][feet]") {
    const auto dem = feet_fixtures::needle_dem();
    const RasterGeometry& g = dem.geometry();
    const auto start = feet_fixtures::needle_start(g);
    const auto out = run(dem, start, 0.0, true);
    // Not Delaunay here, with feet off as well: see the [!shouldfail] case below.
    const Feet feet = check(dem, start, out, 0.0, false);  // spacing >= floor, one foot per node
    REQUIRE(out.max_error == 0.0);
    REQUIRE(out.feet >= 1);
    // The needle got its foot, and then went in anyway: R5's fallback.
    REQUIRE(has_vertex(g, out.vertices, project(g, segments(g, start)[0], kNeedle)));
    REQUIRE(has_vertex(g, out.vertices, kNeedle));
    // R6.2: one foot per node (check() above), so no more feet than nodes.
    REQUIRE(feet.source.size() <= g.rows() * g.cols());
}

TEST_CASE("Known defect, not 20b's: at tolerance 0 the needle fixture is not constrained Delaunay",
          "[refinement][feet]") {
    // Found by the incircle oracle when it was added to check(). Feet off or on,
    // the output has the unconstrained interior edge (27, 2)-(27, 7) in
    // (col, row) with (28, 4) strictly inside the circle of (27, 2) (26, 6)
    // (27, 7) in the LatticeFrame: exact incircle determinant 12500. All four
    // are lattice nodes, nowhere near a constraint, and the quad is convex, so
    // Lawson should have flipped it. At 0.5 m the same fixture is Delaunay.
    // [!shouldfail]: this case goes red when the defect is fixed, so the
    // exemption in F4 above is removed with it.
    const auto dem = feet_fixtures::needle_dem();
    const auto start = feet_fixtures::needle_start(dem.geometry());
    const bool feet = GENERATE(false, true);
    CAPTURE(feet);
    delaunay_oracle(dem.geometry(), run(dem, start, 0.0, feet));
}

// ------------------------------------------------------------------------ F5

TEST_CASE("F5: the default RefineOptions leave feet off", "[refinement][feet]") {
    REQUIRE_FALSE(RefineOptions{}.constraint_feet);
}

TEST_CASE("F5: feet off reproduces increment 20 on the T12 fixtures and a domain start",
          "[refinement][feet]") {
    // Q7's goldens (quality off) and 517e0e3's (quality on); see the header.
    struct Case {
        const char* name;
        double min_angle;
        std::uint64_t golden;
    };
    const auto c = GENERATE(Case{"cone", 0.0, 0x2326c8c0a462ac4dull}, Case{"island", 0.0, 0x39c68151cb4e1186ull},
                            Case{"ring", 0.0, 0x2b9d05a87b859d47ull}, Case{"cone", 25.0, 0x2326c8c0a462ac4dull},
                            Case{"island", 25.0, 0x39c68151cb4e1186ull}, Case{"ring", 25.0, 0x95b958bd85d2c5f9ull});
    CAPTURE(c.name, c.min_angle);
    const std::string name = c.name;
    auto both = [&](const Raster<float>& dem, const Start& start, double tol) {
        const auto before = run20(dem, start, tol, c.min_angle);
        const auto off = run(dem, start, tol, false, c.min_angle);
        REQUIRE(refine_digest::topology_digest(off) == c.golden);
        REQUIRE(refine_digest::digest(off) == refine_digest::digest(before));
        REQUIRE(off.feet == 0);
        REQUIRE(off.feet_refused == 0);
    };
    if (name == "ring") {
        const std::size_t n = 33;
        const Raster<float> dem{refinement_fixtures::geometry(n, n), refinement_fixtures::rough_dem(n, n, 11)};
        both(dem, quality_fixtures::fan(dem.geometry(), quality_fixtures::q7_ring()), 3.0);
    } else {
        const auto dem = name == "island" ? quality_fixtures::island(129, 2.0, 2.0) : quality_fixtures::cone(129, 2.0, 2.0);
        const auto grid = refinement_fixtures::grid_mesh(dem.geometry(), 16);
        both(dem, Start{grid.mesh, grid.edges, grid.masks}, 1.0);
    }
}

TEST_CASE("F5: feet on changes the needle fixture's output", "[refinement][feet]") {
    const auto dem = feet_fixtures::needle_dem();
    const auto start = feet_fixtures::needle_start(dem.geometry());
    REQUIRE(refine_digest::digest(run(dem, start, 0.5, true)) != refine_digest::digest(run(dem, start, 0.5, false)));
}

// ---------------------------------------------------------------- 14b T6, feet on

TEST_CASE("T6 with feet: bit-identical for 1 2 7 and all threads", "[refinement][feet]") {
    const auto dem = feet_fixtures::needle_dem(6.0, 1.0);
    const auto start = feet_fixtures::needle_start(dem.geometry(), 1.0 / 7.77);
    const auto ref = run(dem, start, 0.1, true, 25.0, 1);
    REQUIRE(ref.ok());
    REQUIRE(ref.feet > 0);
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    const auto other = run(dem, start, 0.1, true, 25.0, threads);
    REQUIRE(refine_digest::digest(other) == refine_digest::digest(ref));
    REQUIRE(other.feet == ref.feet);
    REQUIRE(other.feet_refused == ref.feet_refused);
}

// ---------------------------------------------------------------- the helpers, directly
//
// detail::foot_epsilon, foot_of and foot_fits on hand-built inputs. Each case
// names the mutant it was checked against; see the commit that added them.

namespace {

using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::MeshVertex;

// refinement_fixtures::geometry: dx = 10 m, dy = 5 m, so cap = 2.5 m and floor = 0.05 m.
Raster<float> plane(std::size_t rows, std::size_t cols, double per_x, std::optional<float> nodata = std::nullopt) {
    const auto g = refinement_fixtures::geometry(rows, cols);
    std::vector<float> z;
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) z.push_back(static_cast<float>(per_x * static_cast<double>(c) * g.delta_x()));
    return Raster<float>{g, std::move(z), nodata};
}

}  // namespace

TEST_CASE("R3: foot_epsilon is the cap where G is 0, for any tolerance", "[refinement][feet][helpers]") {
    const auto dem = plane(5, 5, 0.0);
    const double tol = GENERATE(0.0, 1e-3, 1.0, 1e6);
    CAPTURE(tol);
    REQUIRE(terrain::refinement::detail::foot_epsilon(dem, LatticeVertex{2, 2}, tol) == cap(dem.geometry()));
    REQUIRE(terrain::refinement::detail::foot_epsilon(dem, LatticeVertex{0, 0}, tol) == cap(dem.geometry()));  // a corner node: one cell
}

TEST_CASE("R3: foot_epsilon is tol / G clamped to [floor, cap]", "[refinement][feet][helpers]") {
    // z = x: G = 1 exactly (10 m per 10 m column, 0 per row).
    const auto dem = plane(5, 5, 1.0);
    const RasterGeometry& g = dem.geometry();
    using terrain::refinement::detail::foot_epsilon;
    REQUIRE(foot_epsilon(dem, LatticeVertex{2, 2}, 1.0) == 1.0);        // inside the band: tol / G
    REQUIRE(foot_epsilon(dem, LatticeVertex{2, 2}, 100.0) == cap(g));   // clamped down
    REQUIRE(foot_epsilon(dem, LatticeVertex{2, 2}, 1e-4) == floor_eps(g));  // clamped up
    REQUIRE(foot_epsilon(dem, LatticeVertex{2, 2}, 0.0) == floor_eps(g));   // tolerance 0 on a slope: the floor
}

TEST_CASE("R3: a cell with a NoData corner is left out of G", "[refinement][feet][helpers]") {
    // Flat, except node (1, 1) is the sentinel: the four cells around it each
    // have that corner, so G takes none of them. If it did, -9999 m over 10 m
    // would put ε at the floor. Node (2, 2)'s cells include (1, 1)'s cell too.
    const float sentinel = -9999.0f;
    auto dem0 = plane(5, 5, 0.0);
    std::vector<float> z(25, 0.0f);
    z[1 * 5 + 1] = sentinel;
    const Raster<float> dem{dem0.geometry(), std::move(z), sentinel};
    using terrain::refinement::detail::foot_epsilon;
    REQUIRE(foot_epsilon(dem, LatticeVertex{2, 2}, 1.0) == cap(dem.geometry()));
    REQUIRE(foot_epsilon(dem, LatticeVertex{1, 1}, 1.0) == cap(dem.geometry()));  // every cell has it: G = 0
}

TEST_CASE("R2 step 2: no foot within eps of a segment end", "[refinement][feet][helpers]") {
    // One triangle (col, row): A (6, 0.8), B (1.8, 0.8), C (1, 3); A-B is
    // constrained and runs along row 0.8, 1 m from row 1. Flat DEM, so eps is
    // the cap, 2.5 m. The node (row 1, col 4) projects 22 m from B: a foot.
    // The node (row 1, col 2) projects 2 m from B, inside eps: no foot (N goes in).
    const auto dem = plane(6, 8, 0.0);
    auto m = LatticeMesh::build(std::vector<MeshVertex>{{6.0, 0.8}, {1.8, 0.8}, {1.0, 3.0}}, {{0, 1, 2}}, {1u}, {{{7u, 0u, 0u}}});
    REQUIRE(m.has_value());
    using terrain::refinement::detail::foot_of;

    const auto far = foot_of(dem, *m, 0, LatticeVertex{1, 4}, 1.0);
    REQUIRE(far.has_value());
    REQUIRE(far->edge == 0u);
    REQUIRE(std::abs(far->at.col - 4.0) <= 1e-12);  // the projection of N
    REQUIRE(far->at.row == 0.8);

    REQUIRE_FALSE(foot_of(dem, *m, 0, LatticeVertex{1, 2}, 1.0).has_value());  // M2

    // The same at A's end: A moved to (4.1, 0.8), the node at col 4 is 1 m from it.
    auto m2 = LatticeMesh::build(std::vector<MeshVertex>{{4.1, 0.8}, {1.8, 0.8}, {1.0, 3.0}}, {{0, 1, 2}}, {1u}, {{{7u, 0u, 0u}}});
    REQUIRE(m2.has_value());
    REQUIRE_FALSE(foot_of(dem, *m2, 0, LatticeVertex{1, 4}, 1.0).has_value());  // M2
    REQUIRE(foot_of(dem, *m2, 0, LatticeVertex{1, 3}, 1.0).has_value());          // the control, 9 m from A
}

TEST_CASE("R2 step 4: foot_fits refuses a foot that folds either side of the edge", "[refinement][feet][helpers]") {
    // t = (A, B, C) and its neighbour u = (B, A, D) across the constrained
    // edge A-B, both slivers: A (0, 2), B (4, 2), C (2, 1.999), D (2, 2.001)
    // in (col, row). A candidate foot a hundredth of a row off A-B folds
    // whichever sliver it is off towards.
    auto m = LatticeMesh::build(std::vector<MeshVertex>{{0.0, 2.0}, {4.0, 2.0}, {2.0, 1.999}, {2.0, 2.001}},
                                {{0, 1, 2}, {1, 0, 3}}, {1u, 1u}, {{{1u, 0u, 0u}, {1u, 0u, 0u}}});
    REQUIRE(m.has_value());
    REQUIRE(m->neighbours(0)[0] == 1u);
    using terrain::refinement::detail::foot_fits;

    REQUIRE(foot_fits(*m, 0, 0, MeshVertex{1.0, 2.0}));         // on the edge: fits
    REQUIRE_FALSE(foot_fits(*m, 0, 0, MeshVertex{1.0, 1.99}));  // M5: past C, t's side folds
    REQUIRE_FALSE(foot_fits(*m, 0, 0, MeshVertex{1.0, 2.01}));  // M8: past D, only u's side folds
    REQUIRE_FALSE(foot_fits(*m, 1, 0, MeshVertex{1.0, 1.99}));  // from u: t is the neighbour now

    // With no neighbour only t's side is asked: past D is then fine.
    auto lone = LatticeMesh::build(std::vector<MeshVertex>{{0.0, 2.0}, {4.0, 2.0}, {2.0, 1.999}}, {{0, 1, 2}}, {1u}, {{{1u, 0u, 0u}}});
    REQUIRE(lone.has_value());
    REQUIRE(foot_fits(*lone, 0, 0, MeshVertex{1.0, 2.01}));
    REQUIRE_FALSE(foot_fits(*lone, 0, 0, MeshVertex{1.0, 1.99}));
}
