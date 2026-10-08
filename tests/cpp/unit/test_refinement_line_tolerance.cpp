// Increment 33 (docs/increments/33-feature-tolerance.md, sections 3, 4.2, 8
// and 9): the tolerance field on its own, before refine sees it.
//
//   test 4   the ramp ToleranceRamp::at: 0, S, the middle, E and beyond; the
//            step S = E; the margin; never decreasing over 10 000 distances.
//   test 6   LineTolerance::distance: the eight named triangle-to-segment
//            cases, the search's stop (M4), and 10 000 random triangles whose
//            indexed distance equals the brute-force minimum exactly, with
//            and without the cap.
//   test 10  LineTolerance::make's refusals, one each, and what it accepts.
//
// Invariant-critical (section 9): tests 4 and 6. Mutants M1 (the ramp from E,
// or no margin), M2 (centroid distance), M3 (a crossing segment with both
// ends outside gets a positive distance) and M4 (stop at best <= 2g) are
// killed here.
//
// Interface used, as section 4.2 writes it (namespace terrain::refinement):
//   ToleranceRamp{near, far, start, end, margin}, .at(distance)
//   UniformTolerance{value}, the TolerancePolicy concept
//   LineTolerance::make(geometry, span<const array<double, 4>>, ramp, why)
//       -> optional<LineTolerance>; lowest(), highest(), at(m, t), distance(m, t)
//
// "Exactly" in test 6 is read in the producer's own arithmetic: the brute
// force is the minimum over one LineTolerance per segment (one segment, so no
// index can lose it), each capped at end + margin as distance() is. The pair
// distance itself is checked against this file's own formula
// (support/line_tolerance_oracle.hpp) to 1e-9 m, its coordinates below 1e3 m.
//
// Guarded: until include/terrain/refinement/line_tolerance.hpp exists, this
// file compiles to one failing case that says so, and the rest of the tree
// builds.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#if __has_include(<terrain/refinement/line_tolerance.hpp>)

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/line_tolerance.hpp>

#include "line_tolerance_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <span>
#include <string>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::mesh::LatticeMesh;
using terrain::mesh::MeshVertex;
using terrain::raster::RasterGeometry;
using terrain::refinement::LineTolerance;
using terrain::refinement::ToleranceRamp;
using terrain::refinement::TolerancePolicy;
using terrain::refinement::UniformTolerance;
using Seg = line_tolerance_oracle::Seg;

static_assert(TolerancePolicy<UniformTolerance>);
static_assert(TolerancePolicy<LineTolerance>);
static_assert(!TolerancePolicy<double>);

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// 1 m cells with world (0, 0) at node (5000, 5000), so every case reads in
// metres and every vertex of this file lies inside the node rectangle.
RasterGeometry unit_geometry() { return RasterGeometry{-5000.0, 5000.0, 1.0, 1.0, 20001, 20001}; }

Point2 world(const RasterGeometry& g, const MeshVertex& v) {
    return Point2{g.x_min() + v.col * g.delta_x(), g.y_max() - v.row * g.delta_y()};
}

MeshVertex vertex_at(const RasterGeometry& g, Point2 p) {
    return MeshVertex{(p.x - g.x_min()) / g.delta_x(), (g.y_max() - p.y) / g.delta_y()};
}

// One mesh of disjoint triangles, each given by its world corners, CCW in world.
LatticeMesh disjoint(const RasterGeometry& g, const std::vector<std::array<Point2, 3>>& tris) {
    std::vector<MeshVertex> v;
    std::vector<terrain::TriangleIndices> t;
    for (const auto& tri : tris) {
        const auto i = static_cast<std::uint32_t>(v.size());
        for (const Point2 p : tri) v.push_back(vertex_at(g, p));
        t.push_back({i, i + 1, i + 2});
    }
    const std::size_t n = t.size();
    auto m = LatticeMesh::build(std::move(v), std::move(t), std::vector<std::uint8_t>(n, 0),
                                std::vector<std::array<std::uint32_t, 3>>(n, {0, 0, 0}));
    REQUIRE(m.has_value());
    return std::move(*m);
}

LineTolerance field(const RasterGeometry& g, const std::vector<Seg>& segs, ToleranceRamp r) {
    std::string why;
    auto f = LineTolerance::make(g, std::span<const Seg>{segs}, r, why);
    INFO(why);
    REQUIRE(f.has_value());
    return std::move(*f);
}

std::string refusal(const RasterGeometry& g, const std::vector<Seg>& segs, ToleranceRamp r) {
    std::string why;
    const auto f = LineTolerance::make(g, std::span<const Seg>{segs}, r, why);
    CHECK_FALSE(f.has_value());
    CHECK_FALSE(why.empty());
    return why;
}

// The ramp the named cases use: 1 m on the line to 20 m at 1 km, no margin.
constexpr ToleranceRamp kWide{1.0, 20.0, 0.0, 1000.0, 0.0};

// The triangle of the named cases: A (0, 0), B (10, 0), C (0, 10) in metres.
const std::array<Point2, 3> kT{Point2{0.0, 0.0}, Point2{10.0, 0.0}, Point2{0.0, 10.0}};

double distance_to(const Seg& s, ToleranceRamp r = kWide) {
    const auto g = unit_geometry();
    const auto m = disjoint(g, {kT});
    return field(g, {s}, r).distance(m, 0);
}

}  // namespace

// ------------------------------------------------------------------- test 4

TEST_CASE("33 test 4: the ramp at 0, S, the middle, E and beyond", "[line_tolerance][ramp]") {
    const ToleranceRamp r{1.0, 20.0, 100.0, 3100.0, 0.0};
    CHECK(r.at(0.0) == 1.0);
    CHECK(r.at(50.0) == 1.0);
    CHECK(r.at(100.0) == 1.0);
    CHECK(r.at(1600.0) == Catch::Approx(10.5).epsilon(1e-12));  // (S + E) / 2: (N + F) / 2
    CHECK(r.at(850.0) == Catch::Approx(1.0 + 19.0 * 0.25).epsilon(1e-12));
    CHECK(r.at(3100.0) == 20.0);
    CHECK(r.at(3101.0) == 20.0);
    CHECK(r.at(1e9) == 20.0);
}

TEST_CASE("33 test 4: Ola's ramp, 1 m on the line to 20 m at 3 km", "[line_tolerance][ramp]") {
    // Section 3's figures: 1 m on the line, 1.63 m at 100 m, 4.2 m at 500 m.
    const ToleranceRamp r{1.0, 20.0, 0.0, 3000.0, 0.0};
    CHECK(r.at(0.0) == 1.0);
    CHECK(r.at(100.0) == Catch::Approx(1.0 + 19.0 / 30.0).epsilon(1e-12));
    CHECK(r.at(500.0) == Catch::Approx(1.0 + 19.0 / 6.0).epsilon(1e-12));
    CHECK(r.at(3000.0) == 20.0);
}

TEST_CASE("33 test 4: the step S = E evaluates nothing between", "[line_tolerance][ramp]") {
    const ToleranceRamp r{0.5, 8.0, 250.0, 250.0, 0.0};
    CHECK(r.at(0.0) == 0.5);
    CHECK(r.at(250.0) == 0.5);  // d <= S: N
    CHECK(r.at(std::nextafter(250.0, kInf)) == 8.0);
    CHECK(r.at(1000.0) == 8.0);
    CHECK(std::isfinite(r.at(250.0)));
}

TEST_CASE("33 test 4: the margin is subtracted from every distance, never below 0", "[line_tolerance][ramp]") {
    // M1: a ramp that drops the margin is 1.0 + 19 * 2 / 3000 above these.
    const ToleranceRamp r{1.0, 20.0, 0.0, 3000.0, 2.0};
    CHECK(r.at(0.0) == 1.0);
    CHECK(r.at(1.0) == 1.0);
    CHECK(r.at(2.0) == 1.0);
    CHECK(r.at(1502.0) == Catch::Approx(10.5).epsilon(1e-12));
    CHECK(r.at(3001.0) == Catch::Approx(20.0 - 19.0 / 3000.0).epsilon(1e-12));
    CHECK(r.at(3002.0) == 20.0);
    const ToleranceRamp s{1.0, 20.0, 100.0, 3100.0, 2.0};
    CHECK(s.at(102.0) == 1.0);
    CHECK(s.at(103.0) == Catch::Approx(1.0 + 19.0 / 3000.0).epsilon(1e-12));
}

TEST_CASE("33 test 4: the ramp interpolates from S, not from E", "[line_tolerance][ramp]") {
    // M1's first form: at S + a quarter of (E - S), a quarter of the way to F.
    const ToleranceRamp r{2.0, 10.0, 400.0, 800.0, 0.0};
    CHECK(r.at(500.0) == Catch::Approx(4.0).epsilon(1e-12));
    CHECK(r.at(700.0) == Catch::Approx(8.0).epsilon(1e-12));
}

TEST_CASE("33 test 4: the ramp never decreases, a sweep of 10 000 distances", "[line_tolerance][ramp]") {
    const ToleranceRamp r = GENERATE(ToleranceRamp{1.0, 20.0, 0.0, 3000.0, 1.0}, ToleranceRamp{0.0, 5.0, 10.0, 10.0, 0.0},
                                     ToleranceRamp{3.0, 3.0, 0.0, 50.0, 0.5}, ToleranceRamp{0.25, 7.0, 33.0, 900.0, 4.0});
    CAPTURE(r.near, r.far, r.start, r.end, r.margin);
    double last = -kInf;
    for (int i = 0; i <= 10000; ++i) {
        const double d = 1.2 * (r.end + r.margin) * i / 10000.0;
        const double t = r.at(d);
        CAPTURE(d, t);
        REQUIRE(t >= last);
        REQUIRE(t >= r.near);
        REQUIRE(t <= r.far);
        REQUIRE(t == Catch::Approx(line_tolerance_oracle::ramp({r.near, r.far, r.start, r.end, r.margin}, d))
                         .epsilon(1e-12)
                         .margin(1e-12));
        last = t;
    }
}

TEST_CASE("33: UniformTolerance is today's single number", "[line_tolerance][policy]") {
    const UniformTolerance u{2.5};
    const auto g = unit_geometry();
    const auto m = disjoint(g, {kT});
    CHECK(u.lowest() == 2.5);
    CHECK(u.highest() == 2.5);
    CHECK(u.at(m, 0) == 2.5);
}

TEST_CASE("33: lowest and highest bound the field", "[line_tolerance][policy]") {
    const auto g = unit_geometry();
    const auto m = disjoint(g, {kT});
    const auto with = field(g, {Seg{0.0, -3.0, 10.0, -3.0}}, kWide);
    CHECK(with.lowest() == 1.0);
    CHECK(with.highest() == 20.0);
    const auto none = field(g, {}, kWide);  // zero segments: every triangle gets far (section 8)
    CHECK(none.lowest() == 20.0);
    CHECK(none.highest() == 20.0);
    CHECK(none.at(m, 0) == 20.0);
}

// ------------------------------------------------------------------- test 6

TEST_CASE("33 test 6: a segment crossing the triangle with both ends outside is at 0", "[line_tolerance][distance]") {
    // M3.
    CHECK(distance_to(Seg{-5.0, 2.0, 15.0, 2.0}) == 0.0);
    CHECK(distance_to(Seg{-3.0, 9.0, 9.0, -3.0}) == 0.0);  // through two edges, no corner
}

TEST_CASE("33 test 6: a segment with an end inside is at 0", "[line_tolerance][distance]") {
    CHECK(distance_to(Seg{2.0, 2.0, 20.0, 20.0}) == 0.0);
    CHECK(distance_to(Seg{1.0, 1.0, 2.0, 1.5}) == 0.0);  // both ends inside
}

TEST_CASE("33 test 6: a segment touching an edge or a corner is at 0", "[line_tolerance][distance]") {
    CHECK(distance_to(Seg{5.0, 0.0, 5.0, -10.0}) == 0.0);    // an end on edge AB
    CHECK(distance_to(Seg{-4.0, 4.0, 4.0, -4.0}) == 0.0);    // through corner A only
    CHECK(distance_to(Seg{10.0, 0.0, 20.0, 0.0}) == 0.0);    // from corner B outward
}

TEST_CASE("33 test 6: a segment along an edge is at 0", "[line_tolerance][distance]") {
    CHECK(distance_to(Seg{2.0, 0.0, 8.0, 0.0}) == 0.0);
    CHECK(distance_to(Seg{-5.0, 0.0, 15.0, 0.0}) == 0.0);
    CHECK(distance_to(Seg{12.0, -2.0, -2.0, 12.0}) == 0.0);  // along BC's line, past both ends
}

TEST_CASE("33 test 6: a parallel segment at a known offset", "[line_tolerance][distance]") {
    CHECK(distance_to(Seg{0.0, -3.0, 10.0, -3.0}) == Catch::Approx(3.0).epsilon(1e-12));
    CHECK(distance_to(Seg{-0.5, 0.0, -0.5, 10.0}) == Catch::Approx(0.5).epsilon(1e-12));
    CHECK(distance_to(Seg{12.0, 2.0, 2.0, 12.0}) == Catch::Approx(4.0 / std::sqrt(2.0)).epsilon(1e-12));
}

TEST_CASE("33 test 6: nearest at a corner, not at the centroid", "[line_tolerance][distance]") {
    // M2: the centroid (10/3, 10/3) is 9.7 m from this segment, corner A 5 m.
    CHECK(distance_to(Seg{-3.0, -4.0, -7.0, -1.0}) == Catch::Approx(5.0).epsilon(1e-12));
    CHECK(distance_to(Seg{13.0, -4.0, 13.0, -40.0}) == Catch::Approx(5.0).epsilon(1e-12));  // corner B
}

TEST_CASE("33 test 6: nearest at an end of the segment", "[line_tolerance][distance]") {
    CHECK(distance_to(Seg{5.0, -2.0, 5.0, -12.0}) == Catch::Approx(2.0).epsilon(1e-12));
    CHECK(distance_to(Seg{-1.5, 5.0, -30.0, 5.0}) == Catch::Approx(1.5).epsilon(1e-12));
}

TEST_CASE("33 test 6: a zero-length segment is a point", "[line_tolerance][distance]") {
    CHECK(distance_to(Seg{13.0, 13.0, 13.0, 13.0}) == Catch::Approx(16.0 / std::sqrt(2.0)).epsilon(1e-12));
    CHECK(distance_to(Seg{-3.0, -4.0, -3.0, -4.0}) == Catch::Approx(5.0).epsilon(1e-12));
    CHECK(distance_to(Seg{2.0, 3.0, 2.0, 3.0}) == 0.0);  // inside
}

TEST_CASE("33 test 6: the distance is capped at end + margin", "[line_tolerance][distance]") {
    const ToleranceRamp r{1.0, 20.0, 0.0, 100.0, 2.0};
    CHECK(distance_to(Seg{0.0, -500.0, 10.0, -500.0}, r) == 102.0);
    CHECK(distance_to(Seg{0.0, -50.0, 10.0, -50.0}, r) == Catch::Approx(50.0).epsilon(1e-12));
}

TEST_CASE("33 test 6: at() is the ramp at the triangle's distance", "[line_tolerance][distance]") {
    const auto g = unit_geometry();
    const auto m = disjoint(g, {kT});
    const ToleranceRamp r{1.0, 20.0, 0.0, 100.0, 1.0};
    const auto f = field(g, {Seg{0.0, -21.0, 10.0, -21.0}}, r);
    CHECK(f.distance(m, 0) == Catch::Approx(21.0).epsilon(1e-12));
    CHECK(f.at(m, 0) == Catch::Approx(1.0 + 19.0 * 20.0 / 100.0).epsilon(1e-12));
}

TEST_CASE("33 test 6: distances are world metres at UTM scale, non-square cells", "[line_tolerance][distance]") {
    // dx 10, dy 5: a row/col swap or a dropped offset moves every answer.
    const RasterGeometry g{500000.0, 7000000.0, 10.0, 5.0, 100, 100};
    // Nodes (row 10, col 10), (20, 10), (10, 20): world (500100, 6999950),
    // (500100, 6999900), (500200, 6999950); CCW in world.
    const std::array<Point2, 3> t{world(g, MeshVertex{10.0, 10.0}), world(g, MeshVertex{10.0, 20.0}),
                                  world(g, MeshVertex{20.0, 10.0})};
    REQUIRE(t[0].x == 500100.0);
    REQUIRE(t[0].y == 6999950.0);
    const auto m = disjoint(g, {t});
    const auto f = field(g, {Seg{500050.0, 6999957.0, 500300.0, 6999957.0}}, kWide);
    CHECK(f.distance(m, 0) == Catch::Approx(7.0).epsilon(1e-9));
    const auto side = field(g, {Seg{500097.0, 6999800.0, 500097.0, 6999990.0}}, kWide);
    CHECK(side.distance(m, 0) == Catch::Approx(3.0).epsilon(1e-9));
}

TEST_CASE("33 test 6: the search stops only when the best is within the distance its box grew by",
          "[line_tolerance][distance][search]") {
    // M4. E = 1600, so the first query grows the triangle's box by g = 100 m.
    // P, a zero-length segment 92 m left of and 92 m below corner A, is 130.1 m
    // away and inside that box; Q, a vertical segment 110 m right of corner B,
    // is 110 m away and outside it. The answer is Q's. A search that stops when
    // the best is within 2g returns P's.
    //
    // Twelve far segments (1 885 m and more away, beyond the cap) set the
    // index's extent; the scene is shifted along x in 64 steps of 7.3 m so that,
    // for any bucket size above about 10 m, some shift puts P's and Q's
    // buckets apart and a false positive cannot hide the mutant. At shift 0
    // and k = ceil(sqrt(14)) = 4 buckets over x in [3115, 7115], a bucket edge
    // lies at x = 5115, between the first box's x_max 5110 and Q at 5120.
    const auto g = unit_geometry();
    const ToleranceRamp r{1.0, 20.0, 0.0, 1600.0, 0.0};
    for (int k = 0; k < 64; ++k) {
        const double sx = 7.3 * k;
        CAPTURE(k, sx);
        const std::array<Point2, 3> tri{Point2{5000.0 + sx, -5000.0}, Point2{5010.0 + sx, -5000.0},
                                        Point2{5000.0 + sx, -4990.0}};
        const auto m = disjoint(g, {tri});
        std::vector<Seg> segs{Seg{4908.0 + sx, -5092.0, 4908.0 + sx, -5092.0},
                              Seg{5120.0 + sx, -5000.0, 5120.0 + sx, -4990.0}};
        for (int i = 0; i < 6; ++i) {
            const double y = -7000.0 + 800.0 * i;
            segs.push_back(Seg{3115.0, y, 3115.0, y});
            segs.push_back(Seg{7115.0, y, 7115.0, y});
        }
        const auto f = field(g, segs, r);
        REQUIRE(f.distance(m, 0) == Catch::Approx(110.0).epsilon(1e-12));
    }
}

namespace {

struct RandomCase {
    std::uint32_t seed;
    std::size_t segments;
    ToleranceRamp ramp;
};

double uniform(std::mt19937& gen, double lo, double hi) {
    return lo + (hi - lo) * static_cast<double>(gen() % 1000001u) / 1000000.0;
}

}  // namespace

TEST_CASE("33 test 6: 10 000 random triangles, the indexed distance is the brute-force minimum exactly",
          "[line_tolerance][distance][random]") {
    // Ten sets of 1 000 triangles in a 400 x 400 cell window of a non-square
    // grid (dx 1.5 m, dy 0.75 m), against 0 to 300 random segments (a tenth
    // of them zero-length, lengths to 60 m). E = 16 m caps most triangles
    // (the cap's path); E = 2 000 m caps none.
    const RandomCase c = GENERATE(RandomCase{1, 0, {1.0, 20.0, 0.0, 2000.0, 0.0}},
                                  RandomCase{2, 1, {1.0, 20.0, 0.0, 2000.0, 0.0}},
                                  RandomCase{3, 5, {1.0, 20.0, 0.0, 16.0, 0.0}},
                                  RandomCase{4, 5, {1.0, 20.0, 0.0, 2000.0, 1.5}},
                                  RandomCase{5, 50, {1.0, 20.0, 0.0, 16.0, 1.5}},
                                  RandomCase{6, 50, {0.5, 5.0, 10.0, 2000.0, 0.0}},
                                  RandomCase{7, 300, {1.0, 20.0, 0.0, 16.0, 0.0}},
                                  RandomCase{8, 300, {1.0, 20.0, 0.0, 2000.0, 0.0}},
                                  RandomCase{9, 120, {1.0, 20.0, 0.0, 40.0, 1.0}},
                                  RandomCase{10, 300, {2.0, 2.0, 0.0, 300.0, 0.0}});
    CAPTURE(c.seed, c.segments, c.ramp.end, c.ramp.margin);
    const RasterGeometry g{1000.0, 2000.0, 1.5, 0.75, 401, 401};
    std::mt19937 gen{c.seed};
    std::vector<Seg> segs;
    for (std::size_t i = 0; i < c.segments; ++i) {
        const double x = uniform(gen, 1000.0, 1600.0), y = uniform(gen, 1700.0, 2000.0);
        const bool point = gen() % 10u == 0u;
        const double len = point ? 0.0 : uniform(gen, 0.0, 60.0), ang = uniform(gen, 0.0, 6.283185307179586);
        segs.push_back(Seg{x, y, x + len * std::cos(ang), y + len * std::sin(ang)});
    }
    std::vector<std::array<Point2, 3>> tris;
    while (tris.size() < 1000) {
        std::array<Point2, 3> t{};
        const double cx = uniform(gen, 1010.0, 1590.0), cy = uniform(gen, 1710.0, 1990.0);
        const double size = uniform(gen, 0.1, 20.0);
        for (auto& p : t) p = Point2{cx + uniform(gen, -size, size) / 2, cy + uniform(gen, -size, size) / 2};
        const double o = line_tolerance_oracle::orient(t[0], t[1], t[2]);
        if (std::abs(o) < 1e-3) continue;
        if (o < 0) std::swap(t[1], t[2]);
        tris.push_back(t);
    }
    const auto m = disjoint(g, tris);
    const auto all = field(g, segs, c.ramp);
    std::vector<LineTolerance> one;
    for (const Seg& s : segs) one.push_back(field(g, {s}, c.ramp));
    const double cap = c.ramp.end + c.ramp.margin;
    std::size_t capped = 0;
    for (std::uint32_t t = 0; t < tris.size(); ++t) {
        double brute = cap;  // no segment: the cap
        for (const auto& f : one) brute = std::min(brute, f.distance(m, t));
        CAPTURE(t);
        REQUIRE(all.distance(m, t) == brute);  // bit for bit
        REQUIRE(all.at(m, t) == c.ramp.at(brute));
        // The pair distance against this file's own formula, world metres.
        const auto tw = std::array<Point2, 3>{world(g, m.vertices()[m.triangles()[t][0]]),
                                              world(g, m.vertices()[m.triangles()[t][1]]),
                                              world(g, m.vertices()[m.triangles()[t][2]])};
        double mine = cap;
        for (const Seg& s : segs)
            mine = std::min(mine, line_tolerance_oracle::triangle_segment(tw[0], tw[1], tw[2], Point2{s[0], s[1]},
                                                                          Point2{s[2], s[3]}));
        REQUIRE(brute == Catch::Approx(mine).margin(1e-9));
        capped += brute == cap ? 1 : 0;
    }
    // Both paths are exercised where they are meant to be.
    if (c.segments == 0) CHECK(capped == tris.size());
    if (c.ramp.end == 2000.0 && c.segments > 0) CHECK(capped == 0);
    if (c.ramp.end == 16.0 && c.segments >= 50) {
        CHECK(capped > 0);
        CHECK(capped < tris.size());
    }
}

// ------------------------------------------------------------------ test 10

TEST_CASE("33 test 10: make refuses each broken bound with a reason", "[line_tolerance][make]") {
    const auto g = unit_geometry();
    const std::vector<Seg> segs{Seg{0.0, -3.0, 10.0, -3.0}};
    CHECK(refusal(g, segs, {-0.1, 20.0, 0.0, 100.0, 0.0}).find("near") != std::string::npos);
    CHECK(refusal(g, segs, {21.0, 20.0, 0.0, 100.0, 0.0}).find("near") != std::string::npos);
    CHECK(refusal(g, segs, {1.0, 20.0, -1.0, 100.0, 0.0}).find("start") != std::string::npos);
    CHECK(refusal(g, segs, {1.0, 20.0, 200.0, 100.0, 0.0}).find("start") != std::string::npos);
    CHECK(refusal(g, segs, {1.0, 20.0, 0.0, 100.0, -0.5}).find("margin") != std::string::npos);
}

TEST_CASE("33 test 10: make refuses every non-finite ramp value", "[line_tolerance][make]") {
    const auto g = unit_geometry();
    const std::vector<Seg> segs{Seg{0.0, -3.0, 10.0, -3.0}};
    const double bad = GENERATE(kNaN, kInf, -kInf);
    CAPTURE(bad);
    refusal(g, segs, {bad, 20.0, 0.0, 100.0, 0.0});
    refusal(g, segs, {1.0, bad, 0.0, 100.0, 0.0});
    refusal(g, segs, {1.0, 20.0, bad, 100.0, 0.0});
    refusal(g, segs, {1.0, 20.0, 0.0, bad, 0.0});
    refusal(g, segs, {1.0, 20.0, 0.0, 100.0, bad});
}

TEST_CASE("33 test 10: make refuses a non-finite coordinate in any place", "[line_tolerance][make]") {
    const auto g = unit_geometry();
    const double bad = GENERATE(kNaN, kInf, -kInf);
    const int place = GENERATE(0, 1, 2, 3);
    CAPTURE(bad, place);
    std::vector<Seg> segs{Seg{0.0, -3.0, 10.0, -3.0}, Seg{1.0, 1.0, 2.0, 2.0}};
    segs[1][static_cast<std::size_t>(place)] = bad;
    CHECK(refusal(g, segs, kWide).find("coordinate") != std::string::npos);
}

TEST_CASE("33 test 10: make accepts the edges of the bounds", "[line_tolerance][make]") {
    const auto g = unit_geometry();
    const std::vector<Seg> segs{Seg{0.0, -3.0, 10.0, -3.0}};
    for (const ToleranceRamp r : {ToleranceRamp{0.0, 0.0, 0.0, 0.0, 0.0},       // N = 0, as --tolerance 0
                                  ToleranceRamp{5.0, 5.0, 10.0, 10.0, 0.0},     // N = F, S = E
                                  ToleranceRamp{1.0, 20.0, 0.0, 3000.0, 1.0}}) {
        std::string why;
        CHECK(LineTolerance::make(g, std::span<const Seg>{segs}, r, why).has_value());
    }
    std::string why;
    CHECK(LineTolerance::make(g, std::span<const Seg>{}, kWide, why).has_value());  // zero segments
}

#else

TEST_CASE("33: include/terrain/refinement/line_tolerance.hpp does not exist yet", "[line_tolerance]") {
    FAIL("increment 33 is not built: no terrain/refinement/line_tolerance.hpp (LineTolerance, "
         "ToleranceRamp, UniformTolerance, TolerancePolicy), so tests 4, 6 and 10 cannot compile");
}

#endif
