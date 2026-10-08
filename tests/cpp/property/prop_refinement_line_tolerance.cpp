// Increment 33 (docs/increments/33-feature-tolerance.md, sections 3, 4.3, 4.4,
// 7 and 9): the tolerance field inside the three refinement entry points.
//
//   test 1  refine(..., options) and refine(..., options, UniformTolerance{t})
//           give the same output, bit for bit; the same for refine_points and
//           refine_strip (G1).
//   test 2  LineTolerance with zero segments equals UniformTolerance{F} (G2).
//   test 3  with segments: (a) N = F equals UniformTolerance{F}; (b) a ramp
//           that holds every triangle at N equals UniformTolerance{N}, the
//           LineTolerance run with options.tolerance = F (G2, G8; kills M6
//           and M7).
//   test 5  the guarantee: every valid DEM node within t of its brute-force
//           distance to the original segments (G3).
//   test 7  laziness changes nothing: a policy whose bounds are [0, inf)
//           gives LineTolerance's output, bit for bit (G6; kills M5).
//   test 8  the strip: a constraint edge along a line, every strip point
//           within S of the line ends within N (G4).
//   test 9  1 and 8 threads give the same mesh with a field (G5).
//
// Invariant-critical (section 9): tests 1, 2, 3, 5, 7 and 8 here.
//
// Interface used, as sections 4.2 and 4.3 write it, with the argument order
// pinned where 4.3 says only "the same for refine_points and refine_strip":
// the policy comes right after the options, as in refine, so
//   refine(dem, start, edges, masks, options, policy)
//   refine_points(store, start, z, valid, edges, masks, options, policy, strip = nullptr)
//   refine_strip(dem, strip, start, z, valid, edges, masks, options, policy)
// and RefineOutcome::max_error_near (4.5).
//
// Fixtures: strip_oracle's 33 x 33 exact grid (dx 2 m, dy 1 m; world x in
// [0, 64], y in [-32, 0]) and its jittered start, whose constraint chains run
// off-grid (the outline), along grid row 16 (mask 2), along column 9 (mask 4)
// and through nodes (mask 8); terrain_dem's Smooth and SmoothWithHole (a 3 x 3
// NoData patch). Constraint feet are on wherever the design says "feet on".
//
// Oracles (QA rules, section D): every LineTolerance output below that is not
// already pinned equal to a uniform one carries the constrained-Delaunay
// oracle (strip_oracle::delaunay_violations, exact incircle in the producer's
// frame) and the tolerance oracle with the ramp
// (line_tolerance_oracle::ramp_findings: every valid DEM node in every closed
// triangle with three valid vertices within t(d(n)) of the plane recomputed
// from the output, d(n) by brute force to the original segments). Test 5 and
// test 8 show the ramp oracle fails on the uniform-F mesh of the same input,
// so it can fail.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/line_tolerance.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>

#include "line_tolerance_oracle.hpp"
#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

using namespace strip_oracle;
using line_tolerance_oracle::Ramp;
using line_tolerance_oracle::Seg;
using terrain::refinement::CheckPoints;
using terrain::refinement::constraint_check_points;
using terrain::refinement::ConstraintCheckPoints;
using terrain::refinement::LineTolerance;
using terrain::refinement::PointRefineOptions;
using terrain::refinement::PointRefineOutcome;
using terrain::refinement::refine;
using terrain::refinement::refine_points;
using terrain::refinement::refine_strip;
using terrain::refinement::RefineOptions;
using terrain::refinement::RefineOutcome;
using terrain::refinement::ToleranceRamp;
using terrain::refinement::TolerancePolicy;
using terrain::refinement::UniformTolerance;

namespace {

// ---------------------------------------------------------------- the runs

using Span2 = std::span<const std::array<std::uint32_t, 2>>;
using SpanM = std::span<const std::uint32_t>;

RefineOptions options(double tol, bool feet = true, unsigned threads = 1, double min_angle = 0.0) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.constraint_feet = feet;
    o.min_angle_deg = min_angle;
    return o;
}

PointRefineOptions point_options(double tol, bool feet = true, unsigned threads = 1) {
    PointRefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.constraint_feet = feet;
    return o;
}

RefineOutcome run(const Raster<float>& dem, const Start& s, const RefineOptions& o) {
    return refine(dem, s.mesh, Span2{s.edges}, SpanM{s.masks}, o);
}

template <TolerancePolicy P>
RefineOutcome run(const Raster<float>& dem, const Start& s, const RefineOptions& o, const P& p) {
    return refine(dem, s.mesh, Span2{s.edges}, SpanM{s.masks}, o, p);
}

// What refine_points and refine_strip start from: phase 1's output.
struct Begin {
    IndexedMesh2 mesh;
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    Edges edges;
    std::vector<std::uint32_t> masks;
};

Begin begin_from(RefineOutcome r) {
    REQUIRE(r.ok());
    const std::size_t nt = r.triangles.size();
    return Begin{IndexedMesh2{std::move(r.vertices), std::move(r.triangles), std::vector<std::uint8_t>(nt, 0)},
                 std::move(r.z), std::move(r.valid), std::move(r.edges), std::move(r.masks)};
}

ConstraintCheckPoints strip_of(const Raster<float>& dem, const Begin& b) {
    return constraint_check_points(dem, b.mesh.vertices(), Span2{b.edges});
}

PointRefineOutcome strip_run(const Raster<float>& dem, const ConstraintCheckPoints& strip, const Begin& b,
                             const PointRefineOptions& o) {
    return refine_strip(dem, strip, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                        Span2{b.edges}, SpanM{b.masks}, o);
}

template <TolerancePolicy P>
PointRefineOutcome strip_run(const Raster<float>& dem, const ConstraintCheckPoints& strip, const Begin& b,
                             const PointRefineOptions& o, const P& p) {
    return refine_strip(dem, strip, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                        Span2{b.edges}, SpanM{b.masks}, o, p);
}

PointRefineOutcome points_run(const CheckPoints& store, const Begin& b, const PointRefineOptions& o) {
    return refine_points(store, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                         Span2{b.edges}, SpanM{b.masks}, o);
}

template <TolerancePolicy P>
PointRefineOutcome points_run(const CheckPoints& store, const Begin& b, const PointRefineOptions& o, const P& p) {
    return refine_points(store, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                         Span2{b.edges}, SpanM{b.masks}, o, p);
}

// ---------------------------------------------------------------- inputs

struct Sources {
    std::vector<Point2> xy;
    std::vector<float> z;
};

// Two points per cell at seeded dyadic offsets in (0, 1), z the DEM's bilinear
// value plus up to +-4 m of noise in 1/64 m steps, so it is a float exactly
// (ES9's scattered()).
Sources scattered(const Raster<float>& dem, std::uint32_t seed) {
    const auto& g = dem.geometry();
    std::mt19937 gen{seed};
    Sources p;
    for (std::size_t r = 0; r + 1 < g.rows(); ++r)
        for (std::size_t c = 0; c + 1 < g.cols(); ++c)
            for (int k = 0; k < 2; ++k) {
                const double col = static_cast<double>(c) + static_cast<double>(1u + gen() % 1023u) / 1024.0;
                const double row = static_cast<double>(r) + static_cast<double>(1u + gen() % 1023u) / 1024.0;
                const double noise = (static_cast<double>(gen() % 513u) - 256.0) / 64.0;
                const auto h = height(dem, Lat{col, row});
                if (!h) continue;
                p.xy.push_back(world(g, col, row));
                p.z.push_back(static_cast<float>(std::round((*h + noise) * 64.0) / 64.0));
            }
    return p;
}

CheckPoints store_of(const RasterGeometry& g, const Sources& src) {
    CheckPoints store{g};
    store.add(src.xy, src.z);
    store.freeze();
    return store;
}

// The jittered start's vertex (i, j) of its 5 x 5 grid.
std::uint32_t at(std::uint32_t i, std::uint32_t j) { return i * 5 + j; }

// A line crossing the whole domain, corner region to corner region.
Seg crossing() { return Seg{-5.0, -3.0, 70.0, -29.0}; }

// The three lines of test 5: one crossing the domain, one along the start's
// mask-2 chain (grid row 16, the chain's own edges), one from outside ending
// inside a start triangle.
std::vector<Seg> three_lines(const Start& s) {
    std::vector<Seg> out{crossing()};
    for (std::uint32_t j = 0; j + 1 < 5; ++j) {
        const Point2 a = s.mesh.vertices()[at(2, j)], b = s.mesh.vertices()[at(2, j + 1)];
        out.push_back(Seg{a.x, a.y, b.x, b.y});
    }
    out.push_back(Seg{75.0, -20.0, 40.3, -11.7});
    return out;
}

LineTolerance field(const RasterGeometry& g, const std::vector<Seg>& segs, const ToleranceRamp& r) {
    std::string why;
    auto f = LineTolerance::make(g, std::span<const Seg>{segs}, r, why);
    INFO(why);
    REQUIRE(f.has_value());
    return std::move(*f);
}

Ramp oracle_ramp(const ToleranceRamp& r) { return Ramp{r.near, r.far, r.start, r.end, 0.0}; }

const char* name(Terrain k) { return k == Terrain::Smooth ? "Smooth" : "SmoothWithHole"; }

// ---------------------------------------------------------------- equality

std::uint64_t bits(double v) { return std::bit_cast<std::uint64_t>(v); }

// Everything refine returned before increment 33, bit for bit; the timings
// are left out (they are not part of the determinism guarantee).
void same(const RefineOutcome& a, const RefineOutcome& b) {
    REQUIRE(a.status == b.status);
    REQUIRE(a.message == b.message);
    REQUIRE(a.vertices.size() == b.vertices.size());
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        CAPTURE(i);
        REQUIRE(bits(a.vertices[i].x) == bits(b.vertices[i].x));
        REQUIRE(bits(a.vertices[i].y) == bits(b.vertices[i].y));
        REQUIRE(bits(a.z[i]) == bits(b.z[i]));
    }
    REQUIRE(a.valid == b.valid);
    REQUIRE(a.triangles == b.triangles);
    REQUIRE(a.edges == b.edges);
    REQUIRE(a.masks == b.masks);
    CHECK(a.rounds == b.rounds);
    CHECK(a.inserted == b.inserted);
    CHECK(a.flips == b.flips);
    CHECK(bits(a.max_error) == bits(b.max_error));
    CHECK(a.uncovered == b.uncovered);
    CHECK(a.carved == b.carved);
    CHECK(a.quality_inserted == b.quality_inserted);
    CHECK(a.quality_skipped == b.quality_skipped);
    CHECK(a.quality_feet == b.quality_feet);
    CHECK(a.quality_no_gain == b.quality_no_gain);
    CHECK(a.quality_line_splits == b.quality_line_splits);
    CHECK(a.feet == b.feet);
    CHECK(a.feet_refused == b.feet_refused);
}

void same(const PointRefineOutcome& a, const PointRefineOutcome& b) {
    same(static_cast<const RefineOutcome&>(a), static_cast<const RefineOutcome&>(b));
    CHECK(a.coincident == b.coincident);
    CHECK(bits(a.coincident_max_error) == bits(b.coincident_max_error));
    CHECK(a.strip_points == b.strip_points);
    CHECK(a.strip_inserted == b.strip_inserted);
    CHECK(bits(a.strip_max_error) == bits(b.strip_max_error));
    CHECK(a.strip_refused == b.strip_refused);
    CHECK(bits(a.strip_refused_max_error) == bits(b.strip_refused_max_error));
    CHECK(a.nodes_inserted == b.nodes_inserted);
    CHECK(a.on_frozen == b.on_frozen);
    CHECK(bits(a.on_frozen_max_error) == bits(b.on_frozen_max_error));
    CHECK(a.feet_fallback == b.feet_fallback);
}

// The two section-D oracles on one output with a field.
template <class Outcome>
void oracles(const Raster<float>& dem, const Outcome& out, const std::vector<Seg>& segs, const Ramp& r) {
    REQUIRE(out.ok());
    const Mesh m = mesh_of(out);
    const auto f = line_tolerance_oracle::ramp_findings(dem, m, std::span<const Seg>{segs}, r);
    CAPTURE(f.nodes, f.worst_excess);
    CHECK(f.not_ccw == 0);
    CHECK(f.over == 0);
    CHECK(f.nodes > 0);
    CHECK(delaunay_violations(dem.geometry(), m) == 0);
}

// A test-only policy (test 7): the field's own at(), with bounds that never
// let the scan skip a query.
struct Unbounded {
    const LineTolerance* f;
    [[nodiscard]] double lowest() const noexcept { return 0.0; }
    [[nodiscard]] double highest() const noexcept { return std::numeric_limits<double>::infinity(); }
    [[nodiscard]] double at(const terrain::mesh::LatticeMesh& m, std::uint32_t t) const { return f->at(m, t); }
};
static_assert(TolerancePolicy<Unbounded>);

// One input of tests 1 to 3: a DEM, the start, and the sources and phase-1
// begin that refine_points and refine_strip need.
struct Scene {
    Raster<float> dem;
    Start start;
};

Scene scene(Terrain kind, std::uint32_t seed) {
    const auto g = exact_geometry();
    return Scene{terrain_dem(g, kind, seed), jittered(g, seed)};
}

}  // namespace

// ------------------------------------------------------------------- test 1

TEST_CASE("33 test 1: refine with UniformTolerance{t} is today's refine, bit for bit",
          "[line_tolerance][property][test1]") {
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const double tol = GENERATE(0.0, 0.5, 3.0);
    const double min_angle = GENERATE(0.0, 20.0);
    CAPTURE(name(kind), tol, min_angle);
    const auto sc = scene(kind, 3);
    const auto o = options(tol, true, 1, min_angle);
    same(run(sc.dem, sc.start, o), run(sc.dem, sc.start, o, UniformTolerance{tol}));
}

TEST_CASE("33 test 1: refine_strip and refine_points with UniformTolerance{t} are today's, bit for bit",
          "[line_tolerance][property][test1]") {
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const double tol = GENERATE(0.0, 0.5, 3.0);
    CAPTURE(name(kind), tol);
    const auto sc = scene(kind, 3);
    const Begin b = begin_from(run(sc.dem, sc.start, options(tol)));
    const auto strip = strip_of(sc.dem, b);
    const auto o = point_options(tol);
    same(strip_run(sc.dem, strip, b, o), strip_run(sc.dem, strip, b, o, UniformTolerance{tol}));
    const auto store = store_of(sc.dem.geometry(), scattered(sc.dem, 3));
    same(points_run(store, b, o), points_run(store, b, o, UniformTolerance{tol}));
}

// ------------------------------------------------------------------- test 2

TEST_CASE("33 test 2: LineTolerance with zero segments equals UniformTolerance{F} on all three entry points",
          "[line_tolerance][property][test2]") {
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const double far = GENERATE(0.5, 3.0);
    CAPTURE(name(kind), far);
    const auto sc = scene(kind, 5);
    const auto& g = sc.dem.geometry();
    const auto none = field(g, {}, ToleranceRamp{0.1, far, 0.0, 50.0, 1.0});
    const auto o = options(far);
    same(run(sc.dem, sc.start, o, none), run(sc.dem, sc.start, o, UniformTolerance{far}));
    const Begin b = begin_from(run(sc.dem, sc.start, o));
    const auto strip = strip_of(sc.dem, b);
    const auto po = point_options(far);
    same(strip_run(sc.dem, strip, b, po, none), strip_run(sc.dem, strip, b, po, UniformTolerance{far}));
    const auto store = store_of(g, scattered(sc.dem, 5));
    same(points_run(store, b, po, none), points_run(store, b, po, UniformTolerance{far}));
}

// ------------------------------------------------------------------- test 3

TEST_CASE("33 test 3 (a): LineTolerance with N = F equals UniformTolerance{F} on all three entry points",
          "[line_tolerance][property][test3]") {
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const double far = GENERATE(0.5, 3.0);
    CAPTURE(name(kind), far);
    const auto sc = scene(kind, 3);
    const auto& g = sc.dem.geometry();
    const auto f = field(g, three_lines(sc.start), ToleranceRamp{far, far, 0.0, 20.0, 0.5});
    const auto o = options(far);
    same(run(sc.dem, sc.start, o, f), run(sc.dem, sc.start, o, UniformTolerance{far}));
    const Begin b = begin_from(run(sc.dem, sc.start, o));
    const auto strip = strip_of(sc.dem, b);
    const auto po = point_options(far);
    same(strip_run(sc.dem, strip, b, po, f), strip_run(sc.dem, strip, b, po, UniformTolerance{far}));
    const auto store = store_of(g, scattered(sc.dem, 3));
    same(points_run(store, b, po, f), points_run(store, b, po, UniformTolerance{far}));
}

TEST_CASE("33 test 3 (b): a ramp holding every triangle at N equals UniformTolerance{N}, options at F",
          "[line_tolerance][property][test3]") {
    // M6 (the foot's epsilon from options.tolerance) and M7 (refine_points
    // comparing with options.tolerance): the LineTolerance run's options say
    // F, so only the field can bring N to the foot's epsilon and to the
    // comparisons. S = 100 m is above the domain's diagonal (71.6 m), so every
    // triangle's distance to the crossing line is below S.
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const double near = GENERATE(0.25, 0.5);
    constexpr double far = 6.0;
    CAPTURE(name(kind), near, far);
    const auto sc = scene(kind, 3);
    const auto& g = sc.dem.geometry();
    const std::vector<Seg> line{crossing()};
    const ToleranceRamp ramp{near, far, 100.0, 200.0, 1.0};
    const auto f = field(g, line, ramp);

    const auto line_run = run(sc.dem, sc.start, options(far), f);
    const auto uniform_run = run(sc.dem, sc.start, options(near), UniformTolerance{near});
    same(line_run, uniform_run);
    CHECK(line_run.feet + line_run.quality_feet > 0);  // the feet fire, so M6 can show
    CHECK(line_run.max_error_near == line_run.max_error);  // every triangle's allowed error is N (4.5)
    oracles(sc.dem, line_run, line, oracle_ramp(ramp));

    // Phase 1 with the field, as _dem_mesh runs it: the strip's DEM rescan
    // reaches only the triangles it writes, so the rest must already hold N.
    const Begin b = begin_from(line_run);
    const auto strip = strip_of(sc.dem, b);
    const auto line_strip = strip_run(sc.dem, strip, b, point_options(far), f);
    same(line_strip, strip_run(sc.dem, strip, b, point_options(near), UniformTolerance{near}));
    CHECK(line_strip.inserted > 0);
    oracles(sc.dem, line_strip, line, oracle_ramp(ramp));

    const auto src = scattered(sc.dem, 3);
    const auto store = store_of(g, src);
    const auto line_points = points_run(store, b, point_options(far), f);
    same(line_points, points_run(store, b, point_options(near), UniformTolerance{near}));
    CHECK(line_points.inserted > 0);
    CHECK(line_points.max_error <= near);
    CHECK(delaunay_violations(g, mesh_of(line_points)) == 0);
}

// ------------------------------------------------------------------- test 5

TEST_CASE("33 test 5: every valid DEM node is within t of its own distance to the original lines",
          "[line_tolerance][property][test5]") {
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const bool feet = GENERATE(false, true);
    const double min_angle = GENERATE(0.0, 20.0);
    const double margin = GENERATE(0.0, 0.5);
    const unsigned threads = GENERATE(1u, 4u);
    CAPTURE(name(kind), feet, min_angle, margin, threads);
    const auto sc = scene(kind, 7);
    const auto lines = three_lines(sc.start);
    // 5 cm on the lines to 3 m at 20 m: the domain spans 0 to about 30 m from them.
    const ToleranceRamp ramp{0.05, 3.0, 0.0, 20.0, margin};
    const auto f = field(sc.dem.geometry(), lines, ramp);
    const auto out = run(sc.dem, sc.start, options(3.0, feet, threads, min_angle), f);
    oracles(sc.dem, out, lines, oracle_ramp(ramp));
    CHECK(out.max_error <= 3.0);
    CHECK(out.max_error_near <= 0.05);

    // The oracle can fail here: the uniform-F mesh of the same input breaks it.
    const auto uniform = run(sc.dem, sc.start, options(3.0, feet, threads, min_angle));
    CHECK(line_tolerance_oracle::ramp_findings(sc.dem, mesh_of(uniform), std::span<const Seg>{lines},
                                               oracle_ramp(ramp))
              .over > 0);
}

// ------------------------------------------------------------------- test 7

TEST_CASE("33 test 7: a policy that is queried for every triangle gives LineTolerance's output",
          "[line_tolerance][property][test7]") {
    // M5: with lowest() = 0 and highest() = inf no triangle is decided by the
    // bounds, so every comparison goes through at().
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const bool feet = GENERATE(false, true);
    CAPTURE(name(kind), feet);
    const auto sc = scene(kind, 11);
    const auto& g = sc.dem.geometry();
    const auto lines = three_lines(sc.start);
    const ToleranceRamp ramp{0.1, 4.0, 2.0, 25.0, 0.5};
    const auto f = field(g, lines, ramp);
    const Unbounded every{&f};

    const auto lazy = run(sc.dem, sc.start, options(4.0, feet), f);
    same(lazy, run(sc.dem, sc.start, options(4.0, feet), every));
    oracles(sc.dem, lazy, lines, oracle_ramp(ramp));

    const Begin b = begin_from(lazy);  // phase 1 with the field, as in test 3 (b)
    const auto strip = strip_of(sc.dem, b);
    const auto lazy_strip = strip_run(sc.dem, strip, b, point_options(4.0, feet), f);
    same(lazy_strip, strip_run(sc.dem, strip, b, point_options(4.0, feet), every));
    oracles(sc.dem, lazy_strip, lines, oracle_ramp(ramp));

    const auto store = store_of(g, scattered(sc.dem, 11));
    same(points_run(store, b, point_options(4.0, feet), f), points_run(store, b, point_options(4.0, feet), every));
}

// ------------------------------------------------------------------- test 8

TEST_CASE("33 test 8: strip points within S of a line along a constraint end within N",
          "[line_tolerance][property][test8]") {
    // The line is the start's mask-2 chain itself (grid row 16), so every
    // strip point on it is at distance 0; points of the crossing chains within
    // S = 4 m of it are held to N too (their triangles are no farther).
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const bool feet = GENERATE(false, true);
    CAPTURE(name(kind), feet);
    const auto sc = scene(kind, 13);
    const auto& g = sc.dem.geometry();
    std::vector<Seg> line;
    for (std::uint32_t j = 0; j + 1 < 5; ++j) {
        const Point2 a = sc.start.mesh.vertices()[at(2, j)], b = sc.start.mesh.vertices()[at(2, j + 1)];
        line.push_back(Seg{a.x, a.y, b.x, b.y});
    }
    constexpr double near = 0.1, far = 8.0;
    const ToleranceRamp ramp{near, far, 4.0, 30.0, 0.0};
    const auto f = field(g, line, ramp);

    // As _dem_mesh runs it: refine with the field, then the strip with it.
    const Begin b = begin_from(run(sc.dem, sc.start, options(far, feet), f));
    const auto strip = strip_of(sc.dem, b);
    const auto out = strip_run(sc.dem, strip, b, point_options(far, feet), f);
    REQUIRE(out.ok());
    const Mesh m = mesh_of(out);

    std::vector<OraclePoint> close;
    for (const OraclePoint& p : ruled_points(sc.dem, b.mesh.vertices(), b.edges))
        if (line_tolerance_oracle::point_lines(world(g, p.at.col, p.at.row), std::span<const Seg>{line}) <= 4.0)
            close.push_back(p);
    REQUIRE(close.size() > 20);
    const auto s = strip_findings(g, close, m, near);
    CAPTURE(s.checked, s.worst);
    CHECK(s.unfiled == 0);
    CHECK(s.over == 0);
    oracles(sc.dem, out, line, oracle_ramp(ramp));

    // The strip check can fail here: the uniform-F strip run leaves points over N.
    const Begin u = begin_from(run(sc.dem, sc.start, options(far, feet)));
    const auto ustrip = strip_of(sc.dem, u);
    const auto uniform = strip_run(sc.dem, ustrip, u, point_options(far, feet));
    std::vector<OraclePoint> uclose;
    for (const OraclePoint& p : ruled_points(sc.dem, u.mesh.vertices(), u.edges))
        if (line_tolerance_oracle::point_lines(world(g, p.at.col, p.at.row), std::span<const Seg>{line}) <= 4.0)
            uclose.push_back(p);
    CHECK(strip_findings(g, uclose, mesh_of(uniform), near).over > 0);
}

// ------------------------------------------------------------------- test 9

TEST_CASE("33 test 9: 1 and 8 threads give the same mesh with a field", "[line_tolerance][property][test9]") {
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    CAPTURE(name(kind));
    const auto sc = scene(kind, 17);
    const auto& g = sc.dem.geometry();
    const auto lines = three_lines(sc.start);
    const ToleranceRamp ramp{0.05, 3.0, 0.0, 20.0, 0.5};
    const auto f = field(g, lines, ramp);
    const auto one = run(sc.dem, sc.start, options(3.0, true, 1, 20.0), f);
    const auto eight = run(sc.dem, sc.start, options(3.0, true, 8, 20.0), f);
    same(one, eight);
    CHECK(bits(one.max_error_near) == bits(eight.max_error_near));
    oracles(sc.dem, eight, lines, oracle_ramp(ramp));

    const Begin b = begin_from(run(sc.dem, sc.start, options(3.0)));
    const auto strip = strip_of(sc.dem, b);
    same(strip_run(sc.dem, strip, b, point_options(3.0, true, 1), f),
         strip_run(sc.dem, strip, b, point_options(3.0, true, 8), f));
    const auto store = store_of(g, scattered(sc.dem, 17));
    same(points_run(store, b, point_options(3.0, true, 1), f), points_run(store, b, point_options(3.0, true, 8), f));
}

// ------------------------------------------------- max_error_near, many blocks

namespace {

// A policy whose at() depends on the slot alone: N = 1 on every third slot,
// otherwise 1.5 to 4.5, so lowest() 1 and highest() 4.5 (the mesh is unread).
struct BySlot {
    [[nodiscard]] double lowest() const noexcept { return 1.0; }
    [[nodiscard]] double highest() const noexcept { return 4.5; }
    [[nodiscard]] double at(const terrain::mesh::LatticeMesh&, std::uint32_t t) const noexcept {
        return t % 3 == 0 ? 1.0 : 1.0 + 0.5 * static_cast<double>(1 + t % 7);
    }
};
static_assert(TolerancePolicy<BySlot>);

// No slot is held to lowest(): every at() is 2.
struct FarSlots {
    [[nodiscard]] double lowest() const noexcept { return 1.0; }
    [[nodiscard]] double highest() const noexcept { return 4.5; }
    [[nodiscard]] double at(const terrain::mesh::LatticeMesh&, std::uint32_t) const noexcept { return 2.0; }
};

// Section 9.1's definition, serially: the largest error over the slots whose
// at() is at most lowest(); 0.0 when there are none.
double brute_near(const BySlot& p, const terrain::mesh::LatticeMesh& m, const std::vector<double>& error) {
    double best = 0.0;
    for (std::uint32_t t = 0; t < error.size(); ++t)
        if (p.at(m, t) <= p.lowest()) best = std::max(best, error[t]);
    return best;
}

}  // namespace

TEST_CASE("33: max_error_near over many blocks is the serial definition, for 1 and 8 threads",
          "[line_tolerance][property][max_error_near]") {
    // Code review round 2: detail::max_error_near works in blocks of 4 096
    // slots, and for_each_block runs a single block inline, so test 9's
    // 33 x 33 grid never starts a second worker. 3 x 4 096 + 17 slots here.
    const BySlot p;
    auto built = terrain::mesh::LatticeMesh::build(
        std::vector<terrain::mesh::MeshVertex>{{0.0, 0.0}, {1.0, 1.0}, {1.0, 0.0}}, {{0, 1, 2}},
        std::vector<std::uint8_t>(1, 0), std::vector<std::array<std::uint32_t, 3>>(1, {0, 0, 0}));
    REQUIRE(built.has_value());
    const terrain::mesh::LatticeMesh& m = *built;

    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    CAPTURE(seed);
    constexpr std::size_t n = 3 * 4096 + 17;
    std::mt19937 gen{seed};
    std::vector<double> error(n), allowed(n);
    std::size_t at_lowest = 0, splits = 0, asked = 0;
    for (std::uint32_t t = 0; t < n; ++t) {
        error[t] = static_cast<double>(gen() % 6001u) / 1000.0;  // 0 to 6 m in mm steps
        allowed[t] = terrain::refinement::detail::allowed_at(p, m, t, error[t]);
        at_lowest += allowed[t] == p.lowest() ? 1 : 0;
        splits += allowed[t] == -std::numeric_limits<double>::infinity() ? 1 : 0;
        asked += allowed[t] > p.lowest() ? 1 : 0;
    }
    REQUIRE(at_lowest > 0);  // all three of allowed_at's values occur
    REQUIRE(splits > 0);
    REQUIRE(asked > 0);
    // The largest near error sits in the last, partial block, so a block
    // writing to another block's slot of the partial maxima shows.
    error[3 * 4096 + 6] = 7.0;  // 12 294 % 3 == 0: a near slot
    allowed[3 * 4096 + 6] = -std::numeric_limits<double>::infinity();
    // And a far slot above it everywhere else, which must not count.
    error[4097] = 9.0;  // 4 097 % 3 == 2: far
    allowed[4097] = -std::numeric_limits<double>::infinity();

    const auto err = [&](std::uint32_t t) { return error[t]; };
    const double expected = brute_near(p, m, error);
    REQUIRE(expected == 7.0);
    const double one = terrain::refinement::detail::max_error_near(p, m, allowed, 1, err);
    const double eight = terrain::refinement::detail::max_error_near(p, m, allowed, 8, err);
    CHECK(bits(one) == bits(expected));
    CHECK(bits(eight) == bits(expected));

    // Per block: each block's own near maximum is found (the last block's
    // 7.0 removed, the result is the brute force over the rest).
    error[3 * 4096 + 6] = 0.0;
    const double rest = brute_near(p, m, error);
    CHECK(bits(terrain::refinement::detail::max_error_near(p, m, allowed, 8, err)) == bits(rest));
    CHECK(rest < 7.0);

    // No slot is near: 0.0, whatever the errors.
    const FarSlots far;
    std::vector<double> far_allowed(n);
    for (std::uint32_t t = 0; t < n; ++t) far_allowed[t] = terrain::refinement::detail::allowed_at(far, m, t, error[t]);
    CHECK(terrain::refinement::detail::max_error_near(far, m, far_allowed, 8, err) == 0.0);
}
