// Increment 34 (docs/increments/34-slope-tolerance.md, sections 3, 4.3, 4.4,
// 7 and 9): the slope tolerance inside the three refinement entry points.
//
//   test 1  on V1c (constraints and feet on), refine(..., options) and
//           refine(..., options, UniformTolerance{t}) give the same output,
//           bit for bit, on all three entry points (G1).
//   test 4  equalities, bit for bit, on all three entry points, on V1c:
//           (a) Sloped<UniformTolerance{10}> with N = F = 10 equals
//               UniformTolerance{10} (G2);
//               also on a 5 x 5 fixture whose two inner errors are adjacent
//               doubles with one product, so only the second key picks
//               today's node (M2);
//           (b) START = END = 0 with N = 2, F = 10 equals UniformTolerance{2}
//               (G6; the slope must reach the foot's epsilon, M6, and
//               refine_points' comparison, M7);
//           (c) Sloped<LineTolerance> with N = F equals LineTolerance (G8).
//   test 5  the guarantee: every valid DEM node within the ramp at the
//           test's own class, on V1 and V1c; with V1's line, within the
//           smaller of that and 33's ramp at its brute-force distance (G3).
//   test 6  (a) check points (4 per cell over V1) and the edge strip on V1c:
//               every point within t of its cell's largest valid corner
//               class (G4, M3, M7);
//           (b) V1n: the 16 check points in the four cells around the NoData
//               node end within N = 2 (M9);
//           (c) V1n: RefineOutcome::slope_nodes reports 4 224 valid nodes and
//               1 495 tightened; the NoData node is in neither.
//   test 7  1 and 8 threads give the same mesh, with the slope alone and with
//           lines (G7). TSan: this suite starts threads.
//   test 8  exactness of `over` (M1): an error one unit in the last place
//           above allowed(class), chosen so that error * weight rounds to
//           exactly 1, still splits.
//   M8      the lazy test never asks the triangle part about a triangle whose
//           nodes already split it (a counting policy).
//
// Invariant-critical (section 9): tests 4, 5, 6 and 8 here.
//
// Interface used, as sections 4.3 and 4.5 write it: SlopeTolerance::make(dem,
// SlopeRamp{near, far, start_deg, end_deg}, threads, why); the policy
// Sloped<P>{triangle, &slope} in the entry points' policy slot (33's
// argument order); RefineOutcome::max_slope_share. CHOSEN HERE where 4.5
// says only "RefineOutcome::slope_nodes (two counts)": members `valid` and
// `tightened`.
//
// Oracles (QA rules, section D): every output with a slope carries the
// constrained-Delaunay oracle (strip_oracle::delaunay_violations) and the
// slope tolerance oracle (slope_oracle::node_findings with each node's own
// allowed error from the test's Horn); the controls show the slope oracles
// fail on the uniform-F mesh of the same input.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/line_tolerance.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>
#include <terrain/refinement/slope_tolerance.hpp>

#include "line_tolerance_oracle.hpp"
#include "slope_oracle.hpp"
#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <string>
#include <utility>
#include <vector>

using namespace strip_oracle;
using slope_oracle::Sources;
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
using terrain::refinement::SlopeRamp;
using terrain::refinement::Sloped;
using terrain::refinement::SlopeTolerance;
using terrain::refinement::ToleranceRamp;
using terrain::refinement::TolerancePolicy;
using terrain::refinement::UniformTolerance;

namespace {

// ---------------------------------------------------------------- the runs

using Span2 = std::span<const std::array<std::uint32_t, 2>>;
using SpanM = std::span<const std::uint32_t>;

RefineOptions options(double tol, bool feet = true, unsigned threads = 1) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.constraint_feet = feet;
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

CheckPoints store_of(const RasterGeometry& g, const Sources& src) {
    CheckPoints store{g};
    store.add(src.xy, src.z);
    store.freeze();
    return store;
}

// ---------------------------------------------------------------- equality

std::uint64_t bits(double v) { return std::bit_cast<std::uint64_t>(v); }

// Everything refine returned before increment 34, bit for bit; the timings
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
    CHECK(a.nodes_inserted == b.nodes_inserted);
    CHECK(a.on_frozen == b.on_frozen);
    CHECK(a.feet_fallback == b.feet_fallback);
}

}  // namespace

// ------------------------------------------------------------------- test 1

TEST_CASE("34 test 1: on V1c, refine with UniformTolerance{t} is today's refine on all three entry points",
          "[slope][property][test1]") {
    const double tol = GENERATE(2.0, 10.0);
    CAPTURE(tol);
    const auto dem = slope_oracle::v1();
    const auto start = slope_oracle::v1c_start();
    const auto o = options(tol);
    const auto plain = run(dem, start, o);
    same(plain, run(dem, start, o, UniformTolerance{tol}));
    CHECK(plain.feet + plain.quality_feet > 0);  // the feet act on the cut edge
    const Begin b = begin_from(plain);
    const auto strip = strip_of(dem, b);
    const auto po = point_options(tol);
    same(strip_run(dem, strip, b, po), strip_run(dem, strip, b, po, UniformTolerance{tol}));
    const auto store = store_of(dem.geometry(), slope_oracle::scattered(dem, 7));
    same(points_run(store, b, po), points_run(store, b, po, UniformTolerance{tol}));
}

namespace {

SlopeTolerance slope_of(const Raster<float>& dem, SlopeRamp r, unsigned threads = 1) {
    std::string why;
    auto s = SlopeTolerance::make(dem, r, threads, why);
    INFO(why);
    REQUIRE(s.has_value());
    return std::move(*s);
}

LineTolerance field(const RasterGeometry& g, const std::vector<line_tolerance_oracle::Seg>& segs,
                    const ToleranceRamp& r) {
    std::string why;
    auto f = LineTolerance::make(g, std::span<const line_tolerance_oracle::Seg>{segs}, r, why);
    INFO(why);
    REQUIRE(f.has_value());
    return std::move(*f);
}

// V1's line: N = 1 on it, F = 10 from 200 m (section 9, *Fixtures*).
constexpr ToleranceRamp kLineRamp{1.0, 10.0, 0.0, 200.0, 0.0};
line_tolerance_oracle::Ramp line_oracle() { return {1.0, 10.0, 0.0, 200.0, 0.0}; }

slope_oracle::Ramp oracle(const SlopeRamp& r) { return {r.near, r.far, r.start_deg, r.end_deg}; }

// The QA rules' section D oracles on one output, each node held to allowed[n].
template <class Outcome>
void oracles(const Raster<float>& dem, const Outcome& out, const std::vector<double>& allowed) {
    REQUIRE(out.ok());
    const Mesh m = mesh_of(out);
    const auto f = slope_oracle::node_findings(dem, m, std::span<const double>{allowed});
    CAPTURE(f.nodes, f.worst_excess);
    CHECK(f.not_ccw == 0);
    CHECK(f.over == 0);
    CHECK(f.nodes > 0);
    CHECK(delaunay_violations(dem.geometry(), m) == 0);
}

// Two policies with their options; equal_on_three compares their runs on all three entry points.
template <class A, class B>
struct Pair {
    const A& a;
    RefineOptions oa;
    PointRefineOptions pa;
    const B& b;
    RefineOptions ob;
    PointRefineOptions pb;
};

template <class A, class B>
void equal_on_three(const Raster<float>& dem, const Start& start, const Pair<A, B>& p) {
    const auto first = run(dem, start, p.oa, p.a);
    same(first, run(dem, start, p.ob, p.b));
    const Begin b = begin_from(first);
    const auto strip = strip_of(dem, b);
    same(strip_run(dem, strip, b, p.pa, p.a), strip_run(dem, strip, b, p.pb, p.b));
    const auto store = store_of(dem.geometry(), slope_oracle::scattered(slope_oracle::v1(), 7));
    same(points_run(store, b, p.pa, p.a), points_run(store, b, p.pb, p.b));
}

}  // namespace

// ------------------------------------------------------------------- test 4

TEST_CASE("34 test 4 (a): Sloped<UniformTolerance{10}> with N = F = 10 is UniformTolerance{10}, bit for bit",
          "[slope][property][test4]") {
    const auto dem = slope_oracle::v1();
    const auto start = slope_oracle::v1c_start();
    const double end = GENERATE(0.0, 30.0);
    CAPTURE(end);
    const auto s = slope_of(dem, SlopeRamp{10.0, 10.0, 0.0, end});
    const Sloped<UniformTolerance> sloped{UniformTolerance{10.0}, &s};
    const UniformTolerance uniform{10.0};
    equal_on_three(dem, start,
                   Pair<Sloped<UniformTolerance>, UniformTolerance>{sloped, options(10.0), point_options(10.0),
                                                                    uniform, options(10.0), point_options(10.0)});
}

TEST_CASE("34 test 4 (b): START = END = 0 with N = 2, F = 10 is UniformTolerance{2}, bit for bit, feet included",
          "[slope][property][test4]") {
    // M6 (the foot's epsilon from the triangle part's 10) and M7 (point_loop
    // without the slope): the Sloped run's triangle part and options say 10,
    // so only the slope can bring 2 to the foot's epsilon and to
    // refine_points' comparison. V1c's cut edge crosses the wall, where three
    // DEM-node vertices of the --tolerance 2 mesh lie between eps(2) and
    // eps(10) of it (README item 4b).
    const auto dem = slope_oracle::v1();
    const auto start = slope_oracle::v1c_start();
    const auto s = slope_of(dem, SlopeRamp{2.0, 10.0, 0.0, 0.0});
    const Sloped<UniformTolerance> sloped{UniformTolerance{10.0}, &s};
    const UniformTolerance uniform{2.0};
    const auto uniform_run = run(dem, start, options(2.0), uniform);
    CHECK(uniform_run.feet > 0);  // the feet fire, so M6 can show
    equal_on_three(dem, start,
                   Pair<Sloped<UniformTolerance>, UniformTolerance>{sloped, options(10.0), point_options(10.0),
                                                                    uniform, options(2.0), point_options(2.0)});
    const auto sloped_run = run(dem, start, options(10.0), sloped);
    CHECK(sloped_run.max_slope_share <= 1.0 + 1e-12);
    CHECK(sloped_run.max_slope_share > 0.5);
}

TEST_CASE("34 test 4 (c): Sloped<LineTolerance> with N = F is LineTolerance alone, bit for bit",
          "[slope][property][test4]") {
    const auto dem = slope_oracle::v1();
    const auto start = slope_oracle::v1c_start();
    const auto lines = slope_oracle::v1_line();
    const auto f = field(dem.geometry(), lines, kLineRamp);
    const auto s = slope_of(dem, SlopeRamp{10.0, 10.0, 30.0, 30.0});
    const Sloped<LineTolerance> sloped{f, &s};
    equal_on_three(dem, start,
                   Pair<Sloped<LineTolerance>, LineTolerance>{sloped, options(10.0), point_options(10.0), f,
                                                              options(10.0), point_options(10.0)});
}

// ------------------------------------------------------------------- test 5

TEST_CASE("34 test 5: every valid DEM node is within the ramp at its own slope", "[slope][property][test5]") {
    const bool cut = GENERATE(false, true);
    const bool feet = GENERATE(false, true);
    const unsigned threads = GENERATE(1u, 4u);
    const auto ramp = GENERATE(SlopeRamp{2.0, 10.0, 30.0, 30.0}, SlopeRamp{2.0, 10.0, 25.0, 35.0});
    CAPTURE(cut, feet, threads, ramp.start_deg, ramp.end_deg);
    const auto dem = slope_oracle::v1();
    const auto start = cut ? slope_oracle::v1c_start() : slope_oracle::outline(dem.geometry());
    const auto s = slope_of(dem, ramp, threads);
    const auto allowed = slope_oracle::slope_allowed(slope_oracle::classes(dem), oracle(ramp));
    const auto out = run(dem, start, options(10.0, feet, threads), Sloped<UniformTolerance>{UniformTolerance{10.0}, &s});
    oracles(dem, out, allowed);
    CHECK(out.max_error <= 10.0);
    CHECK(out.max_slope_share <= 1.0 + 1e-12);

    // The oracle can fail here: the uniform-F mesh of the same input breaks it.
    const auto uniform = run(dem, start, options(10.0, feet, threads));
    CHECK(slope_oracle::node_findings(dem, mesh_of(uniform), std::span<const double>{allowed}).over > 0);
}

TEST_CASE("34 test 5: with V1's line, every node within the smaller of the two ramps", "[slope][property][test5]") {
    const bool cut = GENERATE(false, true);
    const bool feet = GENERATE(false, true);
    CAPTURE(cut, feet);
    const auto dem = slope_oracle::v1();
    const auto& g = dem.geometry();
    const auto start = cut ? slope_oracle::v1c_start() : slope_oracle::outline(g);
    const auto lines = slope_oracle::v1_line();
    const SlopeRamp ramp{2.0, 10.0, 30.0, 30.0};
    const auto s = slope_of(dem, ramp);
    const auto f = field(g, lines, kLineRamp);
    const auto allowed = slope_oracle::slope_and_lines_allowed(
        g, slope_oracle::classes(dem), oracle(ramp), std::span<const line_tolerance_oracle::Seg>{lines}, line_oracle());
    const auto out = run(dem, start, options(10.0, feet), Sloped<LineTolerance>{f, &s});
    oracles(dem, out, allowed);
    CHECK(out.max_slope_share <= 1.0 + 1e-12);

    // Both drivers decide somewhere: neither alone meets the combined oracle.
    const auto slope_alone = run(dem, start, options(10.0, feet), Sloped<UniformTolerance>{UniformTolerance{10.0}, &s});
    const auto lines_alone = run(dem, start, options(10.0, feet), f);
    CHECK(slope_oracle::node_findings(dem, mesh_of(slope_alone), std::span<const double>{allowed}).over > 0);
    CHECK(slope_oracle::node_findings(dem, mesh_of(lines_alone), std::span<const double>{allowed}).over > 0);
}

// ------------------------------------------------------------------- test 6

TEST_CASE("34 test 6 (a): every check point ends within t of its cell's largest valid corner class",
          "[slope][property][test6]") {
    // 33's resampled-path setup on V1: phase 1 with the slope, then the
    // final check with it. 4 points per cell (this file's draw, mt19937; not
    // NumPy's scattered(65, 4, 7), which test_core_slope_tolerance.py uses):
    // 237 points lie in cells straddling 30 degrees nearest a corner below
    // it, 125 of them with noise above 2 m, which M3 would hold to 10.
    const bool feet = GENERATE(false, true);
    CAPTURE(feet);
    const auto dem = slope_oracle::v1();
    const auto& g = dem.geometry();
    const auto start = slope_oracle::outline(g);
    const SlopeRamp ramp{2.0, 10.0, 30.0, 30.0};
    const auto s = slope_of(dem, ramp);
    const Sloped<UniformTolerance> sloped{UniformTolerance{10.0}, &s};
    const auto src = slope_oracle::scattered(dem, 7);
    const auto store = store_of(g, src);
    const auto allowed = slope_oracle::point_allowed(g, slope_oracle::classes(dem), src, oracle(ramp));

    const Begin b = begin_from(run(dem, start, options(10.0, feet), sloped));
    const auto out = points_run(store, b, point_options(10.0, feet), sloped);
    REQUIRE(out.ok());
    const auto f = slope_oracle::point_findings(g, src, std::span<const double>{allowed}, mesh_of(out));
    CAPTURE(f.checked, f.at_vertex, f.worst_excess);
    CHECK(f.checked + f.at_vertex == src.xy.size());
    CHECK(f.checked > src.xy.size() / 2);
    CHECK(f.unheld == 0);
    CHECK(f.over == 0);
    CHECK(out.max_slope_share <= 1.0 + 1e-12);
    CHECK(delaunay_violations(g, mesh_of(out)) == 0);

    // The point oracle can fail here: the uniform-F final check breaks it.
    const Begin u = begin_from(run(dem, start, options(10.0, feet)));
    const auto uniform = points_run(store, u, point_options(10.0, feet));
    CHECK(slope_oracle::point_findings(g, src, std::span<const double>{allowed}, mesh_of(uniform)).over > 0);
}

TEST_CASE("34 test 6 (a): every strip point along V1c's edges ends within its cell's t", "[slope][property][test6]") {
    const bool feet = GENERATE(false, true);
    CAPTURE(feet);
    const auto dem = slope_oracle::v1();
    const auto& g = dem.geometry();
    const auto start = slope_oracle::v1c_start();
    const SlopeRamp ramp{2.0, 10.0, 30.0, 30.0};
    const auto s = slope_of(dem, ramp);
    const Sloped<UniformTolerance> sloped{UniformTolerance{10.0}, &s};
    const auto cls = slope_oracle::classes(dem);

    // The ruled points grouped by their cell's allowed error.
    const auto groups = [&](const Begin& b) {
        std::map<double, std::vector<OraclePoint>> by;
        for (const OraclePoint& p : ruled_points(dem, b.mesh.vertices(), b.edges))
            by[oracle(ramp).of_class(slope_oracle::cell_class(g, cls, p.at.col, p.at.row))].push_back(p);
        return by;
    };

    const Begin b = begin_from(run(dem, start, options(10.0, feet), sloped));
    const auto out = strip_run(dem, strip_of(dem, b), b, point_options(10.0, feet), sloped);
    REQUIRE(out.ok());
    const auto by = groups(b);
    REQUIRE(by.contains(2.0));
    for (const auto& [t, pts] : by) {
        const auto f = strip_findings(g, pts, mesh_of(out), t);
        CAPTURE(t, f.checked, f.worst);
        CHECK(f.unfiled == 0);
        CHECK(f.over == 0);
    }
    oracles(dem, out, slope_oracle::slope_allowed(cls, oracle(ramp)));

    // The strip check can fail here: the uniform-F strip run leaves points over 2.
    const Begin u = begin_from(run(dem, start, options(10.0, feet)));
    const auto uniform = strip_run(dem, strip_of(dem, u), u, point_options(10.0, feet));
    CHECK(strip_findings(g, groups(u).at(2.0), mesh_of(uniform), 2.0).over > 0);
}

TEST_CASE("34 test 6 (b): on V1n the 16 check points around the NoData node end within N = 2",
          "[slope][property][test6]") {
    // Their cells each have the NoData node as a corner and valid corners at
    // classes 81 to 86 (40.5 to 43 degrees). z: V1's bilinear surface (the
    // NoData node's V1 value standing in) plus noise; 8 of the 16 carry
    // noise above 2 m in this file's draw. M9 would hold them to F = 10.
    const auto dem = slope_oracle::v1n();
    const auto& g = dem.geometry();
    const SlopeRamp ramp{2.0, 10.0, 30.0, 30.0};
    const auto s = slope_of(dem, ramp);
    const Sloped<UniformTolerance> sloped{UniformTolerance{10.0}, &s};
    const auto all = slope_oracle::scattered(slope_oracle::v1(), 7);
    Sources around;
    for (std::size_t k = 0; k < all.xy.size(); ++k) {
        const Lat l = lat(g, all.xy[k]);
        if ((l.row >= 31.0 && l.row < 33.0) && (l.col >= 29.0 && l.col < 31.0)) {
            around.xy.push_back(all.xy[k]);
            around.z.push_back(all.z[k]);
        }
    }
    REQUIRE(around.xy.size() == 16);
    const auto store = store_of(g, all);
    const Begin b = begin_from(run(dem, slope_oracle::outline(g), options(10.0), sloped));
    const auto out = points_run(store, b, point_options(10.0), sloped);
    REQUIRE(out.ok());
    const std::vector<double> two(around.xy.size(), 2.0);
    const auto f = slope_oracle::point_findings(g, around, std::span<const double>{two}, mesh_of(out));
    CAPTURE(f.worst_excess);
    CHECK(f.checked + f.at_vertex == 16);
    CHECK(f.unheld == 0);
    CHECK(f.over == 0);
    // Every point of the store, by the oracle's own classes.
    const auto allowed = slope_oracle::point_allowed(g, slope_oracle::classes(dem), all, oracle(ramp));
    CHECK(slope_oracle::point_findings(g, all, std::span<const double>{allowed}, mesh_of(out)).over == 0);
}

TEST_CASE("34 test 6 (c): slope_nodes counts V1n's valid nodes and the tightened ones, NoData in neither",
          "[slope][property][test6]") {
    // README item 4b: 4 224 valid nodes of 4 225, 1 495 held to 2 m by the
    // step at 30. Counted by refine over the final triangles (section 4.5).
    const unsigned threads = GENERATE(1u, 8u);
    CAPTURE(threads);
    const auto dem = slope_oracle::v1n();
    const auto s = slope_of(dem, SlopeRamp{2.0, 10.0, 30.0, 30.0});
    const auto out = run(dem, slope_oracle::outline(dem.geometry()), options(10.0, true, threads),
                         Sloped<UniformTolerance>{UniformTolerance{10.0}, &s});
    REQUIRE(out.ok());
    CHECK(out.slope_nodes.valid == 4224u);
    CHECK(out.slope_nodes.tightened == 1495u);
}

TEST_CASE("34 test 6 (c): slope_nodes counts only the nodes inside the mesh", "[slope][property][test6]") {
    // V1c's domain, a convex pentagon: the nodes in it (closed), counted here.
    const auto dem = slope_oracle::v1();
    const auto& g = dem.geometry();
    const auto ramp = GENERATE(SlopeRamp{2.0, 10.0, 30.0, 30.0}, SlopeRamp{2.0, 10.0, 25.0, 35.0});
    CAPTURE(ramp.start_deg, ramp.end_deg);
    const auto start = slope_oracle::v1c_start();
    const auto& ring = start.mesh.vertices();
    const auto cls = slope_oracle::classes(dem);
    std::size_t valid = 0, tightened = 0;
    for (std::size_t i = 0; i < g.size(); ++i) {
        const Point2 p = g.node({i / g.cols(), i % g.cols()});
        bool in = true;
        for (std::size_t k = 0; k < ring.size(); ++k)
            in = in && line_tolerance_oracle::orient(ring[k], ring[(k + 1) % ring.size()], p) >= 0.0;
        if (!in) continue;
        ++valid;
        tightened += oracle(ramp).of_class(cls[i]) < ramp.far ? 1 : 0;
    }
    REQUIRE(valid < g.size());
    const auto s = slope_of(dem, ramp);
    const auto out = run(dem, start, options(10.0), Sloped<UniformTolerance>{UniformTolerance{10.0}, &s});
    REQUIRE(out.ok());
    CHECK(out.slope_nodes.valid == valid);
    CHECK(out.slope_nodes.tightened == tightened);
}

TEST_CASE("34: without a slope, max_slope_share and slope_nodes are zero (G1)", "[slope][property][test1]") {
    const auto dem = slope_oracle::v1();
    const auto out = run(dem, slope_oracle::v1c_start(), options(10.0));
    REQUIRE(out.ok());
    CHECK(out.max_slope_share == 0.0);
    CHECK(out.slope_nodes.valid == 0u);
    CHECK(out.slope_nodes.tightened == 0u);
}

// ------------------------------------------------------------------- test 7

TEST_CASE("34 test 7: 1 and 8 threads give the same mesh, with the slope alone and with lines",
          "[slope][property][test7][threads]") {
    const bool with_lines = GENERATE(false, true);
    CAPTURE(with_lines);
    const auto dem = slope_oracle::v1();
    const auto& g = dem.geometry();
    const auto start = slope_oracle::v1c_start();
    const SlopeRamp ramp{2.0, 10.0, 25.0, 35.0};
    const auto s1 = slope_of(dem, ramp, 1);
    const auto s8 = slope_of(dem, ramp, 8);
    const auto lines = slope_oracle::v1_line();
    const auto f = field(g, lines, kLineRamp);
    const auto check = [&](const auto& p1, const auto& p8) {
        const auto one = run(dem, start, options(10.0, true, 1), p1);
        const auto eight = run(dem, start, options(10.0, true, 8), p8);
        same(one, eight);
        CHECK(bits(one.max_slope_share) == bits(eight.max_slope_share));
        CHECK(one.slope_nodes.valid == eight.slope_nodes.valid);
        CHECK(one.slope_nodes.tightened == eight.slope_nodes.tightened);
        const Begin b = begin_from(one);
        const auto strip = strip_of(dem, b);
        same(strip_run(dem, strip, b, point_options(10.0, true, 1), p1),
             strip_run(dem, strip, b, point_options(10.0, true, 8), p8));
        const auto store = store_of(g, slope_oracle::scattered(dem, 17));
        same(points_run(store, b, point_options(10.0, true, 1), p1),
             points_run(store, b, point_options(10.0, true, 8), p8));
    };
    if (with_lines)
        check(Sloped<LineTolerance>{f, &s1}, Sloped<LineTolerance>{f, &s8});
    else
        check(Sloped<UniformTolerance>{UniformTolerance{10.0}, &s1}, Sloped<UniformTolerance>{UniformTolerance{10.0}, &s8});
}

// ------------------------------------------------------------------- test 8

namespace {

// A 4 x 4 grid at 10 m, flat 0 but for node (1, 1) at z; one start triangle
// on nodes (0, 0), (3, 0), (0, 3), whose only node strictly inside is (1, 1).
struct OneNode {
    Raster<float> dem;
    Start start;
};

OneNode one_node(float z) {
    const auto g = slope_oracle::geometry(4, 4, 10.0, 10.0);
    std::vector<float> v(16, 0.0f);
    v[1 * 4 + 1] = z;
    return OneNode{Raster<float>{g, std::move(v)}, slope_oracle::fan({{0.0, 0.0}, {0.0, -30.0}, {30.0, 0.0}})};
}

}  // namespace

TEST_CASE("34 test 8: an error one ulp above allowed splits even when error * weight rounds to 1 (M1)",
          "[slope][property][test8]") {
    // N is the double just below a float f, so the node's error f is one unit
    // in the last place above N; f is the first float down from 4 for which
    // f * (1 / N) rounds to exactly 1 (3.99998069), so a test on the ratio
    // would converge. START = END = 0 holds every node to N; F and the
    // triangle part are 1000.
    float f = std::nextafter(4.0f, 0.0f);
    double near = 0.0;
    for (int k = 0; k < 1 << 20; ++k, f = std::nextafter(f, 0.0f)) {
        near = std::nextafter(static_cast<double>(f), 0.0);
        if (static_cast<double>(f) * (1.0 / near) == 1.0) break;
    }
    CAPTURE(f, near);
    REQUIRE(static_cast<double>(f) * (1.0 / near) == 1.0);
    REQUIRE(static_cast<double>(f) > near);
    REQUIRE(std::nextafter(near, 8.0) == static_cast<double>(f));

    const auto one = one_node(f);
    const auto s = slope_of(one.dem, SlopeRamp{near, 1000.0, 0.0, 0.0});
    REQUIRE(s.allowed(0) == near);
    const auto out = run(one.dem, one.start, options(1000.0, false), Sloped<UniformTolerance>{UniformTolerance{1000.0}, &s});
    REQUIRE(out.ok());
    CHECK(out.inserted == 1);
    CHECK(out.max_error == 0.0);

    // At exactly N the node is within its tolerance: nothing to insert.
    const auto exact = slope_of(one.dem, SlopeRamp{static_cast<double>(f), 1000.0, 0.0, 0.0});
    const auto none = run(one.dem, one.start, options(1000.0, false),
                          Sloped<UniformTolerance>{UniformTolerance{1000.0}, &exact});
    REQUIRE(none.ok());
    CHECK(none.inserted == 0);
}

namespace {

// M2's fixture: a 5 x 5 grid at 10 m, one start triangle on nodes (0, 0),
// (4, 0), (0, 4), z 0 but for its corner (0, 4) at -2^-50 and nodes (1, 1) and
// (1, 2) at 1.5. The plane falls by 2^-52 per column, so the two nodes'
// errors are 1.5 + 2^-52 and 1.5 + 2^-51: adjacent doubles, the smaller
// first in the row-major walk.
OneNode two_errors() {
    const auto g = slope_oracle::geometry(5, 5, 10.0, 10.0);
    std::vector<float> v(25, 0.0f);
    v[0 * 5 + 4] = -std::ldexp(1.0f, -50);
    v[1 * 5 + 1] = v[1 * 5 + 2] = 1.5f;
    return OneNode{Raster<float>{g, std::move(v)}, slope_oracle::fan({{0.0, 0.0}, {0.0, -40.0}, {40.0, 0.0}})};
}

}  // namespace

TEST_CASE("34 test 4 (a): two different errors with one product: the larger error is inserted, as today (M2)",
          "[slope][property][test4]") {
    // F is the first of 1.4, 1.39, ... whose weight 1 / F maps both errors to
    // one product, so only the second key, the error, picks today's node.
    const double e1 = 1.5 + std::ldexp(1.0, -52), e2 = 1.5 + std::ldexp(1.0, -51);
    double far = 1.4;
    for (int k = 0; k < 1000 && e1 * (1.0 / far) != e2 * (1.0 / far); ++k) far -= 0.0001;
    CAPTURE(far);
    REQUIRE(e1 * (1.0 / far) == e2 * (1.0 / far));
    const auto fx = two_errors();
    const auto s = slope_of(fx.dem, SlopeRamp{far, far, 30.0, 30.0});
    const auto uniform = run(fx.dem, fx.start, options(far, false), UniformTolerance{far});
    REQUIRE(uniform.ok());
    REQUIRE(uniform.vertices.size() > 3);
    CHECK(uniform.vertices[3] == Point2{20.0, -10.0});  // node (1, 2): the larger error goes in first
    same(run(fx.dem, fx.start, options(far, false), Sloped<UniformTolerance>{UniformTolerance{far}, &s}), uniform);
}

// ------------------------------------------------------------------- M8

namespace {

// A triangle part that never splits anything when asked (at() = 1000), but
// whose bounds put every error above 1 m in its query range; it counts the
// questions about a triangle whose own nodes are already above N = 1 m (the
// scan's result, recomputed here from the mesh it is asked about).
struct CountingAsks {
    const Raster<float>* dem;
    std::atomic<std::size_t>* over_asked;
    [[nodiscard]] double lowest() const noexcept { return 1.0; }
    [[nodiscard]] double highest() const noexcept { return 1000.0; }
    [[nodiscard]] double at(const terrain::mesh::LatticeMesh& m, std::uint32_t t) const {
        const auto r = terrain::refinement::scan(*dem, m, t);
        if (!r.is_void && r.max_error > 1.0) over_asked->fetch_add(1, std::memory_order_relaxed);
        return 1000.0;
    }
};
static_assert(TolerancePolicy<CountingAsks>);

}  // namespace

TEST_CASE("34 (M8): the triangle part is not asked about a triangle its nodes' slope already splits",
          "[slope][property][laziness]") {
    // Every node held to N = 1 (START = END = 0), so every triangle with an
    // error above 1 m is split by the slope; the triangle part's query range
    // is (1, 1000], so laziness asks it nothing during the loop. Feet off:
    // a foot would rightly ask (section 4.3).
    const auto dem = slope_oracle::v1();
    const auto s = slope_of(dem, SlopeRamp{1.0, 1000.0, 0.0, 0.0});
    std::atomic<std::size_t> over_asked{0};
    const CountingAsks counting{&dem, &over_asked};
    const auto out = run(dem, slope_oracle::outline(dem.geometry()), options(1000.0, false, 4),
                         Sloped<CountingAsks>{counting, &s});
    REQUIRE(out.ok());
    CHECK(out.inserted > 100);
    CHECK(over_asked.load() == 0);
    same(out, run(dem, slope_oracle::outline(dem.geometry()), options(1.0, false, 4), UniformTolerance{1.0}));
}
