// Increment 65 (docs/increments/65-delaunay-as-written.md, section 4): the
// refined mesh is Delaunay where it is written. The last swap pass of
// point_loop decides on the written doubles, a triangle flat as written is
// counted and left out of it, and nothing changes where nothing was wrong.
// INVARIANT-CRITICAL: mutation targets T1 (the check removed), T3 (the count
// removed), T5 (skip the check when any triangle is flat), T4 (break after
// the pass, step 7 removed).
// M4 (step 7 removed) survives all tests, 3000 outlines probed; kept as the
// guard that `max_error` describes the written mesh.
//
// Interface (section 3, "Interfaces"): PointRefineOutcome::flat_as_written,
// flips_as_written (a subset of flips) and illegal_as_written, std::size_t.
//
// The Delaunay oracle here is its own, no margin: for each interior edge not
// in out.edges whose two triangles both have positive exact orientation on
// out.vertices (an edge of a flat triangle is not counted), the edge counts
// when DefaultKernel::incircle on out.vertices says Inside. It decides on the
// written world doubles, the frame of the last swap pass. The two Delaunay
// oracles in support/ decide on grid coordinates recovered from the world
// point and cannot see this defect. The tolerance oracle is
// constraint_foot_oracles::node_findings, recomputed from the output.
//
// Platforms: fixture W1's z values are small integers and its vertices are
// decimetres, so the input is the same everywhere; the output can differ by
// rounding (FMA contraction on arm64, none on x86-64 without -march). Every
// figure of C1-C4 was measured the same on macOS arm64 with contraction on
// and with -ffp-contract=off (increment 66, section 4, "Platforms"). T1
// asserts per outline; feet-off flat_as_written is asserted by C1
// (increment 66).
//
// Increment 66 (docs/increments/66-coincident-written-vertices.md, section
// 4): no two vertices written at one point, and none flat as written. C1-C4,
// tagged [coincident]. "Coincident": the output vertices whose written point
// is bit-equal to an earlier vertex's. INVARIANT-CRITICAL: mutation targets
// R1 removed, R1 moved after the foot search, R2 removed.
// Kill record (planted in a scratch copy of the green commit 20253257, the
// test file touched and its object seen rebuilt each time; ctest 1040 tests):
// - "R1 removed" (refine.hpp, the near_constraint branch before the foot
//   search deleted): killed by C2 and C3.
// - "R1 moved after the foot search" (refine.hpp, the same test run only when
//   no foot was found): killed by C2 and C3.
// - "R2 removed" (constraint_points.hpp, the coincidence-radius test beside
//   orient_sign deleted): killed by C2, C4 and CC7
//   (test_refinement_constraint_points.cpp).
// - "the R1-fired flag replaced by !edge" (refine.hpp, the foot search also
//   skipped for a node found on an edge): survives ctest; killed by the
//   Python test test_cli_start_quality_gain.py::
//   test_minus_1_is_20c1s_mesh_on_the_box[gallery] (the mesh topology hash
//   changes; it passes with the flag restored), the only one of 7911 to fail.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>

#include "constraint_foot_oracles.hpp"
#include "quality_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::TriangleIndices;
using terrain::pred::DefaultKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::CheckPoints;
using terrain::refinement::constraint_check_points;
using terrain::refinement::PointRefineOptions;
using terrain::refinement::PointRefineOutcome;
using terrain::refinement::refine;
using terrain::refinement::refine_points;
using terrain::refinement::refine_strip;
using terrain::refinement::RefineOptions;
using quality_fixtures::Start;

namespace cfo = constraint_foot_oracles;

namespace {

// The oracle, on the written doubles (see the head of the file).
std::size_t written_violations(const PointRefineOutcome& out) {
    const auto& v = out.vertices;
    std::set<std::pair<std::uint32_t, std::uint32_t>> constrained;
    for (const auto& e : out.edges) constrained.insert(std::minmax(e[0], e[1]));
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> sides;
    for (std::size_t t = 0; t < out.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            sides[std::minmax(out.triangles[t][k], out.triangles[t][(k + 1) % 3])].push_back(t);
    auto ccw = [&](std::size_t t) {
        const auto& tri = out.triangles[t];
        return DefaultKernel::orient2d(v[tri[0]], v[tri[1]], v[tri[2]]) == Orientation::CounterClockwise;
    };
    std::size_t bad = 0;
    for (const auto& [e, ts] : sides) {
        if (ts.size() > 2) ++bad;  // an edge in three triangles is no triangulation
        if (ts.size() != 2 || constrained.count(e) != 0 || !ccw(ts[0]) || !ccw(ts[1])) continue;
        bool inside = false;
        for (unsigned s = 0; s < 2 && !inside; ++s) {
            const auto& tri = out.triangles[ts[s]];
            std::uint32_t apex = 0;
            for (const auto x : out.triangles[ts[1 - s]])
                if (x != e.first && x != e.second) apex = x;
            inside = DefaultKernel::incircle(v[tri[0]], v[tri[1]], v[tri[2]], v[apex]) == Incircle::Inside;
        }
        if (inside) ++bad;
    }
    return bad;
}

// ------------------------------------------------------------------ fixture W1

constexpr double kX0 = 500000.0, kY0 = 6600000.0, kTol = 0.5;
constexpr int kOutlines = 100;

RasterGeometry w1_geometry() { return RasterGeometry{kX0, kY0, 10.0, 10.0, 12, 12}; }

Raster<float> w1_dem() {
    const auto g = w1_geometry();
    std::vector<float> z(g.rows() * g.cols());
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c)
            z[r * g.cols() + c] = static_cast<float>(static_cast<int>((r * 7 + c * 13) % 11) - 5);
    return Raster<float>{g, std::move(z)};
}

// h(i, lo, hi) in decimetres from the corner, as section 4 defines it.
double h(int k, int i, int lo, int hi) {
    const std::int64_t v = lo + (static_cast<std::int64_t>(k) * (2 * i + 7) * 7919 + i * 104729) % (hi - lo);
    return static_cast<double>(v) / 10.0;
}

// Outline k's four world vertices, counter-clockwise.
std::vector<Point2> outline(int k) {
    return {Point2{kX0 + h(k, 0, 50, 300), kY0 - h(k, 1, 700, 1050)},
            Point2{kX0 + h(k, 2, 800, 1050), kY0 - h(k, 3, 750, 1050)},
            Point2{kX0 + h(k, 4, 700, 1050), kY0 - h(k, 5, 50, 300)},
            Point2{kX0 + h(k, 6, 50, 400), kY0 - h(k, 7, 50, 300)}};
}

// Outline k as a fan {0,1,2}, {0,2,3}, edges (i, i+1) mask 1, plus an
// optional second piece of four vertices in two triangles (T5).
struct Piece {
    std::vector<Point2> xy;
    std::vector<TriangleIndices> tris;
    quality_fixtures::Edges edges;
};

Start w1_start(int k, const Piece* extra = nullptr) {
    std::vector<Point2> xy = outline(k);
    std::vector<TriangleIndices> tris{{0, 1, 2}, {0, 2, 3}};
    quality_fixtures::Edges edges{{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    if (extra != nullptr) {
        const auto base = static_cast<std::uint32_t>(xy.size());
        xy.insert(xy.end(), extra->xy.begin(), extra->xy.end());
        for (const auto& t : extra->tris) tris.push_back({t[0] + base, t[1] + base, t[2] + base});
        for (const auto& e : extra->edges) edges.push_back({e[0] + base, e[1] + base});
    }
    const std::size_t nt = tris.size(), ne = edges.size();
    return Start{IndexedMesh2{std::move(xy), std::move(tris), std::vector<std::uint8_t>(nt, 0)}, std::move(edges),
                 std::vector<std::uint32_t>(ne, 1)};
}

// The CLI's chain: refine (phase 1), constraint_check_points on its output,
// refine_strip; one thread, tolerance and feet as given on both stages.
// `phase1`, if given, gets phase 1's triangles (T4). refine_strip's outcome
// is returned as it is, refused or not (C2 counts a refused run); w1_run
// below requires it ok.
PointRefineOutcome w1_try(const Raster<float>& dem, const Start& s, bool feet, double tol,
                          std::vector<TriangleIndices>* phase1 = nullptr) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = 1;
    o.constraint_feet = feet;
    auto r = refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                    std::span<const std::uint32_t>{s.masks}, o);
    REQUIRE(r.ok());
    const std::size_t nt = r.triangles.size();
    if (phase1 != nullptr) *phase1 = r.triangles;
    const IndexedMesh2 mesh{std::move(r.vertices), std::move(r.triangles), std::vector<std::uint8_t>(nt, 0)};
    const auto strip = constraint_check_points(dem, mesh.vertices(), std::span<const std::array<std::uint32_t, 2>>{r.edges});
    PointRefineOptions po;
    po.tolerance = tol;
    po.threads = 1;
    po.constraint_feet = feet;
    return refine_strip(dem, strip, mesh, std::span<const double>{r.z}, std::span<const std::uint8_t>{r.valid},
                        std::span<const std::array<std::uint32_t, 2>>{r.edges},
                        std::span<const std::uint32_t>{r.masks}, po);
}

PointRefineOutcome w1_run(const Raster<float>& dem, const Start& s, bool feet, double tol,
                          std::vector<TriangleIndices>* phase1 = nullptr) {
    auto out = w1_try(dem, s, feet, tol, phase1);
    INFO("refine_strip refused: " << out.message);
    REQUIRE(out.ok());
    return out;
}

// Rotation-normalised triangles, for "as given".
std::set<std::array<std::uint32_t, 3>> normalised(const std::vector<TriangleIndices>& ts) {
    std::set<std::array<std::uint32_t, 3>> out;
    for (const auto& t : ts) {
        const unsigned m = t[0] < t[1] ? (t[0] < t[2] ? 0 : 2) : (t[1] < t[2] ? 1 : 2);
        out.insert({t[m], t[(m + 1) % 3], t[(m + 2) % 3]});
    }
    return out;
}

// The outlines whose `bad` is true, for a failure message that names k.
struct Tally {
    std::size_t total = 0;
    std::vector<int> ks;
    void add(int k, std::size_t n) {
        total += n;
        if (n != 0) ks.push_back(k);
    }
    [[nodiscard]] std::string names() const {
        std::ostringstream os;
        os << total << " in " << ks.size() << " outlines, k =";
        for (const int k : ks) os << ' ' << k;
        return os.str();
    }
};

// Both oracles and the three counts over all 100 outlines.
struct Sweep {
    Tally delaunay, over, not_ccw, flat, flips_written, illegal;
    // (k, worst - max_error) where max_error < worst, over the triangles not
    // in phase 1's output: refine_strip's max_error covers DEM nodes only in
    // the triangles it wrote (scan_one's `written`), and every triangle not in
    // its start was written by an insertion, a flip or (step 7) the pass.
    std::vector<std::pair<int, double>> max_error_short;
};

// The output restricted to the triangles not in `phase1`.
PointRefineOutcome new_triangles(const PointRefineOutcome& out, const std::vector<TriangleIndices>& phase1) {
    const auto old = normalised(phase1);
    PointRefineOutcome n;
    n.vertices = out.vertices;
    n.z = out.z;
    n.valid = out.valid;
    for (const auto& t : out.triangles)
        if (old.count(*normalised({t}).begin()) == 0) n.triangles.push_back(t);
    return n;
}

Sweep sweep(bool feet) {
    const auto dem = w1_dem();
    Sweep s;
    for (int k = 0; k < kOutlines; ++k) {
        CAPTURE(k);
        std::vector<TriangleIndices> phase1;
        const auto out = w1_run(dem, w1_start(k), feet, kTol, &phase1);
        s.delaunay.add(k, written_violations(out));
        const auto nf = cfo::node_findings(dem, out, kTol);
        s.over.add(k, nf.over);
        s.not_ccw.add(k, nf.not_ccw);
        // max_error is the loop's own; worst is recomputed here. Slack 1e-9 m
        // for the two evaluation orders, on |z| <= 5 m (W1's z range).
        const double worst = cfo::node_findings(dem, new_triangles(out, phase1), kTol).worst;
        if (out.max_error + 1e-9 < worst) s.max_error_short.emplace_back(k, worst - out.max_error);
        s.flat.add(k, out.flat_as_written);
        s.flips_written.add(k, out.flips_as_written);
        s.illegal.add(k, out.illegal_as_written);
    }
    return s;
}

}  // namespace

// ------------------------------------------------------------------ T1

TEST_CASE("T1: with constraint feet on, every W1 outline's output is Delaunay as written, none flat, none left illegal",
          "[refine_points][written_delaunay][T1]") {
    const auto s = sweep(true);
    {
        INFO("Delaunay as written: " << s.delaunay.names());
        CHECK(s.delaunay.total == 0);
    }
    {
        INFO("flat_as_written: " << s.flat.names());
        CHECK(s.flat.total == 0);
    }
    {
        INFO("illegal_as_written: " << s.illegal.names());
        CHECK(s.illegal.total == 0);
    }
    {
        INFO("tolerance oracle, (node, triangle) pairs over " << kTol << " m: " << s.over.names());
        CHECK(s.over.total == 0);
    }
}

// ------------------------------------------------------------------ T2

TEST_CASE("T2: with constraint feet off, the check swaps nothing and the output is Delaunay as written",
          "[refine_points][written_delaunay][T2]") {
    // Pins K1 (a run that needs no fix is bit-identical): no pass runs.
    // Feet-off flat_as_written is asserted by C1 (increment 66).
    const auto s = sweep(false);
    {
        INFO("flips_as_written: " << s.flips_written.names());
        CHECK(s.flips_written.total == 0);
    }
    {
        INFO("illegal_as_written: " << s.illegal.names());
        CHECK(s.illegal.total == 0);
    }
    {
        INFO("Delaunay as written: " << s.delaunay.names());
        CHECK(s.delaunay.total == 0);
    }
    {
        INFO("tolerance oracle, (node, triangle) pairs over " << kTol << " m: " << s.over.names());
        CHECK(s.over.total == 0);
    }
}

// ------------------------------------------------------------------ T3

TEST_CASE("T3: a triangle collinear as written and counter-clockwise on the grid is counted and kept as given",
          "[refine_points][written_delaunay][T3]") {
    const RasterGeometry g{0.0, 30.0, 3.0, 3.0, 11, 11};
    // A, B, C on the line x = y as written (exact orientation 0); D left of
    // A -> C, so (A, C, D) is counter-clockwise.
    const Point2 a{0.5, 0.5}, b{0.51, 0.51}, c{0.52, 0.52}, d{0.4, 0.6};
    REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::Collinear);
    REQUIRE(DefaultKernel::orient2d(a, c, d) == Orientation::CounterClockwise);
    // Precondition, the loop's grid frame ((x - x_min) / dx, ...) scaled back
    // to (col dx, -(row dy)): A, B, C counter-clockwise there (about 9.25e-20).
    auto grid = [&](Point2 p) {
        return Point2{((p.x - g.x_min()) / g.delta_x()) * g.delta_x(), -(((g.y_max() - p.y) / g.delta_y()) * g.delta_y())};
    };
    REQUIRE(DefaultKernel::orient2d(grid(a), grid(b), grid(c)) == Orientation::CounterClockwise);
    const std::vector<TriangleIndices> tris{{0, 1, 2}, {0, 2, 3}};
    const IndexedMesh2 start{std::vector<Point2>{a, b, c, d}, tris, std::vector<std::uint8_t>(2, 0)};
    const quality_fixtures::Edges edges{{0, 1}, {1, 2}, {2, 0}, {0, 3}, {2, 3}};
    const std::vector<std::uint32_t> masks(edges.size(), 1);
    const std::vector<double> z(4, 0.0);
    const std::vector<std::uint8_t> valid(4, 1);
    CheckPoints store{g};
    store.freeze();
    PointRefineOptions o;
    o.tolerance = kTol;
    o.threads = 1;
    const auto out = refine_points(store, start, std::span<const double>{z}, std::span<const std::uint8_t>{valid},
                                   std::span<const std::array<std::uint32_t, 2>>{edges},
                                   std::span<const std::uint32_t>{masks}, o);
    REQUIRE(out.ok());
    CHECK(out.flat_as_written == 1);
    CHECK(out.vertices == std::vector<Point2>{a, b, c, d});
    CHECK(normalised(out.triangles) == normalised(tris));
}

// ------------------------------------------------------------------ T4

TEST_CASE("T4: on T1's outlines every DEM node is within tolerance of the written mesh, and max_error covers it",
          "[refine_points][written_delaunay][T4]") {
    const auto s = sweep(true);
    {
        INFO("tolerance oracle, (node, triangle) pairs over " << kTol << " m: " << s.over.names());
        CHECK(s.over.total == 0);
    }
    {
        INFO("triangles not counter-clockwise in (col, -row): " << s.not_ccw.names());
        CHECK(s.not_ccw.total == 0);
    }
    for (const auto& [k, short_by] : s.max_error_short) UNSCOPED_INFO("k = " << k << ": max_error short by " << short_by);
    CHECK(s.max_error_short.empty());
}

// ------------------------------------------------------------------ T5

TEST_CASE("T5: a triangle flat as written does not stop the pass on the rest of the mesh",
          "[refine_points][written_delaunay][T5]") {
    // Multiples of 1/64, collinear as written, counter-clockwise on the grid
    // (about +1.08e-20); D makes (A, C, D) counter-clockwise.
    const Point2 a{500000.515625, 6599999.4375}, b{500000.53125, 6599999.421875}, c{500000.546875, 6599999.40625},
        d{500001.5, 6599999.75};
    REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::Collinear);
    REQUIRE(DefaultKernel::orient2d(a, c, d) == Orientation::CounterClockwise);
    const Piece piece{{a, b, c, d}, {{0, 1, 2}, {0, 2, 3}}, {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {2, 3}}};
    const auto dem = w1_dem();
    for (const int k : {0, 10}) {
        CAPTURE(k);
        const auto out = w1_run(dem, w1_start(k, &piece), true, kTol);
        CHECK(out.flat_as_written == 1);
        CHECK(out.flips_as_written > 0);
        CHECK(written_violations(out) == 0);
        CHECK(cfo::node_findings(dem, out, kTol).over == 0);
    }
}

// ================================================================== increment 66

namespace {

// The output vertices whose written point is bit-equal to an earlier
// vertex's; a point written three times counts 2.
std::size_t coincident(const std::vector<Point2>& vertices) {
    std::set<std::pair<double, double>> seen;
    std::size_t n = 0;
    for (const auto& v : vertices)
        n += seen.insert({v.x, v.y}).second ? 0 : 1;
    return n;
}

// Outline k = 82's edge from vertex 1 to vertex 2 passes, in exact decimals,
// through the DEM node (col 9, row 4), written (500090.0, 6599960.0), z -3.
constexpr int kK82 = 82;
constexpr Point2 kNode94{500090.0, 6599960.0};
constexpr double kNode94Z = -3.0;

}  // namespace

// ------------------------------------------------------------------ C1

TEST_CASE("C1: at tolerance 0.5, no W1 outline has two vertices written at one point or a triangle flat as written",
          "[refine_points][written_delaunay][coincident][C1]") {
    const auto dem = w1_dem();
    for (const bool feet : {true, false}) {
        CAPTURE(feet);
        Tally pairs, flat, delaunay, over;
        for (int k = 0; k < kOutlines; ++k) {
            CAPTURE(k);
            const auto out = w1_run(dem, w1_start(k), feet, kTol);
            pairs.add(k, coincident(out.vertices));
            flat.add(k, out.flat_as_written);
            delaunay.add(k, written_violations(out));
            over.add(k, cfo::node_findings(dem, out, kTol).over);
        }
        {
            INFO("coincident vertices: " << pairs.names());
            CHECK(pairs.total == 0);
        }
        {
            INFO("flat_as_written: " << flat.names());
            CHECK(flat.total == 0);
        }
        {
            INFO("Delaunay as written: " << delaunay.names());
            CHECK(delaunay.total == 0);
        }
        {
            INFO("tolerance oracle, (node, triangle) pairs over " << kTol << " m: " << over.names());
            CHECK(over.total == 0);
        }
    }
}

// ------------------------------------------------------------------ C2

TEST_CASE("C2: at tolerance 0, the ten W1 outlines that fail on master run, with no coincident vertex and none flat",
          "[refine_points][written_delaunay][coincident][C2]") {
    // The other 2990 outlines of 0..2999 are 0 on master and stay 0 (probed
    // by the design, not asserted here).
    const auto dem = w1_dem();
    for (const bool feet : {true, false}) {
        CAPTURE(feet);
        std::ostringstream failing;
        std::size_t bad = 0, delaunay = 0, over = 0;
        for (const int k : {82, 1172, 1437, 1457, 1922, 2137, 2157, 2622, 2852, 2922}) {
            const auto out = w1_try(dem, w1_start(k), feet, 0.0);
            if (!out.ok()) {
                ++bad;
                failing << " k = " << k << " refused (" << out.message << ");";
                continue;
            }
            const std::size_t pairs = coincident(out.vertices);
            if (pairs != 0 || out.flat_as_written != 0) {
                ++bad;
                failing << " k = " << k << " (coincident " << pairs << ", flat " << out.flat_as_written << ");";
            }
            delaunay += written_violations(out);
            // node_findings' own slack, 1e-9 * max|z| (5 m on W1), makes
            // tolerance 0 checkable.
            over += cfo::node_findings(dem, out, 0.0).over;
        }
        INFO("outlines refused, with a coincident vertex or flat:" << failing.str());
        CHECK(bad == 0);
        CHECK(delaunay == 0);
        CHECK(over == 0);
    }
}

// ------------------------------------------------------------------ C3

TEST_CASE("C3: phase 1 inserts the DEM node a hair off a constraint on the constraint, not its foot",
          "[refine][coincident][C3]") {
    const auto dem = w1_dem();
    const Start s = w1_start(kK82);
    RefineOptions o;
    o.tolerance = kTol;
    o.threads = 1;
    o.constraint_feet = true;
    const auto r = refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                          std::span<const std::uint32_t>{s.masks}, o);
    REQUIRE(r.ok());
    // The vertex nearest the node; 1e-6 m on 10 m cells, coordinates near
    // 6.6e6 m (an ulp there is 9.3e-10 m), checked on W1 only.
    std::size_t near = r.vertices.size();
    double best = 1e-6;
    for (std::size_t i = 0; i < r.vertices.size(); ++i) {
        const double d = std::hypot(r.vertices[i].x - kNode94.x, r.vertices[i].y - kNode94.y);
        if (d < best) {
            best = d;
            near = i;
        }
    }
    REQUIRE(near < r.vertices.size());
    CAPTURE(r.vertices[near].x, r.vertices[near].y, r.z[near]);
    CHECK(r.vertices[near] == kNode94);
    CHECK(r.z[near] == kNode94Z);
    const auto ends = std::count_if(r.edges.begin(), r.edges.end(),
                                    [&](const auto& e) { return e[0] == near || e[1] == near; });
    CHECK(ends == 2);
}

// ------------------------------------------------------------------ C4

TEST_CASE("C4: a crossing whose node lies within the coincidence radius of the edge is that node",
          "[constraint_points][coincident][C4]") {
    const auto dem = w1_dem();
    const std::vector<Point2> ends{Point2{500094.6, 6599920.9}, Point2{500088.6, 6599971.9}};
    const std::vector<std::array<std::uint32_t, 2>> edge{{0, 1}};
    const auto strip = constraint_check_points(dem, std::span<const Point2>{ends},
                                               std::span<const std::array<std::uint32_t, 2>>{edge});
    REQUIRE(strip.edge_count() == 1);
    // Points within 1e-6 lattice units of (9, 4): far wider than the
    // coincidence radius (1e-10 on a 12 x 12 grid), far narrower than a cell.
    std::vector<terrain::refinement::ConstraintPoint> there;
    std::ostringstream at;
    at.precision(17);
    for (const auto& p : strip.on_edge(0))
        if (std::hypot(p.at.col - 9.0, p.at.row - 4.0) < 1e-6) {
            there.push_back(p);
            at << " (" << p.at.col << ", " << p.at.row << ")";
        }
    INFO("points near (9, 4):" << at.str());
    REQUIRE(there.size() == 1);
    CHECK(there[0].at.col == 9.0);
    CHECK(there[0].at.row == 4.0);
    CHECK(there[0].z == kNode94Z);
    CHECK(strip.size() == 11);
    CHECK(strip.duplicates() == 1);
}
