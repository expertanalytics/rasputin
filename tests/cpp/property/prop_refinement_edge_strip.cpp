// Increment 15f-2 (docs/increments/15f-edge-strip.md, D4, "The guarantee, and
// its wording" E1 to E8, and "Tests for @tester" ES1 to ES11): the edge strip
// in the refinement loop, through the C++ API.
//
//   refine_strip(dem, strip, start, z, valid, edges, masks, options)
//       the projected path: the strip's points and the DEM's nodes in every
//       triangle the run writes, one loop;
//   refine_points(store, start, z, valid, edges, masks, options, &strip)
//       the reprojected path: the source points and the strip, one loop;
//
// both returning PointRefineOutcome with strip_points, strip_inserted,
// strip_max_error, strip_refused, strip_refused_max_error and nodes_inserted
// (D4, "Outcome").
//
// Invariant-critical (the record): ES2 (E1 by an independent oracle), ES3
// (E2, the DEM nodes kept) and ES5 (ownership and the sub-edge bookkeeping).
//
// Oracles (support/strip_oracle.hpp; each is shown to fail in
// unit/test_refinement_strip_oracle.cpp). Every property case below carries
// both of the QA rules' section D oracles on the projected path:
//   - the constrained-Delaunay oracle (delaunay_violations), the exact
//     incircle in the producer's frame (col dx, -(row dy));
//   - the tolerance oracle (node_findings): every valid DEM node in every
//     closed output triangle with three valid vertices within tolerance of
//     the plane recomputed from the output. This is E2.
// and the strip's own guarantee:
//   - the strip oracle (ruled_points + strip_findings): every crossing of a
//     START constraint edge with a grid line, and every midpoint between two
//     neighbouring ones (ends counting), generated here, filed against the
//     OUTPUT constraint edge that holds it by its own projection, within
//     tolerance of the output's linear z along that edge. This is E1. It
//     reads neither the store's order nor the loop's records.
// On the reprojected path the tolerance oracle is 15c's J2 over the source
// points (sources_over below), since J2 gives up the DEM-node guarantee there.
//
// CHOSEN HERE, where D4 is silent (see the handback):
//   - the logic_error refusals of D4 name their entry point in what():
//     "refine_strip" or "refine_points", in refine's "refine: ..." style;
//   - refine_strip refuses a tolerance that is not finite and >= 0 with
//     RefineStatus::InvalidTolerance, as refine and refine_points do;
//   - a triangle whose named strip point is refused stays in the run, so the
//     points behind the refused one are still repaired (ES6, "E1 holds for
//     the rest").
//
// Q1 (the every-point form, recursive midpoints) is NOT ruled, so the strip
// oracle checks the ruled points only. If Ola rules yes, ES2 alone changes:
// see the block marked "Q1" in it.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>

#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace strip_oracle;
using terrain::Point2;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::refinement::CheckPoints;
using terrain::refinement::constraint_check_points;
using terrain::refinement::ConstraintCheckPoints;
using terrain::refinement::PointRefineOptions;
using terrain::refinement::PointRefineOutcome;
using terrain::refinement::refine;
using terrain::refinement::refine_points;
using terrain::refinement::refine_strip;
using terrain::refinement::RefineOptions;
using terrain::refinement::RefineStatus;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;

namespace {

// What a strip run starts from: a mesh in world coordinates, its z and valid
// flags, and its constraint edges with their masks.
struct Begin {
    IndexedMesh2 mesh;
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    Edges edges;
    std::vector<std::uint32_t> masks;
};

// Phase 1 on the projected path: refine(dem) from `s`, as _dem_mesh runs it.
Begin after_refine(const Raster<float>& dem, const Start& s, double tol) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = 1;
    auto r = refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                    std::span<const std::uint32_t>{s.masks}, o);
    REQUIRE(r.ok());
    const std::size_t nt = r.triangles.size();
    return Begin{IndexedMesh2{std::move(r.vertices), std::move(r.triangles), std::vector<std::uint8_t>(nt, 0)},
                 std::move(r.z), std::move(r.valid), std::move(r.edges), std::move(r.masks)};
}

// A start that did not come from refine: z by the oracle's own heights.
Begin direct(const Raster<float>& dem, const Start& s) {
    auto [z, valid] = start_z(dem, s.mesh.vertices());
    return Begin{s.mesh, std::move(z), std::move(valid), s.edges, s.masks};
}

ConstraintCheckPoints strip_of(const Raster<float>& dem, const Begin& b) {
    return constraint_check_points(dem, b.mesh.vertices(), std::span<const std::array<std::uint32_t, 2>>{b.edges});
}

PointRefineOutcome run_strip(const Raster<float>& dem, const ConstraintCheckPoints& strip, const Begin& b,
                             double tol, unsigned threads = 1) {
    PointRefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    return refine_strip(dem, strip, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                        std::span<const std::array<std::uint32_t, 2>>{b.edges},
                        std::span<const std::uint32_t>{b.masks}, o);
}

PointRefineOutcome run_points(const CheckPoints& store, const ConstraintCheckPoints* strip, const Begin& b,
                              double tol, unsigned threads = 1) {
    PointRefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    return refine_points(store, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                         std::span<const std::array<std::uint32_t, 2>>{b.edges},
                         std::span<const std::uint32_t>{b.masks}, o, strip);
}

// E1 on `out`, from the start `b`'s constraint edges.
StripFindings e1(const Raster<float>& dem, const Begin& b, const Mesh& out, double tol) {
    return strip_findings(dem.geometry(), ruled_points(dem, b.mesh.vertices(), b.edges), out, tol);
}

// Every finding of the projected path's oracles that must be zero, and the
// outcome's own bookkeeping, for one run.
void projected_oracles(const Raster<float>& dem, const ConstraintCheckPoints& strip, const Begin& b,
                       const PointRefineOutcome& out, double tol) {
    const Mesh m = mesh_of(out);
    const auto& g = dem.geometry();
    REQUIRE(out.ok());
    CHECK(out.strip_points == strip.size());
    CHECK(out.strip_max_error <= tol);
    CHECK(out.max_error <= tol);
    CHECK(out.uncovered == 0);
    CHECK(out.vertices.size() == b.mesh.vertices().size() + out.inserted);
    CHECK(out.strip_inserted + out.nodes_inserted == out.inserted);

    const auto f = e1(dem, b, m, tol);  // E1
    CHECK(f.unfiled == 0);
    CHECK(f.on_void == 0);
    CHECK(f.over == 0);

    const auto n = node_findings(dem, m, tol);  // E2: section D's tolerance oracle
    CHECK(n.not_ccw == 0);
    CHECK(n.over == 0);
    CHECK(delaunay_violations(g, m) == 0);  // section D's constrained-Delaunay oracle

    const auto s = shape_findings(dem, b.mesh.vertices(), b.edges, b.masks, m);  // E5, E6, E8
    CHECK(s.moved_start == 0);
    CHECK(s.stray == 0);
    CHECK(s.wrong_z == 0);
    CHECK(s.coincident == 0);
    CHECK(s.unplaced_edge == 0);
    CHECK(s.wrong_mask == 0);
    CHECK(s.broken_chain == 0);
    CHECK(s.off_node_new <= out.strip_inserted);
}

void require_identical(const PointRefineOutcome& a, const PointRefineOutcome& b) {
    REQUIRE(a.vertices.size() == b.vertices.size());
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        CAPTURE(i);
        REQUIRE(a.vertices[i].x == b.vertices[i].x);
        REQUIRE(a.vertices[i].y == b.vertices[i].y);
        REQUIRE(a.z[i] == b.z[i]);
        REQUIRE(a.valid[i] == b.valid[i]);
    }
    REQUIRE(a.triangles == b.triangles);
    REQUIRE(a.edges == b.edges);
    REQUIRE(a.masks == b.masks);
    CHECK(a.inserted == b.inserted);
    CHECK(a.rounds == b.rounds);
    CHECK(a.strip_inserted == b.strip_inserted);
    CHECK(a.nodes_inserted == b.nodes_inserted);
    CHECK(a.strip_refused == b.strip_refused);
    CHECK(a.max_error == b.max_error);
    CHECK(a.strip_max_error == b.strip_max_error);
}

const char* name(Terrain t) {
    switch (t) {
    case Terrain::Rough: return "rough";
    case Terrain::Smooth: return "smooth";
    case Terrain::SmoothWithHole: return "smooth with a NoData patch";
    case Terrain::Plane: return "plane";
    }
    return "?";
}

}  // namespace

// ---------------------------------------------------------------------------
// ES1: the defect, then its repair
// ---------------------------------------------------------------------------

TEST_CASE("ES1: refine_strip repairs the crossings refine leaves over the tolerance",
          "[edge_strip][ES1]") {
    // The defect half (refine alone, the oracle over 0) is in
    // unit/test_refinement_strip_oracle.cpp, so it runs before refine_strip exists.
    const auto dem = bump_dem();
    const double tol = 1.0;
    const Begin b = after_refine(dem, bump_start(), tol);
    REQUIRE(e1(dem, b, Mesh{{b.mesh.vertices().begin(), b.mesh.vertices().end()}, b.z, b.valid,
                            {b.mesh.triangles().begin(), b.mesh.triangles().end()}, b.edges, b.masks},
               tol).over > 0);
    const auto strip = strip_of(dem, b);
    const auto out = run_strip(dem, strip, b, tol);
    projected_oracles(dem, strip, b, out, tol);
    CHECK(out.strip_inserted > 0);
    CHECK(out.strip_refused == 0);
    // F2: the crossing of column 3 goes in at z 5 beside nodes at 0, so the
    // rescan must have put DEM nodes in the triangles it wrote.
    CHECK(out.nodes_inserted > 0);
}

// ---------------------------------------------------------------------------
// ES2 and ES3: E1 and E2 by independent oracles, over seeded inputs
// ---------------------------------------------------------------------------

TEST_CASE("ES2: every ruled strip point is within tolerance after refine_strip (E1)",
          "[edge_strip][ES2][property]") {
    const auto g = exact_geometry();
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth, Terrain::SmoothWithHole);
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    const double tol = GENERATE(0.0, 0.5, 2.0, 8.0);
    CAPTURE(name(kind), seed, tol);
    const auto dem = terrain_dem(g, kind, seed);
    const Begin b = after_refine(dem, jittered(g, seed), tol);
    const auto strip = strip_of(dem, b);
    const auto out = run_strip(dem, strip, b, tol);
    projected_oracles(dem, strip, b, out, tol);
    CHECK(out.strip_refused == 0);  // no end is within rounding of a grid line here
    if (kind == Terrain::Rough) CHECK(out.strip_inserted > 0);

    // The oracle can fail on this output: the strip's vertices (new and off
    // a node) shifted by twice the tolerance (1 m at tolerance 0).
    Mesh planted = mesh_of(out);
    std::size_t shifted = 0;
    for (std::size_t i = b.mesh.vertices().size(); i < planted.vertices.size(); ++i)
        if (!is_node(lat(g, planted.vertices[i]))) {
            planted.z[i] += std::max(2 * tol, 1.0);
            ++shifted;
        }
    if (shifted > 0) CHECK(e1(dem, b, planted, tol).over > 0);

    // Q1. Ruled no (not ruled, default as designed): only the ruled points
    // are checked. If Ola rules the every-point form, this case adds a dense
    // sample of every output constraint edge, against the bilinear surface,
    // within 1.25 x tolerance (F4), or within tolerance on constraints along
    // grid lines; no other case changes.
}

TEST_CASE("ES3: every valid DEM node is within tolerance after refine_strip (E2)",
          "[edge_strip][ES3][property]") {
    // projected_oracles carries E2 for every ES2 input; this case adds the
    // inputs on which a strip insertion is certain to disturb a DEM node (the
    // two bump sections, by construction), so that a run without the rescan
    // of F2 fails here whatever the seeds do.
    SECTION("the bump, as refine leaves it: the crossings go in beside nodes at 0") {
        const auto dem = bump_dem();
        const double tol = GENERATE(0.0, 0.5, 1.0, 3.0);
        CAPTURE(tol);
        const Begin b = after_refine(dem, bump_start(), tol);
        const auto strip = strip_of(dem, b);
        const auto out = run_strip(dem, strip, b, tol);
        projected_oracles(dem, strip, b, out, tol);
        CHECK(out.nodes_inserted > 0);
    }
    SECTION("the bump side interior: the 10 m nodes lie below it, in a triangle the strip writes") {
        const auto dem = bump_dem();
        const double tol = GENERATE(0.5, 1.0, 3.0);
        CAPTURE(tol);
        const Begin b = direct(dem, bump_start(true));
        const auto strip = strip_of(dem, b);
        const auto out = run_strip(dem, strip, b, tol);
        projected_oracles(dem, strip, b, out, tol);
        CHECK(out.nodes_inserted > 0);
    }
    SECTION("seeded, rough: the rescan has work on these seeds") {
        const auto g = exact_geometry();
        std::size_t nodes = 0;
        for (const std::uint32_t seed : {4u, 5u, 6u, 7u})
            for (const double tol : {1.0, 4.0}) {
                CAPTURE(seed, tol);
                const auto dem = terrain_dem(g, Terrain::Rough, seed);
                const Begin b = after_refine(dem, jittered(g, seed), tol);
                const auto strip = strip_of(dem, b);
                const auto out = run_strip(dem, strip, b, tol);
                projected_oracles(dem, strip, b, out, tol);
                nodes += out.nodes_inserted;
            }
        CHECK(nodes > 0);
    }
}

// ---------------------------------------------------------------------------
// ES4: determinism (E4)
// ---------------------------------------------------------------------------

TEST_CASE("ES4: the output is bit-identical over threads, edge-list order and pair order",
          "[edge_strip][ES4][property]") {
    const auto g = exact_geometry();
    const auto kind = GENERATE(Terrain::Rough, Terrain::SmoothWithHole);
    const std::uint32_t seed = GENERATE(1u, 2u);
    const double tol = GENERATE(0.0, 0.5);
    CAPTURE(name(kind), seed, tol);
    const auto dem = terrain_dem(g, kind, seed);
    const Begin b = after_refine(dem, jittered(g, seed), tol);
    const auto reference = run_strip(dem, strip_of(dem, b), b, tol, 1);
    REQUIRE(reference.ok());
    REQUIRE(reference.inserted > 0);

    for (const unsigned threads : {2u, 8u}) {
        CAPTURE(threads);
        require_identical(reference, run_strip(dem, strip_of(dem, b), b, tol, threads));
    }

    // The edge list permuted (masks with it) and every pair reversed, for the
    // generator and the run alike.
    Begin shuffled = b;
    std::vector<std::size_t> order(b.edges.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::mt19937 gen{seed};
    for (std::size_t i = order.size(); i > 1; --i) std::swap(order[i - 1], order[gen() % i]);
    for (std::size_t i = 0; i < order.size(); ++i) {
        shuffled.edges[i] = {b.edges[order[i]][1], b.edges[order[i]][0]};
        shuffled.masks[i] = b.masks[order[i]];
    }
    for (const unsigned threads : {1u, 8u}) {
        CAPTURE(threads);
        require_identical(reference, run_strip(dem, strip_of(dem, shuffled), shuffled, tol, threads));
    }
}

// ---------------------------------------------------------------------------
// ES5: ownership and the sub-edge bookkeeping
// ---------------------------------------------------------------------------

TEST_CASE("ES5: a boundary edge run high to low is scanned, an interior one repaired once",
          "[edge_strip][ES5]") {
    const auto dem = bump_dem();
    const auto& g = dem.geometry();
    const bool interior = GENERATE(false, true);
    const double tol = GENERATE(0.0, 1.0);
    CAPTURE(interior, tol);
    // The bottom side v3 -> v0 (bump_start): in its only triangle (boundary)
    // or in the upper one (interior) it runs from the higher index to the lower.
    const Begin b = direct(dem, bump_start(interior));
    const auto strip = strip_of(dem, b);
    const auto out = run_strip(dem, strip, b, tol);
    projected_oracles(dem, strip, b, out, tol);

    // The bottom side's own points: seven crossings and eight midpoints, all
    // off-node at row 3.5. Its new vertices are some of them, each once.
    const auto mine = ruled_points(dem, b.mesh.vertices(), Edges{{0, 3}});
    REQUIRE(mine.size() == 15);
    const auto f = strip_findings(g, mine, mesh_of(out), tol);
    CHECK(f.over == 0);
    std::set<Lat> on_side;
    for (std::size_t i = b.mesh.vertices().size(); i < out.vertices.size(); ++i) {
        const Lat p = lat(g, out.vertices[i]);
        if (p.row == 3.5) {
            CHECK(on_side.insert(p).second);
            CHECK(std::any_of(mine.begin(), mine.end(), [&](const OraclePoint& q) { return q.at == p; }));
        }
    }
    CHECK(!on_side.empty());
    CHECK(on_side.size() <= 15);
    // Every edge of this fixture is off-node, so every off-node new vertex is a
    // strip point and every strip point inserted is off-node.
    CHECK(shape_findings(dem, b.mesh.vertices(), b.edges, b.masks, mesh_of(out)).off_node_new == out.strip_inserted);
}

// ---------------------------------------------------------------------------
// ES6: a refused insertion
// ---------------------------------------------------------------------------

namespace {

// The CC7 probe (P0 one lattice ulp before column 1) as a one-triangle mesh:
// a = (nextafter(1, 0), 1.5), b = (4, 3.25), and an apex c. The crossing of
// column 1 is the strip point f = (1, 1.5), one ulp from a, and f - a is
// horizontal. Splitting a-b at f makes the child (a, f, c), which is strictly
// counter-clockwise exactly when c's row is above a's: with c = (3, 2.296875)
// (inside the wedge between a-b and the row through a) foot_fits refuses f;
// with c = (3, 1) it does not. The DEM is a plane and a's z is lifted 10 m,
// so the error along a-b falls linearly from a: f is the worst point on a-b
// and is named in the first round, while the start triangle is the only one.
struct Probe {
    Raster<float> dem;
    Begin begin;
};

Probe ulp_probe(double apex_row) {
    const auto g = exact_geometry(9, 9);
    std::vector<float> v(81);
    for (std::size_t r = 0; r < 9; ++r)
        for (std::size_t c = 0; c < 9; ++c) v[r * 9 + c] = static_cast<float>(3 * c + 2 * r);
    Raster<float> dem{g, std::move(v)};
    Start s;
    s.mesh = IndexedMesh2{{world(g, std::nextafter(1.0, 0.0), 1.5), world(g, 4.0, 3.25), world(g, 3.0, apex_row)},
                          {{0, 1, 2}},
                          {0}};
    s.edges = {{0, 1}};
    s.masks = {1};
    Begin b = direct(dem, s);
    b.z[0] += 10.0;
    return Probe{std::move(dem), std::move(b)};
}

}  // namespace

TEST_CASE("ES6: a strip point foot_fits refuses is counted with its error, and the rest is repaired",
          "[edge_strip][ES6]") {
    const double tol = 0.5;
    SECTION("the apex in the wedge: refused") {
        const auto p = ulp_probe(2.296875);
        const auto& g = p.dem.geometry();
        const auto strip = strip_of(p.dem, p.begin);
        REQUIRE(strip.on_edge(0).size() > 0);
        REQUIRE(strip.on_edge(0)[0].at == terrain::mesh::MeshVertex{1.0, 1.5});
        const auto out = run_strip(p.dem, strip, p.begin, tol);
        REQUIRE(out.ok());
        CHECK(out.strip_refused == 1);
        CHECK(out.strip_refused_max_error == Catch::Approx(10.0).margin(1e-9));
        CHECK(out.strip_max_error <= tol);  // refused points excluded
        CHECK(out.strip_inserted > 0);      // the points behind f still went in
        const Mesh m = mesh_of(out);
        const auto f = e1(p.dem, p.begin, m, tol);
        CHECK(f.unfiled == 0);
        REQUIRE(f.over > 0);
        const Lat a = lat(g, p.begin.mesh.vertices()[0]);
        for (const Lat q : f.over_at) {
            CAPTURE(q.col, q.row);
            CHECK(std::hypot(q.col - a.col, q.row - a.row) <= 1e-9);  // only at f (and its midpoint with a)
        }
        // No vertex at f; every triangle strictly counter-clockwise; E2 and
        // the Delaunay oracle hold.
        for (const Point2 v : out.vertices) CHECK_FALSE(lat(g, v) == (Lat{1.0, 1.5}));
        const auto n = node_findings(p.dem, m, tol);
        CHECK(n.not_ccw == 0);
        CHECK(n.over == 0);
        CHECK(delaunay_violations(g, m) == 0);
    }
    SECTION("the apex above the row of a: f goes in, a hair from a") {
        const auto p = ulp_probe(1.0);
        const auto& g = p.dem.geometry();
        const auto strip = strip_of(p.dem, p.begin);
        const auto out = run_strip(p.dem, strip, p.begin, tol);
        REQUIRE(out.ok());
        CHECK(out.strip_refused == 0);
        const Mesh m = mesh_of(out);
        CHECK(e1(p.dem, p.begin, m, tol).over == 0);
        CHECK(std::any_of(out.vertices.begin(), out.vertices.end(),
                          [&](Point2 v) { return lat(g, v) == Lat{1.0, 1.5}; }));
        const auto n = node_findings(p.dem, m, tol);
        CHECK(n.not_ccw == 0);
        CHECK(n.over == 0);
        CHECK(delaunay_violations(g, m) == 0);
        CHECK(shape_findings(p.dem, p.begin.mesh.vertices(), p.begin.edges, p.begin.masks, m).coincident == 0);
    }
}

// ---------------------------------------------------------------------------
// ES7: a constraint edge with an invalid end
// ---------------------------------------------------------------------------

TEST_CASE("ES7: a void sub-edge is carved along the edge until no strip point is uncovered",
          "[edge_strip][ES7]") {
    // The bump with node (row 4, col 0) NoData: v3 = (0.5, 3.5) has no z, so
    // the bottom side and the left side each have an invalid end.
    const auto g = exact_geometry(9, 9);
    std::vector<float> z(81, 0.0f);
    z[4 * 9 + 3] = 10.0f;
    z[4 * 9 + 4] = 10.0f;
    z[4 * 9 + 0] = -9999.0f;
    const Raster<float> dem{g, std::move(z), -9999.0f};
    const double tol = GENERATE(0.0, 1.0);
    CAPTURE(tol);
    const Begin b = direct(dem, bump_start());
    REQUIRE(b.valid[3] == 0);
    const auto strip = strip_of(dem, b);
    REQUIRE(strip.no_data() > 0);
    const auto out = run_strip(dem, strip, b, tol);
    projected_oracles(dem, strip, b, out, tol);  // uncovered 0 and on_void 0 among them
    CHECK(out.carved > 0);
}

// ---------------------------------------------------------------------------
// ES8: a point of another set on a constrained edge consumes the strip's point
// ---------------------------------------------------------------------------

namespace {

// A square A (0.5, 0.5), B (3.5, 3.5), C (3.5, 0.5), D (0.5, 3.5) cut by the
// constrained diagonal A-B, which passes through nodes (1, 1), (2, 2) and
// (3, 3). Triangle 0 = (A, D, B) holds A-B as B -> A (not its owner);
// triangle 1 = (A, B, C) holds it as A -> B (its owner). The outline is
// constrained too.
Start diagonal_start(const terrain::raster::RasterGeometry& g) {
    Start s;
    s.mesh = IndexedMesh2{{world(g, 0.5, 0.5), world(g, 3.5, 3.5), world(g, 3.5, 0.5), world(g, 0.5, 3.5)},
                          {{0, 3, 1}, {0, 1, 2}},
                          {0, 0}};
    s.edges = {{0, 1}, {0, 2}, {1, 2}, {1, 3}, {0, 3}};
    s.masks = {1, 2, 2, 2, 2};
    return s;
}

Raster<float> node_bump(std::size_t row, std::size_t col, float height) {
    const auto g = exact_geometry(6, 6);
    std::vector<float> z(36, 0.0f);
    z[row * 6 + col] = height;
    return Raster<float>{g, std::move(z)};
}

}  // namespace

TEST_CASE("ES8: a node on a constrained diagonal is inserted once, by whichever set names it",
          "[edge_strip][ES8]") {
    SECTION("projected path: the rescan and the strip both see node (2, 2)") {
        const auto row = GENERATE(as<std::size_t>{}, 1, 2, 3);
        const double tol = GENERATE(0.0, 0.5, 3.0);
        CAPTURE(row, tol);
        // The bump on the diagonal node (row, row). Which set inserts the node
        // depends on the split phase's triangle order, which the design does
        // not fix; either way it goes in once.
        auto dem = node_bump(row, row, 10.0f);
        const auto g = dem.geometry();
        const Begin b = direct(dem, diagonal_start(g));
        const auto strip = strip_of(dem, b);
        const auto out = run_strip(dem, strip, b, tol);
        projected_oracles(dem, strip, b, out, tol);  // coincident 0 among them
        const Lat node{static_cast<double>(row), static_cast<double>(row)};
        CHECK(std::count_if(out.vertices.begin(), out.vertices.end(),
                            [&](Point2 v) { return lat(g, v) == node; }) == 1);
    }
    SECTION("reprojected path: a source point at node (2, 2) wins the tie and consumes the strip's") {
        // Source and strip both name (2, 2) with error 10 (the strip's z is the
        // node's, 10; the source point's z is 10). Ties go to the source set
        // (D4, "Combining"); its split is on the constrained edge, so the strip
        // point at exactly that position is consumed (D4 step 4).
        const double tol = GENERATE(0.0, 1.0);
        CAPTURE(tol);
        const auto dem = node_bump(2, 2, 10.0f);
        const auto g = dem.geometry();
        const Begin b = direct(dem, diagonal_start(g));
        const auto strip = strip_of(dem, b);
        CheckPoints store{g};
        const std::vector<Point2> xy{world(g, 2.0, 2.0)};
        const std::vector<float> zs{10.0f};
        store.add(xy, zs);
        store.freeze();
        const auto out = run_points(store, &strip, b, tol);
        REQUIRE(out.ok());
        const Lat node{2.0, 2.0};
        CHECK(std::count_if(out.vertices.begin(), out.vertices.end(),
                            [&](Point2 v) { return lat(g, v) == node; }) == 1);
        CHECK(out.strip_inserted + 1 == out.inserted);  // the one source point, and strip points
        const Mesh m = mesh_of(out);
        CHECK(shape_findings(dem, b.mesh.vertices(), b.edges, b.masks, m).coincident == 0);
        const auto f = e1(dem, b, m, tol);
        CHECK(f.unfiled == 0);
        CHECK(f.over == 0);
        CHECK(out.strip_max_error <= tol);
        CHECK(delaunay_violations(g, m) == 0);
    }
}

// ---------------------------------------------------------------------------
// ES9: the reprojected path, source points and the strip in one loop
// ---------------------------------------------------------------------------

namespace {

struct Sources {
    std::vector<Point2> xy;
    std::vector<float> z;
};

// Two points per cell at seeded dyadic offsets in (0, 1), z the DEM's bilinear
// value plus noise of up to +-4 m, rounded to 1/64 m so it is a float exactly.
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

// 15c's J2 (RP3's oracle): every source point that is not a start vertex, in
// every closed output triangle with three valid vertices, within tolerance of
// that triangle's plane recomputed from the output.
std::size_t sources_over(const terrain::raster::RasterGeometry& g, const Sources& pts, const Begin& b,
                         const Mesh& out, double tol) {
    std::set<std::pair<double, double>> start_xy;
    for (const Point2 v : b.mesh.vertices()) start_xy.insert({v.x, v.y});
    std::vector<Point2> fp;
    double zmax = 1.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const Lat l = lat(g, out.vertices[i]);
        fp.push_back(Point2{l.col, -l.row});
        zmax = std::max(zmax, std::abs(out.z[i]));
    }
    std::size_t bad = 0;
    for (const auto& t : out.triangles) {
        if (!(out.valid[t[0]] && out.valid[t[1]] && out.valid[t[2]])) continue;
        const Point2 a = fp[t[0]], bb = fp[t[1]], c = fp[t[2]];
        const double two_a = cross(a, bb, c);
        for (std::size_t i = 0; i < pts.xy.size(); ++i) {
            if (start_xy.contains({pts.xy[i].x, pts.xy[i].y})) continue;
            const Lat l = lat(g, pts.xy[i]);
            const Point2 p{l.col, -l.row};
            if (DefaultKernel::orient2d(a, bb, p) == Orientation::Clockwise
                || DefaultKernel::orient2d(bb, c, p) == Orientation::Clockwise
                || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                continue;
            const double plane = (cross(p, bb, c) * out.z[t[0]] + cross(a, p, c) * out.z[t[1]]
                                  + cross(a, bb, p) * out.z[t[2]]) / two_a;
            if (std::abs(plane - static_cast<double>(pts.z[i])) > tol + 1e-9 * zmax) ++bad;
        }
    }
    return bad;
}

}  // namespace

TEST_CASE("ES9: refine_points with the strip keeps J2 at the source points and E1 at the strip's",
          "[edge_strip][ES9][property]") {
    const auto g = exact_geometry();
    const auto kind = GENERATE(Terrain::Smooth, Terrain::SmoothWithHole);
    const std::uint32_t seed = GENERATE(1u, 2u);
    const double tol = GENERATE(0.0, 0.5, 2.0, 8.0);
    CAPTURE(name(kind), seed, tol);
    // The resampled target grid stands in as `dem`: phase 1 refines against
    // it, the strip is generated on it, and the sources are scattered over it.
    const auto dem = terrain_dem(g, kind, seed);
    const Begin b = after_refine(dem, jittered(g, seed), tol);
    const auto strip = strip_of(dem, b);
    const Sources src = scattered(dem, seed);
    CheckPoints store{g};
    store.add(src.xy, src.z);
    store.freeze();

    const auto out = run_points(store, &strip, b, tol);
    REQUIRE(out.ok());
    CHECK(out.strip_points == strip.size());
    CHECK(out.strip_refused == 0);
    CHECK(out.strip_max_error <= tol);
    CHECK(out.max_error <= tol);
    CHECK(out.uncovered == 0);
    CHECK(out.nodes_inserted == 0);  // refine_strip only
    const Mesh m = mesh_of(out);
    CHECK(sources_over(g, src, b, m, tol) == 0);  // E3: 15c's J2
    const auto f = e1(dem, b, m, tol);            // E1, against the target grid
    CHECK(f.unfiled == 0);
    CHECK(f.on_void == 0);
    CHECK(f.over == 0);
    CHECK(node_findings(dem, m, tol).not_ccw == 0);
    CHECK(delaunay_violations(g, m) == 0);

    // Every new vertex is a source point (its own z) or lies strictly inside
    // a start constraint edge (a strip point, z the target grid's there).
    std::map<std::pair<double, double>, float> given;
    for (std::size_t i = 0; i < src.xy.size(); ++i) given[{src.xy[i].x, src.xy[i].y}] = src.z[i];
    const auto s = shape_findings(dem, b.mesh.vertices(), b.edges, b.masks, m);
    CHECK(s.coincident == 0);
    CHECK(s.unplaced_edge == 0);
    CHECK(s.wrong_mask == 0);
    CHECK(s.broken_chain == 0);
    std::size_t from_sources = 0;
    for (std::size_t i = b.mesh.vertices().size(); i < out.vertices.size(); ++i)
        if (const auto it = given.find({out.vertices[i].x, out.vertices[i].y}); it != given.end()) {
            CHECK(out.z[i] == static_cast<double>(it->second));
            ++from_sources;
        }
    CHECK(from_sources + out.strip_inserted == out.inserted);

    // Without the strip: the 15c call, no strip fields.
    const auto plain = run_points(store, nullptr, b, tol);
    REQUIRE(plain.ok());
    CHECK(plain.strip_points == 0);
    CHECK(plain.strip_inserted == 0);
    CHECK(plain.strip_refused == 0);
}

// ---------------------------------------------------------------------------
// ES10: refusals, and a start edge the strip leaves out
// ---------------------------------------------------------------------------

TEST_CASE("ES10: a strip from another geometry or another mesh is a logic_error", "[edge_strip][ES10]") {
    const auto dem = bump_dem();
    const auto& g = dem.geometry();
    const Begin b = direct(dem, bump_start());
    const auto other = Raster<float>{exact_geometry(10, 9), std::vector<float>(90, 0.0f)};
    const auto foreign = strip_of(other, b);                 // same edges, another geometry
    const auto off_mesh = constraint_check_points(dem, b.mesh.vertices(),
                                                  Edges{{0, 3}, {0, 2}});  // 0-2 is not a constraint
    CheckPoints store{g};
    store.freeze();

    SECTION("refine_strip") {
        REQUIRE_THROWS_MATCHES(run_strip(dem, foreign, b, 1.0), std::logic_error,
                               MessageMatches(ContainsSubstring("refine_strip")));
        REQUIRE_THROWS_MATCHES(run_strip(dem, off_mesh, b, 1.0), std::logic_error,
                               MessageMatches(ContainsSubstring("refine_strip")));
    }
    SECTION("refine_points") {
        REQUIRE_THROWS_MATCHES(run_points(store, &foreign, b, 1.0), std::logic_error,
                               MessageMatches(ContainsSubstring("refine_points")));
        REQUIRE_THROWS_MATCHES(run_points(store, &off_mesh, b, 1.0), std::logic_error,
                               MessageMatches(ContainsSubstring("refine_points")));
    }
    SECTION("refine_strip refuses a tolerance that is not finite and >= 0") {
        const auto strip = strip_of(dem, b);
        for (const double tol : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity()}) {
            CAPTURE(tol);
            const auto out = run_strip(dem, strip, b, tol);
            CHECK(out.status == RefineStatus::InvalidTolerance);
            CHECK_FALSE(out.ok());
        }
    }
    SECTION("a start constraint edge with no strip edge is allowed, and gets no strip point") {
        const auto only_bottom = constraint_check_points(dem, b.mesh.vertices(), Edges{{0, 3}});
        const double tol = 1.0;
        const auto out = run_strip(dem, only_bottom, b, tol);
        REQUIRE(out.ok());
        CHECK(out.strip_points == only_bottom.size());
        const Mesh m = mesh_of(out);
        CHECK(strip_findings(g, ruled_points(dem, b.mesh.vertices(), Edges{{0, 3}}), m, tol).over == 0);
        // Off-node new vertices lie on the bottom side only (row 3.5).
        for (std::size_t i = b.mesh.vertices().size(); i < out.vertices.size(); ++i) {
            const Lat p = lat(g, out.vertices[i]);
            if (!is_node(p)) CHECK(p.row == 3.5);
        }
    }
}

// ---------------------------------------------------------------------------
// ES11: nothing to do, nothing done (E7)
// ---------------------------------------------------------------------------

TEST_CASE("ES11: a start whose strip points are all within tolerance comes back unchanged",
          "[edge_strip][ES11]") {
    const auto g = exact_geometry();
    const std::uint32_t seed = GENERATE(1u, 2u);
    const auto [kind, tol] = GENERATE(std::pair{Terrain::Plane, 1e-6}, std::pair{Terrain::Rough, 1e4});
    CAPTURE(seed, name(kind), tol);
    const auto dem = terrain_dem(g, kind, seed);
    const Begin b = after_refine(dem, jittered(g, seed), tol);
    const auto strip = strip_of(dem, b);
    REQUIRE(strip.size() > 0);
    const auto out = run_strip(dem, strip, b, tol);
    REQUIRE(out.ok());
    CHECK(out.inserted == 0);
    CHECK(out.strip_inserted == 0);
    CHECK(out.nodes_inserted == 0);
    CHECK(out.strip_max_error <= tol);
    REQUIRE(out.vertices.size() == b.mesh.vertices().size());
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        CAPTURE(i);
        CHECK(out.vertices[i].x == b.mesh.vertices()[i].x);
        CHECK(out.vertices[i].y == b.mesh.vertices()[i].y);
        CHECK(out.z[i] == b.z[i]);
        CHECK(out.valid[i] == b.valid[i]);
    }
    CHECK(std::equal(out.triangles.begin(), out.triangles.end(), b.mesh.triangles().begin(),
                     b.mesh.triangles().end()));
    CHECK(out.edges == b.edges);
    CHECK(out.masks == b.masks);
}

// ---------------------------------------------------------------------------
// ES12: extreme scale (UTM coordinates, non-square 10 m x 5 m cells)
// ---------------------------------------------------------------------------

TEST_CASE("ES12: E1, E2 and Delaunay hold at UTM scale, where world points round",
          "[edge_strip][ES12][property]") {
    // refinement_fixtures' geometry: x_min 500000, y_max 7000000. Lattice
    // positions no longer map to the world and back exactly, so this is the
    // one case where the oracles read positions to rounding (1e-9 cells).
    const terrain::raster::RasterGeometry g{500000.0, 7000000.0, 10.0, 5.0, kN, kN};
    const double tol = GENERATE(0.0, 0.5);
    CAPTURE(tol);
    const auto dem = terrain_dem(g, Terrain::Smooth, 1);
    const Begin b = after_refine(dem, jittered(g, 1), tol);
    const auto strip = strip_of(dem, b);
    const auto out = run_strip(dem, strip, b, tol);
    projected_oracles(dem, strip, b, out, tol);
}

// ---------------------------------------------------------------------------
// ES13: ends at ulp distance from grid lines (D2 steps 5-6, I1-I3, at the loop)
// ---------------------------------------------------------------------------

namespace {

double ulps(double x, int k) {
    const double to = k > 0 ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    for (int i = 0; i < std::abs(k); ++i) x = std::nextafter(x, to);
    return x;
}

// One triangle (a, b, apex) with a-b constrained, over a 24 x 16 node DEM.
// Asserts what the loop owes such a start: it ends; no two output vertices
// share a position; every triangle is strictly counter-clockwise; E2 and
// Delaunay hold; and every point the strip oracle finds over the tolerance is
// within 1e-9 cells of an end and accounted for as refused.
void ulp_run(Point2 a, Point2 b, Point2 apex, double tol) {
    const auto g = exact_geometry(24, 16);
    std::vector<float> v(24 * 16);
    for (std::size_t r = 0; r < 16; ++r)
        for (std::size_t c = 0; c < 24; ++c)
            v[r * 24 + c] = static_cast<float>(1000 + 7 * c + 13 * r + 4 * ((5 * r + 3 * c) % 11));
    const Raster<float> dem{g, std::move(v)};
    Start s;
    s.mesh = IndexedMesh2{{a, b, apex}, {{0, 1, 2}}, {0}};
    s.edges = {{0, 1}};
    s.masks = {1};
    const Begin begin = direct(dem, s);
    const auto strip = strip_of(dem, begin);
    const auto out = run_strip(dem, strip, begin, tol);
    REQUIRE(out.ok());
    const Mesh m = mesh_of(out);
    const auto sh = shape_findings(dem, begin.mesh.vertices(), begin.edges, begin.masks, m);
    CHECK(sh.coincident == 0);
    CHECK(sh.stray == 0);
    CHECK(sh.broken_chain == 0);
    const auto n = node_findings(dem, m, tol);
    CHECK(n.not_ccw == 0);
    CHECK(n.over == 0);
    CHECK(delaunay_violations(g, m) == 0);
    const auto f = e1(dem, begin, m, tol);
    CHECK(f.unfiled == 0);
    if (f.over > 0) CHECK(out.strip_refused > 0);
    const Lat la = lat(g, a), lb = lat(g, b);
    for (const Lat q : f.over_at) {
        CAPTURE(q.col, q.row);
        CHECK(std::min(std::hypot(q.col - la.col, q.row - la.row), std::hypot(q.col - lb.col, q.row - lb.row))
              <= 1e-9);
    }
}

}  // namespace

TEST_CASE("ES13: ends within ulps of grid lines: the loop ends, with no coincident vertex and no flat triangle",
          "[edge_strip][ES13]") {
    const auto g = exact_geometry(24, 16);
    const double tol = GENERATE(0.0, 0.5);
    CAPTURE(tol);
    SECTION("P1 one ulp past column 10 (the crossing at t == 1 is not a point)") {
        ulp_run(world(g, 0.8932792255671602, 1.5), world(g, std::nextafter(10.0, 11.0), 3.25), world(g, 6.0, 0.5),
                tol);
    }
    SECTION("P0 one ulp before column 1 (the crossing a hair from P0 is a point)") {
        ulp_run(world(g, std::nextafter(1.0, 0.0), 1.5), world(g, 4.0, 3.25), world(g, 3.0, 1.0), tol);
    }
    SECTION("a seeded sweep: ends 1 to 4 ulps either side of grid lines, apexes either side of the wedge") {
        for (const std::uint32_t seed : {11u, 12u, 13u}) {
            std::mt19937 gen{seed};
            for (int k = 0; k < 40; ++k) {
                auto near = [&](double line) { return ulps(line, static_cast<int>(gen() % 9u) - 4); };
                const double ac = near(static_cast<double>(1 + gen() % 8u)), ar = near(static_cast<double>(2 + gen() % 4u));
                const double bc = near(static_cast<double>(12 + gen() % 8u)), br = near(static_cast<double>(8 + gen() % 6u));
                // The apex: left of a -> b in (col, -row), from just off the
                // line to well clear of it, dyadic.
                const double tc = static_cast<double>(1 + gen() % 255u) / 256.0;
                const double lift = static_cast<double>(1 + gen() % 64u) / 16.0;
                const double mc = ac + tc * (bc - ac), mr = ar + tc * (br - ar) - lift;
                const Point2 a = world(g, ac, ar), b = world(g, bc, br);
                const Point2 apex = world(g, std::round(mc * 64) / 64, std::clamp(std::round(mr * 64) / 64, 0.0, 15.0));
                const Lat la{ac, -ar}, lb{bc, -br}, lp{std::round(mc * 64) / 64, -std::clamp(std::round(mr * 64) / 64, 0.0, 15.0)};
                if (DefaultKernel::orient2d(Point2{la.col, la.row}, Point2{lb.col, lb.row}, Point2{lp.col, lp.row})
                    != Orientation::CounterClockwise)
                    continue;
                CAPTURE(seed, k, ac, ar, bc, br);
                ulp_run(a, b, apex, tol);
            }
        }
    }
}
