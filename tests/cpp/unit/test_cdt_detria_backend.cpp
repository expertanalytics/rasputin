// Unit tests for the detria backend: src/cdt/detria_backend.cpp behind
// terrain/cdt/detria_backend.hpp, plus the result channel in
// terrain/cdt/result.hpp.
//
// INVARIANT-CRITICAL. A mutation round applies to this file. It is where the
// increment's topology decisions live: which detria call each ChainRole gets,
// which triangles come back, and which detria error becomes which CdtStatus.
//
// ONE NAMED TEST PER ROW of the degeneracy table in docs/increments/04-cdt.md.
// Every row of that table was MEASURED against lib/detria/detria.hpp at the
// pinned SHA rather than inferred, and these tests are the mechanism a sceptical
// reader re-runs. For the Ok rows the exact interior triangle count is asserted,
// because "it succeeded" is satisfied by a backend that returns the convex hull,
// the hole triangles, or nothing at all.
//
// NO TEST HERE ASSERTS MESSAGE TEXT. `message` carries detria's own diagnosis --
// "Point 4 is exactly on a constrained edge (between points 0 and 1)" -- which
// is a third party's prose and the most useful thing at the point of failure
// precisely because it names indices our enum has no field for. Pinning its
// wording would fail the suite on an upstream improvement to it. What is pinned
// is the empty/non-empty split, which is ours.
//
// Every fixture is hand-made NON-CROSSING where a mesh is expected. Between
// increments 4 and 5 the wrapper is correct on all input and useful only on
// input that happens not to cross; that is a stated risk, not an oversight
// here.

#include <catch2/catch_test_macros.hpp>

#include <cdt_cases.hpp>
#include <mesh_queries.hpp>

#include <terrain/cdt/detria_backend.hpp>
#include <terrain/cdt/result.hpp>
#include <terrain/cdt/triangulate.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <set>
#include <string_view>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::Pslg;
using terrain::PslgBuilder;
using terrain::TriangleIndices;
using terrain::cdt::CdtOptions;
using terrain::cdt::CdtOutcome;
using terrain::cdt::CdtStatus;
using terrain::cdt::DetriaBackend;
using terrain::cdt::describe;
using terrain::pred::DefaultKernel;
using terrain::test::build_fixture;
using terrain::test::edge_uses;
using terrain::test::undirected_key;

namespace {

[[nodiscard]] CdtOutcome run(PslgBuilder b, const CdtOptions& options = {}) {
    const Pslg pslg = build_fixture<DefaultKernel>(std::move(b));
    return terrain::cdt::triangulate<DetriaBackend>(pslg, options);
}

// The set of triangles, each rotated to start at its lowest index, so two
// meshes can be compared without depending on triangle order or on which
// vertex a triangle happens to be listed from. Winding is preserved by the
// rotation, so this is not a laxer comparison than it looks.
[[nodiscard]] std::set<std::array<std::uint32_t, 3>> triangle_set(const IndexedMesh2& mesh) {
    std::set<std::array<std::uint32_t, 3>> out;
    for (const TriangleIndices& t : mesh.triangles()) {
        const std::size_t lowest =
            (t[0] <= t[1] && t[0] <= t[2]) ? 0 : (t[1] <= t[2] ? 1 : 2);
        out.insert({t[lowest], t[(lowest + 1) % 3], t[(lowest + 2) % 3]});
    }
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// The result channel
// ---------------------------------------------------------------------------

// NotRun is a failure value and it is the DEFAULT member initializer, for the
// same reason Chain::role defaults to Breakline: a default-constructed outcome
// that reaches a reader must not be able to claim success.
TEST_CASE("a default-constructed outcome is NotRun and not ok", "[cdt][backend][result]") {
    const CdtOutcome out;

    CHECK(out.status == CdtStatus::NotRun);
    CHECK_FALSE(out.ok());
    CHECK(out.message.empty());
    CHECK(out.mesh.empty());
}

// describe() renders every enumerator, INCLUDING the two no reachable input
// produces. MalformedInput fires only if the Pslg did not come from the builder
// or our wrapper mis-built a span, and BackendFailure only on a detria bug;
// both are self-checks that must still be able to say what they are when they
// appear in a log.
TEST_CASE("describe renders every status, distinctly", "[cdt][backend][result]") {
    constexpr std::array all = {
        CdtStatus::Ok,          CdtStatus::NotRun,         CdtStatus::NotNoded,
        CdtStatus::DegenerateGeometry, CdtStatus::InvalidTopology, CdtStatus::MalformedInput,
        CdtStatus::BackendFailure,
    };

    std::set<std::string_view> rendered;
    for (const CdtStatus s : all) {
        const std::string_view text = describe(s);
        INFO("status " << static_cast<int>(s));
        CHECK_FALSE(text.empty());
        rendered.insert(text);
    }
    CHECK(rendered.size() == all.size());
}

// Ok says nothing; failure says something. getErrorMessage() builds through
// <sstream>, so it is called only on the failure path -- and detria's own
// success string ("The triangulation was successful") is deliberately NOT
// forwarded. This assertion is what kills the mutant that formats a message on
// the happy path, which a profiler would otherwise be the only thing to notice.
TEST_CASE("Ok implies an empty message and failure implies a non-empty one",
          "[cdt][backend][result]") {
    const CdtOutcome good = run(terrain::test::polygon_with_holes_domain());
    REQUIRE(good.ok());
    CHECK(good.message.empty());

    const CdtOutcome bad = run(terrain::test::crossing_breaklines_domain());
    REQUIRE_FALSE(bad.ok());
    CHECK_FALSE(bad.message.empty());
}

// No outcome carries a mesh alongside a failure, and none carries Ok with an
// empty mesh. The second half is not hypothetical: a Triangulation with points
// but no addOutline call SUCCEEDS and yields zero interior triangles, so
// "forgot to add the outline" is a silent success without this. The generic
// entry point asserts both in debug; in a release build this test is the only
// guard.
TEST_CASE("status and mesh never disagree", "[cdt][backend][result]") {
    const CdtOutcome bad = run(terrain::test::hole_outside_outline_domain());
    REQUIRE_FALSE(bad.ok());
    CHECK(bad.mesh.empty());
    // empty() asks about TRIANGLES, so it is satisfied by a mesh that carries
    // the whole vertex array and no triangles. A failure returns no mesh at
    // all, and this is the assertion that says so -- without it, a wrapper that
    // copies the vertices before checking the status passes.
    CHECK(bad.mesh.vertices().empty());

    const CdtOutcome good = run(terrain::test::square_domain());
    REQUIRE(good.ok());
    CHECK_FALSE(good.mesh.empty());
}

// ---------------------------------------------------------------------------
// The degeneracy table, row by row -- the Ok rows
// ---------------------------------------------------------------------------

TEST_CASE("the minimal valid input is one Outer triangle", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::single_triangle_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 1);
}

TEST_CASE("a square is two triangles", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::square_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 2);  // Euler: 2*4 - 4 - 2
}

TEST_CASE("a square with a disjoint hole is eight triangles", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::disjoint_hole_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 8);  // Euler: 2*8 - 8 - 2 + 2
}

// The goal line: a real CDT over a polygon with holes, which is what increment
// 1's plan was written against and what nothing in the tree has produced yet.
//
// This is ALSO THE MULTI-RING MULTI-HOLE ASAN FIXTURE. READ BEFORE SIMPLIFYING:
// detria's addOutline, addHole and setPoints store a span and DO NOT COPY, so a
// wrapper that builds any per-ring buffer of its own hands detria a dangling
// span at triangulate(). With a single ring that dangling span usually still
// points at live memory and the bug passes; it takes several rings and several
// holes for the allocator to reuse the storage. Under the asan+ubsan Debug job
// this test passes by not crashing, and that is a second thing it asserts on
// top of the count.
TEST_CASE("an outer ring with two holes and two breaklines triangulates",
          "[cdt][backend][goal]") {
    const CdtOutcome out = run(terrain::test::polygon_with_holes_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    // n = 4 + 4 + 4 + 3 + 2 referenced, b = 12 on rings, h = 2:
    // 2*17 - 12 - 2 + 4.
    CHECK(out.mesh.triangle_count() == 24);
    CHECK(out.mesh.vertices().size() == 17);
}

// A hole may touch its outer ring, but ONLY IF THE TOUCH SHARES AN INDEX. Two
// coincident vertices are DuplicatePointsFound; see the companion test below.
//
// The rings are no longer disjoint here, so THE EULER HELPER DOES NOT APPLY:
// with n = 6, b = 6, h = 1 the formula predicts 6 and the measured answer is 5,
// because the shared vertex is one vertex on two rings and the domain is not a
// disk with a hole. The count below is explicit for that reason. Do not feed
// this fixture to the general helper.
TEST_CASE("a hole touching its outer ring at a shared index is accepted",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::corner_touching_hole_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 5);
}

TEST_CASE("collinear consecutive vertices within a ring are accepted",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::collinear_ring_vertices_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 3);  // Euler: 2*5 - 5 - 2
}

TEST_CASE("an island in a lake is supported", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::island_in_lake_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    // Not an Euler case either: two outer rings, one nested inside the hole.
    // 8 triangles between the square and the lake, plus the island itself.
    CHECK(out.mesh.triangle_count() == 9);
}

// A vertex outside the outer ring stays in the vertex array -- index k means
// the same point coming out as going in, which is what lets Python hold a
// parallel attribute array -- and is referenced by no interior triangle. The
// mesh explicitly does NOT guarantee that every vertex is referenced, and the
// Euler check counts referenced vertices for exactly this reason.
TEST_CASE("a vertex outside the outer ring is kept but unreferenced",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::outside_vertex_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 2);
    REQUIRE(out.mesh.vertices().size() == 5);
    CHECK(out.mesh.vertices()[4] == Point2{9.0, 9.0});
    for (const TriangleIndices& t : out.mesh.triangles()) {
        CHECK(t[0] != 4u);
        CHECK(t[1] != 4u);
        CHECK(t[2] != 4u);
    }
}

// There is no area threshold anywhere in this project -- increment 2 settled
// that there are no tolerances -- so a sliver is an ordinary domain. "No
// near-zero-area triangle" is not an invariant this project can state; the
// strongest true statement is that every triangle is counterclockwise under an
// exact predicate, which the property suite asserts.
TEST_CASE("a sliver outer ring is an ordinary one-triangle domain",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::sliver_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 1);
}

// Breaklines reach the backend EDGE BY EDGE, never as a polyline, and
// addPolylineAutoDetectType is prohibited: we know each chain's role from its
// ChainRole and letting the library re-derive it geometrically is a
// silent-wrong-answer channel. This is the concrete failure -- a closed-loop
// breakline (a contour, a ring road) auto-detected as a hole would CARVE A VOID
// out of the domain. Fed as individual constrained edges the full interior
// survives, so the count is what tells the two apart.
TEST_CASE("a closed-loop breakline does not carve a hole", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::closed_loop_breakline_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 10);  // 8 were it a hole
}

// ...but only when the loop closes on a SHARED INDEX. The same ring road with a
// coincident closing vertex is the duplicate-coordinate row of the table all
// over again, and the pair of tests is what stops the one above being read as
// "closed-loop breaklines are fine".
TEST_CASE("a closed-loop breakline spelled with a coincident vertex is DegenerateGeometry",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::coincident_closed_breakline_domain());

    CHECK(out.status == CdtStatus::DegenerateGeometry);
}

// ---------------------------------------------------------------------------
// The degeneracy table -- the failing rows, and the status each maps to
// ---------------------------------------------------------------------------
//
// The mapping is grouped by WHAT THE CALLER SHOULD DO, not by where in detria
// the failure arose, which is why two detria enumerators can land on one status
// below. Nothing is lost by the grouping: the offending indices are in
// `message`.

// A Pslg permits coincident vertices; detria's duplicate scan runs over the
// WHOLE point array, unreferenced vertices included. This is the most likely
// failure on real data -- two features digitised to the same corner -- and it
// is not fixable here: a coordinate-keyed dedup needs the snap grid increment 5
// owns.
TEST_CASE("coincident vertices are DegenerateGeometry even when unreferenced",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::duplicate_unreferenced_vertex_domain());

    CHECK(out.status == CdtStatus::DegenerateGeometry);
}

TEST_CASE("a hole touching its outer ring at coincident vertices is DegenerateGeometry",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::coincident_touch_hole_domain());

    // Same geometry as the accepted corner-touching fixture, spelled with a
    // second vertex instead of a shared index. How you spell the touch decides
    // whether it works.
    CHECK(out.status == CdtStatus::DegenerateGeometry);
}

// {0, 0, 1, 2, 3}: a valid Pslg -- increment 2's degeneracy 4 accepts a
// repeated index -- that detria refuses. Passed through as a status rather than
// pre-checked.
TEST_CASE("a repeated consecutive index is DegenerateGeometry", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::repeated_index_ring_domain());

    CHECK(out.status == CdtStatus::DegenerateGeometry);
}

// A T-junction: a vertex lying exactly ON a constraint edge it is not part of.
// parallel_refinement.md assigns this to the noder, so until increment 5 lands
// it is a loud failure and never a mesh that ignores it.
TEST_CASE("a foreign vertex on a constraint edge is NotNoded", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::foreign_vertex_on_constraint_domain());

    CHECK(out.status == CdtStatus::NotNoded);
    CHECK(out.mesh.empty());
}

// Crossing constraints. The load-bearing part is the second assertion: the CDT
// NEVER SILENTLY PRODUCES A MESH THAT IGNORES A CROSSING. It does not repair,
// does not insert intersection points, and returns no mesh.
TEST_CASE("crossing constraint edges are NotNoded and yield no mesh",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::crossing_breaklines_domain());

    CHECK(out.status == CdtStatus::NotNoded);
    CHECK(out.mesh.empty());
    CHECK_FALSE(describe(CdtStatus::NotNoded).empty());
}

TEST_CASE("a hole sharing a whole edge with its outer ring is InvalidTopology",
          "[cdt][backend][degeneracy]") {
    // A hole that shares an edge with its boundary is a notch in the boundary,
    // not a hole.
    const CdtOutcome out = run(terrain::test::hole_sharing_edge_domain());

    CHECK(out.status == CdtStatus::InvalidTopology);
}

// The nesting check increment 3 deferred, landing here as a status instead of a
// diagnostic. detria reports this one as StackedPolylines, not as
// HoleNotInsideOutline -- that enumerator fires on the different arm below --
// and both map to InvalidTopology, which is why the grouping is by action.
TEST_CASE("a hole outside every outline is InvalidTopology", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::hole_outside_outline_domain());

    CHECK(out.status == CdtStatus::InvalidTopology);
}

// The other arm: a hole that CONTAINS the outer ring, so the outermost polyline
// is a hole. Same status, different detria enumerator, and the pair of them is
// what the mapping's grouping has to survive.
TEST_CASE("a hole containing the outer ring is InvalidTopology", "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::hole_containing_outline_domain());

    CHECK(out.status == CdtStatus::InvalidTopology);
}

// A constraint edge OUTSIDE the domain is silently dropped, and that is
// accepted behaviour rather than a bug: a Pslg promises no disjointness and no
// nesting, so a breakline may legitimately lie outside the outer ring or inside
// a hole, and no interior triangle can carry it.
//
// This is the input that bounds the constraint-preservation property. That
// property reads "every constraint edge appears as an edge of one or two output
// triangles", which is FALSE here -- so it is asserted over the property
// suite's generated families, whose features are all in-domain by construction,
// and the exception is pinned here instead of being left as an unstated
// precondition.
TEST_CASE("a constraint edge outside the domain is dropped, not honoured",
          "[cdt][backend][degeneracy]") {
    const CdtOutcome out = run(terrain::test::out_of_domain_breaklines_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    CHECK(out.mesh.triangle_count() == 8);  // the hole domain; the breaklines add nothing

    const auto uses = edge_uses(out.mesh);
    CHECK(uses.find(undirected_key(8, 9)) == uses.end());    // outside the outer ring
    CHECK(uses.find(undirected_key(10, 11)) == uses.end());  // inside the hole
    // Both breakline vertex pairs are still in the vertex array, in order:
    // index identity does not depend on a vertex being used.
    CHECK(out.mesh.vertices().size() == 12);
}

// ---------------------------------------------------------------------------
// The characterisation test: detria closes its own polylines
// ---------------------------------------------------------------------------
//
// READ THIS BEFORE DELETING IT AS "TESTING A THIRD-PARTY LIBRARY". It is a
// CHARACTERISATION TEST pinning behaviour we depend on at a pinned SHA, which
// is a carve-out from the rule that vendored libraries are trusted to their own
// suites.
//
// The dependency: detria's createConstrainedEdges seeds prevVertexIdx with
// polyline.back() when the polyline is open (detria.hpp:3543-3560), so an open
// index span is closed by detria itself. That is why Pslg::indices_of(c) goes
// straight into addOutline with NO SCRATCH BUFFER, no ClosedChainSpans type and
// no closing index materialised anywhere in the project -- and why
// closed_index_buffer_size was deleted. It is a property of the vendored
// SOURCE, not of a documented API contract.
//
// IF THIS TEST FAILS AFTER A RE-PIN: the wrapper must materialise closing
// indices again. The design is one scratch buffer for the whole triangulation,
// sized by a restored closed_index_buffer_size, with sub-spans handed to
// addOutline/addHole -- a five-line change, and knowing where it goes is the
// whole value of this test.
//
// It cannot be written as the design states it -- "the same square given open
// and closed produces identical output" -- because A PSLG CANNOT CARRY A CLOSED
// RING. Validator stage 4 rejects a stored closure (guarantee 5), and the
// closed spelling has no other route to the backend: test targets cannot
// include detria.hpp. What is asserted instead is the observable consequence:
// the closing edge is present, on a non-convex ring where it is a wall of the
// domain rather than a hull edge, and the interior is the closed polygon's.
TEST_CASE("an open ring span is closed by the backend", "[cdt][backend][characterisation]") {
    const CdtOutcome out = run(terrain::test::l_shaped_domain());

    REQUIRE(out.status == CdtStatus::Ok);
    // The L-shaped hexagon, closed: 2*6 - 6 - 2. An unclosed polyline bounds no
    // region at all, and detria's interior classification would not produce
    // this.
    REQUIRE(out.mesh.triangle_count() == 4);

    // The closing edge (v5, v0) is a boundary edge -- used by exactly one
    // triangle -- and it carries its mask bit there.
    const auto uses = edge_uses(out.mesh);
    const auto it = uses.find(undirected_key(5, 0));
    REQUIRE(it != uses.end());
    REQUIRE(it->second.size() == 1);
    CHECK(out.mesh.is_constrained(it->second.front().triangle, it->second.front().slot));
}

// ---------------------------------------------------------------------------
// The Delaunay knob
// ---------------------------------------------------------------------------

// CdtOptions::delaunay is forwarded to triangulate(bool). With it off, detria
// still produces a valid CONSTRAINED triangulation -- same domain, same
// triangle count by Euler -- so the count alone cannot tell the two apart. The
// property suite is where the difference shows, as a visibility-Delaunay
// violation; here we pin that the knob is wired through at all and that
// switching it does change the mesh.
TEST_CASE("delaunay=false still triangulates the same domain", "[cdt][backend][options]") {
    const CdtOutcome on = run(terrain::test::polygon_with_holes_domain(), CdtOptions{true});
    const CdtOutcome off = run(terrain::test::polygon_with_holes_domain(), CdtOptions{false});

    REQUIRE(on.status == CdtStatus::Ok);
    REQUIRE(off.status == CdtStatus::Ok);
    CHECK(on.mesh.triangle_count() == off.mesh.triangle_count());
    CHECK(triangle_set(on.mesh) != triangle_set(off.mesh));
}

// ---------------------------------------------------------------------------
// What the backend promises about the mesh it returns
// ---------------------------------------------------------------------------

// Obligation 1 of the CdtBackend concept, for this backend in its strongest
// form: DetriaBackend appends nothing, so the arrays are equal, not merely
// prefixed. Index k means the same point coming out as going in.
TEST_CASE("the mesh vertex array is the Pslg's, element for element",
          "[cdt][backend][identity]") {
    const Pslg pslg = build_fixture<DefaultKernel>(terrain::test::polygon_with_holes_domain());
    const CdtOutcome out = terrain::cdt::triangulate<DetriaBackend>(pslg);

    REQUIRE(out.ok());
    REQUIRE(out.mesh.vertices().size() == pslg.vertices().size());
    for (std::size_t i = 0; i < pslg.vertices().size(); ++i) {
        INFO("vertex " << i);
        CHECK(out.mesh.vertices()[i] == pslg.vertices()[i]);
    }
}

// Obligation 2: only in-domain triangles. forEachTriangle filters to
// TriangleLocation::Interior, which is what makes the counts above meaningful;
// forEachTriangleOfEveryLocation would add the hole and convex-hull triangles
// and this assertion is what kills that substitution on a fixture whose hull is
// strictly larger than its domain.
TEST_CASE("no hole triangle and no convex-hull triangle is returned",
          "[cdt][backend][identity]") {
    const CdtOutcome out = run(terrain::test::disjoint_hole_domain());

    REQUIRE(out.ok());
    // The hole is a square on indices 4..7; its own two triangles would appear
    // if the location filter were dropped, as would nothing outside the square
    // -- which is why the L-shaped domain is checked too.
    CHECK(out.mesh.triangle_count() == 8);

    const CdtOutcome l = run(terrain::test::l_shaped_domain());
    REQUIRE(l.ok());
    CHECK(l.mesh.triangle_count() == 4);  // 5 over every location: the hull adds one
}
