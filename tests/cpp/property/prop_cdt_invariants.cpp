// Property tests for the mesh a CDT produces.
//
// INVARIANT-CRITICAL, alongside the backend and mask suites. The unit suites
// pin named configurations; this file asserts the laws that must hold for EVERY
// mesh the backend returns, over generated rectangle-with-holes families.
//
// THE PROPERTY SUITE RUNS UNDER DefaultKernel ONLY, AND FastKernel IS FORBIDDEN
// HERE. This is a ruling, not a default. The visibility-Delaunay property
// evaluates incircle as an ORACLE against a mesh detria built with its own
// robust predicates; an oracle less exact than the thing it judges reports
// false failures on precisely the cocircular inputs the property is about.
// FastKernel is header-only and therefore the tempting choice for a suite that
// is slow to link -- that is why this paragraph is here.
//
// No rapidcheck. Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64 gives reproducible input without a second framework; a
// failure prints its seed as the generator index and rerunning that section
// reproduces it exactly.
//
// The generator is cdt_cases.hpp's, not increment 3's valid_chain_specs: a CDT
// property over crossing constraints would assert almost nothing, because a
// crossing is a NotNoded failure with no mesh. See that header for the full
// reason. That reason SURVIVES increment 5c and is worth saying, because the
// noder looks like it removes it: node() would turn a generated crossing into a
// mesh, but it would also change the vertex and chain counts the Euler property
// takes from DomainShape, so the property would be asserting Euler's formula
// against a domain nobody declared. The generator's disjointness is now about
// the ORACLE rather than about the backend's tolerance.
//
// AMENDED AT INCREMENT 5c. Every case is noded before it is triangulated, at
// cdt_cases.hpp's kFixtureSpacing. Measured over all 24 seeds on both arms
// before this was written: every generated domain nodes Ok with its vertex
// count and its chain structure unchanged, so every count below is still the
// count the generator declared.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <cdt_cases.hpp>
#include <mesh_queries.hpp>

#include <terrain/cdt/constrained_edges.hpp>
#include <terrain/cdt/detria_backend.hpp>
#include <terrain/cdt/result.hpp>
#include <terrain/cdt/triangulate.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <cstddef>
#include <cstdint>
#include <format>
#include <map>
#include <random>
#include <set>
#include <span>
#include <utility>
#include <vector>

using terrain::Chain;
using terrain::ChainRole;
using terrain::IndexedMesh2;
using terrain::IndexedRing;
using terrain::NodedPslg;
using terrain::Point2;
using terrain::PointInRing;
using terrain::Segment2;
using terrain::TriangleIndices;
using terrain::cdt::CdtOptions;
using terrain::cdt::CdtOutcome;
using terrain::cdt::CdtStatus;
using terrain::cdt::ConstraintEdgeSet;
using terrain::cdt::DetriaBackend;
using terrain::is_closed;
using terrain::point_in_ring;
using terrain::pred::DefaultKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::test::DomainShape;
using terrain::test::EdgeUse;
using terrain::test::GeneratedDomain;
using terrain::test::apex;
using terrain::test::blocks_visibility;
using terrain::test::edge_uses;
using terrain::test::euler_triangles;
using terrain::test::generated_domain;
using terrain::test::node_fixture;
using terrain::test::undirected_key;

namespace {

constexpr int seed_count = 24;

[[nodiscard]] std::mt19937_64 seeded(int seed) {
    return std::mt19937_64{0x0CD7'0000ULL + static_cast<std::uint64_t>(seed)};
}

// A generated case, triangulated. Every property goes through this, so "the
// generator produced input the backend rejects" fails at the generator rather
// than as a confusing downstream assertion -- and a rejection IS a generator
// bug here: the layout guarantees disjoint features, which is precisely the
// precondition a mesh existing at all depends on before increment 5.
struct Case {
    NodedPslg pslg;
    DomainShape shape;
    CdtOutcome outcome;
};

[[nodiscard]] Case generated(int seed, const CdtOptions& options = {}) {
    std::mt19937_64 rng = seeded(seed);
    GeneratedDomain domain = generated_domain(rng);
    NodedPslg pslg = node_fixture<DefaultKernel>(std::move(domain.builder));
    CdtOutcome outcome = terrain::cdt::triangulate<DetriaBackend>(pslg, options);
    INFO("status " << static_cast<int>(outcome.status) << ": " << outcome.message);
    REQUIRE(outcome.ok());
    return Case{std::move(pslg), domain.shape, std::move(outcome)};
}

// The constraint edges of a Pslg, derived here by hand rather than through
// ConstraintEdgeSet: the mask property must not be checked against the same
// code that computes the mask, or it asserts only that a function equals
// itself.
[[nodiscard]] std::set<std::uint64_t> input_edges(const NodedPslg& p) {
    std::set<std::uint64_t> out;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        const std::span<const std::uint32_t> idx = p.indices_of(c);
        for (std::size_t k = 0; k < p.edge_count(c); ++k) {
            out.insert(undirected_key(idx[k], idx[(k + 1) % idx.size()]));
        }
    }
    return out;
}

[[nodiscard]] Point2 centroid(const IndexedMesh2& mesh, const TriangleIndices& t) {
    return (mesh.vertices()[t[0]] + mesh.vertices()[t[1]] + mesh.vertices()[t[2]]) / 3.0;
}

}  // namespace

// ---------------------------------------------------------------------------
// Constraint preservation
// ---------------------------------------------------------------------------

// Every constraint edge of the input appears as an edge of one or two output
// triangles, and carries its mask bit on every triangle that has it.
//
// This works only because detria NEVER SPLITS A CONSTRAINT EDGE -- it inserts
// no Steiner points, and a vertex lying on a constraint is a hard
// PointOnConstrainedEdge failure rather than a split. So each input edge is
// EXACTLY ONE output edge. A backend that splits edges (poly2tri does) would
// need the mask computed by walking the split chain, and would fail this
// property rather than silently producing a wrong mask; that assumption belongs
// to DetriaBackend and is stated in constrained_edges.hpp.
//
// The generator keeps every constraint inside the domain. An input edge lying
// INSIDE A HOLE or OUTSIDE THE OUTER RING is legal in a Pslg, triangulates Ok,
// and appears in no output triangle at all -- so this property is false as
// stated for such an input, and the fixture that pins that is in the backend
// suite.
TEST_CASE("every constraint edge survives into the mesh", "[cdt][property][constraints]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);
    const auto uses = edge_uses(c.outcome.mesh);

    for (const std::uint64_t key : input_edges(c.pslg)) {
        INFO(std::format("constraint edge key {}", key));
        const auto it = uses.find(key);
        REQUIRE(it != uses.end());
        REQUIRE_FALSE(it->second.empty());
        REQUIRE(it->second.size() <= 2);
        for (const EdgeUse& u : it->second) {
            REQUIRE(c.outcome.mesh.is_constrained(u.triangle, u.slot));
        }
    }
}

// And the converse, which is the half that catches a mask set too widely: a
// triangle edge flagged constrained must be an edge of the input. Without this
// the property above is satisfied by a backend that returns 7 for every
// triangle.
TEST_CASE("no unconstrained edge is flagged", "[cdt][property][constraints]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);
    const std::set<std::uint64_t> expected = input_edges(c.pslg);

    for (std::size_t t = 0; t < c.outcome.mesh.triangle_count(); ++t) {
        const TriangleIndices& tri = c.outcome.mesh.triangles()[t];
        for (std::size_t e = 0; e < 3; ++e) {
            if (!c.outcome.mesh.is_constrained(t, e)) continue;
            INFO(std::format("triangle {} slot {}", t, e));
            REQUIRE(expected.contains(undirected_key(tri[e], tri[(e + 1) % 3])));
        }
    }
}

// ConstraintEdgeSet against an independent derivation of the same set. This is
// not the circular check the comment on input_edges warns about: there the
// question is whether the MESH's mask matches the input, and using the set to
// answer it would compare a function with itself. Here the set itself is the
// subject, and input_edges is the oracle.
TEST_CASE("the constraint edge set is exactly the input's edges",
          "[cdt][property][constraints]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);
    const std::set<std::uint64_t> expected = input_edges(c.pslg);
    const ConstraintEdgeSet edges{c.pslg};

    REQUIRE(edges.size() == expected.size());
    for (const std::uint64_t key : expected) {
        const auto a = static_cast<std::uint32_t>(key >> 32);
        const auto b = static_cast<std::uint32_t>(key & 0xFFFF'FFFFu);
        REQUIRE(edges.contains(a, b));
        REQUIRE(edges.contains(b, a));
    }
}

// ---------------------------------------------------------------------------
// Visibility-Delaunay
// ---------------------------------------------------------------------------

// THE DELAUNAY INVARIANT OF A CDT IS THE VISIBILITY FORM, NOT PLAIN INCIRCLE.
// For every interior non-constrained edge (a, b) with apexes c and d, either
// incircle(a, b, c, d) is not Inside, OR the open segment cd crosses at least
// one constraint edge -- the flip that would fix the incircle violation is the
// flip a constraint forbids.
//
// Plain incircle is FALSE for every CDT that has a constraint, which increment
// 1 recorded and testing.md still states wrongly. A suite asserting it would
// fail on the first hole.
TEST_CASE("every flippable interior edge is locally Delaunay",
          "[cdt][property][delaunay]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);
    const IndexedMesh2& mesh = c.outcome.mesh;
    const std::set<std::uint64_t> constraints = input_edges(c.pslg);

    std::vector<Segment2> constraint_segments;
    for (std::size_t ch = 0; ch < c.pslg.chains().size(); ++ch) {
        for (std::size_t k = 0; k < c.pslg.edge_count(ch); ++k) {
            constraint_segments.push_back(c.pslg.edge(ch, k));
        }
    }

    for (const auto& [key, uses] : edge_uses(mesh)) {
        if (uses.size() != 2 || constraints.contains(key)) continue;

        const std::uint32_t a = mesh.triangles()[uses[0].triangle][uses[0].slot];
        const std::uint32_t b = mesh.triangles()[uses[0].triangle][(uses[0].slot + 1) % 3];
        const Point2& pa = mesh.vertices()[a];
        const Point2& pb = mesh.vertices()[b];
        const Point2& pc = mesh.vertices()[apex(mesh, uses[0])];
        const Point2& pd = mesh.vertices()[apex(mesh, uses[1])];

        if (DefaultKernel::incircle(pa, pb, pc, pd) != Incircle::Inside) continue;

        INFO(std::format("edge key {} apexes ({}, {}) and ({}, {})", key, pc.x, pc.y, pd.x, pd.y));
        bool blocked = false;
        for (const Segment2& s : constraint_segments) {
            if (blocks_visibility<DefaultKernel>(Segment2{pc, pd}, s)) {
                blocked = true;
                break;
            }
        }
        REQUIRE(blocked);
    }
}

// The knob is what proves the property above is EARNED rather than vacuous.
// With delaunay=false detria produces a valid constrained triangulation that
// violates it, so the mutant "pass false unconditionally" dies here -- and if
// this test ever stops finding a violation, the test above has stopped saying
// anything.
TEST_CASE("delaunay=false produces a triangulation that is not Delaunay",
          "[cdt][property][delaunay]") {
    const Case c = generated(0, CdtOptions{false});
    const IndexedMesh2& mesh = c.outcome.mesh;
    const std::set<std::uint64_t> constraints = input_edges(c.pslg);

    std::size_t violations = 0;
    for (const auto& [key, uses] : edge_uses(mesh)) {
        if (uses.size() != 2 || constraints.contains(key)) continue;
        const std::uint32_t a = mesh.triangles()[uses[0].triangle][uses[0].slot];
        const std::uint32_t b = mesh.triangles()[uses[0].triangle][(uses[0].slot + 1) % 3];
        if (DefaultKernel::incircle(mesh.vertices()[a], mesh.vertices()[b],
                                    mesh.vertices()[apex(mesh, uses[0])],
                                    mesh.vertices()[apex(mesh, uses[1])]) == Incircle::Inside) {
            ++violations;
        }
    }
    CHECK(violations > 0);
}

// ---------------------------------------------------------------------------
// Euler
// ---------------------------------------------------------------------------

// triangles == 2n - b - 2 + 2h, with n the number of REFERENCED vertices, b the
// number on any boundary ring, and h the number of hole rings.
//
// Two conditions, both of which a naive test gets wrong. n is not
// vertices().size(): a vertex outside the domain or inside a hole is in the
// array and in no triangle. And THE RINGS MUST BE DISJOINT -- for a hole
// touching its outer ring the formula predicts 6 against a measured 5, because
// the shared vertex is one vertex on two rings. So n, b and h come from the
// fixture's DECLARED structure, and the touching fixture is asserted with an
// explicit count in the backend suite instead of coming through here.
TEST_CASE("the triangle count is Euler's", "[cdt][property][euler]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);

    INFO(std::format("n {} b {} h {}", c.shape.referenced_vertices, c.shape.boundary_vertices,
                     c.shape.holes));
    CHECK(c.outcome.mesh.triangle_count() == euler_triangles(c.shape));
}

// And the declared n is the true one: every vertex of a generated domain is
// referenced, because the generator places nothing outside the domain or inside
// a hole. This is what stops the Euler check from passing on a wrong mesh whose
// vertex bookkeeping is wrong in a compensating way.
TEST_CASE("every generated vertex is referenced by some triangle", "[cdt][property][euler]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);

    std::set<std::uint32_t> referenced;
    for (const TriangleIndices& t : c.outcome.mesh.triangles()) {
        referenced.insert(t.begin(), t.end());
    }
    CHECK(referenced.size() == c.shape.referenced_vertices);
    CHECK(referenced.size() == c.outcome.mesh.vertices().size());
}

// ---------------------------------------------------------------------------
// Orientation, domain and identity
// ---------------------------------------------------------------------------

// Every output triangle is counterclockwise under an EXACT predicate. This is
// the strongest true statement about triangle shape this project can make:
// "no near-zero-area triangle" is not one -- there is no epsilon here, and a
// legitimate CDT over near-collinear terrain constraints produces slivers.
//
// It is also what kills the forgotten cwTriangles argument. forEachTriangle's
// parameter DEFAULTS TO TRUE and the wrapper must pass false; under the default
// v1 and v2 arrive swapped, which both reverses winding and transposes mask
// bits 0 and 2. One line, two invariants -- this property catches the first and
// the constraint-preservation property the second.
TEST_CASE("every triangle is counterclockwise", "[cdt][property][orientation]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);
    const IndexedMesh2& mesh = c.outcome.mesh;

    for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
        const TriangleIndices& tri = mesh.triangles()[t];
        INFO(std::format("triangle {}", t));
        REQUIRE(DefaultKernel::orient2d(mesh.vertices()[tri[0]], mesh.vertices()[tri[1]],
                                        mesh.vertices()[tri[2]]) == Orientation::CounterClockwise);
    }
}

// Only in-domain triangles come back: nothing inside a Hole, nothing outside
// every Outer. The centroid is the witness -- it is strictly interior to its
// triangle, and the generator's features are separated by whole units, so no
// centroid lands near a ring.
TEST_CASE("no triangle lies in a hole or outside the domain", "[cdt][property][domain]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);
    const IndexedMesh2& mesh = c.outcome.mesh;

    for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
        const Point2 p = centroid(mesh, mesh.triangles()[t]);
        INFO(std::format("triangle {} centroid ({}, {})", t, p.x, p.y));
        bool in_some_outer = false;
        for (std::size_t ch = 0; ch < c.pslg.chains().size(); ++ch) {
            const Chain& chain = c.pslg.chains()[ch];
            if (!is_closed(chain.role)) continue;
            const IndexedRing ring = c.pslg.ring(ch);
            const PointInRing where = point_in_ring<DefaultKernel>(ring, p);
            if (chain.role == ChainRole::Hole) {
                REQUIRE(where != PointInRing::Inside);
            } else if (where == PointInRing::Inside) {
                in_some_outer = true;
            }
        }
        REQUIRE(in_some_outer);
    }
}

// Obligation 1 of the CdtBackend concept: the returned vertex array begins with
// pslg.vertices(), element-wise, in order. DetriaBackend appends nothing, so
// for it the arrays are equal -- the seam permits a future backend that inserts
// Steiner points to append them, which is why the concept says "begins with".
// As of 5c the input side of this identity is the NODED graph, and that is not
// a rename: the two index spaces genuinely differ, because node ids are the
// sorted node set's order. Every downstream consumer -- the mask join, cli.py's
// scene -- reads the mesh against the graph that was handed to the backend.
TEST_CASE("the mesh vertex array is the noded input's", "[cdt][property][identity]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);

    REQUIRE(c.outcome.mesh.vertices().size() == c.pslg.vertices().size());
    for (std::size_t i = 0; i < c.pslg.vertices().size(); ++i) {
        INFO(std::format("vertex {}", i));
        REQUIRE(c.outcome.mesh.vertices()[i] == c.pslg.vertices()[i]);
    }
}

// Guarantees 1-3 of IndexedMesh2, asserted on real output rather than on a
// hand-built mesh: the constructor checks them as debug asserts, which are
// compiled out in Release, and this is the one suite that sees meshes the
// backend actually built.
TEST_CASE("the mesh is internally well-formed", "[cdt][property][identity]") {
    const int seed = GENERATE(range(0, seed_count));
    const Case c = generated(seed);
    const IndexedMesh2& mesh = c.outcome.mesh;

    REQUIRE(mesh.triangles().size() == mesh.constrained_edges().size());
    for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
        REQUIRE(mesh.constrained_edges()[t] < 8u);
        for (const std::uint32_t i : mesh.triangles()[t]) {
            REQUIRE(i < mesh.vertices().size());
        }
    }
}

// No status/mesh disagreement, over generated input on both arms: the Ok arm
// above, and a corruption of each case that must fail and must return nothing.
// A generated hole is moved outside the outer ring, which is the one corruption
// that needs no crossing to be invalid.
TEST_CASE("a corrupted domain fails and returns no mesh", "[cdt][property][identity]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng = seeded(seed);
    GeneratedDomain domain = generated_domain(rng);
    domain.builder.add_chain(terrain::test::points(terrain::test::cw_rect(1000.0, 1000.0,
                                                                         1010.0, 1010.0)),
                             ChainRole::Hole);
    // The corruption is a hole outside every outline, which the NODER accepts --
    // NodedPslgBuilder is not the Pslg validator and computes no nesting
    // (noded_pslg_builder.hpp:45-49). That is what keeps InvalidTopology
    // reachable after the signature change, and it is the only reason this
    // property did not go with the deleted backend cases.
    const NodedPslg pslg = node_fixture<DefaultKernel>(std::move(domain.builder));

    const CdtOutcome out = terrain::cdt::triangulate<DetriaBackend>(pslg);

    CHECK(out.status == CdtStatus::InvalidTopology);
    CHECK_FALSE(out.ok());
    CHECK(out.mesh.empty());
    CHECK_FALSE(out.message.empty());
}
