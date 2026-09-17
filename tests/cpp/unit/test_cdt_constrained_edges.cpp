// Unit tests for terrain/cdt/constrained_edges.hpp: the constraint edge set and
// the per-triangle mask.
//
// INVARIANT-CRITICAL. A mutation round applies to this file. It is also the
// cheapest suite in the increment -- no triangulation happens here at all, just
// a hand-built Pslg and hand-built triangles -- and that is exactly why the
// mask convention is pinned here rather than downstream: a transposed bit dies
// on a three-line fixture instead of being inferred from a triangle count.
//
// THE CONVENTION, AND THE ONE IT IS NOT. Bit e of a mask is set iff the edge
// from v[e] to v[(e + 1) % 3] is a constraint edge: bit 0 is v0->v1, bit 1 is
// v1->v2, bit 2 is v2->v0. The rejected alternative is CGAL's "edge e is
// opposite vertex e", which a reader with that background will assume; it maps
// edge e to (v[e+1], v[e+2]), which is our bit e+1. The two are A ROTATION OF
// EACH OTHER, so they agree on every triangle with zero or three constrained
// edges and disagree on every other one.
//
// THEREFORE EVERY MASK ASSERTION BELOW THAT MATTERS USES A TRIANGLE WITH
// EXACTLY ONE CONSTRAINED EDGE, and names the bit by number. A test over a
// fully constrained triangle passes under both conventions and under the
// reversed-bit spelling (2,1,0) as well; it proves nothing about which one we
// implement. Do not "simplify" these fixtures into one.

#include <catch2/catch_test_macros.hpp>

#include <cdt_cases.hpp>
#include <pslg_cases.hpp>

#include <terrain/cdt/constrained_edges.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <cstdint>
#include <utility>
#include <vector>

using terrain::ChainRole;
using terrain::Point2;
using terrain::Pslg;
using terrain::PslgBuilder;
using terrain::TriangleIndices;
using terrain::cdt::ConstraintEdgeSet;
using terrain::cdt::constrained_mask;
using terrain::kEdgeMask01;
using terrain::kEdgeMask12;
using terrain::kEdgeMask20;
using terrain::pred::DefaultKernel;
using terrain::test::build_fixture;
using terrain::test::ccw_rect;
using terrain::test::cw_rect;
using terrain::test::indices;
using terrain::test::points;

namespace {

[[nodiscard]] Pslg built(PslgBuilder b) { return build_fixture<DefaultKernel>(std::move(b)); }

// A square (indices 0-3) and one interior breakline (indices 4-5). Five
// constraint edges: the square's four, closing edge included, and the
// breakline's one. Nothing here is symmetric enough to hide a transposition.
[[nodiscard]] Pslg square_and_breakline() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(std::vector<Point2>{Point2{3.0, 5.0}, Point2{7.0, 5.0}}),
                ChainRole::Breakline);
    return built(std::move(b));
}

}  // namespace

// ---------------------------------------------------------------------------
// key(): symmetry, pinned rather than inferred
// ---------------------------------------------------------------------------

TEST_CASE("key is symmetric in its arguments", "[cdt][constraints]") {
    STATIC_REQUIRE(ConstraintEdgeSet::key(2, 7) == ConstraintEdgeSet::key(7, 2));
    STATIC_REQUIRE(ConstraintEdgeSet::key(0, 1) == ConstraintEdgeSet::key(1, 0));
    STATIC_REQUIRE(ConstraintEdgeSet::key(5, 5) == ((std::uint64_t{5} << 32) | 5));
}

// The packing itself, so "symmetric" cannot be satisfied by a function that
// discards information -- key(a, b) = a + b is symmetric too, and collides.
TEST_CASE("key packs min into the high word and max into the low", "[cdt][constraints]") {
    STATIC_REQUIRE(ConstraintEdgeSet::key(2, 7) == ((std::uint64_t{2} << 32) | 7));
    STATIC_REQUIRE(ConstraintEdgeSet::key(7, 2) == ((std::uint64_t{2} << 32) | 7));
    STATIC_REQUIRE(ConstraintEdgeSet::key(0, 0) == 0);

    constexpr std::uint32_t max = 0xFFFF'FFFFu;
    STATIC_REQUIRE(ConstraintEdgeSet::key(0, max) == std::uint64_t{max});
    STATIC_REQUIRE(ConstraintEdgeSet::key(max, max) ==
                   ((std::uint64_t{max} << 32) | std::uint64_t{max}));
    // Distinct pairs get distinct keys at the word boundary, where an off-by-one
    // shift would collide them.
    STATIC_REQUIRE(ConstraintEdgeSet::key(1, 0) != ConstraintEdgeSet::key(0, 2));
}

// ---------------------------------------------------------------------------
// Which edges a Pslg contributes
// ---------------------------------------------------------------------------

// A ring's CLOSING edge is a constraint edge like any other. It is the edge no
// index pair in chain_indices() spells out, so it is the one a loop written as
// `for k in [0, count - 1)` silently drops.
TEST_CASE("a ring contributes its closing edge", "[cdt][constraints]") {
    const Pslg p = square_and_breakline();
    const ConstraintEdgeSet edges{p};

    CHECK(edges.contains(0, 1));
    CHECK(edges.contains(1, 2));
    CHECK(edges.contains(2, 3));
    CHECK(edges.contains(3, 0));  // the closing edge
    CHECK(edges.contains(0, 3));  // and it is undirected
}

// A breakline is OPEN and contributes edge_count(c) edges, not count. The pair
// (first, last) is not one of them -- that is the whole difference between the
// two roles, and the mutant that closes an open chain dies here.
TEST_CASE("a breakline contributes no closing edge", "[cdt][constraints]") {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(std::vector<Point2>{Point2{2.0, 2.0}, Point2{5.0, 8.0}, Point2{8.0, 2.0}}),
                ChainRole::Breakline);
    const Pslg p = built(std::move(b));
    const ConstraintEdgeSet edges{p};

    CHECK(edges.contains(4, 5));
    CHECK(edges.contains(5, 6));
    CHECK_FALSE(edges.contains(6, 4));
    CHECK(edges.size() == 6);  // four ring edges, two breakline edges
}

TEST_CASE("a diagonal of the square is not a constraint edge", "[cdt][constraints]") {
    const Pslg p = square_and_breakline();
    const ConstraintEdgeSet edges{p};

    CHECK_FALSE(edges.contains(0, 2));
    CHECK_FALSE(edges.contains(1, 3));
    CHECK_FALSE(edges.contains(0, 4));
    CHECK(edges.size() == 5);
}

// A hole contributes its edges exactly as an outline does. The set knows
// nothing about roles -- that distinction belongs to the backend, which decides
// addOutline versus addHole, and reappears in the mesh only as which triangles
// exist.
TEST_CASE("a hole contributes its edges like any other ring", "[cdt][constraints]") {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(cw_rect(3.0, 3.0, 6.0, 6.0)), ChainRole::Hole);
    const ConstraintEdgeSet edges{built(std::move(b))};

    CHECK(edges.size() == 8);
    CHECK(edges.contains(4, 5));
    CHECK(edges.contains(7, 4));  // the hole's closing edge
}

// Sorting dedups for free, and a breakline may legitimately repeat an edge --
// a contour digitised twice, a river segment shared with a road. The count is
// what a mutant that forgets to dedup gets wrong.
TEST_CASE("a repeated constraint edge is stored once", "[cdt][constraints]") {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    const std::uint32_t first =
        b.append_vertices(points(std::vector<Point2>{Point2{3.0, 5.0}, Point2{7.0, 5.0}}));
    const std::vector<std::uint32_t> line = {first, first + 1u};
    b.add_chain(indices(line), ChainRole::Breakline);
    b.add_chain(indices(line), ChainRole::Breakline);
    const ConstraintEdgeSet edges{built(std::move(b))};

    CHECK(edges.contains(first, first + 1u));
    CHECK(edges.size() == 5);
}

// A breakline edge that duplicates a ring edge collapses into it too: the same
// index pair, arrived at through two roles.
TEST_CASE("a breakline duplicating a ring edge is stored once", "[cdt][constraints]") {
    PslgBuilder b;
    b.append_vertices(points(ccw_rect(0.0, 0.0, 10.0, 10.0)));
    const std::vector<std::uint32_t> ring = {0u, 1u, 2u, 3u};
    const std::vector<std::uint32_t> line = {1u, 0u};  // the ring's first edge, reversed
    b.add_chain(indices(ring), ChainRole::Outer);
    b.add_chain(indices(line), ChainRole::Breakline);
    const ConstraintEdgeSet edges{built(std::move(b))};

    CHECK(edges.size() == 4);
}

// ---------------------------------------------------------------------------
// The mask, by bit index
// ---------------------------------------------------------------------------
//
// Vertices 4 and 5 are the breakline's; vertex 2 is a square corner. Exactly
// one edge of each triangle below is constrained, so each assertion picks out
// one bit and no rotation of the convention can satisfy it by accident.

TEST_CASE("one constrained edge in slot 0 sets bit 0", "[cdt][constraints][mask]") {
    const ConstraintEdgeSet edges{square_and_breakline()};

    CHECK(constrained_mask(edges, TriangleIndices{4, 5, 2}) == kEdgeMask01);
}

TEST_CASE("one constrained edge in slot 1 sets bit 1", "[cdt][constraints][mask]") {
    const ConstraintEdgeSet edges{square_and_breakline()};

    CHECK(constrained_mask(edges, TriangleIndices{2, 4, 5}) == kEdgeMask12);
}

TEST_CASE("one constrained edge in slot 2 sets bit 2", "[cdt][constraints][mask]") {
    const ConstraintEdgeSet edges{square_and_breakline()};

    CHECK(constrained_mask(edges, TriangleIndices{5, 2, 4}) == kEdgeMask20);
}

// The mask is undirected: a triangle listing the constrained edge backwards
// carries the same bit. Winding is a separate invariant, asserted on real
// output in the backend suite.
TEST_CASE("the mask does not care which way the constrained edge runs",
          "[cdt][constraints][mask]") {
    const ConstraintEdgeSet edges{square_and_breakline()};

    CHECK(constrained_mask(edges, TriangleIndices{5, 4, 2}) == kEdgeMask01);
}

TEST_CASE("an unconstrained triangle has mask zero", "[cdt][constraints][mask]") {
    const ConstraintEdgeSet edges{square_and_breakline()};

    CHECK(constrained_mask(edges, TriangleIndices{0, 2, 4}) == 0u);
}

// Two bits, so that the single-edge fixtures above are not the only shape the
// mask is ever asked about -- and the pair is asymmetric, unlike the all-three
// case.
TEST_CASE("two constrained edges set exactly their two bits", "[cdt][constraints][mask]") {
    const ConstraintEdgeSet edges{square_and_breakline()};

    // (0,1) is a ring edge, (1,2) is a ring edge, (2,0) is a diagonal.
    CHECK(constrained_mask(edges, TriangleIndices{0, 1, 2}) == (kEdgeMask01 | kEdgeMask12));
}

TEST_CASE("a fully constrained triangle has mask seven", "[cdt][constraints][mask]") {
    PslgBuilder b;
    b.add_chain(points(std::vector<Point2>{Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{0.0, 4.0}}),
                ChainRole::Outer);
    const ConstraintEdgeSet edges{built(std::move(b))};

    CHECK(constrained_mask(edges, TriangleIndices{0, 1, 2}) == 7u);
    // ...and that this assertion is convention-blind is exactly why it is not
    // the only one in this section.
}
