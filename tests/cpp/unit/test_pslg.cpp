// Unit tests for terrain/core/pslg.hpp: the accessors, the span algebra and the
// value semantics of the immutable constraint set.
//
// NOT invariant-critical. No mutation round is spent here and there is a single
// kernel instantiation: this file's failure mode is a typo, not a topology
// decision. What lives here is everything a downstream module reads --
// indices_of, ring, edge, edge_count, closed_index_buffer_size -- plus the
// value semantics that make `const Pslg&` shareable across refinement threads.
//
// The only way to hold a Pslg is to have passed validation, so every fixture
// here goes through the builder. That is not a limitation of the test: it is
// the property the private constructor exists to enforce, and a test that could
// hand-build a corrupt Pslg would be testing a type the program cannot produce.

#include <catch2/catch_test_macros.hpp>

#include <pslg_cases.hpp>
#include <ring_cases.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using terrain::Chain;
using terrain::ChainRole;
using terrain::IndexedRing;
using terrain::Point2;
using terrain::Pslg;
using terrain::PslgBuilder;
using terrain::PslgBuildResult;
using terrain::Segment2;
using terrain::closed_index_buffer_size;
using terrain::is_closed;
using terrain::pred::DefaultKernel;
using terrain::test::ccw_square;
using terrain::test::cw_hole;
using terrain::test::indices;
using terrain::test::open_breakline;
using terrain::test::points;
using terrain::test::render;
using terrain::test::translated;

namespace {

// One outer square, one clockwise hole, one open breakline: the smallest set
// that exercises both edge_count arms and both closed roles.
[[nodiscard]] Pslg three_chain_pslg() {
    PslgBuilder b;
    b.add_chain(points(ccw_square()), ChainRole::Outer);
    b.add_chain(points(cw_hole()), ChainRole::Hole, true);
    b.add_chain(points(open_breakline()), ChainRole::Breakline);
    PslgBuildResult r = std::move(b).build<DefaultKernel>();
    REQUIRE(r.ok());
    return std::move(*r.pslg);
}

}  // namespace

// ---------------------------------------------------------------------------
// Value semantics
// ---------------------------------------------------------------------------

// The private constructor is the whole enforcement mechanism: the only way to
// hold a Pslg is to have passed validation. A hand-built std::vector<Chain>
// with begin/count that do not match the flat buffer has no public path to a
// Pslg, so the "Pslg-shaped thing the validator never saw" is unrepresentable
// rather than merely discouraged.
static_assert(!std::is_default_constructible_v<Pslg>);
static_assert(!std::is_constructible_v<Pslg, std::vector<Point2>, std::vector<std::uint32_t>,
                                       std::vector<Chain>>);

// Freely copyable and movable, which is what makes `const Pslg&` shareable
// across threads with no synchronisation.
static_assert(std::is_copy_constructible_v<Pslg>);
static_assert(std::is_copy_assignable_v<Pslg>);
static_assert(std::is_nothrow_move_constructible_v<Pslg>);
static_assert(std::is_move_assignable_v<Pslg>);

// Every accessor is const and returns a view of const. There is no mutator, no
// non-const accessor, no lazy cache and no mutable member.
static_assert(std::is_same_v<decltype(std::declval<const Pslg&>().vertices()),
                             std::span<const Point2>>);
static_assert(std::is_same_v<decltype(std::declval<const Pslg&>().chains()),
                             std::span<const Chain>>);
static_assert(std::is_same_v<decltype(std::declval<const Pslg&>().chain_indices()),
                             std::span<const std::uint32_t>>);
static_assert(std::is_same_v<decltype(std::declval<const Pslg&>().indices_of(0)),
                             std::span<const std::uint32_t>>);
static_assert(std::is_same_v<decltype(std::declval<const Pslg&>().ring(0)), IndexedRing>);

// ring() is deliberately NOT noexcept even though it cannot throw. IndexedRing's
// constructor rechecks size, closure and the two boundary indices, and every one
// of those was established by the validator -- but the compiler cannot know
// that, and a noexcept that rests on an argument the compiler cannot see is a
// std::terminate waiting for the argument to stop being true.
static_assert(!noexcept(std::declval<const Pslg&>().ring(0)));
static_assert(noexcept(std::declval<const Pslg&>().vertices()));
static_assert(noexcept(std::declval<const Pslg&>().edge_count(0)));

// ---------------------------------------------------------------------------
// The buffers
// ---------------------------------------------------------------------------

TEST_CASE("the vertex buffer is the builder's input, element for element", "[pslg][accessors]") {
    // No dedup, no reordering, no reversal, no insertion. A caller's index k
    // means the same point going out as it did going in, which is what lets
    // Python hold a parallel attribute array.
    const std::vector<Point2> outer = ccw_square();
    const std::vector<Point2> hole = cw_hole();

    PslgBuilder b;
    b.add_chain(points(outer), ChainRole::Outer);
    b.add_chain(points(hole), ChainRole::Hole);
    const PslgBuildResult r = std::move(b).build<DefaultKernel>();
    INFO(render(r));
    REQUIRE(r.ok());
    const Pslg& p = *r.pslg;

    REQUIRE(p.vertices().size() == outer.size() + hole.size());
    for (std::size_t i = 0; i < outer.size(); ++i) CHECK(p.vertices()[i] == outer[i]);
    for (std::size_t i = 0; i < hole.size(); ++i)
        CHECK(p.vertices()[outer.size() + i] == hole[i]);
}

TEST_CASE("the point-taking add_chain does not dedup coincident points", "[pslg][accessors]") {
    // Two chains over the same coordinates give eight vertices, not four. Dedup
    // is a mutation of caller-declared topology and there is no legal key for
    // it until the snap grid exists.
    PslgBuilder b;
    b.add_chain(points(ccw_square()), ChainRole::Outer);
    b.add_chain(points(ccw_square()), ChainRole::Outer);
    const PslgBuildResult r = std::move(b).build<DefaultKernel>();
    REQUIRE(r.ok());
    CHECK(r.pslg->vertices().size() == 8);
    CHECK(r.pslg->vertices()[0] == r.pslg->vertices()[4]);
}

TEST_CASE("append_vertices returns the index of the first new vertex", "[pslg][builder]") {
    const std::vector<Point2> a = ccw_square();
    const std::vector<Point2> b_pts = cw_hole();

    PslgBuilder b;
    CHECK(b.append_vertices(points(a)) == 0u);
    CHECK(b.append_vertices(points(b_pts)) == static_cast<std::uint32_t>(a.size()));
    CHECK(b.append_vertices(points(a)) ==
          static_cast<std::uint32_t>(a.size() + b_pts.size()));
}

// Guarantee 7, stated as the arithmetic the builder exists to do once: the
// caller never writes begin or count, so an off-by-one in them is
// unrepresentable rather than validated.
TEST_CASE("indices_of sub-spans partition chain_indices contiguously", "[pslg][accessors]") {
    const Pslg p = three_chain_pslg();

    std::size_t expected_begin = 0;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        const Chain& ch = p.chains()[c];
        CHECK(ch.begin == expected_begin);
        CHECK(p.indices_of(c).size() == ch.count);
        CHECK(p.indices_of(c).data() == p.chain_indices().data() + ch.begin);
        expected_begin += ch.count;
    }
    CHECK(expected_begin == p.chain_indices().size());
}

TEST_CASE("roles and the is_river bit survive the build", "[pslg][accessors]") {
    const Pslg p = three_chain_pslg();
    REQUIRE(p.chains().size() == 3);
    CHECK(p.chains()[0].role == ChainRole::Outer);
    CHECK(p.chains()[1].role == ChainRole::Hole);
    CHECK(p.chains()[2].role == ChainRole::Breakline);
    CHECK_FALSE(p.chains()[0].is_river);
    CHECK(p.chains()[1].is_river);
    CHECK_FALSE(p.chains()[2].is_river);
}

// ---------------------------------------------------------------------------
// ring()
// ---------------------------------------------------------------------------

TEST_CASE("ring(c) is a zero-copy view of a closed chain", "[pslg][accessors][ring]") {
    const Pslg p = three_chain_pslg();
    const std::vector<Point2> hole = cw_hole();

    const IndexedRing r = p.ring(1);
    REQUIRE(r.size() == hole.size());
    for (std::size_t i = 0; i < hole.size(); ++i) {
        CHECK(r.vertex(i) == hole[i]);
        // The view really is a view: the reference names the Pslg's own buffer.
        CHECK(&r.vertex(i) == &p.vertices()[p.indices_of(1)[i]]);
    }
}

TEST_CASE("ring(c) never throws for a closed chain", "[pslg][accessors][ring]") {
    const Pslg p = three_chain_pslg();
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        if (!is_closed(p.chains()[c].role)) continue;
        CHECK_NOTHROW(p.ring(c));
    }
}

// ---------------------------------------------------------------------------
// edge() and edge_count()
// ---------------------------------------------------------------------------
//
// This is the accessor the noder's broad phase and the CDT wrapper's
// setConstrainedEdge loop both consume, and it exists so neither of them writes
// (k + 1) % count and gets the open case wrong.

TEST_CASE("edge_count is count for a closed chain and count-1 for an open one",
          "[pslg][accessors][edge]") {
    const Pslg p = three_chain_pslg();
    CHECK(p.edge_count(0) == p.chains()[0].count);
    CHECK(p.edge_count(1) == p.chains()[1].count);
    CHECK(p.edge_count(2) == static_cast<std::size_t>(p.chains()[2].count) - 1);
}

TEST_CASE("edge(c, k) walks the chain and wraps only when it is closed",
          "[pslg][accessors][edge]") {
    const Pslg p = three_chain_pslg();

    SECTION("a closed chain: the last edge closes the ring") {
        const std::span<const std::uint32_t> idx = p.indices_of(0);
        for (std::size_t k = 0; k + 1 < idx.size(); ++k) {
            CHECK(p.edge(0, k) == Segment2{p.vertices()[idx[k]], p.vertices()[idx[k + 1]]});
        }
        const Segment2 closing = p.edge(0, p.edge_count(0) - 1);
        CHECK(closing == Segment2{p.vertices()[idx.back()], p.vertices()[idx.front()]});
    }

    SECTION("an open chain: the last edge is the last segment, and does not close") {
        const std::span<const std::uint32_t> idx = p.indices_of(2);
        const Segment2 last = p.edge(2, p.edge_count(2) - 1);
        CHECK(last == Segment2{p.vertices()[idx[idx.size() - 2]], p.vertices()[idx.back()]});
        CHECK_FALSE(last.b == p.vertices()[idx.front()]);
    }
}

TEST_CASE("edge agrees with the ring's own edge for closed chains",
          "[pslg][accessors][edge][ring]") {
    // Two spellings of the same traversal. They are allowed to be separate
    // implementations; they are not allowed to disagree.
    const Pslg p = three_chain_pslg();
    const IndexedRing r = p.ring(1);
    for (std::size_t k = 0; k < p.edge_count(1); ++k) {
        CHECK(p.edge(1, k) == terrain::edge(r, k));
    }
}

// ---------------------------------------------------------------------------
// closed_index_buffer_size
// ---------------------------------------------------------------------------
//
// The only concession core/ makes to a detria-shaped need. detria wants a
// CLOSED contiguous span and indices_of(c) is not closed, so the wrapper builds
// ONE scratch index buffer per triangulation and hands detria sub-spans of it.
// This function returns the exact element count so the wrapper's shape is
// "reserve exactly this, then fill", under which a per-ring std::vector inside
// the loop is visibly wrong rather than merely wrong.
//
// Note what does not exist: there is deliberately no closed_indices_of(c). The
// only way to get a closed span is to allocate, and the only correct place to
// allocate is once, up front.

TEST_CASE("closed_index_buffer_size counts closed chains only", "[pslg][accessors][cdt_seam]") {
    const Pslg p = three_chain_pslg();

    std::size_t expected = 0;
    for (const Chain& c : p.chains()) {
        if (is_closed(c.role)) expected += static_cast<std::size_t>(c.count) + 1;
    }
    CHECK(closed_index_buffer_size(p) == expected);
    CHECK(expected == (4 + 1) + (4 + 1));  // the breakline contributes nothing
}

TEST_CASE("closed_index_buffer_size is zero when no chain is closed",
          "[pslg][accessors][cdt_seam]") {
    // Not reachable through a valid build -- a Pslg has an Outer chain by
    // guarantee 3 -- so the degenerate arm is asserted by construction instead:
    // an outer-only set contributes exactly its own count plus one.
    PslgBuilder b;
    b.add_chain(points(ccw_square()), ChainRole::Outer);
    b.add_chain(points(open_breakline()), ChainRole::Breakline);
    b.add_chain(points(translated(points(open_breakline()), Point2{1.0, 1.0})),
                ChainRole::Breakline);
    const PslgBuildResult r = std::move(b).build<DefaultKernel>();
    REQUIRE(r.ok());
    CHECK(closed_index_buffer_size(*r.pslg) == 5);
}

// ---------------------------------------------------------------------------
// Copy and move
// ---------------------------------------------------------------------------

// Pslg stores no span and no iterator into itself; every view is computed on
// demand from the vectors. That is what makes the defaulted copy correct. A
// cached IndexedRing member would make a copied Pslg point at the ORIGINAL's
// buffers, and that bug survives every test that never copies -- so this one
// copies, and then checks the addresses rather than only the values.
TEST_CASE("a copied Pslg owns its own buffers", "[pslg][value_semantics]") {
    const Pslg original = three_chain_pslg();
    const Pslg copy = original;

    REQUIRE(copy.vertices().size() == original.vertices().size());
    CHECK(copy.vertices().data() != original.vertices().data());
    CHECK(copy.chain_indices().data() != original.chain_indices().data());
    CHECK(copy.chains().data() != original.chains().data());

    for (std::size_t i = 0; i < original.vertices().size(); ++i)
        CHECK(copy.vertices()[i] == original.vertices()[i]);
    for (std::size_t i = 0; i < original.chains().size(); ++i)
        CHECK(copy.chains()[i] == original.chains()[i]);

    // And the views a copy hands out point into the copy.
    const IndexedRing r = copy.ring(0);
    CHECK(&r.vertex(0) == &copy.vertices()[copy.indices_of(0)[0]]);
    CHECK(&r.vertex(0) != &original.vertices()[original.indices_of(0)[0]]);
}

TEST_CASE("a moved-from Pslg hands its buffers over intact", "[pslg][value_semantics]") {
    Pslg source = three_chain_pslg();
    const Point2* const vertex_data = source.vertices().data();
    const std::size_t n = source.vertices().size();

    const Pslg moved = std::move(source);
    CHECK(moved.vertices().data() == vertex_data);  // no reallocation
    CHECK(moved.vertices().size() == n);
}
