// Unit tests for terrain/core/pslg_builder.hpp: the validator, its diagnostics
// and the order it runs its checks in.
//
// INVARIANT-CRITICAL. This suite and prop_pslg_invariants.cpp are the two the
// increment's mutation round is spent on. Every check below is written to fail
// under a specific mutation of the validator, and the nine mutants the design
// names are called out by the test that kills each.
//
// Three things about this file are deliberate and easy to erode:
//
// THE ORDER IS THE SUBJECT, not the individual checks. Finiteness must precede
// winding because orient2d on a NaN coordinate returns Collinear by
// orientation.hpp's documented total behaviour -- so a validator that winds
// first reports DegenerateRing for a NaN ring and sends the reader looking for
// collinear geometry that is not there. Index range and structure must precede
// stored closure and winding because those two stages dereference indices and
// name a first and last vertex. Each of those orderings has its own test, and
// each asserts the ABSENCE of the diagnostic the wrong order would produce.
//
// THE VALIDATOR NEVER EARLY-RETURNS. Inputs at this boundary are wrong in bulk
// -- a whole layer digitised clockwise, a whole file with a shifted index base
// -- and a channel that reports the first failure turns one fix into N round
// trips through a pipeline whose cheapest stage is a DEM decode. So "N
// independently broken chains produce at least N diagnostics" is asserted
// directly, and so is "one diagnostic from every stage in a single build".
//
// NO TEST ASSERTS MESSAGE TEXT. The message is std::formatted at diagnosis time
// and carries the offending values -- the out-of-range index and its position,
// the observed winding -- precisely because those do not fit the two index
// fields. Pinning its wording freezes prose and produces a suite that fails on
// improvements to its own error messages. Tests assert `error`, `chain` and
// `vertex`, and render() is used only to print what actually happened when an
// assertion fails.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <pslg_cases.hpp>
#include <ring_cases.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

using terrain::Chain;
using terrain::ChainRole;
using terrain::Point2;
using terrain::Pslg;
using terrain::PslgBuilder;
using terrain::PslgBuildResult;
using terrain::PslgError;
using terrain::is_closed;
using terrain::kNoChain;
using terrain::kNoVertex;
using terrain::pred::DefaultKernel;
using terrain::pred::FastKernel;
using terrain::test::bowtie_ring;
using terrain::test::ccw_hole;
using terrain::test::ccw_square;
using terrain::test::closed_polyline;
using terrain::test::collinear_spine;
using terrain::test::count_errors;
using terrain::test::cw_hole;
using terrain::test::cw_square;
using terrain::test::find_error;
using terrain::test::has_error;
using terrain::test::indices;
using terrain::test::open_breakline;
using terrain::test::points;
using terrain::test::render;
using terrain::test::translated;
using terrain::test::utm33_enclosing_outer;
using terrain::test::utm33_sliver_hole;

namespace {

const double quiet_nan = std::numeric_limits<double>::quiet_NaN();
constexpr double inf = std::numeric_limits<double>::infinity();

// A builder carrying one valid counterclockwise outer square, so that a test
// about some other chain does not also trip NoOuterChain and have to reason
// about two diagnostics at once.
[[nodiscard]] PslgBuilder with_outer_square() {
    PslgBuilder b;
    b.add_chain(points(ccw_square()), ChainRole::Outer);
    return b;
}

// build() is rvalue-ref-qualified, so the whole suite spells this.
[[nodiscard]] PslgBuildResult build(PslgBuilder b) {
    return std::move(b).build<DefaultKernel>();
}

}  // namespace

// ---------------------------------------------------------------------------
// The shape of the channel
// ---------------------------------------------------------------------------

// build() is consumed. "Build twice and get two Pslgs sharing nothing" is a
// compile error rather than a subtle question about what the second build sees,
// and a successful build costs zero copies of the vertex data.
template <typename B>
concept BuildableFromLvalue = requires(B& b) { b.template build<DefaultKernel>(); };
template <typename B>
concept BuildableFromRvalue = requires(B& b) { std::move(b).template build<DefaultKernel>(); };

static_assert(BuildableFromRvalue<PslgBuilder>);
static_assert(!BuildableFromLvalue<PslgBuilder>,
              "build() must be rvalue-ref-qualified: an lvalue build would silently copy the "
              "whole vertex buffer");

TEST_CASE("ok() and pslg agree in both directions", "[pslg][builder][diagnostics]") {
    SECTION("a valid set yields a Pslg and no diagnostics") {
        const PslgBuildResult r = build(with_outer_square());
        INFO(render(r));
        REQUIRE(r.ok());
        REQUIRE(r.pslg.has_value());
        REQUIRE(r.diagnostics.empty());
    }

    SECTION("a failed build yields no Pslg at all") {
        PslgBuilder b;
        b.add_chain(points(cw_square()), ChainRole::Outer);  // wrong winding
        const PslgBuildResult r = build(std::move(b));
        REQUIRE_FALSE(r.ok());
        REQUIRE_FALSE(r.pslg.has_value());
        REQUIRE_FALSE(r.diagnostics.empty());
    }
}

// add_chain is total: it records and defers. A builder that can fail
// mid-accumulation needs a second failure channel, and two channels is how half
// the failures end up unreported.
TEST_CASE("add_chain accepts input the validator will reject", "[pslg][builder]") {
    const std::vector<Point2> verts = ccw_square();
    PslgBuilder b{verts};
    const std::vector<std::uint32_t> empty{};
    const std::vector<std::uint32_t> single{0};
    const std::vector<std::uint32_t> out_of_range{0, 1, 99};

    REQUIRE_NOTHROW(b.add_chain(indices(empty), ChainRole::Outer));
    REQUIRE_NOTHROW(b.add_chain(indices(single), ChainRole::Breakline));
    REQUIRE_NOTHROW(b.add_chain(indices(out_of_range), ChainRole::Hole));

    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    REQUIRE_FALSE(r.ok());
}

// ---------------------------------------------------------------------------
// One named test per PslgError enumerator
// ---------------------------------------------------------------------------

// Mutant 9: NoOuterChain dropped, or satisfied by a Hole.
TEST_CASE("NoOuterChain: a constraint set needs a bounded domain", "[pslg][builder][domain]") {
    SECTION("an empty builder has no domain") {
        const PslgBuildResult r = build(PslgBuilder{});
        INFO(render(r));
        const auto* d = find_error(r, PslgError::NoOuterChain);
        REQUIRE(d != nullptr);
        CHECK(d->chain == kNoChain);
        CHECK(d->vertex == kNoVertex);
    }

    SECTION("a Hole does not satisfy it") {
        PslgBuilder b;
        b.add_chain(points(cw_hole()), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE_FALSE(r.ok());
        CHECK(has_error(r, PslgError::NoOuterChain));
    }

    SECTION("a Breakline does not satisfy it") {
        PslgBuilder b;
        b.add_chain(points(open_breakline()), ChainRole::Breakline);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE_FALSE(r.ok());
        CHECK(has_error(r, PslgError::NoOuterChain));
    }

    SECTION("one Outer satisfies it, however many other chains there are") {
        PslgBuilder b = with_outer_square();
        b.add_chain(points(cw_hole()), ChainRole::Hole);
        b.add_chain(points(open_breakline()), ChainRole::Breakline);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }
}

// Mutant 6: count >= 3 weakened to >= 2 for closed chains, or >= 2 to >= 1 for
// breaklines. Both the rejection and the acceptance at the boundary are
// asserted, because a mutant that moves the threshold either way has to fail
// one of them.
TEST_CASE("ChainTooShort: the per-role minimum vertex count", "[pslg][builder][structure]") {
    const std::vector<Point2> two{Point2{0.0, 0.0}, Point2{10.0, 0.0}};
    const std::vector<Point2> one{Point2{0.0, 0.0}};
    const std::vector<Point2> none{};

    SECTION("a two-vertex Outer is not a ring") {
        PslgBuilder b;
        b.add_chain(points(two), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::ChainTooShort);
        REQUIRE(d != nullptr);
        CHECK(d->chain == 0u);
        CHECK(d->vertex == kNoVertex);
    }

    SECTION("a two-vertex Hole is not a ring either") {
        PslgBuilder b = with_outer_square();
        b.add_chain(points(two), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE(has_error(r, PslgError::ChainTooShort, 1u));
    }

    SECTION("a one-vertex Breakline is a point constraint, and this pipeline has none") {
        PslgBuilder b = with_outer_square();
        b.add_chain(points(one), ChainRole::Breakline);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE(has_error(r, PslgError::ChainTooShort, 1u));
    }

    SECTION("an empty chain is too short whatever its role") {
        PslgBuilder b = with_outer_square();
        b.add_chain(points(none), ChainRole::Breakline);
        b.add_chain(points(none), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(has_error(r, PslgError::ChainTooShort, 1u));
        CHECK(has_error(r, PslgError::ChainTooShort, 2u));
    }

    SECTION("three vertices is enough for a ring and two for a breakline") {
        const std::vector<Point2> tri{Point2{1.0, 1.0}, Point2{9.0, 1.0}, Point2{1.0, 9.0}};
        PslgBuilder b;
        b.add_chain(points(tri), ChainRole::Outer);
        b.add_chain(points(two), ChainRole::Breakline);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }
}

// Mutant 4: the index range check using <= instead of <. The boundary case --
// an index exactly equal to vertices().size() -- is the only input that
// separates the two, and it is the one a caller writing a half-open loop
// actually produces.
TEST_CASE("IndexOutOfRange: every index is < vertices().size()",
          "[pslg][builder][index_range]") {
    const std::vector<Point2> verts = ccw_square();

    SECTION("an index exactly equal to the vertex count is out of range") {
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 4};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::IndexOutOfRange);
        REQUIRE(d != nullptr);
        CHECK(d->chain == 0u);
    }

    SECTION("an index far out of range is reported once per occurrence") {
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 77, 2, 99};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE(count_errors(r, PslgError::IndexOutOfRange) >= 2);
    }

    SECTION("the last valid index is accepted") {
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 3};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }
}

// Mutant 5: finiteness restricted to referenced vertices. The scan is over the
// WHOLE vertex buffer because the CDT wrapper hands detria's setPoints the
// entire point array -- an unreferenced NaN reaches the backend regardless of
// whether any chain names it.
TEST_CASE("NonFiniteVertex: the scan covers the whole vertex buffer",
          "[pslg][builder][finiteness]") {
    SECTION("a referenced NaN is rejected") {
        std::vector<Point2> verts = ccw_square();
        verts[2] = Point2{quiet_nan, 40.0};
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 3};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::NonFiniteVertex);
        REQUIRE(d != nullptr);
        CHECK(d->chain == kNoChain);
        CHECK(d->vertex == 2u);
    }

    SECTION("an UNREFERENCED NaN is rejected -- setPoints sees it anyway") {
        std::vector<Point2> verts = ccw_square();
        verts.push_back(Point2{5.0, quiet_nan});  // named by no chain
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 3};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::NonFiniteVertex);
        REQUIRE(d != nullptr);
        CHECK(d->chain == kNoChain);
        CHECK(d->vertex == 4u);
        CHECK_FALSE(r.ok());
    }

    SECTION("both infinities are non-finite too") {
        std::vector<Point2> verts = ccw_square();
        verts.push_back(Point2{inf, 0.0});
        verts.push_back(Point2{0.0, -inf});
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 3};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(count_errors(r, PslgError::NonFiniteVertex) == 2);
    }

    SECTION("a finite unreferenced vertex is fine -- it is not junk, it is a future node") {
        std::vector<Point2> verts = ccw_square();
        verts.push_back(Point2{50.0, 50.0});
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 3};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }
}

// Mutant 3, all three directions: the closure check extended to Breakline; the
// closure check skipped for Hole; the closure check comparing indices instead
// of points.
TEST_CASE("StoredClosure: closure is implied, never stored", "[pslg][builder][closure]") {
    SECTION("an Outer that repeats its first index is rejected") {
        const std::vector<Point2> verts = ccw_square();
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 3, 0};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::StoredClosure);
        REQUIRE(d != nullptr);
        CHECK(d->chain == 0u);
        CHECK(d->vertex == 0u);  // the repeated vertex, named by one index twice
    }

    SECTION("a Hole is checked too -- the rule is per closed role, not per Outer") {
        PslgBuilder b = with_outer_square();
        std::vector<Point2> hole = cw_hole();
        hole.push_back(hole.front());
        b.add_chain(points(hole), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE(has_error(r, PslgError::StoredClosure, 1u));
    }

    SECTION("the comparison is on POINTS, not on indices") {
        // Two distinct indices onto coincident vertices close the ring just as
        // surely as one index used twice, and an implementation that compares
        // idx.front() != idx.back() accepts this.
        std::vector<Point2> verts = ccw_square();
        verts.push_back(verts.front());  // a second index onto the same point
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 1, 2, 3, 4};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::StoredClosure);
        REQUIRE(d != nullptr);
        CHECK(d->chain == 0u);
        CHECK((d->vertex == 0u || d->vertex == 4u));
    }

    SECTION("a Breakline is EXEMPT: coincident endpoints make a closed polyline") {
        // A contour, or a ring road that is not a domain boundary. This is
        // legitimate geometry and the only thing standing between it and a
        // rejection is the role check.
        PslgBuilder b = with_outer_square();
        b.add_chain(points(closed_polyline()), ChainRole::Breakline);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
        CHECK(count_errors(r, PslgError::StoredClosure) == 0);
    }

    SECTION("the same index sequence declared Outer IS rejected") {
        // Same points, same order, one field different. Nothing but the role
        // distinguishes the accepted case above from this one.
        PslgBuilder b;
        b.add_chain(points(closed_polyline()), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE(has_error(r, PslgError::StoredClosure, 0u));
    }
}

// Mutant 1: Outer/Hole winding expectations swapped. Four cases, because a
// swap has to be caught by an acceptance as well as by a rejection -- a suite
// that only tested the two rejections would pass under a validator that
// rejected everything.
TEST_CASE("WrongWinding: Outer is counterclockwise, Hole is clockwise",
          "[pslg][builder][winding]") {
    SECTION("a clockwise Outer is rejected and never silently reversed") {
        PslgBuilder b;
        b.add_chain(points(cw_square()), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::WrongWinding);
        REQUIRE(d != nullptr);
        CHECK(d->chain == 0u);
        CHECK(d->vertex == kNoVertex);
        CHECK_FALSE(r.ok());
    }

    SECTION("a counterclockwise Hole is rejected") {
        PslgBuilder b = with_outer_square();
        b.add_chain(points(ccw_hole()), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE(has_error(r, PslgError::WrongWinding, 1u));
    }

    SECTION("a counterclockwise Outer is accepted") {
        const PslgBuildResult r = build(with_outer_square());
        INFO(render(r));
        CHECK(r.ok());
    }

    SECTION("a clockwise Hole is accepted") {
        PslgBuilder b = with_outer_square();
        b.add_chain(points(cw_hole()), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }

    SECTION("a Breakline has no winding contract at all") {
        // Whatever orientation<K> would say about the polyline's vertices, no
        // winding diagnostic may be raised for an open chain.
        PslgBuilder b = with_outer_square();
        b.add_chain(points(cw_hole()), ChainRole::Breakline);      // "clockwise"
        b.add_chain(points(ccw_hole()), ChainRole::Breakline);     // "counterclockwise"
        b.add_chain(points(collinear_spine()), ChainRole::Breakline);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }
}

// Mutant 2: Collinear accepted rather than rejected, for either role.
TEST_CASE("DegenerateRing: a collinear ring is rejected under both closed roles",
          "[pslg][builder][winding][degenerate]") {
    SECTION("as an Outer") {
        PslgBuilder b;
        b.add_chain(points(collinear_spine()), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        const auto* d = find_error(r, PslgError::DegenerateRing);
        REQUIRE(d != nullptr);
        CHECK(d->chain == 0u);
        CHECK(d->vertex == kNoVertex);
    }

    SECTION("as a Hole") {
        PslgBuilder b = with_outer_square();
        b.add_chain(points(collinear_spine()), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        REQUIRE(has_error(r, PslgError::DegenerateRing, 1u));
        CHECK_FALSE(r.ok());
    }

    SECTION("a zero-area ring of three collinear vertices is degenerate, not wrongly wound") {
        PslgBuilder b;
        const std::vector<Point2> spine{Point2{0.0, 0.0}, Point2{1.0, 1.0}, Point2{2.0, 2.0}};
        b.add_chain(points(spine), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(has_error(r, PslgError::DegenerateRing, 0u));
        CHECK(count_errors(r, PslgError::WrongWinding) == 0);
    }
}

// VertexCountOverflow is the one enumerator this suite cannot exercise: making
// it fire needs a vertex buffer of 2^32 Point2, which is 64 GiB, or a flat
// index buffer of 2^32 uint32_t, which is 16 GiB. Neither is a test. What is
// asserted here is that the check does not misfire on ordinary input and that
// the enumerator is reachable through the formatter, so at least a diagnostic
// carrying it renders rather than indexing a table out of bounds.
//
// The gap is structural: the limit lives inside a member template that needs a
// real buffer to reach. If the implementation exposes the comparison as a
// testable predicate, this test should be replaced by one that calls it.
TEST_CASE("VertexCountOverflow: not raised for buffers that fit in uint32",
          "[pslg][builder][overflow]") {
    PslgBuilder b = with_outer_square();
    b.add_chain(points(cw_hole()), ChainRole::Hole);
    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    CHECK(count_errors(r, PslgError::VertexCountOverflow) == 0);

    const terrain::PslgDiagnostic d{PslgError::VertexCountOverflow, kNoChain, kNoVertex, "x"};
    CHECK_FALSE(terrain::describe(d).empty());
}

// ---------------------------------------------------------------------------
// The order of the stages
// ---------------------------------------------------------------------------

// Mutant 7: stage 5 moved before stage 3 -- winding evaluated before
// finiteness. orient2d on a NaN coordinate returns Collinear, so a
// winding-first validator reports this ring as DegenerateRing: a diagnostic
// that sends the reader to look for collinear geometry that is not there, on a
// ring whose vertices are in general position.
//
// The assertion that does the work is the CHECK_FALSE, not the CHECK. A
// validator in either order reports NonFiniteVertex somewhere; only the correct
// order declines to also report DegenerateRing for the chain.
TEST_CASE("finiteness precedes winding: a NaN ring is NonFiniteVertex, not DegenerateRing",
          "[pslg][builder][order][finiteness]") {
    // A TRIANGLE, and that is not incidental. orientation<K> walks prev and
    // next independently when a triple comes back Collinear, so on a ring of
    // four or more vertices the walk can route AROUND a single NaN vertex and
    // return the correct winding -- which lets a validator that reports the
    // NaN but forgets to exclude the chain from stage 5 pass anyway. On a
    // triangle every triple through the extreme vertex involves the NaN, the
    // walk exhausts, and the wrong order is forced into the open.
    //
    // A four-vertex version of this case is below, and it is the weaker of the
    // two on purpose: it is what an implementation that skips the exclusion
    // step gets away with.
    SECTION("a triangular Hole, where the walk cannot route around the NaN") {
        const std::vector<Point2> hole{Point2{20.0, 20.0}, Point2{20.0, 40.0},
                                       Point2{40.0, quiet_nan}};  // otherwise clockwise
        PslgBuilder b = with_outer_square();
        b.add_chain(points(hole), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));

        INFO(render(r));
        REQUIRE_FALSE(r.ok());
        CHECK(has_error(r, PslgError::NonFiniteVertex));
        CHECK(count_errors(r, PslgError::DegenerateRing) == 0);
        CHECK(count_errors(r, PslgError::WrongWinding) == 0);
    }

    SECTION("every vertex of a triangular Hole non-finite") {
        const std::vector<Point2> hole{Point2{quiet_nan, quiet_nan}, Point2{inf, 40.0},
                                       Point2{40.0, -inf}};
        PslgBuilder b = with_outer_square();
        b.add_chain(points(hole), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));

        INFO(render(r));
        CHECK(count_errors(r, PslgError::NonFiniteVertex) == 3);
        CHECK(count_errors(r, PslgError::DegenerateRing) == 0);
    }

    SECTION("a four-vertex Hole") {
        std::vector<Point2> hole = cw_hole();  // correctly wound as a Hole
        hole[1] = Point2{hole[1].x, quiet_nan};

        PslgBuilder b = with_outer_square();
        b.add_chain(points(hole), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));

        INFO(render(r));
        REQUIRE_FALSE(r.ok());
        CHECK(has_error(r, PslgError::NonFiniteVertex));
        CHECK(count_errors(r, PslgError::DegenerateRing) == 0);
        CHECK(count_errors(r, PslgError::WrongWinding) == 0);
    }
}

// The same ordering question for an Outer, where the wrong order produces
// DegenerateRing for a ring whose finite vertices wind counterclockwise.
TEST_CASE("finiteness precedes winding for an Outer too", "[pslg][builder][order][finiteness]") {
    // A triangle, for the reason the Hole case above spells out: on a longer
    // ring the independent prev/next walk can step over the bad vertex.
    std::vector<Point2> outer{Point2{0.0, 0.0}, Point2{100.0, 0.0}, Point2{-inf, 100.0}};

    PslgBuilder b;
    b.add_chain(points(outer), ChainRole::Outer);
    const PslgBuildResult r = build(std::move(b));

    INFO(render(r));
    CHECK(has_error(r, PslgError::NonFiniteVertex));
    CHECK(count_errors(r, PslgError::DegenerateRing) == 0);
}

// Stage 2 before stage 4. The closure check names vertices[idx.front()] and
// vertices[idx.back()], so running it on a chain with an out-of-range index at
// either end is an out-of-bounds read -- caught by the asan job if it is ever
// reordered, and caught here by the absence of the diagnostic that read would
// have produced.
TEST_CASE("index range precedes stored closure", "[pslg][builder][order][index_range]") {
    const std::vector<Point2> verts = ccw_square();
    PslgBuilder b{verts};
    const std::vector<std::uint32_t> idx{0, 1, 2, 9};  // out of range, and last
    b.add_chain(indices(idx), ChainRole::Outer);
    const PslgBuildResult r = build(std::move(b));

    INFO(render(r));
    CHECK(has_error(r, PslgError::IndexOutOfRange, 0u));
    CHECK(count_errors(r, PslgError::StoredClosure) == 0);
    CHECK(count_errors(r, PslgError::WrongWinding) == 0);
    CHECK(count_errors(r, PslgError::DegenerateRing) == 0);
}

// Stage 1 before stage 5. A two-vertex chain declared Outer would reach
// IndexedRing's constructor, which throws std::invalid_argument on fewer than
// three vertices. That the validator does not throw here is the observable form
// of the design's claim that IndexedRing's checks are unreachable inside it.
TEST_CASE("structure precedes winding: a too-short ring never reaches the kernel",
          "[pslg][builder][order][structure]") {
    const std::vector<Point2> two{Point2{0.0, 0.0}, Point2{10.0, 0.0}};
    PslgBuilder b;
    b.add_chain(points(two), ChainRole::Outer);

    PslgBuildResult r{};
    REQUIRE_NOTHROW(r = build(std::move(b)));
    INFO(render(r));
    CHECK(has_error(r, PslgError::ChainTooShort, 0u));
    CHECK(count_errors(r, PslgError::DegenerateRing) == 0);
    CHECK(count_errors(r, PslgError::WrongWinding) == 0);
}

// ---------------------------------------------------------------------------
// Exhaustiveness
// ---------------------------------------------------------------------------

// Mutant 8: the validator early-returning on the first diagnostic. Inputs at
// this boundary are wrong in bulk, so the first broken chain is never the only
// one.
TEST_CASE("three independently broken chains yield at least three diagnostics",
          "[pslg][builder][exhaustive]") {
    SECTION("three chains broken at the same stage") {
        PslgBuilder b;
        b.add_chain(points(cw_square()), ChainRole::Outer);               // wrong winding
        b.add_chain(points(ccw_hole()), ChainRole::Hole);                 // wrong winding
        b.add_chain(points(translated(points(ccw_hole()), Point2{50.0, 0.0})),
                    ChainRole::Hole);                                     // wrong winding
        const PslgBuildResult r = build(std::move(b));

        INFO(render(r));
        REQUIRE(r.diagnostics.size() >= 3);
        CHECK(has_error(r, PslgError::WrongWinding, 0u));
        CHECK(has_error(r, PslgError::WrongWinding, 1u));
        CHECK(has_error(r, PslgError::WrongWinding, 2u));
    }

    SECTION("a later chain is still validated after an earlier one failed") {
        // Chain 0 fails at stage 4, chain 1 at stage 5. An early return anywhere
        // between them loses the second.
        std::vector<Point2> closed_outer = ccw_square();
        closed_outer.push_back(closed_outer.front());

        PslgBuilder b;
        b.add_chain(points(closed_outer), ChainRole::Outer);
        b.add_chain(points(ccw_hole()), ChainRole::Hole);
        const PslgBuildResult r = build(std::move(b));

        INFO(render(r));
        CHECK(has_error(r, PslgError::StoredClosure, 0u));
        CHECK(has_error(r, PslgError::WrongWinding, 1u));
    }
}

// Every stage reports in the same build. Seven of the eight enumerators appear
// here at once; the eighth needs 64 GiB. An early return at ANY stage loses at
// least one of these, so this single case constrains the whole pass rather than
// one boundary within it.
TEST_CASE("one diagnostic from every stage in a single build", "[pslg][builder][exhaustive]") {
    std::vector<Point2> verts = ccw_square();          // 0..3
    verts.push_back(Point2{quiet_nan, quiet_nan});     // 4, referenced by nothing
    verts.push_back(Point2{50.0, 0.0});                // 5, collinear with 0 and 1
    const std::uint32_t nan_index = 4;

    PslgBuilder b{verts};
    const std::vector<std::uint32_t> too_short{0, 1};
    const std::vector<std::uint32_t> out_of_range{0, 1, 2, 88};
    const std::vector<std::uint32_t> stored_closure{0, 1, 2, 3, 0};
    const std::vector<std::uint32_t> clockwise_hole_as_wrong{0, 1, 2, 3};  // CCW, declared Hole
    // Distinct first and last POINTS, so this chain reaches stage 5 rather than
    // being caught by the stored-closure check on its way there.
    const std::vector<std::uint32_t> collinear{0, 5, 1};

    b.add_chain(indices(too_short), ChainRole::Hole);                // 0: ChainTooShort
    b.add_chain(indices(out_of_range), ChainRole::Hole);             // 1: IndexOutOfRange
    b.add_chain(indices(stored_closure), ChainRole::Hole);           // 2: StoredClosure
    b.add_chain(indices(clockwise_hole_as_wrong), ChainRole::Hole);  // 3: WrongWinding
    b.add_chain(indices(collinear), ChainRole::Hole);                // 4: DegenerateRing
    // and no Outer anywhere:                                           NoOuterChain
    const PslgBuildResult r = build(std::move(b));

    INFO(render(r));
    REQUIRE_FALSE(r.ok());
    CHECK(has_error(r, PslgError::ChainTooShort, 0u));
    CHECK(has_error(r, PslgError::IndexOutOfRange, 1u));
    CHECK(has_error(r, PslgError::StoredClosure, 2u));
    CHECK(has_error(r, PslgError::WrongWinding, 3u));
    CHECK(has_error(r, PslgError::DegenerateRing, 4u));
    CHECK(has_error(r, PslgError::NoOuterChain));

    const auto* nan_diag = find_error(r, PslgError::NonFiniteVertex);
    REQUIRE(nan_diag != nullptr);
    CHECK(nan_diag->vertex == nan_index);
}

// ---------------------------------------------------------------------------
// What a valid Pslg does NOT promise
// ---------------------------------------------------------------------------
//
// These are as load-bearing as the rejections. A valid Pslg is a valid NODER
// input, not a valid CDT input: it asserts everything that can be decided
// without constructing a point, and nothing that cannot. Each acceptance below
// is a decision recorded so that it is pinned rather than accidentally reversed
// by someone adding "one more obvious check".

TEST_CASE("no simplicity promise: a self-intersecting ring is accepted",
          "[pslg][builder][negative]") {
    // Intersection CONSTRUCTION rounds, rounding needs the snap grid, and the
    // snap grid is increment 5's vocabulary. The bowtie's winding under the
    // extreme-vertex rule is counterclockwise, so it is a legal Outer.
    PslgBuilder b;
    b.add_chain(points(bowtie_ring()), ChainRole::Outer);
    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    CHECK(r.ok());
}

TEST_CASE("no disjointness promise: overlapping chains are accepted",
          "[pslg][builder][negative]") {
    PslgBuilder b = with_outer_square();
    b.add_chain(points(translated(points(ccw_square()), Point2{50.0, 50.0})), ChainRole::Outer);
    b.add_chain(points(open_breakline()), ChainRole::Breakline);
    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    CHECK(r.ok());
}

// Nesting is neither declared nor computed, and no parent map is stored. On
// un-noded input the computation is not merely expensive, it is wrong: nesting
// by representative vertex is well defined only for non-crossing rings, and
// nothing here promises non-crossing rings. A lake polygon that pokes a metre
// over the DEM border would be rejected as "hole outside any outer ring", and a
// validator that rejects legitimate input is worse than one that misses an
// illegitimate one.
//
// Both of these produce a mesh the caller may not have intended. That is the
// cost of the deferral, taken knowingly, and these two cases are here so the
// deferral is a decision on record rather than an oversight.
TEST_CASE("no nesting promise: a hole inside no outer ring is accepted",
          "[pslg][builder][negative][nesting]") {
    PslgBuilder b = with_outer_square();  // [0, 100]^2
    b.add_chain(points(translated(points(cw_hole()), Point2{1000.0, 1000.0})), ChainRole::Hole);
    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    CHECK(r.ok());
}

TEST_CASE("no nesting promise: a hole inside two outer rings is accepted",
          "[pslg][builder][negative][nesting]") {
    PslgBuilder b = with_outer_square();
    b.add_chain(points(translated(points(ccw_square()), Point2{10.0, 10.0})), ChainRole::Outer);
    b.add_chain(points(cw_hole()), ChainRole::Hole);  // inside both
    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    CHECK(r.ok());
}

TEST_CASE("no distinctness promise: duplicate and coincident vertices are accepted",
          "[pslg][builder][negative]") {
    SECTION("a repeated consecutive vertex gives a zero-length edge, which is harmless") {
        // {A, A, B, C} is a legal ring. Its rotation {A, B, C, A} is an illegal
        // stored closure -- which is safe only because the builder never
        // rotates a chain.
        const std::vector<Point2> ring{Point2{0.0, 0.0}, Point2{0.0, 0.0}, Point2{10.0, 0.0},
                                       Point2{0.0, 10.0}};
        PslgBuilder b;
        b.add_chain(points(ring), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }

    SECTION("a repeated index is a repeated vertex, and is equally accepted") {
        const std::vector<Point2> verts{Point2{0.0, 0.0}, Point2{10.0, 0.0}, Point2{0.0, 10.0}};
        PslgBuilder b{verts};
        const std::vector<std::uint32_t> idx{0, 0, 1, 2};
        b.add_chain(indices(idx), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }

    SECTION("two coincident vertices in the interior of a ring are accepted") {
        const std::vector<Point2> ring{Point2{0.0, 0.0}, Point2{10.0, 0.0}, Point2{5.0, 5.0},
                                       Point2{5.0, 5.0}, Point2{0.0, 10.0}};
        PslgBuilder b;
        b.add_chain(points(ring), ChainRole::Outer);
        const PslgBuildResult r = build(std::move(b));
        INFO(render(r));
        CHECK(r.ok());
    }
}

TEST_CASE("chains may share vertices, which is what the index encoding is for",
          "[pslg][builder][negative]") {
    // A hole touching its outer ring at a node, two breaklines meeting at a
    // junction: append_vertices once, then the index overload.
    std::vector<Point2> verts = ccw_square();
    verts.push_back(Point2{50.0, 50.0});

    PslgBuilder b;
    const std::uint32_t first = b.append_vertices(points(verts));
    CHECK(first == 0u);

    const std::vector<std::uint32_t> outer{0, 1, 2, 3};
    const std::vector<std::uint32_t> line_a{0, 4};  // shares vertex 0 with the outer
    const std::vector<std::uint32_t> line_b{4, 2};  // shares vertex 4 with line_a
    b.add_chain(indices(outer), ChainRole::Outer);
    b.add_chain(indices(line_a), ChainRole::Breakline);
    b.add_chain(indices(line_b), ChainRole::Breakline);

    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    REQUIRE(r.ok());
    CHECK(r.pslg->vertices().size() == 5);
}

// is_river is data, not structure, and is never validated. It is permitted on
// every role -- a wide river or a lake is legitimately an area feature.
TEST_CASE("is_river is carried through untouched on every role", "[pslg][builder][is_river]") {
    PslgBuilder b;
    b.add_chain(points(ccw_square()), ChainRole::Outer, true);
    b.add_chain(points(cw_hole()), ChainRole::Hole, true);
    b.add_chain(points(open_breakline()), ChainRole::Breakline, true);
    b.add_chain(points(translated(points(open_breakline()), Point2{1.0, 1.0})),
                ChainRole::Breakline, false);
    const PslgBuildResult r = build(std::move(b));

    INFO(render(r));
    REQUIRE(r.ok());
    const Pslg& p = *r.pslg;
    CHECK(p.chains()[0].is_river);
    CHECK(p.chains()[1].is_river);
    CHECK(p.chains()[2].is_river);
    CHECK_FALSE(p.chains()[3].is_river);
}

// ---------------------------------------------------------------------------
// Chain's default role
// ---------------------------------------------------------------------------

// Breakline is the role that promises the least: no winding contract, no
// implied closure, no interior. A Chain that reaches a reader without having
// been through the validator must not be able to claim it is a domain boundary,
// and aggregate initialization of a missing enum member would otherwise give
// ChainRole(0).
//
// Note what is NOT asserted: nothing here pins Breakline's underlying value.
// The NSDMI is what keeps this property if the enumerators are ever reordered,
// and a test asserting ChainRole::Breakline == 2 would fail on that reordering
// while the property it is guarding still held.
TEST_CASE("Chain::role defaults to Breakline, the role that promises the least",
          "[pslg][chain][role]") {
    static_assert(Chain{}.role == ChainRole::Breakline);
    static_assert(Chain{.begin = 7, .count = 3}.role == ChainRole::Breakline);
    static_assert(!is_closed(Chain{}.role));

    static_assert(is_closed(ChainRole::Outer));
    static_assert(is_closed(ChainRole::Hole));
    static_assert(!is_closed(ChainRole::Breakline));

    static_assert(Chain{} == Chain{});
    static_assert(Chain{.begin = 0, .count = 3, .role = ChainRole::Outer} != Chain{});

    // A default-constructed Chain is inert rather than merely wrong: a builder
    // that appended it would produce a zero-count breakline, which the
    // validator rejects as ChainTooShort. There is no path from an
    // aggregate-initialised Chain to a Pslg, and this is the runtime half of
    // that statement.
    PslgBuilder b = with_outer_square();
    const std::vector<std::uint32_t> none{};
    b.add_chain(indices(none), Chain{}.role);
    const PslgBuildResult r = build(std::move(b));
    INFO(render(r));
    CHECK(has_error(r, PslgError::ChainTooShort, 1u));
}

// ---------------------------------------------------------------------------
// The one opt-in kernel instantiation
// ---------------------------------------------------------------------------

// The whole template spend of the increment, and it does one job: the same
// input produces DIFFERENT BUILD OUTCOMES under the two kernels. That proves
// the winding decision actually flows through K -- that it has not been written
// kernel-independently, for instance by reaching for signed_area, whose sign
// cancels to noise on exactly this shape -- and it pins the fact increment 4
// depends on: a Pslg is valid WITH RESPECT TO A KERNEL, and the CDT must be
// built with the same one.
//
// There is no Ring-model cross product here. The validator only ever builds
// IndexedRing; PointRing never appears in this increment, and instantiating
// over it would test a code path production does not have.
//
// See utm33_sliver_hole() in pslg_cases.hpp for the construction and for why it
// is independent of floating-point contraction -- a sliver whose edge-vector
// subtractions stay exact demonstrates nothing under -ffp-contract=off, which
// is how three increment-2 tests passed on arm64 and would have reddened the
// ubuntu CI leg.
TEMPLATE_TEST_CASE("the winding decision flows through K", "[pslg][builder][winding][kernel]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    PslgBuilder b;
    b.add_chain(points(utm33_enclosing_outer()), ChainRole::Outer);
    b.add_chain(points(utm33_sliver_hole()), ChainRole::Hole);
    const PslgBuildResult r = std::move(b).build<TestType>();

    INFO(render(r));
    if constexpr (std::is_same_v<TestType, DefaultKernel>) {
        // The sliver really is clockwise. An exact kernel says so and the build
        // succeeds.
        CHECK(r.ok());
    } else {
        // The naive determinant reports counterclockwise, so the same hole is
        // rejected. The outer ring has integer coordinates and no kernel can
        // disagree about it, so this build fails for exactly one reason.
        CHECK_FALSE(r.ok());
        CHECK(has_error(r, PslgError::WrongWinding, 1u));
        CHECK(count_errors(r, PslgError::WrongWinding) == 1);
    }
}
