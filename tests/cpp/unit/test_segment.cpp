// Unit tests for terrain::Segment2 and its one predicate, on_segment.
//
// Two decisions are pinned here.
//
// First, Segment2's equality is *ordered*: {a, b} != {b, a}. Ring edges are
// directed -- the parity rule in point_in_ring and every future constraint
// edge depend on which way an edge runs -- so a segment that compared equal to
// its own reverse would make "the same edge" an ambiguous phrase in the one
// module where it must not be.
//
// Second, on_segment is exact and constructs nothing. It asks the kernel for
// collinearity and then answers betweenness with closed comparisons between
// coordinates the caller supplied. There is no division, no parameter t, and
// no tolerance, so the answer is correct at UTM33 magnitudes where a computed
// intersection point would have rounded. The last test in this file is the
// existence proof that the kernel choice is load-bearing rather than
// decorative.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <point_families.hpp>

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

using terrain::Box2;
using terrain::Point2;
using terrain::Segment2;
using terrain::bbox;
using terrain::is_degenerate;
using terrain::on_segment;
using terrain::reversed;
using terrain::test::utm33_offset;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();
const double quiet_nan = std::numeric_limits<double>::quiet_NaN();

}  // namespace

// ---------------------------------------------------------------------------
// The value type
// ---------------------------------------------------------------------------

TEST_CASE("Segment2 equality is ordered", "[segment][equality]") {
    const Point2 a{0.0, 0.0};
    const Point2 b{1.0, 2.0};

    REQUIRE(Segment2{a, b} == Segment2{a, b});
    REQUIRE_FALSE(Segment2{a, b} == Segment2{b, a});
    REQUIRE(Segment2{a, b} != Segment2{b, a});
}

TEST_CASE("reversed swaps the endpoints and is an involution", "[segment][reversed]") {
    const Segment2 s{Point2{1.0, 2.0}, Point2{3.0, 4.0}};

    REQUIRE(reversed(s) == Segment2{Point2{3.0, 4.0}, Point2{1.0, 2.0}});
    REQUIRE(reversed(reversed(s)) == s);
}

TEST_CASE("a degenerate segment is its own reverse", "[segment][reversed][degenerate]") {
    const Segment2 s{Point2{5.0, 5.0}, Point2{5.0, 5.0}};

    REQUIRE(reversed(s) == s);
    REQUIRE(is_degenerate(s));
}

TEST_CASE("is_degenerate is exact, not approximate", "[segment][degenerate]") {
    const Point2 p = utm33_offset(Point2{0.0, 0.0});
    const Point2 nudged{std::nextafter(p.x, inf), p.y};

    REQUIRE(is_degenerate(Segment2{p, p}));
    REQUIRE_FALSE(is_degenerate(Segment2{p, nudged}));
    REQUIRE_FALSE(is_degenerate(Segment2{Point2{0.0, 0.0}, Point2{0.0, 1e-300}}));
}

// ---------------------------------------------------------------------------
// bbox
// ---------------------------------------------------------------------------

TEST_CASE("bbox spans both endpoints regardless of their order", "[segment][bbox]") {
    const Segment2 s{Point2{3.0, -1.0}, Point2{-2.0, 4.0}};
    const Box2 expected{Point2{-2.0, -1.0}, Point2{3.0, 4.0}};

    REQUIRE(bbox(s) == expected);
    REQUIRE(bbox(reversed(s)) == expected);
}

TEST_CASE("bbox of an axis-aligned segment is degenerate in that axis", "[segment][bbox]") {
    const Box2 horizontal = bbox(Segment2{Point2{0.0, 2.0}, Point2{5.0, 2.0}});

    REQUIRE_FALSE(horizontal.is_empty());
    REQUIRE(horizontal.height() == 0.0);
    REQUIRE(horizontal.width() == 5.0);
}

TEST_CASE("bbox of a degenerate segment is the degenerate box there", "[segment][bbox][degenerate]") {
    const Point2 p{2.0, 3.0};
    const Box2 b = bbox(Segment2{p, p});

    REQUIRE_FALSE(b.is_empty());
    REQUIRE(b.lo() == p);
    REQUIRE(b.hi() == p);
}

TEST_CASE("bbox rejects a non-finite endpoint", "[segment][bbox][nan]") {
    REQUIRE_THROWS_AS(bbox(Segment2{Point2{0.0, 0.0}, Point2{quiet_nan, 1.0}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(bbox(Segment2{Point2{-inf, 0.0}, Point2{1.0, 1.0}}),
                      std::invalid_argument);
}

// ---------------------------------------------------------------------------
// on_segment, over both kernels
// ---------------------------------------------------------------------------

TEMPLATE_TEST_CASE("on_segment includes both endpoints", "[segment][on_segment]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    const Segment2 s{Point2{0.0, 0.0}, Point2{4.0, 2.0}};

    REQUIRE(on_segment<TestType>(s, s.a));
    REQUIRE(on_segment<TestType>(s, s.b));
}

TEMPLATE_TEST_CASE("on_segment accepts interior points and rejects the rest of the line",
                   "[segment][on_segment]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    const Segment2 s{Point2{0.0, 0.0}, Point2{4.0, 2.0}};

    REQUIRE(on_segment<TestType>(s, Point2{2.0, 1.0}));
    REQUIRE(on_segment<TestType>(s, Point2{1.0, 0.5}));

    // Collinear but outside the closed segment, on both sides.
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{-2.0, -1.0}));
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{6.0, 3.0}));
}

TEMPLATE_TEST_CASE("on_segment rejects points off the supporting line", "[segment][on_segment]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    const Segment2 s{Point2{0.0, 0.0}, Point2{4.0, 2.0}};

    REQUIRE_FALSE(on_segment<TestType>(s, Point2{2.0, 1.5}));
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{2.0, 0.5}));
}

TEMPLATE_TEST_CASE("on_segment does not depend on the segment's direction",
                   "[segment][on_segment]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    const Segment2 s{Point2{-3.0, 7.0}, Point2{5.0, -1.0}};
    const Point2 probes[] = {
        Point2{-3.0, 7.0}, Point2{5.0, -1.0}, Point2{1.0, 3.0},
        Point2{9.0, -5.0}, Point2{0.0, 0.0},
    };

    for (const Point2& p : probes) {
        REQUIRE(on_segment<TestType>(s, p) == on_segment<TestType>(reversed(s), p));
    }
}

// A vertical and a horizontal segment are where a slope-based implementation
// divides by zero. This one never divides.
TEMPLATE_TEST_CASE("on_segment handles axis-aligned segments", "[segment][on_segment]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    const Segment2 vertical{Point2{2.0, -1.0}, Point2{2.0, 5.0}};
    const Segment2 horizontal{Point2{-1.0, 2.0}, Point2{5.0, 2.0}};

    REQUIRE(on_segment<TestType>(vertical, Point2{2.0, 0.0}));
    REQUIRE_FALSE(on_segment<TestType>(vertical, Point2{2.0, 6.0}));
    REQUIRE_FALSE(on_segment<TestType>(vertical, Point2{2.5, 0.0}));

    REQUIRE(on_segment<TestType>(horizontal, Point2{0.0, 2.0}));
    REQUIRE_FALSE(on_segment<TestType>(horizontal, Point2{6.0, 2.0}));
    REQUIRE_FALSE(on_segment<TestType>(horizontal, Point2{0.0, 2.5}));
}

// The degenerate case needs no branch in the implementation: a zero-length
// segment makes every point collinear, so the betweenness comparisons collapse
// to p == a on their own. Pinned because an implementation that "fixes" the
// degenerate case with an explicit branch usually gets it subtly different.
TEMPLATE_TEST_CASE("on_segment on a degenerate segment reduces to p == a",
                   "[segment][on_segment][degenerate]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    const Point2 a{3.0, -2.0};
    const Segment2 s{a, a};

    REQUIRE(on_segment<TestType>(s, a));
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{3.0, -2.5}));
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{3.5, -2.0}));
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{0.0, 0.0}));
}

TEMPLATE_TEST_CASE("on_segment is exact one ulp beyond an endpoint at UTM33 magnitudes",
                   "[segment][on_segment][utm33]",
                   terrain::pred::FastKernel, terrain::pred::DefaultKernel) {
    // Horizontal, so collinearity stays exact in double and the ulp step is
    // purely a betweenness question -- both kernels must agree here.
    const Point2 a = utm33_offset(Point2{0.0, 0.0});
    const Point2 b = utm33_offset(Point2{100.0, 0.0});
    const Segment2 s{a, b};

    REQUIRE(on_segment<TestType>(s, b));
    REQUIRE(on_segment<TestType>(s, Point2{std::nextafter(b.x, -inf), b.y}));
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{std::nextafter(b.x, inf), b.y}));
    REQUIRE_FALSE(on_segment<TestType>(s, Point2{b.x, std::nextafter(b.y, inf)}));
}

// The existence proof. Single-instantiation on purpose: the whole point is
// that the two kernels differ.
//
// Construction: the segment runs from a = -m*d to b = n*d for an integer
// direction d = (dx, dy) and integers m, n, so the origin lies exactly on it
// and every coordinate is an exactly representable integer -- nothing is lost
// before the predicate is called. The loss is in orient2d's *subtraction*, not
// in its products: b - a is (m + n)*d, which needs 54 bits, so it rounds, and
// FastKernel's two edge vectors are no longer parallel. Its cross product comes
// out around -4e15 instead of zero and it reports the origin as off the segment.
//
//     d = (33863043, 43150351), m = 117171955, n = 111172354
//     cross of the rounded edge vectors = -3967798950559065, exactly
//
// Siting the error in the subtraction is the load-bearing part, and the reason
// the constants are frozen here rather than searched for. If the edge vectors
// come out exactly parallel -- which is what happens whenever the coordinates
// are small enough for b - a and p - a to be exact -- then a.x*b.y and a.y*b.x
// are the *same real number*, `cross` rounds them identically, and a
// contraction-free build subtracts them to exactly 0.0: FastKernel is
// accidentally right and this REQUIRE_FALSE fails. Such a case only
// demonstrates anything where the compiler contracts `cross` into a single fma,
// and clang on arm64 does while gcc on baseline x86-64 does not. The earlier
// constants here had that defect and failed the ubuntu CI leg.
//
// These are contraction-independent. Each was checked to give a nonzero cross
// under all three forms a compiler may emit for `a.x*b.y - a.y*b.x`:
// fl(fl(a.x*b.y) - fl(a.y*b.x)), fma(a.x, b.y, -fl(a.y*b.x)) and
// fma(-a.y, b.x, fl(a.x*b.y)). Anyone re-tuning a near-degenerate constant in
// this suite has to clear all three; two of them are invisible on any one host.
TEST_CASE("FastKernel misses an exactly-on-segment point that DefaultKernel finds",
          "[segment][on_segment][fast_kernel][wrong]") {
    const Segment2 s{Point2{-3967798950559065.0, -5056010985606205.0},
                     Point2{3764634203913222.0, 4797126096596254.0}};
    const Point2 p{0.0, 0.0};  // the origin: a + m*d, exactly on s

    REQUIRE(on_segment<terrain::pred::DefaultKernel>(s, p));
    REQUIRE_FALSE(on_segment<terrain::pred::FastKernel>(s, p));
}
