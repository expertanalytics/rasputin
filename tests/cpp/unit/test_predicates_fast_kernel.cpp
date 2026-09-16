// Unit tests for terrain::pred::FastKernel -- the unfiltered kernel.
//
// FastKernel evaluates the determinants in plain double arithmetic with no
// filter and no fallback. It is fast, it is correct on well-separated input,
// and it is wrong on input that is merely close to degenerate. Both halves of
// that sentence are tested here: the last test in this file is the existence
// proof that the module needs a filtered kernel at all.

#include <catch2/catch_test_macros.hpp>

#include <point_families.hpp>
#include <exact_reference.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <format>
#include <random>
#include <string>

using terrain::Point2;
using terrain::pred::FastKernel;
using terrain::pred::FilteredKernel;
using terrain::pred::GeometryKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::pred::orientation_of_sign;
using terrain::pred::reversed;
using terrain::test::RefExact;

namespace {

using Filtered = FilteredKernel<RefExact>;

std::string describe(const Point2& a, const Point2& b, const Point2& c) {
    return std::format("a={} b={} c={}", a, b, c);
}

}  // namespace

TEST_CASE("FastKernel satisfies GeometryKernel", "[predicates][fast_kernel]") {
    STATIC_REQUIRE(GeometryKernel<FastKernel>);
}

TEST_CASE("FastKernel orients the unit triangle", "[predicates][fast_kernel]") {
    const Point2 a{0.0, 0.0};
    const Point2 b{1.0, 0.0};
    const Point2 c{0.0, 1.0};

    REQUIRE(FastKernel::orient2d(a, b, c) == Orientation::CounterClockwise);
    REQUIRE(FastKernel::orient2d(a, c, b) == Orientation::Clockwise);
    REQUIRE(FastKernel::orient2d(a, b, Point2{2.0, 0.0}) == Orientation::Collinear);
}

// The kernel and the free `cross` in point.hpp must not be two independent
// definitions of the same determinant that drift apart. On inputs this small
// every subtraction and product is exact, so the two agree bit for bit and the
// comparison is a definition check rather than a numerical one.
TEST_CASE("FastKernel::orient2d has the sign of cross(b - a, c - a)", "[predicates][fast_kernel]") {
    const Point2 cases[][3] = {
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        {{0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}},
        {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}},
        {{-3.0, 2.0}, {4.0, 2.0}, {1.0, -5.0}},
        {{2.5, -0.5}, {2.5, 8.0}, {-1.25, 0.25}},
        {{1.0, 1.0}, {1.0, 1.0}, {4.0, 9.0}},  // duplicate vertex
        {{7.0, 7.0}, {7.0, 7.0}, {7.0, 7.0}},  // fully coincident
    };

    for (const auto& t : cases) {
        INFO(describe(t[0], t[1], t[2]));
        REQUIRE(FastKernel::orient2d(t[0], t[1], t[2]) ==
                orientation_of_sign(terrain::cross(t[1] - t[0], t[2] - t[0])));
    }
}

TEST_CASE("FastKernel::orient2d is antisymmetric on exact input", "[predicates][fast_kernel]") {
    const Point2 cases[][3] = {
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        {{-3.0, 2.0}, {4.0, 2.0}, {1.0, -5.0}},
        {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}},
        {{8.0, 8.0}, {8.0, 8.0}, {0.0, 1.0}},
    };

    for (const auto& t : cases) {
        INFO(describe(t[0], t[1], t[2]));
        REQUIRE(FastKernel::orient2d(t[0], t[1], t[2]) ==
                reversed(FastKernel::orient2d(t[0], t[2], t[1])));
    }
}

TEST_CASE("FastKernel::orient2d is invariant under cyclic rotation on exact input", "[predicates][fast_kernel]") {
    const Point2 cases[][3] = {
        {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}},
        {{-3.0, 2.0}, {4.0, 2.0}, {1.0, -5.0}},
        {{0.0, 0.0}, {1.0, 1.0}, {2.0, 2.0}},
    };

    for (const auto& t : cases) {
        INFO(describe(t[0], t[1], t[2]));
        const Orientation o = FastKernel::orient2d(t[0], t[1], t[2]);
        REQUIRE(FastKernel::orient2d(t[1], t[2], t[0]) == o);
        REQUIRE(FastKernel::orient2d(t[2], t[0], t[1]) == o);
    }
}

TEST_CASE("FastKernel::incircle classifies a point against an exact circle", "[predicates][fast_kernel]") {
    const Point2 a{5.0, 0.0};
    const Point2 b{0.0, 5.0};
    const Point2 c{-5.0, 0.0};
    REQUIRE(FastKernel::orient2d(a, b, c) == Orientation::CounterClockwise);

    REQUIRE(FastKernel::incircle(a, b, c, Point2{0.0, 0.0}) == Incircle::Inside);
    REQUIRE(FastKernel::incircle(a, b, c, Point2{40.0, 40.0}) == Incircle::Outside);
    REQUIRE(FastKernel::incircle(a, b, c, Point2{3.0, -4.0}) == Incircle::Cocircular);
}

// GeometryKernel promises one meaning of `incircle`, so FastKernel owes the
// same totality and argument-order insensitivity FilteredKernel does. A naive
// 3x3 lifted determinant flips sign when a, b, c is clockwise, so this passes
// only if FastKernel normalizes too. If it does not, the two kernels are not
// interchangeable and the concept is not doing its job.
TEST_CASE("FastKernel::incircle is insensitive to the order of the first three points", "[predicates][fast_kernel]") {
    const Point2 a{5.0, 0.0};
    const Point2 b{0.0, 5.0};
    const Point2 c{-5.0, 0.0};
    const Point2 inside{0.0, 0.0};

    REQUIRE(FastKernel::incircle(a, c, b, inside) == Incircle::Inside);
    REQUIRE(FastKernel::incircle(c, b, a, inside) == Incircle::Inside);
    REQUIRE(FastKernel::incircle(b, a, c, inside) == Incircle::Inside);
}

TEST_CASE("FastKernel::incircle treats a collinear triple as a degenerate circle", "[predicates][fast_kernel]") {
    const Point2 a{0.0, 0.0};
    const Point2 b{1.0, 1.0};
    const Point2 c{2.0, 2.0};
    REQUIRE(FastKernel::incircle(a, b, c, Point2{5.0, -9.0}) == Incircle::Cocircular);
}

// ---------------------------------------------------------------------------
// Why the module exists.
// ---------------------------------------------------------------------------

// A hand-checkable single case. With a at the origin the determinant is
//
//     (2^30 + 1)*(2^30 + 1) - 2^30*(2^30 + 2)
//   = (2^60 + 2^31 + 1)     - (2^60 + 2^31)     = 1,
//
// so the triple is counterclockwise. The first product needs 61 significand
// bits; a double keeps 53, and at that magnitude the representable numbers are
// 128 apart, so it rounds to 2^60 + 2^31 -- exactly the second product. The
// naive difference is therefore 0 and the answer comes out Collinear.
TEST_CASE("A counterclockwise triple whose naive determinant vanishes", "[predicates][fast_kernel][adversarial]") {
    const Point2 a{0.0, 0.0};
    const Point2 b{1073741825.0, 1073741824.0};  // (2^30 + 1, 2^30)
    const Point2 c{1073741826.0, 1073741825.0};  // (2^30 + 2, 2^30 + 1)

    REQUIRE(RefExact::orient2d(a, b, c) == Orientation::CounterClockwise);
    REQUIRE(Filtered::orient2d(a, b, c) == Orientation::CounterClockwise);
}

// The existence proof, stated as an existential rather than pinned to one
// input: over a fixed, seeded family of near-degenerate triples at ~2^50
// magnitude, the unfiltered kernel gets at least one answer wrong. Asserting a
// specific wrong answer for a specific input would make the test hostage to
// floating-point contraction -- whether the compiler fuses `x*y - z*w` into an
// fma changes which inputs FastKernel gets wrong, but not that it gets some
// wrong. Measured counts over this family are in the hundreds without
// contraction and in the dozens with it.
TEST_CASE("FastKernel disagrees with FilteredKernel on near-degenerate input", "[predicates][fast_kernel][adversarial]") {
    std::mt19937_64 rng{20240917};
    int disagreements = 0;

    for (int i = 0; i < 2000; ++i) {
        const auto t = terrain::test::near_degenerate_triple(rng);
        REQUIRE(Filtered::orient2d(t.a, t.b, t.c) == t.expected);
        if (FastKernel::orient2d(t.a, t.b, t.c) != t.expected) {
            ++disagreements;
        }
    }

    REQUIRE(disagreements > 0);
}
