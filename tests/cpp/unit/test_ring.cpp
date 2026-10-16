// Unit tests for terrain/core/ring.hpp: IndexedRing, the ring view the PSLG
// hands the algorithms, and orientation over it.
//
// The degeneracy policy is the substance of this file. A ring in this project
// is a sequence of *distinct* vertices with closure implied, never stored.
// Rejecting a stored closure is the highest-value check in the increment: if
// both encodings were accepted, every ring would have two spellings and every
// off-by-one would live in the gap between them. What is accepted, and is
// therefore required to be total here, is everything a real breakline or
// CORINE polygon actually contains -- repeated consecutive vertices, zero-area
// collinear spines, and self-intersections.
//
// Self-intersections are accepted *and not detected*: there is no simplicity
// check in a ring view. The cases that build a ring from a plain vertex list go
// through IndexedRingCase::Holder (ring_cases.hpp), which stores the vertices
// reversed behind a decoy so that only the chain indices give the right ring.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <point_in_ring.hpp>
#include <ring_cases.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using terrain::IndexedRing;
using terrain::Point2;
using terrain::orientation;
using terrain::pred::Orientation;
using terrain::test::ExactRingCases;
using terrain::test::IndexedRingCase;
using terrain::test::PointInRing;
using terrain::test::RingCases;
using terrain::test::bowtie_ring;
using terrain::test::notched_ring;
using terrain::test::unit_square;

namespace {

const double quiet_nan = std::numeric_limits<double>::quiet_NaN();

[[nodiscard]] std::span<const Point2> as_span(const std::vector<Point2>& v) {
    return std::span<const Point2>{v};
}

[[nodiscard]] std::span<const std::uint32_t> as_span(const std::vector<std::uint32_t>& v) {
    return std::span<const std::uint32_t>{v};
}

using Holder = IndexedRingCase::Holder;

}  // namespace

// ---------------------------------------------------------------------------
// Construction: what a ring view refuses to be built from
// ---------------------------------------------------------------------------

// Three vertices that are all equal are closed under the same rule, so they
// are rejected as a closure rather than reaching any algorithm.
TEST_CASE("IndexedRing rejects a ring collapsed to a single point", "[ring][ctor][closure]") {
    const std::vector<Point2> collapsed(3, Point2{2.0, 2.0});
    const Holder h{as_span(collapsed)};

    REQUIRE_THROWS_MATCHES(h.ring(), std::invalid_argument,
                           MessageMatches(ContainsSubstring("closure")));
}

TEST_CASE("IndexedRing accepts the degeneracies that real data contains",
          "[ring][ctor][degenerate]") {
    SECTION("repeated consecutive vertices, i.e. zero-length edges") {
        const std::vector<Point2> pts{
            Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0},
        };
        REQUIRE_NOTHROW(Holder{as_span(pts)}.ring());
    }
    SECTION("an all-collinear, zero-area spine") {
        const std::vector<Point2> pts{Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{2.0, 0.0}};
        REQUIRE_NOTHROW(Holder{as_span(pts)}.ring());
    }
    SECTION("a self-intersecting ring: accepted and not detected") {
        const std::vector<Point2> pts = bowtie_ring();
        REQUIRE_NOTHROW(Holder{as_span(pts)}.ring());
    }
    SECTION("a repeated vertex that is not the stored closure") {
        const std::vector<Point2> pts{
            Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{0.0, 0.0}, Point2{1.0, 1.0},
        };
        REQUIRE_NOTHROW(Holder{as_span(pts)}.ring());
    }
}

// Finiteness is a precondition of the vertex buffer, established once at the
// PSLG and Python boundaries -- not re-litigated by every non-owning view. A
// ring holding a NaN is constructible, and the NaN reaches the caller.
TEST_CASE("IndexedRing does not check finiteness", "[ring][ctor][nan]") {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0}, Point2{quiet_nan, 0.0}, Point2{1.0, 1.0},
    };
    const Holder h{as_span(pts)};

    REQUIRE_NOTHROW(h.ring());
    REQUIRE(std::isnan(h.ring().vertex(1).x));
}

// A ring is a view. Binding one to a temporary vector would dangle on the very
// next line, so the overload is deleted rather than left to a sanitizer.
TEST_CASE("ring views cannot be built from a temporary vector", "[ring][ctor][lifetime]") {
    STATIC_REQUIRE_FALSE(std::is_constructible_v<IndexedRing, std::vector<Point2>&&,
                                                 std::span<const std::uint32_t>>);
    STATIC_REQUIRE_FALSE(std::is_constructible_v<IndexedRing, std::span<const Point2>,
                                                 std::vector<std::uint32_t>&&>);
    STATIC_REQUIRE(std::is_constructible_v<IndexedRing, const std::vector<Point2>&,
                                           const std::vector<std::uint32_t>&>);
}

// ---------------------------------------------------------------------------
// IndexedRing: the indirection, and the checks it does and does not run
// ---------------------------------------------------------------------------

TEST_CASE("IndexedRing resolves vertices through the chain", "[ring][indexed]") {
    const std::vector<Point2> vertices{
        Point2{9.0, 9.0},  // not on this chain
        Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0},
    };
    const std::vector<std::uint32_t> chain{3, 2, 1};
    const IndexedRing r{as_span(vertices), as_span(chain)};

    REQUIRE(r.size() == 3);
    REQUIRE(r.vertex(0) == Point2{1.0, 1.0});
    REQUIRE(r.vertex(1) == Point2{1.0, 0.0});
    REQUIRE(r.vertex(2) == Point2{0.0, 0.0});
}

TEST_CASE("IndexedRing's size is the chain length, not the buffer length", "[ring][indexed]") {
    std::vector<Point2> vertices;
    for (int i = 0; i < 20; ++i) vertices.push_back(Point2{static_cast<double>(i), 0.0});
    vertices[3] = Point2{3.0, 5.0};  // so the chain is not a stored closure
    const std::vector<std::uint32_t> chain{0, 1, 2, 3};
    const IndexedRing r{as_span(vertices), as_span(chain)};

    REQUIRE(r.size() == 4);
}

TEST_CASE("IndexedRing rejects a chain shorter than three", "[ring][indexed][degenerate]") {
    const std::vector<Point2> vertices{Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0}};
    const std::vector<std::uint32_t> two{0, 1};

    REQUIRE_THROWS_AS((IndexedRing{as_span(vertices), as_span(two)}), std::invalid_argument);
}

// Closure is a question about the *points*, not about the indices, so two
// distinct indices into coincident vertices close the ring just as surely as
// the same index used twice.
TEST_CASE("IndexedRing rejects a stored closure by either spelling", "[ring][indexed][closure]") {
    const std::vector<Point2> vertices{
        Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0},
        Point2{0.0, 0.0},  // a duplicate of vertex 0, at a different index
    };

    SECTION("the same index first and last") {
        const std::vector<std::uint32_t> chain{0, 1, 2, 0};
        REQUIRE_THROWS_MATCHES((IndexedRing{as_span(vertices), as_span(chain)}),
                               std::invalid_argument, MessageMatches(ContainsSubstring("closure")));
    }
    SECTION("two indices onto coincident points") {
        const std::vector<std::uint32_t> chain{0, 1, 2, 3};
        REQUIRE_THROWS_MATCHES((IndexedRing{as_span(vertices), as_span(chain)}),
                               std::invalid_argument, MessageMatches(ContainsSubstring("closure")));
    }
}

// Range-checking every index is the PSLG's one-time job, not something every
// zero-copy view redoes. Pinned so that adding the check later is a conscious
// decision rather than a drive-by: the first and last indices stay in range
// here, since the closure check has to dereference those two.
TEST_CASE("IndexedRing does not range-check the chain", "[ring][indexed][policy]") {
    const std::vector<Point2> vertices{Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0}};
    const std::vector<std::uint32_t> chain{0, 99, 2};

    REQUIRE_NOTHROW((IndexedRing{as_span(vertices), as_span(chain)}));
}

// ---------------------------------------------------------------------------
// orientation: exact, and it does not sum
// ---------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("orientation classifies a square both ways round", "[ring][orientation]",
                        RingCases) {
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder ccw{as_span(pts)};
    const std::vector<Point2> reverse = terrain::test::flipped(as_span(pts));
    const typename TestType::Model::Holder cw{as_span(reverse)};

    REQUIRE(orientation<typename TestType::Kernel>(ccw.ring()) == Orientation::CounterClockwise);
    REQUIRE(orientation<typename TestType::Kernel>(cw.ring()) == Orientation::Clockwise);
}

TEMPLATE_LIST_TEST_CASE("orientation is right for a non-convex ring", "[ring][orientation]",
                        RingCases) {
    const std::vector<Point2> pts = notched_ring();
    const typename TestType::Model::Holder h{as_span(pts)};

    REQUIRE(orientation<typename TestType::Kernel>(h.ring()) == Orientation::CounterClockwise);
}

// orientation looks at the extreme vertex -- min y, ties by min x, ties by
// lowest index -- because that vertex is convex in any simple ring. When that
// vertex is repeated, its immediate neighbour is a duplicate of itself and the
// triple is collinear for a reason that says nothing about the ring. The walk
// past such neighbours is what this pins, and a walk that advances prev and
// next in lockstep gets it wrong.
TEMPLATE_LIST_TEST_CASE("orientation looks past a repeated extreme vertex",
                        "[ring][orientation][degenerate]", RingCases) {
    const std::vector<Point2> doubled{
        Point2{0.0, 0.0}, Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{4.0, 4.0},
    };
    const typename TestType::Model::Holder h{as_span(doubled)};

    REQUIRE(orientation<typename TestType::Kernel>(h.ring()) == Orientation::CounterClockwise);
}

// Which extreme vertex the implementation picks is its own business -- every
// convex-hull vertex gives the right answer, so min-y and max-y are both
// defensible. What is *not* negotiable is that the walk past duplicated
// neighbours works wherever it lands, so this ring duplicates every vertex:
// no choice of extreme escapes the degenerate first triple.
TEMPLATE_LIST_TEST_CASE("orientation survives a ring in which every vertex is doubled",
                        "[ring][orientation][degenerate]", RingCases) {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0}, Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{4.0, 0.0},
        Point2{4.0, 4.0}, Point2{4.0, 4.0}, Point2{0.0, 4.0}, Point2{0.0, 4.0},
    };
    const typename TestType::Model::Holder ccw{as_span(pts)};
    const std::vector<Point2> back = terrain::test::flipped(as_span(pts));
    const typename TestType::Model::Holder cw{as_span(back)};

    REQUIRE(orientation<typename TestType::Kernel>(ccw.ring()) == Orientation::CounterClockwise);
    REQUIRE(orientation<typename TestType::Kernel>(cw.ring()) == Orientation::Clockwise);
}

// A collinear run along an edge is the other thing that can make the first
// triple examined useless, depending on where the extreme lands.
TEMPLATE_LIST_TEST_CASE("orientation is right for a ring with collinear runs on every edge",
                        "[ring][orientation][degenerate]", RingCases) {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0}, Point2{2.0, 0.0}, Point2{4.0, 0.0},
        Point2{4.0, 2.0}, Point2{4.0, 4.0},
        Point2{2.0, 4.0}, Point2{0.0, 4.0},
        Point2{0.0, 2.0},
    };
    const typename TestType::Model::Holder ccw{as_span(pts)};
    const std::vector<Point2> back = terrain::test::flipped(as_span(pts));
    const typename TestType::Model::Holder cw{as_span(back)};

    REQUIRE(orientation<typename TestType::Kernel>(ccw.ring()) == Orientation::CounterClockwise);
    REQUIRE(orientation<typename TestType::Kernel>(cw.ring()) == Orientation::Clockwise);
}

TEMPLATE_LIST_TEST_CASE("an all-collinear ring has no orientation",
                        "[ring][orientation][degenerate]", RingCases) {
    const std::vector<Point2> spine{
        Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{3.0, 0.0}, Point2{2.0, 0.0},
    };
    const typename TestType::Model::Holder h{as_span(spine)};

    REQUIRE(orientation<typename TestType::Kernel>(h.ring()) == Orientation::Collinear);
}

TEMPLATE_LIST_TEST_CASE("a vertical all-collinear ring has no orientation",
                        "[ring][orientation][degenerate]", RingCases) {
    const std::vector<Point2> spine{
        Point2{5.0, 0.0}, Point2{5.0, 2.0}, Point2{5.0, 1.0},
    };
    const typename TestType::Model::Holder h{as_span(spine)};

    REQUIRE(orientation<typename TestType::Kernel>(h.ring()) == Orientation::Collinear);
}

// A triangle whose vertices are exactly collinear only at a magnitude where a
// summed shoelace would have cancelled. orientation does not sum, and with an
// exact kernel it is right.
TEMPLATE_LIST_TEST_CASE("orientation is exact on a near-degenerate sliver",
                        "[ring][orientation][sliver]", ExactRingCases) {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0},
        Point2{25304364.0, 25254615.0},
        Point2{715637399591338.0, 714230438918773.0},
    };
    const typename TestType::Model::Holder h{as_span(pts)};

    REQUIRE(orientation<typename TestType::Kernel>(h.ring()) == Orientation::Clockwise);
}

// ---------------------------------------------------------------------------
// point_in_ring (tests/cpp/support/point_in_ring.hpp), the CDT domain oracle
// ---------------------------------------------------------------------------

// The oracle's own direct check, so that a broken oracle cannot pass the
// domain property in prop_cdt_invariants.cpp by agreeing with nothing.
TEMPLATE_LIST_TEST_CASE("point_in_ring classifies a square", "[ring][point_in_ring]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(terrain::test::point_in_ring<K>(r, Point2{0.5, 0.5}) == PointInRing::Inside);
    REQUIRE(terrain::test::point_in_ring<K>(r, Point2{0.25, 0.75}) == PointInRing::Inside);

    REQUIRE(terrain::test::point_in_ring<K>(r, Point2{1.5, 0.5}) == PointInRing::Outside);
    REQUIRE(terrain::test::point_in_ring<K>(r, Point2{-0.5, 0.5}) == PointInRing::Outside);
    REQUIRE(terrain::test::point_in_ring<K>(r, Point2{0.5, 1.5}) == PointInRing::Outside);
    REQUIRE(terrain::test::point_in_ring<K>(r, Point2{0.5, -0.5}) == PointInRing::Outside);
}

// ---------------------------------------------------------------------------
// Where two accepted degeneracies meet
// ---------------------------------------------------------------------------

// A ring whose first vertex is repeated is legal; its cyclic rotation by one
// is not, because the repeat lands on the first and last slot and the closure
// check cannot tell a stored closure apart from a duplicated vertex that
// happens to sit there. The two rules are individually right and together they
// mean "legal ring" is not closed under cyclic rotation.
//
// This is pinned rather than worked around: every consumer that rotates a ring
// -- normalising a chain to start at its lowest vertex, say -- has to know it
// may have to collapse an adjacent duplicate first. Nothing in this increment
// rotates a ring, so nothing here is broken by it.
TEST_CASE("a legal ring's cyclic rotation may be rejected as a closure",
          "[ring][ctor][closure][policy]") {
    const std::vector<Point2> legal{
        Point2{0.0, 0.0}, Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{4.0, 4.0},
    };
    REQUIRE_NOTHROW(Holder{as_span(legal)}.ring());

    const std::vector<Point2> turned{
        Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{4.0, 4.0}, Point2{0.0, 0.0},
    };
    REQUIRE_THROWS_MATCHES(Holder{as_span(turned)}.ring(), std::invalid_argument,
                           MessageMatches(ContainsSubstring("closure")));
}
