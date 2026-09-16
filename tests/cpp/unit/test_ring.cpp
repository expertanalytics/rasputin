// Unit tests for terrain/core/ring.hpp: the Ring concept, its two non-owning
// models, and the algorithms over them.
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
// Self-intersections are accepted *and not detected*. There is no simplicity
// check and no segment-segment intersection anywhere in this increment,
// because intersection construction rounds and rounding needs the snap grid,
// which is increment 5's vocabulary. The even-odd rule gives a non-simple ring
// a total and deterministic classification, and that answer is specified, not
// incidental -- so it is pinned below.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <point_families.hpp>
#include <ring_cases.hpp>

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using terrain::Box2;
using terrain::IndexedRing;
using terrain::Point2;
using terrain::PointInRing;
using terrain::PointRing;
using terrain::Ring;
using terrain::Segment2;
using terrain::all_finite;
using terrain::bounding_box;
using terrain::edge;
using terrain::orientation;
using terrain::point_in_ring;
using terrain::signed_area;
using terrain::pred::DefaultKernel;
using terrain::pred::FastKernel;
using terrain::pred::Orientation;
using terrain::test::ExactRingCases;
using terrain::test::RingCases;
using terrain::test::bowtie_ring;
using terrain::test::notched_ring;
using terrain::test::unit_square;
using terrain::test::utm33_offset;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();
const double quiet_nan = std::numeric_limits<double>::quiet_NaN();

[[nodiscard]] std::span<const Point2> as_span(const std::vector<Point2>& v) {
    return std::span<const Point2>{v};
}

[[nodiscard]] std::span<const std::uint32_t> as_span(const std::vector<std::uint32_t>& v) {
    return std::span<const std::uint32_t>{v};
}

// Probe types for the concept. Declarations only; nothing is ever called.
struct NoVertex {
    [[nodiscard]] std::size_t size() const;
};
struct SizeIsInt {
    [[nodiscard]] int size() const;
    [[nodiscard]] const Point2& vertex(std::size_t) const;
};
struct VertexByValue {
    [[nodiscard]] std::size_t size() const;
    [[nodiscard]] Point2 vertex(std::size_t) const;
};

// A third model, defined here rather than in the header, to prove the
// algorithms really are written against the concept and not against the two
// shipped classes.
class ArrayTriangle {
public:
    explicit constexpr ArrayTriangle(std::array<Point2, 3> v) : v_{v} {}
    [[nodiscard]] std::size_t size() const noexcept { return v_.size(); }
    [[nodiscard]] const Point2& vertex(std::size_t i) const noexcept { return v_[i]; }

private:
    std::array<Point2, 3> v_;
};

}  // namespace

// ---------------------------------------------------------------------------
// The concept
// ---------------------------------------------------------------------------

TEST_CASE("the shipped models satisfy Ring", "[ring][concept]") {
    STATIC_REQUIRE(Ring<PointRing>);
    STATIC_REQUIRE(Ring<IndexedRing>);
    STATIC_REQUIRE(Ring<ArrayTriangle>);
}

// The concept is deliberately narrow: exactly size() and vertex(i), nothing
// about iterators, ownership or storage. A container that merely looks
// ring-shaped is not one.
TEST_CASE("Ring requires exactly size() and vertex(i), with exact types", "[ring][concept]") {
    STATIC_REQUIRE_FALSE(Ring<std::vector<Point2>>);
    STATIC_REQUIRE_FALSE(Ring<std::span<const Point2>>);
    STATIC_REQUIRE_FALSE(Ring<NoVertex>);
    STATIC_REQUIRE_FALSE(Ring<SizeIsInt>);
    // vertex(i) must return a reference: the algorithms iterate, and a
    // by-value model would silently copy on every access.
    STATIC_REQUIRE_FALSE(Ring<VertexByValue>);
}

// ---------------------------------------------------------------------------
// Construction: what a ring view refuses to be built from
// ---------------------------------------------------------------------------

TEST_CASE("PointRing rejects fewer than three vertices", "[ring][ctor][degenerate]") {
    const std::vector<Point2> none;
    const std::vector<Point2> one{Point2{0.0, 0.0}};
    const std::vector<Point2> two{Point2{0.0, 0.0}, Point2{1.0, 0.0}};

    REQUIRE_THROWS_AS(PointRing{as_span(none)}, std::invalid_argument);
    REQUIRE_THROWS_AS(PointRing{as_span(one)}, std::invalid_argument);
    REQUIRE_THROWS_AS(PointRing{as_span(two)}, std::invalid_argument);
}

// The check the whole encoding rests on. GeoJSON stores closed rings; this
// project does not, and the conversion happens once, at the Python boundary.
TEST_CASE("PointRing rejects a stored closure", "[ring][ctor][closure]") {
    const std::vector<Point2> closed{
        Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0}, Point2{0.0, 0.0},
    };

    REQUIRE_THROWS_MATCHES(PointRing{as_span(closed)}, std::invalid_argument,
                           MessageMatches(ContainsSubstring("closure")));
}

// Three vertices that are all equal are closed under the same rule, so they
// are rejected as a closure rather than reaching any algorithm.
TEST_CASE("PointRing rejects a ring collapsed to a single point", "[ring][ctor][closure]") {
    const std::vector<Point2> collapsed(3, Point2{2.0, 2.0});

    REQUIRE_THROWS_AS(PointRing{as_span(collapsed)}, std::invalid_argument);
}

TEST_CASE("PointRing accepts the degeneracies that real data contains", "[ring][ctor][degenerate]") {
    SECTION("repeated consecutive vertices, i.e. zero-length edges") {
        const std::vector<Point2> pts{
            Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0},
        };
        REQUIRE_NOTHROW(PointRing{as_span(pts)});
    }
    SECTION("an all-collinear, zero-area spine") {
        const std::vector<Point2> pts{Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{2.0, 0.0}};
        REQUIRE_NOTHROW(PointRing{as_span(pts)});
    }
    SECTION("a self-intersecting ring: accepted and not detected") {
        const std::vector<Point2> pts = bowtie_ring();
        REQUIRE_NOTHROW(PointRing{as_span(pts)});
    }
    SECTION("a repeated vertex that is not the stored closure") {
        const std::vector<Point2> pts{
            Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{0.0, 0.0}, Point2{1.0, 1.0},
        };
        REQUIRE_NOTHROW(PointRing{as_span(pts)});
    }
}

// Finiteness is a precondition of the vertex buffer, established once at the
// PSLG and Python boundaries -- not re-litigated by every non-owning view. A
// ring holding a NaN is constructible; all_finite is how a caller finds out.
TEST_CASE("PointRing does not check finiteness", "[ring][ctor][nan]") {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0}, Point2{quiet_nan, 0.0}, Point2{1.0, 1.0},
    };

    REQUIRE_NOTHROW(PointRing{as_span(pts)});
    const PointRing r{as_span(pts)};
    REQUIRE_FALSE(all_finite(r));
}

// A ring is a view. Binding one to a temporary vector would dangle on the very
// next line, so the overload is deleted rather than left to a sanitizer.
TEST_CASE("ring views cannot be built from a temporary vector", "[ring][ctor][lifetime]") {
    STATIC_REQUIRE_FALSE(std::is_constructible_v<PointRing, std::vector<Point2>&&>);
    STATIC_REQUIRE(std::is_constructible_v<PointRing, const std::vector<Point2>&>);

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
// edge
// ---------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("edge(r, i) is directed and wraps at the last vertex", "[ring][edge]",
                        RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(edge(r, 0) == Segment2{Point2{0.0, 0.0}, Point2{1.0, 0.0}});
    REQUIRE(edge(r, 3) == Segment2{Point2{0.0, 1.0}, Point2{0.0, 0.0}});

    for (std::size_t i = 0; i < r.size(); ++i) {
        const Segment2 e = edge(r, i);
        REQUIRE(e.a == r.vertex(i));
        REQUIRE(e.b == edge(r, (i + 1) % r.size()).a);
        REQUIRE(terrain::on_segment<K>(e, e.a));
        REQUIRE(terrain::on_segment<K>(e, e.b));
    }
}

TEMPLATE_LIST_TEST_CASE("a repeated vertex produces a degenerate edge", "[ring][edge][degenerate]",
                        RingCases) {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0},
    };
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(terrain::is_degenerate(edge(r, 1)));
    REQUIRE_FALSE(terrain::is_degenerate(edge(r, 0)));
}

// ---------------------------------------------------------------------------
// all_finite and bounding_box
// ---------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("all_finite reports on every vertex", "[ring][all_finite]", RingCases) {
    const std::vector<Point2> clean = unit_square();
    const typename TestType::Model::Holder ok{as_span(clean)};
    REQUIRE(all_finite(ok.ring()));

    const std::vector<Point2> nan_y{Point2{0.0, 0.0}, Point2{1.0, quiet_nan}, Point2{1.0, 1.0}};
    const typename TestType::Model::Holder bad{as_span(nan_y)};
    REQUIRE_FALSE(all_finite(bad.ring()));

    const std::vector<Point2> infinite{Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{inf, 1.0}};
    const typename TestType::Model::Holder unbounded{as_span(infinite)};
    REQUIRE_FALSE(all_finite(unbounded.ring()));
}

TEMPLATE_LIST_TEST_CASE("bounding_box of a ring spans its vertices", "[ring][bbox]", RingCases) {
    const std::vector<Point2> pts = notched_ring();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();
    const Box2 b = bounding_box(r);

    REQUIRE(b == Box2{Point2{0.0, 0.0}, Point2{6.0, 6.0}});
    for (std::size_t i = 0; i < r.size(); ++i) {
        REQUIRE(b.contains(r.vertex(i)));
    }
}

// A ring always has at least three vertices, so its bounding box is never the
// empty one -- but it may be degenerate in an axis.
TEMPLATE_LIST_TEST_CASE("bounding_box of a collinear ring is an axis-degenerate box",
                        "[ring][bbox][degenerate]", RingCases) {
    const std::vector<Point2> pts{Point2{0.0, 5.0}, Point2{2.0, 5.0}, Point2{4.0, 5.0}};
    const typename TestType::Model::Holder h{as_span(pts)};
    const Box2 b = bounding_box(h.ring());

    REQUIRE_FALSE(b.is_empty());
    REQUIRE(b.height() == 0.0);
    REQUIRE(b.width() == 4.0);
}

// ---------------------------------------------------------------------------
// signed_area: approximate, kernel-free, and not a topology oracle
// ---------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("signed_area is positive for a counterclockwise ring", "[ring][area]",
                        RingCases) {
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder ccw{as_span(pts)};
    const std::vector<Point2> reverse = terrain::test::flipped(as_span(pts));
    const typename TestType::Model::Holder cw{as_span(reverse)};

    REQUIRE(signed_area(ccw.ring()) == 1.0);
    REQUIRE(signed_area(cw.ring()) == -1.0);
}

TEMPLATE_LIST_TEST_CASE("signed_area of the notched ring accounts for the notch", "[ring][area]",
                        RingCases) {
    // 6x6 square minus the triangle (6,6)-(3,2)-(0,6), whose area is 12.
    const std::vector<Point2> pts = notched_ring();
    const typename TestType::Model::Holder h{as_span(pts)};

    REQUIRE(signed_area(h.ring()) == 24.0);
}

TEMPLATE_LIST_TEST_CASE("signed_area of a collinear ring is zero", "[ring][area][degenerate]",
                        RingCases) {
    const std::vector<Point2> pts{Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{2.0, 0.0}};
    const typename TestType::Model::Holder h{as_span(pts)};

    REQUIRE(signed_area(h.ring()) == 0.0);
}

// The bowtie encloses two unit-ish lobes of opposite algebraic sign. Its
// signed area is exactly zero, which is not a statement about how much ground
// it covers. First reason the sign of this function must never decide
// anything topological.
TEMPLATE_LIST_TEST_CASE("signed_area of a bowtie cancels to zero", "[ring][area][non_simple]",
                        RingCases) {
    const std::vector<Point2> pts = bowtie_ring();
    const typename TestType::Model::Holder h{as_span(pts)};

    REQUIRE(signed_area(h.ring()) == 0.0);
    REQUIRE(point_in_ring<DefaultKernel>(h.ring(), Point2{1.0, 2.0}) == PointInRing::Inside);
}

// Second, and the reason the header says so in as many words: signed_area is
// plain double. On a sliver whose vertex coordinates are exactly representable
// but whose products are not, the shoelace sum has no correct digits left and
// its sign is simply wrong, while the exact orientation is right. The
// constants come from the near-degenerate family in point_families.hpp and are
// hard-coded so this is a fact rather than a flaky search.
TEST_CASE("sign(signed_area) disagrees with the exact orientation on a sliver",
          "[ring][area][sliver][wrong]") {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0},
        Point2{25304364.0, 25254615.0},
        Point2{715637399591338.0, 714230438918773.0},
    };
    const PointRing r{as_span(pts)};

    REQUIRE(signed_area(r) > 0.0);
    REQUIRE(orientation<DefaultKernel>(r) == Orientation::Clockwise);
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
    REQUIRE(signed_area(ccw.ring()) == 16.0);
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
// point_in_ring
// ---------------------------------------------------------------------------

TEST_CASE("PointInRing's underlying values are the sign convention", "[ring][point_in_ring]") {
    STATIC_REQUIRE(static_cast<int>(PointInRing::Outside) == -1);
    STATIC_REQUIRE(static_cast<int>(PointInRing::Boundary) == 0);
    STATIC_REQUIRE(static_cast<int>(PointInRing::Inside) == 1);
}

TEMPLATE_LIST_TEST_CASE("point_in_ring classifies a square", "[ring][point_in_ring]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(point_in_ring<K>(r, Point2{0.5, 0.5}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{0.25, 0.75}) == PointInRing::Inside);

    REQUIRE(point_in_ring<K>(r, Point2{1.5, 0.5}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{-0.5, 0.5}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{0.5, 1.5}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{0.5, -0.5}) == PointInRing::Outside);
}

TEMPLATE_LIST_TEST_CASE("every vertex and every edge point is Boundary",
                        "[ring][point_in_ring][boundary]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    for (std::size_t i = 0; i < r.size(); ++i) {
        REQUIRE(point_in_ring<K>(r, r.vertex(i)) == PointInRing::Boundary);
    }
    REQUIRE(point_in_ring<K>(r, Point2{0.5, 0.0}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{1.0, 0.5}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{0.5, 1.0}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{0.0, 0.5}) == PointInRing::Boundary);
}

// One ulp off the boundary is not on it. There is no tolerance in this header.
TEMPLATE_LIST_TEST_CASE("Boundary is exact, not approximate", "[ring][point_in_ring][boundary]",
                        RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(point_in_ring<K>(r, Point2{0.5, std::nextafter(0.0, inf)}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{0.5, std::nextafter(0.0, -inf)}) == PointInRing::Outside);
}

// The half-open-in-y rule exists for exactly this geometry: a horizontal ray
// at y = 2 passes through the notch vertex (3, 2), whose two neighbours are
// both above it. Counting that vertex once -- or twice, or not at all --
// changes the parity for every query point on that line.
TEMPLATE_LIST_TEST_CASE("a ray through a local-minimum vertex keeps its parity",
                        "[ring][point_in_ring][half_open]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = notched_ring();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(point_in_ring<K>(r, Point2{-1.0, 2.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{1.0, 2.0}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{5.0, 2.0}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{7.0, 2.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{3.0, 2.0}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{3.0, 1.0}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{3.0, 3.0}) == PointInRing::Outside);
}

// The bottom edge lies *on* the ray. A horizontal edge must contribute nothing
// to the parity, or every point level with it flips.
TEMPLATE_LIST_TEST_CASE("a horizontal edge on the ray contributes no crossing",
                        "[ring][point_in_ring][half_open]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = unit_square();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(point_in_ring<K>(r, Point2{-1.0, 0.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{2.0, 0.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{-1.0, 1.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{2.0, 1.0}) == PointInRing::Outside);
}

// A diamond's left and right vertices sit exactly on the ray through the
// centre, and both are simple crossings rather than local extrema.
TEMPLATE_LIST_TEST_CASE("a ray through two opposite vertices still classifies",
                        "[ring][point_in_ring][half_open]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts{
        Point2{0.0, -1.0}, Point2{1.0, 0.0}, Point2{0.0, 1.0}, Point2{-1.0, 0.0},
    };
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(point_in_ring<K>(r, Point2{-5.0, 0.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{5.0, 0.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{0.0, 0.0}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{-1.0, 0.0}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{1.0, 0.0}) == PointInRing::Boundary);
}

// Zero-length edges contribute nothing by construction -- there is no special
// case for them in the header, and there must not need to be. The classifying
// answer is identical to the same ring with the duplicates removed.
TEMPLATE_LIST_TEST_CASE("repeated vertices do not change the classification",
                        "[ring][point_in_ring][degenerate]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> plain = unit_square();
    const std::vector<Point2> doubled{
        Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 0.0},
        Point2{1.0, 1.0}, Point2{1.0, 1.0}, Point2{0.0, 1.0},
    };
    const typename TestType::Model::Holder a{as_span(plain)};
    const typename TestType::Model::Holder b{as_span(doubled)};

    const Point2 probes[] = {
        Point2{0.5, 0.5}, Point2{1.5, 0.5}, Point2{1.0, 0.5},
        Point2{0.0, 0.0}, Point2{1.0, 1.0}, Point2{0.5, 1.0},
        Point2{-1.0, 0.0}, Point2{0.5, 2.0},
    };
    for (const Point2& p : probes) {
        REQUIRE(point_in_ring<K>(a.ring(), p) == point_in_ring<K>(b.ring(), p));
    }
}

// An all-collinear ring has no interior. Every point on the spine is Boundary
// -- including the stretch covered only by the closing edge -- and everything
// else is Outside.
TEMPLATE_LIST_TEST_CASE("a collinear ring is all boundary and no interior",
                        "[ring][point_in_ring][degenerate]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> spine{Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{4.0, 0.0}};
    const typename TestType::Model::Holder h{as_span(spine)};
    const auto r = h.ring();

    REQUIRE(point_in_ring<K>(r, Point2{0.5, 0.0}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{3.0, 0.0}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{4.0, 0.0}) == PointInRing::Boundary);
    REQUIRE(point_in_ring<K>(r, Point2{4.5, 0.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{-0.5, 0.0}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{2.0, 0.5}) == PointInRing::Outside);
}

// The specified answer for a non-simple ring, pinned so that "even-odd" is a
// contract rather than a description of whatever the loop happens to do. Both
// lobes of the bowtie are Inside; the crossing point is on two edges and is
// therefore Boundary.
TEMPLATE_LIST_TEST_CASE("a self-intersecting ring gets the even-odd classification",
                        "[ring][point_in_ring][non_simple]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> pts = bowtie_ring();
    const typename TestType::Model::Holder h{as_span(pts)};
    const auto r = h.ring();

    REQUIRE(point_in_ring<K>(r, Point2{1.0, 2.0}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{3.0, 2.0}) == PointInRing::Inside);
    REQUIRE(point_in_ring<K>(r, Point2{2.0, 0.5}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{2.0, 3.5}) == PointInRing::Outside);
    REQUIRE(point_in_ring<K>(r, Point2{2.0, 2.0}) == PointInRing::Boundary);
}

// point_in_ring performs no arithmetic of its own: every numeric decision is a
// kernel call or a comparison between two coordinates the caller supplied. The
// observable consequence is that translating the whole problem by an exactly
// representable UTM33 offset -- which leaves every coordinate exact but would
// wreck any internally computed x-intersection -- cannot change an answer.
TEMPLATE_LIST_TEST_CASE("point_in_ring is invariant under an exact UTM33 translation",
                        "[ring][point_in_ring][utm33]", RingCases) {
    using K = typename TestType::Kernel;
    const std::vector<Point2> local = notched_ring();
    const Point2 offset = utm33_offset(Point2{0.0, 0.0});
    const std::vector<Point2> shifted = terrain::test::translated(as_span(local), offset);

    const typename TestType::Model::Holder a{as_span(local)};
    const typename TestType::Model::Holder b{as_span(shifted)};

    const Point2 probes[] = {
        Point2{1.0, 2.0}, Point2{3.0, 3.0}, Point2{3.0, 1.0}, Point2{-1.0, 2.0},
        Point2{3.0, 2.0}, Point2{0.0, 0.0}, Point2{6.0, 6.0}, Point2{5.0, 2.0},
    };
    for (const Point2& p : probes) {
        REQUIRE(point_in_ring<K>(a.ring(), p) ==
                point_in_ring<K>(b.ring(), Point2{p.x + offset.x, p.y + offset.y}));
    }
}

TEST_CASE("the algorithms work through the concept, not the shipped models", "[ring][concept]") {
    const ArrayTriangle t{{Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{0.0, 4.0}}};

    REQUIRE(orientation<DefaultKernel>(t) == Orientation::CounterClockwise);
    REQUIRE(signed_area(t) == 8.0);
    REQUIRE(point_in_ring<DefaultKernel>(t, Point2{1.0, 1.0}) == PointInRing::Inside);
    REQUIRE(point_in_ring<DefaultKernel>(t, Point2{2.0, 2.0}) == PointInRing::Boundary);
    REQUIRE(bounding_box(t) == Box2{Point2{0.0, 0.0}, Point2{4.0, 4.0}});
}

// The existence proof that the kernel parameter on point_in_ring is
// load-bearing. The query point lies exactly on the edge from (0,0) to
// (424845863969746, 748848292265877) -- it is that vertex divided by the
// integer 22350077 -- so the true answer is Boundary. FastKernel's cross
// product is a difference of two independently rounded 75-bit products, misses
// the collinearity, and reports the point as interior. Single instantiation on
// purpose: the point is that the two kernels differ.
TEST_CASE("FastKernel misclassifies a boundary point that DefaultKernel gets right",
          "[ring][point_in_ring][fast_kernel][wrong]") {
    const std::vector<Point2> pts{
        Point2{0.0, 0.0},
        Point2{424845863969746.0, 748848292265877.0},
        Point2{0.0, 748848292265877.0},
    };
    const PointRing r{as_span(pts)};
    const Point2 p{19008698.0, 33505401.0};

    REQUIRE(point_in_ring<DefaultKernel>(r, p) == PointInRing::Boundary);
    REQUIRE(point_in_ring<FastKernel>(r, p) == PointInRing::Inside);
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
    REQUIRE_NOTHROW(PointRing{as_span(legal)});

    const std::vector<Point2> turned{
        Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{4.0, 4.0}, Point2{0.0, 0.0},
    };
    REQUIRE_THROWS_AS(PointRing{as_span(turned)}, std::invalid_argument);
}
