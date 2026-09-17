// Unit tests for terrain::Box2 and the free terrain::bounding_box.
//
// The load-bearing decision in this header is that a default-constructed Box2
// is *empty* -- lo at +inf, hi at -inf -- rather than a degenerate box at the
// origin. That single choice is what makes
//
//     Box2 b; for (const Point2& p : pts) b.expand(p);
//
// a correct fold: an empty input yields an empty box instead of one that
// silently contains (0, 0), and no caller has to special-case the first
// point. Most of this file exists to pin the consequences of that choice --
// empty is contained by everything, empty intersects nothing, and expanding
// by an empty box is a no-op.
//
// There are no tolerances here and none in the header. Every comparison
// Box2 makes is between coordinates the caller supplied.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <point_families.hpp>

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>

#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <span>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <string>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using terrain::Box2;
using terrain::Point2;
using terrain::bounding_box;
using terrain::test::utm33_offset;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();
const double quiet_nan = std::numeric_limits<double>::quiet_NaN();

[[nodiscard]] Box2 unit_box() { return Box2{Point2{0.0, 0.0}, Point2{1.0, 1.0}}; }

[[nodiscard]] std::span<const Point2> as_span(const std::vector<Point2>& v) {
    return std::span<const Point2>{v};
}

}  // namespace

// ---------------------------------------------------------------------------
// The empty box
// ---------------------------------------------------------------------------

TEST_CASE("a default-constructed Box2 is empty, not a box at the origin", "[bbox][empty]") {
    const Box2 b;

    REQUIRE(b.is_empty());
    REQUIRE_FALSE(b.contains(Point2{0.0, 0.0}));
    REQUIRE(b.width() == 0.0);
    REQUIRE(b.height() == 0.0);
}

TEST_CASE("the empty box carries the sentinel corners the fold relies on", "[bbox][empty]") {
    const Box2 b;

    REQUIRE(b.lo() == Point2{inf, inf});
    REQUIRE(b.hi() == Point2{-inf, -inf});
}

TEST_CASE("two empty boxes compare equal", "[bbox][empty]") {
    REQUIRE(Box2{} == Box2{});
    REQUIRE_FALSE(Box2{} == unit_box());
}

// ---------------------------------------------------------------------------
// The two-point constructor validates; expand does not
// ---------------------------------------------------------------------------

TEST_CASE("Box2's two-point constructor accepts a well-ordered pair", "[bbox][ctor]") {
    const Box2 b{Point2{-3.0, 2.0}, Point2{4.5, 9.0}};

    REQUIRE_FALSE(b.is_empty());
    REQUIRE(b.lo() == Point2{-3.0, 2.0});
    REQUIRE(b.hi() == Point2{4.5, 9.0});
    REQUIRE(b.width() == 7.5);
    REQUIRE(b.height() == 7.0);
}

TEST_CASE("Box2 rejects lo > hi in either component independently", "[bbox][ctor]") {
    REQUIRE_THROWS_AS((Box2{Point2{1.0, 0.0}, Point2{0.0, 1.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS((Box2{Point2{0.0, 1.0}, Point2{1.0, 0.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS((Box2{Point2{1.0, 1.0}, Point2{0.0, 0.0}}), std::invalid_argument);
}

// NaN is the case a bare `lo > hi` check misses: every comparison against NaN
// is false, so a range test alone accepts it and the box then answers every
// containment query with a confident `false`. The finiteness test has to be
// explicit, exactly as RasterGeometry::cell_of documents for the same reason.
TEST_CASE("Box2 rejects non-finite corners", "[bbox][ctor][nan]") {
    REQUIRE_THROWS_AS((Box2{Point2{quiet_nan, 0.0}, Point2{1.0, 1.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS((Box2{Point2{0.0, quiet_nan}, Point2{1.0, 1.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS((Box2{Point2{0.0, 0.0}, Point2{quiet_nan, 1.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS((Box2{Point2{0.0, 0.0}, Point2{1.0, quiet_nan}}), std::invalid_argument);
    REQUIRE_THROWS_AS((Box2{Point2{-inf, 0.0}, Point2{1.0, 1.0}}), std::invalid_argument);
    REQUIRE_THROWS_AS((Box2{Point2{0.0, 0.0}, Point2{inf, 1.0}}), std::invalid_argument);
}

// The empty box is reachable only through the default constructor. Handing the
// sentinel corners back to the two-point constructor is an error, not a way to
// spell `Box2{}`.
TEST_CASE("Box2's constructor will not rebuild the empty sentinel", "[bbox][ctor][empty]") {
    REQUIRE_THROWS_AS((Box2{Point2{inf, inf}, Point2{-inf, -inf}}), std::invalid_argument);
}

TEST_CASE("a degenerate Box2 at a single point is not empty", "[bbox][ctor][degenerate]") {
    const Point2 p{7.25, -3.5};
    const Box2 b{p, p};

    REQUIRE_FALSE(b.is_empty());
    REQUIRE(b.width() == 0.0);
    REQUIRE(b.height() == 0.0);
    REQUIRE(b.contains(p));
    REQUIRE(b.center() == p);
}

// A box may be zero-extent in one axis only -- a horizontal or vertical
// sliver. This is what the bounding box of an axis-aligned segment is, so it
// has to be a first-class value rather than a rejected degeneracy.
TEST_CASE("a Box2 may be degenerate in one axis only", "[bbox][ctor][degenerate]") {
    const Box2 horizontal{Point2{0.0, 5.0}, Point2{4.0, 5.0}};

    REQUIRE_FALSE(horizontal.is_empty());
    REQUIRE(horizontal.width() == 4.0);
    REQUIRE(horizontal.height() == 0.0);
    REQUIRE(horizontal.contains(Point2{2.0, 5.0}));
    REQUIRE_FALSE(horizontal.contains(Point2{2.0, 5.5}));
}

// ---------------------------------------------------------------------------
// expand
// ---------------------------------------------------------------------------

// expand is the inner loop of every bounding-box computation in the engine, so
// it validates nothing and is noexcept. Finiteness is the caller's
// precondition, established once at the PSLG and Python boundaries.
TEST_CASE("expand is noexcept in both overloads", "[bbox][expand]") {
    Box2 b;
    const Point2 p{1.0, 1.0};
    const Box2 o = unit_box();

    STATIC_REQUIRE(noexcept(b.expand(p)));
    STATIC_REQUIRE(noexcept(b.expand(o)));
}

TEST_CASE("expanding the empty box by a point yields the degenerate box there", "[bbox][expand]") {
    Box2 b;
    const Point2 p{-2.0, 6.0};
    b.expand(p);

    REQUIRE_FALSE(b.is_empty());
    REQUIRE(b.lo() == p);
    REQUIRE(b.hi() == p);
}

TEST_CASE("expand grows each axis independently", "[bbox][expand]") {
    Box2 b = unit_box();
    b.expand(Point2{0.5, 4.0});

    REQUIRE(b.lo() == Point2{0.0, 0.0});
    REQUIRE(b.hi() == Point2{1.0, 4.0});
}

TEST_CASE("expanding by an already-contained point changes nothing", "[bbox][expand]") {
    Box2 b = unit_box();
    b.expand(Point2{0.5, 0.5});
    b.expand(Point2{0.0, 1.0});  // on the boundary

    REQUIRE(b == unit_box());
}

// Expanding by the empty box must be a no-op, which is what makes a fold over
// sub-boxes correct when some of them are empty. Note that the empty box's own
// corners are non-finite, so this overload's finiteness precondition has to
// admit the empty box as an argument.
TEST_CASE("expanding by an empty box is a no-op in both directions", "[bbox][expand][empty]") {
    Box2 b = unit_box();
    b.expand(Box2{});
    REQUIRE(b == unit_box());

    Box2 e;
    e.expand(unit_box());
    REQUIRE(e == unit_box());

    Box2 both;
    both.expand(Box2{});
    REQUIRE(both.is_empty());
}

TEST_CASE("expand by a box takes the componentwise union", "[bbox][expand]") {
    Box2 b{Point2{0.0, 0.0}, Point2{1.0, 1.0}};
    b.expand(Box2{Point2{-4.0, 0.25}, Point2{0.5, 3.0}});

    REQUIRE(b == Box2{Point2{-4.0, 0.0}, Point2{1.0, 3.0}});
}

// ---------------------------------------------------------------------------
// contains and intersects, closed and exact
// ---------------------------------------------------------------------------

TEST_CASE("contains(Point2) is closed on every edge and corner", "[bbox][contains]") {
    const Box2 b{Point2{0.0, 0.0}, Point2{2.0, 4.0}};

    REQUIRE(b.contains(Point2{0.0, 0.0}));
    REQUIRE(b.contains(Point2{2.0, 4.0}));
    REQUIRE(b.contains(Point2{0.0, 4.0}));
    REQUIRE(b.contains(Point2{2.0, 0.0}));
    REQUIRE(b.contains(Point2{1.0, 0.0}));
    REQUIRE(b.contains(Point2{0.0, 2.0}));
    REQUIRE(b.contains(Point2{1.0, 2.0}));

    REQUIRE_FALSE(b.contains(Point2{-0.5, 2.0}));
    REQUIRE_FALSE(b.contains(Point2{1.0, 4.5}));
}

// One ulp outside is outside. There is no tolerance in this header, so the
// classification flips at the exact representable neighbour of the corner.
TEST_CASE("contains(Point2) is exact to the ulp at UTM33 magnitudes", "[bbox][contains][utm33]") {
    const Point2 lo = utm33_offset(Point2{0.0, 0.0});
    const Point2 hi = utm33_offset(Point2{1000.0, 1000.0});
    const Box2 b{lo, hi};

    REQUIRE(b.contains(hi));
    REQUIRE_FALSE(b.contains(Point2{std::nextafter(hi.x, inf), hi.y}));
    REQUIRE_FALSE(b.contains(Point2{hi.x, std::nextafter(hi.y, inf)}));
    REQUIRE(b.contains(Point2{std::nextafter(hi.x, -inf), hi.y}));
}

TEST_CASE("contains(Point2) rejects a non-finite query", "[bbox][contains][nan]") {
    const Box2 b = unit_box();

    REQUIRE_FALSE(b.contains(Point2{quiet_nan, 0.5}));
    REQUIRE_FALSE(b.contains(Point2{0.5, quiet_nan}));
    REQUIRE_FALSE(b.contains(Point2{inf, 0.5}));
}

TEST_CASE("the empty box contains nothing", "[bbox][contains][empty]") {
    const Box2 e;

    REQUIRE_FALSE(e.contains(Point2{0.0, 0.0}));
    REQUIRE_FALSE(e.contains(Point2{1e9, -1e9}));
    REQUIRE_FALSE(e.contains(unit_box()));
}

TEST_CASE("contains(Box2) is closed and reflexive", "[bbox][contains]") {
    const Box2 outer{Point2{0.0, 0.0}, Point2{4.0, 4.0}};

    REQUIRE(outer.contains(outer));
    REQUIRE(outer.contains(Box2{Point2{1.0, 1.0}, Point2{2.0, 2.0}}));
    REQUIRE(outer.contains(Box2{Point2{0.0, 0.0}, Point2{4.0, 4.0}}));
    REQUIRE(outer.contains(Box2{Point2{0.0, 2.0}, Point2{0.0, 2.0}}));  // on the edge
    REQUIRE_FALSE(outer.contains(Box2{Point2{-0.5, 1.0}, Point2{2.0, 2.0}}));
    REQUIRE_FALSE(outer.contains(Box2{Point2{1.0, 1.0}, Point2{5.0, 2.0}}));
}

// The empty box is the identity of the union, so it is a subset of everything
// -- including of itself.
TEST_CASE("the empty box is contained by every box", "[bbox][contains][empty]") {
    REQUIRE(unit_box().contains(Box2{}));
    REQUIRE(Box2{}.contains(Box2{}));
}

TEST_CASE("intersects is closed: boxes touching at an edge or corner intersect", "[bbox][intersects]") {
    const Box2 b{Point2{0.0, 0.0}, Point2{1.0, 1.0}};

    REQUIRE(b.intersects(Box2{Point2{1.0, 0.0}, Point2{2.0, 1.0}}));   // shared edge
    REQUIRE(b.intersects(Box2{Point2{1.0, 1.0}, Point2{2.0, 2.0}}));   // shared corner
    REQUIRE(b.intersects(Box2{Point2{1.0, 1.0}, Point2{1.0, 1.0}}));   // degenerate at the corner
    REQUIRE(b.intersects(b));
    REQUIRE(b.intersects(Box2{Point2{0.25, 0.25}, Point2{0.5, 0.5}})); // containment implies it
}

TEST_CASE("intersects is symmetric and false for disjoint boxes", "[bbox][intersects]") {
    const Box2 a{Point2{0.0, 0.0}, Point2{1.0, 1.0}};
    const Box2 b{Point2{1.5, 0.0}, Point2{2.0, 1.0}};
    const Box2 c{Point2{0.0, 2.0}, Point2{1.0, 3.0}};  // disjoint in y only

    REQUIRE_FALSE(a.intersects(b));
    REQUIRE_FALSE(b.intersects(a));
    REQUIRE_FALSE(a.intersects(c));
    REQUIRE_FALSE(c.intersects(a));
}

TEST_CASE("the empty box intersects nothing, including itself", "[bbox][intersects][empty]") {
    REQUIRE_FALSE(Box2{}.intersects(unit_box()));
    REQUIRE_FALSE(unit_box().intersects(Box2{}));
    REQUIRE_FALSE(Box2{}.intersects(Box2{}));
}

// ---------------------------------------------------------------------------
// center
// ---------------------------------------------------------------------------

TEST_CASE("center is the midpoint of a non-empty box", "[bbox][center]") {
    REQUIRE(Box2{Point2{0.0, 0.0}, Point2{2.0, 4.0}}.center() == Point2{1.0, 2.0});
    REQUIRE(Box2{Point2{-3.0, -1.0}, Point2{1.0, 1.0}}.center() == Point2{-1.0, 0.0});
}

// `(lo + hi) / 2` overflows to infinity on a box this wide while
// `lo / 2 + hi / 2` does not. The corners are finite and the box is legal, so
// a finite centre is the only defensible answer.
TEST_CASE("center does not overflow on an extreme-scale box", "[bbox][center][extreme]") {
    // Deliberately asymmetric: lo + hi overflows to +inf here, where a box
    // centred on the origin would have cancelled and hidden the bug.
    const Box2 b{Point2{1.0e308, 1.0e308}, Point2{1.5e308, 1.5e308}};
    const Point2 c = b.center();

    REQUIRE(std::isfinite(c.x));
    REQUIRE(std::isfinite(c.y));
    REQUIRE(c == Point2{1.25e308, 1.25e308});
    REQUIRE(b.contains(c));

    const Box2 symmetric{Point2{-1e308, -1e308}, Point2{1e308, 1e308}};
    REQUIRE(symmetric.center() == Point2{0.0, 0.0});
}

TEST_CASE("center of a box lies inside it", "[bbox][center]") {
    const Box2 b{utm33_offset(Point2{0.0, 0.0}), utm33_offset(Point2{1000.0, 250.0})};

    REQUIRE(b.contains(b.center()));
}

// ---------------------------------------------------------------------------
// bounding_box over a span
// ---------------------------------------------------------------------------

TEST_CASE("bounding_box of an empty span is the empty box", "[bbox][bounding_box]") {
    const std::vector<Point2> none;

    REQUIRE(bounding_box(as_span(none)).is_empty());
}

TEST_CASE("bounding_box of one point is degenerate there", "[bbox][bounding_box]") {
    const std::vector<Point2> one{Point2{3.0, -4.0}};
    const Box2 b = bounding_box(as_span(one));

    REQUIRE_FALSE(b.is_empty());
    REQUIRE(b.lo() == Point2{3.0, -4.0});
    REQUIRE(b.hi() == Point2{3.0, -4.0});
}

TEST_CASE("bounding_box spans every point", "[bbox][bounding_box]") {
    const std::vector<Point2> pts{
        Point2{1.0, 1.0}, Point2{-2.0, 5.0}, Point2{4.0, 0.0}, Point2{0.0, -3.0},
    };
    const Box2 b = bounding_box(as_span(pts));

    REQUIRE(b == Box2{Point2{-2.0, -3.0}, Point2{4.0, 5.0}});
    for (const Point2& p : pts) {
        REQUIRE(b.contains(p));
    }
}

// Unlike expand, the free function is a boundary: it is what a caller reaches
// for with data of unknown provenance, so it validates.
TEST_CASE("bounding_box rejects a non-finite point", "[bbox][bounding_box][nan]") {
    const std::vector<Point2> with_nan{Point2{0.0, 0.0}, Point2{quiet_nan, 1.0}};
    const std::vector<Point2> with_inf{Point2{0.0, 0.0}, Point2{1.0, inf}};

    REQUIRE_THROWS_AS(bounding_box(as_span(with_nan)), std::invalid_argument);
    REQUIRE_THROWS_AS(bounding_box(as_span(with_inf)), std::invalid_argument);
}

TEST_CASE("bounding_box is exact at UTM33 magnitudes", "[bbox][bounding_box][utm33]") {
    const std::vector<Point2> pts{
        utm33_offset(Point2{0.0, 0.0}),
        utm33_offset(Point2{0.001, 0.002}),
        utm33_offset(Point2{-0.001, 0.0}),
    };
    const Box2 b = bounding_box(as_span(pts));

    REQUIRE(b.lo() == utm33_offset(Point2{-0.001, 0.0}));
    REQUIRE(b.hi() == Point2{utm33_offset(Point2{0.001, 0.0}).x,
                             utm33_offset(Point2{0.0, 0.002}).y});
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

TEST_CASE("Box2 formats through std::format", "[bbox][format]") {
    const Box2 b{Point2{0.0, 1.0}, Point2{2.0, 3.0}};
    const std::string s = std::format("{}", b);

    REQUIRE_THAT(s, ContainsSubstring("Box2"));
    REQUIRE_THAT(s, ContainsSubstring(std::format("{}", b.lo())));
    REQUIRE_THAT(s, ContainsSubstring(std::format("{}", b.hi())));
}

// The empty box's corners are infinities. Printing them would hand a reader
// "Box2(Point2(inf, inf), Point2(-inf, -inf))", which reads as a bug report
// rather than as the identity element it is.
TEST_CASE("the empty Box2 renders distinctly and without infinities", "[bbox][format][empty]") {
    const std::string s = std::format("{}", Box2{});

    REQUIRE_THAT(s, ContainsSubstring("Box2"));
    REQUIRE_THAT(s, !ContainsSubstring("inf"));
    REQUIRE(s != std::format("{}", unit_box()));
}

// Same reasoning as point.hpp: a spec that parses and is then discarded is a
// silent lie about the output.
TEST_CASE("Box2's formatter rejects a spec it does not honour", "[bbox][format]") {
    const Box2 b = unit_box();

    REQUIRE_THROWS_AS(std::vformat("{:>24}", std::make_format_args(b)), std::format_error);
    REQUIRE_THROWS_AS(std::vformat("{:.3f}", std::make_format_args(b)), std::format_error);
}

// ---------------------------------------------------------------------------
// expand as a fold
// ---------------------------------------------------------------------------

// min and max are associative and commutative, so the fold is order
// independent -- exactly, not approximately, since no arithmetic happens. This
// is what lets a future parallel bounding-box computation split the point
// buffer across threads and combine the partial boxes in whatever order they
// finish in.
TEST_CASE("expand as a fold is order independent", "[bbox][expand][property]") {
    const int seed = GENERATE(range(0, 32));
    std::mt19937_64 rng{0xb0c50000ULL + static_cast<unsigned long long>(seed)};
    std::uniform_real_distribution<double> coord{-1e4, 1e4};

    std::vector<Point2> pts;
    pts.reserve(24);
    for (int i = 0; i < 24; ++i) pts.push_back(Point2{coord(rng), coord(rng)});

    Box2 reference;
    for (const Point2& p : pts) reference.expand(p);

    for (int trial = 0; trial < 8; ++trial) {
        std::shuffle(pts.begin(), pts.end(), rng);
        Box2 b;
        for (const Point2& p : pts) b.expand(p);
        INFO(std::format("seed={} trial={} box={}", seed, trial, b));
        REQUIRE(b == reference);
    }
}

// The same law for the box overload, which is how partial results combine.
TEST_CASE("folding sub-boxes matches folding their points", "[bbox][expand][property]") {
    const int seed = GENERATE(range(0, 32));
    std::mt19937_64 rng{0xb0c51000ULL + static_cast<unsigned long long>(seed)};
    std::uniform_real_distribution<double> coord{-1e4, 1e4};

    std::vector<Point2> pts;
    pts.reserve(30);
    for (int i = 0; i < 30; ++i) pts.push_back(Point2{coord(rng), coord(rng)});

    Box2 whole;
    for (const Point2& p : pts) whole.expand(p);

    // Three chunks, one of them deliberately left empty.
    Box2 first, second, third;
    for (std::size_t i = 0; i < 12; ++i) first.expand(pts[i]);
    for (std::size_t i = 12; i < 30; ++i) second.expand(pts[i]);
    REQUIRE(third.is_empty());

    Box2 combined;
    combined.expand(third);
    combined.expand(first);
    combined.expand(second);

    INFO(std::format("seed={}", seed));
    REQUIRE(combined == whole);
}
