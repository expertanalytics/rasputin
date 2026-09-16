// Unit tests for terrain::pred::DetriaExact -- the exact-arithmetic backend
// behind the default kernel, implemented on top of the vendored detria
// predicates.
//
// Three things are under test here, and only one of them is arithmetic:
//
//   1. The enum mapping. detria answers with its own enumerations,
//      `math::Orientation{CW=0, CCW=1, Collinear=2}` and
//      `math::CircleLocation{Inside=0, Outside=1, Cocircular=2}`. Neither the
//      values nor their order match ours (Clockwise=-1, Collinear=0,
//      CounterClockwise=1; Outside=-1, Cocircular=0, Inside=1), so the
//      translation must be an explicit switch and a `static_cast` is silently
//      wrong. A cast would turn detria's CW into our Collinear and its Inside
//      into our Cocircular -- both of which are plausible-looking answers that
//      a triangulation would act on. The mapping is therefore pinned value by
//      value, in both directions.
//
//   2. Agreement with `RefExact` across the adversarial families. RefExact is
//      an independently written oracle over integer coordinates, so this is the
//      strongest check available: two implementations sharing a bug would have
//      to share it by coincidence.
//
//   3. The shape of the type: stateless, static, header-only at the call site,
//      and -- load-bearing for the I/O and build boundaries -- not dragging
//      detria.hpp into every translation unit that asks a geometric question.
//
// Note on calling convention: every `incircle_ccw` call below passes a triple
// whose counterclockwise orientation is established first, by construction or
// via `as_ccw`. That is not tidiness. `incircle_ccw` inherits detria's Debug
// assertion, which raises SIGTRAP; a test that violated the precondition would
// kill the asan+ubsan CI job instead of failing.

#include <catch2/catch_test_macros.hpp>

#include <ccw_normalize.hpp>
#include <exact_reference.hpp>
#include <point_families.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/detria_exact.hpp>
#include <terrain/predicates/exact.hpp>
#include <terrain/predicates/orientation.hpp>

// The vendored detria.hpp is an implementation detail of one translation unit,
// src/predicates/detria_exact.cpp. It is a 4500-line header that pulls in
// <iostream>, <sstream> and <csignal>, and whose Debug assertions raise
// signals; nothing in the engine, and nothing in this suite, should compile
// against it just to ask which way three points turn.
//
// If detria_exact.hpp ever includes it, the guard macro from the pinned version
// is defined here and this stops the build with a readable message rather than
// a mysterious growth in compile time. The stronger, structural half of the
// same rule is in CMake: the test targets never get lib/ on their include path,
// so such an include would also simply fail to resolve.
#ifdef DETRIA_HPP_INCLUDED
#error "terrain/predicates/detria_exact.hpp must not include detria.hpp -- see src/predicates/detria_exact.cpp"
#endif

#include <array>
#include <cstddef>
#include <format>
#include <random>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

using terrain::Point2;
using terrain::pred::DetriaExact;
using terrain::pred::ExactPredicates;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::test::as_ccw;
using terrain::test::RefExact;
using terrain::test::utm33_offset;

namespace {

std::string describe(const Point2& a, const Point2& b, const Point2& c) {
    return std::format("a={} b={} c={}", a, b, c);
}

std::string describe(const Point2& a, const Point2& b, const Point2& c, const Point2& d) {
    return std::format("a={} b={} c={} d={}", a, b, c, d);
}

}  // namespace

// ---------------------------------------------------------------------------
// The type's shape.
// ---------------------------------------------------------------------------

TEST_CASE("DetriaExact models ExactPredicates", "[predicates][detria]") {
    STATIC_REQUIRE(ExactPredicates<DetriaExact>);
}

// The same structural purity check the kernels get. `ExactPredicates` is
// spelled with qualified static calls so a model must be usable without an
// instance, but the concept cannot pin emptiness: a backend could be
// static-callable and still carry a mutable cache, which is what would make the
// concurrency test below start failing intermittently.
TEST_CASE("DetriaExact is stateless and its predicates are static", "[predicates][detria][purity]") {
    STATIC_REQUIRE(std::is_empty_v<DetriaExact>);
    STATIC_REQUIRE(std::is_function_v<std::remove_pointer_t<decltype(&DetriaExact::orient2d)>>);
    STATIC_REQUIRE(std::is_function_v<std::remove_pointer_t<decltype(&DetriaExact::incircle_ccw)>>);
}

// ---------------------------------------------------------------------------
// The enum mapping, value by value.
// ---------------------------------------------------------------------------

// detria's CW is 0 and its Collinear is 2. `static_cast<Orientation>` on those
// yields Collinear and an enumerator we do not define -- so a clockwise triple
// would be reported as degenerate, and a genuinely degenerate one as a value no
// switch in the engine handles. Both of these are named here so the failure
// says which direction of the mapping broke.
TEST_CASE("DetriaExact::orient2d maps every detria orientation to ours", "[predicates][detria][mapping]") {
    const Point2 a{0, 0};
    const Point2 b{4, 0};
    const Point2 ccw{0, 3};
    const Point2 cw{0, -3};
    const Point2 collinear{9, 0};

    SECTION("counterclockwise") {
        REQUIRE(DetriaExact::orient2d(a, b, ccw) == Orientation::CounterClockwise);
    }
    SECTION("clockwise: detria CW is 0, which a cast would read as Collinear") {
        REQUIRE(DetriaExact::orient2d(a, b, cw) == Orientation::Clockwise);
    }
    SECTION("collinear: detria Collinear is 2, which is not an Orientation at all") {
        REQUIRE(DetriaExact::orient2d(a, b, collinear) == Orientation::Collinear);
    }
    SECTION("the three answers are pairwise distinct") {
        // Catches a mapping that collapses two detria values onto one of ours
        // -- a missing `case` falling through to a shared `default`, say --
        // which the three assertions above would not all detect on their own.
        const Orientation l = DetriaExact::orient2d(a, b, ccw);
        const Orientation r = DetriaExact::orient2d(a, b, cw);
        const Orientation z = DetriaExact::orient2d(a, b, collinear);
        REQUIRE(l != r);
        REQUIRE(l != z);
        REQUIRE(r != z);
    }
}

// The degenerate orientations detria still has to classify: a repeated vertex
// and a fully coincident triple. Both are Collinear, and both are ordinary
// inputs for a snapped PSLG rather than exotica.
TEST_CASE("DetriaExact::orient2d classifies coincident vertices as Collinear", "[predicates][detria][mapping][adversarial]") {
    REQUIRE(DetriaExact::orient2d(Point2{5, 5}, Point2{5, 5}, Point2{9, 1}) == Orientation::Collinear);
    REQUIRE(DetriaExact::orient2d(Point2{5, 5}, Point2{9, 1}, Point2{5, 5}) == Orientation::Collinear);
    REQUIRE(DetriaExact::orient2d(Point2{9, 1}, Point2{5, 5}, Point2{5, 5}) == Orientation::Collinear);
    REQUIRE(DetriaExact::orient2d(Point2{5, 5}, Point2{5, 5}, Point2{5, 5}) == Orientation::Collinear);

    const Point2 p = utm33_offset({0, 0});
    REQUIRE(DetriaExact::orient2d(p, p, utm33_offset({1, 1})) == Orientation::Collinear);
}

// detria's CircleLocation is Inside=0, Outside=1, Cocircular=2 against our
// Outside=-1, Cocircular=0, Inside=1. A cast maps Inside to Cocircular and
// Outside to Inside -- that is, it reports a point inside the circumcircle as
// on it and a point outside it as inside, which in a Delaunay flip loop is the
// difference between terminating and not.
//
// Every triple here is counterclockwise by construction; see the file header.
TEST_CASE("DetriaExact::incircle_ccw maps every detria circle location to ours", "[predicates][detria][mapping]") {
    const Point2 a{5, 0};
    const Point2 b{0, 5};
    const Point2 c{-5, 0};
    REQUIRE(RefExact::orient2d(a, b, c) == Orientation::CounterClockwise);

    SECTION("inside: detria Inside is 0, which a cast would read as Cocircular") {
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{0, 0}) == Incircle::Inside);
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{1, 1}) == Incircle::Inside);
    }
    SECTION("outside: detria Outside is 1, which a cast would read as Inside") {
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{400, 400}) == Incircle::Outside);
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{6, 0}) == Incircle::Outside);
    }
    SECTION("cocircular: detria Cocircular is 2, which is not an Incircle at all") {
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{3, -4}) == Incircle::Cocircular);
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{-4, -3}) == Incircle::Cocircular);
    }
    SECTION("a fourth point coincident with a vertex is Cocircular") {
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, a) == Incircle::Cocircular);
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, b) == Incircle::Cocircular);
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, c) == Incircle::Cocircular);
    }
    SECTION("the three answers are pairwise distinct") {
        const Incircle in = DetriaExact::incircle_ccw(a, b, c, Point2{0, 0});
        const Incircle out = DetriaExact::incircle_ccw(a, b, c, Point2{400, 400});
        const Incircle on = DetriaExact::incircle_ccw(a, b, c, Point2{3, -4});
        REQUIRE(in != out);
        REQUIRE(in != on);
        REQUIRE(out != on);
    }
}

// The sign convention travels with the mapping. detria's incircle determinant
// is positive for Inside given counterclockwise input, as ours is -- but a
// backend written against the opposite convention would still pass every
// "pairwise distinct" check while inverting the predicate. Pinning Inside to a
// point demonstrably enclosed by the circle, and Outside to one demonstrably
// not, is what fixes the polarity.
TEST_CASE("DetriaExact::incircle_ccw has the same polarity as our Incircle", "[predicates][detria][mapping]") {
    // Unit circle through three lattice points, with the centre as the fourth.
    const Point2 a{1, 0};
    const Point2 b{0, 1};
    const Point2 c{-1, 0};
    REQUIRE(RefExact::orient2d(a, b, c) == Orientation::CounterClockwise);

    REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{0, 0}) == Incircle::Inside);
    REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{0, -1}) == Incircle::Cocircular);
    REQUIRE(DetriaExact::incircle_ccw(a, b, c, Point2{0, -2}) == Incircle::Outside);
}

// ---------------------------------------------------------------------------
// Agreement with the independent oracle.
// ---------------------------------------------------------------------------

TEST_CASE("DetriaExact::orient2d agrees with RefExact on the adversarial corpus", "[predicates][detria][adversarial]") {
    const std::vector<std::array<Point2, 3>> corpus = {
        {Point2{0, 0}, Point2{1, 0}, Point2{0, 1}},
        {Point2{0, 0}, Point2{0, 1}, Point2{1, 0}},
        {Point2{-3, 2}, Point2{4, 2}, Point2{1, -5}},
        {Point2{0, 0}, Point2{1, 1}, Point2{2, 2}},
        {Point2{0, 0}, Point2{1000000, 1000000}, Point2{3000000, 3000000}},
        {Point2{5, 5}, Point2{5, 5}, Point2{9, 1}},
        {Point2{5, 5}, Point2{5, 5}, Point2{5, 5}},
        // A sliver: two points a third of a unit apart in y across a span of
        // 2^25 in x.
        {Point2{0, 0}, Point2{16777216, 1}, Point2{33554432, 1}},
        {utm33_offset({0, 0}), utm33_offset({1, 0}), utm33_offset({0, 1})},
        {utm33_offset({0, 0}), utm33_offset({4, 4}), utm33_offset({9, 9})},
        // A narrow corridor at UTM33 magnitude: half a million metres long,
        // one metre wide.
        {utm33_offset({0, 0}), utm33_offset({1000000, 0}), utm33_offset({500000, 1})},
    };

    for (const auto& t : corpus) {
        INFO(describe(t[0], t[1], t[2]));
        REQUIRE(DetriaExact::orient2d(t[0], t[1], t[2]) == RefExact::orient2d(t[0], t[1], t[2]));
    }
}

TEST_CASE("DetriaExact::orient2d agrees with RefExact near degeneracy", "[predicates][detria][adversarial]") {
    std::mt19937_64 rng{424242};
    for (int i = 0; i < 500; ++i) {
        const auto t = terrain::test::near_degenerate_triple(rng);
        INFO(describe(t.a, t.b, t.c));
        REQUIRE(DetriaExact::orient2d(t.a, t.b, t.c) == t.expected);
        REQUIRE(DetriaExact::orient2d(t.a, t.b, t.c) == RefExact::orient2d(t.a, t.b, t.c));
    }
}

TEST_CASE("DetriaExact::orient2d confirms exact collinearity at UTM33 magnitude", "[predicates][detria][adversarial]") {
    std::mt19937_64 rng{99};
    for (int i = 0; i < 200; ++i) {
        const auto t = terrain::test::exactly_collinear_triple(rng);
        INFO(describe(t.a, t.b, t.c));
        REQUIRE(DetriaExact::orient2d(t.a, t.b, t.c) == Orientation::Collinear);
    }

    // A long collinear run on the integer lattice -- what a snap grid produces,
    // and the input a naive determinant is least able to confirm.
    for (int k = 2; k <= 64; ++k) {
        const Point2 a = utm33_offset({0, 0});
        const Point2 b = utm33_offset({3, 7});
        const Point2 c = utm33_offset({3.0 * k, 7.0 * k});
        INFO(describe(a, b, c));
        REQUIRE(DetriaExact::orient2d(a, b, c) == Orientation::Collinear);
    }
}

// A sub-millimetre feature inside a coordinate system spanning hundreds of
// kilometres: the everyday case for a TIN built from a snapped breakline
// network in UTM33. RefExact cannot referee this one -- its coordinates are not
// integers -- so the expected answer is known by construction instead: 2^-14 is
// exactly representable, so adding it to an integer easting loses nothing and
// the triangle really is a right angle of that size.
TEST_CASE("DetriaExact::orient2d resolves a sub-millimetre triangle at UTM33 magnitude", "[predicates][detria][adversarial]") {
    const double tiny = 0.00006103515625;  // about 0.06 mm
    const Point2 a = utm33_offset({0, 0});
    const Point2 b{a.x + tiny, a.y};
    const Point2 c{a.x, a.y + tiny};

    REQUIRE(DetriaExact::orient2d(a, b, c) == Orientation::CounterClockwise);
    REQUIRE(DetriaExact::orient2d(a, c, b) == Orientation::Clockwise);
    REQUIRE(DetriaExact::orient2d(a, b, Point2{a.x + 2 * tiny, a.y}) == Orientation::Collinear);
}

TEST_CASE("DetriaExact::incircle_ccw agrees with RefExact on an exact circle", "[predicates][detria][adversarial]") {
    const auto circle = terrain::test::integer_circle_radius_5();
    const std::array<Point2, 4> probes = {Point2{0, 0}, Point2{400, 400}, Point2{1, 1}, Point2{5, 1}};

    for (std::size_t i = 0; i < circle.size(); ++i) {
        const auto triple = as_ccw(circle[i], circle[(i + 3) % circle.size()],
                                   circle[(i + 6) % circle.size()]);
        REQUIRE(triple.has_value());
        const auto& [a, b, c] = *triple;

        // The fourth point taken from the circle itself: exactly cocircular,
        // Pythagorean rather than approximate, so a Cocircular answer here is
        // right for the right reason.
        const Point2 on = circle[(i + 9) % circle.size()];
        INFO(describe(a, b, c, on));
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, on) == Incircle::Cocircular);
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, on) == RefExact::incircle_ccw(a, b, c, on));

        for (const Point2& d : probes) {
            INFO(describe(a, b, c, d));
            REQUIRE(DetriaExact::incircle_ccw(a, b, c, d) == RefExact::incircle_ccw(a, b, c, d));
        }
    }
}

TEST_CASE("DetriaExact::incircle_ccw agrees with RefExact at UTM33 magnitude", "[predicates][detria][adversarial]") {
    const auto circle = terrain::test::integer_circle_radius_5();

    for (std::size_t i = 0; i < circle.size(); ++i) {
        const auto triple = as_ccw(utm33_offset(circle[i]),
                                   utm33_offset(circle[(i + 3) % circle.size()]),
                                   utm33_offset(circle[(i + 6) % circle.size()]));
        REQUIRE(triple.has_value());
        const auto& [a, b, c] = *triple;

        for (const Point2& d : {utm33_offset(circle[(i + 9) % circle.size()]),
                                utm33_offset({0, 0}),
                                utm33_offset({40, 40}),
                                utm33_offset({1, 1})}) {
            INFO(describe(a, b, c, d));
            REQUIRE(DetriaExact::incircle_ccw(a, b, c, d) == RefExact::incircle_ccw(a, b, c, d));
        }
    }
}

// Clusters of points *nearly* on a circle, at a scale where the naive lifted
// determinant has no significant digits left. The oracle settles each case; the
// point of the family is that the answers are definite, not that they are
// predictable by inspection.
TEST_CASE("DetriaExact::incircle_ccw agrees with RefExact on near-cocircular clusters", "[predicates][detria][adversarial]") {
    const auto circle = terrain::test::integer_circle_radius_5();
    std::mt19937_64 rng{31337};
    std::uniform_int_distribution<int> nudge{-1, 1};
    std::uniform_int_distribution<int> scale{1000, 4000};

    for (int trial = 0; trial < 200; ++trial) {
        const int s = scale(rng);
        const auto blow_up = [s, &rng, &nudge](const Point2& p) {
            return Point2{p.x * s + nudge(rng), p.y * s + nudge(rng)};
        };

        const auto triple = as_ccw(blow_up(circle[0]), blow_up(circle[4]), blow_up(circle[8]));
        if (!triple.has_value()) {
            continue;  // Perturbed into collinearity; no CCW ordering exists.
        }
        const auto& [a, b, c] = *triple;
        const Point2 d = blow_up(circle[2]);

        INFO(describe(a, b, c, d));
        REQUIRE(DetriaExact::incircle_ccw(a, b, c, d) == RefExact::incircle_ccw(a, b, c, d));
    }
}

// ---------------------------------------------------------------------------
// Purity under concurrency.
// ---------------------------------------------------------------------------

// detria's robust predicates evaluate floating-point expansions; the ones we
// call must do so on the stack, with no shared scratch buffer and no lazily
// initialised table. If a future version of the vendored header -- or our
// wrapper -- acquires either, the engine's per-subdomain parallel refinement
// would corrupt silently rather than fail, so the property is asserted rather
// than assumed.
TEST_CASE("DetriaExact is callable concurrently", "[predicates][detria][purity]") {
    std::mt19937_64 rng{2718};
    std::vector<terrain::test::KnownTriple> work;
    work.reserve(256);
    for (int i = 0; i < 256; ++i) {
        work.push_back(i % 2 == 0 ? terrain::test::near_degenerate_triple(rng)
                                  : terrain::test::exactly_collinear_triple(rng));
    }

    constexpr int thread_count = 4;
    std::vector<std::thread> threads;
    std::vector<int> mismatches(thread_count, 0);

    for (int t = 0; t < thread_count; ++t) {
        threads.emplace_back([&work, &mismatches, t] {
            for (int repeat = 0; repeat < 8; ++repeat) {
                for (const auto& triple : work) {
                    if (DetriaExact::orient2d(triple.a, triple.b, triple.c) != triple.expected) {
                        ++mismatches[static_cast<std::size_t>(t)];
                    }
                }
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }

    for (const int m : mismatches) {
        REQUIRE(m == 0);
    }
}
