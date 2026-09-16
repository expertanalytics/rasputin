// Unit tests for terrain::pred::FilteredKernel -- the kernel the rest of the
// engine is expected to use.
//
// The filter is a performance device. It may decide *how* an answer is reached
// and never *what* the answer is, so the headline test in this file is that
// FilteredKernel<E> agrees with E on every input it is given. Everything else
// here either supports that claim (the filter is really being exercised, in
// both directions) or covers the two things the filter is not allowed to
// inherit from its backend: a precondition, and an ordering sensitivity.

#include <catch2/catch_test_macros.hpp>

#include <exact_reference.hpp>
#include <point_families.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/exact.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <array>
#include <format>
#include <limits>
#include <random>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

using terrain::Point2;
using terrain::pred::ExactPredicates;
using terrain::pred::FilteredKernel;
using terrain::pred::GeometryKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::pred::incircle_bound_a;
using terrain::pred::orient2d_bound_a;
using terrain::pred::reversed;
using terrain::test::CcwCheckingExact;
using terrain::test::CountingExact;
using terrain::test::RefExact;
using terrain::test::utm33_offset;

namespace {

using Filtered = FilteredKernel<RefExact>;
using Counted = FilteredKernel<CountingExact<RefExact>>;
using CcwChecked = FilteredKernel<CcwCheckingExact<RefExact>>;

std::string describe(const Point2& a, const Point2& b, const Point2& c) {
    return std::format("a={} b={} c={}", a, b, c);
}

std::string describe(const Point2& a, const Point2& b, const Point2& c, const Point2& d) {
    return std::format("a={} b={} c={} d={}", a, b, c, d);
}

// Triples spanning every family the kernel must survive: general position,
// duplicated and coincident vertices, exact collinearity at unit and at UTM33
// magnitude, sliver-thin triangles, and a sub-millimetre feature embedded in a
// coordinate system spanning hundreds of kilometres.
const std::vector<std::array<Point2, 3>>& orientation_corpus() {
    static const std::vector<std::array<Point2, 3>> corpus = {
        {Point2{0, 0}, Point2{1, 0}, Point2{0, 1}},
        {Point2{0, 0}, Point2{0, 1}, Point2{1, 0}},
        {Point2{-3, 2}, Point2{4, 2}, Point2{1, -5}},
        {Point2{0, 0}, Point2{1, 1}, Point2{2, 2}},
        {Point2{0, 0}, Point2{1000000, 1000000}, Point2{3000000, 3000000}},
        {Point2{5, 5}, Point2{5, 5}, Point2{9, 1}},
        {Point2{5, 5}, Point2{5, 5}, Point2{5, 5}},
        {Point2{0, 0}, Point2{16777216, 1}, Point2{33554432, 1}},
        {utm33_offset({0, 0}), utm33_offset({1, 0}), utm33_offset({0, 1})},
        {utm33_offset({0, 0}), utm33_offset({4, 4}), utm33_offset({9, 9})},
        {utm33_offset({0, 0}), utm33_offset({1000000, 0}), utm33_offset({500000, 1})},
    };
    return corpus;
}

}  // namespace

// ---------------------------------------------------------------------------
// Concepts and the static filter bounds.
// ---------------------------------------------------------------------------

TEST_CASE("RefExact models ExactPredicates and FilteredKernel models GeometryKernel", "[predicates][filtered_kernel]") {
    STATIC_REQUIRE(ExactPredicates<RefExact>);
    STATIC_REQUIRE(GeometryKernel<Filtered>);
}

// The bounds are Shewchuk's static error bounds for the two determinants,
// expressed in half-ulp units. They are part of the interface because the
// filter's correctness argument rests on them, and a typo in a magic constant
// is invisible in behaviour until it silently disables the exact fallback.
TEST_CASE("Static filter bounds have their derived values", "[predicates][filtered_kernel]") {
    constexpr double eps = std::numeric_limits<double>::epsilon() / 2;
    STATIC_REQUIRE(orient2d_bound_a == (3.0 + 16.0 * eps) * eps);
    STATIC_REQUIRE(incircle_bound_a == (10.0 + 96.0 * eps) * eps);
    STATIC_REQUIRE(orient2d_bound_a > 0.0);
    STATIC_REQUIRE(incircle_bound_a > orient2d_bound_a);
}

// ---------------------------------------------------------------------------
// The headline invariant: the filter never changes an answer.
// ---------------------------------------------------------------------------

TEST_CASE("FilteredKernel::orient2d agrees with its exact backend", "[predicates][filtered_kernel]") {
    for (const auto& t : orientation_corpus()) {
        INFO(describe(t[0], t[1], t[2]));
        REQUIRE(Filtered::orient2d(t[0], t[1], t[2]) == RefExact::orient2d(t[0], t[1], t[2]));
    }
}

TEST_CASE("FilteredKernel::orient2d agrees with the exact backend near degeneracy", "[predicates][filtered_kernel][adversarial]") {
    std::mt19937_64 rng{424242};
    for (int i = 0; i < 500; ++i) {
        const auto t = terrain::test::near_degenerate_triple(rng);
        INFO(describe(t.a, t.b, t.c));
        REQUIRE(Filtered::orient2d(t.a, t.b, t.c) == t.expected);
        REQUIRE(Filtered::orient2d(t.a, t.b, t.c) == RefExact::orient2d(t.a, t.b, t.c));
    }
}

// ---------------------------------------------------------------------------
// Symmetries, which must hold on degenerate input too.
// ---------------------------------------------------------------------------

// Two callers asking about the same triple in different vertex orders must not
// disagree. This is unconditional: it has to hold on collinear input as well,
// which is exactly where an implementation that reorders its arguments to help
// the filter would break it.
TEST_CASE("FilteredKernel::orient2d is antisymmetric and cyclic, degeneracies included", "[predicates][filtered_kernel]") {
    for (const auto& t : orientation_corpus()) {
        INFO(describe(t[0], t[1], t[2]));
        const Orientation o = Filtered::orient2d(t[0], t[1], t[2]);

        REQUIRE(Filtered::orient2d(t[0], t[2], t[1]) == reversed(o));
        REQUIRE(Filtered::orient2d(t[1], t[0], t[2]) == reversed(o));
        REQUIRE(Filtered::orient2d(t[2], t[1], t[0]) == reversed(o));

        REQUIRE(Filtered::orient2d(t[1], t[2], t[0]) == o);
        REQUIRE(Filtered::orient2d(t[2], t[0], t[1]) == o);
    }
}

// ---------------------------------------------------------------------------
// Degeneracy is exact, not approximate.
// ---------------------------------------------------------------------------

TEST_CASE("Points on an exactly representable line are Collinear", "[predicates][filtered_kernel][adversarial]") {
    // A long collinear run on the integer lattice -- what a snap grid produces.
    const Point2 a{0.0, 0.0};
    const Point2 b{3.0, 7.0};
    for (int k = 2; k <= 64; ++k) {
        const Point2 c{3.0 * k, 7.0 * k};
        INFO(describe(a, b, c));
        REQUIRE(Filtered::orient2d(a, b, c) == Orientation::Collinear);
    }
}

TEST_CASE("Collinearity is exact at UTM33 magnitude", "[predicates][filtered_kernel][adversarial]") {
    std::mt19937_64 rng{99};
    for (int i = 0; i < 200; ++i) {
        const auto t = terrain::test::exactly_collinear_triple(rng);
        INFO(describe(t.a, t.b, t.c));
        REQUIRE(Filtered::orient2d(t.a, t.b, t.c) == Orientation::Collinear);
    }

    const Point2 a = utm33_offset({0, 0});
    const Point2 b = utm33_offset({1, 1});
    const Point2 c = utm33_offset({123456, 123456});
    REQUIRE(Filtered::orient2d(a, b, c) == Orientation::Collinear);
}

TEST_CASE("Four points on an exactly representable circle are Cocircular", "[predicates][filtered_kernel][adversarial]") {
    const auto circle = terrain::test::integer_circle_radius_5();
    for (std::size_t i = 0; i < circle.size(); ++i) {
        const Point2 a = circle[i];
        const Point2 b = circle[(i + 3) % circle.size()];
        const Point2 c = circle[(i + 6) % circle.size()];
        const Point2 d = circle[(i + 9) % circle.size()];
        INFO(describe(a, b, c, d));
        REQUIRE(Filtered::incircle(a, b, c, d) == Incircle::Cocircular);
    }
}

TEST_CASE("Cocircularity is exact at UTM33 magnitude", "[predicates][filtered_kernel][adversarial]") {
    const auto circle = terrain::test::integer_circle_radius_5();
    for (std::size_t i = 0; i < circle.size(); ++i) {
        const Point2 a = utm33_offset(circle[i]);
        const Point2 b = utm33_offset(circle[(i + 3) % circle.size()]);
        const Point2 c = utm33_offset(circle[(i + 6) % circle.size()]);
        const Point2 d = utm33_offset(circle[(i + 9) % circle.size()]);
        INFO(describe(a, b, c, d));
        REQUIRE(Filtered::incircle(a, b, c, d) == Incircle::Cocircular);
        REQUIRE(Filtered::incircle(a, b, c, utm33_offset({0, 0})) == Incircle::Inside);
        REQUIRE(Filtered::incircle(a, b, c, utm33_offset({40, 40})) == Incircle::Outside);
    }
}

// A sub-millimetre triangle embedded in a coordinate system spanning hundreds
// of kilometres: the everyday case for a TIN built from a snapped breakline
// network in UTM33, and the one naive determinants have no digits left for.
TEST_CASE("FilteredKernel resolves a sub-millimetre triangle at UTM33 magnitude", "[predicates][filtered_kernel][adversarial]") {
    // 2^-14 m is about 0.06 mm, and is exactly representable, so adding it to
    // an integer easting/northing loses nothing.
    const double tiny = 0.00006103515625;
    const Point2 a = utm33_offset({0, 0});
    const Point2 b{a.x + tiny, a.y};
    const Point2 c{a.x, a.y + tiny};

    REQUIRE(Filtered::orient2d(a, b, c) == Orientation::CounterClockwise);
    REQUIRE(Filtered::orient2d(a, c, b) == Orientation::Clockwise);
}

// ---------------------------------------------------------------------------
// Filter discipline: the fallback fires exactly when it should.
// ---------------------------------------------------------------------------

// Without these two tests, a FilteredKernel that ignored its filter and always
// called the backend, and one whose filter swallowed every case, would both
// pass every other test in this file -- the first slowly, the second wrongly.
TEST_CASE("Well-separated input never reaches the exact backend", "[predicates][filtered_kernel][filter]") {
    CountingExact<RefExact>::reset();

    REQUIRE(Counted::orient2d(Point2{0, 0}, Point2{1, 0}, Point2{0, 1}) ==
            Orientation::CounterClockwise);
    REQUIRE(Counted::orient2d(Point2{-3, 2}, Point2{4, 2}, Point2{1, -5}) ==
            Orientation::Clockwise);
    REQUIRE(CountingExact<RefExact>::orient2d_calls == 0);

    // Reset before the incircle block so its two counters are read against a
    // known zero rather than against whatever the orient2d block left behind.
    // `orient2d_calls` is the load-bearing one here: `incircle` normalizes by
    // calling `FilteredKernel::orient2d`, never `E::orient2d`, so on
    // well-separated input the filter answers the orientation too and the
    // backend is not touched by either predicate. Asserting only
    // `incircle_calls` would leave that silent -- and would make the test's
    // meaning depend on the fact that it happened to assert before the
    // incircle calls were made.
    CountingExact<RefExact>::reset();

    const Point2 a{5, 0};
    const Point2 b{0, 5};
    const Point2 c{-5, 0};
    REQUIRE(Counted::incircle(a, b, c, Point2{0, 0}) == Incircle::Inside);
    REQUIRE(Counted::incircle(a, b, c, Point2{400, 400}) == Incircle::Outside);
    REQUIRE(CountingExact<RefExact>::incircle_calls == 0);
    REQUIRE(CountingExact<RefExact>::orient2d_calls == 0);
}

TEST_CASE("A collinear triple reaches the exact backend", "[predicates][filtered_kernel][filter]") {
    CountingExact<RefExact>::reset();

    REQUIRE(Counted::orient2d(Point2{0, 0}, Point2{1, 1}, Point2{2, 2}) == Orientation::Collinear);
    REQUIRE(CountingExact<RefExact>::orient2d_calls > 0);
}

TEST_CASE("A cocircular quadruple reaches the exact backend", "[predicates][filtered_kernel][filter]") {
    CountingExact<RefExact>::reset();

    REQUIRE(Counted::incircle(Point2{5, 0}, Point2{0, 5}, Point2{-5, 0}, Point2{3, -4}) ==
            Incircle::Cocircular);
    REQUIRE(CountingExact<RefExact>::incircle_calls > 0);
}

// ---------------------------------------------------------------------------
// incircle is total: no counterclockwise precondition survives to the backend.
// ---------------------------------------------------------------------------

// The highest-value test in this increment, and the reason it is written as an
// ordinary unconditional TEST_CASE rather than guarded on NDEBUG: it has to run
// in the Debug build.
//
// detria's `incircle` asserts under `#ifndef NDEBUG` that its first three
// points are counterclockwise, and its assert handler raises SIGTRAP rather
// than throwing or returning. Our asan+ubsan CI job is a Debug build, so once
// DetriaExact is wired up, an un-normalized call from FilteredKernel::incircle
// would not fail a test -- it would kill the job on a signal, with no
// assertion text and nothing to reproduce from a Release build.
//
// FilteredKernel::incircle is therefore specified as total and
// precondition-free: it normalizes first, and `incircle_ccw` is only ever
// reached with a genuinely counterclockwise triple. CcwCheckingExact is a
// stand-in for detria's assertion that records the violation instead of
// raising, so the failure is legible today, before the real backend exists.
TEST_CASE("incircle never passes a non-counterclockwise triple to the backend", "[predicates][filtered_kernel][incircle]") {
    const Point2 a{5, 0};
    const Point2 b{0, 5};
    const Point2 c{-5, 0};
    const Point2 inside{0, 0};
    const Point2 outside{400, 400};
    const Point2 on{3, -4};

    CcwCheckingExact<RefExact>::reset();

    // Clockwise first three points.
    REQUIRE(CcwChecked::incircle(a, c, b, inside) == Incircle::Inside);
    REQUIRE(CcwChecked::incircle(a, c, b, outside) == Incircle::Outside);
    REQUIRE(CcwChecked::incircle(a, c, b, on) == Incircle::Cocircular);

    REQUIRE(CcwCheckingExact<RefExact>::violations == 0);

    // Collinear first three points: a degenerate circle, by definition, so
    // every fourth point is Cocircular. This is not merely permission to skip
    // the backend, it is the specified behaviour -- `incircle` returns
    // Cocircular from the orientation alone and never calls `incircle_ccw`,
    // because there is no counterclockwise triple it could legally pass. The
    // call count is the only way to observe that from outside, so it is
    // asserted rather than described; a reset here puts it on a known zero.
    CcwCheckingExact<RefExact>::reset();

    REQUIRE(CcwChecked::incircle(Point2{0, 0}, Point2{1, 1}, Point2{2, 2}, Point2{99, -7}) ==
            Incircle::Cocircular);
    REQUIRE(CcwChecked::incircle(Point2{0, 0}, Point2{2, 2}, Point2{1, 1}, Point2{0, 0}) ==
            Incircle::Cocircular);

    // Coincident points are collinear, hence cocircular, hence still not a
    // legal argument for the backend.
    REQUIRE(CcwChecked::incircle(Point2{7, 7}, Point2{7, 7}, Point2{7, 7}, Point2{1, 2}) ==
            Incircle::Cocircular);

    REQUIRE(CcwCheckingExact<RefExact>::incircle_calls == 0);
    REQUIRE(CcwCheckingExact<RefExact>::violations == 0);
}

// The consequence the CDT's local-Delaunay test and `flip` depend on: a
// triangle's three vertices may arrive in any rotation or reflection, and the
// answer about a fourth point may not depend on which.
TEST_CASE("incircle is insensitive to the order of the first three points", "[predicates][filtered_kernel][incircle]") {
    const auto circle = terrain::test::integer_circle_radius_5();
    const Point2 a = circle[0];
    const Point2 b = circle[4];
    const Point2 c = circle[8];

    const std::array<Point2, 3> permutations[] = {
        {a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a},
    };

    for (const Point2& d : {Point2{0, 0}, Point2{400, 400}, circle[2], Point2{1, 1}}) {
        const Incircle expected = Filtered::incircle(a, b, c, d);
        for (const auto& p : permutations) {
            INFO(describe(p[0], p[1], p[2], d));
            REQUIRE(Filtered::incircle(p[0], p[1], p[2], d) == expected);
        }
    }
}

// ---------------------------------------------------------------------------
// Purity.
// ---------------------------------------------------------------------------

// Static member functions over an empty type, so there is no instance state to
// share, no lazy initialisation to race on, and nothing for a caller to have to
// construct. Checked structurally rather than by comment: taking the address of
// a non-static member function yields a pointer-to-member, which is not a
// pointer to function, so this fails to compile if the signatures drift.
//
// `GeometryKernel` and `ExactPredicates` are spelled with qualified static
// calls -- `{ K::orient2d(a, b, c) }`, not `{ k.orient2d(a, b, c) }` -- so the
// static form is required of every model, not just of this one. This test
// stays because the concept does not pin emptiness: a model could be
// static-callable and still carry mutable static state, which is what would
// make the concurrency test below start failing.
TEST_CASE("FilteredKernel is stateless and its predicates are static", "[predicates][filtered_kernel][purity]") {
    STATIC_REQUIRE(std::is_empty_v<Filtered>);
    STATIC_REQUIRE(std::is_function_v<std::remove_pointer_t<decltype(&Filtered::orient2d)>>);
    STATIC_REQUIRE(std::is_function_v<std::remove_pointer_t<decltype(&Filtered::incircle)>>);
}

TEST_CASE("FilteredKernel is callable concurrently", "[predicates][filtered_kernel][purity]") {
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
                    if (Filtered::orient2d(triple.a, triple.b, triple.c) != triple.expected) {
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
