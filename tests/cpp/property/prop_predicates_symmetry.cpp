// Property tests for the symmetries of terrain::pred::FilteredKernel.
//
// The unit suites pin named configurations; this file asserts the same laws
// over generated input, which is where argument-order bugs that survive a
// hand-picked corpus tend to show up. Two callers asking about the same triple
// in different vertex orders must never disagree -- if they can, a CDT can
// enter an infinite flip loop with every individual predicate call looking
// defensible in isolation.
//
// No rapidcheck: Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64 gives reproducible generated input without standing up a
// second framework alongside a brand new module. A failure prints its seed via
// the generator index, and rerunning that section reproduces it exactly.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <exact_reference.hpp>
#include <point_families.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <array>
#include <cstdint>
#include <format>
#include <random>
#include <string>

using terrain::Point2;
using terrain::pred::FilteredKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::pred::reversed;
using terrain::test::CcwCheckingExact;
using terrain::test::RefExact;
using terrain::test::utm33_offset;

namespace {

using Filtered = FilteredKernel<RefExact>;
using CcwChecked = FilteredKernel<CcwCheckingExact<RefExact>>;

constexpr int seed_count = 32;

// A coarse integer lattice. The range is deliberately tiny: on a 9x9 grid a
// random triple is collinear or duplicated often enough that the degenerate
// cases are sampled heavily rather than as an afterthought.
[[nodiscard]] Point2 lattice_point(std::mt19937_64& rng) {
    std::uniform_int_distribution<int> coord{-4, 4};
    return Point2{static_cast<double>(coord(rng)), static_cast<double>(coord(rng))};
}

[[nodiscard]] Point2 utm_lattice_point(std::mt19937_64& rng) {
    return utm33_offset(lattice_point(rng));
}

std::string describe(const Point2& a, const Point2& b, const Point2& c) {
    return std::format("a={} b={} c={}", a, b, c);
}

std::string describe(const Point2& a, const Point2& b, const Point2& c, const Point2& d) {
    return std::format("a={} b={} c={} d={}", a, b, c, d);
}

}  // namespace

TEST_CASE("prop: orient2d is antisymmetric under a transposition", "[predicates][property][orientation]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed)};
    const bool at_utm_scale = seed % 2 == 0;

    for (int i = 0; i < 200; ++i) {
        const Point2 a = at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng);
        const Point2 b = at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng);
        const Point2 c = at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng);

        INFO(describe(a, b, c));
        const Orientation o = Filtered::orient2d(a, b, c);
        REQUIRE(Filtered::orient2d(a, c, b) == reversed(o));
        REQUIRE(Filtered::orient2d(b, a, c) == reversed(o));
        REQUIRE(Filtered::orient2d(c, b, a) == reversed(o));
    }
}

TEST_CASE("prop: orient2d is invariant under cyclic rotation", "[predicates][property][orientation]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 1000};
    const bool at_utm_scale = seed % 2 == 0;

    for (int i = 0; i < 200; ++i) {
        const Point2 a = at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng);
        const Point2 b = at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng);
        const Point2 c = at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng);

        INFO(describe(a, b, c));
        const Orientation o = Filtered::orient2d(a, b, c);
        REQUIRE(Filtered::orient2d(b, c, a) == o);
        REQUIRE(Filtered::orient2d(c, a, b) == o);
    }
}

TEST_CASE("prop: the filter never changes an orient2d answer", "[predicates][property][orientation]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 2000};

    for (int i = 0; i < 100; ++i) {
        // Three scales in one run: unit lattice, UTM33 lattice, and the
        // near-degenerate family where the naive determinant has no correct
        // digits at all.
        const int family = i % 3;
        Point2 a{};
        Point2 b{};
        Point2 c{};
        if (family == 0) {
            a = lattice_point(rng);
            b = lattice_point(rng);
            c = lattice_point(rng);
        } else if (family == 1) {
            a = utm_lattice_point(rng);
            b = utm_lattice_point(rng);
            c = utm_lattice_point(rng);
        } else {
            const auto t = terrain::test::near_degenerate_triple(rng);
            a = t.a;
            b = t.b;
            c = t.c;
        }

        INFO(describe(a, b, c));
        REQUIRE(Filtered::orient2d(a, b, c) == RefExact::orient2d(a, b, c));
    }
}

TEST_CASE("prop: incircle does not depend on the order of the first three points", "[predicates][property][incircle]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 3000};
    const bool at_utm_scale = seed % 2 == 0;

    for (int i = 0; i < 100; ++i) {
        const auto pick = [&] { return at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng); };
        const Point2 a = pick();
        const Point2 b = pick();
        const Point2 c = pick();
        const Point2 d = pick();

        INFO(describe(a, b, c, d));
        const std::array<Point2, 3> permutations[] = {
            {a, b, c}, {a, c, b}, {b, a, c}, {b, c, a}, {c, a, b}, {c, b, a},
        };
        const Incircle expected = Filtered::incircle(a, b, c, d);
        for (const auto& p : permutations) {
            REQUIRE(Filtered::incircle(p[0], p[1], p[2], d) == expected);
        }
    }
}

// The generated counterpart of the named unit test: over a lattice this coarse,
// collinear and clockwise triples are common, so this exercises the
// normalization on a large sample. It must keep running in Debug builds -- see
// the comment on CcwCheckingExact, which stands in for detria's SIGTRAP-raising
// counterclockwise assertion.
TEST_CASE("prop: incircle never hands a non-counterclockwise triple to the backend", "[predicates][property][incircle]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 4000};

    CcwCheckingExact<RefExact>::reset();

    for (int i = 0; i < 100; ++i) {
        const Point2 a = lattice_point(rng);
        const Point2 b = lattice_point(rng);
        const Point2 c = lattice_point(rng);
        const Point2 d = lattice_point(rng);

        INFO(describe(a, b, c, d));
        const int calls_before = CcwCheckingExact<RefExact>::incircle_calls;
        const Incircle result = CcwChecked::incircle(a, b, c, d);
        if (Filtered::orient2d(a, b, c) == Orientation::Collinear) {
            // The generated counterpart of the short-circuit assertion in the
            // unit suite: a degenerate circle is answered from the orientation
            // alone, so the backend is not consulted at all.
            REQUIRE(result == Incircle::Cocircular);
            REQUIRE(CcwCheckingExact<RefExact>::incircle_calls == calls_before);
        }
        REQUIRE(CcwCheckingExact<RefExact>::violations == 0);
    }
}
