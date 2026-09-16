// Property tests for the detria-backed exact backend.
//
// The unit suite pins named configurations and the enum mapping. This file
// asserts the one law that subsumes all of them over generated input:
// `DetriaExact` and `RefExact` are two independently written exact
// implementations of the same two determinants, so on every input they admit in
// common they must return the same classification. A disagreement is a bug in
// one of them, and which one is a separate question the named tests answer.
//
// The generators stay inside RefExact's domain -- integer-valued coordinates
// within its 2^52 (orient2d) and 2^25 (incircle) bounds -- because outside it
// the oracle throws rather than referees. That is a real limit on this file's
// reach and not a flaw in the oracle: an oracle that guessed would be worse
// than none.
//
// No rapidcheck: Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64, as in prop_predicates_symmetry.cpp. A failure prints its
// seed via the generator index, and rerunning that section reproduces it
// exactly.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ccw_normalize.hpp>
#include <exact_reference.hpp>
#include <point_families.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/detria_exact.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <cstdint>
#include <format>
#include <random>
#include <string>

using terrain::Point2;
using terrain::pred::DefaultKernel;
using terrain::pred::DetriaExact;
using terrain::pred::FilteredKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::test::as_ccw;
using terrain::test::RefExact;
using terrain::test::utm33_offset;

namespace {

using Reference = FilteredKernel<RefExact>;

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

// ---------------------------------------------------------------------------
// Backend against oracle.
// ---------------------------------------------------------------------------

TEST_CASE("prop: DetriaExact::orient2d agrees with RefExact", "[predicates][property][detria]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 5000};

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
        REQUIRE(DetriaExact::orient2d(a, b, c) == RefExact::orient2d(a, b, c));
    }
}

// Every triple is normalized through `as_ccw` -- which consults the oracle, not
// the subject -- before it reaches `incircle_ccw`. Calling the backend directly
// means the precondition is the caller's to establish, and with detria behind
// it the penalty for not doing so is SIGTRAP rather than a wrong answer.
// Collinear triples are skipped because no counterclockwise ordering of them
// exists; `FilteredKernel::incircle` reaches the same conclusion and answers
// Cocircular without consulting the backend at all.
TEST_CASE("prop: DetriaExact::incircle_ccw agrees with RefExact", "[predicates][property][detria]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 6000};
    const bool at_utm_scale = seed % 2 == 0;

    int normalized = 0;
    for (int i = 0; i < 100; ++i) {
        const auto pick = [&] { return at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng); };
        const Point2 a = pick();
        const Point2 b = pick();
        const Point2 c = pick();
        const Point2 d = pick();

        const auto triple = as_ccw(a, b, c);
        if (!triple.has_value()) {
            continue;
        }
        ++normalized;

        INFO(describe(triple->a, triple->b, triple->c, d));
        REQUIRE(DetriaExact::incircle_ccw(triple->a, triple->b, triple->c, d) ==
                RefExact::incircle_ccw(triple->a, triple->b, triple->c, d));
    }

    // On a lattice this coarse a large fraction of triples are degenerate, so
    // without this the test could quietly assert nothing.
    REQUIRE(normalized > 0);
}

// ---------------------------------------------------------------------------
// Kernel against kernel: swapping the backend moves no answer.
// ---------------------------------------------------------------------------

TEST_CASE("prop: DefaultKernel::orient2d agrees with FilteredKernel<RefExact>", "[predicates][property][default_kernel]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 7000};
    const bool at_utm_scale = seed % 2 == 0;

    for (int i = 0; i < 200; ++i) {
        const auto pick = [&] { return at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng); };
        const Point2 a = pick();
        const Point2 b = pick();
        const Point2 c = pick();

        INFO(describe(a, b, c));
        REQUIRE(DefaultKernel::orient2d(a, b, c) == Reference::orient2d(a, b, c));
    }
}

// Unlike the `incircle_ccw` property above, nothing is normalized here and
// nothing is skipped: `FilteredKernel::incircle` is specified as total, so
// clockwise, collinear and coincident triples are all in scope and are common
// on a 9x9 lattice. That totality is exactly what keeps the call from reaching
// detria's counterclockwise assertion.
//
// Be aware of the failure mode: if the normalization regresses, this test does
// not fail, it takes the process down on SIGTRAP mid-run. The forked probe in
// test_predicates_default_kernel.cpp is the place that turns that into a
// readable report; keep it in mind before reading a silent ctest timeout here
// as an unrelated infrastructure problem.
TEST_CASE("prop: DefaultKernel::incircle agrees with FilteredKernel<RefExact>, degeneracies included", "[predicates][property][default_kernel][incircle][signal]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 8000};
    const bool at_utm_scale = seed % 2 == 0;

    int degenerate_seen = 0;
    for (int i = 0; i < 100; ++i) {
        const auto pick = [&] { return at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng); };
        const Point2 a = pick();
        const Point2 b = pick();
        const Point2 c = pick();
        const Point2 d = pick();

        if (RefExact::orient2d(a, b, c) != Orientation::CounterClockwise) {
            ++degenerate_seen;
        }

        INFO(describe(a, b, c, d));
        REQUIRE(DefaultKernel::incircle(a, b, c, d) == Reference::incircle(a, b, c, d));
    }

    // The clockwise and collinear triples are the whole point of this test, so
    // a run that happened to contain none would be reporting a pass it did not
    // earn.
    REQUIRE(degenerate_seen > 0);
}

// The law the CDT's local-Delaunay test and `flip` depend on, asserted against
// the real backend: a triangle's vertices may arrive in any rotation or
// reflection, and the answer about a fourth point may not depend on which.
// Under the default kernel this exercises both branches of the normalization on
// every sample.
TEST_CASE("prop: DefaultKernel::incircle does not depend on the order of the first three points", "[predicates][property][default_kernel][incircle]") {
    const int seed = GENERATE(range(0, seed_count));
    std::mt19937_64 rng{static_cast<std::uint64_t>(seed) + 9000};
    const bool at_utm_scale = seed % 2 == 0;

    for (int i = 0; i < 100; ++i) {
        const auto pick = [&] { return at_utm_scale ? utm_lattice_point(rng) : lattice_point(rng); };
        const Point2 a = pick();
        const Point2 b = pick();
        const Point2 c = pick();
        const Point2 d = pick();

        INFO(describe(a, b, c, d));
        const Incircle expected = DefaultKernel::incircle(a, b, c, d);
        REQUIRE(DefaultKernel::incircle(a, c, b, d) == expected);
        REQUIRE(DefaultKernel::incircle(b, a, c, d) == expected);
        REQUIRE(DefaultKernel::incircle(b, c, a, d) == expected);
        REQUIRE(DefaultKernel::incircle(c, a, b, d) == expected);
        REQUIRE(DefaultKernel::incircle(c, b, a, d) == expected);
    }
}
