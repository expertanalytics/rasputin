// Property tests for the structural laws of orientation over a ring.
//
// The unit suites pin named configurations; this file asserts the laws that
// must hold for *every* ring, over generated input: that a ring's orientation
// does not depend on which vertex the caller happened to list first, and that
// reversing a ring negates its orientation and nothing else.
//
// Every invariant runs as a TEMPLATE_LIST_TEST_CASE over
// {FastKernel, DefaultKernel} x IndexedRing. The IndexedRing holder stores the
// vertices reversed behind a decoy (ring_cases.hpp), so these cases also prove
// the zero-copy path reads through the chain.
//
// No rapidcheck. Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64 gives reproducible input without a second framework; a
// failure prints its seed as the generator index and rerunning that section
// reproduces it exactly.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <ring_cases.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include <cstddef>
#include <cstdint>
#include <format>
#include <random>
#include <span>
#include <vector>

using terrain::Point2;
using terrain::orientation;
using terrain::pred::Orientation;
using terrain::test::RingCases;
using terrain::test::flipped;
using terrain::test::rotated;
using terrain::test::star_ring;

namespace {

constexpr int seed_count = 24;

[[nodiscard]] std::mt19937_64 seeded(int seed) {
    return std::mt19937_64{0x5eed0000ULL + static_cast<std::uint64_t>(seed)};
}

[[nodiscard]] std::span<const Point2> as_span(const std::vector<Point2>& v) {
    return std::span<const Point2>{v};
}

[[nodiscard]] std::size_t ring_size(std::mt19937_64& rng) {
    std::uniform_int_distribution<std::size_t> n{3, 12};
    return n(rng);
}

}  // namespace

// ---------------------------------------------------------------------------
// orientation
// ---------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("orientation is invariant under cyclic rotation",
                        "[ring][property][orientation]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::size_t n = ring_size(rng);
    const std::vector<Point2> base = star_ring(rng, n, Point2{0.0, 0.0}, 1.0, 3.0);
    const typename TestType::Model::Holder h0{as_span(base)};
    const Orientation expected = orientation<K>(h0.ring());

    REQUIRE(expected == Orientation::CounterClockwise);  // the generator emits CCW
    for (std::size_t k = 1; k < n; ++k) {
        const std::vector<Point2> turned = rotated(as_span(base), k);
        const typename TestType::Model::Holder hk{as_span(turned)};
        INFO(std::format("seed={} n={} k={}", seed, n, k));
        REQUIRE(orientation<K>(hk.ring()) == expected);
    }
}

TEMPLATE_LIST_TEST_CASE("orientation negates exactly under reversal",
                        "[ring][property][orientation]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::vector<Point2> base = star_ring(rng, ring_size(rng), Point2{0.0, 0.0}, 1.0, 3.0);
    const std::vector<Point2> back = flipped(as_span(base));
    const typename TestType::Model::Holder forward{as_span(base)};
    const typename TestType::Model::Holder reverse{as_span(back)};

    INFO(std::format("seed={}", seed));
    REQUIRE(orientation<K>(reverse.ring()) ==
            terrain::pred::reversed(orientation<K>(forward.ring())));
}

// The star generator builds counterclockwise rings by construction
// (ring_cases.hpp), and that construction is the oracle: orientation must
// report it, and report its reversal as clockwise.
TEMPLATE_LIST_TEST_CASE("orientation reports a generated star ring's construction",
                        "[ring][property][orientation]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::vector<Point2> base = star_ring(rng, ring_size(rng), Point2{0.0, 0.0}, 1.0, 3.0);
    const std::vector<Point2> back = flipped(as_span(base));
    const typename TestType::Model::Holder forward{as_span(base)};
    const typename TestType::Model::Holder reverse{as_span(back)};

    INFO(std::format("seed={}", seed));
    REQUIRE(orientation<K>(forward.ring()) == Orientation::CounterClockwise);
    REQUIRE(orientation<K>(reverse.ring()) == Orientation::Clockwise);
}
