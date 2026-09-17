// Property tests for the structural laws of terrain/core/ring.hpp.
//
// The unit suites pin named configurations; this file asserts the laws that
// must hold for *every* ring, over generated input. These are the laws a
// downstream CDT will lean on without ever restating them: that a ring's
// classification does not depend on which vertex the caller happened to list
// first, that reversing a ring negates its orientation and nothing else, and
// that a point outside the bounding box is outside the ring.
//
// Every invariant runs as a TEMPLATE_LIST_TEST_CASE over
// {FastKernel, DefaultKernel} x {PointRing, IndexedRing}. The IndexedRing
// instantiations are the point of the matrix: they are what proves increment
// 3's zero-copy path compiles and behaves against these algorithms, and they
// are what nothing else in the tree would instantiate on day one.
//
// No rapidcheck. Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64 gives reproducible input without a second framework; a
// failure prints its seed as the generator index and rerunning that section
// reproduces it exactly.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <point_families.hpp>
#include <ring_cases.hpp>

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/raster/geometry.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <random>
#include <span>
#include <string>
#include <vector>

using Catch::Matchers::WithinRel;
using terrain::Box2;
using terrain::Point2;
using terrain::PointInRing;
using terrain::PointRing;
using terrain::bounding_box;
using terrain::orientation;
using terrain::point_in_ring;
using terrain::signed_area;
using terrain::pred::DefaultKernel;
using terrain::pred::FastKernel;
using terrain::pred::Orientation;
using terrain::raster::RasterGeometry;
using terrain::test::RingCases;
using terrain::test::describe;
using terrain::test::flipped;
using terrain::test::integer_star_ring;
using terrain::test::rotated;
using terrain::test::star_ring;
using terrain::test::utm33_offset;

namespace {

constexpr int seed_count = 24;

[[nodiscard]] std::mt19937_64 seeded(int seed) {
    return std::mt19937_64{0x5eed0000ULL + static_cast<std::uint64_t>(seed)};
}

[[nodiscard]] std::span<const Point2> as_span(const std::vector<Point2>& v) {
    return std::span<const Point2>{v};
}

// Probe points for a ring generated about the origin with radii in [1, 3]:
// deep interior, the annulus where the boundary actually is, and far outside.
[[nodiscard]] std::vector<Point2> probes(std::mt19937_64& rng, std::size_t n) {
    std::uniform_real_distribution<double> coord{-5.0, 5.0};
    std::vector<Point2> out;
    out.reserve(n + 1);
    out.push_back(Point2{0.0, 0.0});  // the star centre: known interior
    for (std::size_t i = 0; i < n; ++i) out.push_back(Point2{coord(rng), coord(rng)});
    return out;
}

[[nodiscard]] std::size_t ring_size(std::mt19937_64& rng) {
    std::uniform_int_distribution<std::size_t> n{3, 12};
    return n(rng);
}

}  // namespace

// ---------------------------------------------------------------------------
// point_in_ring does not depend on the encoding
// ---------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("point_in_ring is invariant under cyclic rotation",
                        "[ring][property][point_in_ring]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::size_t n = ring_size(rng);
    const std::vector<Point2> base = star_ring(rng, n, Point2{0.0, 0.0}, 1.0, 3.0);
    const std::vector<Point2> qs = probes(rng, 12);

    const typename TestType::Model::Holder h0{as_span(base)};
    for (std::size_t k = 1; k < n; ++k) {
        const std::vector<Point2> turned = rotated(as_span(base), k);
        const typename TestType::Model::Holder hk{as_span(turned)};
        for (const Point2& p : qs) {
            const auto a = point_in_ring<K>(h0.ring(), p);
            const auto b = point_in_ring<K>(hk.ring(), p);
            INFO(std::format("seed={} n={} k={} p={} base={} rotated={}", seed, n, k, p,
                             describe(a), describe(b)));
            REQUIRE(a == b);
        }
    }
}

// Reversing a ring reverses every edge. The parity rule must not care: the
// even-odd classification is a property of the closed curve, not of the
// direction it is traced in. A rule that leaned on edge direction rather than
// on the half-open y test would fail exactly here.
TEMPLATE_LIST_TEST_CASE("point_in_ring is invariant under reversal",
                        "[ring][property][point_in_ring]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::vector<Point2> base = star_ring(rng, ring_size(rng), Point2{0.0, 0.0}, 1.0, 3.0);
    const std::vector<Point2> back = flipped(as_span(base));
    const std::vector<Point2> qs = probes(rng, 12);

    const typename TestType::Model::Holder forward{as_span(base)};
    const typename TestType::Model::Holder reverse{as_span(back)};
    for (const Point2& p : qs) {
        const auto a = point_in_ring<K>(forward.ring(), p);
        const auto b = point_in_ring<K>(reverse.ring(), p);
        INFO(std::format("seed={} p={} forward={} reversed={}", seed, p, describe(a), describe(b)));
        REQUIRE(a == b);
    }
}

TEMPLATE_LIST_TEST_CASE("every vertex of every ring is Boundary",
                        "[ring][property][point_in_ring][boundary]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::vector<Point2> base = star_ring(rng, ring_size(rng), Point2{0.0, 0.0}, 1.0, 3.0);
    const typename TestType::Model::Holder h{as_span(base)};
    const auto r = h.ring();

    for (std::size_t i = 0; i < r.size(); ++i) {
        INFO(std::format("seed={} i={} v={}", seed, i, r.vertex(i)));
        REQUIRE(point_in_ring<K>(r, r.vertex(i)) == PointInRing::Boundary);
    }
}

// The centre of a star-shaped ring is interior by construction, so this is an
// oracle rather than a self-consistency check.
TEMPLATE_LIST_TEST_CASE("the centre of a star-shaped ring is Inside",
                        "[ring][property][point_in_ring]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const Point2 c{0.0, 0.0};
    const std::vector<Point2> base = star_ring(rng, ring_size(rng), c, 1.0, 3.0);
    const typename TestType::Model::Holder h{as_span(base)};

    INFO(std::format("seed={}", seed));
    REQUIRE(point_in_ring<K>(h.ring(), c) == PointInRing::Inside);
}

// ---------------------------------------------------------------------------
// bounding_box bounds the classification
// ---------------------------------------------------------------------------

TEMPLATE_LIST_TEST_CASE("nothing outside the bounding box is Inside or Boundary",
                        "[ring][property][bbox]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::vector<Point2> base = star_ring(rng, ring_size(rng), Point2{0.0, 0.0}, 1.0, 3.0);
    const typename TestType::Model::Holder h{as_span(base)};
    const auto r = h.ring();
    const Box2 b = bounding_box(r);

    for (std::size_t i = 0; i < r.size(); ++i) {
        REQUIRE(b.contains(r.vertex(i)));
    }
    for (const Point2& p : probes(rng, 32)) {
        if (b.contains(p)) continue;
        INFO(std::format("seed={} p={} box={}", seed, p, b));
        REQUIRE(point_in_ring<K>(r, p) == PointInRing::Outside);
    }
}

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

// signed_area's sign is not an oracle in general -- that is the whole reason
// orientation exists and does not sum. On rings whose area is well separated
// from zero, though, the two must agree, and a disagreement here would mean
// one of them is simply broken rather than merely imprecise.
TEMPLATE_LIST_TEST_CASE("orientation agrees with sign(signed_area) away from degeneracy",
                        "[ring][property][orientation][area]", RingCases) {
    using K = typename TestType::Kernel;
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::vector<Point2> base = star_ring(rng, ring_size(rng), Point2{0.0, 0.0}, 1.0, 3.0);
    const std::vector<Point2> back = flipped(as_span(base));
    const typename TestType::Model::Holder forward{as_span(base)};
    const typename TestType::Model::Holder reverse{as_span(back)};

    const double area = signed_area(forward.ring());
    INFO(std::format("seed={} area={}", seed, area));
    REQUIRE(std::abs(area) > 1e-3);  // the generator is supposed to stay well clear of zero

    REQUIRE(orientation<K>(forward.ring()) == Orientation::CounterClockwise);
    REQUIRE(area > 0.0);
    REQUIRE(orientation<K>(reverse.ring()) == Orientation::Clockwise);
    REQUIRE(signed_area(reverse.ring()) < 0.0);
}

// ---------------------------------------------------------------------------
// signed_area's magnitude invariances
// ---------------------------------------------------------------------------

// On an integer lattice every product and every partial sum in the translated
// shoelace is exact, so the invariance is bit for bit. Asserting it exactly is
// what distinguishes a genuine ordering bug from ordinary drift.
TEMPLATE_LIST_TEST_CASE("|signed_area| is bit-exact under rotation and reversal on a lattice",
                        "[ring][property][area]", RingCases) {
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::size_t n = ring_size(rng);
    const std::vector<Point2> base = integer_star_ring(rng, n);
    const typename TestType::Model::Holder h0{as_span(base)};
    const double expected = std::abs(signed_area(h0.ring()));

    INFO(std::format("seed={} n={} area={}", seed, n, expected));
    REQUIRE(expected > 0.0);

    for (std::size_t k = 1; k < n; ++k) {
        const std::vector<Point2> turned = rotated(as_span(base), k);
        const typename TestType::Model::Holder hk{as_span(turned)};
        REQUIRE(std::abs(signed_area(hk.ring())) == expected);
    }

    const std::vector<Point2> back = flipped(as_span(base));
    const typename TestType::Model::Holder hr{as_span(back)};
    REQUIRE(std::abs(signed_area(hr.ring())) == expected);
}

// At UTM33 magnitudes the invariance is no longer exact, and the tolerance
// here is a statement about double arithmetic, not about the header: shifting
// the translation origin from one vertex to another changes the rounding.
// Translating to vertex(0) is what keeps the error at this scale instead of
// the 1e-6 absolute a non-translated shoelace would produce.
TEMPLATE_LIST_TEST_CASE("|signed_area| survives rotation at UTM33 magnitudes",
                        "[ring][property][area][utm33]", RingCases) {
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::size_t n = ring_size(rng);
    const Point2 c = utm33_offset(Point2{0.0, 0.0});
    const std::vector<Point2> base = star_ring(rng, n, c, 100.0, 1000.0);
    const typename TestType::Model::Holder h0{as_span(base)};
    const double expected = std::abs(signed_area(h0.ring()));

    INFO(std::format("seed={} n={} area={}", seed, n, expected));
    for (std::size_t k = 1; k < n; ++k) {
        const std::vector<Point2> turned = rotated(as_span(base), k);
        const typename TestType::Model::Holder hk{as_span(turned)};
        REQUIRE_THAT(std::abs(signed_area(hk.ring())), WithinRel(expected, 1e-9));
    }
}

// ---------------------------------------------------------------------------
// The two kernels
// ---------------------------------------------------------------------------

// On well-separated input the filter never reaches its exact backend, so this
// is also a check that FilteredKernel's fast path answers the same question
// FastKernel does rather than a subtly different one.
TEST_CASE("FastKernel and DefaultKernel agree on well-separated input",
          "[ring][property][point_in_ring][kernels]") {
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    const std::vector<Point2> base = star_ring(rng, ring_size(rng), Point2{0.0, 0.0}, 1.0, 3.0);
    const PointRing r{as_span(base)};

    for (const Point2& p : probes(rng, 48)) {
        const auto fast = point_in_ring<FastKernel>(r, p);
        const auto exact = point_in_ring<DefaultKernel>(r, p);
        INFO(std::format("seed={} p={} fast={} exact={}", seed, p, describe(fast),
                         describe(exact)));
        REQUIRE(fast == exact);
    }
}

// ---------------------------------------------------------------------------
// The raster boundary
// ---------------------------------------------------------------------------

// The first place these two modules meet. RasterGeometry::contains_strict is
// the epsilon-tolerant "well inside the DEM" test that boundary draping uses;
// point_in_ring on the ring through the four corner nodes is the exact one.
// The tolerant test must be the stricter of the two, or a point the draping
// code believes is safely interior can be outside the polygon it is meshing.
TEST_CASE("contains_strict implies Inside on the raster's corner ring",
          "[ring][property][raster]") {
    const RasterGeometry geometry{utm33_offset(Point2{0.0, 0.0}).x,
                                  utm33_offset(Point2{0.0, 0.0}).y,
                                  10.0, 10.0, 64, 48};
    const std::vector<Point2> corners{
        geometry.node({0, 0}),
        geometry.node({0, geometry.cols() - 1}),
        geometry.node({geometry.rows() - 1, geometry.cols() - 1}),
        geometry.node({geometry.rows() - 1, 0}),
    };
    const PointRing r{as_span(corners)};

    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);
    std::uniform_real_distribution<double> fx{-0.1, 1.1};

    const double w = geometry.x_max() - geometry.x_min();
    const double h = geometry.y_max() - geometry.y_min();

    int interior = 0;
    for (int i = 0; i < 64; ++i) {
        const Point2 p{geometry.x_min() + fx(rng) * w, geometry.y_min() + fx(rng) * h};
        if (!geometry.contains_strict(p)) continue;
        ++interior;
        INFO(std::format("seed={} p={}", seed, p));
        REQUIRE(point_in_ring<DefaultKernel>(r, p) == PointInRing::Inside);
    }
    REQUIRE(interior > 0);  // the sampler must actually be hitting the interior
}

// The corners themselves lie on the ring, and contains_strict rejects them.
TEST_CASE("the raster's corner nodes are on the corner ring and not strictly inside",
          "[ring][raster][boundary]") {
    const RasterGeometry geometry{utm33_offset(Point2{0.0, 0.0}).x,
                                  utm33_offset(Point2{0.0, 0.0}).y,
                                  10.0, 10.0, 64, 48};
    const std::vector<Point2> corners{
        geometry.node({0, 0}),
        geometry.node({0, geometry.cols() - 1}),
        geometry.node({geometry.rows() - 1, geometry.cols() - 1}),
        geometry.node({geometry.rows() - 1, 0}),
    };
    const PointRing r{as_span(corners)};

    for (const Point2& c : corners) {
        REQUIRE(point_in_ring<DefaultKernel>(r, c) == PointInRing::Boundary);
        REQUIRE_FALSE(geometry.contains_strict(c));
    }
}
