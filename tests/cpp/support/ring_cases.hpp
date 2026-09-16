#pragma once

// Test-only scaffolding for the ring suites: the (kernel x ring model) matrix
// the structural invariants run over, and deterministic generators for rings
// whose classification is known by construction.
//
// The matrix matters. PointRing is the contiguous model and is what every
// hand-written test reaches for, so it is the one that gets exercised by
// accident; IndexedRing is the zero-copy view the PSLG will hand the
// algorithms, and nothing would instantiate it on day one unless a test
// insisted. The holders below store an IndexedRing's points *reversed behind a
// decoy vertex*, so an algorithm that quietly reads the vertex buffer in
// storage order instead of going through vertex(i) fails rather than passes.
//
// Every generator takes an explicit std::mt19937_64, matching point_families:
// a test seeds once and the whole sequence is reproducible.

#include <terrain/core/point.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <span>
#include <string>
#include <tuple>
#include <vector>

namespace terrain::test {

// ---------------------------------------------------------------------------
// Ring models, behind a uniform "build one from a vertex list" interface
// ---------------------------------------------------------------------------

struct PointRingCase {
    static constexpr const char* name = "PointRing";

    class Holder {
    public:
        explicit Holder(std::span<const Point2> pts) : pts_{pts.begin(), pts.end()} {}

        [[nodiscard]] PointRing ring() const { return PointRing{std::span<const Point2>{pts_}}; }

    private:
        std::vector<Point2> pts_;
    };
};

struct IndexedRingCase {
    static constexpr const char* name = "IndexedRing";

    class Holder {
    public:
        explicit Holder(std::span<const Point2> pts) {
            // A decoy at slot 0 and the vertices stored backwards: the identity
            // mapping is wrong for every index, so an implementation that
            // ignores chain_indices cannot accidentally agree.
            storage_.push_back(Point2{987654.0, -987654.0});
            storage_.insert(storage_.end(), pts.rbegin(), pts.rend());
            const auto n = static_cast<std::uint32_t>(pts.size());
            for (std::uint32_t i = 0; i < n; ++i) index_.push_back(n - i);
        }

        [[nodiscard]] IndexedRing ring() const {
            return IndexedRing{std::span<const Point2>{storage_},
                               std::span<const std::uint32_t>{index_}};
        }

    private:
        std::vector<Point2> storage_;
        std::vector<std::uint32_t> index_;
    };
};

template <typename K, typename RC>
struct RingCase {
    using Kernel = K;
    using Model = RC;
};

using FastPointCase = RingCase<pred::FastKernel, PointRingCase>;
using FastIndexedCase = RingCase<pred::FastKernel, IndexedRingCase>;
using DefaultPointCase = RingCase<pred::DefaultKernel, PointRingCase>;
using DefaultIndexedCase = RingCase<pred::DefaultKernel, IndexedRingCase>;

// The four instantiations every structural invariant is required to hold for.
using RingCases =
    std::tuple<FastPointCase, FastIndexedCase, DefaultPointCase, DefaultIndexedCase>;

// The subset that is exact. Anything asserting a *correct* answer on
// near-degenerate input belongs here, not in RingCases.
using ExactRingCases = std::tuple<DefaultPointCase, DefaultIndexedCase>;

// ---------------------------------------------------------------------------
// Vertex-list transforms
// ---------------------------------------------------------------------------

[[nodiscard]] inline std::vector<Point2> rotated(std::span<const Point2> pts, std::size_t k) {
    std::vector<Point2> out;
    out.reserve(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i) out.push_back(pts[(i + k) % pts.size()]);
    return out;
}

[[nodiscard]] inline std::vector<Point2> flipped(std::span<const Point2> pts) {
    return std::vector<Point2>{pts.rbegin(), pts.rend()};
}

[[nodiscard]] inline std::vector<Point2> translated(std::span<const Point2> pts,
                                                    const Point2& offset) {
    std::vector<Point2> out;
    out.reserve(pts.size());
    for (const Point2& p : pts) out.push_back(p + offset);
    return out;
}

// ---------------------------------------------------------------------------
// Rings whose classification is known by construction
// ---------------------------------------------------------------------------

// A star-shaped ring about `center`: n vertices at strictly increasing angles,
// radii in [r_min, r_max]. Star-shaped polygons are simple by construction --
// no rejection loop, no chance of a generated self-intersection sneaking a
// wrong expectation into a property test -- and `center` is strictly interior
// by construction too, which gives every generated ring one point whose
// classification is known without an oracle. Vertices come out
// counterclockwise.
[[nodiscard]] inline std::vector<Point2> star_ring(std::mt19937_64& rng, std::size_t n,
                                                   const Point2& center, double r_min,
                                                   double r_max) {
    std::uniform_real_distribution<double> jitter{0.1, 0.9};
    std::uniform_real_distribution<double> radius{r_min, r_max};

    const double step = 2.0 * std::numbers::pi / static_cast<double>(n);
    std::vector<Point2> pts;
    pts.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double theta = step * (static_cast<double>(i) + jitter(rng));
        const double r = radius(rng);
        pts.push_back(Point2{center.x + r * std::cos(theta), center.y + r * std::sin(theta)});
    }
    return pts;
}

// The same construction rounded onto an integer lattice. Radii are kept large
// enough relative to the rounding error that the angular order, and therefore
// the star-shapedness, survives. Integer coordinates make every step of the
// shoelace sum exact, so magnitude invariances can be asserted bit for bit
// rather than within a tolerance -- which is the only way to tell an ordering
// bug apart from ordinary floating-point drift.
[[nodiscard]] inline std::vector<Point2> integer_star_ring(std::mt19937_64& rng, std::size_t n) {
    std::vector<Point2> pts = star_ring(rng, n, Point2{0.0, 0.0}, 50.0, 200.0);
    for (Point2& p : pts) {
        p.x = std::round(p.x);
        p.y = std::round(p.y);
    }
    return pts;
}

// The axis-aligned unit square, counterclockwise, with no stored closure.
[[nodiscard]] inline std::vector<Point2> unit_square() {
    return {Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0}, Point2{0.0, 1.0}};
}

// A ring with a V-shaped notch cut into its top edge. The notch vertex sits at
// y = 2 with both neighbours above it, so a horizontal ray at y = 2 touches a
// local minimum -- the configuration the half-open-in-y rule exists for.
[[nodiscard]] inline std::vector<Point2> notched_ring() {
    return {
        Point2{0.0, 0.0}, Point2{6.0, 0.0}, Point2{6.0, 6.0},
        Point2{3.0, 2.0},  // the notch
        Point2{0.0, 6.0},
    };
}

// A bowtie: the diagonals cross at (2, 2). Self-intersecting rings are
// accepted and not detected in this increment, and the even-odd rule gives
// them a total, deterministic classification.
[[nodiscard]] inline std::vector<Point2> bowtie_ring() {
    return {Point2{0.0, 0.0}, Point2{4.0, 4.0}, Point2{4.0, 0.0}, Point2{0.0, 4.0}};
}

[[nodiscard]] inline std::string describe(PointInRing c) {
    switch (c) {
        case PointInRing::Outside: return "Outside";
        case PointInRing::Boundary: return "Boundary";
        case PointInRing::Inside: return "Inside";
    }
    return "?";
}

}  // namespace terrain::test
