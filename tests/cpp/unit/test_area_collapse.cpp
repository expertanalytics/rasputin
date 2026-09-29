// Increment 22, PR 2 (docs/increments/22-auto-catchment.md, "The reduction"
// and "The red suites", PR 2 C++): the fine outline reduced by
// area-preserving segment collapse (Kronenfeld et al. 2020) to a horizontal
// tolerance, simple and with the seed kept inside.
//
// The design names this file tests/cpp/unit/area_collapse.cpp; it carries the
// test_ prefix every other suite in tests/cpp/unit/ has.
//
// Interface assumed (the design fixes the header, the function, the statuses
// and which counts exist; the spelling is chosen here and stated in the
// handback):
//
//   #include <terrain/vector_simplify/area_collapse.hpp>
//   namespace terrain::vector_simplify {
//   enum class ReduceStatus : std::uint8_t {
//       Ok, InvalidTolerance, NotCounterClockwise, TooFewVertices };
//   struct ReduceCounts {
//       std::size_t collinear;          // vertices dropped by the collinear pass
//       std::size_t collapses;          // B, C -> E replacements made
//       std::size_t rejected_crossing, rejected_seed, rejected_tolerance;
//   };
//   struct ReduceOutcome {
//       std::vector<Point2> ring;       // open: the first vertex not repeated
//       ReduceStatus status;
//       ReduceCounts counts;
//   };
//   template <pred::GeometryKernel K>
//   ReduceOutcome reduce_ring(std::span<const Point2> ring, double tolerance,
//                             std::span<const Point2> keep);
//   }
//
// The input ring is open too. On a status other than Ok the ring's content is
// not pinned.
//
// THE CHECKS ARE THE DESIGN'S GUARANTEES, computed here from the output and
// the input alone: area to 1e-9 relative; simple by a brute-force pairwise
// classify<DefaultKernel> (non-adjacent edges Disjoint, adjacent ones
// Touching); keep-points strictly inside by an exact winding number; every
// input vertex within the tolerance of the output ring and every output
// vertex within it of the input ring; the same input twice gives equal bits.
// Plus the bookkeeping identity: input size - collinear - collapses = output
// size (each collapse replaces two vertices by one).

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/vector_simplify/area_collapse.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <numeric>
#include <random>
#include <span>
#include <vector>

using terrain::Point2;
using terrain::Segment2;
using terrain::noding::classify;
using terrain::noding::SegmentRelation;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::vector_simplify::reduce_ring;
using terrain::vector_simplify::ReduceOutcome;
using terrain::vector_simplify::ReduceStatus;

namespace {

using Ring = std::vector<Point2>;

ReduceOutcome reduce(const Ring& ring, double tolerance, const Ring& keep = {}) {
    return reduce_ring<DefaultKernel>(std::span<const Point2>{ring}, tolerance,
                                      std::span<const Point2>{keep});
}

double area(const Ring& r) {
    double twice = 0.0;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const Point2& a = r[i];
        const Point2& b = r[(i + 1) % r.size()];
        twice += a.x * b.y - b.x * a.y;
    }
    return 0.5 * twice;
}

Segment2 edge(const Ring& r, std::size_t i) { return Segment2{r[i], r[(i + 1) % r.size()]}; }

// Brute force, exact: adjacent edges share only their vertex, others nothing.
bool simple(const Ring& r) {
    const std::size_t n = r.size();
    if (n < 3)
        return false;
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j) {
            const bool adjacent = j == i + 1 || (i == 0 && j == n - 1);
            const SegmentRelation rel = classify<DefaultKernel>(edge(r, i), edge(r, j));
            if (adjacent ? rel != SegmentRelation::Touching : rel != SegmentRelation::Disjoint)
                return false;
        }
    return true;
}

// Strictly inside: on no edge, and an exact winding number of 1.
bool strictly_inside(const Point2& p, const Ring& r) {
    int winding = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const Point2& a = r[i];
        const Point2& b = r[(i + 1) % r.size()];
        if (classify<DefaultKernel>(Segment2{p, p}, Segment2{a, b}) != SegmentRelation::Disjoint)
            return false;
        const Orientation o = DefaultKernel::orient2d(a, b, p);
        if (a.y <= p.y && b.y > p.y && o == Orientation::CounterClockwise)
            ++winding;
        if (a.y > p.y && b.y <= p.y && o == Orientation::Clockwise)
            --winding;
    }
    return winding == 1;
}

double distance(const Point2& p, const Point2& a, const Point2& b) {
    const double dx = b.x - a.x, dy = b.y - a.y;
    const double len2 = dx * dx + dy * dy;
    double t = len2 > 0.0 ? ((p.x - a.x) * dx + (p.y - a.y) * dy) / len2 : 0.0;
    t = std::clamp(t, 0.0, 1.0);
    return std::hypot(p.x - (a.x + t * dx), p.y - (a.y + t * dy));
}

double distance(const Point2& p, const Ring& r) {
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < r.size(); ++i)
        best = std::min(best, distance(p, r[i], r[(i + 1) % r.size()]));
    return best;
}

// Every guarantee of the design, on one call.
void check_guarantees(const Ring& fine, double tolerance, const Ring& keep, const ReduceOutcome& out) {
    REQUIRE(out.status == ReduceStatus::Ok);
    const Ring& red = out.ring;
    REQUIRE(red.size() >= 3);
    CHECK(fine.size() - out.counts.collinear - out.counts.collapses == red.size());
    const double a0 = area(fine);
    CHECK(std::abs(area(red) - a0) <= 1e-9 * a0);
    CHECK(area(red) > 0.0); // still counter-clockwise
    CHECK(simple(red));
    for (const Point2& k : keep)
        CHECK(strictly_inside(k, red));
    const double slack = 1e-9 * std::max(1.0, tolerance);
    for (const Point2& p : fine)
        CHECK(distance(p, red) <= tolerance + slack);
    for (const Point2& p : red)
        CHECK(distance(p, fine) <= tolerance + slack);
    for (const Point2& p : red)
        CHECK((std::isfinite(p.x) && std::isfinite(p.y)));
}

bool same_bits(const Ring& a, const Ring& b) {
    return a.size() == b.size()
           && std::equal(a.begin(), a.end(), b.begin(), [](const Point2& p, const Point2& q) {
                  return std::bit_cast<std::uint64_t>(p.x) == std::bit_cast<std::uint64_t>(q.x)
                         && std::bit_cast<std::uint64_t>(p.y) == std::bit_cast<std::uint64_t>(q.y);
              });
}

// The fine outline of a w x h block of nodes (cell `cell`, first node at
// `origin`), as the tracer draws it: a vertex at the midpoint of every lattice
// edge leaving the block, counter-clockwise. 2 (w + h) vertices, of which all
// but the 8 corner-cutting ones are collinear; area (w h - 1/2) cells.
Ring block(std::size_t w, std::size_t h, double cell = 10.0, Point2 origin = {1000.0, 2000.0}) {
    Ring r;
    const auto at = [&](double i, double j) { return Point2{origin.x + i * cell, origin.y + j * cell}; };
    const double W = static_cast<double>(w), H = static_cast<double>(h);
    for (std::size_t i = 0; i < w; ++i)
        r.push_back(at(static_cast<double>(i), -0.5));
    for (std::size_t j = 0; j < h; ++j)
        r.push_back(at(W - 0.5, static_cast<double>(j)));
    for (std::size_t i = w; i-- > 0;)
        r.push_back(at(static_cast<double>(i), H - 0.5));
    for (std::size_t j = h; j-- > 0;)
        r.push_back(at(-0.5, static_cast<double>(j)));
    return r;
}

// The same up to where the cycle starts.
bool same_cycle(const Ring& a, const Ring& b) {
    if (a.size() != b.size())
        return false;
    for (std::size_t s = 0; s < a.size(); ++s) {
        bool all = true;
        for (std::size_t i = 0; i < a.size() && all; ++i)
            all = a[(s + i) % a.size()] == b[i];
        if (all)
            return true;
    }
    return a.empty();
}

Ring corners_of(const Ring& r) {
    Ring out;
    for (std::size_t i = 0; i < r.size(); ++i) {
        const Point2& a = r[(i + r.size() - 1) % r.size()];
        const Point2& b = r[(i + 1) % r.size()];
        if (DefaultKernel::orient2d(a, r[i], b) != Orientation::Collinear)
            out.push_back(r[i]);
    }
    return out;
}

// A random star polygon round `centre`: sorted distinct angles, radii on a
// half-metre grid, so it is simple and the centre is strictly inside.
Ring star(std::mt19937& rng, std::size_t n, Point2 centre = {5000.0, 7000.0}) {
    std::uniform_int_distribution<int> radius(100, 300);
    std::vector<int> ticks(3600);
    std::iota(ticks.begin(), ticks.end(), 0);
    std::shuffle(ticks.begin(), ticks.end(), rng);
    ticks.resize(n);
    std::sort(ticks.begin(), ticks.end());
    Ring r;
    for (int t : ticks) {
        const double angle = 2.0 * std::numbers::pi * t / 3600.0;
        const double rho = 0.5 * radius(rng);
        r.push_back(Point2{std::round(centre.x + rho * std::cos(angle)),
                           std::round(centre.y + rho * std::sin(angle))});
    }
    // Rounding may make three consecutive vertices collinear or two equal on
    // a small star; the fixture keeps only simple, counter-clockwise ones.
    return r;
}

} // namespace

TEST_CASE("status: an invalid tolerance, a clockwise ring, fewer than 4 vertices",
          "[vector_simplify][area_collapse][status]") {
    const Ring ring = block(5, 4);
    CHECK(reduce(ring, -1.0).status == ReduceStatus::InvalidTolerance);
    CHECK(reduce(ring, std::numeric_limits<double>::quiet_NaN()).status
          == ReduceStatus::InvalidTolerance);
    CHECK(reduce(ring, std::numeric_limits<double>::infinity()).status
          == ReduceStatus::InvalidTolerance);
    CHECK(reduce(ring, -std::numeric_limits<double>::infinity()).status
          == ReduceStatus::InvalidTolerance);

    Ring clockwise(ring.rbegin(), ring.rend());
    CHECK(reduce(clockwise, 10.0).status == ReduceStatus::NotCounterClockwise);

    const Ring triangle{{0.0, 0.0}, {10.0, 0.0}, {0.0, 10.0}};
    CHECK(reduce(triangle, 10.0).status == ReduceStatus::TooFewVertices);
    CHECK(reduce(Ring{{0.0, 0.0}, {10.0, 0.0}}, 10.0).status == ReduceStatus::TooFewVertices);
    CHECK(reduce(Ring{}, 10.0).status == ReduceStatus::TooFewVertices);

    const Ring square{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {0.0, 10.0}};
    CHECK(reduce(square, 10.0).status == ReduceStatus::Ok);
}

TEST_CASE("tolerance 0 drops only the exactly collinear vertices", "[vector_simplify][area_collapse][zero]") {
    const Ring ring = block(20, 10);
    REQUIRE(ring.size() == 60);
    const ReduceOutcome out = reduce(ring, 0.0);
    REQUIRE(out.status == ReduceStatus::Ok);
    CHECK(out.counts.collinear == 52);
    CHECK(out.counts.collapses == 0);
    CHECK(same_cycle(out.ring, corners_of(ring)));
    CHECK(out.ring.size() == 8);
    CHECK(area(out.ring) == area(ring));

    // A ring with no collinear vertex comes back as it went in.
    const Ring square{{0.0, 0.0}, {10.0, 0.0}, {10.0, 10.0}, {5.0, 12.5}, {0.0, 10.0}};
    const ReduceOutcome same = reduce(square, 0.0);
    REQUIRE(same.status == ReduceStatus::Ok);
    CHECK(same_cycle(same.ring, square));
    CHECK(same.counts.collinear == 0);
}

TEST_CASE("a traced rectangle of nodes reduces to at most 8 vertices at one cell, area unchanged",
          "[vector_simplify][area_collapse][rectangle]") {
    constexpr std::array<std::array<std::size_t, 2>, 4> shapes{{{20, 10}, {3, 17}, {2, 2}, {40, 40}}};
    for (const auto& [w, h] : shapes) {
        const Ring ring = block(w, h);
        const ReduceOutcome out = reduce(ring, 10.0);
        INFO(w << " x " << h);
        check_guarantees(ring, 10.0, {}, out);
        CHECK(out.ring.size() <= 8);
        CHECK(out.ring.size() >= 4);
    }
}

TEST_CASE("a keep-point half a cell inside a tab a collapse would cut off stays inside",
          "[vector_simplify][area_collapse][keep]") {
    // A 200 m square with a tab one cell (10 m) wide and 1.5 cells tall on its
    // top side, densified at 10 m so there is plenty to collapse. At a
    // tolerance of 5 cells the tab is well within reach of a collapse; the
    // keep-point is half a cell inside the tab on every side.
    Ring ring;
    for (int x = 0; x < 200; x += 10)
        ring.push_back({static_cast<double>(x), 0.0});
    for (int y = 0; y < 200; y += 10)
        ring.push_back({200.0, static_cast<double>(y)});
    for (int x = 200; x > 110; x -= 10)
        ring.push_back({static_cast<double>(x), 200.0});
    ring.push_back({110.0, 200.0});
    ring.push_back({110.0, 215.0});
    ring.push_back({100.0, 215.0});
    ring.push_back({100.0, 200.0});
    for (int x = 90; x > 0; x -= 10)
        ring.push_back({static_cast<double>(x), 200.0});
    for (int y = 200; y > 0; y -= 10)
        ring.push_back({0.0, static_cast<double>(y)});
    REQUIRE(simple(ring));
    REQUIRE(area(ring) > 0.0);

    const Ring keep{{105.0, 210.0}};
    REQUIRE(strictly_inside(keep[0], ring));
    const ReduceOutcome out = reduce(ring, 50.0, keep);
    check_guarantees(ring, 50.0, keep, out);
}

TEST_CASE("a thin corridor at a tolerance of five cells does not cross itself",
          "[vector_simplify][area_collapse][corridor]") {
    // Two 100 m squares joined by a corridor 300 m long whose sides are 1.5
    // cells (15 m) apart, each side a zigzag of 2.5 m so every vertex is a
    // candidate. A collapse on one side that reached the other would cross.
    Ring ring{{0.0, 0.0}, {100.0, 0.0}};
    for (int k = 0; k <= 30; ++k)
        ring.push_back({100.0 + 10.0 * k, 42.5 + 2.5 * (k % 2)});
    ring.push_back({400.0, 0.0});
    ring.push_back({500.0, 0.0});
    ring.push_back({500.0, 100.0});
    ring.push_back({400.0, 100.0});
    for (int k = 30; k >= 0; --k)
        ring.push_back({100.0 + 10.0 * k, 57.5 - 2.5 * (k % 2)});
    ring.push_back({100.0, 100.0});
    ring.push_back({0.0, 100.0});
    REQUIRE(simple(ring));
    REQUIRE(area(ring) > 0.0);

    const Ring keep{{50.0, 50.0}, {450.0, 50.0}};
    const ReduceOutcome out = reduce(ring, 50.0, keep);
    check_guarantees(ring, 50.0, keep, out);
    CHECK(out.ring.size() < ring.size());
}

TEST_CASE("the same input twice gives equal bits", "[vector_simplify][area_collapse][determinism]") {
    std::mt19937 rng{22201};
    for (int trial = 0; trial < 20; ++trial) {
        const Ring ring = star(rng, 40);
        if (!simple(ring) || !(area(ring) > 0.0))
            continue;
        const Ring keep{{5000.0, 7000.0}};
        const ReduceOutcome a = reduce(ring, 20.0, keep);
        const ReduceOutcome b = reduce(ring, 20.0, keep);
        CHECK(a.status == b.status);
        CHECK(same_bits(a.ring, b.ring));
        CHECK(a.counts.collapses == b.counts.collapses);
    }
    const Ring rect = block(30, 7);
    CHECK(same_bits(reduce(rect, 10.0).ring, reduce(rect, 10.0).ring));
}

TEST_CASE("random star polygons keep every guarantee", "[vector_simplify][area_collapse][property]") {
    std::mt19937 rng{22202};
    std::uniform_int_distribution<std::size_t> size(8, 80);
    const double tolerances[] = {0.0, 2.0, 10.0, 30.0, 100.0};
    int checked = 0, reduced = 0;
    for (int trial = 0; trial < 300; ++trial) {
        const Ring ring = star(rng, size(rng));
        if (!simple(ring) || !(area(ring) > 0.0))
            continue;
        const Ring keep{{5000.0, 7000.0}};
        if (!strictly_inside(keep[0], ring))
            continue;
        const double tolerance = tolerances[trial % 5];
        const ReduceOutcome out = reduce(ring, tolerance, keep);
        INFO("trial " << trial << ", " << ring.size() << " vertices, tolerance " << tolerance);
        check_guarantees(ring, tolerance, keep, out);
        if (tolerance == 0.0)
            CHECK(out.counts.collapses == 0);
        ++checked;
        reduced += out.counts.collapses > 0 ? 1 : 0;
    }
    // The fixture must exercise the loop, not only the collinear pass.
    CHECK(checked > 200);
    CHECK(reduced > 100);
}
