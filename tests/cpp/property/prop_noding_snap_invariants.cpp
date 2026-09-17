// Property tests for increment 5a: the laws that must hold for every point,
// every segment pair and every node set, rather than for the named
// configurations the unit suites pin.
//
// INVARIANT-CRITICAL, mutation round (docs/increments/05-noder.md, "The
// invariant-critical suite (5a)").
//
// No rapidcheck. Catch2's GENERATE over a fixed seed range plus a seeded
// std::mt19937_64 gives reproducible input without a second framework; a
// failure prints its seed as the generator index and rerunning that section
// reproduces it exactly. Same arrangement as prop_ring_invariants.cpp.
//
// FastKernel IS FORBIDDEN IN THIS FILE, and the ban has two reasons. The
// oracle below is built from the kernel, so an oracle less exact than its
// subject reports false failures on precisely the degenerate inputs the
// property is about; and snapped data is where the filter's fall-throughs live
// -- 38 % of grid-collinear triples at 0.1 m -- so FastKernel is at its worst
// on exactly this increment's inputs. DefaultKernel only, everywhere.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/noding/node_set.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

using terrain::GridPoint;
using terrain::Point2;
using terrain::Segment2;
using terrain::SnapGrid;
using terrain::on_segment;
using terrain::reversed;
using terrain::noding::NodeSet;
using terrain::noding::SegmentRelation;
using terrain::noding::classify;
using terrain::noding::crossing_point;
using terrain::noding::segment_meets_cell;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;

namespace {

constexpr int seed_count = 24;

constexpr double decimal_spacing = 0.1;
constexpr double dyadic_spacing = 0.0625;  // 2^-4

// UTM33 easting/northing, the magnitudes the engine actually runs at.
constexpr double base_x = 500000.0;
constexpr double base_y = 7900000.0;

[[nodiscard]] std::mt19937_64 seeded(int seed) {
    return std::mt19937_64{0x5eed5a00ULL + static_cast<std::uint64_t>(seed)};
}

// One ulp of x, as a positive number. The displacement and containment bounds
// are geometric statements about half a cell; the slack below is the arithmetic
// of evaluating them at UTM33 magnitudes, where fl(ix * spacing) is up to half
// an ulp of the coordinate away from the exact lattice coordinate. It is not a
// tolerance on the geometry -- it is the width of the representation.
[[nodiscard]] double ulp_of(double x) {
    const double a = std::abs(x);
    return std::nextafter(a, std::numeric_limits<double>::infinity()) - a;
}

[[nodiscard]] Point2 random_point(std::mt19937_64& rng, double half_extent) {
    std::uniform_real_distribution<double> d{-half_extent, half_extent};
    return Point2{base_x + d(rng), base_y + d(rng)};
}

// Endpoints on a coarse integer lattice, so exact collinearity, shared
// endpoints and exact overlaps occur often enough for the classification
// property to exercise every arm rather than just Disjoint.
[[nodiscard]] Point2 lattice_point(std::mt19937_64& rng, int extent) {
    std::uniform_int_distribution<int> d{-extent, extent};
    return Point2{base_x + static_cast<double>(d(rng)), base_y + static_cast<double>(d(rng))};
}

[[nodiscard]] Segment2 lattice_segment(std::mt19937_64& rng, int extent) {
    return Segment2{lattice_point(rng, extent), lattice_point(rng, extent)};
}

// The oracle for classify, written out longhand from the four orientations.
//
// IT MAY NOT CALL on_segment<K>: classify's endpoint arm calls it, so an oracle
// that shared it would be tautological on exactly the exact-incidence cases the
// property is about. Betweenness is spelled out here in the four closed
// comparisons, with the collinearity it is conditioned on coming from this
// function's own orient2d call. The only thing shared with the subject is K
// itself, which is what is being tested rather than a component of the test.
[[nodiscard]] bool oracle_between(const Segment2& s, const Point2& p) {
    if (DefaultKernel::orient2d(s.a, s.b, p) != Orientation::Collinear) {
        return false;
    }
    return std::min(s.a.x, s.b.x) <= p.x && p.x <= std::max(s.a.x, s.b.x) &&
           std::min(s.a.y, s.b.y) <= p.y && p.y <= std::max(s.a.y, s.b.y);
}

[[nodiscard]] SegmentRelation oracle_classify(const Segment2& s, const Segment2& t) {
    const Orientation o1 = DefaultKernel::orient2d(s.a, s.b, t.a);
    const Orientation o2 = DefaultKernel::orient2d(s.a, s.b, t.b);
    const Orientation o3 = DefaultKernel::orient2d(t.a, t.b, s.a);
    const Orientation o4 = DefaultKernel::orient2d(t.a, t.b, s.b);

    const bool all_collinear = o1 == Orientation::Collinear && o2 == Orientation::Collinear &&
                               o3 == Orientation::Collinear && o4 == Orientation::Collinear;

    if (all_collinear) {
        const double xlo = std::max(std::min(s.a.x, s.b.x), std::min(t.a.x, t.b.x));
        const double xhi = std::min(std::max(s.a.x, s.b.x), std::max(t.a.x, t.b.x));
        const double ylo = std::max(std::min(s.a.y, s.b.y), std::min(t.a.y, t.b.y));
        const double yhi = std::min(std::max(s.a.y, s.b.y), std::max(t.a.y, t.b.y));

        if (xlo > xhi || ylo > yhi) {
            return SegmentRelation::Disjoint;
        }
        if (xlo == xhi && ylo == yhi) {
            return SegmentRelation::Touching;
        }
        return SegmentRelation::Overlapping;
    }

    if (o1 != o2 && o3 != o4) {
        const bool any_collinear = o1 == Orientation::Collinear || o2 == Orientation::Collinear ||
                                   o3 == Orientation::Collinear || o4 == Orientation::Collinear;
        return any_collinear ? SegmentRelation::Touching : SegmentRelation::Crossing;
    }

    if (oracle_between(s, t.a) || oracle_between(s, t.b) || oracle_between(t, s.a) ||
        oracle_between(t, s.b)) {
        return SegmentRelation::Touching;
    }
    return SegmentRelation::Disjoint;
}

// Squared distance from p to the closed segment s, by clamped projection. Used
// by the containment property only, at a dyadic spacing and at local
// magnitudes, where every product below is exact in double.
[[nodiscard]] double squared_distance(const Segment2& s, const Point2& p) {
    const Point2 d = s.b - s.a;
    const double len2 = dot(d, d);
    if (len2 == 0.0) {
        return dot(p - s.a, p - s.a);
    }
    const double u = std::clamp(dot(p - s.a, d) / len2, 0.0, 1.0);
    const Point2 q{s.a.x + u * d.x, s.a.y + u * d.y};
    return dot(p - q, p - q);
}

}  // namespace

// Displacement: snapping moves a point by at most half a cell per axis, which
// is the numeric content of testing.md's `sqrt(2)*h` bullet -- a vertex moves by
// at most half a cell diagonal, h/sqrt(2), and per axis that is h/2.
TEST_CASE("snapping displaces a point by at most half a cell per axis", "[noding][property][snap]") {
    const int seed = GENERATE(range(0, seed_count));
    const double spacing = GENERATE(decimal_spacing, dyadic_spacing, 1.0);
    auto rng = seeded(seed);
    const SnapGrid grid{spacing};
    const double half = spacing * 0.5;

    for (int i = 0; i < 64; ++i) {
        const Point2 p = random_point(rng, 1000.0);
        REQUIRE(grid.can_snap(p));

        const Point2 s = grid.snapped(p);
        INFO("seed " << seed << " spacing " << spacing << " p = (" << p.x << ", " << p.y << ")");
        REQUIRE(std::abs(s.x - p.x) <= half + ulp_of(p.x));
        REQUIRE(std::abs(s.y - p.y) <= half + ulp_of(p.y));
    }
}

TEST_CASE("snapping is idempotent", "[noding][property][snap]") {
    const int seed = GENERATE(range(0, seed_count));
    const double spacing = GENERATE(decimal_spacing, dyadic_spacing);
    auto rng = seeded(seed);
    const SnapGrid grid{spacing};

    std::uniform_int_distribution<std::int64_t> index{-100000000, 100000000};

    for (int i = 0; i < 64; ++i) {
        const GridPoint g{index(rng), index(rng)};
        INFO("seed " << seed << " spacing " << spacing << " g = (" << g.ix << ", " << g.iy << ")");
        REQUIRE(grid.snap(grid.world(g)) == g);

        const Point2 p = random_point(rng, 1000.0);
        const Point2 once = grid.snapped(p);
        const Point2 twice = grid.snapped(once);
        REQUIRE(once.x == twice.x);  // bitwise
        REQUIRE(once.y == twice.y);
    }
}

// THE PROPERTY DEDUP RESTS ON, in both directions: two points snap to the same
// GridPoint if and only if their snapped world coordinates are bit-identical
// doubles. Pairs are generated both at random and deliberately astride a cell
// boundary, because the boundary is the only place the two directions can come
// apart.
TEST_CASE("the grid key and the snapped coordinate agree on coincidence", "[noding][property][snap]") {
    const int seed = GENERATE(range(0, seed_count));
    const double spacing = GENERATE(decimal_spacing, dyadic_spacing);
    auto rng = seeded(seed);
    const SnapGrid grid{spacing};

    std::uniform_real_distribution<double> jitter{-spacing, spacing};
    std::uniform_int_distribution<int> steps{-4, 4};

    for (int i = 0; i < 64; ++i) {
        const Point2 p = random_point(rng, 1000.0);

        // A near neighbour, and a point placed on the boundary between p's cell
        // and the next one, nudged a few ulps either way.
        const GridPoint g = grid.snap(p);
        double boundary = grid.cell_max(g).x;
        for (int k = 0; k < std::abs(steps(rng)); ++k) {
            boundary = std::nextafter(boundary, std::numeric_limits<double>::infinity());
        }

        const Point2 others[] = {
            Point2{p.x + jitter(rng), p.y + jitter(rng)},
            Point2{boundary, p.y},
            grid.world(g),
        };

        for (const Point2& q : others) {
            INFO("seed " << seed << " spacing " << spacing << " p = (" << p.x << ", " << p.y
                         << ") q = (" << q.x << ", " << q.y << ")");
            const bool same_key = grid.snap(p) == grid.snap(q);
            const bool same_coord =
                grid.snapped(p).x == grid.snapped(q).x && grid.snapped(p).y == grid.snapped(q).y;
            REQUIRE(same_key == same_coord);
        }
    }
}

// Mutant 10, and the shuffle is what kills it. Node ids are in lexicographic
// grid order, which is a function of the SET and of nothing else; first
// appearance order would make the numbering depend on scheduling, and every
// downstream mesh index with it.
TEST_CASE("NodeSet does not depend on the order it was given", "[noding][property][node_set]") {
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    std::uniform_int_distribution<std::int64_t> index{-40, 40};
    std::vector<GridPoint> points;
    points.reserve(200);
    for (int i = 0; i < 200; ++i) {
        points.push_back(GridPoint{index(rng), index(rng)});  // duplicates are the point
    }

    const NodeSet reference{points};

    for (int round = 0; round < 4; ++round) {
        std::shuffle(points.begin(), points.end(), rng);
        const NodeSet shuffled{points};

        INFO("seed " << seed << " round " << round);
        REQUIRE(shuffled.size() == reference.size());
        REQUIRE(std::ranges::equal(shuffled.points(), reference.points()));

        for (const GridPoint& g : reference.points()) {
            REQUIRE(shuffled.id_of(g) == reference.id_of(g));
        }
    }
}

// Totality: classify has an answer for every pair, and it is the answer the
// longhand case analysis gives. Degenerate pairs -- zero-length, collinear,
// shared endpoints -- are common in this generator by construction, so the
// property covers the arms rather than just the generic crossing.
TEST_CASE("classify agrees with the longhand oracle on every pair", "[noding][property][intersect]") {
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);

    for (int i = 0; i < 128; ++i) {
        const Segment2 s = lattice_segment(rng, 6);
        const Segment2 t = lattice_segment(rng, 6);

        INFO("seed " << seed << " s = (" << s.a.x << ", " << s.a.y << ")->(" << s.b.x << ", "
                     << s.b.y << ") t = (" << t.a.x << ", " << t.a.y << ")->(" << t.b.x << ", "
                     << t.b.y << ")");

        const SegmentRelation r = classify<DefaultKernel>(s, t);
        REQUIRE(r == oracle_classify(s, t));

        // The same relation whichever way round, and whichever way along.
        REQUIRE(classify<DefaultKernel>(t, s) == r);
        REQUIRE(classify<DefaultKernel>(reversed(s), reversed(t)) == r);
    }
}

// The clamp, as a universal invariant rather than as a found fixture (mutant 9):
// on EVERY crossing, the constructed node's closed cell meets the intersection
// of the two segments' coordinate ranges, which is the region the true
// intersection provably lies in.
TEST_CASE("a constructed crossing lands in the coordinate ranges it must", "[noding][property][intersect]") {
    const int seed = GENERATE(range(0, seed_count));
    const double spacing = GENERATE(decimal_spacing, dyadic_spacing);
    auto rng = seeded(seed);
    const SnapGrid grid{spacing};

    int crossings = 0;
    for (int i = 0; i < 256; ++i) {
        const Segment2 s = lattice_segment(rng, 6);
        const Segment2 t = lattice_segment(rng, 6);
        if (classify<DefaultKernel>(s, t) != SegmentRelation::Crossing) {
            continue;
        }
        ++crossings;

        const GridPoint g = crossing_point<DefaultKernel>(grid, s, t);
        const Point2 lo = grid.cell_min(g);
        const Point2 hi = grid.cell_max(g);

        const double xlo = std::max(std::min(s.a.x, s.b.x), std::min(t.a.x, t.b.x));
        const double xhi = std::min(std::max(s.a.x, s.b.x), std::max(t.a.x, t.b.x));
        const double ylo = std::max(std::min(s.a.y, s.b.y), std::min(t.a.y, t.b.y));
        const double yhi = std::min(std::max(s.a.y, s.b.y), std::max(t.a.y, t.b.y));

        INFO("seed " << seed << " spacing " << spacing << " crossing " << crossings);
        REQUIRE(lo.x <= xhi);
        REQUIRE(hi.x >= xlo);
        REQUIRE(lo.y <= yhi);
        REQUIRE(hi.y >= ylo);
    }

    // A property that silently generated no crossings would pass while testing
    // nothing.
    REQUIRE(crossings > 0);
}

// HOT-PIXEL DOMINANCE, and the shape of this test is the finding.
//
// on_segment<DefaultKernel>(s, p) implies segment_meets_cell(grid, s, snap(p)):
// exact incidence is strictly stronger than hot-pixel proximity, so the
// hot-pixel predicate never misses an incidence.
//
// THE CONVERSE IS NOT ASSERTED, AND MUST NOT BE. segment_meets_cell is true for
// a great many cells whose lattice point is not exactly on the segment, and
// that gap is precisely what makes snap rounding work where exact incidence
// does not. What is asserted about the gap is only that it is NON-EMPTY on this
// generator -- a run in which the two predicates agreed everywhere would mean
// the generator had stopped producing the case the increment exists for.
TEST_CASE("exact incidence implies hot-pixel incidence, and not the reverse",
          "[noding][property][hot_pixel]") {
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);
    const SnapGrid grid{decimal_spacing};

    std::uniform_int_distribution<std::int64_t> index{-2000, 2000};
    std::uniform_int_distribution<std::int64_t> step{1, 40};
    std::uniform_int_distribution<int> multiple{1, 8};

    const GridPoint origin = grid.snap(Point2{base_x, base_y});
    int gap = 0;

    for (int i = 0; i < 128; ++i) {
        // A host whose endpoints and one interior node are collinear IN INDEX
        // SPACE, plus an unrelated point. The first is where the two predicates
        // come apart; the second is where they agree trivially.
        const GridPoint a{origin.ix + index(rng), origin.iy + index(rng)};
        const std::int64_t dx = step(rng);
        const std::int64_t dy = step(rng);
        const int k = multiple(rng);
        const GridPoint mid{a.ix + k * dx, a.iy + k * dy};
        const GridPoint b{a.ix + 2 * k * dx, a.iy + 2 * k * dy};

        const Segment2 host{grid.world(a), grid.world(b)};

        const GridPoint probes[] = {a, mid, b, GridPoint{a.ix + index(rng), a.iy + index(rng)}};
        for (const GridPoint& g : probes) {
            const Point2 p = grid.world(g);
            INFO("seed " << seed << " g = (" << g.ix << ", " << g.iy << ")");

            const bool exact = on_segment<DefaultKernel>(host, p);
            const bool hot = segment_meets_cell<DefaultKernel>(grid, host, g);
            if (exact) {
                REQUIRE(hot);
            }
            if (hot && !exact) {
                ++gap;
            }
        }
    }

    REQUIRE(gap > 0);
}

// HOT-PIXEL CONTAINMENT, in two forms, because only one of them can be stated
// without arithmetic of the test's own.
//
// The exact form is the one the predicate's first step guarantees and it is
// asserted with no tolerance at all: if the cell meets the segment then the
// segment's coordinate range reaches the cell on both axes.
//
// The geometric form -- world(g) is within spacing/sqrt(2) of s, which is the
// bound 5b's deformation inherits -- needs a distance, and a distance needs
// division and a square. It is asserted at a DYADIC spacing and at local
// magnitudes, where every product in squared_distance is exact in double, with
// a slack of one ulp of the bound rather than an epsilon chosen to make it pass.
TEST_CASE("a cell that meets a segment is within half a cell diagonal of it",
          "[noding][property][hot_pixel]") {
    const int seed = GENERATE(range(0, seed_count));
    auto rng = seeded(seed);
    const SnapGrid grid{dyadic_spacing};

    std::uniform_int_distribution<std::int64_t> index{-4096, 4096};
    const double bound2 = dyadic_spacing * dyadic_spacing * 0.5;  // (spacing/sqrt(2))^2, exact
    int met = 0;

    for (int i = 0; i < 128; ++i) {
        const Segment2 s{grid.world(GridPoint{index(rng), index(rng)}),
                         grid.world(GridPoint{index(rng), index(rng)})};

        // Probe the neighbourhood of the segment, where the answer is
        // interesting, rather than the whole domain, where it is always false.
        const GridPoint anchor = grid.snap(s.a);
        std::uniform_int_distribution<std::int64_t> near{-3, 3};

        for (int j = 0; j < 8; ++j) {
            const GridPoint g{anchor.ix + near(rng), anchor.iy + near(rng)};
            if (!segment_meets_cell<DefaultKernel>(grid, s, g)) {
                continue;
            }
            ++met;

            INFO("seed " << seed << " g = (" << g.ix << ", " << g.iy << ")");

            // Exact, per axis, no tolerance.
            REQUIRE(std::min(s.a.x, s.b.x) <= grid.cell_max(g).x);
            REQUIRE(std::max(s.a.x, s.b.x) >= grid.cell_min(g).x);
            REQUIRE(std::min(s.a.y, s.b.y) <= grid.cell_max(g).y);
            REQUIRE(std::max(s.a.y, s.b.y) >= grid.cell_min(g).y);

            // Geometric, at a dyadic spacing where the arithmetic is exact.
            const double d2 = squared_distance(s, grid.world(g));
            REQUIRE(d2 <= std::nextafter(bound2, std::numeric_limits<double>::infinity()));
        }
    }

    REQUIRE(met > 0);
}
