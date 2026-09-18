// Unit tests for terrain/noding/intersect.hpp -- classify, crossing_point and
// the hot-pixel predicate segment_meets_cell.
//
// INVARIANT-CRITICAL, mutation round (docs/increments/05-noder.md, "The
// invariant-critical suite (5a)").
//
// Three rules from the design govern what this file is allowed to assert, and
// each is repeated at the fixture it constrains:
//
//  1. No fixture may assert that a `Touching` answer came from classify's
//     ENDPOINT ARM under DefaultKernel, because none does -- under an exact
//     kernel that arm's only reachable answer is Disjoint. Its coverage here is
//     the zero-length-segment-off-the-other's-line fixture, and that fixture
//     says so.
//  2. `Overlapping` may be asserted only where the collinear arm provably
//     fires: axis-aligned pairs, or any direction at a dyadic spacing. A 0.1 m
//     diagonal overlap may assert only `Crossing || Disjoint`.
//  3. The hot-pixel block asserts `on_segment => segment_meets_cell` and never
//     the converse. The gap between the two predicates is this increment's
//     finding, not a defect to be papered over.

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <string>

using terrain::GridPoint;
using terrain::Point2;
using terrain::Segment2;
using terrain::SnapGrid;
using terrain::on_segment;
using terrain::reversed;
using terrain::noding::SegmentRelation;
using terrain::noding::classify;
using terrain::noding::crossing_point;
using terrain::noding::segment_meets_cell;
using terrain::pred::DefaultKernel;
using terrain::pred::FastKernel;

namespace {

constexpr double decimal_spacing = 0.1;
constexpr double dyadic_spacing = 0.0625;  // 2^-4

[[nodiscard]] Segment2 seg(double ax, double ay, double bx, double by) {
    return Segment2{Point2{ax, ay}, Point2{bx, by}};
}

[[nodiscard]] SegmentRelation rel(const Segment2& s, const Segment2& t) {
    return classify<DefaultKernel>(s, t);
}

[[nodiscard]] std::string describe(const Segment2& s) {
    return "(" + std::to_string(s.a.x) + ", " + std::to_string(s.a.y) + ") -> (" +
           std::to_string(s.b.x) + ", " + std::to_string(s.b.y) + ")";
}

// classify answers a question about two unordered closed segments, so it must
// be blind to argument order and to each segment's direction. A relation that
// was not would make the broad phase's pair ordering load-bearing, which is
// exactly the kind of hidden coupling 5b cannot afford. Every named fixture
// below goes through this, so the symmetry is asserted 5 ways per fixture
// rather than once in a test of its own.
void require_relation(const Segment2& s, const Segment2& t, SegmentRelation expected) {
    INFO("s = " << describe(s) << "\nt = " << describe(t));
    REQUIRE(classify<DefaultKernel>(s, t) == expected);
    REQUIRE(classify<DefaultKernel>(t, s) == expected);
    REQUIRE(classify<DefaultKernel>(reversed(s), t) == expected);
    REQUIRE(classify<DefaultKernel>(s, reversed(t)) == expected);
    REQUIRE(classify<DefaultKernel>(reversed(s), reversed(t)) == expected);
}

// The segment_meets_cell fixtures get the same treatment on the one symmetry
// the predicate has: reversing the segment cannot change whether it meets a
// cell.
void require_meets(const SnapGrid& grid, const Segment2& s, const GridPoint& g, bool expected) {
    INFO("s = " << describe(s) << "\ng = (" << g.ix << ", " << g.iy << ")");
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, s, g) == expected);
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, reversed(s), g) == expected);
}

}  // namespace

// ---------------------------------------------------------------------------
// classify: one named test per relation, plus the pairs that separate the arms
// ---------------------------------------------------------------------------

// The single most likely defect in the collinear arm, and the most common pair
// in the whole input: two consecutive edges of a ring. An arm that returns
// Overlapping on sight would merge every ring edge with its neighbour. Mutant 7.
TEST_CASE("collinear segments sharing one endpoint touch, not overlap", "[noding][intersect]") {
    require_relation(seg(0, 0, 1, 0), seg(1, 0, 2, 0), SegmentRelation::Touching);
    require_relation(seg(0, 0, 0, 1), seg(0, 1, 0, 2), SegmentRelation::Touching);
    require_relation(seg(0, 0, 3, 3), seg(3, 3, 5, 5), SegmentRelation::Touching);
}

// Two segments meeting at a shared endpoint but not collinear. The Collinear
// check inside the four-orientation arm is what decides Touching over Crossing
// here; mutant 6 swaps them.
TEST_CASE("segments meeting at a shared endpoint touch, not cross", "[noding][intersect]") {
    require_relation(seg(0, 0, 1, 1), seg(0, 0, 1, -1), SegmentRelation::Touching);
    require_relation(seg(-2, 0, 0, 0), seg(0, 0, 0, 5), SegmentRelation::Touching);
}

// A T-junction in exact arithmetic: one endpoint interior to the other segment.
// This is the four-orientation arm's Touching and it is assertable, unlike the
// endpoint arm's.
TEST_CASE("a T-junction touches", "[noding][intersect]") {
    require_relation(seg(0, 0, 4, 0), seg(2, 0, 2, 3), SegmentRelation::Touching);
    require_relation(seg(0, 0, 4, 4), seg(1, 1, 5, 0), SegmentRelation::Touching);
}

TEST_CASE("segments meeting interior to both cross", "[noding][intersect]") {
    require_relation(seg(-2, -2, 2, 2), seg(-2, 2, 2, -2), SegmentRelation::Crossing);
    require_relation(seg(0, 0, 4, 0), seg(1, -1, 1, 1), SegmentRelation::Crossing);
}

// All four containment arrangements. Axis-aligned, so one coordinate is
// literally equal across the points and the collinear arm provably fires --
// rule 2 above.
TEST_CASE("collinear segments sharing positive length overlap", "[noding][intersect]") {
    SECTION("partial, each way") {
        require_relation(seg(0, 0, 4, 0), seg(2, 0, 6, 0), SegmentRelation::Overlapping);
        require_relation(seg(2, 0, 6, 0), seg(0, 0, 4, 0), SegmentRelation::Overlapping);
    }
    SECTION("one strictly inside the other") {
        require_relation(seg(0, 0, 10, 0), seg(3, 0, 7, 0), SegmentRelation::Overlapping);
    }
    SECTION("identical") {
        require_relation(seg(0, 0, 4, 0), seg(0, 0, 4, 0), SegmentRelation::Overlapping);
    }
    SECTION("sharing an endpoint and extending inwards") {
        require_relation(seg(0, 0, 4, 0), seg(0, 0, 2, 0), SegmentRelation::Overlapping);
    }
}

TEST_CASE("collinear segments sharing no point are disjoint", "[noding][intersect]") {
    require_relation(seg(0, 0, 1, 0), seg(2, 0, 3, 0), SegmentRelation::Disjoint);
    require_relation(seg(0, 0, 1, 1), seg(2, 2, 3, 3), SegmentRelation::Disjoint);
}

// The pair that passes an interval test and must fail the orientation test. An
// implementation that reached the interval intersection without first requiring
// all four orientations Collinear would report Overlapping here.
TEST_CASE("parallel non-collinear segments with overlapping ranges are disjoint", "[noding][intersect]") {
    require_relation(seg(0, 0, 4, 0), seg(1, 1, 5, 1), SegmentRelation::Disjoint);
    require_relation(seg(0, 0, 0, 4), seg(1, 1, 1, 5), SegmentRelation::Disjoint);
}

// The both-axes ruling, and mutant 8. A collinear arm that intersects only the
// x-intervals sees [0,0] against [0,0] for every vertical pair -- a single
// point -- and answers Touching for the overlap and Touching for the disjoint
// pair. Both lines below are wrong under that mutant; nothing else in the file
// is.
TEST_CASE("the collinear arm intersects both axes", "[noding][intersect]") {
    SECTION("vertical") {
        require_relation(seg(0, 0, 0, 4), seg(0, 2, 0, 6), SegmentRelation::Overlapping);
        require_relation(seg(0, 0, 0, 1), seg(0, 2, 0, 3), SegmentRelation::Disjoint);
        require_relation(seg(0, 0, 0, 1), seg(0, 1, 0, 2), SegmentRelation::Touching);
    }
    SECTION("horizontal") {
        require_relation(seg(0, 0, 4, 0), seg(2, 0, 6, 0), SegmentRelation::Overlapping);
        require_relation(seg(0, 0, 1, 0), seg(2, 0, 3, 0), SegmentRelation::Disjoint);
    }
}

// MANDATORY, because there is no branch in the production code to point at: a
// Pslg permits a repeated consecutive index, so a zero-length Segment2 reaches
// classify, and the design's ruling is that it gets no special case. All four
// answers below fall out of the general arms.
TEST_CASE("degenerate segments are classified without a special case", "[noding][intersect]") {
    SECTION("a point on the segment touches -- via the collinear arm") {
        require_relation(seg(2, 0, 2, 0), seg(0, 0, 4, 0), SegmentRelation::Touching);
        require_relation(seg(0, 0, 0, 0), seg(0, 0, 4, 0), SegmentRelation::Touching);  // at an endpoint
        require_relation(seg(2, 2, 2, 2), seg(0, 0, 4, 4), SegmentRelation::Touching);  // diagonal host
    }

    SECTION("a point on the segment's line but past its end is disjoint") {
        // Collinear arm, empty interval intersection.
        require_relation(seg(9, 0, 9, 0), seg(0, 0, 4, 0), SegmentRelation::Disjoint);
    }

    SECTION("a point off the segment's line is disjoint -- THE ENDPOINT ARM") {
        // This fixture is the endpoint arm's coverage under DefaultKernel, and
        // the only one in the file that reaches it. o3 == o4 != Collinear, both
        // decided arms are skipped, on_segment is false four times, and the arm
        // answers Disjoint. What it does NOT do under an exact kernel is return
        // Touching -- see rule 1 at the top of this file, and mutant 18, which
        // no fixture here kills and which the design accepts as defended by
        // case analysis rather than by test.
        require_relation(seg(2, 5, 2, 5), seg(0, 0, 4, 0), SegmentRelation::Disjoint);
    }

    SECTION("two zero-length segments") {
        require_relation(seg(3, 7, 3, 7), seg(3, 7, 3, 7), SegmentRelation::Touching);
        require_relation(seg(3, 7, 3, 7), seg(3, 8, 3, 8), SegmentRelation::Disjoint);
    }
}

// ---------------------------------------------------------------------------
// Ruling 4's two sides: where Overlapping may be asserted, and where it may not
// ---------------------------------------------------------------------------

TEST_CASE("a diagonal overlap at a dyadic spacing overlaps", "[noding][intersect]") {
    // world(g) = g * 2^-4 is an exact affine scaling, so grid-collinearity
    // survives into world coordinates and the collinear arm provably fires.
    // Indices: base + k*(37, 19) for k = 0, 2, 3, 6 at UTM33 magnitudes.
    const SnapGrid grid{dyadic_spacing};
    const GridPoint base{8000000, 126400000};
    const auto at = [&](std::int64_t k) {
        return grid.world(GridPoint{base.ix + k * 37, base.iy + k * 19});
    };

    require_relation(Segment2{at(0), at(3)}, Segment2{at(2), at(6)}, SegmentRelation::Overlapping);
}

TEST_CASE("an axis-aligned overlap at a decimal spacing overlaps", "[noding][intersect]") {
    // One coordinate is literally equal across all four points, so every
    // orientation is exactly zero however the other coordinate rounded.
    const SnapGrid grid{decimal_spacing};
    const auto at = [&](std::int64_t k) { return grid.world(GridPoint{5000000 + k, 79000000}); };

    require_relation(Segment2{at(0), at(20)}, Segment2{at(10), at(35)}, SegmentRelation::Overlapping);
}

// THE FIXTURE THE PROVISIONAL DESIGN WOULD NOT HAVE WRITTEN, and the one that
// stops a later reader "fixing" the collinear arm.
//
// These four points are exactly collinear IN INDEX SPACE. They are not
// collinear in world coordinates, because world(g) = g * 0.1 is not an affine
// map -- fl(ix * 0.1) perturbs each coordinate independently by up to half an
// ulp. So the collinear arm does not fire, the pair falls through to the
// four-orientation arm, and the answer is decided by which way an ulp-level
// perturbation fell. Both Crossing and Disjoint are reachable; asserting either
// one is asserting a coin flip that would pass on the author's machine.
//
// Measured for this fixture: DefaultKernel answers Disjoint. That is recorded
// here as an observation, NOT as the assertion.
TEST_CASE("a diagonal overlap at 0.1 m asserts only what is not a coin flip", "[noding][intersect]") {
    const SnapGrid grid{decimal_spacing};
    const GridPoint base{5000000, 79000000};
    const auto at = [&](std::int64_t k) {
        return grid.world(GridPoint{base.ix + k * 37, base.iy + k * 19});
    };

    const Segment2 s{at(0), at(3)};
    const Segment2 t{at(2), at(6)};
    const SegmentRelation r = rel(s, t);

    REQUIRE((r == SegmentRelation::Crossing || r == SegmentRelation::Disjoint));
    REQUIRE(rel(t, s) == r);  // symmetric whichever way it fell
}

// ---------------------------------------------------------------------------
// crossing_point: the only constructive function in the project
// ---------------------------------------------------------------------------

TEST_CASE("crossing_point returns the lattice point a crossing lands on", "[noding][intersect]") {
    const SnapGrid grid{decimal_spacing};
    const Segment2 s = seg(-2, -2, 2, 2);
    const Segment2 t = seg(-2, 2, 2, -2);

    REQUIRE(rel(s, t) == SegmentRelation::Crossing);  // the documented precondition
    REQUIRE(crossing_point<DefaultKernel>(grid, s, t) == GridPoint{0, 0});
    REQUIRE(crossing_point<DefaultKernel>(grid, t, s) == GridPoint{0, 0});

    // Off the origin, and at a dyadic spacing where the arithmetic is exact.
    const SnapGrid dyadic{dyadic_spacing};
    const Segment2 u = seg(0, 0, 8, 8);
    const Segment2 v = seg(0, 8, 8, 0);
    REQUIRE(crossing_point<DefaultKernel>(dyadic, u, v) == GridPoint{64, 64});  // (4, 4) / 2^-4
}

// THE CLAMP, and mutant 9. The invariant holds for EVERY crossing, but that is
// NOT enough to kill the mutant: a well-conditioned crossing lands in range with
// or without the clamp, and measured, the clamp-removed mutant survives every
// case above. The last fixture in this case is the one that kills it, found by a
// seeded scan. The design's contrary ruling is corrected in `05-noder.md`.
//
// Stated on the node's CELL rather than on its coordinate, and that is exact
// rather than a tolerance: the clamped point provably lies in the intersection
// of the two coordinate ranges, and snapping moves it to a lattice point whose
// closed cell still contains it. Asserting world(g) itself lies in the range
// would be wrong -- snapping is allowed to move it by half a cell, and that is
// the standing contract of snap rounding rather than a defect (risk 4).
TEST_CASE("crossing_point lands in the intersection of the coordinate ranges", "[noding][intersect]") {
    const SnapGrid grid{decimal_spacing};

    const Segment2 cases[][2] = {
        {seg(-2, -2, 2, 2), seg(-2, 2, 2, -2)},
        {seg(0, 0, 10, 1), seg(0, 1, 10, 0)},
        {seg(500000.0, 7900000.0, 500100.0, 7900100.0), seg(500000.0, 7900100.0, 500100.0, 7900000.0)},
        // Near-parallel: den is near zero and the unclamped u is off by orders
        // of magnitude.
        {seg(-20037508.0, -20037508.0, 20037508.0, 20037508.0),
         seg(-2650361.1187479496, -2650361.118747951, 2505995.0066200756, 2505995.0066200765)},
        // THE ONE THAT KILLS MUTANT 9 ON ITS OWN. Same near-parallel family,
        // and here the computed u is 1.5 -- half a segment beyond s.b -- so
        // without the clamp the constructed point lands 3.3e7 m outside the
        // intersection of the two coordinate ranges, on a pair whose true
        // crossing is interior to both.
        //
        // PROVENANCE, and it corrects the increment file: 05-noder.md rules
        // that such a fixture "may not exist at reasonable search effort" and
        // that the universal range invariant kills mutant 9 without one.
        // Measured here, it does not -- with only the cases above the mutant
        // survives the whole suite, because a well-conditioned crossing lands
        // in range with or without the clamp. This pair came out of a 4e5-sample
        // scan over near-parallel crossings at Web Mercator magnitudes, seeded,
        // in under a minute.
        {seg(-20037508.0, -20037508.0, 20037508.0, 20037508.0),
         seg(3102738.5570855252, 3102738.557085525, 6767829.455649719, 6767829.45564972)},
    };

    for (const auto& pair : cases) {
        const Segment2& s = pair[0];
        const Segment2& t = pair[1];
        INFO("s = " << describe(s) << "\nt = " << describe(t));
        REQUIRE(rel(s, t) == SegmentRelation::Crossing);

        const GridPoint g = crossing_point<DefaultKernel>(grid, s, t);
        const Point2 lo = grid.cell_min(g);
        const Point2 hi = grid.cell_max(g);

        const double xlo = std::max(std::min(s.a.x, s.b.x), std::min(t.a.x, t.b.x));
        const double xhi = std::min(std::max(s.a.x, s.b.x), std::max(t.a.x, t.b.x));
        const double ylo = std::max(std::min(s.a.y, s.b.y), std::min(t.a.y, t.b.y));
        const double yhi = std::min(std::max(s.a.y, s.b.y), std::max(t.a.y, t.b.y));

        REQUIRE(lo.x <= xhi);
        REQUIRE(hi.x >= xlo);
        REQUIRE(lo.y <= yhi);
        REQUIRE(hi.y >= ylo);
    }
}

// ---------------------------------------------------------------------------
// segment_meets_cell: the hot-pixel predicate
// ---------------------------------------------------------------------------

// IF ONE TEST FROM THIS INCREMENT IS READ IN FIVE YEARS, IT SHOULD BE THIS ONE.
//
// A grid-collinear diagonal T-junction at 0.1 m at UTM33 magnitudes. The node G
// is the exact index-space midpoint of the host's two endpoints, so in index
// space it is on the host by construction. In world coordinates it is not:
// fl(ix * 0.1) perturbs each of the three points independently, and the node
// ends up 4.254e-10 m off the host's line -- 0.457 ulp at this northing, where
// ulp(7.9e6) = 9.313e-10, and 4.25e-9 of a cell. Recomputed with exact rationals
// over the actual doubles.
//
// on_segment<DefaultKernel> asks EXACT INCIDENCE and correctly answers false.
// segment_meets_cell asks HOT-PIXEL PROXIMITY -- does the segment pass through
// the cell -- and answers true. Snap rounding is defined on the second question.
// A design that detects T-junctions through the first reports Disjoint for most
// real ones, never splits the host, and then agrees with its own error when the
// verification pass asks the same question.
//
// classify is asserted here too, to pin what the finding costs: the exact
// classification of host against branch is NOT Touching. (Measured: Disjoint.
// Asserted as "not Touching", because Crossing is reachable for a neighbouring
// fixture and which one comes out is an ulp-level accident.) Mutant 16 is
// segment_meets_cell implemented as on_segment, and this test kills it alone.
TEST_CASE("a grid-collinear diagonal T-junction is missed by on_segment and found by the hot pixel",
          "[noding][intersect][hot_pixel]") {
    const SnapGrid grid{decimal_spacing};

    const GridPoint a{5000000, 79000000};
    const GridPoint g{a.ix + 137, a.iy + 61};    // the index-space midpoint
    const GridPoint b{a.ix + 274, a.iy + 122};   // ... of these two
    const GridPoint branch_end{g.ix - 61, g.iy + 137};

    const Segment2 host{grid.world(a), grid.world(b)};
    const Segment2 branch{grid.world(g), grid.world(branch_end)};

    REQUIRE_FALSE(on_segment<DefaultKernel>(host, grid.world(g)));
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, host, g));

    REQUIRE(rel(host, branch) != SegmentRelation::Touching);
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, branch, g));  // its own endpoint's cell
}

TEST_CASE("segment_meets_cell separates by the box faces and by the segment's line",
          "[noding][intersect][hot_pixel]") {
    const SnapGrid unit{1.0};
    const GridPoint g{0, 0};  // the closed cell [-0.5, 0.5] x [-0.5, 0.5]

    SECTION("crossing the cell's interior through no lattice point") {
        require_meets(unit, seg(-2.0, 0.2, 2.0, 0.3), g, true);
        require_meets(unit, seg(-0.4, -0.4, 0.4, 0.45), g, true);
    }

    SECTION("grazing exactly one corner") {
        // The line x + y = 1 passes through (0.5, 0.5) and leaves the other
        // three corners strictly on one side: one orientation is Collinear and
        // three are equal and non-Collinear. The condition is "all four the
        // same NON-Collinear value", so the answer is true.
        //
        // Written as a strict sign test -- all > 0 or all < 0 -- this becomes
        // false and the predicate loses exactly the incidences it exists to
        // find. Mutant 13, the most innocent-looking diff in 5a.
        require_meets(unit, seg(0.0, 1.0, 1.0, 0.0), g, true);
        require_meets(unit, seg(-1.0, 0.0, 0.0, -1.0), g, true);  // the opposite corner
    }

    SECTION("a far-away segment whose line passes through the cell") {
        // Mutant 14: drop the bbox half and this is true. Nobody writes this
        // fixture by accident, because a hand-picked fixture puts the cell near
        // the segment.
        require_meets(unit, seg(100.0, 100.0, 200.0, 200.0), g, false);
        require_meets(unit, seg(-500.0, -500.0, -100.0, -100.0), g, false);
    }

    SECTION("a long diagonal whose bbox contains the cell but whose line misses it") {
        // Mutant 15: drop the orientation half and this is true, because for a
        // long diagonal the bounding box is most of the domain.
        require_meets(unit, seg(-10.0, -8.0, 10.0, 12.0), g, false);
        require_meets(unit, seg(-10.0, 12.0, 10.0, -8.0), g, false);
    }

    SECTION("axis-aligned along a cell edge and one half-cell out") {
        require_meets(unit, seg(-5.0, 0.5, 5.0, 0.5), g, true);   // along the top edge
        require_meets(unit, seg(-5.0, -0.5, 5.0, -0.5), g, true);  // along the bottom edge
        require_meets(unit, seg(0.5, -5.0, 0.5, 5.0), g, true);    // along the right edge
        require_meets(unit, seg(-5.0, 1.5, 5.0, 1.5), g, false);   // one half-cell out
        require_meets(unit, seg(1.5, -5.0, 1.5, 5.0), g, false);
    }

    SECTION("zero-length segments") {
        require_meets(unit, seg(0.0, 0.0, 0.0, 0.0), g, true);      // the cell centre
        require_meets(unit, seg(0.4999, 0.4999, 0.4999, 0.4999), g, true);   // just inside a corner
        require_meets(unit, seg(0.5, 0.5, 0.5, 0.5), g, true);      // exactly on a corner
        require_meets(unit, seg(0.5001, 0.5, 0.5001, 0.5), g, false);  // just outside
        require_meets(unit, seg(3.0, 3.0, 3.0, 3.0), g, false);
    }
}

// Step 4 of 5b's split pass relies on this to not lose endpoints: a segment
// always meets the cells of its own endpoints, whatever the spacing and
// whatever the rounding did.
TEST_CASE("a segment meets the cells of its own endpoints", "[noding][intersect][hot_pixel]") {
    const SnapGrid grid{decimal_spacing};

    const Point2 raw[][2] = {
        {Point2{500000.04, 7900000.06}, Point2{500123.456, 7900987.654}},
        {Point2{-17.37, 42.4242}, Point2{0.05, -0.05}},
        {Point2{500000.0, 7900000.0}, Point2{500000.0, 7900000.0}},  // zero-length
    };

    for (const auto& pair : raw) {
        const GridPoint ga = grid.snap(pair[0]);
        const GridPoint gb = grid.snap(pair[1]);
        const Segment2 s{grid.world(ga), grid.world(gb)};
        INFO("s = " << describe(s));
        REQUIRE(segment_meets_cell<DefaultKernel>(grid, s, ga));
        REQUIRE(segment_meets_cell<DefaultKernel>(grid, s, gb));
    }
}

// ---------------------------------------------------------------------------
// The one template cross product in 5a
// ---------------------------------------------------------------------------

// THE ONLY TEMPLATE_TEST_CASE IN THIS INCREMENT, and it does one job: prove
// that the classification actually flows through K -- that nobody has written
// the four determinants inline in double. Same job and same justification as
// increment 3's single sliver-winding cross product.
//
// DEVIATION FROM THE DESIGN, recorded here rather than silently weakened.
// docs/increments/05-noder.md specifies this fixture as "a near-parallel
// crossing at UTM33 magnitudes where FastKernel mis-signs one orientation and
// reports Disjoint". At UTM33 magnitudes that pair does not exist, and the
// reason is arithmetic: FastKernel's error in orient2d is about eps * |dx * dy|,
// while the smallest nonzero determinant the coordinates can express is about
// |dx| * ulp(y) ~ eps * |dx| * |y|. The first exceeds the second only when a
// point's offset from the segment is comparable to the absolute coordinate
// magnitude -- at UTM33 northing 7.9e6 over even a 100 km domain that ratio is
// 0.013. Measured: over 1.8e6 near-line triples about a 100 km UTM33 segment,
// FastKernel and the exact orientation disagreed on ZERO of them.
//
// The fixture below is therefore at Web Mercator full-extent magnitudes
// (+/-2.0037e7 with a 4e7 span), where the ratio reaches 2 and FastKernel does
// go wrong. It does not go wrong in the direction the design predicted either:
// it answers Collinear rather than a flipped sign, which drops a genuine
// crossing into classify's ENDPOINT ARM -- exactly the case the design says
// that arm exists for. The arm then answers Touching without floating-point
// contraction and Disjoint with it (both measured, by emulating the two ways a
// compiler may fuse the determinant), so the assertion is the robust one they
// share: under FastKernel the answer is NOT Crossing, and under DefaultKernel
// it is.
TEMPLATE_TEST_CASE("a near-parallel crossing is decided by the kernel", "[noding][intersect]",
                   FastKernel, DefaultKernel) {
    const Segment2 s = seg(-20037508.0, -20037508.0, 20037508.0, 20037508.0);
    const Segment2 t = seg(-2650361.1187479496, -2650361.118747951,
                           2505995.0066200756, 2505995.0066200765);

    const SegmentRelation r = classify<TestType>(s, t);

    if constexpr (std::same_as<TestType, DefaultKernel>) {
        REQUIRE(r == SegmentRelation::Crossing);
    } else {
        REQUIRE(r != SegmentRelation::Crossing);
    }
}
