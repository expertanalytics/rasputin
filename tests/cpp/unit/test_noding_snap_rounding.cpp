// Unit tests for terrain/core/snap_grid.hpp -- increment 5a's arithmetic.
//
// INVARIANT-CRITICAL, mutation round (docs/increments/05-noder.md, "The
// invariant-critical suite (5a)"). The snapping contract is three lines of
// arithmetic and every defect in it is silent: a wrong rounding mode, a cached
// reciprocal or a cell corner computed in world space all produce plausible
// numbers that differ from the right ones in the last bit, and the noder's
// topology is decided by exactly those bits.
//
// Every fixture below whose value could not be derived by hand was found by a
// scan and the comment says which; nothing here is a number somebody thought
// looked adversarial.
//
// Deliberately NOT asserted anywhere in this file: that a snapped coordinate is
// exactly `ix * spacing` as a real number, or that a cell corner is exactly
// `(ix +/- 1/2) * spacing`. Commit e412a43 settled that they are not, at a
// non-dyadic spacing, and a test asserting it "approximately" would be a
// tolerance, which this project does not have. The dyadic fixtures are where
// exactness is asserted, because there it is true.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/core/snap_grid.hpp>

#include <cmath>
#include <cstdint>
#include <limits>

using terrain::GridPoint;
using terrain::Point2;
using terrain::SnapGrid;
using terrain::is_valid_spacing;
using terrain::kMaxGridIndex;

namespace {

// The two spacings the whole increment is argued about: a decimal one, where
// world() is not an affine map, and a dyadic one, where it is.
constexpr double decimal_spacing = 0.1;
constexpr double dyadic_spacing = 0.0625;  // 2^-4

}  // namespace

TEST_CASE("is_valid_spacing accepts exactly the finite positive spacings", "[noding][snap_grid]") {
    STATIC_REQUIRE(is_valid_spacing(0.1));
    STATIC_REQUIRE(is_valid_spacing(dyadic_spacing));
    STATIC_REQUIRE(is_valid_spacing(1.0));

    REQUIRE_FALSE(is_valid_spacing(0.0));

    // -0.0 does NOT separate the two spellings: `> 0.0` is false for it and
    // `!= 0.0` is also false, because `-0.0 == 0.0` in IEEE. Kept because the
    // answer is right and a reader will wonder. Mutant 5 is killed by exactly
    // three of the lines below -- -1.0, -0.1 and -inf -- and by no other. NaN
    // and +inf kill nothing: `<= max()` is false for both under either
    // spelling. Measured by building the mutant and running these arguments.
    REQUIRE_FALSE(is_valid_spacing(-0.0));

    REQUIRE_FALSE(is_valid_spacing(-1.0));
    REQUIRE_FALSE(is_valid_spacing(-0.1));
    REQUIRE_FALSE(is_valid_spacing(std::numeric_limits<double>::infinity()));
    REQUIRE_FALSE(is_valid_spacing(-std::numeric_limits<double>::infinity()));
    REQUIRE_FALSE(is_valid_spacing(std::numeric_limits<double>::quiet_NaN()));

    // Finite and strictly positive is the whole contract, so a subnormal is
    // admitted. Whether a grid that fine is a sane thing to ask for is the
    // driver's policy question and is answered in 5b, not here.
    REQUIRE(is_valid_spacing(std::numeric_limits<double>::denorm_min()));
    REQUIRE(is_valid_spacing(std::numeric_limits<double>::min()));
}

TEST_CASE("SnapGrid reports the spacing it was built with", "[noding][snap_grid]") {
    STATIC_REQUIRE(SnapGrid{decimal_spacing}.spacing() == decimal_spacing);
    STATIC_REQUIRE(SnapGrid{dyadic_spacing}.spacing() == dyadic_spacing);
}

// Ruling 1. llround rounds half away from zero unconditionally; rint and
// nearbyint honour the dynamic rounding mode and round half to even under the
// default one, so they disagree with llround at exactly the ties -- and only
// there. This is the fixture that separates them, and it is the only one that
// does. Mutant 1.
TEST_CASE("snap rounds ties away from zero on both sides of the origin", "[noding][snap_grid]") {
    const SnapGrid unit{1.0};

    // Under round-half-to-even these would be 0, 2, 2 and their negatives.
    REQUIRE(unit.snap(Point2{0.5, -0.5}) == GridPoint{1, -1});
    REQUIRE(unit.snap(Point2{1.5, -1.5}) == GridPoint{2, -2});
    REQUIRE(unit.snap(Point2{2.5, -2.5}) == GridPoint{3, -3});

    // The same tie at the decimal spacing: 0.05 / 0.1 is exactly 0.5 in double
    // (checked, not assumed -- it is asserted on the next line), so the tie
    // survives the division and the rounding mode is observable here too.
    REQUIRE(0.05 / decimal_spacing == 0.5);
    const SnapGrid grid{decimal_spacing};
    REQUIRE(grid.snap(Point2{0.05, -0.05}) == GridPoint{1, -1});
}

// Mutant 2: static_cast<std::int64_t>(x / s) truncates towards zero, which
// agrees with llround on everything below and disagrees on all of this.
TEST_CASE("snap rounds rather than truncates", "[noding][snap_grid]") {
    const SnapGrid unit{1.0};

    REQUIRE(unit.snap(Point2{0.7, -0.7}) == GridPoint{1, -1});
    REQUIRE(unit.snap(Point2{0.999, -0.999}) == GridPoint{1, -1});
    REQUIRE(unit.snap(Point2{-3.6, 3.6}) == GridPoint{-4, 4});
    REQUIRE(unit.snap(Point2{0.4, -0.4}) == GridPoint{0, 0});
}

// Ruling 2, mutant 3. `p.x * (1.0 / spacing)` is two roundings where
// `p.x / spacing` is one, and the two disagree on a set of inputs that is easy
// to find by search and impossible to predict by reading.
//
// PROVENANCE, and the comment is part of the test: this coordinate was not
// hand-picked. It came out of a scan over the doubles adjacent to half-index
// coordinates (k + 1/2) * 0.1 in the UTM33 easting band -- the neighbourhood of
// a tie is the only place the two spellings can disagree, because they differ
// by at most one ulp of the quotient. At 500865.14999999997 the divide gives a
// quotient just below the tie and the reciprocal gives one just above it.
TEST_CASE("snap divides by the spacing rather than multiplying by its reciprocal", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};
    const double x = 500865.14999999997;

    REQUIRE(grid.snap(Point2{x, 0.0}).ix == 5008651);

    // Stated as the mutant rather than as a second expected value, so a reader
    // sees what the line above is defending against.
    REQUIRE(std::llround(x * (1.0 / decimal_spacing)) == 5008652);
    REQUIRE(std::llround(x / decimal_spacing) == 5008651);
}

// Ruling 3. The signature is the enforcement -- there is no argument through
// which a bounding box or a raster origin could arrive -- so this property is
// weak by construction. It is asserted anyway because mutant 12 is live in 5b's
// driver and a reader has to be able to find the statement of the rule
// somewhere.
TEST_CASE("snap is a pure function of the point and the spacing", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};
    const Point2 p{500123.456, 7900987.654};
    const Point2 far{-999999.5, 12345678.25};

    const GridPoint alone = grid.snap(p);

    // The same point, snapped with wildly different traffic either side of it.
    (void)grid.snap(far);
    (void)grid.snap(Point2{0.0, 0.0});
    REQUIRE(grid.snap(p) == alone);

    // And through a second, independently constructed grid of the same spacing.
    const SnapGrid other{decimal_spacing};
    REQUIRE(other.snap(p) == alone);
}

TEST_CASE("snapped is world of snap", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};
    const Point2 points[] = {
        Point2{0.0, 0.0},
        Point2{500000.04, 7900000.06},
        Point2{-17.37, 42.4242},
        Point2{0.05, -0.05},
    };

    for (const Point2& p : points) {
        INFO("p = (" << p.x << ", " << p.y << ")");
        const Point2 w = grid.world(grid.snap(p));
        const Point2 s = grid.snapped(p);
        REQUIRE(s.x == w.x);  // bitwise: both are pure functions of the index
        REQUIRE(s.y == w.y);
    }
}

// can_snap is the driver's one-pass admission check and the only place
// |p| / spacing is compared against kMaxGridIndex. Mutant 4 is this bound raised
// or the comparison written strictly, so the boundary is pinned on each axis
// independently and in all four quadrants, in the style increment 3 used for
// sizes_fit_u32.
TEST_CASE("can_snap admits exactly the indices world stays injective on", "[noding][snap_grid]") {
    const SnapGrid unit{1.0};  // spacing 1.0: the bound is a coordinate, exactly
    const double limit = static_cast<double>(kMaxGridIndex);
    const double just_over = std::nextafter(limit, std::numeric_limits<double>::infinity());
    const double just_under = std::nextafter(limit, 0.0);

    SECTION("at, just under and just over the bound, per axis and per sign") {
        for (const double sx : {1.0, -1.0}) {
            for (const double sy : {1.0, -1.0}) {
                INFO("quadrant (" << sx << ", " << sy << ")");
                REQUIRE(unit.can_snap(Point2{sx * limit, sy * limit}));
                REQUIRE(unit.can_snap(Point2{sx * just_under, sy * just_under}));

                // One axis over at a time: the check is per axis, and an
                // implementation testing only x passes the y line below only by
                // accident.
                REQUIRE_FALSE(unit.can_snap(Point2{sx * just_over, sy * limit}));
                REQUIRE_FALSE(unit.can_snap(Point2{sx * limit, sy * just_over}));
                REQUIRE_FALSE(unit.can_snap(Point2{sx * just_over, sy * just_over}));
            }
        }
    }

    SECTION("the bound scales with the spacing") {
        const SnapGrid grid{decimal_spacing};
        REQUIRE(grid.can_snap(Point2{0.0, 0.0}));
        REQUIRE(grid.can_snap(Point2{500000.0, 7900000.0}));
        REQUIRE_FALSE(grid.can_snap(Point2{limit, 0.0}));       // limit/0.1 = 10x the bound
        REQUIRE_FALSE(grid.can_snap(Point2{0.0, -limit}));
        REQUIRE(grid.can_snap(Point2{limit * decimal_spacing, 0.0}));
    }
}

// Injectivity, and where kMaxGridIndex comes from. The first half says the
// bound is safe; the second says it is a number with a reason rather than a
// number, by exhibiting the collision it holds a factor of two of margin
// against. The over-large points are constructed directly -- can_snap would
// reject them and snap must never see them.
TEST_CASE("world is injective below kMaxGridIndex and collides far above it", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};
    const SnapGrid unit{1.0};

    const GridPoint at{kMaxGridIndex, kMaxGridIndex};
    const GridPoint next{kMaxGridIndex + 1, kMaxGridIndex + 1};

    REQUIRE(grid.world(at).x != grid.world(next).x);
    REQUIRE(grid.world(at).y != grid.world(next).y);
    REQUIRE(unit.world(at).x != unit.world(next).x);

    // 2^53: two adjacent grid points, two identical doubles, and dedup would
    // then leave two nodes at one coordinate. kMaxGridIndex is 2^51.
    constexpr std::int64_t beyond = std::int64_t{1} << 53;
    const GridPoint far{beyond, beyond};
    const GridPoint far_next{beyond + 1, beyond + 1};

    REQUIRE(grid.world(far).x == grid.world(far_next).x);
    REQUIRE(unit.world(far).x == unit.world(far_next).x);
    REQUIRE(far != far_next);  // distinct keys, identical coordinates: the failure
}

TEST_CASE("snap and world round-trip on the index range", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};

    const GridPoint cases[] = {
        GridPoint{0, 0},
        GridPoint{1, -1},
        GridPoint{5000000, 79000000},
        GridPoint{-5000000, -79000000},
        GridPoint{123456789, -987654321},
    };

    for (const GridPoint& g : cases) {
        INFO("g = (" << g.ix << ", " << g.iy << ")");
        REQUIRE(grid.snap(grid.world(g)) == g);
    }
}

// Ruling 5, and mutant 17. Cell corners are one rounding computed in index
// space: (2*ix +/- 1) * (spacing * 0.5), never world(g) +/- spacing * 0.5.
TEST_CASE("cell corners bracket the lattice point and abut the neighbouring cell", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};

    const GridPoint g{5000000, 79000000};
    const Point2 lo = grid.cell_min(g);
    const Point2 hi = grid.cell_max(g);
    const Point2 w = grid.world(g);

    REQUIRE(lo.x < w.x);
    REQUIRE(w.x < hi.x);
    REQUIRE(lo.y < w.y);
    REQUIRE(w.y < hi.y);

    // The cells partition the plane rather than being a set of squares with
    // gaps between them, and the statement of that is bitwise equality.
    const GridPoint east{g.ix + 1, g.iy};
    const GridPoint north{g.ix, g.iy + 1};
    REQUIRE(grid.cell_max(g).x == grid.cell_min(east).x);
    REQUIRE(grid.cell_max(g).y == grid.cell_min(north).y);
}

// PROVENANCE: found by scanning grid indices across the UTM33 easting band at
// 0.1 m and comparing the two spellings bit for bit -- not hand-picked, and not
// findable by reading. At index 4999003 the index-space corner is the double
// 499900.25 and the world-space spelling is 499900.25000000006, one step away.
//
// Note this mutant is behaviourally INVISIBLE at a dyadic spacing, which is why
// it has to be hunted at 0.1 m. Mutant 17.
TEST_CASE("cell corners are computed in index space, not in world space", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};
    const double half = decimal_spacing * 0.5;

    const GridPoint g{4999003, 6};
    REQUIRE(grid.cell_min(g).x == 499900.25);
    REQUIRE(grid.world(g).x - half == 499900.25000000006);
    REQUIRE(grid.cell_min(g).x != grid.world(g).x - half);

    REQUIRE(grid.cell_max(g).y == 0.65);
    REQUIRE(grid.world(g).y + half == 0.6500000000000001);
    REQUIRE(grid.cell_max(g).y != grid.world(g).y + half);
}

// The payoff "Choosing the spacing" claims, asserted rather than described: at
// a dyadic spacing the corner is an exact product of an exact integer and an
// exactly halved spacing, so it is the real number and not an approximation to
// it. 8000000 * 2^-4 = 500000 exactly; (2*8000000 - 1) * 2^-5 = 499999.96875
// exactly, and every digit of that literal is significant.
TEST_CASE("at a dyadic spacing the cell corners are exact", "[noding][snap_grid]") {
    const SnapGrid grid{dyadic_spacing};
    const GridPoint g{8000000, 126400000};

    REQUIRE(grid.world(g).x == 500000.0);
    REQUIRE(grid.cell_min(g).x == 499999.96875);
    REQUIRE(grid.cell_max(g).x == 500000.03125);

    REQUIRE(grid.world(g).y == 7900000.0);
    REQUIRE(grid.cell_min(g).y == 7899999.96875);
    REQUIRE(grid.cell_max(g).y == 7900000.03125);
}

// The one property that matters, stated as equality rather than as accuracy:
// two points snap to the same GridPoint if and only if their snapped world
// coordinates are bit-identical doubles. Dedup is legal because of this line.
TEST_CASE("coincidence has an exact key", "[noding][snap_grid]") {
    const SnapGrid grid{decimal_spacing};

    // Two distinct doubles inside one cell.
    const Point2 a{500000.04, 7900000.045};
    const Point2 b{500000.0400000001, 7900000.035};
    REQUIRE(a != b);
    REQUIRE(grid.snap(a) == grid.snap(b));
    REQUIRE(grid.snapped(a).x == grid.snapped(b).x);
    REQUIRE(grid.snapped(a).y == grid.snapped(b).y);

    // PROVENANCE: the adjacent-double pair below straddles the boundary between
    // cells 5000000 and 5000001 at 0.1 m; found by walking nextafter across the
    // half-index coordinate, not by hand. One ulp apart, two cells, and they
    // must stay two nodes.
    const Point2 under{500000.04999999993, 0.0};
    const Point2 over{500000.05, 0.0};
    REQUIRE(std::nextafter(under.x, std::numeric_limits<double>::infinity()) == over.x);
    REQUIRE(grid.snap(under).ix == 5000000);
    REQUIRE(grid.snap(over).ix == 5000001);
    REQUIRE(grid.snapped(under).x != grid.snapped(over).x);
}
