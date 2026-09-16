// The one and only translation unit in this repository that includes
// detria.hpp. See lib/detria/README.md for why that is a rule and how it is
// enforced.

#include <terrain/predicates/detria_exact.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/orientation.hpp>

#include <detria.hpp>

namespace terrain::pred {

namespace {

// `true` selects detria's adaptive, exactly-rounded evaluation (Shewchuk's
// expansion arithmetic); `false` selects the plain double determinant. This is
// an exact backend, so the choice is not a tuning knob: with `Robust = false`
// the backend would answer the same question the filter already answered in
// double arithmetic, and the fallback would be decorative. The suite kills the
// mutation -- the near-degenerate and exactly-collinear families disagree with
// RefExact immediately.
constexpr bool robust = true;

// detria's predicates are templated on any vector type exposing `.x` and `.y`,
// so `terrain::Point2` is passed straight through: no adapter struct, no copy,
// no conversion. `decltype(Vec2::x)` is `double`, so detria takes its
// floating-point path and no integer widening happens.
using Vec2 = Point2;

// The two enum translations. Both are explicit switches, and neither may become
// a `static_cast`.
//
// detria's values are `math::Orientation{CW=0, CCW=1, Collinear=2}` and
// `math::CircleLocation{Inside=0, Outside=1, Cocircular=2}`; ours are
// `Orientation{Clockwise=-1, Collinear=0, CounterClockwise=1}` and
// `Incircle{Outside=-1, Cocircular=0, Inside=1}`. The values differ *and* the
// orders differ, so a cast is not merely unportable, it is wrong in a way that
// still type-checks: it would report a clockwise triple as collinear, a point
// inside a circumcircle as on it, and a point outside as inside -- the last of
// which is the difference between a Delaunay flip loop terminating and not. The
// unit suite mutation-tests exactly these substitutions.
[[nodiscard]] Orientation from_detria(detria::math::Orientation o) noexcept {
    switch (o) {
        case detria::math::Orientation::CW:
            return Orientation::Clockwise;
        case detria::math::Orientation::CCW:
            return Orientation::CounterClockwise;
        case detria::math::Orientation::Collinear:
            return Orientation::Collinear;
    }
    // Unreachable: the enumeration is closed and detria returns nothing else.
    // Collinear is the answer that does not claim a turn, matching
    // `orientation_of_sign`'s treatment of an uninterpretable determinant.
    return Orientation::Collinear;
}

[[nodiscard]] Incircle from_detria(detria::math::CircleLocation l) noexcept {
    switch (l) {
        case detria::math::CircleLocation::Inside:
            return Incircle::Inside;
        case detria::math::CircleLocation::Outside:
            return Incircle::Outside;
        case detria::math::CircleLocation::Cocircular:
            return Incircle::Cocircular;
    }
    return Incircle::Cocircular;
}

}  // namespace

// A note for the next reader, so this is not re-derived:
//
// detria's `orient2d` and `incircle` carry Debug-only `detail::detriaAssert`
// calls on integer-overflow conditions (detria.hpp lines ~1150, ~1256, ~1275
// and ~1310 at the pinned SHA), and `detriaAssert` raises SIGTRAP. Reading
// `checkOverflowMultiply` it looks as though it could fire spuriously for
// negative `double` operands, because it compares against
// `std::numeric_limits<Scalar>::min()`, which for a floating-point type is the
// smallest *positive* normal rather than the most negative value. It does not
// fire: the whole body is guarded by `if constexpr
// (std::is_floating_point_v<Scalar>) return false;`, so for our double
// coordinates the checks compile to a constant `false` and the asserts are
// vacuous. @tester confirmed empirically as well -- several thousand calls with
// negative coordinates in a Debug build with assertions live, nothing raised.
//
// The `incircle` assertion at line 1237 is the live one, and it is real: it
// checks the counterclockwise precondition. `FilteredKernel::incircle`
// normalizes before ever reaching here.

Orientation DetriaExact::orient2d(const Point2& a, const Point2& b, const Point2& c) noexcept {
    return from_detria(detria::math::orient2d<robust, Vec2>(a, b, c));
}

Incircle DetriaExact::incircle_ccw(const Point2& a, const Point2& b, const Point2& c,
                                   const Point2& d) noexcept {
    return from_detria(detria::math::incircle<robust, Vec2>(a, b, c, d));
}

}  // namespace terrain::pred
