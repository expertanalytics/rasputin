#pragma once

// The vocabulary types of the predicates module: the classifications every
// predicate returns, and the total functions that turn a determinant's sign
// into one.
//
// This header is deliberately at the bottom of the module's dependency order.
// It knows nothing about points and nothing about arithmetic policy -- that is
// why `is_left_turn` takes an `Orientation` rather than three points. Anything
// that needs a coordinate type or a kernel belongs in exact.hpp or kernel.hpp.

namespace terrain::pred {

// The underlying values are part of the interface: `static_cast<int>` on an
// enumerator yields the sign of the determinant that produced it, so callers
// may compare, multiply or negate signs without a lookup table.
enum class Orientation : int {
    Clockwise = -1,
    Collinear = 0,
    CounterClockwise = 1,
};

enum class Incircle : int {
    Outside = -1,
    Cocircular = 0,
    Inside = 1,
};

// The single funnel through which a floating-point determinant becomes a
// classification.
//
// Signed zero is the hazard. `-0.0 == 0.0` is true and `-0.0 < 0.0` is false,
// so the comparison form below maps both zeros to Collinear, while a
// `std::signbit`-based test would split the collinear case in two and let two
// callers whose determinants differ only in the sign of a zero disagree about
// whether three points are collinear. point.hpp's hash specialisation carries
// the same warning for coordinates.
//
// NaN maps to Collinear because both comparisons are false. This is totality,
// not support: finiteness is a precondition checked at the PSLG and Python
// boundaries, so a NaN determinant is an unreachable state. Collinear is the
// only choice that does not report unreachable input as a definite turn.
[[nodiscard]] constexpr Orientation orientation_of_sign(double determinant) noexcept {
    if (determinant > 0.0) {
        return Orientation::CounterClockwise;
    }
    if (determinant < 0.0) {
        return Orientation::Clockwise;
    }
    return Orientation::Collinear;
}

// Same contract as `orientation_of_sign`, including the treatment of signed
// zero and NaN. Positive means the fourth point lies inside the circle through
// the first three, given that those three are counterclockwise.
[[nodiscard]] constexpr Incircle incircle_of_sign(double determinant) noexcept {
    if (determinant > 0.0) {
        return Incircle::Inside;
    }
    if (determinant < 0.0) {
        return Incircle::Outside;
    }
    return Incircle::Cocircular;
}

// The orientation of the same three points visited in the opposite order.
[[nodiscard]] constexpr Orientation reversed(Orientation o) noexcept {
    return static_cast<Orientation>(-static_cast<int>(o));
}

// There is deliberately no `reversed(Incircle)` overload. Reversing the vertex
// order of a triangle does not change the circle through its vertices, so the
// only thing such an overload could be used for is negating the result of a
// reordered `incircle` call -- which is precisely the bug documented in
// kernel.hpp, where an inverted answer is returned for clockwise input. A
// caller with a genuine need to negate an `Incircle` can cast, conspicuously,
// where a reviewer will see it.

// Collinear is not a left turn. A degenerate triple must not be reported as
// turning either way, or a hull walk keeps or drops collinear points depending
// on which side the test happens to be written from.
[[nodiscard]] constexpr bool is_left_turn(Orientation o) noexcept {
    return o == Orientation::CounterClockwise;
}

[[nodiscard]] constexpr bool is_collinear(Orientation o) noexcept {
    return o == Orientation::Collinear;
}

}  // namespace terrain::pred
