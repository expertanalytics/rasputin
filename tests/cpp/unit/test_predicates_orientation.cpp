// Unit tests for terrain/predicates/orientation.hpp: the classification types
// every predicate in the module returns.
//
// These are small but load-bearing. The enumerators' underlying values are part
// of the interface -- callers cast to int and use the result as a sign -- and
// `orientation_of_sign` is the single funnel through which every floating-point
// determinant becomes a classification, which makes its treatment of signed
// zero a correctness property of the whole module rather than a detail.

#include <catch2/catch_test_macros.hpp>

#include <terrain/predicates/orientation.hpp>

#include <limits>
#include <type_traits>

using terrain::pred::Incircle;
using terrain::pred::Orientation;
using terrain::pred::incircle_of_sign;
using terrain::pred::is_collinear;
using terrain::pred::is_left_turn;
using terrain::pred::orientation_of_sign;
using terrain::pred::reversed;

namespace {
constexpr double nan_value = std::numeric_limits<double>::quiet_NaN();
constexpr double positive_zero = 0.0;
constexpr double negative_zero = -0.0;
}  // namespace

TEST_CASE("Orientation enumerators carry their sign as the underlying value", "[predicates][orientation]") {
    STATIC_REQUIRE(std::is_same_v<std::underlying_type_t<Orientation>, int>);
    STATIC_REQUIRE(static_cast<int>(Orientation::Clockwise) == -1);
    STATIC_REQUIRE(static_cast<int>(Orientation::Collinear) == 0);
    STATIC_REQUIRE(static_cast<int>(Orientation::CounterClockwise) == 1);
}

TEST_CASE("Incircle enumerators carry their sign as the underlying value", "[predicates][orientation]") {
    STATIC_REQUIRE(std::is_same_v<std::underlying_type_t<Incircle>, int>);
    STATIC_REQUIRE(static_cast<int>(Incircle::Outside) == -1);
    STATIC_REQUIRE(static_cast<int>(Incircle::Cocircular) == 0);
    STATIC_REQUIRE(static_cast<int>(Incircle::Inside) == 1);
}

TEST_CASE("orientation_of_sign maps a determinant to a classification", "[predicates][orientation]") {
    STATIC_REQUIRE(orientation_of_sign(1.0) == Orientation::CounterClockwise);
    STATIC_REQUIRE(orientation_of_sign(-1.0) == Orientation::Clockwise);
    STATIC_REQUIRE(orientation_of_sign(0.0) == Orientation::Collinear);

    // Subnormal magnitudes are still definite signs, not noise. The decision
    // about whether a determinant this small is trustworthy belongs to the
    // filter, not to the classification function.
    STATIC_REQUIRE(orientation_of_sign(std::numeric_limits<double>::denorm_min()) ==
                   Orientation::CounterClockwise);
    STATIC_REQUIRE(orientation_of_sign(-std::numeric_limits<double>::denorm_min()) ==
                   Orientation::Clockwise);
}

// Signed zero is the hazard here. `-0.0 < 0.0` is false and `-0.0 == 0.0` is
// true, so a sign test written as `v < 0 ? ... : v > 0 ? ...` happens to be
// right while one written via `std::signbit` splits the collinear case in two
// and makes `orientation_of_sign(+0.0) != orientation_of_sign(-0.0)`. Two
// callers whose determinants differ only in the sign of a zero would then
// disagree about whether three points are collinear. point.hpp's hash
// specialisation documents the same hazard for coordinates.
TEST_CASE("orientation_of_sign does not split the collinear case on signed zero", "[predicates][orientation]") {
    STATIC_REQUIRE(orientation_of_sign(positive_zero) == Orientation::Collinear);
    STATIC_REQUIRE(orientation_of_sign(negative_zero) == Orientation::Collinear);
    STATIC_REQUIRE(orientation_of_sign(positive_zero) == orientation_of_sign(negative_zero));
}

TEST_CASE("incircle_of_sign does not split the cocircular case on signed zero", "[predicates][orientation]") {
    STATIC_REQUIRE(incircle_of_sign(1.0) == Incircle::Inside);
    STATIC_REQUIRE(incircle_of_sign(-1.0) == Incircle::Outside);
    STATIC_REQUIRE(incircle_of_sign(positive_zero) == Incircle::Cocircular);
    STATIC_REQUIRE(incircle_of_sign(negative_zero) == Incircle::Cocircular);
    STATIC_REQUIRE(incircle_of_sign(positive_zero) == incircle_of_sign(negative_zero));
}

// Finiteness is a precondition validated at the PSLG and Python boundaries, not
// here; a NaN determinant is unreachable if those hold. The classification is
// nonetheless total, because a function that is total is one fewer thing to
// reason about, and Collinear/Cocircular is the choice that keeps a NaN from
// being reported as a definite turn.
TEST_CASE("orientation_of_sign is total on NaN", "[predicates][orientation]") {
    STATIC_REQUIRE(orientation_of_sign(nan_value) == Orientation::Collinear);
}

TEST_CASE("reversed negates the orientation", "[predicates][orientation]") {
    STATIC_REQUIRE(reversed(Orientation::CounterClockwise) == Orientation::Clockwise);
    STATIC_REQUIRE(reversed(Orientation::Clockwise) == Orientation::CounterClockwise);
    STATIC_REQUIRE(reversed(Orientation::Collinear) == Orientation::Collinear);
}

TEST_CASE("reversed agrees with negating the underlying value", "[predicates][orientation]") {
    constexpr Orientation all[] = {
        Orientation::Clockwise,
        Orientation::Collinear,
        Orientation::CounterClockwise,
    };
    for (const Orientation o : all) {
        REQUIRE(reversed(o) == static_cast<Orientation>(-static_cast<int>(o)));
    }
}

TEST_CASE("reversed is an involution", "[predicates][orientation]") {
    STATIC_REQUIRE(reversed(reversed(Orientation::CounterClockwise)) == Orientation::CounterClockwise);
    STATIC_REQUIRE(reversed(reversed(Orientation::Clockwise)) == Orientation::Clockwise);
    STATIC_REQUIRE(reversed(reversed(Orientation::Collinear)) == Orientation::Collinear);
}

TEST_CASE("is_left_turn is true only for a counterclockwise orientation", "[predicates][orientation]") {
    STATIC_REQUIRE(is_left_turn(Orientation::CounterClockwise));
    STATIC_REQUIRE_FALSE(is_left_turn(Orientation::Clockwise));
    // Collinear is not a left turn: a degenerate triple must not be reported as
    // turning either way, or a convex-hull walk will keep collinear points on
    // the hull depending on which side the test is written from.
    STATIC_REQUIRE_FALSE(is_left_turn(Orientation::Collinear));
}

TEST_CASE("is_collinear is true only for the degenerate orientation", "[predicates][orientation]") {
    STATIC_REQUIRE(is_collinear(Orientation::Collinear));
    STATIC_REQUIRE_FALSE(is_collinear(Orientation::CounterClockwise));
    STATIC_REQUIRE_FALSE(is_collinear(Orientation::Clockwise));
}
