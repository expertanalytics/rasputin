#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>

#include <format>
#include <string>
#include <type_traits>

using terrain::Point2;

TEST_CASE("Point2 aggregate initializes and defaults to origin", "[point][point2]") {
    constexpr Point2 p{1.0, 2.0};
    STATIC_REQUIRE(p.x == 1.0);
    STATIC_REQUIRE(p.y == 2.0);

    constexpr Point2 origin{};
    STATIC_REQUIRE(origin.x == 0.0);
    STATIC_REQUIRE(origin.y == 0.0);
}

TEST_CASE("Point2 equality is coordinate-wise", "[point][point2]") {
    STATIC_REQUIRE(Point2{1.0, 2.0} == Point2{1.0, 2.0});
    STATIC_REQUIRE(Point2{1.0, 2.0} != Point2{1.0, 3.0});
    STATIC_REQUIRE(Point2{1.0, 2.0} != Point2{2.0, 2.0});
}

TEST_CASE("Point2 supports vector arithmetic", "[point][point2]") {
    STATIC_REQUIRE(Point2{1.0, 2.0} + Point2{3.0, 4.0} == Point2{4.0, 6.0});
    STATIC_REQUIRE(Point2{3.0, 4.0} - Point2{1.0, 2.0} == Point2{2.0, 2.0});
    STATIC_REQUIRE(2.0 * Point2{1.0, 2.0} == Point2{2.0, 4.0});
    STATIC_REQUIRE(Point2{1.0, 2.0} * 2.0 == Point2{2.0, 4.0});
    STATIC_REQUIRE(Point2{2.0, 4.0} / 2.0 == Point2{1.0, 2.0});
    STATIC_REQUIRE(-Point2{1.0, 2.0} == Point2{-1.0, -2.0});
}

TEST_CASE("Point2 dot and 2D cross behave as vector algebra requires", "[point][point2]") {
    STATIC_REQUIRE(dot(Point2{1.0, 0.0}, Point2{0.0, 1.0}) == 0.0);
    STATIC_REQUIRE(dot(Point2{1.0, 2.0}, Point2{3.0, 4.0}) == 11.0);
    STATIC_REQUIRE(cross(Point2{1.0, 0.0}, Point2{0.0, 1.0}) == 1.0);
    STATIC_REQUIRE(cross(Point2{0.0, 1.0}, Point2{1.0, 0.0}) == -1.0);
    STATIC_REQUIRE(cross(Point2{2.0, 3.0}, Point2{2.0, 3.0}) == 0.0);
}

TEST_CASE("Point2 formats via std::format", "[point][point2]") {
    REQUIRE(std::format("{}", Point2{1.0, 2.5}) == "Point2(1, 2.5)");
    REQUIRE(std::format("{}", Point2{}) == "Point2(0, 0)");
}

TEST_CASE("Point2 is a trivially copyable, standard-layout aggregate", "[point][traits]") {
    STATIC_REQUIRE(std::is_trivially_copyable_v<Point2>);
    STATIC_REQUIRE(std::is_standard_layout_v<Point2>);
    STATIC_REQUIRE(std::is_aggregate_v<Point2>);
    STATIC_REQUIRE(sizeof(Point2) == 2 * sizeof(double));
}

TEST_CASE("the Point2 formatter rejects a spec it does not honour", "[point][format]") {
    const Point2 p{1.0, 2.0};
    REQUIRE(std::format("{}", p) == "Point2(1, 2)");

    // A width spec used to parse successfully and then be silently dropped,
    // so "{:>24}" produced unpadded output with no diagnostic.
    REQUIRE_THROWS_AS(std::vformat("{:>24}", std::make_format_args(p)), std::format_error);
}
