#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>

#include <format>
#include <string>
#include <type_traits>
#include <unordered_set>

using terrain::Point2;
using terrain::Point3;

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

TEST_CASE("Point2 hash deduplicates identical coordinates in unordered_set", "[point][point2]") {
    std::unordered_set<Point2> s;
    s.insert({1.0, 2.0});
    s.insert({1.0, 2.0});
    s.insert({3.0, 4.0});
    REQUIRE(s.size() == 2);
    REQUIRE(s.contains(Point2{1.0, 2.0}));
    REQUIRE(s.contains(Point2{3.0, 4.0}));
    REQUIRE_FALSE(s.contains(Point2{2.0, 1.0}));
}

TEST_CASE("Point3 aggregate initializes and equality is coordinate-wise", "[point][point3]") {
    constexpr Point3 p{1.0, 2.0, 3.0};
    STATIC_REQUIRE(p.x == 1.0);
    STATIC_REQUIRE(p.y == 2.0);
    STATIC_REQUIRE(p.z == 3.0);

    STATIC_REQUIRE(Point3{1, 2, 3} == Point3{1, 2, 3});
    STATIC_REQUIRE(Point3{1, 2, 3} != Point3{1, 2, 4});
    STATIC_REQUIRE(Point3{1, 2, 3} != Point3{1, 3, 3});
    STATIC_REQUIRE(Point3{1, 2, 3} != Point3{2, 2, 3});
}

TEST_CASE("Point3 supports vector arithmetic", "[point][point3]") {
    STATIC_REQUIRE(Point3{1, 2, 3} + Point3{4, 5, 6} == Point3{5, 7, 9});
    STATIC_REQUIRE(Point3{4, 5, 6} - Point3{1, 2, 3} == Point3{3, 3, 3});
    STATIC_REQUIRE(2.0 * Point3{1, 2, 3} == Point3{2, 4, 6});
    STATIC_REQUIRE(Point3{2, 4, 6} / 2.0 == Point3{1, 2, 3});
    STATIC_REQUIRE(-Point3{1, 2, 3} == Point3{-1, -2, -3});
}

TEST_CASE("Point3 dot and 3D cross follow right-hand rule", "[point][point3]") {
    STATIC_REQUIRE(dot(Point3{1, 0, 0}, Point3{0, 1, 0}) == 0.0);
    STATIC_REQUIRE(dot(Point3{1, 2, 3}, Point3{4, 5, 6}) == 32.0);

    STATIC_REQUIRE(cross(Point3{1, 0, 0}, Point3{0, 1, 0}) == Point3{0, 0, 1});
    STATIC_REQUIRE(cross(Point3{0, 1, 0}, Point3{0, 0, 1}) == Point3{1, 0, 0});
    STATIC_REQUIRE(cross(Point3{0, 0, 1}, Point3{1, 0, 0}) == Point3{0, 1, 0});
    STATIC_REQUIRE(cross(Point3{1, 2, 3}, Point3{1, 2, 3}) == Point3{0, 0, 0});
}

TEST_CASE("Point3 formats via std::format", "[point][point3]") {
    REQUIRE(std::format("{}", Point3{1.0, 2.0, 3.5}) == "Point3(1, 2, 3.5)");
    REQUIRE(std::format("{}", Point3{}) == "Point3(0, 0, 0)");
}

TEST_CASE("Point3 hash deduplicates identical coordinates in unordered_set", "[point][point3]") {
    std::unordered_set<Point3> s;
    s.insert({1, 2, 3});
    s.insert({1, 2, 3});
    s.insert({4, 5, 6});
    REQUIRE(s.size() == 2);
    REQUIRE(s.contains(Point3{1, 2, 3}));
    REQUIRE(s.contains(Point3{4, 5, 6}));
    REQUIRE_FALSE(s.contains(Point3{3, 2, 1}));
}

TEST_CASE("Point types are trivially copyable, standard-layout aggregates", "[point][traits]") {
    STATIC_REQUIRE(std::is_trivially_copyable_v<Point2>);
    STATIC_REQUIRE(std::is_standard_layout_v<Point2>);
    STATIC_REQUIRE(std::is_aggregate_v<Point2>);
    STATIC_REQUIRE(sizeof(Point2) == 2 * sizeof(double));

    STATIC_REQUIRE(std::is_trivially_copyable_v<Point3>);
    STATIC_REQUIRE(std::is_standard_layout_v<Point3>);
    STATIC_REQUIRE(std::is_aggregate_v<Point3>);
    STATIC_REQUIRE(sizeof(Point3) == 3 * sizeof(double));
}

TEST_CASE("Point hashing is consistent with equality for signed zero", "[point][hash]") {
    // Both libc++ and libstdc++ normalize -0.0 in std::hash<double>. That is
    // required here, not incidental: operator== is defaulted, so +0.0 and -0.0
    // compare equal, and equal keys hashing differently would corrupt any
    // unordered container.
    const Point2 pos{0.0, 0.0};
    const Point2 neg{-0.0, 0.0};
    REQUIRE(pos == neg);
    REQUIRE(std::hash<Point2>{}(pos) == std::hash<Point2>{}(neg));
    REQUIRE(std::unordered_set<Point2>{pos, neg}.size() == 1);
}

TEST_CASE("Point formatters reject a spec they do not honour", "[point][format]") {
    const Point2 p{1.0, 2.0};
    const Point3 q{1.0, 2.0, 3.0};
    REQUIRE(std::format("{}", p) == "Point2(1, 2)");
    REQUIRE(std::format("{}", q) == "Point3(1, 2, 3)");

    // A width spec used to parse successfully and then be silently dropped,
    // so "{:>24}" produced unpadded output with no diagnostic.
    REQUIRE_THROWS_AS(std::vformat("{:>24}", std::make_format_args(p)), std::format_error);
    REQUIRE_THROWS_AS(std::vformat("{:>24}", std::make_format_args(q)), std::format_error);
}
