// Increment 34 (docs/increments/34-slope-tolerance.md, sections 3, 4.2, 7 and
// 9), test 2: terrain::raster::steepness, every node's class against the
// test's own Horn in long double (support/slope_oracle.hpp).
//
//   G5  class(n) is the smallest half degree at or above Horn's slope of n,
//       missing neighbours filled as section 3 says; at a class boundary,
//       within rounding, the higher class.
//
// Fixtures (section 9): V1 (65 x 65 at 10 m), V2 (33 x 41 at 10 m by 5 m:
// dx != dy, so a Horn that uses dx for both axes, M4, changes classes on the
// wall), P1 (a plane at 37.3 degrees facing 30 degrees off the row on V2's
// spacing, 21 x 21, one NoData node at its centre), a one-row grid (the named
// departure: a 40-degree plane facing 30 degrees reads 36.005 degrees), and
// the half-weight case of section 3 (NoData above and below a node).
// Invariant-critical (section 9). Starts threads (1 and 8 must agree).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#if __has_include(<terrain/raster/steepness.hpp>)

#include <terrain/raster/raster.hpp>
#include <terrain/raster/steepness.hpp>

#include "slope_oracle.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace slope_oracle;
using terrain::raster::kSteepnessClasses;
using terrain::raster::steepness;

namespace {

// Every node's class within test 2's acceptance; the number of nodes checked.
std::size_t check_classes(const Raster<float>& dem, unsigned threads) {
    const auto got = steepness(dem, threads);
    const auto want = horn(dem);
    REQUIRE(got.size() == dem.geometry().size());
    std::size_t bad = 0;
    for (std::size_t i = 0; i < got.size(); ++i) {
        const bool ok = std::isnan(want[i]) ? got[i] == terrain::raster::kNoDataClass : class_ok(got[i], want[i]);
        if (!ok && ++bad <= 5) {
            CAPTURE(i / dem.geometry().cols(), i % dem.geometry().cols(), static_cast<double>(want[i]),
                    static_cast<int>(class_of(want[i])), static_cast<int>(got[i]));
            CHECK(ok);
        }
    }
    CHECK(bad == 0);
    return got.size();
}

}  // namespace

TEST_CASE("34 test 2: the constants are section 4.2's", "[slope][steepness][test2]") {
    STATIC_REQUIRE(kSteepnessClasses == 181);
    STATIC_REQUIRE(terrain::raster::kNoDataClass == 255);
    STATIC_REQUIRE(terrain::raster::kNoDataClass == slope_oracle::kNoDataClass);
}

TEST_CASE("34 test 2: every class on V1 and V2 is the smallest half degree at or above Horn's slope",
          "[slope][steepness][test2]") {
    const bool second = GENERATE(false, true);
    const auto dem = second ? v2() : v1();
    CAPTURE(second ? "V2" : "V1");
    CHECK(check_classes(dem, 1) == dem.geometry().size());
}

TEST_CASE("34 test 2: V2's classes are not Horn with dx for both axes (M4)", "[slope][steepness][test2]") {
    // The fixture can show M4: with dx for dy, 733 of V2's 1 353 classes change.
    const auto dem = v2();
    const auto got = steepness(dem, 1);
    const auto m4 = classes(dem, true);
    std::size_t differ = 0;
    for (std::size_t i = 0; i < got.size(); ++i) differ += got[i] != m4[i] ? 1 : 0;
    CHECK(differ > 700);
}

TEST_CASE("34 test 2: classes are rounded up, never down (M5)", "[slope][steepness][test2]") {
    // On V1 nearly every slope lies strictly inside a class, so rounding down
    // would put it one class lower; never lower is the rule.
    const auto dem = v1();
    const auto got = steepness(dem, 1);
    const auto want = horn(dem);
    std::size_t inside = 0;
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (near_boundary(want[i]) || want[i] == 0.0L) continue;
        ++inside;
        CHECK(static_cast<long double>(got[i]) >= 2.0L * want[i]);
    }
    CHECK(inside > 4000);
}

TEST_CASE("34 test 2: P1, a plane at 37.3 degrees with a NoData node, is class 75 at every valid node",
          "[slope][steepness][test2]") {
    // The border, the four corners and the eight nodes around the NoData one
    // included (item 1b: 37.300 degrees at every valid node); round 1's rule
    // read 12.7 degrees at a corner of this plane.
    const std::array<CellIndex, 1> hole{CellIndex{10, 10}};
    const auto dem = plane(21, 21, 10.0, 5.0, 37.3, 30.0, hole);
    const auto got = steepness(dem, 1);
    REQUIRE(got.size() == 21u * 21u);
    for (std::size_t r = 0; r < 21; ++r)
        for (std::size_t c = 0; c < 21; ++c) {
            CAPTURE(r, c);
            CHECK(static_cast<int>(got[r * 21 + c]) == (r == 10 && c == 10 ? 255 : 75));
        }
}

TEST_CASE("34 test 2: a one-row grid loses the cross-row part of the slope (the named departure)",
          "[slope][steepness][test2]") {
    // A 40-degree plane facing 30 degrees off the row: gy = 0, so Horn reads
    // atan(tan 40 cos 30) = 36.005 degrees, class 73, at every node.
    const auto dem = plane(1, 21, 10.0, 10.0, 40.0, 30.0);
    const auto got = steepness(dem, 1);
    REQUIRE(got.size() == 21u);
    for (std::size_t c = 0; c < 21; ++c) {
        CAPTURE(c);
        CHECK(static_cast<int>(got[c]) == 73);
    }
    check_classes(dem, 1);
}

TEST_CASE("34 test 2: NoData above and below a node leaves that axis to the corners, at half weight",
          "[slope][steepness][test2]") {
    // Section 3: a 40-degree plane facing north reads 22.760 degrees there
    // (atan(tan 40 / 2)), class 46; its neighbours read the plane, class 80.
    const std::array<CellIndex, 2> holes{CellIndex{9, 10}, CellIndex{11, 10}};
    const auto dem = plane(21, 21, 10.0, 10.0, 40.0, 90.0, holes);
    const auto got = steepness(dem, 1);
    CHECK(static_cast<int>(got[10 * 21 + 10]) == 46);
    CHECK(static_cast<int>(got[9 * 21 + 10]) == 255);
    CHECK(static_cast<int>(got[11 * 21 + 10]) == 255);
    check_classes(dem, 1);
}

TEST_CASE("34 test 2: a slope exactly on a class boundary takes the higher class (pinned)",
          "[slope][steepness][test2]") {
    // z = x on a 10 m grid: gx = 1 exactly, so Horn's slope is 45 degrees,
    // the boundary between classes 90 and 91. "At a class boundary, within
    // rounding, the higher class": 91 at every node, the border included.
    const auto dem = plane(9, 9, 10.0, 10.0, 45.0, 0.0);
    for (std::size_t i = 0; i < 81; ++i) REQUIRE(dem.value_at({i / 9, i % 9}) == static_cast<float>(10 * (i % 9)));
    const auto got = steepness(dem, 1);
    for (std::size_t i = 0; i < got.size(); ++i) {
        CAPTURE(i);
        CHECK(static_cast<int>(got[i]) == 91);
    }
}

TEST_CASE("34 test 2: a flat DEM is class 0 everywhere, NoData 255", "[slope][steepness][test2]") {
    std::vector<float> z(7 * 5, 12.5f);
    z[2 * 5 + 2] = kNoData;
    const Raster<float> dem{geometry(7, 5, 10.0, 5.0), std::move(z), kNoData};
    const auto got = steepness(dem, 1);
    for (std::size_t i = 0; i < got.size(); ++i) {
        CAPTURE(i);
        CHECK(static_cast<int>(got[i]) == (i == 12 ? 255 : 0));
    }
}

TEST_CASE("34 test 2: V1n's NoData node is class 255 and its neighbours keep a plane-exact fill",
          "[slope][steepness][test2]") {
    const auto dem = v1n();
    const auto got = steepness(dem, 1);
    CHECK(static_cast<int>(got[kV1nHole.row * 65 + kV1nHole.col]) == 255);
    check_classes(dem, 1);
}

TEST_CASE("34 test 2: 1 and 8 threads give the same classes", "[slope][steepness][test2][threads]") {
    const int which = GENERATE(0, 1, 2);
    const auto dem = which == 0 ? v1() : which == 1 ? v2() : v1n();
    CAPTURE(which);
    CHECK(steepness(dem, 1) == steepness(dem, 8));
    // A grid tall enough to give each of eight threads several rows.
    const auto tall = valley(400, 65, 10.0, 10.0);
    CHECK(steepness(tall, 1) == steepness(tall, 8));
    check_classes(tall, 8);
}

#else

TEST_CASE("34: include/terrain/raster/steepness.hpp does not exist yet", "[slope][steepness]") {
    FAIL("increment 34 is not built: no terrain/raster/steepness.hpp, so terrain::raster::steepness, "
         "kSteepnessClasses and kNoDataClass do not exist and test 2 cannot compile");
}

#endif
