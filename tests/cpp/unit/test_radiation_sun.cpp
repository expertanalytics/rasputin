// Increment 46 PR 2, the terrain-radiation design's TR1: time and the sun's
// position (docs/increments/46-hbv-snow.md, sections 4.4 and 7, "PR 2, TR1";
// docs/increments/terrain-radiation.md, R2 and "Verification plan", item 1).
//
//   J1-J3  Julian day from UTC milliseconds and back, exact to the millisecond
//          (I2), before and after the Unix epoch;
//   N1-N3  NREL's published test case (Reda and Andreas 2004: 2003-10-17
//          12:30:30 UTC-7, Golden CO) through the split and the reference
//          path (I1); the split agrees with the reference path bit for bit;
//   S1     the term tables equal SPA's, term by term (I1);
//   D1-D2  the default delta T (Espenak and Meeus's polynomials, the /100
//          parabola outside 1986-2050) and the legacy defect it replaces;
//   A1-A3  the refraction guard: eq. 42 at and above -1.8 degrees, held below;
//   T1     sunrise, transit and sunset of the NREL case on the reference path;
//   G1     ranges over seeded random inputs.
//
// The expected term tables below are pvlib 0.16.1's (pvlib/spa.py, BSD-3), an
// independent transcription of the paper's tables; the legacy header
// (legacy/rasputin/solar_position.h at the legacy-archive tag) agrees with
// them on all 195 L/B/R rows and all 63 nutation rows.
//
// Every angle is in degrees, every time UTC. Absolute bounds: angles of at
// most 360 degrees, times between the years -2000 and 6000.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain/radiation/solar_position.hpp>
#include <terrain/radiation/time.hpp>

#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <span>
#include <utility>

namespace {

using namespace terrain::radiation;
using std::chrono::milliseconds;
using Catch::Matchers::WithinAbs;

// UTC milliseconds since the Unix epoch.
utc_ms at(std::int64_t ms) { return utc_ms{milliseconds{ms}}; }

// NREL's test case. 12:30:30 at UTC-7 is 19:30:30 UTC.
constexpr std::int64_t kNrelMs = 1066419030000;
constexpr std::int64_t kNrelDay = 1066348800000;  // 2003-10-17T00:00Z
constexpr double kLat = 39.742476;
constexpr double kLon = -105.1786;
constexpr double kHeight = 1830.14;
constexpr double kDeltaT = 67.0;
constexpr double kPressure = 820.0;  // hPa
constexpr double kTemperature = 11.0;  // degrees C

// The smallest angle between two azimuths, in degrees.
double azimuth_gap(double a, double b) {
    double d = std::fmod(std::abs(a - b), 360.0);
    return d > 180.0 ? 360.0 - d : d;
}

// Eq. 42 itself (Bennett's formula scaled by pressure and temperature), the
// refraction in degrees at true elevation e0, with no guard.
double bennett(double e0, double p_hpa, double t_c) {
    constexpr double deg = 3.14159265358979323846 / 180.0;
    return (p_hpa / 1010.0) * (283.0 / (273.0 + t_c)) * 1.02 /
           (60.0 * std::tan((e0 + 10.3 / (e0 + 5.11)) * deg));
}

double refraction(double e0, double p_hpa, double t_c) {
    return apparent_elevation(e0, p_hpa, t_c) - e0;
}

// The millisecond in [lo, hi] where f changes sign, f(lo) and f(hi) of
// opposite signs: bisection on whole milliseconds.
template <class F>
std::int64_t crossing(F f, std::int64_t lo, std::int64_t hi) {
    const bool lo_positive = f(lo) > 0.0;
    REQUIRE(lo_positive != (f(hi) > 0.0));
    while (hi - lo > 1) {
        const std::int64_t mid = lo + (hi - lo) / 2;
        if ((f(mid) > 0.0) == lo_positive)
            lo = mid;
        else
            hi = mid;
    }
    return lo;
}

struct Term {
    double a, b, c;
};
struct Nutation {
    std::array<int, 5> y;
    double a, b, c, d;
};

constexpr std::array<Term, 64> kL0{{
    {175347046.0, 0.0, 0.0},
    {3341656.0, 4.6692568, 6283.07585},
    {34894.0, 4.6261, 12566.1517},
    {3497.0, 2.7441, 5753.3849},
    {3418.0, 2.8289, 3.5231},
    {3136.0, 3.6277, 77713.7715},
    {2676.0, 4.4181, 7860.4194},
    {2343.0, 6.1352, 3930.2097},
    {1324.0, 0.7425, 11506.7698},
    {1273.0, 2.0371, 529.691},
    {1199.0, 1.1096, 1577.3435},
    {990.0, 5.233, 5884.927},
    {902.0, 2.045, 26.298},
    {857.0, 3.508, 398.149},
    {780.0, 1.179, 5223.694},
    {753.0, 2.533, 5507.553},
    {505.0, 4.583, 18849.228},
    {492.0, 4.205, 775.523},
    {357.0, 2.92, 0.067},
    {317.0, 5.849, 11790.629},
    {284.0, 1.899, 796.298},
    {271.0, 0.315, 10977.079},
    {243.0, 0.345, 5486.778},
    {206.0, 4.806, 2544.314},
    {205.0, 1.869, 5573.143},
    {202.0, 2.458, 6069.777},
    {156.0, 0.833, 213.299},
    {132.0, 3.411, 2942.463},
    {126.0, 1.083, 20.775},
    {115.0, 0.645, 0.98},
    {103.0, 0.636, 4694.003},
    {102.0, 0.976, 15720.839},
    {102.0, 4.267, 7.114},
    {99.0, 6.21, 2146.17},
    {98.0, 0.68, 155.42},
    {86.0, 5.98, 161000.69},
    {85.0, 1.3, 6275.96},
    {85.0, 3.67, 71430.7},
    {80.0, 1.81, 17260.15},
    {79.0, 3.04, 12036.46},
    {75.0, 1.76, 5088.63},
    {74.0, 3.5, 3154.69},
    {74.0, 4.68, 801.82},
    {70.0, 0.83, 9437.76},
    {62.0, 3.98, 8827.39},
    {61.0, 1.82, 7084.9},
    {57.0, 2.78, 6286.6},
    {56.0, 4.39, 14143.5},
    {56.0, 3.47, 6279.55},
    {52.0, 0.19, 12139.55},
    {52.0, 1.33, 1748.02},
    {51.0, 0.28, 5856.48},
    {49.0, 0.49, 1194.45},
    {41.0, 5.37, 8429.24},
    {41.0, 2.4, 19651.05},
    {39.0, 6.17, 10447.39},
    {37.0, 6.04, 10213.29},
    {37.0, 2.57, 1059.38},
    {36.0, 1.71, 2352.87},
    {36.0, 1.78, 6812.77},
    {33.0, 0.59, 17789.85},
    {30.0, 0.44, 83996.85},
    {30.0, 2.74, 1349.87},
    {25.0, 3.16, 4690.48},
}};
constexpr std::array<Term, 34> kL1{{
    {628331966747.0, 0.0, 0.0},
    {206059.0, 2.678235, 6283.07585},
    {4303.0, 2.6351, 12566.1517},
    {425.0, 1.59, 3.523},
    {119.0, 5.796, 26.298},
    {109.0, 2.966, 1577.344},
    {93.0, 2.59, 18849.23},
    {72.0, 1.14, 529.69},
    {68.0, 1.87, 398.15},
    {67.0, 4.41, 5507.55},
    {59.0, 2.89, 5223.69},
    {56.0, 2.17, 155.42},
    {45.0, 0.4, 796.3},
    {36.0, 0.47, 775.52},
    {29.0, 2.65, 7.11},
    {21.0, 5.34, 0.98},
    {19.0, 1.85, 5486.78},
    {19.0, 4.97, 213.3},
    {17.0, 2.99, 6275.96},
    {16.0, 0.03, 2544.31},
    {16.0, 1.43, 2146.17},
    {15.0, 1.21, 10977.08},
    {12.0, 2.83, 1748.02},
    {12.0, 3.26, 5088.63},
    {12.0, 5.27, 1194.45},
    {12.0, 2.08, 4694.0},
    {11.0, 0.77, 553.57},
    {10.0, 1.3, 6286.6},
    {10.0, 4.24, 1349.87},
    {9.0, 2.7, 242.73},
    {9.0, 5.64, 951.72},
    {8.0, 5.3, 2352.87},
    {6.0, 2.65, 9437.76},
    {6.0, 4.67, 4690.48},
}};
constexpr std::array<Term, 20> kL2{{
    {52919.0, 0.0, 0.0},
    {8720.0, 1.0721, 6283.0758},
    {309.0, 0.867, 12566.152},
    {27.0, 0.05, 3.52},
    {16.0, 5.19, 26.3},
    {16.0, 3.68, 155.42},
    {10.0, 0.76, 18849.23},
    {9.0, 2.06, 77713.77},
    {7.0, 0.83, 775.52},
    {5.0, 4.66, 1577.34},
    {4.0, 1.03, 7.11},
    {4.0, 3.44, 5573.14},
    {3.0, 5.14, 796.3},
    {3.0, 6.05, 5507.55},
    {3.0, 1.19, 242.73},
    {3.0, 6.12, 529.69},
    {3.0, 0.31, 398.15},
    {3.0, 2.28, 553.57},
    {2.0, 4.38, 5223.69},
    {2.0, 3.75, 0.98},
}};
constexpr std::array<Term, 7> kL3{{
    {289.0, 5.844, 6283.076},
    {35.0, 0.0, 0.0},
    {17.0, 5.49, 12566.15},
    {3.0, 5.2, 155.42},
    {1.0, 4.72, 3.52},
    {1.0, 5.3, 18849.23},
    {1.0, 5.97, 242.73},
}};
constexpr std::array<Term, 3> kL4{{
    {114.0, 3.142, 0.0},
    {8.0, 4.13, 6283.08},
    {1.0, 3.84, 12566.15},
}};
constexpr std::array<Term, 1> kL5{{
    {1.0, 3.14, 0.0},
}};
constexpr std::array<Term, 5> kB0{{
    {280.0, 3.199, 84334.662},
    {102.0, 5.422, 5507.553},
    {80.0, 3.88, 5223.69},
    {44.0, 3.7, 2352.87},
    {32.0, 4.0, 1577.34},
}};
constexpr std::array<Term, 2> kB1{{
    {9.0, 3.9, 5507.55},
    {6.0, 1.73, 5223.69},
}};
constexpr std::array<Term, 40> kR0{{
    {100013989.0, 0.0, 0.0},
    {1670700.0, 3.0984635, 6283.07585},
    {13956.0, 3.05525, 12566.1517},
    {3084.0, 5.1985, 77713.7715},
    {1628.0, 1.1739, 5753.3849},
    {1576.0, 2.8469, 7860.4194},
    {925.0, 5.453, 11506.77},
    {542.0, 4.564, 3930.21},
    {472.0, 3.661, 5884.927},
    {346.0, 0.964, 5507.553},
    {329.0, 5.9, 5223.694},
    {307.0, 0.299, 5573.143},
    {243.0, 4.273, 11790.629},
    {212.0, 5.847, 1577.344},
    {186.0, 5.022, 10977.079},
    {175.0, 3.012, 18849.228},
    {110.0, 5.055, 5486.778},
    {98.0, 0.89, 6069.78},
    {86.0, 5.69, 15720.84},
    {86.0, 1.27, 161000.69},
    {65.0, 0.27, 17260.15},
    {63.0, 0.92, 529.69},
    {57.0, 2.01, 83996.85},
    {56.0, 5.24, 71430.7},
    {49.0, 3.25, 2544.31},
    {47.0, 2.58, 775.52},
    {45.0, 5.54, 9437.76},
    {43.0, 6.01, 6275.96},
    {39.0, 5.36, 4694.0},
    {38.0, 2.39, 8827.39},
    {37.0, 0.83, 19651.05},
    {37.0, 4.9, 12139.55},
    {36.0, 1.67, 12036.46},
    {35.0, 1.84, 2942.46},
    {33.0, 0.24, 7084.9},
    {32.0, 0.18, 5088.63},
    {32.0, 1.78, 398.15},
    {28.0, 1.21, 6286.6},
    {28.0, 1.9, 6279.55},
    {26.0, 4.59, 10447.39},
}};
constexpr std::array<Term, 10> kR1{{
    {103019.0, 1.10749, 6283.07585},
    {1721.0, 1.0644, 12566.1517},
    {702.0, 3.142, 0.0},
    {32.0, 1.02, 18849.23},
    {31.0, 2.84, 5507.55},
    {25.0, 1.32, 5223.69},
    {18.0, 1.42, 1577.34},
    {10.0, 5.91, 10977.08},
    {9.0, 1.42, 6275.96},
    {9.0, 0.27, 5486.78},
}};
constexpr std::array<Term, 6> kR2{{
    {4359.0, 5.7846, 6283.0758},
    {124.0, 5.579, 12566.152},
    {12.0, 3.14, 0.0},
    {9.0, 3.63, 77713.77},
    {6.0, 1.87, 5573.14},
    {3.0, 5.47, 18849.23},
}};
constexpr std::array<Term, 2> kR3{{
    {145.0, 4.273, 6283.076},
    {7.0, 3.92, 12566.15},
}};
constexpr std::array<Term, 1> kR4{{
    {4.0, 2.56, 6283.08},
}};
constexpr std::array<Nutation, 63> kNutation{{
    {{0, 0, 0, 0, 1}, -171996.0, -174.2, 92025.0, 8.9},
    {{-2, 0, 0, 2, 2}, -13187.0, -1.6, 5736.0, -3.1},
    {{0, 0, 0, 2, 2}, -2274.0, -0.2, 977.0, -0.5},
    {{0, 0, 0, 0, 2}, 2062.0, 0.2, -895.0, 0.5},
    {{0, 1, 0, 0, 0}, 1426.0, -3.4, 54.0, -0.1},
    {{0, 0, 1, 0, 0}, 712.0, 0.1, -7.0, 0.0},
    {{-2, 1, 0, 2, 2}, -517.0, 1.2, 224.0, -0.6},
    {{0, 0, 0, 2, 1}, -386.0, -0.4, 200.0, 0.0},
    {{0, 0, 1, 2, 2}, -301.0, 0.0, 129.0, -0.1},
    {{-2, -1, 0, 2, 2}, 217.0, -0.5, -95.0, 0.3},
    {{-2, 0, 1, 0, 0}, -158.0, 0.0, 0.0, 0.0},
    {{-2, 0, 0, 2, 1}, 129.0, 0.1, -70.0, 0.0},
    {{0, 0, -1, 2, 2}, 123.0, 0.0, -53.0, 0.0},
    {{2, 0, 0, 0, 0}, 63.0, 0.0, 0.0, 0.0},
    {{0, 0, 1, 0, 1}, 63.0, 0.1, -33.0, 0.0},
    {{2, 0, -1, 2, 2}, -59.0, 0.0, 26.0, 0.0},
    {{0, 0, -1, 0, 1}, -58.0, -0.1, 32.0, 0.0},
    {{0, 0, 1, 2, 1}, -51.0, 0.0, 27.0, 0.0},
    {{-2, 0, 2, 0, 0}, 48.0, 0.0, 0.0, 0.0},
    {{0, 0, -2, 2, 1}, 46.0, 0.0, -24.0, 0.0},
    {{2, 0, 0, 2, 2}, -38.0, 0.0, 16.0, 0.0},
    {{0, 0, 2, 2, 2}, -31.0, 0.0, 13.0, 0.0},
    {{0, 0, 2, 0, 0}, 29.0, 0.0, 0.0, 0.0},
    {{-2, 0, 1, 2, 2}, 29.0, 0.0, -12.0, 0.0},
    {{0, 0, 0, 2, 0}, 26.0, 0.0, 0.0, 0.0},
    {{-2, 0, 0, 2, 0}, -22.0, 0.0, 0.0, 0.0},
    {{0, 0, -1, 2, 1}, 21.0, 0.0, -10.0, 0.0},
    {{0, 2, 0, 0, 0}, 17.0, -0.1, 0.0, 0.0},
    {{2, 0, -1, 0, 1}, 16.0, 0.0, -8.0, 0.0},
    {{-2, 2, 0, 2, 2}, -16.0, 0.1, 7.0, 0.0},
    {{0, 1, 0, 0, 1}, -15.0, 0.0, 9.0, 0.0},
    {{-2, 0, 1, 0, 1}, -13.0, 0.0, 7.0, 0.0},
    {{0, -1, 0, 0, 1}, -12.0, 0.0, 6.0, 0.0},
    {{0, 0, 2, -2, 0}, 11.0, 0.0, 0.0, 0.0},
    {{2, 0, -1, 2, 1}, -10.0, 0.0, 5.0, 0.0},
    {{2, 0, 1, 2, 2}, -8.0, 0.0, 3.0, 0.0},
    {{0, 1, 0, 2, 2}, 7.0, 0.0, -3.0, 0.0},
    {{-2, 1, 1, 0, 0}, -7.0, 0.0, 0.0, 0.0},
    {{0, -1, 0, 2, 2}, -7.0, 0.0, 3.0, 0.0},
    {{2, 0, 0, 2, 1}, -7.0, 0.0, 3.0, 0.0},
    {{2, 0, 1, 0, 0}, 6.0, 0.0, 0.0, 0.0},
    {{-2, 0, 2, 2, 2}, 6.0, 0.0, -3.0, 0.0},
    {{-2, 0, 1, 2, 1}, 6.0, 0.0, -3.0, 0.0},
    {{2, 0, -2, 0, 1}, -6.0, 0.0, 3.0, 0.0},
    {{2, 0, 0, 0, 1}, -6.0, 0.0, 3.0, 0.0},
    {{0, -1, 1, 0, 0}, 5.0, 0.0, 0.0, 0.0},
    {{-2, -1, 0, 2, 1}, -5.0, 0.0, 3.0, 0.0},
    {{-2, 0, 0, 0, 1}, -5.0, 0.0, 3.0, 0.0},
    {{0, 0, 2, 2, 1}, -5.0, 0.0, 3.0, 0.0},
    {{-2, 0, 2, 0, 1}, 4.0, 0.0, 0.0, 0.0},
    {{-2, 1, 0, 2, 1}, 4.0, 0.0, 0.0, 0.0},
    {{0, 0, 1, -2, 0}, 4.0, 0.0, 0.0, 0.0},
    {{-1, 0, 1, 0, 0}, -4.0, 0.0, 0.0, 0.0},
    {{-2, 1, 0, 0, 0}, -4.0, 0.0, 0.0, 0.0},
    {{1, 0, 0, 0, 0}, -4.0, 0.0, 0.0, 0.0},
    {{0, 0, 1, 2, 0}, 3.0, 0.0, 0.0, 0.0},
    {{0, 0, -2, 2, 2}, -3.0, 0.0, 0.0, 0.0},
    {{-1, -1, 1, 0, 0}, -3.0, 0.0, 0.0, 0.0},
    {{0, 1, 1, 0, 0}, -3.0, 0.0, 0.0, 0.0},
    {{0, -1, 1, 2, 2}, -3.0, 0.0, 0.0, 0.0},
    {{2, -1, -1, 2, 2}, -3.0, 0.0, 0.0, 0.0},
    {{0, 0, 3, 2, 2}, -3.0, 0.0, 0.0, 0.0},
    {{2, -1, 0, 2, 2}, -3.0, 0.0, 0.0, 0.0},
}};

template <std::size_t N>
void require_terms(std::span<const spa::PeriodicTerm> got, const std::array<Term, N>& want,
                   const char* group) {
    INFO(group);
    REQUIRE(got.size() == N);
    for (std::size_t i = 0; i < N; ++i) {
        INFO("row " << i);
        CHECK(got[i].a == want[i].a);
        CHECK(got[i].b == want[i].b);
        CHECK(got[i].c == want[i].c);
    }
}

}  // namespace

// --- J: Julian day ----------------------------------------------------------

TEST_CASE("J1 Julian day at the epoch, J2000 and the NREL case", "[radiation][time]") {
    CHECK(julian_day(at(0)) == 2440587.5);
    CHECK(julian_day(at(946728000000)) == 2451545.0);  // 2000-01-01T12:00Z
    // Published to six decimals: 2452930.312847.
    CHECK_THAT(julian_day(at(kNrelMs)), WithinAbs(2452930.312847, 5e-7));
}

TEST_CASE("J2 Julian day before the epoch: one millisecond before is not floored to a second",
          "[radiation][time]") {
    // At JD 2.44e6 a double resolves 4.7e-10 days (40 microseconds).
    CHECK_THAT(julian_day(at(-1)), WithinAbs(2440587.5 - 1.0 / 86400000.0, 1e-9));
    CHECK_THAT(julian_day(at(-86400000)), WithinAbs(2440586.5, 1e-9));
}

TEST_CASE("J3 UTC milliseconds round-trip through the Julian day exactly", "[radiation][time]") {
    std::mt19937_64 rng(462);
    // Years 1800 to 2200 (ms since the epoch), plus the edges of the epoch.
    std::uniform_int_distribution<std::int64_t> ms(-5364662400000, 7258118400000);
    for (std::int64_t t : {std::int64_t{-1}, std::int64_t{0}, std::int64_t{1}, kNrelMs}) {
        INFO("ms " << t);
        CHECK(utc_from_julian_day(julian_day(at(t))) == at(t));
    }
    for (int i = 0; i < 10000; ++i) {
        const std::int64_t t = ms(rng);
        INFO("ms " << t);
        REQUIRE(utc_from_julian_day(julian_day(at(t))) == at(t));
    }
}

// --- N: NREL's test case ----------------------------------------------------

TEST_CASE("N1 the geocentric sun of the NREL case", "[radiation][spa]") {
    const GeocentricSun g = geocentric_sun(at(kNrelMs), kDeltaT);
    // Published: alpha 202.22741, delta -9.31434, R 0.9965422974 AU,
    // local hour angle H = nu + longitude - alpha = 11.105900 (spa.c: 11.105902).
    CHECK_THAT(g.alpha_deg, WithinAbs(202.22741, 5e-6));
    CHECK_THAT(g.delta_deg, WithinAbs(-9.31434, 5e-6));
    CHECK_THAT(g.r_au, WithinAbs(0.9965422974, 5e-11));
    const double h = std::fmod(g.nu_deg + kLon - g.alpha_deg + 720.0, 360.0);
    CHECK_THAT(h, WithinAbs(11.10590, 5e-6));
}

TEST_CASE("N2 the reference path reproduces the published zenith and azimuth",
          "[radiation][spa]") {
    const SolarPosition s =
        solar_position(at(kNrelMs), kLat, kLon, kHeight, kDeltaT, kPressure, kTemperature);
    // Published: topocentric zenith 50.11162 (refracted), azimuth 194.34024
    // from north, eastward. e0 without refraction, 39.872046, is pvlib's SPA.
    CHECK_THAT(90.0 - s.e_deg, WithinAbs(50.11162, 5e-6));
    CHECK(azimuth_gap(s.azimuth_deg, 194.34024) <= 5e-6);
    CHECK_THAT(s.e0_deg, WithinAbs(39.872046, 5e-7));
    CHECK_THAT(s.e_deg - s.e0_deg, WithinAbs(bennett(s.e0_deg, kPressure, kTemperature), 1e-12));
}

TEST_CASE("N3 the split (geocentric, then topocentric) equals the reference path bit for bit",
          "[radiation][spa]") {
    std::mt19937_64 rng(463);
    std::uniform_int_distribution<std::int64_t> ms(-2208988800000, 4102444800000);  // 1900-2100
    std::uniform_real_distribution<double> lat(-89.0, 89.0), lon(-180.0, 180.0), h(0.0, 3000.0);
    for (int i = 0; i < 1000; ++i) {
        const utc_ms t = at(ms(rng));
        const double la = lat(rng), lo = lon(rng), z = h(rng);
        const TopocentricSun split = topocentric_sun(geocentric_sun(t, 69.0), la, lo, z);
        const SolarPosition ref = solar_position(t, la, lo, z, 69.0, 1010.0, 10.0);
        INFO("case " << i);
        REQUIRE(split.e0_deg == ref.e0_deg);
        REQUIRE(split.azimuth_deg == ref.azimuth_deg);
    }
}

// --- S: the term tables -----------------------------------------------------

TEST_CASE("S1 the periodic-term tables equal SPA's, term by term", "[radiation][spa][tables]") {
    require_terms(spa::l_terms(0), kL0, "L0");
    require_terms(spa::l_terms(1), kL1, "L1");
    require_terms(spa::l_terms(2), kL2, "L2");
    require_terms(spa::l_terms(3), kL3, "L3");
    require_terms(spa::l_terms(4), kL4, "L4");
    require_terms(spa::l_terms(5), kL5, "L5");
    require_terms(spa::b_terms(0), kB0, "B0");
    require_terms(spa::b_terms(1), kB1, "B1");
    require_terms(spa::r_terms(0), kR0, "R0");
    require_terms(spa::r_terms(1), kR1, "R1");
    require_terms(spa::r_terms(2), kR2, "R2");
    require_terms(spa::r_terms(3), kR3, "R3");
    require_terms(spa::r_terms(4), kR4, "R4");

    const std::span<const spa::NutationTerm> nut = spa::nutation_terms();
    REQUIRE(nut.size() == kNutation.size());
    for (std::size_t i = 0; i < kNutation.size(); ++i) {
        INFO("nutation row " << i);
        CHECK(nut[i].y == kNutation[i].y);
        CHECK(nut[i].a == kNutation[i].a);
        CHECK(nut[i].b == kNutation[i].b);
        CHECK(nut[i].c == kNutation[i].c);
        CHECK(nut[i].d == kNutation[i].d);
    }
}

// --- D: delta T -------------------------------------------------------------

TEST_CASE("D1 the default delta T follows Espenak and Meeus, /100 parabola outside 1986-2050",
          "[radiation][time][delta_t]") {
    // Mid-January dates, where the decimal year y = year + (month - 0.5)/12 and
    // the whole year give values 0.04 s apart at most; the bound covers both.
    CHECK_THAT(delta_t_seconds(at(1768478400000)), WithinAbs(75.10, 0.1));   // 2026-01-15
    CHECK_THAT(delta_t_seconds(at(947937600000)), WithinAbs(63.87, 0.1));    // 2000-01-15
    CHECK_THAT(delta_t_seconds(at(790171200000)), WithinAbs(60.83, 0.1));    // 1995-01-15
    // 2200-01-15: -20 + 32 ((y - 1820)/100)^2. The legacy's /20 gives 11532 s.
    CHECK_THAT(delta_t_seconds(at(7259371200000)), WithinAbs(442.18, 0.15));
}

TEST_CASE("D2 delta T at 2026 is not the legacy's 3375 s: Oslo's sun", "[radiation][time][delta_t]") {
    const utc_ms t = at(1768478400000);  // 2026-01-15T12:00Z, Oslo 59.91 N 10.75 E
    CHECK(delta_t_seconds(t) < 100.0);
    // The legacy's measured outputs: delta T = 69 s gives e0 8.72791, azimuth
    // 187.92892; its defective 3374.9 s gives 8.73784, 187.89079.
    const SolarPosition s69 = solar_position(t, 59.91, 10.75, 0.0, 69.0, 1010.0, 10.0);
    CHECK_THAT(s69.e0_deg, WithinAbs(8.72791, 5e-6));
    CHECK_THAT(s69.azimuth_deg, WithinAbs(187.92892, 5e-6));
    const SolarPosition sd =
        solar_position(t, 59.91, 10.75, 0.0, delta_t_seconds(t), 1010.0, 10.0);
    CHECK_THAT(sd.e0_deg, WithinAbs(8.72791, 1e-4));
    CHECK_THAT(sd.azimuth_deg, WithinAbs(187.92892, 1e-4));
    CHECK(std::abs(sd.e0_deg - 8.73784) > 5e-3);
    CHECK(std::abs(sd.azimuth_deg - 187.89079) > 2e-2);
}

// --- A: refraction ----------------------------------------------------------

TEST_CASE("A1 eq. 42 at 1010 hPa and 10 C, at and above -1.8 degrees", "[radiation][refraction]") {
    // The design's table (terrain-radiation.md, "What was measured", P3).
    const std::array<std::pair<double, double>, 7> table{{
        {-1.0, 0.647}, {-0.5667, 0.573}, {0.0, 0.483}, {1.0, 0.362},
        {3.0, 0.228},  {10.0, 0.090},    {-1.8, 0.742},
    }};
    for (auto [e0, r] : table) {
        INFO("e0 " << e0);
        CHECK_THAT(refraction(e0, 1010.0, 10.0), WithinAbs(r, 5e-4));
        CHECK_THAT(refraction(e0, 1010.0, 10.0), WithinAbs(bennett(e0, 1010.0, 10.0), 1e-12));
    }
    // Pressure and temperature scale it: P/1010 and 283/(273 + T).
    CHECK_THAT(refraction(5.0, 505.0, 10.0), WithinAbs(0.5 * refraction(5.0, 1010.0, 10.0), 1e-12));
    CHECK_THAT(refraction(5.0, 1010.0, -10.0),
               WithinAbs(283.0 / 263.0 * refraction(5.0, 1010.0, 10.0), 1e-12));
}

TEST_CASE("A2 below -1.8 degrees the refraction is held at its -1.8 degree value",
          "[radiation][refraction]") {
    const double held = bennett(-1.8, 820.0, 11.0);
    // -5.11 is the formula's own pole; -3 and -4 are where it turns
    // non-physical (0.518, 0.184 at 1010 hPa).
    for (double e0 : {-1.8000001, -2.0, -3.0, -4.0, -5.11, -10.0, -45.0, -90.0}) {
        INFO("e0 " << e0);
        const double r = refraction(e0, 820.0, 11.0);
        REQUIRE(std::isfinite(r));
        CHECK_THAT(r, WithinAbs(held, 1e-12));
    }
}

TEST_CASE("A3 the guarded refraction never decreases as the sun rises, -90 to 90 degrees",
          "[radiation][refraction]") {
    // Eq. 42 crosses zero at e0 = 89.8916 (where e0 + 10.3/(e0 + 5.11) = 90)
    // and is -3.2e-5 degrees at the zenith: still non-increasing.
    double previous = refraction(-90.0, 1010.0, 10.0);
    for (int k = -9000; k <= 9000; ++k) {
        const double e0 = k / 100.0;
        const double r = refraction(e0, 1010.0, 10.0);
        INFO("e0 " << e0);
        REQUIRE(r <= previous + 1e-15);
        previous = r;
    }
}

// --- T: sunrise, transit, sunset --------------------------------------------

TEST_CASE("T1 sunrise, transit and sunset of the NREL case on the reference path",
          "[radiation][spa]") {
    auto e0 = [](std::int64_t ms) {
        return solar_position(at(ms), kLat, kLon, kHeight, kDeltaT, kPressure, kTemperature)
                   .e0_deg +
               0.8333;  // SPA's h0': the sun's radius plus standard refraction
    };
    auto azimuth = [](std::int64_t ms) {
        return solar_position(at(ms), kLat, kLon, kHeight, kDeltaT, kPressure, kTemperature)
                   .azimuth_deg -
               180.0;
    };
    constexpr std::int64_t hour = 3600000;
    // Published, local (UTC-7): sunrise 06:12:43, transit 11:46:04. The
    // reference path's crossings lie 1.3 s and 1.0 s after them (pvlib's SPA,
    // bisected the same way).
    const std::int64_t rise = crossing(e0, kNrelDay + 12 * hour, kNrelDay + 15 * hour);
    const std::int64_t transit = crossing(azimuth, kNrelDay + 17 * hour, kNrelDay + 20 * hour);
    CHECK(std::abs(rise - (kNrelDay + 47563000)) <= 3000);
    CHECK(std::abs(transit - (kNrelDay + 67564000)) <= 2000);
    // SPA's published sunset, 17:20:19 local, is the UT day's: 16 October's
    // evening (00:20:19 UTC on the 17th); the reference path crosses h0' at
    // 17:20:18.40 (pvlib's SPA, bisected the same way).
    const std::int64_t set = crossing(e0, kNrelDay - hour, kNrelDay + 2 * hour);
    CHECK(std::abs(set - (kNrelDay + 1219000)) <= 2000);
}

// --- G: ranges --------------------------------------------------------------

TEST_CASE("G1 elevations and azimuths stay in range over the years -2000 to 6000",
          "[radiation][spa]") {
    std::mt19937_64 rng(4646);
    // -2000-01-01 to 6000-12-31, proleptic Gregorian, in ms since the epoch.
    std::uniform_int_distribution<std::int64_t> ms(-125281123200000, 127206028800000);
    std::uniform_real_distribution<double> lat(-90.0, 90.0), lon(-180.0, 180.0), h(-400.0, 9000.0);
    for (int i = 0; i < 2000; ++i) {
        const utc_ms t = at(ms(rng));
        const SolarPosition s = solar_position(t, lat(rng), lon(rng), h(rng),
                                               delta_t_seconds(t), 1010.0, 10.0);
        INFO("case " << i);
        REQUIRE(std::isfinite(s.e0_deg));
        REQUIRE(s.e0_deg >= -90.0);
        REQUIRE(s.e0_deg <= 90.0);
        REQUIRE(s.azimuth_deg >= 0.0);
        REQUIRE(s.azimuth_deg < 360.0);
        REQUIRE(std::isfinite(s.e_deg));
    }
}
