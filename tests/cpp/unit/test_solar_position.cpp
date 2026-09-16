#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <rasputin/solar_position.h>

#include <chrono>

using Catch::Matchers::WithinAbs;

using namespace std::chrono;

using rasputin::solar_position::jd_from_cal;
using rasputin::solar_position::jd_from_clock;
using rasputin::solar_position::delta_t_calculator::coarse_timestamp_calc;

namespace {

// Julian Day is a fractional day count, so 1e-9 d is ~86 us -- far tighter
// than the second-level resolution these conversions actually carry.
constexpr double jd_tol = 1e-9;

// These entry points take a system_clock time_point. Every one of them is UTC
// by construction: no zone lookup, no local time, no tzdb.
constexpr system_clock::time_point at(sys_days day) {
    return system_clock::time_point{day.time_since_epoch()};
}

} // namespace

TEST_CASE("jd_from_clock anchors the Unix epoch at JD 2440587.5", "[solar][julian]") {
    REQUIRE_THAT(jd_from_clock(at(sys_days{January / 1 / 1970})),
                 WithinAbs(2440587.5, jd_tol));
}

TEST_CASE("jd_from_clock reproduces the J2000.0 epoch", "[solar][julian]") {
    // J2000.0 is 2000-01-01T12:00:00Z == JD 2451545.0, the standard astronomical epoch.
    REQUIRE_THAT(jd_from_clock(at(sys_days{January / 1 / 2000}) + hours{12}),
                 WithinAbs(2451545.0, jd_tol));
}

TEST_CASE("jd_from_clock handles pre-epoch time points", "[solar][julian][edge]") {
    // Negative durations must not truncate toward zero.
    REQUIRE_THAT(jd_from_clock(at(sys_days{December / 31 / 1969})),
                 WithinAbs(2440586.5, jd_tol));
}

TEST_CASE("jd_from_clock advances exactly one Julian Day per 86400 s", "[solar][julian]") {
    const auto base = at(sys_days{March / 14 / 2021});
    const double jd0 = jd_from_clock(base);

    REQUIRE_THAT(jd_from_clock(base + days{1}) - jd0, WithinAbs(1.0, jd_tol));
    REQUIRE_THAT(jd_from_clock(base + hours{12}) - jd0, WithinAbs(0.5, jd_tol));
    REQUIRE_THAT(jd_from_clock(base + seconds{86400}) - jd0, WithinAbs(1.0, jd_tol));
}

TEST_CASE("jd_from_cal agrees with jd_from_clock on the same instant", "[solar][julian]") {
    // Fractional day 1.5 == the 12:00 mark on day 1.
    REQUIRE_THAT(jd_from_cal(2000, 1, 1.5), WithinAbs(2451545.0, jd_tol));
    REQUIRE_THAT(jd_from_cal(2000, 1, 1.5),
                 WithinAbs(jd_from_clock(at(sys_days{January / 1 / 2000}) + hours{12}), jd_tol));

    REQUIRE_THAT(jd_from_cal(1970, 1, 1.5),
                 WithinAbs(jd_from_clock(at(sys_days{January / 1 / 1970}) + hours{12}), jd_tol));
}

TEST_CASE("jd_from_cal takes the Julian-calendar branch before the 1582 reform", "[solar][julian][edge]") {
    // jd_temp < 2299160 skips the Gregorian century correction.
    REQUIRE_THAT(jd_from_cal(1582, 10, 3.5), WithinAbs(2299159.0, jd_tol));
}

TEST_CASE("jd_from_cal normalises January and February into the prior year", "[solar][julian]") {
    // month < 3 shifts to month+12 of year-1; the result must stay continuous
    // across the year boundary.
    const double dec31 = jd_from_cal(1999, 12, 31.5);
    const double jan01 = jd_from_cal(2000, 1, 1.5);
    REQUIRE_THAT(jan01 - dec31, WithinAbs(1.0, jd_tol));
}

TEST_CASE("coarse_timestamp_calc selects the tabulated delta-T bucket", "[solar][deltat]") {
    const auto dt = coarse_timestamp_calc();

    // 1970 falls in the 1955-2009 table at index (1970-1955)/5 == 3.
    REQUIRE_THAT(dt(at(sys_days{January / 1 / 1970})), WithinAbs(40.2, 1e-12));
}

TEST_CASE("coarse_timestamp_calc falls back to the parabolic model outside 1955-2009", "[solar][deltat][edge]") {
    const auto dt = coarse_timestamp_calc();

    // -20 + 32*((1900-1820)/20)^2 == 492.
    REQUIRE_THAT(dt(at(sys_days{January / 1 / 1900})), WithinAbs(492.0, 1e-12));
}

TEST_CASE("coarse_timestamp_calc rounds the year rather than flooring it", "[solar][deltat][characterisation]") {
    const auto dt = coarse_timestamp_calc();

    // CHARACTERISATION, not endorsement. The year is derived as
    // round(elapsed / mean_tropical_year), so a date late in calendar 2004
    // resolves to 2005 and picks bucket 10 (64.7) instead of bucket 9 (63.8).
    // A port that switches to floor() will change this value -- deliberately.
    REQUIRE_THAT(dt(at(sys_days{December / 1 / 2004})), WithinAbs(64.7, 1e-12));
}
