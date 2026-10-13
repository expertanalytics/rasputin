// Time for the sun's position: UTC as Unix time in milliseconds, the Julian
// day, and a default for delta T (TT - UT). Only <chrono>; nothing here leaves
// UTC, so no time-zone database is involved. Unix time has no leap seconds;
// UTC stands in for UT1, which it follows within 0.9 s (0.004 degrees of the
// Earth's rotation).

#pragma once

#include <chrono>
#include <cmath>
#include <cstdint>

namespace terrain::radiation {

using utc_ms = std::chrono::sys_time<std::chrono::milliseconds>;

inline constexpr double kUnixEpochJulianDay = 2440587.5;
inline constexpr double kMsPerDay = 86400000.0;

// The Julian day (UT) of a UTC instant. A double near 2.4e6 resolves 40
// microseconds, so the millisecond survives the round trip below.
inline double julian_day(utc_ms t) {
    return kUnixEpochJulianDay + static_cast<double>(t.time_since_epoch().count()) / kMsPerDay;
}

// The UTC instant of a Julian day, rounded to the nearest millisecond.
inline utc_ms utc_from_julian_day(double jd) {
    const auto ms = std::llround((jd - kUnixEpochJulianDay) * kMsPerDay);
    return utc_ms{std::chrono::milliseconds{ms}};
}

// Espenak and Meeus's delta T in seconds at the decimal year
// y = year + (month - 0.5) / 12: their 1986-2005 and 2005-2050 polynomials,
// and the long-term parabola -20 + 32 ((y - 1820) / 100)^2 outside them.
// A consumer with IERS values passes those instead.
inline double delta_t_seconds(utc_ms t) {
    using namespace std::chrono;
    const year_month_day ymd{floor<days>(t)};
    const double y = static_cast<int>(ymd.year()) + (static_cast<unsigned>(ymd.month()) - 0.5) / 12.0;
    const double u = y - 2000.0;
    if (y >= 1986.0 && y < 2005.0)
        return 63.86 + u * (0.3345 + u * (-0.060374 + u * (0.0017275 + u * (0.000651814 + u * 0.00002373599))));
    if (y >= 2005.0 && y < 2050.0)
        return 62.92 + u * (0.32217 + u * 0.005589);
    const double c = (y - 1820.0) / 100.0;
    return -20.0 + 32.0 * c * c;
}

}  // namespace terrain::radiation
