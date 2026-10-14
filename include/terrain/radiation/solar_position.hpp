// The sun's position by NREL's Solar Position Algorithm (Reda and Andreas
// 2004, "Solar position algorithm for solar radiation applications"), in its
// order and notation; the numbers in comments are its equations. Valid for the
// years -2000 to 6000, about 0.0003 degrees.
//
// Split in two: `geocentric_sun` depends on time only and carries the cost
// (about 200 periodic terms); `topocentric_sun` is the cheap part per place.
// `solar_position` runs both, then refraction: the reference path, point by
// point. Every angle is in degrees, every time UTC.

#pragma once

#include <terrain/radiation/spa_tables.hpp>
#include <terrain/radiation/time.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>

namespace terrain::radiation {

namespace spa {

// The periodic terms of L, B and R by power of jme, and the nutation terms.
inline std::span<const PeriodicTerm> l_terms(std::size_t k) {
    using namespace tables;
    static constexpr std::array<std::span<const PeriodicTerm>, 6> groups{kL0, kL1, kL2, kL3, kL4, kL5};
    return groups.at(k);
}
inline std::span<const PeriodicTerm> b_terms(std::size_t k) {
    static constexpr std::array<std::span<const PeriodicTerm>, 2> groups{tables::kB0, tables::kB1};
    return groups.at(k);
}
inline std::span<const PeriodicTerm> r_terms(std::size_t k) {
    using namespace tables;
    static constexpr std::array<std::span<const PeriodicTerm>, 5> groups{kR0, kR1, kR2, kR3, kR4};
    return groups.at(k);
}
inline std::span<const NutationTerm> nutation_terms() { return tables::kNutation; }

inline constexpr double kRad = std::numbers::pi / 180.0;

// x mod 360 in [0, 360).
inline double wrap360(double x) {
    double r = std::fmod(x, 360.0);
    if (r < 0.0)
        r += 360.0;
    return r >= 360.0 ? 0.0 : r;
}

// Eqs. 9-12: sum_k (sum_i A_i cos(B_i + C_i jme)) jme^k / 1e8.
template <std::size_t N, class Groups>
double series(Groups groups, double jme) {
    double total = 0.0;
    for (std::size_t k = N; k-- > 0;) {
        double s = 0.0;
        for (const PeriodicTerm& t : groups(k))
            s += t.a * std::cos(t.b + t.c * jme);
        total = total * jme + s;
    }
    return total / 1e8;
}

}  // namespace spa

// The sun seen from the Earth's centre: right ascension alpha, declination
// delta, apparent sidereal time at Greenwich nu, and the distance R in AU.
struct GeocentricSun {
    double alpha_deg, delta_deg, nu_deg, r_au;
};

// The sun seen from a place: elevation without refraction, and azimuth from
// true north, clockwise.
struct TopocentricSun {
    double e0_deg, azimuth_deg;
};

// The reference path's result: e0 without refraction, e with it.
struct SolarPosition {
    double e0_deg, e_deg, azimuth_deg;
};

inline GeocentricSun geocentric_sun(utc_ms t, double delta_t_s) {
    using namespace spa;
    const double jd = julian_day(t);                       // 4
    const double jde = jd + delta_t_s / 86400.0;           // 5
    const double jc = (jd - 2451545.0) / 36525.0;          // 6
    const double jce = (jde - 2451545.0) / 36525.0;        // 7
    const double jme = jce / 10.0;                         // 8
    const double l = wrap360(series<6>(l_terms, jme) / kRad);  // 9-11
    const double b = series<2>(b_terms, jme) / kRad;
    const double r = series<5>(r_terms, jme);
    const double theta = wrap360(l + 180.0);  // 13
    const double beta = -b;                   // 14
    const double jce2 = jce * jce, jce3 = jce2 * jce;
    const std::array<double, 5> x{
        297.85036 + 445267.111480 * jce - 0.0019142 * jce2 + jce3 / 189474.0,  // 15
        357.52772 + 35999.050340 * jce - 0.0001603 * jce2 - jce3 / 300000.0,   // 16
        134.96298 + 477198.867398 * jce + 0.0086972 * jce2 + jce3 / 56250.0,   // 17
        93.27191 + 483202.017538 * jce - 0.0036825 * jce2 + jce3 / 327270.0,   // 18
        125.04452 - 1934.136261 * jce + 0.0020708 * jce2 + jce3 / 450000.0,    // 19
    };
    double dpsi = 0.0, deps = 0.0;  // 20-23
    for (const NutationTerm& n : nutation_terms()) {
        double arg = 0.0;
        for (std::size_t j = 0; j < 5; ++j)
            arg += n.y[j] * x[j];
        dpsi += (n.a + n.b * jce) * std::sin(arg * kRad);
        deps += (n.c + n.d * jce) * std::cos(arg * kRad);
    }
    dpsi /= 36000000.0;
    deps /= 36000000.0;
    const double u = jme / 10.0;  // 24
    const double eps0 =
        84381.448 +
        u * (-4680.93 + u * (-1.55 + u * (1999.25 + u * (-51.38 + u * (-249.67 + u * (-39.05 + u * (7.12 + u * (27.87 + u * (5.79 + u * 2.45)))))))));
    const double eps = eps0 / 3600.0 + deps;  // 25
    const double dtau = -20.4898 / (3600.0 * r);  // 26
    const double lambda = theta + dpsi + dtau;    // 27
    const double nu0 = wrap360(280.46061837 + 360.98564736629 * (jd - 2451545.0) +
                               0.000387933 * jc * jc - jc * jc * jc / 38710000.0);  // 28
    const double nu = nu0 + dpsi * std::cos(eps * kRad);                          // 29
    const double sl = std::sin(lambda * kRad), ce = std::cos(eps * kRad), se = std::sin(eps * kRad);
    const double alpha = wrap360(std::atan2(sl * ce - std::tan(beta * kRad) * se, std::cos(lambda * kRad)) / kRad);  // 30
    const double delta = std::asin(std::sin(beta * kRad) * ce + std::cos(beta * kRad) * se * sl) / kRad;            // 31
    return {alpha, delta, nu, r};
}

inline TopocentricSun topocentric_sun(const GeocentricSun& g, double lat_deg, double lon_deg, double height_m) {
    using spa::kRad;
    const double h = spa::wrap360(g.nu_deg + lon_deg - g.alpha_deg) * kRad;  // 32
    const double xi = 8.794 / (3600.0 * g.r_au) * kRad;                      // 33
    const double phi = lat_deg * kRad;
    const double u = std::atan(0.99664719 * std::tan(phi));                        // 34
    const double x = std::cos(u) + height_m / 6378140.0 * std::cos(phi);           // 35
    const double y = 0.99664719 * std::sin(u) + height_m / 6378140.0 * std::sin(phi);  // 36
    const double delta = g.delta_deg * kRad;
    const double den = std::cos(delta) - x * std::sin(xi) * std::cos(h);
    const double dalpha = std::atan2(-x * std::sin(xi) * std::sin(h), den);                          // 37
    const double delta_p = std::atan2((std::sin(delta) - y * std::sin(xi)) * std::cos(dalpha), den);  // 39
    const double h_p = h - dalpha;                                                                    // 40
    const double s = std::sin(phi) * std::sin(delta_p) + std::cos(phi) * std::cos(delta_p) * std::cos(h_p);
    const double e0 = std::asin(std::clamp(s, -1.0, 1.0)) / kRad;  // 41
    const double gamma =
        std::atan2(std::sin(h_p), std::cos(h_p) * std::sin(phi) - std::tan(delta_p) * std::cos(phi)) / kRad;  // 45
    return {e0, spa::wrap360(gamma + 180.0)};                                                                  // 46
}

// Eq. 42 at pressure P (hPa) and temperature T (degrees C), added to e0. Below
// e0 = -1.8 degrees the formula turns non-physical (a pole at -5.11), so the
// refraction is held at its -1.8 degree value there; a summit's horizon dips to
// about -1.5 degrees, so SPA's own cut at -0.8333 is too high.
inline double apparent_elevation(double e0_deg, double pressure_hpa, double temperature_c) {
    const double e = std::max(e0_deg, -1.8);
    return e0_deg + (pressure_hpa / 1010.0) * (283.0 / (273.0 + temperature_c)) * 1.02 /
                        (60.0 * std::tan((e + 10.3 / (e + 5.11)) * spa::kRad));
}

inline SolarPosition solar_position(utc_ms t, double lat_deg, double lon_deg, double height_m, double delta_t_s,
                                    double pressure_hpa, double temperature_c) {
    const TopocentricSun s = topocentric_sun(geocentric_sun(t, delta_t_s), lat_deg, lon_deg, height_m);
    return {s.e0_deg, apparent_elevation(s.e0_deg, pressure_hpa, temperature_c), s.azimuth_deg};
}

}  // namespace terrain::radiation
