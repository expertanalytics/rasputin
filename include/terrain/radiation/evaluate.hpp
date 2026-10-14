// The instant evaluator: per receiver and UTC time, whether the sun is above
// the receiver's stored horizon and in front of its plane (lit), the cosine of
// the incidence angle (cos_i), the beam factor b = lit max(0, cos_i) that
// multiplies direct normal irradiance, and the sine of the apparent sun
// elevation (sin_e). The geocentric sun is computed once per time; per
// receiver it is seen along the receiver's east-north-up basis, without the
// parallax of at most 0.0024 degrees that the reference path keeps, and
// refracted as the reference path refracts it. Horizons are int16 in 0.01
// degrees, direction k at true azimuth k 360 / N.

#pragma once

#include <terrain/parallel_util/chunks.hpp>
#include <terrain/radiation/receivers.hpp>
#include <terrain/radiation/solar_position.hpp>
#include <terrain/radiation/time.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace terrain::radiation {

// The upper limb's lift over the sun's centre, degrees: the lit test counts
// the sun as risen when its upper edge clears the horizon.
inline constexpr double kUpperLimbDeg = 0.2667;

// The horizon at a true azimuth, degrees: the stored value on a direction,
// linear between the two neighbouring directions otherwise, direction N - 1
// next to direction 0.
inline double horizon_at(std::span<const std::int16_t> row, double true_azimuth_deg) noexcept {
    const std::size_t n = row.size();
    double x = std::fmod(true_azimuth_deg, 360.0) * static_cast<double>(n) / 360.0;
    if (x < 0.0)
        x += static_cast<double>(n);
    const double k = std::floor(x), f = x - k;
    const std::size_t k0 = static_cast<std::size_t>(k) % n, k1 = (k0 + 1) % n;
    const double h0 = row[k0] / 100.0, h1 = row[k1] / 100.0;
    return h0 + f * (h1 - h0);
}

struct SampleFactors {
    bool lit;
    double cos_i, b, sin_e;
};

namespace detail {

// The factors for an apparent sun given as its elevation, its true azimuth and
// that azimuth's sine and cosine.
inline SampleFactors factors(std::span<const std::int16_t> row, std::span<const double, 3> n, double e_deg,
                             double a_deg, double sin_a, double cos_a, double limb_deg) noexcept {
    const double e = e_deg * spa::kRad, ce = std::cos(e), se = std::sin(e);
    const double cos_i = n[0] * ce * sin_a + n[1] * ce * cos_a + n[2] * se;
    const bool lit = e_deg + limb_deg > horizon_at(row, a_deg) && cos_i > 0.0;
    return {lit, cos_i, lit ? std::max(0.0, cos_i) : 0.0, se};
}

}  // namespace detail

// One receiver and one apparent sun: lit is strict (a sun exactly on the
// horizon is not lit) and needs n.s > 0; cos_i and sin_e are unclipped.
inline SampleFactors sample_factors(std::span<const std::int16_t> row, std::span<const double, 3> normal_enu,
                                    double apparent_elevation_deg, double true_azimuth_deg,
                                    double limb_deg) noexcept {
    const double a = true_azimuth_deg * spa::kRad;
    return detail::factors(row, normal_enu, apparent_elevation_deg, true_azimuth_deg, std::sin(a), std::cos(a),
                           limb_deg);
}

// The limb (0 for the sun's centre) and the atmosphere refraction is
// computed for; a pressure of 0 is no refraction.
struct SunSettings {
    double limb_deg = 0.0;
    double pressure_hpa = 1010.0;
    double temperature_c = 10.0;
};

// Non-owning views: per receiver a geodetic latitude and longitude, a unit
// normal (east, north, up) and n_directions horizon entries, [i N + k].
struct LitReceivers {
    std::span<const double> lat_deg, lon_deg;
    std::span<const double> normal;
    std::span<const std::int16_t> horizon;
    std::size_t n_directions;
};

// T n values each, time-major: [t n + i].
struct InstantFactors {
    std::vector<std::uint8_t> lit;
    std::vector<float> cos_i, b, sin_e;
};

namespace detail {

inline void validate(const LitReceivers& r) {
    const std::size_t n = r.lat_deg.size();
    if (r.n_directions == 0)
        throw std::invalid_argument("n_directions must be at least 1");
    if (r.lon_deg.size() != n)
        throw std::invalid_argument("lon_deg needs one value per receiver, as lat_deg");
    if (r.normal.size() != 3 * n)
        throw std::invalid_argument("normal needs three values (east, north, up) per receiver");
    if (r.horizon.size() != n * r.n_directions)
        throw std::invalid_argument("horizon needs n_directions values per receiver");
    for (std::size_t i = 0; i < n; ++i) {
        const auto v = r.normal.subspan(3 * i, 3);
        const bool finite = std::isfinite(r.lat_deg[i]) && std::isfinite(r.lon_deg[i]) &&
                            std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
        if (!finite)
            throw std::invalid_argument("receiver " + std::to_string(i) + " has a number that is not finite");
        if (std::abs(r.lat_deg[i]) > kMaxLatitudeDeg)
            throw std::invalid_argument("receiver " + std::to_string(i) +
                                        " is at a latitude beyond 89.9 degrees, where east and north are undefined");
    }
}

inline double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) noexcept {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

}  // namespace detail

// The factors for every receiver at every time, receivers over `threads`
// threads (0: the hardware's count); the result does not depend on the count.
// Refusals (std::invalid_argument) name the field or the receiver.
inline InstantFactors evaluate_instant(const LitReceivers& r, std::span<const std::int64_t> times_utc_ms,
                                       const SunSettings& s, unsigned threads) {
    detail::validate(r);
    const std::size_t n = r.lat_deg.size(), nt = times_utc_ms.size(), N = r.n_directions;
    // The sun in Earth-fixed axes (x toward 0 N 0 E, z toward the pole), with
    // G = nu - alpha the Greenwich hour angle.
    std::vector<std::array<double, 3>> sun(nt);
    for (std::size_t t = 0; t < nt; ++t) {
        const utc_ms at{std::chrono::milliseconds{times_utc_ms[t]}};
        const GeocentricSun g = geocentric_sun(at, delta_t_seconds(at));
        const double d = g.delta_deg * spa::kRad, gh = (g.nu_deg - g.alpha_deg) * spa::kRad;
        sun[t] = {std::cos(d) * std::cos(gh), -std::cos(d) * std::sin(gh), std::sin(d)};
    }
    InstantFactors out{std::vector<std::uint8_t>(nt * n), std::vector<float>(nt * n), std::vector<float>(nt * n),
                       std::vector<float>(nt * n)};
    parallel_util::for_each_block(n, threads, {}, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            const EnuBasis basis = enu_basis(r.lat_deg[i], r.lon_deg[i]);
            const std::span<const double, 3> normal{r.normal.subspan(3 * i, 3)};
            const auto row = r.horizon.subspan(i * N, N);
            for (std::size_t t = 0; t < nt; ++t) {
                const double se = detail::dot(sun[t], basis.east), sn = detail::dot(sun[t], basis.north);
                const double e0 = std::asin(std::clamp(detail::dot(sun[t], basis.up), -1.0, 1.0)) / spa::kRad;
                const double h = std::hypot(se, sn);
                const double sin_a = h > 0.0 ? se / h : 0.0, cos_a = h > 0.0 ? sn / h : 1.0;
                const double a = spa::wrap360(std::atan2(se, sn) / spa::kRad);
                const SampleFactors f = detail::factors(row, normal, apparent_elevation(e0, s.pressure_hpa, s.temperature_c),
                                                        a, sin_a, cos_a, s.limb_deg);
                const std::size_t o = t * n + i;
                out.lit[o] = f.lit ? 1 : 0;
                out.cos_i[o] = static_cast<float>(f.cos_i);
                out.b[o] = static_cast<float>(f.b);
                out.sin_e[o] = static_cast<float>(f.sin_e);
            }
        }
    });
    return out;
}

}  // namespace terrain::radiation
