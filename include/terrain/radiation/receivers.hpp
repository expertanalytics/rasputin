// The receivers of terrain radiation and the local frame at each: the core's
// view of numbers Python computed. Each receiver has a position (x, y, z_obs)
// in the computation frame, a geodetic latitude and longitude, a 2x2 matrix J
// taking a frame displacement (dx, dy) to local true (east, north) metres, a
// unit normal in true (east, north, up), and an exclusion radius in metres.
// The core never learns which CRS the numbers came from.

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>

namespace terrain::radiation {

// Non-owning structure-of-arrays views; jacobian holds J00 J01 J10 J11 per
// receiver, normal E N U per receiver.
struct ReceiverSet {
    std::span<const double> x, y, z_obs;
    std::span<const double> lat_deg, lon_deg;
    std::span<const double> jacobian;
    std::span<const double> normal;
    std::span<const double> exclusion_radius;

    std::size_t size() const noexcept { return x.size(); }
};

enum class ReceiverStatus { ok, size_mismatch, non_finite, polar, negative_exclusion_radius };

struct ReceiverCheck {
    ReceiverStatus status;
    std::size_t index;  // the first failing receiver; 0 when ok or on a size mismatch
};

// East and north are undefined at the poles; beyond this latitude a receiver
// is refused (Svalbard, at most 81 degrees, is well inside).
inline constexpr double kMaxLatitudeDeg = 89.9;

// The first failing receiver. Within one receiver a non-finite number is
// reported before the pole, and the pole before a negative exclusion radius.
inline ReceiverCheck check_receivers(const ReceiverSet& r) noexcept {
    const std::size_t n = r.size();
    if (r.y.size() != n || r.z_obs.size() != n || r.lat_deg.size() != n || r.lon_deg.size() != n ||
        r.jacobian.size() != 4 * n || r.normal.size() != 3 * n || r.exclusion_radius.size() != n)
        return {ReceiverStatus::size_mismatch, 0};
    const auto finite = [](std::span<const double> s) {
        for (const double v : s)
            if (!std::isfinite(v))
                return false;
        return true;
    };
    for (std::size_t i = 0; i < n; ++i) {
        const double scalars[] = {r.x[i], r.y[i], r.z_obs[i], r.lat_deg[i], r.lon_deg[i], r.exclusion_radius[i]};
        if (!finite(scalars) || !finite(r.jacobian.subspan(4 * i, 4)) || !finite(r.normal.subspan(3 * i, 3)))
            return {ReceiverStatus::non_finite, i};
        if (std::abs(r.lat_deg[i]) > kMaxLatitudeDeg)
            return {ReceiverStatus::polar, i};
        if (r.exclusion_radius[i] < 0.0)
            return {ReceiverStatus::negative_exclusion_radius, i};
    }
    return {ReceiverStatus::ok, 0};
}

// check_receivers as a refusal: std::invalid_argument naming the first failing
// receiver, in the words the binding gives Python.
inline void require_receivers(const ReceiverSet& r) {
    const ReceiverCheck c = check_receivers(r);
    const std::string at = "receiver " + std::to_string(c.index);
    switch (c.status) {
    case ReceiverStatus::ok:
        return;
    case ReceiverStatus::size_mismatch:
        throw std::invalid_argument("the receiver arrays differ in length: x, y, z_obs, lat_deg, lon_deg and "
                                    "exclusion_radius need n values, jacobian 4 n and normal 3 n");
    case ReceiverStatus::non_finite:
        throw std::invalid_argument(at + " has a number that is not finite");
    case ReceiverStatus::polar:
        throw std::invalid_argument(at + " is at a latitude beyond 89.9 degrees, where east and north are undefined");
    case ReceiverStatus::negative_exclusion_radius:
        throw std::invalid_argument(at + " has a negative exclusion radius");
    }
}

// Unit east, north and up at a geodetic latitude and longitude, in Earth-fixed
// axes: x toward (0 N, 0 E), y toward (0 N, 90 E), z toward the north pole.
// Up is the ellipsoid normal, so it depends on the geodetic latitude alone.
struct EnuBasis {
    std::array<double, 3> east, north, up;
};

inline EnuBasis enu_basis(double lat_deg, double lon_deg) noexcept {
    constexpr double rad = std::numbers::pi / 180.0;
    const double sp = std::sin(lat_deg * rad), cp = std::cos(lat_deg * rad);
    const double sl = std::sin(lon_deg * rad), cl = std::cos(lon_deg * rad);
    return {{-sl, cl, 0.0}, {-sp * cl, -sp * sl, cp}, {cp * cl, cp * sl, sp}};
}

// The frame displacement (dx, dy) of one ground metre toward a true azimuth
// (degrees clockwise from true north): the solution of J (dx, dy) =
// (sin A, cos A). J is a receiver's four entries, row by row.
inline std::array<double, 2> frame_step(std::span<const double, 4> j, double true_azimuth_deg) noexcept {
    constexpr double rad = std::numbers::pi / 180.0;
    const double e = std::sin(true_azimuth_deg * rad), n = std::cos(true_azimuth_deg * rad);
    const double det = j[0] * j[3] - j[1] * j[2];
    return {(j[3] * e - j[1] * n) / det, (j[0] * n - j[2] * e) / det};
}

}  // namespace terrain::radiation
