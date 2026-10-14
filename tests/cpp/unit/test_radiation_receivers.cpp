// Increment 46 PR 4, the terrain-radiation design's TR2: the receiver frame in
// the core (docs/increments/46-hbv-snow.md, section 10.8 "PR 4, TR2";
// docs/increments/terrain-radiation.md, R3, R10, I3 and its degeneracy
// policy).
//
//   E1-E3  the east-north-up basis at a receiver's geodetic latitude and
//          longitude: orthonormal, right-handed, up the ellipsoid normal;
//   F1-F3  frame_step, the frame displacement for one ground metre toward a
//          true azimuth, through J; the sign of the convergence pinned at
//          Solund and Vardo (true azimuth = grid azimuth + gamma, I3);
//   C1-C6  check_receivers on the SoA views: sizes, non-finite numbers, the
//          pole, a negative exclusion radius, the first failing receiver;
//          empty and single-receiver sets.
//
// Interface assumed (header include/terrain/radiation/receivers.hpp, namespace
// terrain::radiation; names chosen here and listed in the handback):
//
//   struct ReceiverSet {                       // non-owning SoA views
//       std::span<const double> x, y, z_obs;   // frame, n each
//       std::span<const double> lat_deg, lon_deg;
//       std::span<const double> jacobian;      // 4n: per receiver J00 J01 J10 J11,
//                                              // (east, north) = J (dx, dy)
//       std::span<const double> normal;        // 3n: E, N, U per receiver
//       std::span<const double> exclusion_radius;  // n, metres
//       std::size_t size() const noexcept;     // x.size()
//   };
//   enum class ReceiverStatus { ok, size_mismatch, non_finite, polar,
//                               negative_exclusion_radius };
//   struct ReceiverCheck { ReceiverStatus status; std::size_t index; };
//   ReceiverCheck check_receivers(const ReceiverSet&) noexcept;
//   struct EnuBasis { std::array<double, 3> east, north, up; };
//   EnuBasis enu_basis(double lat_deg, double lon_deg) noexcept;
//   std::array<double, 2> frame_step(std::span<const double, 4> j,
//                                    double true_azimuth_deg) noexcept;
//
// The basis is in Earth-fixed Cartesian axes: x toward (0 N, 0 E), y toward
// (0 N, 90 E), z toward the north pole. The J values below are pyproj 3.8.0
// (PROJ 9.8.1) get_factors for EPSG:25833, J = (1/m) [[cos g, sin g],
// [-sin g, cos g]]; tests/python/test_radiation_frame.py checks that form
// against a finite-difference transform.
//
// Absolute bounds: unit vectors and J entries of order 1; 1e-12 is about
// 4500 ulp at 1, checked at latitudes up to 89.9 degrees.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain/radiation/receivers.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <vector>

namespace {

using namespace terrain::radiation;
using Catch::Matchers::WithinAbs;

constexpr double kDeg = std::numbers::pi / 180.0;
constexpr double kTol = 1e-12;

using Vec3 = std::array<double, 3>;

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// The point at geodetic (lat, lon, h) on GRS80, Earth-fixed metres.
Vec3 ecef(double lat_deg, double lon_deg, double h) {
    constexpr double a = 6378137.0;
    constexpr double f = 1.0 / 298.257222101;
    constexpr double e2 = f * (2.0 - f);
    const double p = lat_deg * kDeg, l = lon_deg * kDeg;
    const double n = a / std::sqrt(1.0 - e2 * std::sin(p) * std::sin(p));
    return {(n + h) * std::cos(p) * std::cos(l), (n + h) * std::cos(p) * std::sin(l),
            (n * (1.0 - e2) + h) * std::sin(p)};
}

// J = (1/m) Rot(gamma) for a projected frame: grid north (0, 1) has true
// azimuth gamma.
std::array<double, 4> rot_j(double gamma_deg, double m) {
    const double c = std::cos(gamma_deg * kDeg), s = std::sin(gamma_deg * kDeg);
    return {c / m, s / m, -s / m, c / m};
}

// pyproj's meridian_convergence and point scale in EPSG:25833 (probe in the
// commit message).
constexpr double kSolundGamma = -8.906046, kSolundM = 1.00325523;  // 4.85 E, 61.08 N
constexpr double kVardoGamma = 15.218465, kVardoM = 1.00397453;    // 31.11 E, 70.37 N

double azimuth_of(double east, double north) {
    double a = std::atan2(east, north) / kDeg;
    return a < 0.0 ? a + 360.0 : a;
}

// One valid receiver's worth of arrays, n receivers long; tests break one.
struct Arrays {
    std::vector<double> x, y, z, lat, lon, j, normal, r;
    explicit Arrays(std::size_t n)
        : x(n, 500000.0), y(n, 6800000.0), z(n, 1000.0), lat(n, 61.0), lon(n, 8.5), j(4 * n, 0.0),
          normal(3 * n, 0.0), r(n, 10.0) {
        for (std::size_t i = 0; i < n; ++i) {
            j[4 * i] = j[4 * i + 3] = 1.0;
            normal[3 * i + 2] = 1.0;
        }
    }
    ReceiverSet view() const {
        return ReceiverSet{.x = x, .y = y, .z_obs = z, .lat_deg = lat, .lon_deg = lon,
                           .jacobian = j, .normal = normal, .exclusion_radius = r};
    }
};

}  // namespace

// --- E: the east-north-up basis -------------------------------------------

TEST_CASE("E1: at (0 N, 0 E) the basis is the Earth-fixed axes", "[radiation][receivers][enu]") {
    const EnuBasis b = enu_basis(0.0, 0.0);
    const Vec3 east{0, 1, 0}, north{0, 0, 1}, up{1, 0, 0};
    for (std::size_t k = 0; k < 3; ++k) {
        CHECK_THAT(b.east[k], WithinAbs(east[k], kTol));
        CHECK_THAT(b.north[k], WithinAbs(north[k], kTol));
        CHECK_THAT(b.up[k], WithinAbs(up[k], kTol));
    }
}

TEST_CASE("E2: the basis is orthonormal and right-handed (east x north = up)",
          "[radiation][receivers][enu]") {
    const std::array<std::array<double, 2>, 7> sites{{{61.3694, 8.4852},  // Bygdin
                                                      {70.37, 31.11},     // Vardo
                                                      {-33.9, 18.4},
                                                      {0.0, -179.9},
                                                      {89.9, 45.0},       // the largest latitude taken
                                                      {-89.9, -120.0},
                                                      {45.0, 180.0}}};
    for (const auto& [lat, lon] : sites) {
        CAPTURE(lat, lon);
        const EnuBasis b = enu_basis(lat, lon);
        CHECK_THAT(dot(b.east, b.east), WithinAbs(1.0, kTol));
        CHECK_THAT(dot(b.north, b.north), WithinAbs(1.0, kTol));
        CHECK_THAT(dot(b.up, b.up), WithinAbs(1.0, kTol));
        CHECK_THAT(dot(b.east, b.north), WithinAbs(0.0, kTol));
        CHECK_THAT(dot(b.east, b.up), WithinAbs(0.0, kTol));
        CHECK_THAT(dot(b.north, b.up), WithinAbs(0.0, kTol));
        const Vec3 u = cross(b.east, b.north);
        for (std::size_t k = 0; k < 3; ++k)
            CHECK_THAT(u[k], WithinAbs(b.up[k], kTol));
        CHECK_THAT(b.east[2], WithinAbs(0.0, kTol));  // east is horizontal and has no polar part
        CHECK(b.north[2] >= 0.0);                     // north leans toward the north pole
    }
}

TEST_CASE("E3: up is the ellipsoid normal (geodetic, not geocentric) at Bygdin",
          "[radiation][receivers][enu]") {
    // Bygdin's centre: geodetic and geocentric latitude differ by about 0.17
    // degrees there, so up's polar part separates the two by 2.6e-3.
    const double lat = 61.3694, lon = 8.4852;
    const EnuBasis b = enu_basis(lat, lon);
    CHECK_THAT(b.up[2], WithinAbs(std::sin(lat * kDeg), kTol));
    // One metre up the ellipsoid normal, in Earth-fixed metres; the cancellation
    // in the difference costs about 1e-9 at 6.4e6 m.
    const Vec3 p0 = ecef(lat, lon, 1500.0), p1 = ecef(lat, lon, 1501.0);
    for (std::size_t k = 0; k < 3; ++k)
        CHECK_THAT(p1[k] - p0[k], WithinAbs(b.up[k], 1e-8));
    // North points along the meridian toward rising latitude.
    const Vec3 q0 = ecef(lat - 1e-6, lon, 1500.0), q1 = ecef(lat + 1e-6, lon, 1500.0);
    const Vec3 dq{q1[0] - q0[0], q1[1] - q0[1], q1[2] - q0[2]};
    const double len = std::sqrt(dot(dq, dq));
    for (std::size_t k = 0; k < 3; ++k)
        CHECK_THAT(dq[k] / len, WithinAbs(b.north[k], 1e-7));
}

// --- F: frame_step, the sign of the convergence through J ------------------

TEST_CASE("F1: Solund, true north is grid azimuth +8.906 and grid north is true -8.906",
          "[radiation][receivers][frame_step]") {
    const auto j = rot_j(kSolundGamma, kSolundM);
    // One ground metre toward true north: grid azimuth = true - gamma.
    const auto n = frame_step(j, 0.0);
    CHECK_THAT(azimuth_of(n[0], n[1]), WithinAbs(8.906046, 1e-6));
    CHECK_THAT(std::hypot(n[0], n[1]), WithinAbs(kSolundM, 1e-12));  // m frame units per ground metre
    // One ground metre toward true azimuth gamma is pure grid north.
    const auto g = frame_step(j, kSolundGamma);
    CHECK_THAT(g[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(g[1], WithinAbs(kSolundM, 1e-12));
}

TEST_CASE("F2: Vardo, grid north is true +15.218", "[radiation][receivers][frame_step]") {
    const auto j = rot_j(kVardoGamma, kVardoM);
    const auto g = frame_step(j, 15.218465);
    CHECK_THAT(g[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(g[1], WithinAbs(kVardoM, 1e-12));
    const auto e = frame_step(j, 90.0);  // true east is grid azimuth 90 - 15.218
    CHECK_THAT(azimuth_of(e[0], e[1]), WithinAbs(90.0 - 15.218465, 1e-6));
}

TEST_CASE("F3: frame_step inverts J for every true azimuth, anisotropic J included",
          "[radiation][receivers][frame_step]") {
    // A lattice-frame J (diagonal, unequal) and a general one: J applied to the
    // step gives (sin A, cos A), one ground metre toward A.
    const std::array<std::array<double, 4>, 3> js{{rot_j(kSolundGamma, kSolundM),
                                                   {53.6, 0.0, 0.0, 111.4},
                                                   {0.9, 0.3, -0.2, 1.1}}};
    for (const auto& j : js) {
        for (int k = 0; k < 64; ++k) {
            const double a = k * 360.0 / 64.0;
            CAPTURE(j, a);
            const auto s = frame_step(j, a);
            CHECK_THAT(j[0] * s[0] + j[1] * s[1], WithinAbs(std::sin(a * kDeg), kTol));
            CHECK_THAT(j[2] * s[0] + j[3] * s[1], WithinAbs(std::cos(a * kDeg), kTol));
        }
    }
}

// --- C: check_receivers ------------------------------------------------------

TEST_CASE("C1: valid sets pass, empty and single-receiver included", "[radiation][receivers][check]") {
    for (std::size_t n : {0u, 1u, 5u}) {
        CAPTURE(n);
        const Arrays a(n);
        CHECK(a.view().size() == n);
        CHECK(check_receivers(a.view()).status == ReceiverStatus::ok);
    }
    Arrays v(1);  // a vertical receiver (n_U = 0) is allowed
    v.normal = {1.0, 0.0, 0.0};
    CHECK(check_receivers(v.view()).status == ReceiverStatus::ok);
    Arrays z(1);  // a zero exclusion radius is allowed
    z.r[0] = 0.0;
    CHECK(check_receivers(z.view()).status == ReceiverStatus::ok);
}

TEST_CASE("C2: every array of the wrong length is a size mismatch", "[radiation][receivers][check]") {
    for (int which = 0; which < 7; ++which) {
        CAPTURE(which);
        Arrays a(3);
        std::vector<double>* arr[] = {&a.y, &a.z, &a.lat, &a.lon, &a.j, &a.normal, &a.r};
        arr[which]->pop_back();
        CHECK(check_receivers(a.view()).status == ReceiverStatus::size_mismatch);
    }
}

TEST_CASE("C3: a NaN or infinite number names its receiver", "[radiation][receivers][check]") {
    const double bad[] = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
                          -std::numeric_limits<double>::infinity()};
    for (const double v : bad) {
        for (int which = 0; which < 8; ++which) {
            CAPTURE(v, which);
            Arrays a(4);
            std::vector<double>* arr[] = {&a.x, &a.y, &a.z, &a.lat, &a.lon, &a.j, &a.normal, &a.r};
            const std::size_t stride = which == 5 ? 4 : which == 6 ? 3 : 1;
            (*arr[which])[2 * stride + stride - 1] = v;  // receiver 2's last component
            const auto c = check_receivers(a.view());
            CHECK(c.status == ReceiverStatus::non_finite);
            CHECK(c.index == 2);
        }
    }
}

TEST_CASE("C4: |latitude| above 89.9 is refused, 89.9 itself is taken", "[radiation][receivers][check]") {
    for (const double lat : {89.9, -89.9}) {
        Arrays a(2);
        a.lat[1] = lat;
        CHECK(check_receivers(a.view()).status == ReceiverStatus::ok);
    }
    for (const double lat : {89.900001, -89.900001, 90.0, -90.0}) {
        CAPTURE(lat);
        Arrays a(2);
        a.lat[1] = lat;
        const auto c = check_receivers(a.view());
        CHECK(c.status == ReceiverStatus::polar);
        CHECK(c.index == 1);
    }
}

TEST_CASE("C5: a negative exclusion radius is refused", "[radiation][receivers][check]") {
    Arrays a(3);
    a.r[0] = -1e-9;
    const auto c = check_receivers(a.view());
    CHECK(c.status == ReceiverStatus::negative_exclusion_radius);
    CHECK(c.index == 0);
}

TEST_CASE("C6: the first failing receiver is reported; within one, non-finite before the pole",
          "[radiation][receivers][check]") {
    Arrays a(5);
    a.r[4] = -1.0;
    a.lat[3] = 90.5;
    a.lat[1] = 91.0;
    a.x[1] = std::numeric_limits<double>::quiet_NaN();
    const auto c = check_receivers(a.view());
    CHECK(c.status == ReceiverStatus::non_finite);
    CHECK(c.index == 1);
    a.x[1] = 0.0;
    const auto d = check_receivers(a.view());
    CHECK(d.status == ReceiverStatus::polar);
    CHECK(d.index == 1);
}
