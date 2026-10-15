// Increment 46 PR 6, the terrain-radiation design's TR5: the `instant`
// evaluator (docs/increments/46-hbv-snow.md, section 10.8 "PR 6", tests L1,
// L2, I7, I8 and S; docs/increments/terrain-radiation.md, R5, I7, I8, its
// degeneracy policy and verification 3's flat plane).
//
//   H1-H3  the horizon at a true azimuth: the stored value on a direction,
//          linear between neighbours, wrapping past 360 degrees;
//   L1     the lit test on the V-valley's stored row: strict at 30.00, the
//          interpolated 26.5333 at azimuth 120, self-shading, the limb;
//   L2     the beam factor and the sine of the elevation at 20 and 40 degrees;
//   I7     the evaluator against the reference path (solar_position) at
//          Bygdin's centre, Solund and Vardo, every 15 minutes of 2027-04-15,
//          05-15 and 06-21 with the sun above -1 degree: apparent elevation
//          and true azimuth within 0.003 degrees; lit on the V-valley row
//          agreeing with the reference sun's lit away from the horizon;
//   I8     every output equal bit for bit at 1 and 8 threads;
//   S      the flat plane's sunrise: a vertical east-facing receiver, horizon
//          0 in every direction, the upper limb; unlit 4 s before the
//          reference path's crossing, lit 4 s after, on 04-15 and 06-21;
//   V      refusals.
//
// Interface under test (header include/terrain/radiation/evaluate.hpp,
// namespace terrain::radiation):
//
//   inline constexpr double kUpperLimbDeg = 0.2667;
//   // Degrees: entry / 100.0 on a stored direction, linear between the two
//   // neighbouring directions otherwise (direction N - 1 next to 0); row has
//   // N entries, direction k at true azimuth k 360 / N.
//   double horizon_at(std::span<const std::int16_t> row, double true_azimuth_deg) noexcept;
//   struct SampleFactors { bool lit; double cos_i, b, sin_e; };
//   // One receiver and one apparent sun: lit = e_app + limb > horizon_at(A)
//   // and n.s_app > 0; cos_i = n.s_app, unclipped; b = lit max(0, cos_i);
//   // sin_e = sin e_app, unclipped. s_app = (cos e sin A, cos e cos A, sin e).
//   SampleFactors sample_factors(std::span<const std::int16_t> row, std::span<const double, 3> normal_enu,
//                                double apparent_elevation_deg, double true_azimuth_deg,
//                                double limb_deg) noexcept;
//   struct SunSettings { double limb_deg = 0.0; double pressure_hpa = 1010.0; double temperature_c = 10.0; };
//   struct LitReceivers {
//       std::span<const double> lat_deg, lon_deg;
//       std::span<const double> normal;         // E N U per receiver
//       std::span<const std::int16_t> horizon;  // n N, [i N + k], 0.01 degrees
//       std::size_t n_directions;
//   };
//   struct InstantFactors {                     // T n each, time-major: [t n + i]
//       std::vector<std::uint8_t> lit;
//       std::vector<float> cos_i, b, sin_e;
//   };
//   // The sun per time from geocentric_sun at delta_t_seconds(t), seen along
//   // each receiver's ENU basis (no parallax), refracted by
//   // apparent_elevation at the settings' pressure and temperature.
//   InstantFactors evaluate_instant(const LitReceivers&, std::span<const std::int64_t> times_utc_ms,
//                                   const SunSettings&, unsigned threads);
//
// evaluate_instant throws std::invalid_argument naming "n_directions" for 0,
// "horizon", "normal" or "lon_deg" for a length that does not match lat_deg,
// and "receiver i" for a latitude or longitude that is not finite or beyond
// 89.9 degrees, or a normal that is not finite.
//
// Absolute bounds: angles in degrees; 0.003 is the design's I7 bound (the
// parallax dropped is at most 0.0024). The float32 outputs carry about 1e-7
// relative, about 1e-5 degrees after asin or atan2 at the largest elevation
// met here (52 degrees); checked at Vardo, the furthest north (70.37 N).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain/radiation/evaluate.hpp>
#include <terrain/radiation/solar_position.hpp>
#include <terrain/radiation/time.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace terrain::radiation;
using Catch::Matchers::WithinAbs;

constexpr double kDeg = std::numbers::pi / 180.0;
constexpr std::size_t kN = 64;
constexpr double kP = 1010.0, kT = 10.0;  // the standard atmosphere of SPA's refraction
constexpr std::array<double, 3> kUp{0.0, 0.0, 1.0}, kEast{1.0, 0.0, 0.0}, kNorth{0.0, 1.0, 0.0};
constexpr std::int64_t kMinute = 60'000;

// The V-valley's row as the radiation design stores it: tan h = |sin A| tan 30.
std::vector<std::int16_t> valley_row() {
    std::vector<std::int16_t> row(kN);
    for (std::size_t k = 0; k < kN; ++k) {
        const double a = static_cast<double>(k) * 360.0 / kN;
        row[k] = static_cast<std::int16_t>(std::lround(100.0 * std::atan(std::abs(std::sin(a * kDeg)) * std::tan(30.0 * kDeg)) / kDeg));
    }
    return row;
}

std::vector<std::int16_t> constant_row(std::int16_t v) { return std::vector<std::int16_t>(kN, v); }

std::int64_t utc_day_ms(std::chrono::year_month_day d) {
    return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::sys_days{d}.time_since_epoch()).count();
}

utc_ms at(std::int64_t ms) { return utc_ms{std::chrono::milliseconds{ms}}; }

SolarPosition reference(std::int64_t t, double lat, double lon, double h) {
    return solar_position(at(t), lat, lon, h, delta_t_seconds(at(t)), kP, kT);
}

double azimuth_gap(double a, double b) {
    const double d = std::fmod(std::abs(a - b), 360.0);
    return std::min(d, 360.0 - d);
}

// Receivers that own their arrays.
struct Receivers {
    std::vector<double> lat, lon, normal;
    std::vector<std::int16_t> horizon;
    void add(double la, double lo, std::array<double, 3> n, const std::vector<std::int16_t>& row) {
        lat.push_back(la);
        lon.push_back(lo);
        normal.insert(normal.end(), n.begin(), n.end());
        horizon.insert(horizon.end(), row.begin(), row.end());
    }
    LitReceivers set() const { return LitReceivers{lat, lon, normal, horizon, kN}; }
    std::size_t size() const { return lat.size(); }
};

struct Site {
    const char* name;
    double lat, lon, height;
};
constexpr std::array<Site, 3> kSites{{
    {"Bygdin", 61.3694, 8.4852, 1500.0},
    {"Solund", 61.08, 4.85, 0.0},
    {"Vardo", 70.37, 31.11, 0.0},
}};

std::vector<std::chrono::year_month_day> days() {
    using namespace std::chrono;
    return {2027y / April / 15, 2027y / May / 15, 2027y / June / 21};
}

}  // namespace

TEST_CASE("H1-H3 the horizon at a true azimuth", "[radiation][evaluate]") {
    const auto row = valley_row();
    REQUIRE(row[16] == 3000);
    REQUIRE(row[21] == 2698);
    REQUIRE(row[22] == 2564);

    SECTION("H1 on a stored direction, the stored value exactly") {
        CHECK(horizon_at(row, 90.0) == 30.0);
        CHECK(horizon_at(row, 0.0) == 0.0);
        CHECK(horizon_at(row, 118.125) == 26.98);
    }
    SECTION("H2 between two directions, linear: 26.5333 at 120 degrees") {
        CHECK_THAT(horizon_at(row, 120.0), WithinAbs(26.98 + (25.64 - 26.98) / 3.0, 1e-9));
    }
    SECTION("H3 past the last direction, toward direction 0") {
        std::vector<std::int16_t> ramp(kN, 0);
        ramp[kN - 1] = 400;  // 4 degrees at 354.375; 0 at 0
        CHECK_THAT(horizon_at(ramp, 357.1875), WithinAbs(2.0, 1e-9));
        CHECK_THAT(horizon_at(ramp, 2.8125), WithinAbs(0.0, 1e-12));  // between 0 and 1, both 0
    }
}

TEST_CASE("L1 the lit test on the V-valley floor", "[radiation][evaluate]") {
    const auto row = valley_row();
    const auto lit = [&](double e, double a, double limb = 0.0) { return sample_factors(row, kUp, e, a, limb).lit; };

    SECTION("sun at azimuth 90: lit is strict at the stored 30.00") {
        CHECK_FALSE(lit(29.99, 90.0));
        CHECK(lit(30.01, 90.0));
        CHECK_FALSE(lit(30.0, 90.0));
    }
    SECTION("sun at azimuth 120: the interpolated 26.5333") {
        CHECK_FALSE(lit(26.52, 120.0));
        CHECK(lit(26.55, 120.0));
    }
    SECTION("the upper limb raises the sun by 0.2667") {
        CHECK_FALSE(lit(29.8, 90.0));
        CHECK(lit(29.8, 90.0, kUpperLimbDeg));
    }
    SECTION("self-shading: a 30 degree slope facing north, open sky, sun south at 20 degrees") {
        const std::array<double, 3> north_slope{0.0, std::sin(30.0 * kDeg), std::cos(30.0 * kDeg)};
        const SampleFactors f = sample_factors(constant_row(-9000), north_slope, 20.0, 180.0, 0.0);
        CHECK_FALSE(f.lit);
        CHECK(f.b == 0.0);
        // cos_i unclipped: -0.5 cos 20 + cos 30 sin 20.
        CHECK_THAT(f.cos_i, WithinAbs(-0.5 * std::cos(20.0 * kDeg) + std::cos(30.0 * kDeg) * std::sin(20.0 * kDeg), 1e-12));
    }
    SECTION("a lit sloping receiver: b is cos_i") {
        const std::array<double, 3> south_slope{0.0, -std::sin(30.0 * kDeg), std::cos(30.0 * kDeg)};
        const SampleFactors f = sample_factors(constant_row(-9000), south_slope, 30.0, 180.0, 0.0);
        CHECK(f.lit);
        CHECK(f.b == f.cos_i);
        // Normal (0, -sin 30, cos 30), sun (0, -cos 30, sin 30): cos_i = 2 sin 30 cos 30 = 0.8660.
        CHECK_THAT(f.cos_i, WithinAbs(2.0 * std::sin(30.0 * kDeg) * std::cos(30.0 * kDeg), 1e-12));
    }
}

TEST_CASE("L2 the beam factor and the sine of the elevation on the valley floor", "[radiation][evaluate]") {
    const auto row = valley_row();
    const SampleFactors low = sample_factors(row, kUp, 20.0, 90.0, 0.0);
    CHECK_FALSE(low.lit);
    CHECK(low.b == 0.0);
    CHECK_THAT(low.sin_e, WithinAbs(0.3420, 5e-5));
    CHECK_THAT(800.0 * (low.b - std::max(0.0, low.sin_e)), WithinAbs(-273.6, 0.05));

    const SampleFactors high = sample_factors(row, kUp, 40.0, 90.0, 0.0);
    CHECK(high.lit);
    CHECK_THAT(high.b, WithinAbs(0.6428, 5e-5));
    CHECK(high.b == high.sin_e);  // so the combiner's anomaly is exactly 0

    SECTION("open sky, a flat receiver: b is max(0, sin_e) at every elevation") {
        for (double e = -5.0; e <= 89.0; e += 0.5) {
            const SampleFactors f = sample_factors(constant_row(-9000), kUp, e, 123.0, 0.0);
            INFO("e " << e);
            CHECK(f.b == std::max(0.0, f.sin_e));
        }
    }
}

TEST_CASE("I7 the evaluator against the reference sun", "[radiation][evaluate]") {
    const auto open = constant_row(-9000);
    const auto valley = valley_row();
    for (const Site& site : kSites) {
        // Up, a vertical wall facing east and one facing north: asin(sin_e) is
        // the apparent elevation, atan2(cos_i east, cos_i north) the azimuth.
        Receivers r;
        r.add(site.lat, site.lon, kUp, open);
        r.add(site.lat, site.lon, kEast, open);
        r.add(site.lat, site.lon, kNorth, open);
        r.add(site.lat, site.lon, kUp, valley);
        for (const auto& day : days()) {
            std::vector<std::int64_t> times;
            for (std::int64_t q = 0; q < 96; ++q) {
                const std::int64_t t = utc_day_ms(day) + q * 15 * kMinute;
                if (reference(t, site.lat, site.lon, site.height).e_deg > -1.0)
                    times.push_back(t);
            }
            REQUIRE(times.size() > 40);
            const InstantFactors f = evaluate_instant(r.set(), times, SunSettings{0.0, kP, kT}, 1);
            REQUIRE(f.lit.size() == times.size() * r.size());
            REQUIRE(f.cos_i.size() == f.lit.size());
            REQUIRE(f.b.size() == f.lit.size());
            REQUIRE(f.sin_e.size() == f.lit.size());
            std::size_t lit_count = 0, shaded_by_valley = 0;
            for (std::size_t s = 0; s < times.size(); ++s) {
                const SolarPosition ref = reference(times[s], site.lat, site.lon, site.height);
                const std::size_t o = s * r.size();
                const double e = std::asin(static_cast<double>(f.sin_e[o])) / kDeg;
                const double a = std::fmod(std::atan2(static_cast<double>(f.cos_i[o + 1]),
                                                      static_cast<double>(f.cos_i[o + 2])) / kDeg + 360.0, 360.0);
                INFO(site.name << " sample " << s << " ref e " << ref.e_deg << " A " << ref.azimuth_deg);
                CHECK_THAT(e, WithinAbs(ref.e_deg, 0.003));
                CHECK(azimuth_gap(a, ref.azimuth_deg) <= 0.003);
                // The sine of the elevation is the up receiver's cos_i.
                CHECK(f.sin_e[o] == f.cos_i[o]);

                const double h = horizon_at(valley, ref.azimuth_deg);
                const bool ref_lit = sample_factors(valley, kUp, ref.e_deg, ref.azimuth_deg, 0.0).lit;
                if (std::abs(ref.e_deg - h) > 0.003) {
                    CHECK((f.lit[o + 3] != 0) == ref_lit);
                }
                lit_count += f.lit[o + 3];
                shaded_by_valley += (ref.e_deg > 0.0 && !ref_lit);
            }
            // The probe can fail: the valley both lights and shades the floor.
            CHECK(lit_count > 0);
            CHECK(shaded_by_valley > 0);
        }
    }
}

TEST_CASE("I8 every output is the same at 1 and 8 threads", "[radiation][evaluate]") {
    // 997 receivers (not a multiple of any chunk) on a spread of latitudes,
    // longitudes, normals and horizon rows; 96 samples of 2027-04-15.
    Receivers r;
    for (std::size_t i = 0; i < 997; ++i) {
        const double u = static_cast<double>(i);
        const double slope = std::fmod(u * 7.3, 80.0) * kDeg, aspect = std::fmod(u * 37.1, 360.0) * kDeg;
        std::vector<std::int16_t> row(kN);
        for (std::size_t k = 0; k < kN; ++k)
            row[k] = static_cast<std::int16_t>(std::lround(1500.0 * std::sin(u * 0.37 + static_cast<double>(k) * 0.29)));
        r.add(58.0 + std::fmod(u * 0.013, 13.0), 4.0 + std::fmod(u * 0.029, 27.0),
              {std::sin(slope) * std::sin(aspect), std::sin(slope) * std::cos(aspect), std::cos(slope)}, row);
    }
    std::vector<std::int64_t> times;
    for (std::int64_t q = 0; q < 96; ++q)
        times.push_back(utc_day_ms(std::chrono::year_month_day{std::chrono::year{2027}, std::chrono::April, std::chrono::day{15}}) + q * 15 * kMinute);

    const InstantFactors one = evaluate_instant(r.set(), times, SunSettings{}, 1);
    const InstantFactors eight = evaluate_instant(r.set(), times, SunSettings{}, 8);
    CHECK(one.lit == eight.lit);
    CHECK(one.cos_i == eight.cos_i);
    CHECK(one.b == eight.b);
    CHECK(one.sin_e == eight.sin_e);
    std::size_t lit = 0, unlit_above = 0;
    for (std::size_t j = 0; j < one.lit.size(); ++j) {
        lit += one.lit[j];
        unlit_above += (!one.lit[j] && one.sin_e[j] > 0.2f);
    }
    CHECK(lit > 0);
    CHECK(unlit_above > 0);
}

TEST_CASE("S the flat plane's sunrise, upper limb, vertical east-facing receiver", "[radiation][evaluate]") {
    using namespace std::chrono;
    const Site b = kSites[0];
    Receivers r;
    r.add(b.lat, b.lon, kEast, constant_row(0));  // -atan(c r_i) = -0.0000391 degrees, stored as 0
    struct Case {
        year_month_day day;
        std::int64_t noaa_ms, lo_ms, hi_ms;  // NOAA's sunrise (probe S) and a bracket
    };
    const std::array<Case, 2> cases{{
        {2027y / April / 15, (4 * 60 + 6) * kMinute, 3 * 60 * kMinute, 5 * 60 * kMinute},
        {2027y / June / 21, (1 * 60 + 44) * kMinute + 30'000, 30 * kMinute, 3 * 60 * kMinute},
    }};
    for (const Case& c : cases) {
        const std::int64_t d0 = utc_day_ms(c.day);
        const auto upper = [&](std::int64_t t) { return reference(t, b.lat, b.lon, b.height).e_deg + kUpperLimbDeg; };
        std::int64_t lo = d0 + c.lo_ms, hi = d0 + c.hi_ms;
        REQUIRE(upper(lo) < 0.0);
        REQUIRE(upper(hi) > 0.0);
        while (hi - lo > 1) {
            const std::int64_t mid = lo + (hi - lo) / 2;
            (upper(mid) > 0.0 ? hi : lo) = mid;
        }
        INFO("crossing at " << (hi - d0) / 1000 << " s after midnight UTC");
        // NOAA's 0.833 degree convention against SPA's refraction at 1010 hPa
        // and 10 C differs by about 0.05 degrees: under a minute either day.
        CHECK(std::abs(hi - (d0 + c.noaa_ms)) <= 90'000);
        const std::vector<std::int64_t> times{hi - 4000, hi + 4000};
        const InstantFactors f = evaluate_instant(r.set(), times, SunSettings{kUpperLimbDeg, kP, kT}, 1);
        CHECK(f.lit[0] == 0);
        CHECK(f.lit[1] == 1);
        CHECK(f.cos_i[1] > 0.0f);  // the east wall faces the rising sun
    }
}

TEST_CASE("No times and no receivers give empty arrays", "[radiation][evaluate]") {
    Receivers r;
    r.add(61.0, 8.0, kUp, constant_row(0));
    const InstantFactors none = evaluate_instant(r.set(), {}, SunSettings{}, 1);
    CHECK(none.lit.empty());
    CHECK(none.b.empty());
    const Receivers empty;
    const std::vector<std::int64_t> t{0};
    const InstantFactors nobody = evaluate_instant(empty.set(), t, SunSettings{}, 1);
    CHECK(nobody.lit.empty());
    CHECK(nobody.sin_e.empty());
}

TEST_CASE("V refusals name the field or the receiver", "[radiation][evaluate]") {
    const std::vector<std::int64_t> t{0};
    const auto refused = [&](const LitReceivers& set, const std::string& says) {
        try {
            (void)evaluate_instant(set, t, SunSettings{}, 1);
        } catch (const std::invalid_argument& e) {
            INFO(e.what());
            CHECK(std::string(e.what()).find(says) != std::string::npos);
            return;
        }
        FAIL("not refused: " << says);
    };
    Receivers r;
    r.add(61.0, 8.0, kUp, constant_row(0));
    r.add(62.0, 9.0, kUp, constant_row(0));

    LitReceivers zero = r.set();
    zero.n_directions = 0;
    refused(zero, "n_directions");

    LitReceivers short_h = r.set();
    short_h.horizon = short_h.horizon.first(kN);
    refused(short_h, "horizon");

    LitReceivers short_n = r.set();
    short_n.normal = short_n.normal.first(3);
    refused(short_n, "normal");

    LitReceivers short_lon = r.set();
    short_lon.lon_deg = short_lon.lon_deg.first(1);
    refused(short_lon, "lon_deg");

    Receivers polar = r;
    polar.lat[1] = 89.95;
    refused(polar.set(), "receiver 1");

    Receivers nan = r;
    nan.normal[4] = std::numeric_limits<double>::quiet_NaN();
    refused(nan.set(), "receiver 1");

    Receivers inf = r;
    inf.lon[0] = std::numeric_limits<double>::infinity();
    refused(inf.set(), "receiver 0");
}
