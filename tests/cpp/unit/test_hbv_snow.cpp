// Increment 46 PR 1 (docs/increments/46-hbv-snow.md, sections 4.2 and 7):
// the HBV snow kernel, terrain::snow in include/terrain/snow/hbv_snow.hpp.
//
//   H1-H12  hand-computed single steps, each named in section 7's mutant table
//           as the fixture that kills one mutant of the per-bin update;
//   P1-P5   invariants I1-I5 over seeded random forcing (fixed seeds, so the
//           suite is deterministic): mass, bounds, flat is plain, N = 1 is
//           plain HBV-96, threads do not change a bit;
//   V1      every refusal in the home section 4.2 gives it, the state
//           compared bit for bit before and after a refused step;
//   V2      the distribution's renormalisation.
//
// The bins of a state are read through sp() and sw(): spans of n*N values,
// triangle-major (triangle i, bin k at i*N + k). A state with water in the pack
// (sw > 0) is made by one priming step of rain at T = tx = ts = 0, which adds
// the rain to sw, melts nothing and refreezes nothing (H10's path).
//
// Invariant-critical (section 7): mutation-tested. Starts threads (P5: 1 and 8).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/snow/hbv_snow.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using terrain::snow::HbvSnowParameters;
using terrain::snow::SnowDistribution;
using terrain::snow::SnowState;
using terrain::snow::step;
using terrain::snow::StepForcing;
using terrain::snow::StepStatus;
using terrain::snow::SunCorrection;
using terrain::snow::SunCorrectionParameters;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double kDay = 86400.0;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// Section 5's defaults: tx = ts = 0, cx = 3, lw = 0.1, cfr = 0.05.
constexpr HbvSnowParameters kHbv{.tx = 0.0, .ts = 0.0, .cx = 3.0, .lw = 0.1, .cfr = 0.05};

constexpr SunCorrectionParameters sun(SunCorrection mode) {
    return {.mode = mode, .r = 0.015, .d = 0.02};
}

SnowDistribution one_bin() {
    const std::array f{1.0};
    return SnowDistribution{f};
}

SnowDistribution four_bins() {
    const std::array f{0.5, 0.8, 1.2, 1.5};
    return SnowDistribution{f};
}

// Section 7's D1 factors, the log-normal at CV 0.4 over ten bins.
SnowDistribution ten_bins() {
    const std::array f{0.4974, 0.6287, 0.7228, 0.8080, 0.8930,
                       0.9838, 1.0873, 1.2154, 1.3973, 1.7664};
    std::array<double, 10> g{};
    double sum = 0.0;
    for (double v : f)
        sum += v;
    std::ranges::transform(f, g.begin(), [&](double v) { return v * 10.0 / sum; });
    return SnowDistribution{g};
}

// One step with the same forcing on every triangle of the state.
StepStatus step_all(SnowState& s, std::size_t n, const HbvSnowParameters& p,
                    const SunCorrectionParameters& c, double t, double precip,
                    double anomaly, double dt = kDay, unsigned threads = 1) {
    const std::vector<double> temperature(n, t), precipitation(n, precip);
    const std::vector<double> di =
        c.mode == SunCorrection::none ? std::vector<double>{} : std::vector<double>(n, anomaly);
    return step(p, c, s, StepForcing{temperature, precipitation, di, dt}, threads);
}

// One bin holding sp, with sw of water from a priming rain step (sw <= lw*sp).
SnowState one_bin_state(double sp, double sw) {
    auto s = SnowState::uniform(1, one_bin(), sp);
    if (sw > 0.0)
        REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 0.0, sw, 0.0) == StepStatus::ok);
    REQUIRE_THAT(s.sp()[0], WithinAbs(sp, 1e-12));
    REQUIRE_THAT(s.sw()[0], WithinAbs(sw, 1e-12));
    return s;
}

struct Snapshot {
    std::vector<double> sp, sw, swe, sca, outflow;
};

std::vector<double> copy(std::span<const double> v) { return {v.begin(), v.end()}; }

Snapshot snapshot(const SnowState& s) {
    return {copy(s.sp()), copy(s.sw()), copy(s.swe()), copy(s.sca()), copy(s.outflow())};
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() &&
           std::ranges::equal(a, b, [](double x, double y) {
               return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y);
           });
}

bool same_bits(const Snapshot& a, const Snapshot& b) {
    return same_bits(a.sp, b.sp) && same_bits(a.sw, b.sw) && same_bits(a.swe, b.swe) &&
           same_bits(a.sca, b.sca) && same_bits(a.outflow, b.outflow);
}

void check_one_bin(const SnowState& s, double sp, double sw, double outflow, double sca) {
    CHECK_THAT(s.sp()[0], WithinAbs(sp, 1e-12));
    CHECK_THAT(s.sw()[0], WithinAbs(sw, 1e-12));
    CHECK_THAT(s.outflow()[0], WithinAbs(outflow, 1e-12));
    CHECK_THAT(s.swe()[0], WithinAbs(sp + sw, 1e-12));
    CHECK(s.sca()[0] == sca);
}

}  // namespace

// ------------------------------------------------------------- H1-H12

TEST_CASE("H1 melt; H2 refreeze from H1's state", "[snow][hbv]") {
    auto s = one_bin_state(100.0, 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 5.0, 0.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 85.0, 8.5, 6.5, 1.0);
    CHECK_THAT(s.swe()[0], WithinAbs(93.5, 1e-12));

    // T = -2: potential refreeze cfr*cx*dt*(ts - T) = 0.05*3*1*2 = 0.3 < sw.
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), -2.0, 0.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 85.3, 8.2, 0.0, 1.0);
    CHECK_THAT(s.swe()[0], WithinAbs(93.5, 1e-12));
}

TEST_CASE("H3 snowfall spreads over the bins by their factors", "[snow][hbv]") {
    auto s = SnowState::uniform(1, four_bins(), 100.0);
    const auto before = copy(s.sp());
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), -1.0, 10.0, 0.0) == StepStatus::ok);
    const std::array rise{5.0, 8.0, 12.0, 15.0};
    for (std::size_t k = 0; k < 4; ++k)
        CHECK_THAT(s.sp()[k] - before[k], WithinAbs(rise[k], 1e-12));
    CHECK_THAT(s.swe()[0], WithinAbs(110.0, 1e-12));
    CHECK(s.outflow()[0] == 0.0);
}

TEST_CASE("H4 rain on snow", "[snow][hbv]") {
    // Rain 10 into sw, melt 3*2 = 6: sp 94, sw 16 capped at 9.4, outflow 6.6.
    auto s = one_bin_state(100.0, 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 2.0, 10.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 94.0, 9.4, 6.6, 1.0);
    CHECK_THAT(s.swe()[0], WithinAbs(103.4, 1e-12));
}

TEST_CASE("H5 rain on bare ground leaves at once", "[snow][hbv]") {
    auto s = SnowState::uniform(1, one_bin(), 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 2.0, 10.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 0.0, 0.0, 10.0, 0.0);
}

TEST_CASE("H5b rain on a bare bin below ts is not kept to refreeze", "[snow][hbv]") {
    // tx = -1 so T = -0.5 rains; T < ts, so water kept in the empty bin would
    // refreeze min(0.05*3*0.5, 10) = 0.075 and give sca 1.
    auto p = kHbv;
    p.tx = -1.0;
    auto s = SnowState::uniform(1, one_bin(), 0.0);
    REQUIRE(step_all(s, 1, p, sun(SunCorrection::none), -0.5, 10.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 0.0, 0.0, 10.0, 0.0);
}

TEST_CASE("H6 melt-out within a step", "[snow][hbv]") {
    // Potential melt 3*10 = 30 > sp: the melt is capped at sp.
    auto s = one_bin_state(10.0, 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 10.0, 0.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 0.0, 0.0, 10.0, 0.0);
}

TEST_CASE("H7 four bins melt each by its own snow", "[snow][hbv]") {
    auto s = SnowState::uniform(1, four_bins(), 100.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 20.0, 0.0, 0.0) == StepStatus::ok);
    const std::array sp{0.0, 20.0, 60.0, 90.0};
    const std::array sw{0.0, 2.0, 6.0, 9.0};
    for (std::size_t k = 0; k < 4; ++k) {
        CHECK_THAT(s.sp()[k], WithinAbs(sp[k], 1e-12));
        CHECK_THAT(s.sw()[k], WithinAbs(sw[k], 1e-12));
    }
    CHECK(s.sca()[0] == 0.75);
    CHECK_THAT(s.swe()[0], WithinAbs(46.75, 1e-12));
    CHECK_THAT(s.outflow()[0], WithinAbs(53.25, 1e-12));
}

TEST_CASE("H8 the corrections at their break-even; none ignores the anomaly", "[snow][hbv]") {
    // T = 4, dI = +100: melt_factor m = 3 + 1.5 = 4.5, melt 18; threshold
    // theta = -2, melt 3*6 = 18; none melts 3*4 = 12 whatever r, d and dI say.
    const auto [mode, melt] = GENERATE(table<SunCorrection, double>({
        {SunCorrection::melt_factor, 18.0},
        {SunCorrection::threshold, 18.0},
        {SunCorrection::none, 12.0},
    }));
    CAPTURE(static_cast<int>(mode));
    auto s = one_bin_state(100.0, 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(mode), 4.0, 0.0, 100.0) == StepStatus::ok);
    CHECK_THAT(s.sp()[0], WithinAbs(100.0 - melt, 1e-12));
    CHECK_THAT(s.swe()[0] + s.outflow()[0], WithinAbs(100.0, 1e-12));
}

TEST_CASE("H8b a shaded slope with water in the pack: no melt, no refreeze", "[snow][hbv]") {
    // dI = -300, T = 4. melt_factor: m = max(0, 3 - 4.5) = 0. threshold:
    // theta = +6 > T, and T > ts, so no refreeze. Both: nothing moves.
    const auto mode = GENERATE(SunCorrection::melt_factor, SunCorrection::threshold);
    CAPTURE(static_cast<int>(mode));
    auto s = one_bin_state(100.0, 5.0);
    REQUIRE(step_all(s, 1, kHbv, sun(mode), 4.0, 0.0, -300.0) == StepStatus::ok);
    check_one_bin(s, 100.0, 5.0, 0.0, 1.0);
    CHECK_THAT(s.swe()[0], WithinAbs(105.0, 1e-12));
}

TEST_CASE("H9 a sub-day step scales the melt", "[snow][hbv]") {
    auto s = one_bin_state(100.0, 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 5.0, 0.0, 0.0, 3600.0) ==
            StepStatus::ok);
    CHECK_THAT(s.sp()[0], WithinAbs(100.0 - 0.625, 1e-12));
    CHECK_THAT(s.sw()[0], WithinAbs(0.625, 1e-12));
}

TEST_CASE("H10 precipitation at T = tx is rain", "[snow][hbv]") {
    auto s = one_bin_state(100.0, 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), 0.0, 10.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 100.0, 10.0, 0.0, 1.0);
    CHECK_THAT(s.swe()[0], WithinAbs(110.0, 1e-12));
}

TEST_CASE("H11 refreeze capped by the water in the pack", "[snow][hbv]") {
    // ts = 0 (kHbv), T = -10: potential refreeze 0.05*3*1*10 = 1.5 > sw = 1.
    REQUIRE(kHbv.ts == 0.0);
    auto s = one_bin_state(100.0, 1.0);
    REQUIRE(step_all(s, 1, kHbv, sun(SunCorrection::none), -10.0, 0.0, 0.0) == StepStatus::ok);
    check_one_bin(s, 101.0, 0.0, 0.0, 1.0);
    CHECK(s.sw()[0] >= 0.0);
}

TEST_CASE("H12 melt below ts in the threshold form, and no refreeze after it", "[snow][hbv]") {
    // T = -1, dI = +100: threshold theta = -2 < T, melt 3: sp 97, sw 3 (cap 9.7).
    // melt_factor keeps theta = ts > T: no melt; refreeze min(0.15, 0) = 0.
    const auto [mode, sp, sw] = GENERATE(table<SunCorrection, double, double>({
        {SunCorrection::threshold, 97.0, 3.0},
        {SunCorrection::melt_factor, 100.0, 0.0},
    }));
    CAPTURE(static_cast<int>(mode));
    auto s = one_bin_state(100.0, 0.0);
    REQUIRE(step_all(s, 1, kHbv, sun(mode), -1.0, 0.0, 100.0) == StepStatus::ok);
    check_one_bin(s, sp, sw, 0.0, 1.0);
}

TEST_CASE("a state of zero triangles steps on empty forcing", "[snow][hbv]") {
    auto s = SnowState::uniform(0, ten_bins(), 600.0);
    CHECK(s.swe().empty());
    CHECK(s.sp().empty());
    CHECK(step_all(s, 0, kHbv, sun(SunCorrection::threshold), 3.0, 1.0, 50.0) == StepStatus::ok);
}

TEST_CASE("uniform spreads the initial SWE by the factors, all ice", "[snow][hbv]") {
    const auto s = SnowState::uniform(3, four_bins(), 100.0);
    const std::array sp{50.0, 80.0, 120.0, 150.0};
    REQUIRE(s.sp().size() == 12);
    REQUIRE(s.sw().size() == 12);
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t k = 0; k < 4; ++k) {
            CHECK_THAT(s.sp()[i * 4 + k], WithinAbs(sp[k], 1e-12));
            CHECK(s.sw()[i * 4 + k] == 0.0);
        }
        CHECK_THAT(s.swe()[i], WithinAbs(100.0, 1e-12));
        CHECK(s.sca()[i] == 1.0);
        CHECK(s.outflow()[i] == 0.0);
    }
}

// ------------------------------------------------------------- P1-P5

namespace {

// Random forcing at section 7's ranges: T in [-30, 30], P in [0, 100] on
// three steps in ten (else 0), dI in [-600, 600]; daily steps. Swe drifts but
// stays at the 1e4 mm scale the mass tolerance is stated for (section 5),
// which run() checks.
struct Forcing {
    std::vector<double> t, p, di;
};

Forcing draw(std::mt19937_64& rng, std::size_t n, bool with_anomaly) {
    std::uniform_real_distribution<double> temperature(-30.0, 30.0), rain(0.0, 100.0),
        anomaly(-600.0, 600.0), coin(0.0, 1.0);
    Forcing f{std::vector<double>(n), std::vector<double>(n), {}};
    for (std::size_t i = 0; i < n; ++i) {
        f.t[i] = temperature(rng);
        f.p[i] = coin(rng) < 0.3 ? rain(rng) : 0.0;
    }
    if (with_anomaly) {
        f.di.resize(n);
        for (auto& v : f.di)
            v = anomaly(rng);
    }
    return f;
}

// I1 and I2 after every step of a seeded run; returns the largest swe seen.
double check_mass_and_bounds(const SnowDistribution& dist, SunCorrection mode,
                             std::uint64_t seed, std::size_t n, std::size_t steps,
                             double initial) {
    std::mt19937_64 rng{seed};
    auto s = SnowState::uniform(n, dist, initial);
    const std::size_t bins = dist.factors().size();
    const auto lw = kHbv.lw;
    double largest = initial;
    for (std::size_t j = 0; j < steps; ++j) {
        const auto f = draw(rng, n, mode != SunCorrection::none);
        const auto before = copy(s.swe());
        REQUIRE(step(kHbv, sun(mode), s, StepForcing{f.t, f.p, f.di, kDay}, 1) == StepStatus::ok);
        for (std::size_t i = 0; i < n; ++i) {
            const double expected = before[i] + f.p[i] - s.outflow()[i];
            const double tol = 1e-9 * std::max(1.0, before[i] + f.p[i]);
            if (std::abs(s.swe()[i] - expected) > tol) {
                CAPTURE(seed, j, i, before[i], f.p[i], s.swe()[i], s.outflow()[i]);
                FAIL("I1: mass not conserved");
            }
            if (!(s.sca()[i] >= 0.0 && s.sca()[i] <= 1.0 && s.outflow()[i] >= 0.0)) {
                CAPTURE(seed, j, i, s.sca()[i], s.outflow()[i]);
                FAIL("I2: sca or outflow out of bounds");
            }
            for (std::size_t k = 0; k < bins; ++k) {
                const double sp = s.sp()[i * bins + k], sw = s.sw()[i * bins + k];
                if (!(sp >= 0.0 && sw >= 0.0 && sw <= lw * sp + 1e-12)) {
                    CAPTURE(seed, j, i, k, sp, sw);
                    FAIL("I2: a bin out of bounds");
                }
            }
            largest = std::max(largest, s.swe()[i]);
        }
    }
    return largest;
}

}  // namespace

TEST_CASE("P1, P2 mass and bounds over 1e4 random steps, every mode", "[snow][hbv][property]") {
    const auto mode =
        GENERATE(SunCorrection::none, SunCorrection::melt_factor, SunCorrection::threshold);
    const auto seed = GENERATE(as<std::uint64_t>{}, 46, 4646);
    CAPTURE(static_cast<int>(mode), seed);
    // Checked at: swe start 1e4 mm, the largest swe the run reaches (asserted
    // below to stay within 2e4, the scale section 5 states).
    const double largest = check_mass_and_bounds(ten_bins(), mode, seed, 8, 10000, 1e4);
    CHECK(largest <= 2e4);
}

TEST_CASE("P3 flat is plain: dI = 0 gives bit-identical states in all modes",
          "[snow][hbv][property]") {
    const std::size_t n = 16;
    std::mt19937_64 rng{46};
    auto none = SnowState::uniform(n, ten_bins(), 600.0);
    auto mf = SnowState::uniform(n, ten_bins(), 600.0);
    auto th = SnowState::uniform(n, ten_bins(), 600.0);
    const std::vector<double> zero(n, 0.0);
    for (int j = 0; j < 2000; ++j) {
        const auto f = draw(rng, n, false);
        REQUIRE(step(kHbv, sun(SunCorrection::none), none, StepForcing{f.t, f.p, {}, kDay}, 1) ==
                StepStatus::ok);
        REQUIRE(step(kHbv, sun(SunCorrection::melt_factor), mf, StepForcing{f.t, f.p, zero, kDay},
                     1) == StepStatus::ok);
        REQUIRE(step(kHbv, sun(SunCorrection::threshold), th, StepForcing{f.t, f.p, zero, kDay},
                     1) == StepStatus::ok);
        if (!same_bits(snapshot(none), snapshot(mf)) || !same_bits(snapshot(none), snapshot(th))) {
            CAPTURE(j);
            FAIL("I3: a flat triangle's state differs between modes");
        }
    }
}

namespace {

// Plain single-point HBV-96, section 4.2's steps 1-7 with one bin, mode none.
struct Point {
    double sp, sw, out;
};

Point hbv96(Point s, const HbvSnowParameters& p, double t, double precip, double dt) {
    const double days = dt / kDay;
    const double snow = t < p.tx ? precip : 0.0;
    const double rain = precip - snow;
    s.out = 0.0;
    s.sp += snow;
    if (s.sp > 0.0)
        s.sw += rain;
    else
        s.out += rain;
    if (t > p.ts) {
        const double melt = std::min(p.cx * days * (t - p.ts), s.sp);
        s.sp -= melt;
        s.sw += melt;
    } else if (t < p.ts) {
        const double refreeze = std::min(p.cfr * p.cx * days * (p.ts - t), s.sw);
        s.sp += refreeze;
        s.sw -= refreeze;
    }
    s.out += std::max(0.0, s.sw - p.lw * s.sp);
    s.sw = std::min(s.sw, p.lw * s.sp);
    return s;
}

}  // namespace

TEST_CASE("P4 one bin is plain HBV-96 for one point", "[snow][hbv][property]") {
    const auto seed = GENERATE(as<std::uint64_t>{}, 46, 4646);
    std::mt19937_64 rng{seed};
    const std::size_t n = 4;
    auto s = SnowState::uniform(n, one_bin(), 600.0);
    std::vector<Point> ref(n, Point{600.0, 0.0, 0.0});
    std::uniform_real_distribution<double> step_hours(1.0, 24.0);
    for (int j = 0; j < 10000; ++j) {
        const auto f = draw(rng, n, false);
        const double dt = 3600.0 * std::floor(step_hours(rng));
        REQUIRE(step(kHbv, sun(SunCorrection::none), s, StepForcing{f.t, f.p, {}, dt}, 1) ==
                StepStatus::ok);
        for (std::size_t i = 0; i < n; ++i) {
            ref[i] = hbv96(ref[i], kHbv, f.t[i], f.p[i], dt);
            const double tol = 1e-12 * std::max(1.0, ref[i].sp + ref[i].sw + f.p[i]);
            const bool sp_ok = std::abs(s.sp()[i] - ref[i].sp) <= tol;
            const bool sw_ok = std::abs(s.sw()[i] - ref[i].sw) <= tol;
            const bool out_ok = std::abs(s.outflow()[i] - ref[i].out) <= tol;
            const bool sca_ok = s.sca()[i] == (ref[i].sp > 0.0 ? 1.0 : 0.0);
            if (!(sp_ok && sw_ok && out_ok && sca_ok)) {
                CAPTURE(seed, j, i, s.sp()[i], ref[i].sp, s.sw()[i], ref[i].sw, s.outflow()[i],
                        ref[i].out, s.sca()[i]);
                FAIL("I4: one bin departs from HBV-96");
            }
            ref[i].sp = s.sp()[i];  // follow the kernel, so one rounding never compounds
            ref[i].sw = s.sw()[i];
        }
    }
}

TEST_CASE("P5 the step is independent of threads: 1 and 8 bit-identical",
          "[snow][hbv][property]") {
    const auto mode =
        GENERATE(SunCorrection::none, SunCorrection::melt_factor, SunCorrection::threshold);
    CAPTURE(static_cast<int>(mode));
    const std::size_t n = 20000;  // enough blocks that eight workers start
    std::mt19937_64 rng{46};
    auto one = SnowState::uniform(n, ten_bins(), 600.0);
    auto eight = SnowState::uniform(n, ten_bins(), 600.0);
    for (int j = 0; j < 20; ++j) {
        const auto f = draw(rng, n, mode != SunCorrection::none);
        const StepForcing forcing{f.t, f.p, f.di, kDay};
        REQUIRE(step(kHbv, sun(mode), one, forcing, 1) == StepStatus::ok);
        REQUIRE(step(kHbv, sun(mode), eight, forcing, 8) == StepStatus::ok);
        REQUIRE(same_bits(snapshot(one), snapshot(eight)));
    }
}

// ------------------------------------------------------------- V1, V2

TEST_CASE("V1 SnowDistribution refuses bad factors, naming the rule", "[snow][hbv][refusal]") {
    using V = std::vector<double>;
    const auto [factors, reason] = GENERATE(table<V, std::string>({
        {V{}, "between 1 and 64"},
        {V(65, 1.0), "between 1 and 64"},
        {V{0.5, kNaN, 1.2, 1.5}, "finite"},
        {V{0.5, 0.8, 1.2, kInf}, "finite"},
        {V{-0.1, 0.5, 1.2, 2.4}, "negative"},
        {V{1.5, 1.2, 0.8, 0.5}, "ascending"},
        {V{0.5, 0.8, 1.2, 1.5 + 4e-6}, "mean"},
    }));
    CAPTURE(factors, reason);
    REQUIRE_THROWS_MATCHES(SnowDistribution{factors}, std::invalid_argument,
                           MessageMatches(ContainsSubstring(reason)));
}

TEST_CASE("V1 SnowState::uniform refuses a negative or non-finite initial SWE",
          "[snow][hbv][refusal]") {
    const double swe = GENERATE(-1.0, kNaN, kInf);
    CAPTURE(swe);
    REQUIRE_THROWS_MATCHES(SnowState::uniform(3, four_bins(), swe), std::invalid_argument,
                           MessageMatches(ContainsSubstring("initial")));
}

namespace {

// Three triangles, four bins, one ok step behind them so every slot (sp, sw,
// outflow) is non-trivial; the refused step's bad value is on the last
// triangle, so a kernel that writes before it validates is seen.
SnowState used_state() {
    auto s = SnowState::uniform(3, four_bins(), 100.0);
    REQUIRE(step_all(s, 3, kHbv, sun(SunCorrection::none), 1.0, 5.0, 0.0) == StepStatus::ok);
    REQUIRE(s.outflow()[2] > 0.0);
    REQUIRE(s.sw()[11] > 0.0);
    return s;
}

struct Case {
    std::string name;
    StepStatus want = StepStatus::ok;
    HbvSnowParameters p = kHbv;
    SunCorrectionParameters c = sun(SunCorrection::melt_factor);
    std::vector<double> t = {1.0, 1.0, 1.0}, precip = {2.0, 2.0, 2.0}, di = {50.0, 50.0, 50.0};
    double dt = kDay;
};

template <class Edit>
Case make(std::string name, StepStatus want, Edit edit) {
    Case c;
    c.name = std::move(name);
    c.want = want;
    edit(c);
    return c;
}

std::vector<Case> refusals() {
    using S = StepStatus;
    using P = HbvSnowParameters;
    std::vector<Case> cs;
    const auto param = [&](std::string name, double P::* field, double v) {
        cs.push_back(make(std::move(name), S::invalid_parameters, [&](Case& c) { c.p.*field = v; }));
    };
    param("cx -0.1", &P::cx, -0.1);
    param("lw -0.1", &P::lw, -0.1);
    param("cfr -0.1", &P::cfr, -0.1);
    for (auto [name, field] : {std::pair{"tx", &P::tx}, std::pair{"ts", &P::ts},
                               std::pair{"cx", &P::cx}, std::pair{"lw", &P::lw},
                               std::pair{"cfr", &P::cfr}})
        param(std::string{name} + " NaN", field, kNaN);
    // r and d are checked whatever the mode, so a NaN refuses under none too.
    for (auto mode : {SunCorrection::none, SunCorrection::melt_factor}) {
        const bool off = mode == SunCorrection::none;
        cs.push_back(make("r NaN", S::invalid_parameters, [&](Case& c) {
            c.c = {.mode = mode, .r = kNaN, .d = 0.02};
            if (off)
                c.di.clear();
        }));
        cs.push_back(make("d NaN", S::invalid_parameters, [&](Case& c) {
            c.c = {.mode = mode, .r = 0.015, .d = kNaN};
            if (off)
                c.di.clear();
        }));
    }
    cs.push_back(make("mode outside the enum", S::invalid_parameters,
                      [](Case& c) { c.c.mode = static_cast<SunCorrection>(7); }));

    cs.push_back(make("temperature n - 1", S::size_mismatch, [](Case& c) { c.t.pop_back(); }));
    cs.push_back(make("temperature n + 1", S::size_mismatch, [](Case& c) { c.t.push_back(1.0); }));
    cs.push_back(make("precipitation n - 1", S::size_mismatch, [](Case& c) { c.precip.pop_back(); }));
    cs.push_back(
        make("precipitation n + 1", S::size_mismatch, [](Case& c) { c.precip.push_back(2.0); }));
    cs.push_back(make("dI of length n under none", S::size_mismatch,
                      [](Case& c) { c.c = sun(SunCorrection::none); }));
    cs.push_back(make("dI empty under melt_factor", S::size_mismatch, [](Case& c) { c.di.clear(); }));
    cs.push_back(make("dI empty under threshold", S::size_mismatch, [](Case& c) {
        c.c = sun(SunCorrection::threshold);
        c.di.clear();
    }));
    cs.push_back(make("dI of length n - 1", S::size_mismatch, [](Case& c) { c.di.pop_back(); }));

    for (double bad : {kNaN, kInf, -kInf}) {
        const std::string v = std::isnan(bad) ? "NaN" : (bad > 0 ? "+inf" : "-inf");
        cs.push_back(make("temperature " + v, S::non_finite_forcing, [&](Case& c) { c.t[2] = bad; }));
        cs.push_back(
            make("precipitation " + v, S::non_finite_forcing, [&](Case& c) { c.precip[2] = bad; }));
        cs.push_back(make("dI " + v, S::non_finite_forcing, [&](Case& c) { c.di[2] = bad; }));
    }
    cs.push_back(make("precipitation -1e-12", S::negative_precipitation,
                      [](Case& c) { c.precip[2] = -1e-12; }));
    for (double dt : {0.0, -3600.0, kNaN, kInf})
        cs.push_back(make("dt " + std::to_string(dt), S::bad_step_length, [&](Case& c) { c.dt = dt; }));
    return cs;
}

}  // namespace

TEST_CASE("V1 a refused step returns its status and leaves the state bit for bit",
          "[snow][hbv][refusal]") {
    for (const auto& c : refusals()) {
        CAPTURE(c.name);
        auto s = used_state();
        const auto before = snapshot(s);
        CHECK(step(c.p, c.c, s, StepForcing{c.t, c.precip, c.di, c.dt}, 1) == c.want);
        CHECK(same_bits(before, snapshot(s)));
    }
}

TEST_CASE("V1 the control: the refusal cases' good forcing is accepted", "[snow][hbv][refusal]") {
    auto s = used_state();
    const Case good;
    CHECK(step(good.p, good.c, s, StepForcing{good.t, good.precip, good.di, good.dt}, 1) ==
          StepStatus::ok);
}

TEST_CASE("V2 factors within 1e-9 of mean 1 are accepted and renormalised",
          "[snow][hbv][refusal]") {
    // Mean off by 1e-10: accepted, and returned with mean 1 to 1e-15.
    const std::array f{0.5, 0.8, 1.2, 1.5 + 4e-10};
    const SnowDistribution dist{f};
    const auto g = dist.factors();
    REQUIRE(g.size() == 4);
    double sum = 0.0;
    for (double v : g)
        sum += v;
    CHECK_THAT(sum / 4.0, WithinAbs(1.0, 1e-15));
    CHECK(std::ranges::is_sorted(g));

    // P1's mass check at 1e-9 * max(1 mm, swe_before + P), P up to 100 mm.
    const auto mode = GENERATE(SunCorrection::none, SunCorrection::threshold);
    CAPTURE(static_cast<int>(mode));
    check_mass_and_bounds(dist, mode, 46, 8, 2000, 600.0);
}
