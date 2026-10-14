// Increment 46 PR 8 (docs/increments/46-hbv-snow.md, section 12): linear bins,
// a continuous snow-covered area in the active bin, in
// include/terrain/snow/hbv_snow.hpp.
//
//   L0-L10  section 12.7's fixtures, one bin, factors (1.0), deep (1.5), so the
//           shallow edge is 0.5; each expected value is printed by
//           docs/increments/46-probes/probe_sca.py into probe_sca.out;
//   P6-P8   invariants I1 and I2 (extended), I7 and I8 over seeded random
//           forcing on random linear distributions (fixed seeds: deterministic);
//   V3      the distribution's refusals of a bad deep edge.
//
// I6 (a flat distribution is the kernel before this PR, bit for bit) is held
// by test_hbv_snow.cpp's H1-H12 and P1-P5, unedited (section 12.7, P9).
//
// Interface under test (section 12.3; the accessors are named here):
// - SnowDistribution(factors, deep): deep empty means flat; deep() is each
//   bin's deep-edge factor divided by the factors' mean, equal to factors()
//   for a flat distribution.
// - SnowState::deep(): n*N values, triangle-major as sp(), the solid depth at
//   each bin's deep edge; uniform() sets it to deep_k * initial swe.
//
// Invariant-critical (section 12.7's mutation targets): mutation-tested.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/snow/hbv_snow.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>
#include <string>
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
constexpr double kHour = 3600.0;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

// Section 12.7's fixture parameters: tx = ts = 0, cx = 3, lw = 0.1, cfr = 0.05.
constexpr HbvSnowParameters kHbv{.tx = 0.0, .ts = 0.0, .cx = 3.0, .lw = 0.1, .cfr = 0.05};
constexpr SunCorrectionParameters kNone{.mode = SunCorrection::none, .r = 0.015, .d = 0.02};

constexpr SunCorrectionParameters sun(SunCorrection mode) {
    return {.mode = mode, .r = 0.015, .d = 0.02};
}

HbvSnowParameters with(double lw, double cfr) {
    auto p = kHbv;
    p.lw = lw;
    p.cfr = cfr;
    return p;
}

// One bin, factor 1, deep edge 1.5 (shallow edge 0.5, cv 0.2887).
SnowDistribution linear_bin() {
    const std::array f{1.0}, d{1.5};
    return SnowDistribution{f, d};
}

// One step of mode none on a one-triangle state.
void step1(SnowState& s, const HbvSnowParameters& p, double t, double precip, double dt = kDay) {
    const std::array temperature{t}, precipitation{precip};
    REQUIRE(step(p, kNone, s, StepForcing{temperature, precipitation, {}, dt}, 1) ==
            StepStatus::ok);
}

struct Expect {
    double sp, deep, sw, outflow, sca;
};

// Absolute 1e-12 at the fixtures' scale (values up to 150 mm).
void check(const SnowState& s, const Expect& e) {
    CHECK_THAT(s.sp()[0], WithinAbs(e.sp, 1e-12));
    CHECK_THAT(s.deep()[0], WithinAbs(e.deep, 1e-12));
    CHECK_THAT(s.sw()[0], WithinAbs(e.sw, 1e-12));
    CHECK_THAT(s.outflow()[0], WithinAbs(e.outflow, 1e-12));
    CHECK_THAT(s.swe()[0], WithinAbs(e.sp + e.sw, 1e-12));
    CHECK_THAT(s.sca()[0], WithinAbs(e.sca, 1e-12));
}

// L0 and then L2: 60 mm of melt potential takes it from covered to partly bare.
SnowState l2_state(const HbvSnowParameters& p) {
    auto s = SnowState::uniform(1, linear_bin(), 100.0);
    step1(s, p, 20.0, 0.0);
    return s;
}

}  // namespace

// ------------------------------------------------------------- L0-L10

TEST_CASE("L0 a fresh linear bin is covered, its cover clamped to 1", "[snow][hbv][linear]") {
    const auto s = SnowState::uniform(1, linear_bin(), 100.0);
    REQUIRE(s.deep().size() == 1);
    check(s, {.sp = 100.0, .deep = 150.0, .sw = 0.0, .outflow = 0.0, .sca = 1.0});
    CHECK(s.sca()[0] == 1.0);  // not 2*100/150 = 1.33
}

TEST_CASE("L1 melt within the shallow edge: today's melt, still covered", "[snow][hbv][linear]") {
    auto s = SnowState::uniform(1, linear_bin(), 100.0);
    step1(s, kHbv, 5.0, 0.0);  // M = 15 <= lo = 50
    check(s, {.sp = 85.0, .deep = 135.0, .sw = 8.5, .outflow = 6.5, .sca = 1.0});
    CHECK(s.sca()[0] == 1.0);
}

TEST_CASE("L2-L4 the front moves: covered, partly bare, bare", "[snow][hbv][linear]") {
    auto s = l2_state(kHbv);  // M = 60 > lo = 50: deep 90, sp = 90^2 / (4 * 50)
    check(s, {.sp = 40.5, .deep = 90.0, .sw = 4.05, .outflow = 55.45, .sca = 0.9});
    CHECK_THAT(s.swe()[0], WithinAbs(44.55, 1e-12));

    step1(s, kHbv, 10.0, 0.0);  // L3, M = 30 from partly bare: sp = 40.5 * (60/90)^2
    check(s, {.sp = 18.0, .deep = 60.0, .sw = 1.8, .outflow = 24.75, .sca = 0.6});

    step1(s, kHbv, 20.0, 0.0);  // L4, M = 60 >= deep: melts out
    check(s, {.sp = 0.0, .deep = 0.0, .sw = 0.0, .outflow = 19.8, .sca = 0.0});
    CHECK(s.sp()[0] == 0.0);
    CHECK(s.deep()[0] == 0.0);
    CHECK(s.sca()[0] == 0.0);
}

TEST_CASE("L5 rain on a partly bare bin enters by its cover", "[snow][hbv][linear]") {
    const auto p = with(2.0, 0.05);
    auto s = l2_state(p);  // lw = 2 keeps the melt: sp 40.5, deep 90, sw 59.5, cover 0.9
    REQUIRE_THAT(s.sw()[0], WithinAbs(59.5, 1e-12));
    step1(s, p, 0.0, 10.0);  // T = tx = ts: rain, no melt, no refreeze
    check(s, {.sp = 40.5, .deep = 90.0, .sw = 68.5, .outflow = 1.0, .sca = 0.9});
    CHECK_THAT(s.swe()[0], WithinAbs(109.0, 1e-12));
}

TEST_CASE("L6 refreeze on a partly bare bin, by its cover, keeps the cover",
          "[snow][hbv][linear]") {
    const auto p = with(2.0, 0.05);
    auto s = l2_state(p);
    step1(s, p, -10.0, 0.0);  // F = 0.9 * 0.05 * 3 * 10 = 1.35; deep 90 * 41.85 / 40.5
    check(s, {.sp = 41.85, .deep = 93.0, .sw = 58.15, .outflow = 0.0, .sca = 0.9});
    CHECK_THAT(s.swe()[0], WithinAbs(100.0, 1e-12));
}

TEST_CASE("L7 snow on a partly bare bin spreads the old snow as a wedge; L7c melts it",
          "[snow][hbv][linear]") {
    const auto p = with(0.1, 0.0);
    auto s = l2_state(p);
    step1(s, p, 10.0, 0.0);  // L3 with cfr = 0: sp 18, deep 60, sw 1.8
    step1(s, p, -1.0, 10.0);  // wedge deep = 36, then sp += 10, deep += 15
    check(s, {.sp = 28.0, .deep = 51.0, .sw = 1.8, .outflow = 0.0, .sca = 1.0});
    CHECK_THAT(s.swe()[0], WithinAbs(29.8, 1e-12));
    CHECK(s.sca()[0] == 1.0);

    // L7c: M = 10 > lo = 5, covered to the front: deep 41, sp = 41^2 / (4 * 23).
    step1(s, p, 10.0 / 3.0, 0.0);
    const double sp = 41.0 * 41.0 / 92.0;
    check(s, {.sp = sp, .deep = 41.0, .sw = 0.1 * sp, .outflow = 29.8 - 1.1 * sp,
              .sca = 41.0 / 46.0});
    CHECK_THAT(sp, WithinAbs(18.271739, 1e-6));  // probe_sca.out's L7c
}

TEST_CASE("L8 snow on a bare bin by both edges' factors; L8b melts it", "[snow][hbv][linear]") {
    const auto p = with(0.1, 0.0);
    auto s = SnowState::uniform(1, linear_bin(), 0.0);
    check(s, {.sp = 0.0, .deep = 0.0, .sw = 0.0, .outflow = 0.0, .sca = 0.0});
    step1(s, p, -1.0, 10.0);
    check(s, {.sp = 10.0, .deep = 15.0, .sw = 0.0, .outflow = 0.0, .sca = 1.0});

    step1(s, p, 2.5, 0.0);  // M = 7.5 > lo = 5: deep 7.5, sp = 7.5^2 / (4 * 5)
    check(s, {.sp = 2.8125, .deep = 7.5, .sw = 0.28125, .outflow = 6.90625, .sca = 0.75});
}

TEST_CASE("L9 at T = ts = theta a partly bare bin is left alone", "[snow][hbv][linear]") {
    auto s = l2_state(kHbv);
    step1(s, kHbv, 0.0, 0.0);
    check(s, {.sp = 40.5, .deep = 90.0, .sw = 4.05, .outflow = 0.0, .sca = 0.9});
}

namespace {

// Clamp((deep0 - melt) / (deep0 - lo0), 0, 1): I8's depletion curve.
double depletion(double deep0, double lo0, double melt) {
    return std::clamp((deep0 - melt) / (deep0 - lo0), 0.0, 1.0);
}

double cover(double sp, double deep) {
    if (sp == 0.0)
        return 0.0;
    return 2.0 * sp >= deep ? 1.0 : 2.0 * sp / deep;
}

}  // namespace

TEST_CASE("L10 400 hourly melt steps follow the depletion curve", "[snow][hbv][linear]") {
    // probe_sca.py's L10: T = 0.5 + (k mod 7), cfr = 0; 174.625 mm of melt
    // potential in all, through covered (to 50 mm), partly bare and bare (150 mm).
    // Scale: 150 mm, at most 6.7e-16 in the probe.
    const auto p = with(0.1, 0.0);
    auto s = SnowState::uniform(1, linear_bin(), 100.0);
    double melt = 0.0;
    bool seen_covered = false, seen_partly = false;
    for (int k = 0; k < 400; ++k) {
        const double t = 0.5 + static_cast<double>(k % 7);
        step1(s, p, t, 0.0, kHour);
        melt += 3.0 / 24.0 * t;
        const double want = depletion(150.0, 50.0, melt);
        seen_covered = seen_covered || want == 1.0;
        seen_partly = seen_partly || (want > 0.0 && want < 1.0);
        if (std::abs(s.sca()[0] - want) > 1e-12) {
            CAPTURE(k, melt, s.sca()[0], want);
            FAIL("I8: the cover departs from the depletion curve");
        }
    }
    CHECK_THAT(melt, WithinAbs(174.625, 1e-9));
    CHECK((seen_covered && seen_partly));
    CHECK(s.sca()[0] == 0.0);
    CHECK(s.sp()[0] == 0.0);
}

// ------------------------------------------------------------- P6-P8

namespace {

// A random linear distribution: N in 1..8, ascending factors in [0.2, 2]
// before the division by their mean, deep_k = f_k * (1 + u), u in [0, 1]
// with the ends drawn one time in eight each.
SnowDistribution random_linear(std::mt19937_64& rng) {
    std::uniform_int_distribution<std::size_t> bins(1, 8);
    std::uniform_real_distribution<double> factor(0.2, 2.0), u(0.0, 1.0), coin(0.0, 1.0);
    std::vector<double> f(bins(rng));
    for (auto& v : f)
        v = factor(rng);
    std::ranges::sort(f);
    double sum = 0.0;
    for (double v : f)
        sum += v;
    for (auto& v : f)
        v *= static_cast<double>(f.size()) / sum;
    std::vector<double> d(f.size());
    for (std::size_t k = 0; k < f.size(); ++k) {
        const double c = coin(rng);
        const double spread = c < 0.125 ? 0.0 : c < 0.25 ? 1.0 : u(rng);
        d[k] = f[k] * (1.0 + spread);
    }
    return SnowDistribution{f, d};
}

struct Forcing {
    std::vector<double> t, p, di;
};

// Spring-like forcing, so bins are often partly bare: T in [-10, 15], P in
// [0, 20] on one step in five (snow below tx = 0, else rain), dI in [-300, 300].
Forcing draw(std::mt19937_64& rng, std::size_t n, bool with_anomaly, bool snow) {
    std::uniform_real_distribution<double> temperature(-10.0, 15.0), precip(0.0, 20.0),
        anomaly(-300.0, 300.0), coin(0.0, 1.0);
    Forcing f{std::vector<double>(n), std::vector<double>(n, 0.0), {}};
    for (std::size_t i = 0; i < n; ++i) {
        f.t[i] = temperature(rng);
        if (coin(rng) < 0.2 && (snow || f.t[i] >= kHbv.tx))
            f.p[i] = precip(rng);
    }
    if (with_anomaly) {
        f.di.resize(n);
        for (auto& v : f.di)
            v = anomaly(rng);
    }
    return f;
}

std::vector<double> copy(std::span<const double> v) { return {v.begin(), v.end()}; }

}  // namespace

TEST_CASE("P6 mass and bounds (I1, I2 extended) on random linear distributions",
          "[snow][hbv][linear][property]") {
    const auto mode =
        GENERATE(SunCorrection::none, SunCorrection::melt_factor, SunCorrection::threshold);
    const auto seed = GENERATE(range(std::uint64_t{0}, std::uint64_t{8}));
    CAPTURE(static_cast<int>(mode), seed);
    std::mt19937_64 rng{4600 + seed};
    const auto dist = random_linear(rng);
    const std::size_t n = 8, bins = dist.factors().size();
    std::uniform_real_distribution<double> initial(0.0, 1000.0);
    auto s = SnowState::uniform(n, dist, initial(rng));
    // What the run went through, so a pass is not a run that never reached a
    // partly bare bin with rain, snow, melt or refreeze on it.
    int partly_rain = 0, partly_snow = 0, partly_melt = 0, partly_freeze = 0;
    for (int j = 0; j < 3000; ++j) {
        const auto f = draw(rng, n, mode != SunCorrection::none, true);
        const auto before = copy(s.swe());
        const auto sp0 = copy(s.sp()), deep0 = copy(s.deep()), sw0 = copy(s.sw());
        REQUIRE(step(kHbv, sun(mode), s, StepForcing{f.t, f.p, f.di, kHour * 6}, 1) ==
                StepStatus::ok);
        for (std::size_t i = 0; i < n; ++i) {
            // Scale: swe up to about 2e3 mm (1e3 start, 20 mm snowfalls).
            const double expected = before[i] + f.p[i] - s.outflow()[i];
            if (std::abs(s.swe()[i] - expected) > 1e-9 * std::max(1.0, before[i] + f.p[i])) {
                CAPTURE(j, i, before[i], f.p[i], s.swe()[i], s.outflow()[i]);
                FAIL("I1: mass not conserved");
            }
            double covers = 0.0;
            for (std::size_t k = 0; k < bins; ++k) {
                const std::size_t at = i * bins + k;
                const double sp = s.sp()[at], deep = s.deep()[at], sw = s.sw()[at];
                if (!(sp >= 0.0 && sp <= deep && (sp == 0.0) == (deep == 0.0) && sw >= 0.0 &&
                      sw <= kHbv.lw * sp + 1e-12)) {
                    CAPTURE(j, i, k, sp, deep, sw);
                    FAIL("I2: a bin out of bounds");
                }
                covers += cover(sp, deep);
                const bool partly = sp0[at] > 0.0 && 2.0 * sp0[at] < deep0[at];
                if (partly && f.p[i] > 0.0)
                    ++(f.t[i] < kHbv.tx ? partly_snow : partly_rain);
                else if (partly && sp < sp0[at] && f.p[i] == 0.0)
                    ++partly_melt;
                else if (partly && sw0[at] > 0.0 && sp > sp0[at])
                    ++partly_freeze;
            }
            const double sca = s.sca()[i];
            if (!(sca >= 0.0 && sca <= 1.0 && s.outflow()[i] >= 0.0 &&
                  std::abs(sca - covers / static_cast<double>(bins)) <= 1e-12)) {
                CAPTURE(j, i, sca, covers, s.outflow()[i]);
                FAIL("I2: sca is not the bins' mean cover in [0, 1], or outflow is negative");
            }
        }
    }
    CAPTURE(partly_rain, partly_snow, partly_melt, partly_freeze);
    CHECK((partly_rain > 0 && partly_snow > 0 && partly_melt > 0 && partly_freeze > 0));
}

TEST_CASE("P7 without snowfall no bin's cover rises (I7)", "[snow][hbv][linear][property]") {
    const auto mode =
        GENERATE(SunCorrection::none, SunCorrection::melt_factor, SunCorrection::threshold);
    const auto seed = GENERATE(range(std::uint64_t{0}, std::uint64_t{8}));
    CAPTURE(static_cast<int>(mode), seed);
    std::mt19937_64 rng{4700 + seed};
    const auto dist = random_linear(rng);
    const std::size_t n = 8, bins = dist.factors().size();
    auto s = SnowState::uniform(n, dist, 600.0);
    int fell = 0;
    for (int j = 0; j < 3000; ++j) {
        // Rain only above tx, never snow: melt, refreeze and rain on the pack.
        const auto f = draw(rng, n, mode != SunCorrection::none, false);
        const auto sp0 = copy(s.sp()), deep0 = copy(s.deep()), sca0 = copy(s.sca());
        REQUIRE(step(kHbv, sun(mode), s, StepForcing{f.t, f.p, f.di, kHour * 6}, 1) ==
                StepStatus::ok);
        for (std::size_t i = 0; i < n; ++i) {
            // Scale: covers in [0, 1]; the probe's largest rise was 2.2e-16.
            if (s.sca()[i] > sca0[i] + 1e-12) {
                CAPTURE(j, i, sca0[i], s.sca()[i]);
                FAIL("I7: sca rose in a step without snowfall");
            }
            fell += s.sca()[i] < sca0[i] && s.sca()[i] > 0.0 ? 1 : 0;
            for (std::size_t k = 0; k < bins; ++k) {
                const std::size_t at = i * bins + k;
                const double was = cover(sp0[at], deep0[at]), now = cover(s.sp()[at], s.deep()[at]);
                if (now > was + 1e-12) {
                    CAPTURE(j, i, k, sp0[at], deep0[at], s.sp()[at], s.deep()[at]);
                    FAIL("I7: a bin's cover rose in a step without snowfall");
                }
            }
        }
    }
    CHECK(fell > 0);  // the cover did fall through partial values
}

TEST_CASE("P8 one bin follows the depletion curve under random melt (I8)",
          "[snow][hbv][linear][property]") {
    const auto seed = GENERATE(range(std::uint64_t{0}, std::uint64_t{16}));
    CAPTURE(seed);
    std::mt19937_64 rng{4800 + seed};
    // Deep in [1.1, 2], not (1, 2]: the curve divides by deep0 - lo0 =
    // 2 (deep - 1) swe, so near-flat bins amplify the rounding of the summed
    // melt past 1e-12. Scale: swe 10 to 1000 mm, 400 steps.
    std::uniform_real_distribution<double> deep_factor(1.1, 2.0), swe(10.0, 1000.0);
    const double d = deep_factor(rng), w = swe(rng);
    const std::array f{1.0}, h{d};
    auto s = SnowState::uniform(1, SnowDistribution{f, h}, w);
    const double deep0 = d * w, lo0 = (2.0 - d) * w;
    // Mean melt potential 1.25 deep0 / 400 per daily step: bare before the end.
    std::uniform_real_distribution<double> temperature(1e-3, 2.0 * 1.25 * deep0 / (400.0 * 3.0));
    const auto p = with(0.1, 0.0);
    double melt = 0.0;
    for (int j = 0; j < 400; ++j) {
        const double t = temperature(rng);
        step1(s, p, t, 0.0);
        melt += 3.0 * t;
        const double want = depletion(deep0, lo0, melt);
        if (std::abs(s.sca()[0] - want) > 1e-12) {
            CAPTURE(d, w, j, melt, s.sca()[0], want);
            FAIL("I8: the cover departs from the depletion curve");
        }
    }
}

// ------------------------------------------------------------- V3

TEST_CASE("V3 SnowDistribution refuses a bad deep edge", "[snow][hbv][linear][refusal]") {
    using V = std::vector<double>;
    const V f{0.5, 0.8, 1.2, 1.5};
    const auto [deep, reason] = GENERATE_COPY(table<V, std::string>({
        {V{0.5, 0.8, 1.2}, "deep"},
        {V{0.5, 0.8, 1.2, 1.5, 3.0}, "deep"},
        {V{0.5, kNaN, 1.2, 1.5}, "deep"},
        {V{0.5, 0.8, 1.2, kInf}, "deep"},
        {V{0.5, 0.8, 0.99 * 1.2, 1.5}, "a bin's deep edge must lie between its mean and twice its mean"},
        {V{0.5, 0.8, 1.2, 2.01 * 1.5}, "a bin's deep edge must lie between its mean and twice its mean"},
    }));
    CAPTURE(deep, reason);
    REQUIRE_THROWS_MATCHES((SnowDistribution{f, deep}), std::invalid_argument,
                           MessageMatches(ContainsSubstring(reason)));
}

TEST_CASE("V3 a deep edge at the mean or at twice the mean is accepted",
          "[snow][hbv][linear][refusal]") {
    const std::array f{0.0, 0.5, 1.0, 2.5};  // mean 1; a zero factor allows only deep 0
    const std::array d{0.0, 0.5, 2.0, 5.0};
    const SnowDistribution dist{f, d};
    REQUIRE(dist.deep().size() == 4);
    for (std::size_t k = 0; k < 4; ++k)
        CHECK(dist.deep()[k] == d[k]);
}

TEST_CASE("V3 deep is divided by the factors' mean, as the factors are",
          "[snow][hbv][linear][refusal]") {
    const std::array f{0.5, 1.5 + 2e-10}, d{1.0, 3.0 + 4e-10};  // mean 1 + 1e-10
    const SnowDistribution dist{f, d};
    const double mean = (f[0] + f[1]) / 2.0;
    for (std::size_t k = 0; k < 2; ++k) {
        CHECK_THAT(dist.deep()[k], WithinAbs(d[k] / mean, 1e-15));
        CHECK(dist.deep()[k] <= 2.0 * dist.factors()[k]);
    }
}

TEST_CASE("V3 an empty deep is flat: deep() equals factors()", "[snow][hbv][linear][refusal]") {
    const std::array f{0.5, 0.8, 1.2, 1.5};
    const SnowDistribution plain{f}, empty{f, std::span<const double>{}};
    for (const auto* dist : {&plain, &empty}) {
        REQUIRE(dist->deep().size() == 4);
        for (std::size_t k = 0; k < 4; ++k)
            CHECK(dist->deep()[k] == dist->factors()[k]);
    }
    const auto s = SnowState::uniform(2, plain, 100.0);
    REQUIRE(s.deep().size() == 8);
    for (std::size_t k = 0; k < 8; ++k)
        CHECK(s.deep()[k] == s.sp()[k]);
}
