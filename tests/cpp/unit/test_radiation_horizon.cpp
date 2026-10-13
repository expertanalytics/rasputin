// Increment 46 PR 5, the terrain-radiation design's TR3a: the horizon kernel
// (docs/increments/46-hbv-snow.md, section 10.8 "PR 5, TR3a", its five pins
// and its 13 mutation targets; docs/increments/terrain-radiation.md, R4, I4,
// I5, its degeneracy policy and verification 3). Invariant-critical.
//
//   Q1-Q2   the int16 table entry: 100 h rounded half away from zero (pin 5),
//           -9000 for no valid sample;
//   P1-P3   the max-pyramid: every block's maximum bounds its nodes and is
//           bounded by its nodes plus one row and column; NaN and the NoData
//           sentinel are not nodes; the top block and global_max agree;
//   F1-F7   section 10.8's mutation fixtures (heights depending on x only,
//           receiver on a node, J the identity, direction k = 16, east),
//           each with the reference h and reach of probe_review3.py, checked
//           on the exhaustive search and on the pruned one at 1 and 8
//           threads, and I4 in all 64 directions;
//   S1-S3   the far raster's first sample at D_near + delta_far (pin 3), the
//           sample at D_near itself read from the near raster, and a far
//           raster of None equal bit for bit to one never reached (pin 4);
//   T1-T6   verification 3's analytic terrains: the V-valley with c = 0 and
//           c > 0, the V-valley turned by J (gamma = +-15 degrees, and a
//           scale m), the square crater, the isolated cone at 20, 50 and
//           100 km, the flat plane;
//   D1-D5   degeneracies: a receiver off the raster (-9000, reach 0), r_i = 0,
//           a NoData sentinel as a coverage gap, no receivers, refusals;
//   I4      pruned equals exhaustive, horizon and reach, bit for bit, at 1 and
//           8 threads, on random terrains with NaN holes and 300 receivers;
//   I5      the horizon never rises when c grows, never falls when D_far grows.
//
// Kill map (section 10.8's targets; the fixture is the test that names it):
//   a < 0 branch a/d_lo for a/d_hi ................ F1
//   a >= 0 branch a/d_hi for a/d_lo ............... F2
//   curvature term c d_hi for c d_lo .............. F1, F2
//   early stop on the current raster only (pin 2)  F2
//   search starting at r_i + delta_near ........... F4
//   search starting at 0 .......................... F4
//   near raster read past D_near .................. F5, F1, F2
//   far raster read from r_i ...................... F5, F2, F4, F6, F7
//   march not ended at the first missing sample ... F6
//   reach rounded, not floored (pin 1) ............ F6
//   reach not recorded (255) ...................... F6
//   None read as zeros to the asked D_far (pin 4) . F7, S3
//   int16 truncated or floored (pin 5) ............ Q1, T1 h(45), F7
//
// Interface assumed (header include/terrain/radiation/horizon.hpp, namespace
// terrain::radiation; names chosen here and listed in the handback):
//
//   struct HorizonParams {
//       std::size_t n_directions;  // N; direction k at true azimuth k 360/N
//       double curvature;          // c = (1 - k)/(2 R), per metre, >= 0
//       double d_near;             // metres, > 0
//       double d_far;              // metres, >= d_near; ignored without a far raster
//   };
//   struct HorizonTable {
//       std::vector<std::int16_t> horizon;  // n N, receiver-major: [i N + k]
//       std::vector<std::uint8_t> reach;    // n N, whole km, 255 = nothing missing
//   };
//   std::int16_t horizon_centidegrees(double tan_h) noexcept;  // -inf: -9000
//   class MaxPyramid {
//     public:
//       explicit MaxPyramid(const raster::RasterView<float>&);
//       std::size_t levels() const noexcept;  // the last level is one block
//       // The largest valid node of block (row, col) at `level`, whose nodes
//       // are rows [row 2^level, (row + 1) 2^level) and the same columns, or
//       // that block plus the next row and column; -inf when it has none.
//       float block_max(std::size_t level, std::size_t row, std::size_t col) const noexcept;
//       float global_max() const noexcept;  // -inf when no node is valid
//   };
//   HorizonTable horizon_exhaustive(const ReceiverSet&, const raster::RasterView<float>& near,
//                                   const std::optional<raster::RasterView<float>>& far,
//                                   const HorizonParams&);
//   HorizonTable horizon_search(const ReceiverSet&, const raster::RasterView<float>& near,
//                               const std::optional<raster::RasterView<float>>& far,
//                               const HorizonParams&, unsigned threads);
//
// Both searches throw std::invalid_argument, naming the field, for N = 0, a
// curvature that is negative or not finite, a d_near that is not positive and
// finite, a far raster with a d_far below d_near or not finite, and naming
// "receiver i" for a receiver check_receivers refuses.
//
// Absolute bounds: table entries are in 0.01 degrees, so a closed form is met
// when |entry - 100 h| <= 0.5 + 1e-3; the 1e-3 covers float32 heights (about
// 1e-4 centidegrees at the 577 m ridge, the largest checked), and no closed
// form below is within 2.7e-3 of a half (probe in the commit message).

#include <catch2/catch_test_macros.hpp>

#include <terrain/radiation/horizon.hpp>
#include <terrain/radiation/receivers.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/view.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace terrain::radiation;
using terrain::raster::RasterGeometry;
using terrain::raster::RasterView;

constexpr double kDeg = std::numbers::pi / 180.0;
constexpr double kC = (1.0 - 0.13) / (2.0 * 6371000.0);  // 6.828e-8 per metre
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr float kNegInf = -std::numeric_limits<float>::infinity();
constexpr std::size_t kEast = 16;  // of 64

// A raster that owns its heights; view() is the kernel's input.
struct Grid {
    RasterGeometry g;
    std::vector<float> z;
    std::optional<float> nodata;
    RasterView<float> view() const { return RasterView<float>{g, z.data(), nodata}; }
};

// Heights depending on x only: profile[c] on every row; node (row 0, col 0)
// at (0, y_max).
Grid x_only(const std::vector<double>& profile, double delta, std::size_t rows, double y_max) {
    Grid out{RasterGeometry{0.0, y_max, delta, delta, profile.size(), rows}, {}, std::nullopt};
    out.z.reserve(profile.size() * rows);
    for (std::size_t r = 0; r < rows; ++r)
        for (const double v : profile)
            out.z.push_back(static_cast<float>(v));
    return out;
}

std::vector<double> profile(std::size_t n, double base, std::initializer_list<std::pair<std::size_t, double>> nodes) {
    std::vector<double> p(n, base);
    for (const auto& [i, v] : nodes)
        p[i] = v;
    return p;
}

// Section 10.8's rasters: the near one 10 m with 201 rows (y from 1000 to
// -1000), the far one 100 m with 41 rows (2000 to -2000); the receiver at
// (0, 0) is a node of both.
Grid near_x(const std::vector<double>& p) { return x_only(p, 10.0, 201, 1000.0); }
Grid far_x(const std::vector<double>& p) { return x_only(p, 100.0, 41, 2000.0); }

// A square grid centred on (0, 0): node (r, c) at (-half + c d, half - r d).
template <typename F>
Grid centred(double half, double d, F&& height) {
    const auto n = static_cast<std::size_t>(std::lround(2.0 * half / d)) + 1;
    Grid out{RasterGeometry{-half, half, d, d, n, n}, std::vector<float>(n * n), std::nullopt};
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            out.z[r * n + c] = static_cast<float>(height(-half + static_cast<double>(c) * d,
                                                         half - static_cast<double>(r) * d));
    return out;
}

std::array<double, 4> rot_j(double gamma_deg, double m) {
    const double c = std::cos(gamma_deg * kDeg), s = std::sin(gamma_deg * kDeg);
    return {c / m, s / m, -s / m, c / m};
}

// Receivers that own their arrays; latitude, longitude and the normal are
// valid placeholders the search does not read.
struct Receivers {
    std::vector<double> x, y, z, lat, lon, j, normal, r;
    void add(double xi, double yi, double zi, double ri, std::array<double, 4> ji = {1.0, 0.0, 0.0, 1.0}) {
        x.push_back(xi);
        y.push_back(yi);
        z.push_back(zi);
        r.push_back(ri);
        lat.push_back(61.3694);
        lon.push_back(8.4852);
        j.insert(j.end(), ji.begin(), ji.end());
        normal.insert(normal.end(), {0.0, 0.0, 1.0});
    }
    ReceiverSet set() const {
        return ReceiverSet{.x = x, .y = y, .z_obs = z, .lat_deg = lat, .lon_deg = lon,
                           .jacobian = j, .normal = normal, .exclusion_radius = r};
    }
};

Receivers one(double z_obs, double r_i, std::array<double, 4> j = {1.0, 0.0, 0.0, 1.0}) {
    Receivers out;
    out.add(0.0, 0.0, z_obs, r_i, j);
    return out;
}

std::optional<RasterView<float>> view_of(const std::optional<Grid>& g) {
    return g ? std::optional<RasterView<float>>{g->view()} : std::nullopt;
}

// The exhaustive table, after checking (I4) that the pruned search at 1 and
// 8 threads gives the same table bit for bit.
HorizonTable searched(const Receivers& rc, const Grid& near, const std::optional<Grid>& far,
                      const HorizonParams& p) {
    HorizonTable ref = horizon_exhaustive(rc.set(), near.view(), view_of(far), p);
    REQUIRE(ref.horizon.size() == rc.x.size() * p.n_directions);
    REQUIRE(ref.reach.size() == rc.x.size() * p.n_directions);
    for (const unsigned threads : {1u, 8u}) {
        INFO("pruned search, threads " << threads);
        const HorizonTable got = horizon_search(rc.set(), near.view(), view_of(far), p, threads);
        CHECK(got.horizon == ref.horizon);
        CHECK(got.reach == ref.reach);
    }
    return ref;
}

// |entry - 100 h| within half a centidegree (and the float32 margin).
bool meets(std::int16_t entry, double tan_h) {
    return std::abs(static_cast<double>(entry) - 100.0 * std::atan(tan_h) / kDeg) <= 0.5 + 1e-3;
}

HorizonParams params(double c, double d_near, double d_far, std::size_t n = 64) {
    return HorizonParams{.n_directions = n, .curvature = c, .d_near = d_near, .d_far = d_far};
}

// splitmix64: the same numbers on every standard library, unlike
// std::uniform_real_distribution.
struct Random {
    std::uint64_t s;
    double uniform(double lo, double hi) {
        std::uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        z ^= z >> 31;
        return lo + (hi - lo) * static_cast<double>(z >> 11) * 0x1.0p-53;
    }
};

// Rugged terrain: a few random bumps and ridges, noise, and NaN holes.
Grid random_terrain(Random& rnd, double half, double d, std::size_t holes, double relief) {
    std::array<std::array<double, 4>, 6> bumps{};
    for (auto& b : bumps)
        b = {rnd.uniform(-half, half), rnd.uniform(-half, half), rnd.uniform(0.1, 1.0) * relief,
             rnd.uniform(0.05, 0.4) * half};
    Grid g = centred(half, d, [&](double x, double y) {
        double z = 0.0;
        for (const auto& [bx, by, h, w] : bumps)
            z += h * std::exp(-((x - bx) * (x - bx) + (y - by) * (y - by)) / (w * w));
        return z + 0.3 * relief * std::sin(x / (7.0 * d)) * std::cos(y / (11.0 * d));
    });
    for (auto& v : g.z)
        v += static_cast<float>(rnd.uniform(-2.0, 2.0));
    for (std::size_t k = 0; k < holes; ++k)
        g.z[static_cast<std::size_t>(rnd.uniform(0.0, static_cast<double>(g.z.size())))] = static_cast<float>(kNaN);
    return g;
}

}  // namespace

// --- Q: the int16 table entry ----------------------------------------------

TEST_CASE("Q1: the entry is 100 h rounded half away from zero (pin 5)", "[radiation][horizon][quantise]") {
    const auto at = [](double h_deg) { return horizon_centidegrees(std::tan(h_deg * kDeg)); };
    CHECK(at(22.2077) == 2221);   // the V-valley's h(45); truncation gives 2220
    CHECK(at(-0.3647) == -36);    // F7; flooring gives -37
    CHECK(at(-0.3667) == -37);    // truncation gives -36
    CHECK(at(0.3647) == 36);
    CHECK(at(30.0) == 3000);
    CHECK(at(0.0) == 0);
    CHECK(at(-0.0000391) == 0);   // the flat plane at r_i = 10 m
    CHECK(at(89.9999) == 9000);
    CHECK(at(-89.9999) == -9000);
    CHECK(horizon_centidegrees(std::numeric_limits<double>::infinity()) == 9000);
}

TEST_CASE("Q2: no valid sample is -90 degrees", "[radiation][horizon][quantise]") {
    CHECK(horizon_centidegrees(-std::numeric_limits<double>::infinity()) == -9000);
}

// --- P: the max-pyramid ------------------------------------------------------

namespace {

// Largest valid node in rows [r0, r1) x cols [c0, c1), clipped; -inf if none.
float max_over(const Grid& g, std::size_t r0, std::size_t r1, std::size_t c0, std::size_t c1) {
    const auto v = g.view();
    float m = kNegInf;
    for (std::size_t r = r0; r < std::min(r1, g.g.rows()); ++r)
        for (std::size_t c = c0; c < std::min(c1, g.g.cols()); ++c)
            if (!v.is_nodata({r, c}))
                m = std::max(m, v.value_at({r, c}));
    return m;
}

void check_pyramid(const Grid& g) {
    const MaxPyramid p{g.view()};
    REQUIRE(p.levels() >= 1);
    CHECK(p.global_max() == max_over(g, 0, g.g.rows(), 0, g.g.cols()));
    CHECK(p.block_max(p.levels() - 1, 0, 0) == p.global_max());
    for (std::size_t level = 0; level < p.levels(); ++level) {
        const std::size_t s = std::size_t{1} << level;
        for (std::size_t r = 0; r * s < g.g.rows(); ++r)
            for (std::size_t c = 0; c * s < g.g.cols(); ++c) {
                INFO("level " << level << " block " << r << ", " << c);
                const float m = p.block_max(level, r, c);
                CHECK(m >= max_over(g, r * s, (r + 1) * s, c * s, (c + 1) * s));
                CHECK(m <= max_over(g, r * s, (r + 1) * s + 1, c * s, (c + 1) * s + 1));
            }
    }
}

}  // namespace

TEST_CASE("P1: every block's maximum bounds its nodes, on odd and even shapes", "[radiation][horizon][pyramid]") {
    Random rnd{11};
    for (const auto& [rows, cols] : std::array<std::array<std::size_t, 2>, 5>{{{1, 1}, {1, 7}, {5, 3}, {16, 16}, {37, 23}}}) {
        Grid g{RasterGeometry{0.0, 0.0, 1.0, 1.0, cols, rows}, std::vector<float>(rows * cols), std::nullopt};
        for (auto& v : g.z)
            v = static_cast<float>(rnd.uniform(-500.0, 2500.0));
        INFO(rows << " x " << cols);
        check_pyramid(g);
    }
}

TEST_CASE("P2: NaN and the NoData sentinel are not nodes of the pyramid", "[radiation][horizon][pyramid]") {
    Grid g{RasterGeometry{0.0, 0.0, 1.0, 1.0, 9, 6}, std::vector<float>(54, 100.0f), -9999.0f};
    g.z[13] = static_cast<float>(kNaN);
    g.z[40] = -9999.0f;
    g.z[7] = 300.0f;
    check_pyramid(g);
    CHECK(MaxPyramid{g.view()}.global_max() == 300.0f);

    Grid high_sentinel{RasterGeometry{0.0, 0.0, 1.0, 1.0, 4, 4}, std::vector<float>(16, 5.0f), 32767.0f};
    high_sentinel.z[5] = 32767.0f;
    CHECK(MaxPyramid{high_sentinel.view()}.global_max() == 5.0f);
}

TEST_CASE("P3: a raster with no valid node has a maximum of -inf", "[radiation][horizon][pyramid]") {
    const Grid g{RasterGeometry{0.0, 0.0, 1.0, 1.0, 3, 2}, std::vector<float>(6, static_cast<float>(kNaN)), std::nullopt};
    const MaxPyramid p{g.view()};
    CHECK(p.global_max() == kNegInf);
    CHECK(p.block_max(p.levels() - 1, 0, 0) == kNegInf);
}

// --- F: section 10.8's mutation fixtures -------------------------------------

TEST_CASE("F1: a receiver above all terrain (the a < 0 bound, curvature at d_lo)", "[radiation][horizon][mutation]") {
    const Grid near = near_x(profile(2101, 0.0, {{1000, 2000.0}}));
    const Grid far = far_x(profile(2002, 1000.0, {}));
    const HorizonTable t = searched(one(3000.0, 10.0), near, far, params(kC, 20000.0, 200000.0));
    CHECK(t.horizon[kEast] == -134);  // -1.3388 degrees
    CHECK(t.reach[kEast] == 255);
}

TEST_CASE("F2: a far peak at 25 km (the a >= 0 bound, the early stop over both rasters)",
          "[radiation][horizon][mutation]") {
    const Grid near = near_x(profile(2101, 0.0, {}));
    const Grid far = far_x(profile(2002, 0.0, {{250, 100.0}}));
    const HorizonTable t = searched(one(0.0, 10.0), near, far, params(kC, 20000.0, 200000.0));
    CHECK(t.horizon[kEast] == 13);  // 0.1314 degrees
    CHECK(t.reach[kEast] == 255);
}

TEST_CASE("F4: a spike at r_i (the search starts at r_i, not at 0 or a step later)",
          "[radiation][horizon][mutation]") {
    const Grid near = near_x(profile(401, 0.0, {{1, 100.0}}));
    const HorizonTable t = searched(one(20.0, 7.5), near, std::nullopt, params(0.0, 2000.0, 2000.0));
    CHECK(t.horizon[kEast] == 8223);  // 82.2348 degrees: (75 - 20)/7.5 at the first sample
    CHECK(t.reach[kEast] == 255);
}

TEST_CASE("F5: near and far rasters disagree (the switch at D_near)", "[radiation][horizon][mutation]") {
    // Near: a 1000 m node at 3 km, past D_near. Far: a 1000 m node at 1.5 km,
    // inside D_near, and a 400 m node at 4 km. Only the last is searched.
    const Grid near = near_x(profile(401, 0.0, {{300, 1000.0}}));
    const Grid far = far_x(profile(61, 0.0, {{15, 1000.0}, {40, 400.0}}));
    const HorizonTable t = searched(one(0.0, 10.0), near, far, params(0.0, 2000.0, 5000.0));
    CHECK(t.horizon[kEast] == 571);  // 5.7106 degrees = atan(400/4000)
    CHECK(t.reach[kEast] == 255);
}

TEST_CASE("F6: a coverage gap at 3.7 km ends the march, reach floored (pin 1)", "[radiation][horizon][mutation]") {
    const Grid near = near_x(profile(2101, 0.0, {{370, kNaN}, {500, 1000.0}}));
    const HorizonTable t = searched(one(0.0, 10.0), near, std::nullopt, params(0.0, 20000.0, 20000.0));
    CHECK(t.horizon[kEast] == 0);  // the 1000 m node at 5 km lies behind the gap
    CHECK(t.reach[kEast] == 3);    // the first missing sample near 3.69 km
}

TEST_CASE("F7: a far raster of None ends the search at D_near (pin 4)", "[radiation][horizon][mutation]") {
    const Grid near = near_x(profile(2101, 0.0, {}));
    const HorizonTable t = searched(one(100.0, 10.0), near, std::nullopt, params(kC, 20000.0, 60000.0));
    CHECK(t.horizon[kEast] == -36);  // -0.3647 degrees: -100/20000 - c 20000 at D_near
    CHECK(t.reach[kEast] == 255);
}

// --- S: the switch between rasters ------------------------------------------

TEST_CASE("S1: the far raster's first sample is at D_near + delta_far (pin 3)", "[radiation][horizon][switch]") {
    // A 100 m far node at 2100 m = D_near + 2 delta_far. Samples counted from
    // r_i (2060, 2110) or anywhere off the node read less of it.
    const Grid near = near_x(profile(401, 0.0, {}));
    const Grid far = far_x(profile(61, 0.0, {{21, 100.0}}));
    const HorizonTable t = searched(one(0.0, 10.0), near, far, params(0.0, 2000.0, 5000.0));
    CHECK(t.horizon[kEast] == 273);  // atan(100/2100) = 2.7263 degrees
}

TEST_CASE("S2: the sample at D_near is the near raster's", "[radiation][horizon][switch]") {
    // A 100 m near node at exactly D_near = 2000 m; the far raster is flat there.
    const Grid near = near_x(profile(401, 0.0, {{200, 100.0}}));
    const Grid far = far_x(profile(61, 0.0, {}));
    const HorizonTable t = searched(one(0.0, 10.0), near, far, params(0.0, 2000.0, 5000.0));
    CHECK(t.horizon[kEast] == 286);  // atan(100/2000) = 2.8624 degrees
}

TEST_CASE("S3: a far raster of None equals, bit for bit, one never reached", "[radiation][horizon][switch]") {
    Random rnd{3};
    const Grid near = random_terrain(rnd, 1500.0, 10.0, 0, 300.0);
    const Grid far = random_terrain(rnd, 5000.0, 100.0, 0, 3000.0);
    Receivers rc;
    for (int i = 0; i < 5; ++i)
        rc.add(rnd.uniform(-500.0, 500.0), rnd.uniform(-500.0, 500.0), rnd.uniform(0.0, 200.0), 10.0);
    const HorizonTable none = searched(rc, near, std::nullopt, params(kC, 1000.0, 4000.0));
    const HorizonTable unreached = searched(rc, near, far, params(kC, 1000.0, 1000.0));
    CHECK(none.horizon == unreached.horizon);
    CHECK(none.reach == unreached.reach);
    // ... and a far raster, once reached, is read: a 5000 m plateau beyond
    // D_near stands above everything the near raster holds.
    const Grid plateau = centred(5000.0, 100.0, [](double, double) { return 5000.0; });
    CHECK(searched(rc, near, plateau, params(kC, 1000.0, 4000.0)).horizon != none.horizon);
}

// --- T: verification 3's analytic terrains ----------------------------------

namespace {

const double kTan30 = std::tan(30.0 * kDeg);

// The V-valley: axis on the grid line x = 0, slopes of 30 degrees, ridges at
// half-width 1000 m and flat beyond, on a 10 m grid to +-1500 m.
Grid v_valley() {
    return centred(1500.0, 10.0, [](double x, double) { return std::min(std::abs(x), 1000.0) * kTan30; });
}

}  // namespace

TEST_CASE("T1: the V-valley with c = 0: tan h = |sin A| tan 30", "[radiation][horizon][analytic]") {
    const HorizonTable t = searched(one(0.0, 5.0), v_valley(), std::nullopt, params(0.0, 1200.0, 1200.0));
    for (std::size_t k = 0; k < 64; ++k) {
        INFO("k " << k << " entry " << t.horizon[k]);
        CHECK(meets(t.horizon[k], std::abs(std::sin(static_cast<double>(k) * 5.625 * kDeg)) * kTan30));
        CHECK(t.reach[k] == 255);
    }
    CHECK(t.horizon[8] == 2221);  // h(45) = 22.2077 degrees
    CHECK(t.horizon[16] == 3000);
    CHECK(t.horizon[0] == 0);
    CHECK(t.horizon[32] == 0);
}

TEST_CASE("T2: the V-valley with curvature: tan h = |sin A| tan 30 - c r_i", "[radiation][horizon][analytic]") {
    const HorizonTable t = searched(one(0.0, 5.0), v_valley(), std::nullopt, params(kC, 1200.0, 1200.0));
    for (std::size_t k = 0; k < 64; ++k) {
        INFO("k " << k << " entry " << t.horizon[k]);
        CHECK(meets(t.horizon[k], std::abs(std::sin(static_cast<double>(k) * 5.625 * kDeg)) * kTan30 - 5.0 * kC));
    }
    CHECK(t.horizon[0] == 0);  // -0.000020 degrees
}

TEST_CASE("T3: the V-valley turned by J: rays are J^-1 (sin A, cos A), one ground metre per unit",
          "[radiation][horizon][analytic]") {
    // A ray's frame step u = J^-1 (sin A, cos A); the valley rises |u_x| tan 30
    // per ground metre, so tan h = |u_x| tan 30 whatever J is.
    for (const auto& [gamma, m] : std::array<std::array<double, 2>, 3>{{{15.0, 1.0}, {-15.0, 1.0}, {15.0, 0.9996}}}) {
        INFO("gamma " << gamma << " m " << m);
        const auto j = rot_j(gamma, m);
        const HorizonTable t = searched(one(0.0, 5.0, j), v_valley(), std::nullopt, params(0.0, 1200.0, 1200.0));
        for (std::size_t k = 0; k < 64; ++k) {
            const auto u = frame_step(j, static_cast<double>(k) * 5.625);
            INFO("k " << k << " entry " << t.horizon[k]);
            CHECK(meets(t.horizon[k], std::abs(u[0]) * kTan30));
        }
    }
    // The sign of gamma shows (frame_step is the convention TR2 pinned).
    CHECK(searched(one(0.0, 5.0, rot_j(15.0, 1.0)), v_valley(), std::nullopt, params(0.0, 1200.0, 1200.0)).horizon[8] == 1610);
    CHECK(searched(one(0.0, 5.0, rot_j(-15.0, 1.0)), v_valley(), std::nullopt, params(0.0, 1200.0, 1200.0)).horizon[8] == 2657);
    CHECK(searched(one(0.0, 5.0, rot_j(15.0, 0.9996)), v_valley(), std::nullopt, params(0.0, 1200.0, 1200.0)).horizon[16] == 2914);
}

TEST_CASE("T4: the square crater toward its sides: tan h = (H - c d^2)/d", "[radiation][horizon][analytic]") {
    // Floor to W = 1 km, ramp to H = 500 m at W + L = 2 km, flat beyond; the
    // samples 5 + 5 j hit the rim at 2000 m along the axes.
    const Grid crater = centred(2500.0, 10.0, [](double x, double y) {
        return std::clamp((std::max(std::abs(x), std::abs(y)) - 1000.0) / 1000.0, 0.0, 1.0) * 500.0;
    });
    for (const double c : {0.0, kC}) {
        INFO("c " << c);
        const HorizonTable t = searched(one(0.0, 5.0), crater, std::nullopt, params(c, 2400.0, 2400.0));
        for (const std::size_t k : {0u, 16u, 32u, 48u}) {
            INFO("k " << k);
            CHECK(meets(t.horizon[k], (500.0 - c * 2000.0 * 2000.0) / 2000.0));
        }
        CHECK(t.horizon[16] == (c == 0.0 ? 1404 : 1403));
        CHECK(std::all_of(t.reach.begin(), t.reach.end(), [](std::uint8_t r) { return r == 255; }));
    }
}

TEST_CASE("T5: an isolated cone at 20, 50 and 100 km: tan h = (H - c r0^2)/r0", "[radiation][horizon][analytic]") {
    // H = 1000 m, base radius 1 km, apex on a node at a sample distance: on the
    // near raster at 20 km (10 + 5 j), on the far one at 50 and 100 km
    // (20 000 + 50 j). The drops c r0^2 are 27.3, 170.7 and 682.8 m.
    const auto cone = [](double r0, double d, std::size_t n) {
        std::vector<double> p(n);
        for (std::size_t i = 0; i < n; ++i)
            p[i] = std::max(0.0, 1000.0 - std::abs(static_cast<double>(i) * d - r0));
        return p;
    };
    const HorizonTable t20 = searched(one(0.0, 10.0), near_x(cone(20000.0, 10.0, 2101)), std::nullopt,
                                      params(kC, 20000.0, 20000.0));
    CHECK(t20.horizon[kEast] == 278);  // 2.7844 degrees
    for (const auto& [r0, expected] : std::array<std::pair<double, int>, 2>{{{50000.0, 95}, {100000.0, 18}}}) {
        INFO("r0 " << r0);
        const HorizonTable t = searched(one(0.0, 10.0), near_x(profile(2101, 0.0, {})), far_x(cone(r0, 100.0, 1102)),
                                        params(kC, 20000.0, 110000.0));
        CHECK(meets(t.horizon[kEast], (1000.0 - kC * r0 * r0) / r0));
        CHECK(t.horizon[kEast] == expected);  // 0.9502 and 0.1818 degrees
        CHECK(t.reach[kEast] == 255);
    }
}

TEST_CASE("T6: the flat plane: h = -atan(c r_i), stored as 0", "[radiation][horizon][analytic]") {
    const Grid flat = centred(1500.0, 10.0, [](double, double) { return 0.0; });
    const HorizonTable t = searched(one(0.0, 10.0), flat, std::nullopt, params(kC, 1200.0, 1200.0));
    CHECK(std::all_of(t.horizon.begin(), t.horizon.end(), [](std::int16_t h) { return h == 0; }));
    CHECK(std::all_of(t.reach.begin(), t.reach.end(), [](std::uint8_t r) { return r == 255; }));
}

// --- D: degeneracies ---------------------------------------------------------

TEST_CASE("D1: a receiver off the raster has no sample: -9000 and reach 0", "[radiation][horizon][degenerate]") {
    const Grid flat = centred(500.0, 10.0, [](double, double) { return 0.0; });
    Receivers rc;
    rc.add(5000.0, 5000.0, 0.0, 10.0);
    const HorizonTable t = searched(rc, flat, std::nullopt, params(kC, 1000.0, 1000.0));
    CHECK(std::all_of(t.horizon.begin(), t.horizon.end(), [](std::int16_t h) { return h == -9000; }));
    CHECK(std::all_of(t.reach.begin(), t.reach.end(), [](std::uint8_t r) { return r == 0; }));
}

TEST_CASE("D2: r_i = 0 takes no height at distance 0", "[radiation][horizon][degenerate]") {
    // The receiver 1 m below a flat plane: the first sample with a distance is
    // at 5 m, tan h = 1/5.
    const Grid flat = centred(500.0, 10.0, [](double, double) { return 0.0; });
    const HorizonTable t = searched(one(-1.0, 0.0), flat, std::nullopt, params(0.0, 400.0, 400.0));
    CHECK(std::all_of(t.horizon.begin(), t.horizon.end(), [](std::int16_t h) { return h == 1131; }));  // 11.3099
    CHECK(std::all_of(t.reach.begin(), t.reach.end(), [](std::uint8_t r) { return r == 255; }));
}

TEST_CASE("D3: the NoData sentinel is a coverage gap, as NaN is", "[radiation][horizon][degenerate]") {
    Grid near = near_x(profile(2101, 0.0, {{370, -9999.0}, {500, 1000.0}}));
    near.nodata = -9999.0f;
    const HorizonTable t = searched(one(0.0, 10.0), near, std::nullopt, params(0.0, 20000.0, 20000.0));
    CHECK(t.horizon[kEast] == 0);
    CHECK(t.reach[kEast] == 3);
}

TEST_CASE("D4: no receivers, and one direction", "[radiation][horizon][degenerate]") {
    const Grid flat = centred(500.0, 10.0, [](double, double) { return 0.0; });
    const HorizonTable none = searched(Receivers{}, flat, std::nullopt, params(kC, 400.0, 400.0));
    CHECK(none.horizon.empty());
    CHECK(none.reach.empty());
    const Grid wall = centred(500.0, 10.0, [](double, double y) { return y >= 100.0 ? 100.0 : 0.0; });
    const HorizonTable north = searched(one(0.0, 10.0), wall, std::nullopt, params(0.0, 400.0, 400.0, 1));
    REQUIRE(north.horizon.size() == 1);
    CHECK(north.horizon[0] == 4500);  // atan(100/100): the wall's first node, at y = 100 m, d = 100 m
}

TEST_CASE("D5: refusals name what is wrong", "[radiation][horizon][refusal]") {
    const Grid flat = centred(500.0, 10.0, [](double, double) { return 0.0; });
    const auto refused = [&](const Receivers& rc, const std::optional<Grid>& far, const HorizonParams& p,
                             const std::string& word) {
        INFO("expects '" << word << "'");
        for (const int path : {0, 1}) {
            try {
                if (path == 0)
                    (void)horizon_exhaustive(rc.set(), flat.view(), view_of(far), p);
                else
                    (void)horizon_search(rc.set(), flat.view(), view_of(far), p, 1);
                FAIL("not refused");
            } catch (const std::invalid_argument& e) {
                CHECK(std::string{e.what()}.find(word) != std::string::npos);
            }
        }
    };
    const Receivers ok = one(0.0, 10.0);
    const double inf = std::numeric_limits<double>::infinity();
    refused(ok, std::nullopt, params(kC, 400.0, 400.0, 0), "n_directions");
    refused(ok, std::nullopt, params(-1e-9, 400.0, 400.0), "curvature");
    refused(ok, std::nullopt, params(kNaN, 400.0, 400.0), "curvature");
    refused(ok, std::nullopt, params(kC, 0.0, 400.0), "d_near");
    refused(ok, std::nullopt, params(kC, inf, inf), "d_near");
    refused(ok, flat, params(kC, 400.0, 300.0), "d_far");
    refused(ok, flat, params(kC, 400.0, kNaN), "d_far");
    Receivers bad;
    bad.add(0.0, 0.0, 0.0, 10.0);
    bad.add(0.0, 0.0, kNaN, 10.0);
    refused(bad, std::nullopt, params(kC, 400.0, 400.0), "receiver 1");
    Receivers negative = one(0.0, -1.0);
    refused(negative, std::nullopt, params(kC, 400.0, 400.0), "receiver 0");
}

// --- I4, I5: properties on random terrains ------------------------------------

namespace {

struct Scene {
    Grid near, far;
    Receivers rc;
};

// 300 receivers (more than for_each_block runs inline) over rugged near and
// far rasters with NaN holes; heights from below the terrain to above all of
// it, exclusion radii 0 to 25 m, J turned by up to 15 degrees with scales
// 0.9996 to 1.0004.
Scene random_scene(std::uint64_t seed) {
    Random rnd{seed};
    Scene s{random_terrain(rnd, 2500.0, 10.0, 60, 400.0), random_terrain(rnd, 5000.0, 100.0, 20, 1500.0), {}};
    for (int i = 0; i < 300; ++i) {
        const double z = i % 10 == 0 ? rnd.uniform(2000.0, 4000.0) : rnd.uniform(-100.0, 600.0);
        const double r_i = i % 7 == 0 ? 0.0 : rnd.uniform(0.0, 25.0);
        s.rc.add(rnd.uniform(-1200.0, 1200.0), rnd.uniform(-1200.0, 1200.0), z, r_i,
                 rot_j(rnd.uniform(-15.0, 15.0), rnd.uniform(0.9996, 1.0004)));
    }
    return s;
}

}  // namespace

TEST_CASE("I4: the pruned search equals the exhaustive one, bit for bit, at 1 and 8 threads",
          "[radiation][horizon][invariant]") {
    for (const std::uint64_t seed : {1u, 2u}) {
        const Scene s = random_scene(seed);
        for (const double c : {0.0, kC, 1e-4}) {
            INFO("seed " << seed << " c " << c);
            (void)searched(s.rc, s.near, s.far, params(c, 1000.0, 4500.0));
            (void)searched(s.rc, s.near, std::nullopt, params(c, 1000.0, 4500.0));
            (void)searched(s.rc, s.near, s.far, params(c, 700.0, 6000.0, 7));
        }
    }
}

TEST_CASE("I5: the horizon never rises as c grows, never falls as D_far grows", "[radiation][horizon][invariant]") {
    const Scene s = random_scene(5);
    const HorizonTable flat_earth = searched(s.rc, s.near, s.far, params(0.0, 1000.0, 4500.0));
    const HorizonTable curved = searched(s.rc, s.near, s.far, params(1e-5, 1000.0, 4500.0));
    const HorizonTable shorter = searched(s.rc, s.near, s.far, params(1e-5, 1000.0, 2500.0));
    std::size_t rose = 0, fell = 0;
    for (std::size_t i = 0; i < curved.horizon.size(); ++i) {
        rose += curved.horizon[i] > flat_earth.horizon[i];
        fell += shorter.horizon[i] > curved.horizon[i];
    }
    CHECK(rose == 0);
    CHECK(fell == 0);
    CHECK(curved.horizon != flat_earth.horizon);  // the fixture moves something
}
