// Each pixel's area once (docs/increments/137-dem-stats.md, invariant 3 and
// test 9 of section 5): over every triangle of a triangulation, the weights a
// pixel gives sum to its area inside the triangulated region; with invariants
// 2 (every row sums to 1), 5 (threads change nothing) and 8 (the span) on the
// same triangulations, and the refusals of 3.4.
//
// Interface (3.4):
//
//   #include <terrain/statistics/height_shares.hpp>
//   namespace terrain::statistics {
//   struct HeightShares { std::vector<std::uint32_t> offsets; std::vector<std::int32_t> first;
//                         std::vector<float> shares; std::vector<float> nodata; };
//   template <class T> HeightShares height_shares(std::span<const Point2> corners,
//       std::span<const double> z, const raster::RasterView<T>& dem,
//       std::span<const double> row_area, double class_m, int half, unsigned threads);
//   }
//
// THE PROPERTY: 300 random triangulations of a rectangle in index space
// (pixel (r, c) is [c, c + 1] x [r, r + 1]): a jittered lattice split into
// two triangles per quad, with off-node fans (a quad replaced by four
// triangles about a point inside it). For 20 pixels each, that pixel is
// NoData (NaN, or the view's sentinel on every other trial) and every other
// pixel 0, with every corner z 0. NoData, not a high pixel: departure (c)
// folds a high pixel's pieces centred outside into class 0. Then the sum over
// triangles of (plan area x NoData share) is the pixel's area inside the
// rectangle, |[c, c + 1] x [r, r + 1] n rectangle|, computed here as the
// product of the two interval overlaps.
//
// Bound: the NoData shares are float32, so each triangle's term carries up to
// 2^-24 of its plan area; the bound is 1e-6 times the plan area of the
// triangles the pixel touches (at least the pixel's own area inside), on
// rectangles of up to 12 x 12 pixels. A pixel outside the rectangle gives
// exactly 0 to every triangle. Every row sums to 1 within (count + 1) 2^-24;
// a triangle with data has exactly one class, index half (all residuals are
// 0), and one wholly inside the NoData pixel has count 0 and NoData 1.
//
// Degenerate geometry (tester.md 3A): points snapped onto pixel lines and
// pixel corners, fan points 1e-7 from a quad edge (slivers), sub-millimetre
// rectangles (1e-5 pixels, 0.1 mm at 10 m pixels), and rectangles at column
// 1e5 (1000 km of 10 m pixels from the origin). Duplicate vertices, collinear
// corners and non-finite input are refusals (the unit cases at the end);
// cocircular points and holes do not arise: the call takes triangles, not a
// point set or a domain.
//
// INVARIANT-CRITICAL (section 5): the overlap weights; the kill record is in
// the red commit's handback.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain/core/point.hpp>
#include <terrain/land_cover/overlap.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/view.hpp>
#include <terrain/statistics/height_shares.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using Catch::Matchers::WithinAbs;
using terrain::Point2;
using terrain::raster::RasterGeometry;
using terrain::raster::RasterView;
namespace st = terrain::statistics;

constexpr double kW = 1.0;
constexpr int kHalf = 10;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kSentinel = -9999.0;

struct Rect {
    double x0, y0, x1, y1;
};

struct Scene {
    Rect rect;
    std::vector<Point2> corners;  // 3 per triangle
    std::size_t rows, cols;
};

double shoelace(const Point2& a, const Point2& b, const Point2& c) {
    return 0.5 * std::fabs((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

double overlap(double lo, double hi, double a, double b) {
    return std::max(0.0, std::min(hi, b) - std::max(lo, a));
}

// A jittered (m + 1) x (n + 1) lattice over `rect`: corners fixed, boundary
// points moved along their edge only, interior points by up to 0.2 of a
// lattice step in each axis, which keeps every quad strictly convex. Every
// fourth interior point is snapped onto a pixel line, every ninth onto a
// pixel corner (when the snap stays within its jitter box).
Scene triangulate(std::mt19937_64& rng, Rect rect, std::size_t rows, std::size_t cols) {
    std::uniform_int_distribution<int> count(2, 6);
    const auto m = static_cast<std::size_t>(count(rng)), n = static_cast<std::size_t>(count(rng));
    const double sx = (rect.x1 - rect.x0) / static_cast<double>(m);
    const double sy = (rect.y1 - rect.y0) / static_cast<double>(n);
    std::uniform_real_distribution<double> jitter(-0.2, 0.2);
    std::vector<Point2> p((m + 1) * (n + 1));
    std::size_t counter = 0;
    for (std::size_t j = 0; j <= n; ++j)
        for (std::size_t i = 0; i <= m; ++i) {
            double x = i == m ? rect.x1 : rect.x0 + static_cast<double>(i) * sx;
            double y = j == n ? rect.y1 : rect.y0 + static_cast<double>(j) * sy;
            const bool edge_x = i == 0 || i == m, edge_y = j == 0 || j == n;
            if (!edge_x)
                x += jitter(rng) * sx;
            if (!edge_y)
                y += jitter(rng) * sy;
            if (!edge_x && !edge_y) {
                ++counter;
                const double gx = std::round(x), gy = std::round(y);
                const double bx = rect.x0 + static_cast<double>(i) * sx, by = rect.y0 + static_cast<double>(j) * sy;
                if (counter % 4 == 0 && std::fabs(gx - bx) < 0.2 * sx)
                    x = gx;
                if (counter % 9 == 0 && std::fabs(gx - bx) < 0.2 * sx && std::fabs(gy - by) < 0.2 * sy) {
                    x = gx;
                    y = gy;
                }
            }
            p[j * (m + 1) + i] = Point2{x, y};
        }
    Scene s{rect, {}, rows, cols};
    std::uniform_int_distribution<std::size_t> pick_i(0, m - 1), pick_j(0, n - 1);
    const std::size_t fan_i = pick_i(rng), fan_j = pick_j(rng);
    const std::size_t sliver_i = pick_i(rng), sliver_j = pick_j(rng);
    std::bernoulli_distribution diagonal(0.5);
    for (std::size_t j = 0; j < n; ++j)
        for (std::size_t i = 0; i < m; ++i) {
            const Point2 a = p[j * (m + 1) + i], b = p[j * (m + 1) + i + 1];
            const Point2 c = p[(j + 1) * (m + 1) + i + 1], d = p[(j + 1) * (m + 1) + i];
            const bool fan = (i == fan_i && j == fan_j) || (i == sliver_i && j == sliver_j);
            if (fan) {
                // The fan point: inside the quad, or 1e-7 of the way from the
                // quad's centre to the midpoint of its edge a-b (a sliver).
                const Point2 centre{(a.x + b.x + c.x + d.x) / 4, (a.y + b.y + c.y + d.y) / 4};
                const Point2 mid{(a.x + b.x) / 2, (a.y + b.y) / 2};
                const double t = (i == sliver_i && j == sliver_j) ? 1.0 - 1e-7 : 0.37;
                const Point2 f{centre.x + t * (mid.x - centre.x), centre.y + t * (mid.y - centre.y)};
                for (const auto& [u, v] : {std::pair{a, b}, std::pair{b, c}, std::pair{c, d}, std::pair{d, a}})
                    s.corners.insert(s.corners.end(), {u, v, f});
            } else if (diagonal(rng)) {
                s.corners.insert(s.corners.end(), {a, b, c, a, c, d});
            } else {
                s.corners.insert(s.corners.end(), {a, b, d, b, c, d});
            }
        }
    return s;
}

Rect random_rect(std::mt19937_64& rng, int kind, double col_offset) {
    std::uniform_real_distribution<double> at(0.1, 2.9), size(3.0, 12.0), tiny(0.0, 1.0);
    if (kind == 1) {
        // Sub-millimetre: 1e-5 pixels across, inside one pixel or across a line.
        const double x0 = col_offset + 3.0 + tiny(rng) * 2.0, y0 = 3.0 + tiny(rng) * 2.0;
        return {x0, y0, x0 + 1e-5, y0 + 1e-5 * (1.0 + tiny(rng))};
    }
    const double x0 = col_offset + at(rng), y0 = at(rng);
    return {x0, y0, x0 + size(rng), y0 + size(rng)};
}

}  // namespace

TEST_CASE("each pixel's area once over a triangulation, 300 cases", "[statistics][height_shares][property]") {
    std::mt19937_64 rng(137);
    std::size_t checked = 0;
    for (int trial = 0; trial < 300; ++trial) {
        const int kind = trial % 10 == 9 ? 1 : 0;
        const double offset = trial % 7 == 3 ? 1.0e5 : 0.0;
        const std::size_t rows = 16, cols = offset > 0 ? 100016 : 16;
        const Scene s = triangulate(rng, random_rect(rng, kind, offset), rows, cols);
        const std::size_t T = s.corners.size() / 3;
        const std::vector<double> z(s.corners.size(), 0.0), row_area(rows, 1.0);
        std::vector<double> area(T);
        for (std::size_t t = 0; t < T; ++t)
            area[t] = shoelace(s.corners[3 * t], s.corners[3 * t + 1], s.corners[3 * t + 2]);
        std::vector<double> dem(rows * cols, 0.0);
        const bool sentinel = trial % 2 == 1;
        const std::optional<double> nodata = sentinel ? std::optional<double>{kSentinel} : std::nullopt;
        const RasterView<double> view{RasterGeometry{0.0, 0.0, 1.0, 1.0, cols, rows}, dem.data(), nodata};
        const auto c_lo = static_cast<std::size_t>(std::max(0.0, std::floor(s.rect.x0) - 1.0));
        const auto r_lo = static_cast<std::size_t>(std::max(0.0, std::floor(s.rect.y0) - 1.0));
        const auto c_hi = std::min(cols - 1, static_cast<std::size_t>(std::ceil(s.rect.x1)));
        const auto r_hi = std::min(rows - 1, static_cast<std::size_t>(std::ceil(s.rect.y1)));
        std::uniform_int_distribution<std::size_t> pc(c_lo, c_hi), pr(r_lo, r_hi);
        for (int k = 0; k < 20; ++k) {
            const std::size_t r = pr(rng), c = pc(rng);
            dem[r * cols + c] = sentinel ? kSentinel : kNaN;
            const unsigned threads = static_cast<unsigned>(1 + k % 3);
            const st::HeightShares out = st::height_shares<double>(s.corners, z, view, row_area, kW, kHalf, threads);
            if (k == 0) {
                const st::HeightShares serial = st::height_shares<double>(s.corners, z, view, row_area, kW, kHalf, 1);
                const st::HeightShares wide = st::height_shares<double>(s.corners, z, view, row_area, kW, kHalf, 7);
                for (const st::HeightShares* other : {&serial, &wide}) {
                    REQUIRE(other->offsets == out.offsets);
                    REQUIRE(other->first == out.first);
                    REQUIRE(other->shares == out.shares);
                    REQUIRE(other->nodata == out.nodata);
                }
            }
            dem[r * cols + c] = 0.0;
            REQUIRE(out.offsets.size() == T + 1);
            REQUIRE(out.first.size() == T);
            REQUIRE(out.nodata.size() == T);
            REQUIRE(out.offsets.front() == 0u);
            REQUIRE(out.offsets.back() == out.shares.size());
            const double want = overlap(s.rect.x0, s.rect.x1, static_cast<double>(c), static_cast<double>(c) + 1.0) *
                                overlap(s.rect.y0, s.rect.y1, static_cast<double>(r), static_cast<double>(r) + 1.0);
            double got = 0.0, touched = 0.0;
            for (std::size_t t = 0; t < T; ++t) {
                const std::uint32_t count = out.offsets[t + 1] - out.offsets[t];
                INFO("trial " << trial << " triangle " << t << " count " << count << " first " << out.first[t]);
                double sum = out.nodata[t];
                for (std::uint32_t i = out.offsets[t]; i < out.offsets[t + 1]; ++i) {
                    REQUIRE(out.shares[i] >= 0.0f);
                    sum += out.shares[i];
                }
                REQUIRE(std::fabs(sum - 1.0) <= static_cast<double>(count + 1) * std::ldexp(1.0, -24));
                if (count == 0) {
                    REQUIRE(out.first[t] == 0);
                    REQUIRE(out.nodata[t] == 1.0f);
                } else {
                    REQUIRE(count == 1u);
                    REQUIRE(out.first[t] == kHalf);
                }
                got += area[t] * static_cast<double>(out.nodata[t]);
                if (out.nodata[t] > 0.0f)
                    touched += area[t];
            }
            INFO("trial " << trial << " pixel (" << r << ", " << c << ") got " << got << " want " << want);
            if (want == 0.0)
                REQUIRE(got == 0.0);
            else
                REQUIRE(std::fabs(got - want) <= 1e-6 * std::max(want, touched));
            ++checked;
        }
    }
    REQUIRE(checked == 300 * 20);
}

TEST_CASE("a whole pixel inside a triangle weighs exactly its row area", "[statistics][height_shares]") {
    // Pixel (1, 1) inside the triangle (0, 0), (4, 0), (0, 4), NoData: rows
    // 0 to 3 hold 3.5, 2.5, 1.5 and 0.5 of it, row 1 weighing 2.
    std::vector<double> dem(16, 0.0);
    dem[1 * 4 + 1] = kNaN;
    const RasterView<double> view{RasterGeometry{0.0, 0.0, 1.0, 1.0, 4, 4}, dem.data(), std::nullopt};
    const std::vector<Point2> corners{{0.0, 0.0}, {4.0, 0.0}, {0.0, 4.0}};
    const std::vector<double> z(3, 0.0), row_area{1.0, 2.0, 1.0, 1.0};
    const auto out = st::height_shares<double>(corners, z, view, row_area, kW, kHalf, 1);
    const double total = 3.5 + 2.5 * 2.0 + 1.5 + 0.5;
    // Bound: float32 rounding of a share, 6e-8.
    REQUIRE(out.offsets == std::vector<std::uint32_t>{0, 1});
    REQUIRE(out.first == std::vector<std::int32_t>{kHalf});
    REQUIRE_THAT(out.nodata[0], WithinAbs(2.0 / total, 1e-7));
    REQUIRE_THAT(out.shares[0], WithinAbs((total - 2.0) / total, 1e-7));
}

TEST_CASE("a float DEM and its sentinel", "[statistics][height_shares]") {
    // Test 4's six marked nodes in float, one of them the sentinel.
    std::vector<float> dem(36, 0.0f);
    dem[0] = 3.0f;
    dem[1] = -9999.0f;
    const RasterView<float> view{RasterGeometry{0.0, 0.0, 1.0, 1.0, 6, 6}, dem.data(), std::optional<float>{-9999.0f}};
    const std::vector<Point2> corners{{0.0, 0.0}, {6.0, 0.0}, {0.0, 6.0}};
    const std::vector<double> z(3, 0.0), row_area(6, 1.0);
    const auto out = st::height_shares<float>(corners, z, view, row_area, kW, kHalf, 1);
    REQUIRE(out.first == std::vector<std::int32_t>{kHalf});
    REQUIRE(out.offsets == std::vector<std::uint32_t>{0, 4});
    REQUIRE_THAT(out.shares[0], WithinAbs(16.0 / 18.0, 1e-7));
    REQUIRE_THAT(out.shares[3], WithinAbs(1.0 / 18.0, 1e-7));
    REQUIRE_THAT(out.nodata[0], WithinAbs(1.0 / 18.0, 1e-7));
}

TEST_CASE("refusals before any work", "[statistics][height_shares]") {
    std::vector<double> dem(16, 0.0);
    const RasterView<double> view{RasterGeometry{0.0, 0.0, 1.0, 1.0, 4, 4}, dem.data(), std::nullopt};
    const std::vector<double> rows4(4, 1.0);
    const std::vector<Point2> good{{0.5, 0.5}, {1.5, 0.5}, {1.0, 1.5}};
    const std::vector<double> z3(3, 0.0);
    const auto call = [&](const std::vector<Point2>& c, const std::vector<double>& z, const std::vector<double>& ra,
                          double class_m = kW, int half = kHalf) {
        return st::height_shares<double>(c, z, view, ra, class_m, half, 1);
    };
    REQUIRE_THROWS_AS(call({{0.5, 0.5}, {1.5, 0.5}}, {0.0, 0.0}, rows4), std::invalid_argument);
    REQUIRE_THROWS_AS(call(good, {0.0, 0.0}, rows4), std::invalid_argument);
    REQUIRE_THROWS_AS(call(good, z3, std::vector<double>(3, 1.0)), std::invalid_argument);
    const double inf = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS(call({{0.5, kNaN}, {1.5, 0.5}, {1.0, 1.5}}, z3, rows4), std::invalid_argument);
    REQUIRE_THROWS_AS(call({{inf, 0.5}, {1.5, 0.5}, {1.0, 1.5}}, z3, rows4), std::invalid_argument);
    REQUIRE_THROWS_AS(call(good, {0.0, kNaN, 0.0}, rows4), std::invalid_argument);
    REQUIRE_THROWS_AS(call(good, {0.0, -inf, 0.0}, rows4), std::invalid_argument);
    // Collinear corners, and a duplicate vertex: zero area.
    REQUIRE_THROWS_AS(call({{0.5, 0.5}, {1.5, 1.5}, {2.5, 2.5}}, z3, rows4), std::invalid_argument);
    REQUIRE_THROWS_AS(call({{0.5, 0.5}, {0.5, 0.5}, {1.0, 1.5}}, z3, rows4), std::invalid_argument);
    REQUIRE_THROWS_AS(call({{0.5, 0.5}, {1.5, 0.5}, {1.0, 4.5}}, z3, rows4), terrain::land_cover::OutsideWindow);
    for (const double w : {0.0, -1.0, kNaN, inf})
        REQUIRE_THROWS_AS(call(good, z3, rows4, w), std::invalid_argument);
    for (const int half : {-1, 128})
        REQUIRE_THROWS_AS(call(good, z3, rows4, kW, half), std::invalid_argument);
    REQUIRE_NOTHROW(call(good, z3, rows4, kW, 0));
    REQUIRE_NOTHROW(call(good, z3, rows4, kW, 127));
    // No triangles: one offset, nothing else.
    const auto empty = call({}, {}, rows4);
    REQUIRE(empty.offsets == std::vector<std::uint32_t>{0});
    REQUIRE(empty.first.empty());
    REQUIRE(empty.shares.empty());
    REQUIRE(empty.nodata.empty());
}
