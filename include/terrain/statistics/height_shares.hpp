#pragma once

// Each triangle's share of area per height class about its own plane.
//
// Everything is in the DEM's index space, as in land_cover/overlap.hpp: x is
// the column and y the row, pixel (r, c) is the unit square [c, c + 1] x
// [r, r + 1] and node (r, c) its centre. Per triangle, overlap.hpp's clip
// gives each pixel's piece; a piece weighs its index area times its row's
// pixel area. Its residual is the pixel's height minus the triangle's plane
// (through its three corners with their z) at the pixel centre, for every
// piece of the pixel: that is the point a refinement tolerance bounds.
// Class k of width w holds (k - 1/2) w <= residual < (k + 1/2) w, clamped to
// -half ... half and stored as index k + half.
//
// A triangle's classes run from the lowest to the highest class of the valid
// nodes it holds (each edge function over its edge's length >= -1e-9, index
// units: inside or on an edge, with room for the rounding of a corner placed
// on a node). The pieces of pixels centred outside, whose residuals
// extrapolate the plane, fold into the nearer end of that span. A triangle
// holding no valid node has one class, that of its mean residual weighted by
// piece weight: it holds no node a tolerance checked. NoData pixels (the
// view's sentinel, or NaN) weigh on the NoData share and on no class; every
// share is over the whole weight, NoData included.
//
// Triangles are independent: blocks run in parallel and rows are joined in
// triangle order, so the result does not depend on the thread count. Pure: no
// I/O, no Python.

#include <terrain/core/point.hpp>
#include <terrain/land_cover/overlap.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/raster/view.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <thread>
#include <vector>

namespace terrain::statistics {

// A sparse table, one row per triangle: triangle t's shares are
// [offsets[t], offsets[t + 1]), of class indices first[t], first[t] + 1, ...;
// first[t] is 0 for a row without shares. nodata[t] is its NoData share.
struct HeightShares {
    std::vector<std::uint32_t> offsets;
    std::vector<std::int32_t> first;
    std::vector<float> shares;
    std::vector<float> nodata;
};

// The largest half count: a .ply list of shares holds at most 255.
constexpr int kMaxHalf = 127;

namespace detail {

// The class of residual `rho`, exact against the edges (j + 1/2) w: the
// quotient's floor, corrected by the two edge products, then clamped.
[[nodiscard]] inline int class_of(double rho, double w, int half) noexcept {
    const double q = std::floor(rho / w + 0.5);
    if (!(q <= half + 1.0))
        return half;
    if (q < -half - 1.0)
        return -half;
    int k = static_cast<int>(q);
    while (rho < (k - 0.5) * w)
        --k;
    while (rho >= (k + 0.5) * w)
        ++k;
    return std::clamp(k, -half, half);
}

}  // namespace detail

// The height shares of every triangle (three corners in `corners`, their
// heights in `z`) over `dem`, `row_area` the area of one pixel of each row.
// Refusals, all before any work and for the whole batch: corners or z not
// three per triangle, row_area not one per row, a non-finite corner or z,
// class_m not finite and above 0, half outside 0 ... kMaxHalf, a zero-area
// triangle (std::invalid_argument), a corner past the closed window
// (land_cover::OutsideWindow). threads 0 means every core.
template <class T>
[[nodiscard]] HeightShares height_shares(std::span<const Point2> corners, std::span<const double> z,
                                         const raster::RasterView<T>& dem, std::span<const double> row_area,
                                         double class_m, int half, unsigned threads) {
    const std::size_t rows = dem.geometry().rows(), cols = dem.geometry().cols();
    if (z.size() != corners.size())
        throw std::invalid_argument("height_shares: z must hold three heights per triangle");
    if (row_area.size() != rows)
        throw std::invalid_argument("height_shares: row_area must hold one pixel area per row");
    if (!std::isfinite(class_m) || !(class_m > 0.0))
        throw std::invalid_argument("height_shares: class_m must be finite and above 0");
    if (half < 0 || half > kMaxHalf)
        throw std::invalid_argument("height_shares: half must lie in 0 ... 127");
    for (const double h : z)
        if (!std::isfinite(h))
            throw std::invalid_argument("height_shares: a corner height is not finite");
    land_cover::detail::check_corners(corners, rows, cols, "height_shares");

    const std::size_t n = corners.size() / 3;
    const auto classes = static_cast<std::size_t>(2 * half + 1);
    HeightShares result;
    result.first.resize(n);
    result.nodata.resize(n);
    std::vector<std::uint32_t> counts(n);
    const unsigned workers = threads == 0 ? std::max(1u, std::thread::hardware_concurrency()) : threads;
    const std::size_t block = parallel_util::default_block(n, workers);
    std::vector<std::vector<float>> per_block(n / block + (n % block != 0 ? 1 : 0));
    const parallel_util::BlockSchedule schedule{};
    parallel_util::for_each_block(n, threads, {block, schedule.inline_below}, [&](std::size_t begin, std::size_t end) {
        std::vector<double> weight(classes);
        std::vector<float>& out = per_block[begin / block];
        for (std::size_t t = begin; t < end; ++t) {
            const auto at = land_cover::detail::canonical_order(corners, t);
            const std::array<Point2, 3> tri{corners[at[0]], corners[at[1]], corners[at[2]]};
            const Point2 d1 = tri[1] - tri[0], d2 = tri[2] - tri[0];
            const double dz1 = z[at[1]] - z[at[0]], dz2 = z[at[2]] - z[at[0]], det = cross(d1, d2);
            const double gx = (dz1 * d2.y - dz2 * d1.y) / det, gy = (d1.x * dz2 - d2.x * dz1) / det;
            const auto plane = [&](const Point2& p) { return z[at[0]] + gx * (p.x - tri[0].x) + gy * (p.y - tri[0].y); };
            std::array<double, 3> length{};
            for (std::size_t e = 0; e < 3; ++e)
                length[e] = std::hypot(tri[(e + 1) % 3].x - tri[e].x, tri[(e + 1) % 3].y - tri[e].y);
            const auto held = [&](const Point2& p) {
                for (std::size_t e = 0; e < 3; ++e)
                    if (cross(tri[(e + 1) % 3] - tri[e], p - tri[e]) / length[e] < -1e-9)
                        return false;
                return true;
            };
            std::fill(weight.begin(), weight.end(), 0.0);
            double total = 0.0, missing = 0.0, valid = 0.0, moment = 0.0;
            int lo = half + 1, hi = -half - 1;
            land_cover::detail::for_each_piece(tri, rows, cols, [&](std::size_t r, std::size_t c, const land_cover::detail::Polygon* piece) {
                const double area = land_cover::detail::piece_area(piece);
                if (!(area > 0.0))
                    return;
                const double w = area * row_area[r];
                total += w;
                if (dem.is_nodata(raster::CellIndex{r, c})) {
                    missing += w;
                    return;
                }
                const Point2 node{static_cast<double>(c) + 0.5, static_cast<double>(r) + 0.5};
                const double rho = static_cast<double>(dem.row(r)[c]) - plane(node);
                const int k = detail::class_of(rho, class_m, half);
                weight[static_cast<std::size_t>(k + half)] += w;
                valid += w;
                moment += w * rho;
                if (piece == nullptr || held(node)) {
                    lo = std::min(lo, k);
                    hi = std::max(hi, k);
                }
            });
            if (lo > hi && valid > 0.0) {
                lo = hi = detail::class_of(moment / valid, class_m, half);
                std::fill(weight.begin(), weight.end(), 0.0);
                weight[static_cast<std::size_t>(lo + half)] = valid;
            }
            result.nodata[t] = static_cast<float>(missing / total);
            if (lo > hi)
                continue;  // no valid piece: count 0, first 0
            const auto a = static_cast<std::size_t>(lo + half), b = static_cast<std::size_t>(hi + half);
            for (std::size_t k = 0; k < a; ++k)
                weight[a] += weight[k];
            for (std::size_t k = b + 1; k < classes; ++k)
                weight[b] += weight[k];
            result.first[t] = static_cast<std::int32_t>(a);
            counts[t] = static_cast<std::uint32_t>(b - a + 1);
            for (std::size_t k = a; k <= b; ++k)
                out.push_back(static_cast<float>(weight[k] / total));
        }
    });

    result.offsets.resize(n + 1);
    std::uint64_t sum = 0;
    for (std::size_t t = 0; t < n; ++t) {
        result.offsets[t] = static_cast<std::uint32_t>(sum);
        sum += counts[t];
        if (sum > std::numeric_limits<std::uint32_t>::max())
            throw std::length_error("height_shares: more than 2^32 - 1 shares");
    }
    result.offsets[n] = static_cast<std::uint32_t>(sum);
    result.shares.reserve(sum);
    for (const std::vector<float>& shares : per_block)
        result.shares.insert(result.shares.end(), shares.begin(), shares.end());
    return result;
}

}  // namespace terrain::statistics
