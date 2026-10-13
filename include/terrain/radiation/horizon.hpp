// The horizon of each receiver in N true azimuths, searched on the DEM with
// Earth curvature and terrestrial refraction (docs/increments/terrain-radiation.md,
// R4; increment 46, section 10.8's five pins).
//
// The reference definition, which horizon_exhaustive evaluates sample by
// sample: receiver i looks along true azimuth A_k = k 360/N, whose frame step
// u = J_i^-1 (sin A_k, cos A_k) is one ground metre. Near samples lie at
// d = r_i + j delta_near while d <= D_near (delta_near half the near cell),
// read bilinearly from the near raster; a sample at d = 0 has no height and is
// not taken. Far samples lie at D_near + j delta_far, j >= 1, while d <= D_far,
// on the far raster (pin 3); without one, D_far is D_near (pin 4). The march
// ends at the first sample off its raster or touching a NoData node, and
// `reach` is that sample's distance floored to whole km, 255 when none is
// missing (pin 1). tan h = max (z - z_obs - c d^2) / d over the samples before
// it, -inf when there are none, and the entry is 100 h in degrees rounded half
// away from zero (pin 5), -9000 for -inf.
//
// horizon_search gives the same table bit for bit (I4): it finds `reach` first
// from coverage alone, skipping whole pyramid blocks with no NoData node, then
// marches the heights, skipping a stretch inside a block whose maximum Z cannot
// beat the best tangent so far, and stopping when the global maxima of the
// rasters still ahead cannot (pin 2). Every sample it takes is computed by the
// same functions as in the reference, so only samples that cannot raise the
// maximum are left out.
//
// Locality: each ray reads only the rasters and their pyramids, both
// read-only; receivers are independent, so the parallel-for over them writes
// disjoint slots and the table does not depend on the thread count.

#pragma once

#include <terrain/parallel_util/chunks.hpp>
#include <terrain/radiation/receivers.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/raster/view.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace terrain::radiation {

struct HorizonParams {
    std::size_t n_directions;  // N; direction k at true azimuth k 360/N
    double curvature;          // c = (1 - k)/(2 R), per metre, >= 0
    double d_near;             // metres, > 0
    double d_far;              // metres, >= d_near; ignored without a far raster
};

struct HorizonTable {
    std::vector<std::int16_t> horizon;  // n N, receiver-major: [i N + k], 0.01 degrees
    std::vector<std::uint8_t> reach;    // n N, whole km, 255 when nothing was missing
};

inline std::int16_t horizon_centidegrees(double tan_h) noexcept {
    if (tan_h == -std::numeric_limits<double>::infinity())
        return -9000;
    return static_cast<std::int16_t>(std::round(100.0 * std::atan(tan_h) * 180.0 / std::numbers::pi));
}

// Maxima over 2^L x 2^L node blocks. Block (row, col) at level L covers the
// nodes of rows [row 2^L, (row + 1) 2^L] and the same columns, both ends
// included, so it holds every node a bilinear sample in one of its cells
// reads. NaN and the NoData sentinel are not nodes; `full` says a block has
// none of them. Level 0 is read from the raster; levels 1 up are stored, the
// last being one block.
class MaxPyramid {
public:
    explicit MaxPyramid(const raster::RasterView<float>& v) : v_{v} {
        const std::size_t rows = v.geometry().rows(), cols = v.geometry().cols();
        while ((std::size_t{1} << (levels_ - 1)) < std::max(rows, cols))
            ++levels_;
        for (std::size_t level = 1; level < levels_; ++level) {
            const std::size_t nr = blocks(rows, level), nc = blocks(cols, level);
            std::vector<float> m(nr * nc);
            std::vector<std::uint8_t> f(nr * nc);
            for (std::size_t r = 0; r < nr; ++r)
                for (std::size_t c = 0; c < nc; ++c) {
                    auto [best, all] = level == 1 ? nodes(2 * r, 2 * c, 2) : std::pair{kNone, true};
                    if (level > 1) {
                        for (std::size_t a = 2 * r; a < std::min(2 * r + 2, blocks(rows, level - 1)); ++a)
                            for (std::size_t b = 2 * c; b < std::min(2 * c + 2, blocks(cols, level - 1)); ++b) {
                                const std::size_t k = a * blocks(cols, level - 1) + b;
                                best = std::max(best, max_[level - 2][k]);
                                all = all && full_[level - 2][k] != 0;
                            }
                    }
                    m[r * nc + c] = best;
                    f[r * nc + c] = all ? 1 : 0;
                }
            max_.push_back(std::move(m));
            full_.push_back(std::move(f));
        }
    }

    std::size_t levels() const noexcept { return levels_; }

    float block_max(std::size_t level, std::size_t row, std::size_t col) const noexcept {
        if (level == 0)
            return nodes(row, col, 1).first;
        return max_[level - 1][row * blocks(v_.geometry().cols(), level) + col];
    }

    bool block_full(std::size_t level, std::size_t row, std::size_t col) const noexcept {
        return full_[level - 1][row * blocks(v_.geometry().cols(), level) + col] != 0;
    }

    float global_max() const noexcept { return block_max(levels_ - 1, 0, 0); }

private:
    static constexpr float kNone = -std::numeric_limits<float>::infinity();

    static std::size_t blocks(std::size_t n, std::size_t level) noexcept {
        return (n + (std::size_t{1} << level) - 1) >> level;
    }

    // The largest valid node of rows and columns [r0, r0 + span], clipped, and
    // whether every one is valid.
    std::pair<float, bool> nodes(std::size_t r0, std::size_t c0, std::size_t span) const noexcept {
        float best = kNone;
        bool all = true;
        const auto& g = v_.geometry();
        for (std::size_t r = r0; r <= std::min(r0 + span, g.rows() - 1); ++r)
            for (std::size_t c = c0; c <= std::min(c0 + span, g.cols() - 1); ++c) {
                if (v_.is_nodata({r, c}))
                    all = false;
                else
                    best = std::max(best, v_.value_at({r, c}));
            }
        return {best, all};
    }

    raster::RasterView<float> v_;
    std::size_t levels_ = 1;
    std::vector<std::vector<float>> max_;
    std::vector<std::vector<std::uint8_t>> full_;
};

namespace detail {

// Samples d(j) = base + j step for j in [first, end) on one raster.
struct Segment {
    const raster::RasterView<float>* raster;
    const MaxPyramid* pyramid;
    double base, step;
    std::size_t first, end;
    double d(std::size_t j) const noexcept { return base + static_cast<double>(j) * step; }
};

struct Ray {
    Point2 o, u;
    double z_obs, c;
    Point2 at(double d) const noexcept { return {o.x + d * u.x, o.y + d * u.y}; }
    double tangent(double z, double d) const noexcept { return (z - z_obs - c * d * d) / d; }
    // An upper bound on tangent(z, d) for z <= z_max and d in [d_lo, d_hi]
    // (R4). The slack covers the rounding of a bilinear value and of both
    // formulas, so the bound holds for the computed numbers too.
    double bound(float z_max, double d_lo, double d_hi) const noexcept {
        if (z_max == -std::numeric_limits<float>::infinity())
            return -std::numeric_limits<double>::infinity();
        const double a = z_max - z_obs + 1e-6 + 1e-9 * (std::abs(z_max) + std::abs(z_obs));
        return (a >= 0.0 ? a / d_lo : a / d_hi) - c * d_lo;
    }
};

// The number of j >= 0 with base + j step <= limit.
inline std::size_t count_upto(double base, double step, double limit) noexcept {
    if (!(base <= limit))
        return 0;
    auto n = static_cast<std::size_t>(std::floor((limit - base) / step));
    while (base + static_cast<double>(n + 1) * step <= limit)
        ++n;
    while (n > 0 && base + static_cast<double>(n) * step > limit)
        --n;
    return n + 1;
}

// Whole km, held below 255, which means nothing was missing.
inline std::uint8_t km(double d) noexcept {
    return static_cast<std::uint8_t>(std::min(254.0, std::floor(d / 1000.0)));
}

struct Search {
    const ReceiverSet& rc;
    std::array<const raster::RasterView<float>*, 2> raster;
    std::array<const MaxPyramid*, 2> pyramid;
    double c, d_near, d_far;
    std::size_t n;

    std::array<Segment, 2> segments(std::size_t i) const noexcept {
        const double r_i = rc.exclusion_radius[i];
        const auto half = [](const raster::RasterView<float>* v) {
            return 0.5 * std::min(v->geometry().delta_x(), v->geometry().delta_y());
        };
        Segment near{raster[0], pyramid[0], r_i, half(raster[0]), r_i == 0.0 ? 1u : 0u, 0};
        near.end = std::max(near.first, count_upto(r_i, near.step, d_near));
        Segment far{raster[1], pyramid[1], d_near, raster[1] ? half(raster[1]) : 1.0, 1, 0};
        far.end = raster[1] ? std::max<std::size_t>(1, count_upto(d_near, far.step, d_far)) : 1;
        return {near, far};
    }

    Ray ray(std::size_t i, std::size_t k) const noexcept {
        const double azimuth = static_cast<double>(k) * 360.0 / static_cast<double>(n);
        const auto u = frame_step(rc.jacobian.subspan(4 * i).first<4>(), azimuth);
        return {{rc.x[i], rc.y[i]}, {u[0], u[1]}, rc.z_obs[i], c};
    }

    // The reference: every sample, in order, until the first missing one.
    void exhaustive(std::size_t i, std::size_t k, double& best, std::uint8_t& reach) const {
        const Ray r = ray(i, k);
        for (const Segment& s : segments(i))
            for (std::size_t j = s.first; j < s.end; ++j) {
                const double d = s.d(j);
                const auto z = raster::bilinear(*s.raster, r.at(d));
                if (!z) {
                    reach = km(d);
                    return;
                }
                best = std::max(best, r.tangent(*z, d));
            }
    }

    // The last sample in [j, stop) whose cell lies in the level-L block of
    // `cell`. Cell indices are monotone along the ray, so every sample between
    // j and it is in the block too.
    static std::size_t last_in_block(const Segment& s, const Ray& r, std::size_t level, raster::CellIndex cell,
                                     std::size_t j, std::size_t stop) noexcept {
        const auto& g = s.raster->geometry();
        const std::size_t R = cell.row >> level, C = cell.col >> level, w = std::size_t{1} << level;
        const auto leave = [](double v, double lo, double hi, double o) {
            return v > 0.0 ? (hi - o) / v : v < 0.0 ? (lo - o) / v : std::numeric_limits<double>::infinity();
        };
        const double x0 = g.x_min() + static_cast<double>(C * w) * g.delta_x();
        const double x1 = g.x_min() + static_cast<double>(std::min((C + 1) * w, g.cols() - 1)) * g.delta_x();
        const double y1 = g.y_max() - static_cast<double>(R * w) * g.delta_y();
        const double y0 = g.y_max() - static_cast<double>(std::min((R + 1) * w, g.rows() - 1)) * g.delta_y();
        const double d_exit = std::min(leave(r.u.x, x0, x1, r.o.x), leave(r.u.y, y0, y1, r.o.y));
        const double jd = std::floor((d_exit - s.base) / s.step);
        std::size_t last = jd >= static_cast<double>(stop - 1) ? stop - 1
                         : jd <= static_cast<double>(j)       ? j
                                                              : static_cast<std::size_t>(jd);
        for (; last > j; --last) {
            const auto c = g.bilinear_cell_of(r.at(s.d(last)));
            if (c && c->row >> level == R && c->col >> level == C)
                break;
        }
        return last;
    }

    void pruned(std::size_t i, std::size_t k, double& best, std::uint8_t& reach) const {
        const Ray r = ray(i, k);
        const auto segs = segments(i);
        // Coverage first (pin 1): the first missing sample ends the march.
        std::array<std::size_t, 2> stop{segs[0].end, segs[1].end};
        for (std::size_t si = 0; si < 2 && reach == 255; ++si) {
            const Segment& s = segs[si];
            for (std::size_t j = s.first; j < s.end;) {
                const auto cell = s.raster->geometry().bilinear_cell_of(r.at(s.d(j)));
                std::size_t level = s.pyramid->levels() - 1;
                while (cell && level > 0 && !s.pyramid->block_full(level, cell->row >> level, cell->col >> level))
                    --level;
                if (cell && level > 0) {
                    j = last_in_block(s, r, level, *cell, j, s.end) + 1;
                } else if (!raster::bilinear(*s.raster, r.at(s.d(j)))) {
                    reach = km(s.d(j));
                    stop[si] = j;
                    if (si == 0)
                        stop[1] = segs[1].first;
                    break;
                } else {
                    ++j;
                }
            }
        }
        // The heights, up to the last sample before the missing one.
        const double d_last = stop[1] > segs[1].first ? segs[1].d(stop[1] - 1)
                            : stop[0] > segs[0].first ? segs[0].d(stop[0] - 1) : 0.0;
        const float far_max =
            stop[1] > segs[1].first ? segs[1].pyramid->global_max() : -std::numeric_limits<float>::infinity();
        const std::array<float, 2> ahead{std::max(segs[0].pyramid->global_max(), far_max), far_max};
        for (std::size_t si = 0; si < 2; ++si) {
            const Segment& s = segs[si];
            for (std::size_t j = s.first; j < stop[si];) {
                const double d = s.d(j);
                if (r.bound(ahead[si], d, d_last) < best)
                    return;  // pin 2: nothing ahead can beat the best
                const auto cell = s.raster->geometry().bilinear_cell_of(r.at(d));
                bool skipped = false;
                for (std::size_t level = s.pyramid->levels() - 1; level > 0 && !skipped; --level) {
                    const float z = s.pyramid->block_max(level, cell->row >> level, cell->col >> level);
                    if (r.bound(z, d, d) >= best && z - r.z_obs >= 0.0)
                        continue;  // with a >= 0 the bound does not depend on d_hi
                    const std::size_t last = last_in_block(s, r, level, *cell, j, stop[si]);
                    if (r.bound(z, d, s.d(last)) < best) {
                        j = last + 1;
                        skipped = true;
                    }
                }
                if (!skipped) {
                    best = std::max(best, r.tangent(*raster::bilinear(*s.raster, r.at(d)), d));
                    ++j;
                }
            }
        }
    }
};

// Refusals name the field, or the receiver check_receivers refuses.
inline void validate(const ReceiverSet& rc, bool far, const HorizonParams& p) {
    if (p.n_directions == 0)
        throw std::invalid_argument("n_directions must be at least 1");
    if (!(p.curvature >= 0.0) || !std::isfinite(p.curvature))
        throw std::invalid_argument("curvature must be finite and not negative");
    if (!(p.d_near > 0.0) || !std::isfinite(p.d_near))
        throw std::invalid_argument("d_near must be positive and finite");
    if (far && (!(p.d_far >= p.d_near) || !std::isfinite(p.d_far)))
        throw std::invalid_argument("d_far must be finite and at least d_near");
    require_receivers(rc);
}

inline Search prepare(const ReceiverSet& rc, const raster::RasterView<float>& near,
                      const std::optional<raster::RasterView<float>>& far, const HorizonParams& p,
                      const MaxPyramid* near_pyramid, const MaxPyramid* far_pyramid) {
    return Search{rc, {&near, far ? &*far : nullptr}, {near_pyramid, far_pyramid}, p.curvature, p.d_near, p.d_far,
                  p.n_directions};
}

template <bool Pruned>
HorizonTable run(const Search& s, unsigned threads) {
    const std::size_t n = s.rc.size(), N = s.n;
    HorizonTable t{std::vector<std::int16_t>(n * N), std::vector<std::uint8_t>(n * N)};
    parallel_util::for_each_block(n, threads, {}, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i)
            for (std::size_t k = 0; k < N; ++k) {
                double best = -std::numeric_limits<double>::infinity();
                std::uint8_t reach = 255;
                if constexpr (Pruned)
                    s.pruned(i, k, best, reach);
                else
                    s.exhaustive(i, k, best, reach);
                t.horizon[i * N + k] = horizon_centidegrees(best);
                t.reach[i * N + k] = reach;
            }
    });
    return t;
}

}  // namespace detail

// The reference definition, sample by sample, on the calling thread.
inline HorizonTable horizon_exhaustive(const ReceiverSet& rc, const raster::RasterView<float>& near,
                                       const std::optional<raster::RasterView<float>>& far, const HorizonParams& p) {
    detail::validate(rc, far.has_value(), p);
    return detail::run<false>(detail::prepare(rc, near, far, p, nullptr, nullptr), 1);
}

// The same table through the max-pyramids, receivers over `threads` threads
// (0: the hardware's count).
inline HorizonTable horizon_search(const ReceiverSet& rc, const raster::RasterView<float>& near,
                                   const std::optional<raster::RasterView<float>>& far, const HorizonParams& p,
                                   unsigned threads) {
    detail::validate(rc, far.has_value(), p);
    const MaxPyramid near_pyramid{near};
    const std::optional<MaxPyramid> far_pyramid = far ? std::optional<MaxPyramid>{MaxPyramid{*far}} : std::nullopt;
    return detail::run<true>(
        detail::prepare(rc, near, far, p, &near_pyramid, far_pyramid ? &*far_pyramid : nullptr), threads);
}

}  // namespace terrain::radiation
