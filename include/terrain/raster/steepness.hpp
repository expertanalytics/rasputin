#pragma once

// Horn's slope of every DEM node as a class: half degrees, rounded up, 0 to
// 180, and kNoDataClass for a NoData node
// (docs/increments/34-slope-tolerance.md, sections 3 and 4.2).
//
// A missing neighbour (outside the grid, or NoData) is filled so that a plane
// reads its own slope: an edge neighbour by its reflection through the node,
// 2 z(node) - z(opposite), when the opposite one is there, else by z(node); a
// corner neighbour by its reflection when the opposite corner is there, else
// by z(row neighbour) + z(column neighbour) - z(node), from the edge
// neighbours as filled.
//
// The class is the count of half degrees whose tan^2 is at most gx^2 + gy^2,
// plus one (0 for a flat node), from a table, with no atan per node. Each
// tan^2 is lowered by a relative 1e-12 (about 1e-11 degrees), so a slope at a
// class boundary, within rounding, takes the higher class (G5). Each byte is
// a function of the node's 3 by 3 neighbourhood only, so the result does not
// depend on `threads`.

#include <terrain/parallel_util/chunks.hpp>
#include <terrain/raster/raster.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

namespace terrain::raster {

inline constexpr std::size_t kSteepnessClasses = 181;  // half degrees, 0 to 90
inline constexpr std::uint8_t kNoDataClass = 255;     // a NoData node (section 4.5)

// tan^2 of 0.5 to 89.5 degrees, each lowered by a relative 1e-12; ascending.
[[nodiscard]] inline const std::array<double, kSteepnessClasses - 2>& steepness_bounds() {
    static const auto table = [] {
        std::array<double, kSteepnessClasses - 2> t{};
        for (std::size_t k = 0; k < t.size(); ++k) {
            const double a = std::tan(static_cast<double>(k + 1) * std::numbers::pi / 360.0);
            t[k] = a * a * (1.0 - 1e-12);
        }
        return t;
    }();
    return table;
}

template <RasterSource R>
[[nodiscard]] std::vector<std::uint8_t> steepness(const R& dem, unsigned threads) {
    using T = typename R::value_type;
    const RasterGeometry& g = dem.geometry();
    const auto rows = static_cast<std::ptrdiff_t>(g.rows()), cols = static_cast<std::ptrdiff_t>(g.cols());
    const double dx8 = 8.0 * g.delta_x(), dy8 = 8.0 * g.delta_y();
    const auto& bounds = steepness_bounds();
    const std::optional<T>& nd = dem.nodata();
    std::vector<std::uint8_t> out(g.size());
    parallel_util::for_each_block(g.rows(), threads, parallel_util::BlockSchedule{}, [&](std::size_t r0, std::size_t r1) {
        for (auto r = static_cast<std::ptrdiff_t>(r0); r < static_cast<std::ptrdiff_t>(r1); ++r) {
            const std::array<std::span<const T>, 3> line{r > 0 ? dem.row(static_cast<std::size_t>(r - 1)) : std::span<const T>{},
                                                         dem.row(static_cast<std::size_t>(r)),
                                                         r + 1 < rows ? dem.row(static_cast<std::size_t>(r + 1)) : std::span<const T>{}};
            for (std::ptrdiff_t c = 0; c < cols; ++c) {
                // The height at (r + dr, c + dc), NaN where it is missing.
                const auto at = [&](std::ptrdiff_t dr, std::ptrdiff_t dc) {
                    const auto& s = line[static_cast<std::size_t>(dr + 1)];
                    const std::ptrdiff_t j = c + dc;
                    if (s.empty() || j < 0 || j >= cols)
                        return std::numeric_limits<double>::quiet_NaN();
                    const T v = s[static_cast<std::size_t>(j)];
                    return v != v || (nd && v == *nd) ? std::numeric_limits<double>::quiet_NaN() : static_cast<double>(v);
                };
                const double z = at(0, 0);
                std::uint8_t& to = out[static_cast<std::size_t>(r * cols + c)];
                if (z != z) {
                    to = kNoDataClass;
                    continue;
                }
                const auto edge = [&](std::ptrdiff_t dr, std::ptrdiff_t dc) {
                    const double v = at(dr, dc), o = v == v ? v : at(-dr, -dc);
                    return v == v ? v : o == o ? 2.0 * z - o : z;
                };
                const auto corner = [&](std::ptrdiff_t dr, std::ptrdiff_t dc) {
                    const double v = at(dr, dc), o = v == v ? v : at(-dr, -dc);
                    return v == v ? v : o == o ? 2.0 * z - o : edge(dr, 0) + edge(0, dc) - z;
                };
                const double a = corner(-1, -1), b = edge(-1, 0), cc = corner(-1, 1), d = edge(0, -1);
                const double f = edge(0, 1), gg = corner(1, -1), h = edge(1, 0), i = corner(1, 1);
                const double gx = ((cc + 2.0 * f + i) - (a + 2.0 * d + gg)) / dx8;
                const double gy = ((gg + 2.0 * h + i) - (a + 2.0 * b + cc)) / dy8;
                const double s2 = gx * gx + gy * gy;
                to = s2 == 0.0 ? 0
                               : static_cast<std::uint8_t>(1 + std::upper_bound(bounds.begin(), bounds.end(), s2)
                                                               - bounds.begin());
            }
        }
    });
    return out;
}

}  // namespace terrain::raster
