#pragma once

#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>

namespace terrain::raster {

// Bilinear interpolation over the cell containing p.
//
// Returns nullopt rather than a plausible number when p is outside the
// domain, when the raster is too small to hold a 2x2 neighbourhood, or when
// any of the four corners is NoData. The legacy returned a value in all three
// cases -- reading out of bounds in the second.
template <RasterSource R>
[[nodiscard]] std::optional<double> bilinear(const R& raster, const Point2& p) noexcept {
    const RasterGeometry& g = raster.geometry();

    const auto cell = g.bilinear_cell_of(p);
    if (!cell)
        return std::nullopt;

    const CellIndex c00{cell->row, cell->col};
    const CellIndex c01{cell->row, cell->col + 1};
    const CellIndex c10{cell->row + 1, cell->col};
    const CellIndex c11{cell->row + 1, cell->col + 1};

    if (raster.is_nodata(c00) || raster.is_nodata(c01)
        || raster.is_nodata(c10) || raster.is_nodata(c11))
        return std::nullopt;

    const Point2 upper_left = g.node(c00);
    const double tx = (p.x - upper_left.x) / g.delta_x();
    const double ty = (upper_left.y - p.y) / g.delta_y();

    const double z00 = static_cast<double>(raster.value_at(c00));
    const double z01 = static_cast<double>(raster.value_at(c01));
    const double z10 = static_cast<double>(raster.value_at(c10));
    const double z11 = static_cast<double>(raster.value_at(c11));

    return z00 * (1.0 - tx) * (1.0 - ty)
         + z01 * tx * (1.0 - ty)
         + z10 * (1.0 - tx) * ty
         + z11 * tx * ty;
}

// bilinear over every point: z[i] is the sample and valid[i] says whether
// there is one. Where valid[i] is false, z[i] is 0.0 -- never NaN, so a
// caller that ignores the flag still cannot leak a NaN into a file.
// Spans of different lengths are a caller bug and throw before any write.
template <RasterSource R>
void bilinear_batch(const R& raster, std::span<const Point2> points,
                    std::span<double> z, std::span<bool> valid) {
    if (z.size() != points.size() || valid.size() != points.size())
        throw std::invalid_argument("bilinear_batch: z and valid must match points in length");
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto v = bilinear(raster, points[i]);
        valid[i] = v.has_value();
        z[i] = v.value_or(0.0);
    }
}

} // namespace terrain::raster
