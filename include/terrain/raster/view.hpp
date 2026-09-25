#pragma once

#include <terrain/raster/geometry.hpp>

#include <cstddef>
#include <optional>
#include <span>

namespace terrain::raster {

// Non-owning, row-major view over a caller's contiguous buffer -- in practice a
// numpy array the bindings keep alive. Same semantics as Raster<T>, no copy:
// whoever builds a view owns the buffer's lifetime and guarantees it holds
// geometry.size() elements.
template <typename T>
class RasterView {
public:
    using value_type = T;

    RasterView(RasterGeometry geometry, const T* data, std::optional<T> nodata) noexcept
        : geometry_{geometry}, data_{data}, nodata_{nodata} {}

    [[nodiscard]] const RasterGeometry& geometry() const noexcept { return geometry_; }

    // Unchecked, like Raster::value_at.
    [[nodiscard]] T value_at(const CellIndex& c) const noexcept {
        return data_[geometry_.linear_index(c)];
    }

    [[nodiscard]] bool is_nodata(const CellIndex& c) const noexcept {
        const T v = value_at(c);
        // v != v is NaN, which counts as NoData with or without a sentinel.
        return v != v || (nodata_.has_value() && v == *nodata_);
    }

    [[nodiscard]] std::span<const T> row(std::size_t i) const noexcept {
        return {data_ + i * geometry_.cols(), geometry_.cols()};
    }

private:
    RasterGeometry geometry_;
    const T* data_;
    std::optional<T> nodata_;
};

} // namespace terrain::raster
