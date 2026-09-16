#pragma once

#include <terrain/raster/geometry.hpp>

#include <concepts>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace terrain::raster {

// Everything that samples a raster goes through this contract, so an owning
// grid and a zero-copy view over a numpy buffer stay interchangeable.
template <typename R>
concept RasterSource = requires(const R& r, CellIndex c) {
    typename R::value_type;
    // same_as, not convertible_to: an implementation returning by value
    // satisfies convertible_to, and sample.hpp binds the result to a
    // reference. Lifetime extension saves that today, but the concept should
    // state the requirement rather than rely on the call site's shape.
    { r.geometry() } -> std::same_as<const RasterGeometry&>;
    { r.value_at(c) } -> std::convertible_to<double>;
    { r.is_nodata(c) } -> std::same_as<bool>;
};

// Owning raster. The zero-copy RasterView over a Python buffer lands with the
// reader; this one exists so the C++ core is testable with no I/O at all.
template <typename T>
class Raster {
public:
    using value_type = T;

    Raster(RasterGeometry geometry, std::vector<T> data,
           std::optional<T> nodata = std::nullopt)
        : geometry_{geometry}, data_{std::move(data)}, nodata_{nodata} {
        if (data_.size() != geometry_.size())
            throw std::invalid_argument("Raster: data size does not match geometry");
    }

    [[nodiscard]] const RasterGeometry& geometry() const noexcept { return geometry_; }

    // Deliberately not named `at`: in C++ that means bounds-checked and
    // throwing, and this is neither. Unchecked is the right default on a
    // sampling hot path, but the name has to say so -- especially in a module
    // whose reason for existing is an out-of-bounds read.
    [[nodiscard]] T value_at(const CellIndex& c) const noexcept {
        return data_[geometry_.linear_index(c)];
    }

    // NoData is new: the legacy had no sentinel handling anywhere, so a void
    // in the DEM silently interpolated as terrain. It is a constructor
    // argument rather than a setter so a Raster is never half-configured.
    [[nodiscard]] const std::optional<T>& nodata() const noexcept { return nodata_; }

    [[nodiscard]] bool is_nodata(const CellIndex& c) const noexcept {
        const T v = value_at(c);
        // v != v catches NaN, which must never propagate silently into a mesh.
        return v != v || (nodata_.has_value() && v == *nodata_);
    }

private:
    RasterGeometry geometry_;
    std::vector<T> data_;
    std::optional<T> nodata_;
};

} // namespace terrain::raster
