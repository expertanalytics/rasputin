// A zero-copy raster view over a NumPy array, for any pybind11 module built
// against the installed core headers. Installed as
// <rasputin_bindings/numpy_raster.hpp>; it lives under bindings/ so the core
// headers under include/terrain/ stay free of pybind11.
//
// `Tag` is never used in the body. Each module declares its own tag, so its
// view is a distinct C++ type: pybind11 finds a type by its typeid name, and
// one struct bound in two modules would let either accept the other's view.

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <terrain/raster/geometry.hpp>
#include <terrain/raster/view.hpp>

#include <cstddef>
#include <format>
#include <optional>
#include <utility>
#include <variant>

// The variant picks the cell type once, at construction; `array` is the numpy
// object the view points into, held so the buffer lives as long as the view
// does. Destroyed by pybind11 with the GIL held.
template <class Tag>
struct TaggedRasterView {
    pybind11::object array;
    std::variant<terrain::raster::RasterView<float>, terrain::raster::RasterView<double>> view;
};

namespace rasputin_bindings {

// Exact dtype and C order, or nothing: a conversion would be a silent copy the
// view then outlives, so anything else is a TypeError and the caller converts.
template <class Tag, typename T>
[[nodiscard]] std::optional<TaggedRasterView<Tag>> try_view(const pybind11::object& array,
                                                            const terrain::raster::RasterGeometry* g,
                                                            std::optional<double> nodata) {
    using Exact = pybind11::array_t<T, pybind11::array::c_style>;
    if (!pybind11::isinstance<Exact>(array))
        return std::nullopt;
    const auto a = array.cast<Exact>();
    if (a.ndim() != 2)
        throw pybind11::value_error(std::format("raster_view: array must be 2-D, got {} dimensions",
                                                a.ndim()));
    // Rows and columns come from the array's own shape, never from an argument.
    const terrain::raster::RasterGeometry geometry{
        g->x_min(), g->y_max(), g->delta_x(), g->delta_y(),
        static_cast<std::size_t>(a.shape(1)), static_cast<std::size_t>(a.shape(0))};
    // One conversion of the sentinel to T. decode_dem has already refused a
    // sentinel the cell type cannot hold, so it is exact.
    const std::optional<T> sentinel =
        nodata ? std::optional<T>{static_cast<T>(*nodata)} : std::nullopt;
    return TaggedRasterView<Tag>{array, terrain::raster::RasterView<T>{geometry, a.data(), sentinel}};
}

// The body of a module's `raster_view`: a float32 or float64 view, or a TypeError.
template <class Tag>
[[nodiscard]] TaggedRasterView<Tag> raster_view(const pybind11::object& array, double x_min,
                                                double y_max, double delta_x, double delta_y,
                                                std::optional<double> nodata) {
    // Placeholder cols/rows of 1: only the affine half is read from it,
    // and the RasterGeometry constructor validates the deltas here.
    const terrain::raster::RasterGeometry affine{x_min, y_max, delta_x, delta_y, 1, 1};
    if (auto v = try_view<Tag, float>(array, &affine, nodata))
        return std::move(*v);
    if (auto v = try_view<Tag, double>(array, &affine, nodata))
        return std::move(*v);
    throw pybind11::type_error("raster_view: array must be a C-contiguous float32 or float64 "
                               "numpy array; convert it explicitly rather than have it copied");
}

}  // namespace rasputin_bindings
