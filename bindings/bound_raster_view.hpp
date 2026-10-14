// The bound RasterView, shared by the binding files: core.cpp builds it
// (raster_view), radiation.cpp reads it (horizon_search).

#pragma once

#include <pybind11/pybind11.h>

#include <terrain/raster/view.hpp>

#include <variant>

// The variant picks the cell type once, at construction; `array` is the numpy
// object the view points into, held so the buffer lives as long as the view
// does. Destroyed by pybind11 with the GIL held.
struct BoundRasterView {
    pybind11::object array;
    std::variant<terrain::raster::RasterView<float>, terrain::raster::RasterView<double>> view;
};
