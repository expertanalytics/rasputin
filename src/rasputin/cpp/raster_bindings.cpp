#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "cgal_raster_data.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<double>>);

void bind_raster(py::module& m) {
    rasputin::bind_raster<rasputin::RasterData>(m);
}
