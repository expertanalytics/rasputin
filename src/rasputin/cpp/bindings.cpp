#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "cgal_raster_data.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<double>>);

void bind_vectors(py::module&);
void bind_mesh(py::module&);
void bind_polygons(py::module&);
void bind_solars(py::module&);
void bind_shades(py::module&);
void bind_slopes(py::module&);

PYBIND11_MODULE(triangulate_dem, m) {
    bind_vectors(m);
    rasputin::bind_raster<rasputin::RasterData>(m);
    bind_mesh(m);
    bind_polygons(m);
    bind_slopes(m);
    bind_shades(m);
    bind_solars(m);
}
