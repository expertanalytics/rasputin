#include <pybind11/pybind11.h>

namespace py = pybind11;

void bind_vectors(py::module&);
void bind_raster(py::module& m);
void bind_meshes(py::module& m);
void bind_polygons(py::module&);
void bind_shades(py::module&);
void bind_solars(py::module&);
void bind_slopes(py::module&);

PYBIND11_MODULE(triangulate_dem, m) {
    bind_vectors(m);
    bind_raster(m);
    bind_meshes(m);
    bind_polygons(m);
    bind_slopes(m);
    bind_shades(m);
    bind_solars(m);
}
