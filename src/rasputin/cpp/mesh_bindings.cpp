#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cgal_mesh.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(rasputin::point3_vector);
PYBIND11_MAKE_OPAQUE(rasputin::point2_vector);
PYBIND11_MAKE_OPAQUE(rasputin::face_vector);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<double>>);

void bind_meshes(py::module& m) {
    rasputin::bind_mesh<rasputin::Mesh, rasputin::RasterData<float>, rasputin::RasterData<double>>(m);
}
