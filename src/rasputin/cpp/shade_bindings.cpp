#include <pybind11/pybind11.h>

#include "cgal_shade.h"

void bind_shades(py::module& m) {
    rasputin::bind_shades<rasputin::Mesh>(m);
}
