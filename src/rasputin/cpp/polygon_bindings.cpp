#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "cgal_polygon.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(CGAL::Polygon);
PYBIND11_MAKE_OPAQUE(CGAL::SimplePolygon);
PYBIND11_MAKE_OPAQUE(CGAL::MultiPolygon);

void bind_polygons(py::module& m) {
    rasputin::bind_polygons<rasputin::PolygonUtils<CGAL::Polygon, CGAL::SimplePolygon>>(m);
}

