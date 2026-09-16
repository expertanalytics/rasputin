#include <pybind11/operators.h>
#include <pybind11/pybind11.h>

#include <terrain/core/point.hpp>

#include <format>
#include <functional>

namespace py = pybind11;

using terrain::Point2;
using terrain::Point3;

PYBIND11_MODULE(_core, m) {
    m.doc() = "C++20 geometry primitives for the rasputin terrain engine.";

    py::class_<Point2>(m, "Point2", R"doc(
A 2D point with double-precision coordinates.

Immutable: the coordinates are read-only because the type is hashable, and a
mutable hashable value silently corrupts any set or dict holding it.

Equality is exact, component-wise floating-point comparison, so NaN coordinates
never compare equal and +0.0 equals -0.0. Hashing is consistent with that.
Coordinate de-duplication that must tolerate floating-point drift should snap
-round upstream rather than rely on this.
)doc")
        .def(py::init([](double x, double y) { return Point2{x, y}; }),
             py::arg("x") = 0.0, py::arg("y") = 0.0,
             "Construct from coordinates, defaulting to the origin.")
        .def_readonly("x", &Point2::x, "The x coordinate (read-only).")
        .def_readonly("y", &Point2::y, "The y coordinate (read-only).")
        .def(py::self == py::self, "Exact component-wise equality.")
        .def(py::self != py::self, "Exact component-wise inequality.")
        .def(py::self + py::self, "Component-wise sum.")
        .def(py::self - py::self, "Component-wise difference.")
        .def(py::self * double(), "Scale by a scalar.")
        .def(double() * py::self, "Scale by a scalar.")
        .def(py::self / double(), "Divide by a scalar. Division by zero yields inf or nan.")
        .def(-py::self, "Component-wise negation.")
        // Mirrors the std::formatter specialization in point.hpp, so the C++
        // and Python renderings of a point cannot drift apart.
        .def("__repr__", [](const Point2& p) { return std::format("{}", p); },
             "Render as Point2(x, y), matching the C++ std::formatter output.")
        .def("__hash__", [](const Point2& p) { return std::hash<Point2>{}(p); },
             "Hash of the coordinate pair, consistent with equality.");

    py::class_<Point3>(m, "Point3", R"doc(
A 3D point with double-precision coordinates.

Immutable and hashable on the same terms as Point2: exact component-wise
equality, read-only coordinates, and a hash consistent with equality.
)doc")
        .def(py::init([](double x, double y, double z) { return Point3{x, y, z}; }),
             py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0,
             "Construct from coordinates, defaulting to the origin.")
        .def_readonly("x", &Point3::x, "The x coordinate (read-only).")
        .def_readonly("y", &Point3::y, "The y coordinate (read-only).")
        .def_readonly("z", &Point3::z, "The z coordinate (read-only).")
        .def(py::self == py::self, "Exact component-wise equality.")
        .def(py::self != py::self, "Exact component-wise inequality.")
        .def(py::self + py::self, "Component-wise sum.")
        .def(py::self - py::self, "Component-wise difference.")
        .def(py::self * double(), "Scale by a scalar.")
        .def(double() * py::self, "Scale by a scalar.")
        .def(py::self / double(), "Divide by a scalar. Division by zero yields inf or nan.")
        .def(-py::self, "Component-wise negation.")
        .def("__repr__", [](const Point3& p) { return std::format("{}", p); },
             "Render as Point3(x, y, z), matching the C++ std::formatter output.")
        .def("__hash__", [](const Point3& p) { return std::hash<Point3>{}(p); },
             "Hash of the coordinate triple, consistent with equality.");

    // Overloaded rather than method-bound: dot and cross are symmetric binary
    // operations, and mixing dimensions is a TypeError rather than a coercion.
    m.def("dot", py::overload_cast<const Point2&, const Point2&>(&terrain::dot),
          py::arg("a"), py::arg("b"),
          R"doc(
Dot product of two points of the same dimension.

Accepts two Point2 or two Point3. Mixing dimensions raises TypeError rather
than coercing. Runs in constant time.
)doc");
    m.def("dot", py::overload_cast<const Point3&, const Point3&>(&terrain::dot),
          py::arg("a"), py::arg("b"), "Dot product of two Point3.");

    m.def("cross", py::overload_cast<const Point2&, const Point2&>(&terrain::cross),
          py::arg("a"), py::arg("b"),
          R"doc(
Cross product, whose result type depends on the dimension.

For two Point2 this returns a float: the scalar z-component of the 3D cross
product, i.e. a.x*b.y - a.y*b.x. Its sign gives the orientation of the turn
from a to b. For two Point3 it returns a Point3 orthogonal to both.

Mixing dimensions raises TypeError. Runs in constant time.
)doc");
    m.def("cross", py::overload_cast<const Point3&, const Point3&>(&terrain::cross),
          py::arg("a"), py::arg("b"),
          "Cross product of two Point3, returning a Point3 orthogonal to both.");
}
