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

    py::class_<Point2>(m, "Point2")
        .def(py::init([](double x, double y) { return Point2{x, y}; }),
             py::arg("x") = 0.0, py::arg("y") = 0.0)
        .def_readwrite("x", &Point2::x)
        .def_readwrite("y", &Point2::y)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def(-py::self)
        // Mirrors the std::formatter specialization in point.hpp, so the C++
        // and Python renderings of a point cannot drift apart.
        .def("__repr__", [](const Point2& p) { return std::format("{}", p); })
        .def("__hash__", [](const Point2& p) { return std::hash<Point2>{}(p); });

    py::class_<Point3>(m, "Point3")
        .def(py::init([](double x, double y, double z) { return Point3{x, y, z}; }),
             py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0)
        .def_readwrite("x", &Point3::x)
        .def_readwrite("y", &Point3::y)
        .def_readwrite("z", &Point3::z)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def(-py::self)
        .def("__repr__", [](const Point3& p) { return std::format("{}", p); })
        .def("__hash__", [](const Point3& p) { return std::hash<Point3>{}(p); });

    // Overloaded rather than method-bound: dot and cross are symmetric binary
    // operations, and mixing dimensions is a TypeError rather than a coercion.
    m.def("dot", py::overload_cast<const Point2&, const Point2&>(&terrain::dot),
          py::arg("a"), py::arg("b"));
    m.def("dot", py::overload_cast<const Point3&, const Point3&>(&terrain::dot),
          py::arg("a"), py::arg("b"));

    // 2D cross is the scalar z-component; 3D cross is a vector.
    m.def("cross", py::overload_cast<const Point2&, const Point2&>(&terrain::cross),
          py::arg("a"), py::arg("b"));
    m.def("cross", py::overload_cast<const Point3&, const Point3&>(&terrain::cross),
          py::arg("a"), py::arg("b"));
}
