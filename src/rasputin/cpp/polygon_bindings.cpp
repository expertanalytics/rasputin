#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Boolean_set_operations_2.h>

#include "polygon.h"

namespace py = pybind11;

template<typename P0, typename P1>
CGAL::MultiPolygon difference_polygons(const P0& polygon0, const P1& polygon1) {
    CGAL::MultiPolygon difference_result;
    CGAL::difference(polygon0, polygon1, std::back_inserter(difference_result));

    return difference_result;
}

template<typename P0, typename P1>
CGAL::MultiPolygon intersect_polygons(const P0& polygon0, const P1& polygon1) {
    CGAL::MultiPolygon intersection_result;
    CGAL::intersection(polygon0, polygon1, std::back_inserter(intersection_result));

    return intersection_result;
}

template<typename P0, typename P1>
CGAL::MultiPolygon join_polygons(const P0& polygon0, const P1& polygon1) {
    CGAL::Polygon joined;
    CGAL::MultiPolygon join_result;
    if (CGAL::join(polygon0, polygon1, joined))
        join_result.push_back(std::move(joined));
    else {
        join_result.push_back(CGAL::Polygon(polygon0));
        join_result.push_back(CGAL::Polygon(polygon1));
    }
    return join_result;
}

CGAL::MultiPolygon join_multipolygons(const CGAL::MultiPolygon& polygon0, const CGAL::MultiPolygon& polygon1) {
    CGAL::MultiPolygon join_result;
    CGAL::join(polygon0.begin(), polygon0.end(), polygon1.begin(), polygon1.end(), std::back_inserter(join_result));

    return join_result;
}

CGAL::SimplePolygon polygon_from_numpy(py::array_t<double>& buf) {
        auto info = buf.request();
        if (info.ndim < 1 or info.ndim > 2)
            throw py::type_error("Can only convert one- and two-dimensional arrays.");

        // Make sure total size can be preserved
        auto size = info.shape[0];
        if (info.ndim == 1 and size % 2)
            throw py::type_error("Size of one-dimensional array must be divisible by 2.");

        if (info.ndim == 2 and info.shape[1] != 2)
            throw py::type_error("Second dimension does not have size equal to 2.");

        if (info.ndim == 2)
            size *= 2;

        CGAL::SimplePolygon exterior;

        // Copy values
        double* p = static_cast<double*>(info.ptr);
        double* end = p + size;
        for (; p < end; p += 2)
            exterior.push_back(CGAL::Point2(*p, *(p+1)));
        return exterior;
}


void bind_polygons(py::module& m) {
    py::class_<CGAL::SimplePolygon, std::unique_ptr<CGAL::SimplePolygon>>(m, "simple_polygon")
        .def(py::init(&polygon_from_numpy))
        .def("num_vertices", &CGAL::SimplePolygon::size)
        .def("array", [] (const CGAL::SimplePolygon& self) {
                py::array_t<double> result(self.size() * 2);
                result.resize(py::array::ShapeContainer({static_cast<long int>(self.size()), 2}), true);
                auto info = result.request();
                double* data = static_cast<double*>(info.ptr);
                for (auto v = self.vertices_begin(); v != self.vertices_end(); ++v) {
                    data[0] = v->x();
                    data[1] = v->y();
                    data += 2;
                }
                return result;
                }, py::return_value_policy::move)
        .def("join", &join_polygons<CGAL::SimplePolygon, CGAL::SimplePolygon>)
        .def("join", &join_polygons<CGAL::SimplePolygon, CGAL::Polygon>)
        .def("difference", &difference_polygons<CGAL::SimplePolygon, CGAL::SimplePolygon>)
        .def("difference", &difference_polygons<CGAL::SimplePolygon, CGAL::Polygon>)
        .def("intersection", &intersect_polygons<CGAL::SimplePolygon, CGAL::SimplePolygon>)
        .def("intersection", &intersect_polygons<CGAL::SimplePolygon, CGAL::Polygon>);

    py::class_<CGAL::Polygon, std::unique_ptr<CGAL::Polygon>>(m, "polygon")
        .def(py::init([] (py::array_t<double>& buf) {
            const CGAL::SimplePolygon exterior = polygon_from_numpy(buf);
            return CGAL::Polygon(exterior);}))
        .def("holes", [] (const CGAL::Polygon& self) {
                    py::list result;
                    for (auto h = self.holes_begin(); h != self.holes_end(); ++h)
                        result.append(*h);
                    return result;
                    })
        .def("extract_boundaries", [] (const CGAL::Polygon& self) {return CGAL::extract_boundaries(self);},
                py::return_value_policy::take_ownership)
        .def("exterior", [] (const CGAL::Polygon& self) {return self.outer_boundary();})
        .def("join", &join_polygons<CGAL::Polygon, CGAL::SimplePolygon>)
        .def("join", &join_polygons<CGAL::Polygon, CGAL::Polygon>)
        .def("difference", &difference_polygons<CGAL::Polygon, CGAL::SimplePolygon>)
        .def("difference", &difference_polygons<CGAL::Polygon, CGAL::Polygon>)
        .def("intersection", &intersect_polygons<CGAL::Polygon, CGAL::SimplePolygon>)
        .def("intersection", &intersect_polygons<CGAL::Polygon, CGAL::Polygon>);

    py::class_<CGAL::MultiPolygon, std::unique_ptr<CGAL::MultiPolygon>>(m, "multi_polygon")
        .def(py::init(
            [] (const CGAL::Polygon& polygon) {
                CGAL::MultiPolygon self;
                self.push_back(CGAL::Polygon(polygon));
                return self;
            }))
        .def(py::init(
            [] (const CGAL::SimplePolygon& polygon) {
                CGAL::MultiPolygon self;
                self.push_back(CGAL::Polygon(polygon));
                return self;
            }))
        .def("num_parts", &CGAL::MultiPolygon::size)
        .def("parts",
            [] (const CGAL::MultiPolygon& self) {
                py::list result;
                for (auto p = self.begin(); p != self.end(); ++p)
                    result.append(*p);
                return result;
            })
        .def("extract_boundaries", [] (const CGAL::MultiPolygon& self) {return CGAL::extract_boundaries(self);},
                py::return_value_policy::take_ownership)
        .def("join",
            [] (const CGAL::MultiPolygon a, CGAL::SimplePolygon b) {
                return join_multipolygons(a, CGAL::MultiPolygon({static_cast<CGAL::Polygon>(b)}));
            })
        .def("join",
            [] (const CGAL::MultiPolygon a, CGAL::Polygon b) {
                return join_multipolygons(a, CGAL::MultiPolygon({b}));
            })
        .def("join",
            [] (const CGAL::MultiPolygon a, CGAL::MultiPolygon b) {
                return join_multipolygons(a, b);
            })
        .def("__getitem__", [] (const CGAL::MultiPolygon& self, int idx) {
                    return self.at(idx);
                    }, py::return_value_policy::reference_internal);
}

