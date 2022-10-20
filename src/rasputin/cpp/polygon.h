#pragma once
#include <concepts>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace rasputin {
namespace traits {
template<typename Polygon, typename SimplePolygon>
struct PolygonUtils {
    static_assert("PolygonUtils trait not implemented");
};

template<typename P, typename Polygon, typename SimplePolygon>
concept PolygonType = std::same_as<P, Polygon> || std::same_as<P, SimplePolygon>;
} // namespace traits

template<typename Polygon, typename SimplePolygon>
struct PolygonUtils {
    using P = traits::PolygonUtils<Polygon, SimplePolygon>;
    using Polygon_t = Polygon;
    using SimplePolygon_t = SimplePolygon;
    using MultiPolygon_t = std::vector<Polygon>;

    template<typename Point2>
    static bool point_inside_polygon(const Point2 &x, const SimplePolygon &polygon) {
        return P::point_inside_polygon(x, polygon);
    }

    template<typename Point2>
    static bool point_inside_polygon(const Point2 &x, const Polygon &polygon) {
        return P::point_inside_polygon(x, polygon);
    }

    template<typename Point2>
    static bool point_inside_polygon(const Point2 &x, const MultiPolygon_t &polygon) {
        for (const auto& part: polygon)
            if (P::point_inside_polygon(x, part))
                return true;
        return false;
    }

    static std::vector<SimplePolygon> extract_boundaries(const SimplePolygon &polygon) {
        std::vector<SimplePolygon> ret;
        ret.emplace_back(polygon);

        return ret;
    }

    static py::array_t<double> polygon_array(const SimplePolygon& polygon) {
        return P::polygon_array(polygon);
    }

    static std::vector<SimplePolygon> extract_boundaries(const Polygon &polygon) {
        return P::extract_boundaries(polygon);
    }

    static std::vector<SimplePolygon> extract_boundaries(const MultiPolygon_t &polygon) {
        return P::extract_boundaries(polygon);
    }

    static SimplePolygon polygon_from_numpy(py::array_t<double>& buf) {
        return P::polygon_from_numpy(buf);
    }

    static MultiPolygon_t join_multipolygons(const MultiPolygon_t& polygon0, const MultiPolygon_t& polygon1) {
        return P::join_multipolygons(polygon0, polygon1);
    }

    template<typename P0, typename P1>
    requires traits::PolygonType<P0, Polygon_t, SimplePolygon_t> && traits::PolygonType<P1, Polygon_t, SimplePolygon_t>
    static MultiPolygon_t join_polygons(const P0& polygon0, const P1& polygon1) {
        return P::join_polygons(polygon0, polygon1);
    }

    template<typename P0, typename P1>
    requires traits::PolygonType<P0, Polygon_t, SimplePolygon_t> && traits::PolygonType<P1, Polygon_t, SimplePolygon_t>
    static MultiPolygon_t intersect_polygons(const P0& polygon0, const P1& polygon1) {
        return P::intersect_polygons(polygon0, polygon1);
    }

    template<typename P0, typename P1>
    requires traits::PolygonType<P0, Polygon_t, SimplePolygon_t> && traits::PolygonType<P1, Polygon_t, SimplePolygon_t>
    static MultiPolygon_t difference_polygons(const P0& polygon0, const P1& polygon1) {
        return P::difference_polygons(polygon0, polygon1);
    }

    static SimplePolygon exterior(const Polygon& polygon) {
        return P::exterior(polygon);
    }

    static py::list holes(const Polygon& polygon) {
        return P::holes(polygon);
    }

    static std::size_t num_vertices(const SimplePolygon& polygon) {
        return P::num_vertices(polygon);
    }
};

template<typename P>
concept PolygonUtilsType =
    requires {
        typename P::Polygon_t;
        typename P::SimplePolygon_t;
        typename P::MultiPolygon_t;
    }
    && std::same_as<typename P::MultiPolygon_t, std::vector<typename P::Polygon_t>>;

template<PolygonUtilsType P>
void bind_polygon_class(py::module& m) {
    using Polygon = typename P::Polygon_t;
    using SimplePolygon = typename P::SimplePolygon_t;

    py::class_<Polygon, std::unique_ptr<Polygon>>(m, "polygon")
        .def(py::init([] (py::array_t<double>& buf) {
            const SimplePolygon exterior = P::polygon_from_numpy(buf);
            return Polygon(exterior);}))
        .def("holes", &P::holes)
        .def("extract_boundaries", [] (const Polygon& self) { return P::extract_boundaries(self); },
            py::return_value_policy::take_ownership)
        .def("exterior", &P::exterior, py::return_value_policy::take_ownership)
        .def("join", &P::template join_polygons<Polygon, SimplePolygon>)
        .def("join", &P::template join_polygons<Polygon, Polygon>)
        .def("difference", &P::template difference_polygons<Polygon, SimplePolygon>)
        .def("difference", &P::template difference_polygons<Polygon, Polygon>)
        .def("intersection", &P::template intersect_polygons<Polygon, SimplePolygon>)
        .def("intersection", &P::template intersect_polygons<Polygon, Polygon>);
}

template<PolygonUtilsType P>
void bind_simple_polygon_class(py::module& m) {
    using Polygon = typename P::Polygon_t;
    using SimplePolygon = typename P::SimplePolygon_t;

    py::class_<SimplePolygon, std::unique_ptr<SimplePolygon>>(m, "simple_polygon")
        .def(py::init(&P::polygon_from_numpy))
        .def("num_vertices", &P::num_vertices)
        .def("array",
            [] (const SimplePolygon& self) { return P::polygon_array(self); },
            py::return_value_policy::move)
        .def("join", &P::template join_polygons<SimplePolygon, SimplePolygon>)
        .def("join", &P::template join_polygons<SimplePolygon, Polygon>)
        .def("difference", &P::template difference_polygons<SimplePolygon, SimplePolygon>)
        .def("difference", &P::template difference_polygons<SimplePolygon, Polygon>)
        .def("intersection", &P::template intersect_polygons<SimplePolygon, SimplePolygon>)
        .def("intersection", &P::template intersect_polygons<SimplePolygon, Polygon>);
}

template<PolygonUtilsType P>
void bind_multi_polygon_class(py::module& m) {
    using Polygon = typename P::Polygon_t;
    using SimplePolygon = typename P::SimplePolygon_t;
    using MultiPolygon = typename P::MultiPolygon_t;

    py::class_<MultiPolygon, std::unique_ptr<MultiPolygon>>(m, "multi_polygon")
        .def(py::init(
            [] (const Polygon& polygon) {
                MultiPolygon self;
                self.push_back(Polygon(polygon));
                return self;
            }))
        .def(py::init(
            [] (const SimplePolygon& polygon) {
                MultiPolygon self;
                self.push_back(Polygon(polygon));
                return self;
            }))
        .def("num_parts", &MultiPolygon::size)
        .def("parts",
            [] (const MultiPolygon& self) {
                py::list result;
                for (auto p = self.begin(); p != self.end(); ++p)
                    result.append(*p);
                return result;
            })
        .def("extract_boundaries", [] (const MultiPolygon& self) { return P::extract_boundaries(self); },
            py::return_value_policy::take_ownership)
        .def("join",
            [] (const MultiPolygon a, SimplePolygon b) {
                return P::join_multipolygons(a, MultiPolygon({static_cast<Polygon>(b)}));
            })
        .def("join",
            [] (const MultiPolygon a, Polygon b) {
                return P::join_multipolygons(a, MultiPolygon({b}));
            })
        .def("join",
            [] (const MultiPolygon a, MultiPolygon b) {
                return P::join_multipolygons(a, b);
            })
        .def("__getitem__", [] (const MultiPolygon& self, int idx) {
                return self.at(idx);
            }, py::return_value_policy::reference_internal);
}

template<PolygonUtilsType P>
void bind_polygons(py::module& m) {
    bind_polygon_class<P>(m);
    bind_simple_polygon_class<P>(m);
    bind_multi_polygon_class<P>(m);
}
} // namespace rasputin

