#include <concepts>
#include <memory>

#include <boost/geometry/geometry.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "triangulate_dem.h"

namespace py = pybind11;

namespace bg = boost::geometry;
using Point = bg::model::point<double, 3, bg::cs::cartesian>;
namespace boost::geometry::traits {
    using shared_type = std::shared_ptr<const Point>;
    BOOST_GEOMETRY_DETAIL_SPECIALIZE_POINT_TRAITS(shared_type, 3, double, cs::cartesian)
    template<>
    struct access<shared_type, 0> {
        static inline double get(shared_type const& p) { return p->get<0>(); }
    };
    template<>
    struct access<shared_type, 1> {
        static inline double get(shared_type const& p) { return p->get<1>(); }
    };
    template<>
    struct access<shared_type, 2> {
        static inline double get(shared_type const& p) { return p->get<2>(); }
    };
}

namespace rasputin {
template<typename T, std::size_t n=3>
requires std::floating_point<T>
py::buffer_info vecarray_buffer(typename RasterData<Point, T>::MultiPoint &v) {
    return py::buffer_info(
        &v[0],                                        /* Pointer to buffer */
        sizeof(T),                                    /* Size of one scalar */
        py::format_descriptor<T>::format(),           /* Python struct-style format descriptor */
        bg::dimension<Point>(),                       /* Number of dimensions */
        std::vector<std::size_t> { v.size(), n },     /* Buffer dimensions */
        { sizeof(T) * n, sizeof(T) }                  /* Strides (in bytes) for each index */
    );
}

template<typename T, std::size_t n>
requires std::floating_point<T>
std::vector<std::array<T, n>>* vecarray_from_numpy(py::array_t<T> buf) {
    auto info = buf.request();
    if (info.ndim < 1 or info.ndim > 2)
        throw py::type_error("Can only convert one- and two-dimensional arrays.");

    // Make sure total size can be preserved
    auto size = info.shape[0];
    if (info.ndim == 1 and size % n)
        throw py::type_error("Size of one-dimensional array must be divisible by!" + std::to_string(n) + ".");

    if (info.ndim == 2 and info.shape[1] != n)
        throw py::type_error("Second dimension does not have size equal to " + std::to_string(n) + ".");

    if (info.ndim == 2)
        size *= n;

    auto vec = std::unique_ptr<std::vector<std::array<T, n>>> (new std::vector<std::array<T, n>> ());
    vec->reserve(size / n);

    // Copy values
    double* p = static_cast<double*>(info.ptr);
    double* end = p + size;
    for (; p < end; p += n) {
        auto a = std::array<T, n> {};
        for (auto k = 0; k < n; k++)
            a[k] = *(p + k);
        vec->push_back(a);
    }

    return vec.release();
}

template<typename T, std::size_t n>
requires std::floating_point<T>
py::buffer_info vector_buffer(std::vector<T> &v) {
    return py::buffer_info(&v[0], sizeof(T), py::format_descriptor<T>::format(), 1, { v.size() }, { sizeof(T) });
}

template<typename FT, typename Point>
void bind_raster_data(py::module& m, const std::string& pyname) {
    using RasterData = rasputin::RasterData<Point, FT>;
    py::class_<RasterData, std::unique_ptr<RasterData>>(m, pyname.c_str(), py::buffer_protocol())
    .def(py::init([] (py::array_t<FT>& data_array, double x_min, double y_max, double delta_x, double delta_y) {
            auto buffer = data_array.request();
            int m = buffer.shape[0], n = buffer.shape[1];

            return RasterData(x_min, y_max, delta_x, delta_y, n, m, static_cast<FT*>(buffer.ptr) );
        }), py::return_value_policy::take_ownership,  py::keep_alive<1, 2>(),
            py::arg("data_array").noconvert(), py::arg("x_min"), py::arg("y_max"), py::arg("delta_x"), py::arg("delta_y"))
    .def_buffer([] (RasterData& self) {
            return py::buffer_info(
                (void*)(&(self.get_points()[0])),
                sizeof(FT),
                py::format_descriptor<FT>::format(),
                bg::dimension<Point>(),
                std::vector<std::size_t> { self.get_num_y(), self.get_num_y() },
                { sizeof(FT) * self.get_num_x(), sizeof(FT) }
            );
            })
    .def_property_readonly("x_max", &RasterData::get_xmax)
    .def_property_readonly("y_min", &RasterData::get_ymin)
    .def("__getitem__", [](RasterData& self, std::pair<int, int> idx) {
                auto [i, j] = idx;
                return self.get_points()[self.get_num_x() * i + j]; })
    .def("get_indices", &RasterData::get_indices)
    .def("exterior", &RasterData::exterior)
    .def("contains", &RasterData::contains)
    .def("get_interpolated_value_at_point", &RasterData::get_interpolated_value_at_point);
}

template<typename T, typename Point>
void bind_raster_list(py::module &m, const std::string& pyname) {
    using RasterData = rasputin::RasterData<Point, T>;
    py::class_<std::vector<RasterData>, std::unique_ptr<std::vector<RasterData>>> (m, pyname.c_str())
        .def(py::init( [] () {std::vector<RasterData> self; return self;}))
        .def("add_raster",
            [] (std::vector<RasterData>& self, RasterData raster_data) {
                self.push_back(raster_data);
            }, py::keep_alive<1,2>())
        .def("__getitem__",
            [] (std::vector<RasterData>& self, int index) {
                return self.at(index);
            }, py::return_value_policy::reference_internal);
}
}

PYBIND11_MAKE_OPAQUE(rasputin::RasterData<Point, double>::MultiPoint)
PYBIND11_MODULE(triangulate_dem, m) {
    py::bind_vector<rasputin::RasterData<Point, double>::MultiPoint>(m, "points", py::buffer_protocol())
        .def_buffer(&rasputin::vecarray_buffer<double, bg::dimension<Point>::value>);
    py::bind_vector<std::vector<int>>(m, "int_vector");
    py::bind_vector<std::vector<std::vector<int>>>(m, "shadow_vector");
    rasputin::bind_raster_data<double, Point>(m, "raster_data");
    rasputin::bind_raster_list<double, Point>(m, "raster_list");
};
