#pragma once

#include <concepts>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "types.h"

namespace py = pybind11;

namespace rasputin {
namespace traits {
template<
    typename FT,
    typename PointList,
    typename SimplePolygon,
    typename MultiPolygon
>
requires std::same_as<FT, double> || std::same_as<FT, float>
class RasterBase {
    public:
    double x_min;
    double delta_x;
    std::size_t num_points_x;

    double y_max;
    double delta_y;
    std::size_t num_points_y;

    FT* data;

    RasterBase(
        double x_min,
        double y_max,
        double delta_x,
        double delta_y,
        std::size_t num_points_x,
        std::size_t num_points_y,
        FT* data
    ) :
        x_min(x_min),
        delta_x(delta_x),
        num_points_x(num_points_x),
        y_max(y_max),
        delta_y(delta_y),
        num_points_y(num_points_y),
        data(data) {
    }

    double get_x_max() const {
        return x_min + (num_points_x - 1)*delta_x;
    }
    double get_y_min() const {
        return y_max - (num_points_y - 1)*delta_y;
    }

    // For every point inside the raster rectangle we identify indices (i, j) of the upper-left vertex of the cell containing the point
    std::pair<int, int> get_indices(double x, double y) const {
        int i = std::min<int>(std::max<int>(static_cast<int>((y_max-y) /delta_y), 0), num_points_y - 1);
        int j = std::min<int>(std::max<int>(static_cast<int>((x-x_min) /delta_x), 0), num_points_x - 1);

        return std::make_pair(i,j);
    }

    // Interpolate data using using a bilinear interpolation rule on each cell
    FT get_interpolated_value_at_point(double x, double y) const {
        // Determine indices of the cell containing (x, y)
        auto [i, j] = get_indices(x, y);

        // Determine the cell corners
        //     (x0, y0) -- upper left
        //     (x1, y1) -- lower right
        double x_0 = x_min + j * delta_x,
               y_0 = y_max - i * delta_y,
               x_1 = x_min + (j+1) * delta_x,
               y_1 = y_max - (i+1) * delta_y;

        // Using bilinear interpolation on the celll
        double h = data[(i + 0)*num_points_x + j + 0] * (x_1 - x)/delta_x * (y - y_1)/delta_y   // (x0, y0)
                 + data[(i + 0)*num_points_x + j + 1] * (x - x_0)/delta_x * (y - y_1)/delta_y   // (x1, y0)
                 + data[(i + 1)*num_points_x + j + 0] * (x_1 - x)/delta_x * (y_0 - y)/delta_y   // (x0, y1)
                 + data[(i + 1)*num_points_x + j + 1] * (x - x_0)/delta_x * (y_0 - y)/delta_y;  // (x1, y1)

        return h;
    }

    PointList raster_points() const;
    SimplePolygon exterior() const;

    template<typename P>
    MultiPolygon compute_intersection(const P& polygon) const;

    bool contains(double x, double y) const;
};
} // namespace traits

template<template<typename> class RasterData, typename FT>
void bind_raster_list(py::module &m, const std::string& pyname) {
    py::class_<std::vector<RasterData<FT>>, std::unique_ptr<std::vector<RasterData<FT>>>> (m, pyname.c_str())
        .def(py::init( [] () {std::vector<RasterData<FT>> self; return self;}))
        .def("add_raster",
            [] (std::vector<RasterData<FT>>& self, RasterData<FT> raster_data) {
                self.push_back(raster_data);
            }, py::keep_alive<1,2>())
        .def("__getitem__",
            [] (std::vector<RasterData<FT>>& self, int index) {
                return self.at(index);
            }, py::return_value_policy::reference_internal);
}

template<template<typename> class RasterData, typename FT>
void bind_rasterdata(py::module &m, const std::string& pyname) {
    py::class_<RasterData<FT>, std::unique_ptr<RasterData<FT>>>(m, pyname.c_str(), py::buffer_protocol())
    .def(py::init([] (py::array_t<FT>& data_array, double x_min, double y_max, double delta_x, double delta_y) {
            auto buffer = data_array.request();
            int m = buffer.shape[0], n = buffer.shape[1];

            return RasterData<FT>(x_min, y_max, delta_x, delta_y, n, m, static_cast<FT*>(buffer.ptr) );
        }), py::return_value_policy::take_ownership,  py::keep_alive<1, 2>(),
            py::arg("data_array").noconvert(), py::arg("x_min"), py::arg("y_max"), py::arg("delta_x"), py::arg("delta_y"))
    .def_buffer([] (RasterData<FT>& self) {
            return py::buffer_info(
                self.data,
                sizeof(FT),
                py::format_descriptor<FT>::format(),
                2,
                std::vector<std::size_t> { self.num_points_y, self.num_points_x },
                { sizeof(FT) * self.num_points_x, sizeof(FT) }
            );
            })
    .def_readwrite("x_min", &RasterData<FT>::x_min)
    .def_readwrite("y_max", &RasterData<FT>::y_max)
    .def_readwrite("delta_x", &RasterData<FT>::delta_x)
    .def_readwrite("delta_y", &RasterData<FT>::delta_y)
    .def_readonly("num_points_x", &RasterData<FT>::num_points_x)
    .def_readonly("num_points_y", &RasterData<FT>::num_points_y)
    .def_property_readonly("x_max", &RasterData<FT>::get_x_max)
    .def_property_readonly("y_min", &RasterData<FT>::get_y_min)
    .def("__getitem__", [](RasterData<FT>& self, std::pair<int, int> idx) {
                auto [i, j] = idx;
                return self.data[self.num_points_x * i + j]; })
    .def("get_indices", &RasterData<FT>::get_indices)
    .def("exterior", &RasterData<FT>::exterior, py::return_value_policy::take_ownership)
    .def("contains", &RasterData<FT>::contains)
    .def("get_interpolated_value_at_point", &RasterData<FT>::get_interpolated_value_at_point);
}

template<template<typename> class RasterData>
void bind_raster(py::module& m) {
    bind_rasterdata<RasterData, float>(m, "raster_data_float");
    bind_rasterdata<RasterData, double>(m, "raster_data_double");

    bind_raster_list<RasterData, float>(m, "raster_list_float");
    bind_raster_list<RasterData, double>(m, "raster_list_double");
}
} // namespace rasputin
