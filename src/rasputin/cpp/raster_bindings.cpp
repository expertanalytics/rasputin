#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "raster_data.h"

namespace py = pybind11;

template<typename T>
void bind_raster_list(py::module &m, const std::string& pyname) {
    py::class_<std::vector<rasputin::RasterData<T>>, std::unique_ptr<std::vector<rasputin::RasterData<T>>>> (m, pyname.c_str())
        .def(py::init( [] () {std::vector<rasputin::RasterData<T>> self; return self;}))
        .def("add_raster",
            [] (std::vector<rasputin::RasterData<T>>& self, rasputin::RasterData<T> raster_data) {
                self.push_back(raster_data);
            }, py::keep_alive<1,2>())
        .def("__getitem__",
            [] (std::vector<rasputin::RasterData<T>>& self, int index) {
                return self.at(index);
            }, py::return_value_policy::reference_internal);
}

template<typename FT>
void bind_rasterdata(py::module &m, const std::string& pyname) {
    py::class_<rasputin::RasterData<FT>, std::unique_ptr<rasputin::RasterData<FT>>>(m, pyname.c_str(), py::buffer_protocol())
    .def(py::init([] (py::array_t<FT>& data_array, double x_min, double y_max, double delta_x, double delta_y) {
            auto buffer = data_array.request();
            int m = buffer.shape[0], n = buffer.shape[1];

            return rasputin::RasterData<FT>(x_min, y_max, delta_x, delta_y, n, m, static_cast<FT*>(buffer.ptr) );
        }), py::return_value_policy::take_ownership,  py::keep_alive<1, 2>(),
            py::arg("data_array").noconvert(), py::arg("x_min"), py::arg("y_max"), py::arg("delta_x"), py::arg("delta_y"))
    .def_buffer([] (rasputin::RasterData<FT>& self) {
            return py::buffer_info(
                self.data,
                sizeof(FT),
                py::format_descriptor<FT>::format(),
                2,
                std::vector<std::size_t> { self.num_points_y, self.num_points_x },
                { sizeof(FT) * self.num_points_x, sizeof(FT) }
            );
            })
    .def_readwrite("x_min", &rasputin::RasterData<FT>::x_min)
    .def_readwrite("y_max", &rasputin::RasterData<FT>::y_max)
    .def_readwrite("delta_x", &rasputin::RasterData<FT>::delta_x)
    .def_readwrite("delta_y", &rasputin::RasterData<FT>::delta_y)
    .def_readonly("num_points_x", &rasputin::RasterData<FT>::num_points_x)
    .def_readonly("num_points_y", &rasputin::RasterData<FT>::num_points_y)
    .def_property_readonly("x_max", &rasputin::RasterData<FT>::get_x_max)
    .def_property_readonly("y_min", &rasputin::RasterData<FT>::get_y_min)
    .def("__getitem__", [](rasputin::RasterData<FT>& self, std::pair<int, int> idx) {
                auto [i, j] = idx;
                return self.data[self.num_points_x * i + j]; })
    .def("get_indices", &rasputin::RasterData<FT>::get_indices)
    .def("exterior", &rasputin::RasterData<FT>::exterior, py::return_value_policy::take_ownership)
    .def("contains", &rasputin::RasterData<FT>::contains)
    .def("get_interpolated_value_at_point", &rasputin::RasterData<FT>::get_interpolated_value_at_point);
}

void bind_raster(py::module& m) {
    bind_rasterdata<float>(m, "raster_data_float");
    bind_rasterdata<double>(m, "raster_data_double");

    bind_raster_list<float>(m, "raster_list_float");
    bind_raster_list<double>(m, "raster_list_double");
}
