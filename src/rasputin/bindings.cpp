#include "triangulate_dem.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <memory.h>

// Surface mesh simplication policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>

namespace py = pybind11;
namespace SMS = CGAL::Surface_mesh_simplification;

PYBIND11_MAKE_OPAQUE(rasputin::PointList);
PYBIND11_MAKE_OPAQUE(rasputin::PointList2D);
PYBIND11_MAKE_OPAQUE(rasputin::FaceList);
PYBIND11_MAKE_OPAQUE(rasputin::ScalarList);

template<typename T, std::size_t n>
py::buffer_info vecarray_buffer(std::vector<std::array<T, n>> &v) {
    return py::buffer_info(
        &v[0],                                        /* Pointer to buffer */
        sizeof(T),                                    /* Size of one scalar */
        py::format_descriptor<T>::format(),           /* Python struct-style format descriptor */
        2,                                            /* Number of dimensions */
        std::vector<std::size_t> { v.size(), n },     /* Buffer dimensions */
        { sizeof(T) * n, sizeof(T) }                  /* Strides (in bytes) for each index */
    );
}

template<typename T, std::size_t n>
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

template<typename T>
py::buffer_info vector_buffer(std::vector<T> &v) {
    return py::buffer_info(&v[0], sizeof(T), py::format_descriptor<T>::format(), 1, { v.size() }, { sizeof(T) });
}


template<typename FT>
void bind_rasterdata(py::module &m, const std::string& pyname) {
    py::class_<rasputin::RasterData<FT>, std::shared_ptr<rasputin::RasterData<FT>>>(m, pyname.c_str(), py::buffer_protocol())
    .def(py::init([] (py::array_t<FT>& data_array, double x_min, double y_max, double delta_x, double delta_y) {
            auto buffer = data_array.request();
            int m = buffer.shape[0], n = buffer.shape[1];

            return rasputin::RasterData<FT>(x_min, y_max, delta_x, delta_y, n, m, (FT*) buffer.ptr );
        }))
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
                return self.data[self.num_points_x * i, j]; })
    .def("get_indices", &rasputin::RasterData<FT>::get_indices)
    .def("get_interpolated_value_at_point", &rasputin::RasterData<FT>::get_interpolated_value_at_point);

    m.def("lindstrom_turk_by_size",
          [] (const rasputin::RasterData<FT>& raster_data, const rasputin::PointList2D& boundary_vertices, size_t result_mesh_size) {
              return rasputin::tin_from_raster(raster_data, boundary_vertices,
                                        SMS::Count_stop_predicate<CGAL::Mesh>(result_mesh_size),
                                        SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                        SMS::LindstromTurk_cost<CGAL::Mesh>());
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold.")
     .def("lindstrom_turk_by_size",
          [] (const rasputin::RasterData<FT>& raster_data, size_t result_mesh_size) {
              return rasputin::tin_from_raster(raster_data,
                                        SMS::Count_stop_predicate<CGAL::Mesh>(result_mesh_size),
                                        SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                        SMS::LindstromTurk_cost<CGAL::Mesh>());
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold.")
        .def("lindstrom_turk_by_ratio",
             [] (const rasputin::RasterData<FT>& raster_data, const rasputin::PointList2D& boundary_vertices, double ratio) {
                 return rasputin::tin_from_raster(raster_data, boundary_vertices,
                                           SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                           SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                           SMS::LindstromTurk_cost<CGAL::Mesh>());
             },
            "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold.")
        .def("lindstrom_turk_by_ratio",
             [] (const rasputin::RasterData<FT>& raster_data, double ratio) {
                 return rasputin::tin_from_raster(raster_data,
                                           SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                           SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                           SMS::LindstromTurk_cost<CGAL::Mesh>());
             },
            "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold.");
}

PYBIND11_MODULE(triangulate_dem, m) {
    py::bind_vector<rasputin::PointList>(m, "PointVector", py::buffer_protocol())
      .def_buffer(&vecarray_buffer<double, 3>)
      .def("from_numpy", &vecarray_from_numpy<double, 3>);
    py::bind_vector<rasputin::PointList2D>(m, "PointVector2D", py::buffer_protocol())
      .def_buffer(&vecarray_buffer<double, 2>)
      .def("from_numpy", &vecarray_from_numpy<double, 2>);
    py::bind_vector<rasputin::FaceList>(m, "FaceVector", py::buffer_protocol())
      .def_buffer(&vecarray_buffer<int, 3>)
      .def("from_numpy", &vecarray_from_numpy<int, 3>);
    py::bind_vector<rasputin::ScalarList>(m, "ScalarVector", py::buffer_protocol())
      .def_buffer(&vector_buffer<double>);

    bind_rasterdata<float>(m, "RasterData_float");
    bind_rasterdata<double>(m, "RasterData_double");

    py::bind_vector<std::vector<int>>(m, "IntVector");
    py::bind_vector<std::vector<std::vector<int>>>(m, "ShadowVector");

    m.def("compute_shadow", &rasputin::compute_shadow, "Compute shadows for given sun ray direction.")
     .def("compute_shadows", &rasputin::compute_shadows, "Compute shadows for a series of times and ray directions.")
     .def("surface_normals", &rasputin::surface_normals, "Compute surface normals for all faces in the mesh.")
     .def("point_normals", &rasputin::point_normals, "Compute surface normals for all vertices in the mesh.")
     .def("orient_tin", &rasputin::orient_tin, "Orient all triangles in the TIN and returns their surface normals.")
     .def("extract_lakes", &rasputin::extract_lakes, "Extract lakes as separate face list.")
     .def("compute_slopes", &rasputin::compute_slopes,"Compute slopes (i.e. angles relative to xy plane) for the all the vectors in list.")
     .def("compute_aspects", &rasputin::compute_aspects, "Compute aspects for the all the vectors in list.")
     .def("extract_avalanche_expositions", &rasputin::extract_avalanche_expositions, "Extract avalanche exposed cells.");
}
