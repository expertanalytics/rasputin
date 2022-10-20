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

#include "triangulate_dem.h"
#include "mesh.h"

namespace py = pybind11;
namespace SMS = CGAL::Surface_mesh_simplification;

void bind_raster(py::module&);
void bind_mesh(py::module&);
void bind_polygons(py::module&);
void bind_solars(py::module&);
void bind_shades(py::module&);
void bind_slopes(py::module&);

PYBIND11_MAKE_OPAQUE(rasputin::point3_vector);
PYBIND11_MAKE_OPAQUE(rasputin::point2_vector);
PYBIND11_MAKE_OPAQUE(rasputin::face_vector);
PYBIND11_MAKE_OPAQUE(rasputin::double_vector);
PYBIND11_MAKE_OPAQUE(rasputin::index_vector);
PYBIND11_MAKE_OPAQUE(CGAL::MultiPolygon);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<double>>);

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

PYBIND11_MODULE(triangulate_dem, m) {
    py::bind_vector<rasputin::point3_vector>(m, "point3_vector", py::buffer_protocol())
      .def_buffer(&vecarray_buffer<double, 3>)
      .def("from_numpy", &vecarray_from_numpy<double, 3>);
    py::bind_vector<rasputin::point2_vector>(m, "point2_vector", py::buffer_protocol())
      .def_buffer(&vecarray_buffer<double, 2>)
      .def("from_numpy", &vecarray_from_numpy<double, 2>);
    py::bind_vector<rasputin::face_vector>(m, "face_vector", py::buffer_protocol())
      .def_buffer(&vecarray_buffer<int, 3>)
      .def("from_numpy", &vecarray_from_numpy<int, 3>);
    py::bind_vector<rasputin::index_vector>(m, "index_vector", py::buffer_protocol())
      .def_buffer(&vecarray_buffer<unsigned int, 2>)
      .def("from_numpy", &vecarray_from_numpy<unsigned int, 2>);
    py::bind_vector<rasputin::double_vector >(m, "double_vector", py::buffer_protocol())
      .def_buffer(&vector_buffer<double>);

    py::bind_vector<std::vector<int>>(m, "int_vector");
    py::bind_vector<std::vector<std::vector<int>>>(m, "shadow_vector");

    bind_raster(m);

    bind_mesh(m);

    bind_polygons(m);

    py::class_<rasputin::Mesh, std::unique_ptr<rasputin::Mesh>>(m, "Mesh")
        .def("lindstrom_turk_by_ratio",
            [] (const rasputin::Mesh& self, double ratio) {
                return self.coarsen(SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                    SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                    SMS::LindstromTurk_cost<CGAL::Mesh>());
            }, py::return_value_policy::take_ownership,
            "Simplify the mesh.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold.")
        .def("lindstrom_turk_by_size",
            [] (const rasputin::Mesh& self, int max_size) {
                return self.coarsen(SMS::Count_stop_predicate<CGAL::Mesh>(max_size),
                                    SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                    SMS::LindstromTurk_cost<CGAL::Mesh>());
            }, py::return_value_policy::take_ownership,
            "Simplify the mesh.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold.")
        .def("copy", &rasputin::Mesh::copy, py::return_value_policy::take_ownership)
        .def("extract_sub_mesh", &rasputin::Mesh::extract_sub_mesh, py::return_value_policy::take_ownership)

        .def_property_readonly("num_vertices", &rasputin::Mesh::num_vertices)
        .def_property_readonly("num_edges", &rasputin::Mesh::num_edges)
        .def_property_readonly("num_faces", &rasputin::Mesh::num_faces)

        .def_property_readonly("points", &rasputin::Mesh::get_points, py::return_value_policy::reference_internal)
        .def_property_readonly("faces", &rasputin::Mesh::get_faces, py::return_value_policy::reference_internal);

    bind_slopes(m);
    bind_shades(m);
    bind_solars(m);
}
