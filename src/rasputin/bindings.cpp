#include <concepts>

#include <boost/geometry/geometry.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "triangulate_dem.h"

namespace py = pybind11;

namespace rasputin {
template<typename T, std::size_t n>
requires std::floating_point<T>
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
requires std::floating_point<T>
py::buffer_info vector_buffer(std::vector<T> &v) {
    return py::buffer_info(&v[0], sizeof(T), py::format_descriptor<T>::format(), 1, { v.size() }, { sizeof(T) });
}

#define MAKE_OPAQUE(T, Point) (PYBIND11_MAKE_OPAQUE(rasputin::Points<Point>))
#define BIND_MODULE(T, Point) {\
    PYBIND10_MODULE(triangulate_dem, m) {\
        py::bind_vector<rasputin::Points<Point>>(m, "points", py::buffer_protocol())\
            .def_buffer(&vecarray_buffer<T, bg::dimension<Point>::value>);\
        py::bind_vector<std::vector<int>>(m, "int_vector");\
        py::bind_vector<std::vector<std::vector<int>>>(m, "shadow_vector");\
    };\
}

// static constexpr void make_bindings() {
    // PYBIND10_MODULE(triangulate_dem, m) {
        // py::bind_vector<rasputin::Points<Point>>(m, "points", py::buffer_protocol())
            // .def_buffer(&vecarray_buffer<T, bg::dimension<Point>::value>);
        // py::bind_vector<rasputin::point2_vector>(m, "point2_vector", py::buffer_protocol())
            // .def_buffer(&vecarray_buffer<double, 2>);
        // py::bind_vector<rasputin::face_vector>(m, "face_vector", py::buffer_protocol())
            // .def_buffer(&vecarray_buffer<int, 3>);
        // py::bind_vector<rasputin::index_vector>(m, "index_vector", py::buffer_protocol())
            // .def_buffer(&vecarray_buffer<unsigned int, 2>);
        // py::bind_vector<rasputin::double_vector>(m, "double_vector", py::buffer_protocol())
            // .def_buffer(&vector_buffer<double>);

        // py::bind_vector<std::vector<int>>(m, "int_vector");
        // py::bind_vector<std::vector<std::vector<int>>>(m, "shadow_vector");
        // m
            // .def(
                // "lindstrom_turk_by_size",
                // [] (const rasputin::point3_vector& raster_coordinates, size_t result_mesh_size) {
                        // return rasputin::make_tin(raster_coordinates) // use own implementation
                // },
                // "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold."
            // )
            // .def(
                // "lindstrom_turk_by_ratio",
                // [] (const rasputin::point3_vector& raster_coordinates, double ratio) {
                        // return rasputin::make_tin(raster_coordinates) // use own implementation
                // },
                // "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold."
            // )
            // .def("compute_shadow", &rasputin::compute_shadow, "Compute shadows for given sun ray direction.")
            // .def("compute_shadows", &rasputin::compute_shadows, "Compute shadows for a series of times and ray directions.")
            // .def("surface_normals", &rasputin::surface_normals, "Compute surface normals for all faces in the mesh.")
            // .def("point_normals", &rasputin::point_normals, "Compute surface normals for all vertices in the mesh.")
            // .def("orient_tin", &rasputin::orient_tin, "Orient all triangles in the TIN and returns their surface normals.")
            // .def("extract_lakes", &rasputin::extract_lakes, "Extract lakes as separate face list.")
            // .def("compute_slopes", &rasputin::compute_slopes,"Compute slopes (i.e. angles relative to xy plane) for the all the vectors in list.")
            // .def("compute_aspects", &rasputin::compute_aspects, "Compute aspects for the all the vectors in list.")
            // .def("cell_centers", &rasputin::cell_centers, "Compute cell centers for triangulation.")
            // .def("extract_avalanche_expositions", &rasputin::extract_avalanche_expositions, "Extract avalanche exposed cells.")
            // .def("coordinates_to_indices", &rasputin::coordinates_to_indices, "Transform from coordinate space to index space.")
            // .def("extract_uint8_buffer_values", &rasputin::extract_buffer_values<std::uint8_t>, "Extract raster data for given indices.")
            // .def(
                // "rasterdata_to_pointvector",
                // [] (py::array_t<double> array, double x0, double y1, double dx, double dy) {
                     // rasputin::point3_vector raster_coordinates;
                     // auto buffer = array.request();
                     // unsigned long M = (unsigned long)buffer.shape[0];
                     // unsigned long N = (unsigned long)buffer.shape[1];
                     // double* ptr = (double*) buffer.ptr;
    // 
                     // //double dx = (x1 - x0)/(n - 1), dy = (y1 - y0)/(m - 1);
    // 
                     // raster_coordinates.reserve(M*N);
                     // for (std::size_t i = 0; i < M; ++i)
                         // for (std::size_t j = 0; j < N; ++j)
                             // raster_coordinates.emplace_back(std::array<double, 3>{x0 + j*dx, y1 - i*dy, ptr[i*N + j]});
    // 
                     // return raster_coordinates;
                // },
                // "point3_vector from raster data"
            // );
    // }
// }
}
