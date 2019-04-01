#include "triangulate_dem.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

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

template<typename T>
py::buffer_info vector_buffer(std::vector<T> &v) {
    return py::buffer_info(&v[0], sizeof(T), py::format_descriptor<T>::format(), 1, { v.size() }, { sizeof(T) });
}


PYBIND11_MODULE(triangulate_dem, m) {
    py::bind_vector<rasputin::PointList>(m, "PointVector", py::buffer_protocol())
        .def_buffer(&vecarray_buffer<double, 3>);
    py::bind_vector<rasputin::FaceList>(m, "FaceVector", py::buffer_protocol())
        .def_buffer(&vecarray_buffer<int, 3>);
    py::bind_vector<rasputin::ScalarList>(m, "ScalarVector", py::buffer_protocol())
        .def_buffer(&vector_buffer<double>);

    py::bind_vector<std::vector<int>>(m, "IntVector");
    py::bind_vector<std::vector<std::vector<int>>>(m, "ShadowVector");

    m.def("lindstrom_turk_by_size",
          [] (const rasputin::PointList& raster_coordinates, size_t result_mesh_size) {
              return rasputin::make_tin(raster_coordinates,
                                        SMS::Count_stop_predicate<CGAL::Mesh>(result_mesh_size),
                                        SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                        SMS::LindstromTurk_cost<CGAL::Mesh>());
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold.")
        .def("lindstrom_turk_by_ratio",
             [] (const rasputin::PointList& raster_coordinates, double ratio) {
                 return rasputin::make_tin(raster_coordinates,
                                           SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                           SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                           SMS::LindstromTurk_cost<CGAL::Mesh>());
             },
            "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold.")
        .def("compute_shadow", &rasputin::compute_shadow, "Compute shadows for given sun ray direction.")
        .def("compute_shadows", &rasputin::compute_shadows, "Compute shadows for a series of times and ray directions.")
        .def("surface_normals", &rasputin::surface_normals, "Compute surface normals for all faces in the mesh.")
        .def("point_normals", &rasputin::point_normals, "Compute surface normals for all vertices in the mesh.")
        .def("orient_tin", &rasputin::orient_tin, "Orients all triangles in the TIN and returns their surface normals.")
        .def("extract_lakes", &rasputin::extract_lakes, "Extract lakes as separate face list.")
        .def("compute_slopes", &rasputin::compute_slopes,"computes slopes (i.e. angles relative to xy plane) for the all the vectors in list.")
        .def("compute_aspect", &rasputin::compute_aspect,"computes aspects for the all the vectors in list.")
        .def("rasterdata_to_pointvector",
             [] (py::array_t<double> array, double x0, double y0, double x1, double y1) {
                 auto buffer = array.request();
                 int m = buffer.shape[0], n = buffer.shape[1];
                 double* ptr = (double*) buffer.ptr;

                 double dx = (x1 - x0)/(n - 1), dy = (y1 - y0)/(m - 1);

                 rasputin::PointList raster_coordinates;
                 raster_coordinates.reserve(m*n);
                 for (std::size_t i = 0; i < m; ++i)
                     for (std::size_t j = 0; j < n; ++j)
                         raster_coordinates.push_back(std::array<double, 3>{x0 + j*dx, y1 - i*dy, ptr[i*n + j]});
                 return raster_coordinates;
            }, "Pointvector from raster data");
}
