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
PYBIND11_MAKE_OPAQUE(rasputin::point2_vector);

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

rasputin::PointList rasterdata_to_pointvector(py::array_t<double> array, double x0, double y0, double x1, double y1) {
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
            }

rasputin::PointList interpolate_onto_polygon(py::array_t<double> array, double x0, double y0, double x1, double y1, rasputin::PointList polygon) {
    auto buffer = array.request();
    int m = buffer.shape[0], n = buffer.shape[1];
    double dx = (x1 - x0)/(n - 1), dy = (y1 - y0)/(m - 1);
    double* ptr = (double*) buffer.ptr;

    auto get_index= [x0,y1,dx,dy](double x, double y) -> std::pair<int, int> {
        int j = static_cast<int>((x-x0) /dx);
        int i = static_cast<int>((y1-y) /dy);

        return std::make_pair(i,j);
    };

    rasputin::PointList interpolated_points;
    for (std::size_t vertex_number = 0; vertex_number < polygon.size(); ++vertex_number) {
        rasputin::Point vertex = polygon[vertex_number];

        rasputin::Point second_vertex = polygon[(vertex_number + 1) % polygon.size()];

        // Sample with the same resolution along the polygon edges
        double edge_len_x = second_vertex[0] - vertex[0];
        double edge_len_y = second_vertex[1] - vertex[1];

        double edge_len = pow(pow(edge_len_x, 2) + pow(edge_len_y, 2), 0.5);

        // Not including next vertex
        std::size_t num_samples = static_cast<int>(std::max<double>(std::fabs(edge_len_x/dx), std::fabs(edge_len_y/dy)));
        double edge_dx = edge_len_x / num_samples;
        double edge_dy = edge_len_y / num_samples;

        for (std::size_t k=0; k < num_samples; ++k) {
            double x = vertex[0] + k * edge_dx;
            double y = vertex[1] + k * edge_dy;
            auto indices = get_index(x, y);

            int i = indices.first, j = indices.second;

            // Bilinear interpolation
            // h = h_i_j0 * (x_i1 - x)/dx
            //   + h_i_j1 * (x - x_i0)/dx
            //
            //   = h_i0_j0 * (x_i1 -x)/dx * (yi1 - y)/dy
            //   + h_i1_j0 * (x_i1 -x)/dx * (y - yi0)/dy
            //   + ...
            //   yi0 = y1 -i * dy
            //   yi1 = y1 -(i+1) * dy
            double x_j0 = x0 + (j+0) * dx;
            double x_j1 = x0 + (j+1) * dx;
            double y_i0 = y1 - (i+0) * dy;
            double y_i1 = y1 - (i+1) * dy;
            /* double h = ptr[(i + 0)*n + j+0] * ((j+1)*dx  +x0 - x)/dx * -(y1 - (i+1) * dy - y)/dy */
            /*          + ptr[(i + 1)*n + j+0] * ((j+1)*dx  +x0 - x)/dx * -(y - y1 + (i+0) * dy)/dy */
            /*          + ptr[(i + 0)*n + j+1] * (x - (j+0)*dx - x0)/dx * -(y1 - (i+1) * dy - y)/dy */
            /*          + ptr[(i + 0)*n + j+1] * (x - (j+0)*dx - x0)/dx * -(y - y1 + (i+0) * dy)/dy; */

            double h = ptr[(i + 0)*n + j + 0] * (x_j1 - x)/dx * (y - y_i1)/dy
                     + ptr[(i + 1)*n + j + 0] * (x_j1 - x)/dx * (y_i0 - y)/dy
                     + ptr[(i + 0)*n + j + 1] * (x - x_j0)/dx * (y - y_i1)/dy
                     + ptr[(i + 1)*n + j + 1] * (x - x_j0)/dx * (y_i0 - y)/dy;

            interpolated_points.push_back(std::array<double, 3>{x, y, h});


        }

    }

    return interpolated_points;

}


PYBIND11_MODULE(triangulate_dem, m) {
    py::bind_vector<rasputin::PointList>(m, "PointVector", py::buffer_protocol())
        .def_buffer(&vecarray_buffer<double, 3>);
    py::bind_vector<rasputin::point2_vector>(m, "point2_vector", py::buffer_protocol())
        .def_buffer(&vecarray_buffer<double, 2>);
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
        .def("orient_tin", &rasputin::orient_tin, "Orient all triangles in the TIN and returns their surface normals.")
        .def("extract_lakes", &rasputin::extract_lakes, "Extract lakes as separate face list.")
        .def("compute_slopes", &rasputin::compute_slopes,"Compute slopes (i.e. angles relative to xy plane) for the all the vectors in list.")
        .def("compute_aspects", &rasputin::compute_aspects, "Compute aspects for the all the vectors in list.")
        .def("extract_avalanche_expositions", &rasputin::extract_avalanche_expositions, "Extract avalanche exposed cells.")
        .def("constrained2", &interpolate_onto_polygon)
        .def("constrained3",
             [] (py::array_t<double> array, double x0, double y0, double x1, double y1, rasputin::PointList boundary, double ratio) {
             rasputin::PointList pts = rasterdata_to_pointvector(array, x0, y0, x1, y1);
             rasputin::PointList boundary_pts = interpolate_onto_polygon(array, x0, y0, x1, y1, boundary);
             return rasputin::constrained(pts, boundary, boundary_pts,
                                          SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                          SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                          SMS::LindstromTurk_cost<CGAL::Mesh>());
             })
        .def("rasterdata_to_pointvector", &rasterdata_to_pointvector, "Pointvector from raster data");
}
