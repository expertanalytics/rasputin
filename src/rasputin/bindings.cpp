#include "triangulate_dem.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>

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


PYBIND11_MODULE(triangulate_dem, m) {
    py::bind_vector<rasputin::PointList>(m, "PointVector");
    py::bind_vector<rasputin::FaceList>(m, "FaceVector");
    py::bind_vector<std::vector<int>>(m, "IntVector");
    py::bind_vector<std::vector<std::vector<int>>>(m, "ShadowVector");
    m.def("lindstrom_turk_by_size",
          [] (const rasputin::PointList& raster_coordinates, size_t result_mesh_size) {
              return rasputin::make_tin(raster_coordinates,
                                        SMS::Count_stop_predicate<CGAL::Mesh>(result_mesh_size),
                                        SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                        SMS::LindstromTurk_cost<CGAL::Mesh>());
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold."
    ).def("lindstrom_turk_by_ratio",
          [] (const rasputin::PointList& raster_coordinates, double ratio) {
              return rasputin::make_tin(raster_coordinates,
                                        SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                        SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                        SMS::LindstromTurk_cost<CGAL::Mesh>());
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold."
    ).def("compute_shadow",
          &rasputin::compute_shadow,
          "Compute shadows for given sun ray direction."
    ).def("compute_shadows", &rasputin::compute_shadows,
          "Compute shadows for a series of times and ray directions."
    ).def("surface_normals", &rasputin::surface_normals,
            "Compute surface normals for all faces in the mesh.");
}
