#include "triangulate_dem.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

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
    m.def("lindstrom_turk_by_size",
          [] (const rasputin::PointList& raster_coordinates, size_t result_mesh_size) {
              rasputin::PointList tin_coordinates;
              rasputin::FaceList tin_faces;
              rasputin::make_tin(raster_coordinates,
                                 SMS::Count_stop_predicate<CGAL::Mesh>(result_mesh_size),
                                 SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                 SMS::LindstromTurk_cost<CGAL::Mesh>(),
                                 tin_coordinates,
                                 tin_faces);
              return std::make_tuple(tin_coordinates, tin_faces);
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold."
    ).def("lindstrom_turk_by_ratio",
          [] (const rasputin::PointList& raster_coordinates, double ratio) {
              rasputin::PointList tin_coordinates;
              rasputin::FaceList tin_faces;
              rasputin::make_tin(raster_coordinates,
                                 SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                 SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                 SMS::LindstromTurk_cost<CGAL::Mesh>(),
                                 tin_coordinates,
                                 tin_faces);
              return std::make_tuple(tin_coordinates, tin_faces);
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold."
    );
}
