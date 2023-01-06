#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

// Surface mesh simplication policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

#include "cgal_mesh.h"
#include "cgal_raster_data.h"

namespace py = pybind11;
namespace SMS = CGAL::Surface_mesh_simplification;

PYBIND11_MAKE_OPAQUE(rasputin::point3_vector);
PYBIND11_MAKE_OPAQUE(rasputin::point2_vector);
PYBIND11_MAKE_OPAQUE(rasputin::face_vector);
// PYBIND11_MAKE_OPAQUE(CGAL::MultiPolygon);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<rasputin::RasterData<double>>);

void bind_mesh(py::module& m) {
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

    rasputin::bind_make_mesh<
        rasputin::RasterData,
        float,
        CGAL::SimplePolygon,
        CGAL::PointList,
        rasputin::Mesh,
        CGAL::Mesh,
        CGAL::DelaunayConstraints
    >(m);
    rasputin::bind_make_mesh<
        rasputin::RasterData,
        double,
        CGAL::SimplePolygon,
        CGAL::PointList,
        rasputin::Mesh,
        CGAL::Mesh,
        CGAL::DelaunayConstraints
    >(m);
    rasputin::bind_make_mesh<
        rasputin::RasterData,
        float,
        CGAL::Polygon,
        CGAL::PointList,
        rasputin::Mesh,
        CGAL::Mesh,
        CGAL::DelaunayConstraints
    >(m);
    rasputin::bind_make_mesh<
        rasputin::RasterData,
        double,
        CGAL::Polygon,
        CGAL::PointList,
        rasputin::Mesh,
        CGAL::Mesh,
        CGAL::DelaunayConstraints
    >(m);

    // bind_make_mesh<std::vector<rasputin::RasterData<float>>, CGAL::MultiPolygon>(m);
    // bind_make_mesh<std::vector<rasputin::RasterData<double>>, CGAL::MultiPolygon>(m);
    // bind_make_mesh<rasputin::RasterData<float>, CGAL::MultiPolygon>(m);
    // bind_make_mesh<rasputin::RasterData<double>, CGAL::MultiPolygon>(m);

     m.def("construct_mesh",
        [] (const rasputin::point3_vector& points, const rasputin::face_vector & faces, const std::string proj4_str) {
            rasputin::VertexIndexMap index_map;
            rasputin::FaceDescrMap face_map;
            return rasputin::Mesh(rasputin::Mesh::construct_mesh(points, faces, index_map, face_map), proj4_str);
        }, py::return_value_policy::take_ownership);
}
