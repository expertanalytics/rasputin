#include <vector>

#include <pybind11/pybind11.h>

#include "raster_data.h"
#include "mesh.h"

namespace py = pybind11;

template<typename R, typename P>
void bind_make_mesh(py::module &m) {
    m.def("make_mesh",
        [] (const R& raster_data, const P polygon, const std::string proj4_str) {
            return rasputin::mesh_from_raster(raster_data, polygon, proj4_str);
        }, py::return_value_policy::take_ownership)
    .def("make_mesh",
        [] (const R& raster_data, const std::string proj4_str) {
            return rasputin::mesh_from_raster(raster_data, proj4_str);
        }, py::return_value_policy::take_ownership);
}

void bind_mesh(py::module& m) {
    bind_make_mesh<std::vector<rasputin::RasterData<float>>, CGAL::SimplePolygon>(m);
    bind_make_mesh<std::vector<rasputin::RasterData<double>>, CGAL::SimplePolygon>(m);
    bind_make_mesh<std::vector<rasputin::RasterData<float>>, CGAL::Polygon>(m);
    bind_make_mesh<std::vector<rasputin::RasterData<double>>, CGAL::Polygon>(m);
    // bind_make_mesh<std::vector<rasputin::RasterData<float>>, CGAL::MultiPolygon>(m);
    // bind_make_mesh<std::vector<rasputin::RasterData<double>>, CGAL::MultiPolygon>(m);
    bind_make_mesh<rasputin::RasterData<float>, CGAL::SimplePolygon>(m);
    bind_make_mesh<rasputin::RasterData<double>, CGAL::SimplePolygon>(m);
    bind_make_mesh<rasputin::RasterData<float>, CGAL::Polygon>(m);
    bind_make_mesh<rasputin::RasterData<double>, CGAL::Polygon>(m);
    // bind_make_mesh<rasputin::RasterData<float>, CGAL::MultiPolygon>(m);
    // bind_make_mesh<rasputin::RasterData<double>, CGAL::MultiPolygon>(m);

     m.def("construct_mesh",
            [] (const rasputin::point3_vector& points, const rasputin::face_vector & faces, const std::string proj4_str) {
                rasputin::VertexIndexMap index_map;
                rasputin::FaceDescrMap face_map;
                return rasputin::Mesh(rasputin::construct_mesh(points, faces, index_map, face_map), proj4_str);
         }, py::return_value_policy::take_ownership);
}
