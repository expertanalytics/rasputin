#include <pybind11/pybind11.h>

#include "slope.h"

namespace py = pybind11;

void bind_slopes(py::module& m) {
     m.def("surface_normals", &rasputin::surface_normals,
          "Compute surface normals for all faces in the mesh.",
          py::return_value_policy::take_ownership)
     .def("point_normals", &rasputin::point_normals, "Compute surface normals for all vertices in the mesh.")
     .def("orient_tin", &rasputin::orient_tin, "Orient all triangles in the TIN and returns their surface normals.")
     .def("extract_lakes", &rasputin::extract_lakes, "Extract lakes as separate face list.")
     .def("compute_slopes", &rasputin::compute_slopes,"Compute slopes (i.e. angles relative to xy plane) for the all the vectors in list.")
     .def("compute_aspects", &rasputin::compute_aspects, "Compute aspects for the all the vectors in list.")
     .def("extract_avalanche_expositions", &rasputin::extract_avalanche_expositions, "Extract avalanche exposed cells.")
     .def("consolidate", &rasputin::consolidate, "Make a stand alone consolidated tin.")
     .def("cell_centers", &rasputin::cell_centers, "Compute cell centers for triangulation.")
     .def("coordinates_to_indices", &rasputin::coordinates_to_indices, "Transform from coordinate space to index space.")
     .def("extract_uint8_buffer_values", &rasputin::extract_buffer_values<std::uint8_t>, "Extract raster data for given indices.");

}
