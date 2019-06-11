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


template<typename P0, typename P1>
CGAL::MultiPolygon difference_polygons(const P0& polygon0, const P1& polygon1) {
    CGAL::MultiPolygon difference_result;
    CGAL::difference(polygon0, polygon1, std::back_inserter(difference_result));

    return difference_result;
}

template<typename P0, typename P1>
CGAL::MultiPolygon intersect_polygons(const P0& polygon0, const P1& polygon1) {
    CGAL::MultiPolygon intersection_result;
    CGAL::intersection(polygon0, polygon1, std::back_inserter(intersection_result));

    return intersection_result;
}

template<typename P0, typename P1>
CGAL::MultiPolygon join_polygons(const P0& polygon0, const P1& polygon1) {
    CGAL::Polygon joined;
    CGAL::MultiPolygon join_result;
    if (CGAL::join(polygon0, polygon1, joined))
        join_result.push_back(std::move(joined));
    else {
        join_result.push_back(CGAL::Polygon(polygon0));
        join_result.push_back(CGAL::Polygon(polygon1));
    }
    return join_result;
}


CGAL::MultiPolygon join_multipolygons(const CGAL::MultiPolygon& polygon0, const CGAL::MultiPolygon& polygon1) {
    CGAL::MultiPolygon join_result;
    CGAL::join(polygon0.begin(), polygon0.end(), polygon1.begin(), polygon1.end(), std::back_inserter(join_result));

    return join_result;
}



CGAL::SimplePolygon polygon_from_numpy(py::array_t<double>& buf) {
        auto info = buf.request();
        if (info.ndim < 1 or info.ndim > 2)
            throw py::type_error("Can only convert one- and two-dimensional arrays.");

        // Make sure total size can be preserved
        auto size = info.shape[0];
        if (info.ndim == 1 and size % 2)
            throw py::type_error("Size of one-dimensional array must be divisible by 2.");

        if (info.ndim == 2 and info.shape[1] != 2)
            throw py::type_error("Second dimension does not have size equal to 2.");

        if (info.ndim == 2)
            size *= 2;

        CGAL::SimplePolygon exterior;

        // Copy values
        double* p = static_cast<double*>(info.ptr);
        double* end = p + size;
        for (; p < end; p += 2)
            exterior.push_back(CGAL::Point2(*p, *(p+1)));
        return exterior;
}

template<typename FT>
void bind_rasterdata(py::module &m, const std::string& pyname) {
    py::class_<rasputin::RasterData<FT>, std::unique_ptr<rasputin::RasterData<FT>>>(m, pyname.c_str(), py::buffer_protocol())
    .def(py::init([] (py::array_t<FT>& data_array, double x_min, double y_max, double delta_x, double delta_y) {
            auto buffer = data_array.request();
            int m = buffer.shape[0], n = buffer.shape[1];

            return rasputin::RasterData<FT>(x_min, y_max, delta_x, delta_y, n, m, static_cast<FT*>(buffer.ptr) );
        }), py::return_value_policy::take_ownership,  py::keep_alive<1, 2>(),
            py::arg("data_array").noconvert(), py::arg("x_min"), py::arg("y_max"), py::arg("delta_x"), py::arg("delta_y"))
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
                return self.data[self.num_points_x * i + j]; })
    .def("get_indices", &rasputin::RasterData<FT>::get_indices)
    .def("get_interpolated_value_at_point", &rasputin::RasterData<FT>::get_interpolated_value_at_point);
}

template<typename R>
void bind_tin_from_raster(py::module &m) {
    m.def("lindstrom_turk_by_size",
          [] (const R& raster_data, size_t result_mesh_size) {
              return rasputin::tin_from_raster(raster_data,
                                        SMS::Count_stop_predicate<CGAL::Mesh>(result_mesh_size),
                                        SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                        SMS::LindstromTurk_cost<CGAL::Mesh>());
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold.")
     .def("lindstrom_turk_by_ratio",
        [] (const R& raster_data, double ratio) {
            return rasputin::tin_from_raster(raster_data,
                                           SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                           SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                           SMS::LindstromTurk_cost<CGAL::Mesh>());
             },
            "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold.");
}

template<typename R, typename P>
void bind_tin_from_raster(py::module &m) {
    m.def("lindstrom_turk_by_size",
          [] (const R& raster_data, const P& boundary_polygon, size_t result_mesh_size) {
              return rasputin::tin_from_raster(raster_data, boundary_polygon,
                                        SMS::Count_stop_predicate<CGAL::Mesh>(result_mesh_size),
                                        SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                        SMS::LindstromTurk_cost<CGAL::Mesh>());
          },
          "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the size threshold.")
        .def("lindstrom_turk_by_ratio",
             [] (const R& raster_data, const P& boundary_polygon, double ratio) {
                 return rasputin::tin_from_raster(raster_data, boundary_polygon,
                                           SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                           SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                           SMS::LindstromTurk_cost<CGAL::Mesh>());
             },
            "Construct a TIN based on the points provided.\n\nThe LindstromTurk cost and placement strategy is used, and simplification process stops when the number of undirected edges drops below the ratio threshold.");
}

template<typename R, typename P>
void bind_make_mesh(py::module &m) {
        m.def("make_mesh",
            [] (const R& raster_data, const P polygon) {
                return rasputin::mesh_from_raster(raster_data, polygon);
            })
        .def("make_mesh",
            [] (const R& raster_data) {
                return rasputin::mesh_from_raster(raster_data);
            });
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

    py::class_<CGAL::SimplePolygon, std::unique_ptr<CGAL::SimplePolygon>>(m, "simple_polygon")
        .def(py::init(&polygon_from_numpy))
        .def("num_vertices", &CGAL::SimplePolygon::size)
        .def("array", [] (const CGAL::SimplePolygon& self) {
                py::array_t<double> result(self.size() * 2);
                result.resize(py::array::ShapeContainer({static_cast<long int>(self.size()), 2}), true);
                auto info = result.request();
                double* data = static_cast<double*>(info.ptr);
                for (auto v = self.vertices_begin(); v != self.vertices_end(); ++v) {
                    data[0] = v->x();
                    data[1] = v->y();
                    data += 2;
                }
                return result;
                }, py::return_value_policy::move)
        .def("join", &join_polygons<CGAL::SimplePolygon, CGAL::SimplePolygon>)
        .def("join", &join_polygons<CGAL::SimplePolygon, CGAL::Polygon>)
        .def("difference", &difference_polygons<CGAL::SimplePolygon, CGAL::SimplePolygon>)
        .def("difference", &difference_polygons<CGAL::SimplePolygon, CGAL::Polygon>)
        .def("intersection", &intersect_polygons<CGAL::SimplePolygon, CGAL::SimplePolygon>)
        .def("intersection", &intersect_polygons<CGAL::SimplePolygon, CGAL::Polygon>);

    py::class_<CGAL::Polygon, std::unique_ptr<CGAL::Polygon>>(m, "Polygon")
        .def(py::init([] (py::array_t<double>& buf) {
            const CGAL::SimplePolygon exterior = polygon_from_numpy(buf);
            return CGAL::Polygon(exterior);}))
        .def("holes", [] (const CGAL::Polygon& self) {
                    py::list result;
                    for (auto h = self.holes_begin(); h != self.holes_end(); ++h)
                        result.append(*h);
                    return result;
                    })
        .def("exterior", [] (const CGAL::Polygon& self) {return self.outer_boundary();})
        .def("join", &join_polygons<CGAL::Polygon, CGAL::SimplePolygon>)
        .def("join", &join_polygons<CGAL::Polygon, CGAL::Polygon>)
        .def("difference", &difference_polygons<CGAL::Polygon, CGAL::SimplePolygon>)
        .def("difference", &difference_polygons<CGAL::Polygon, CGAL::Polygon>)
        .def("intersection", &intersect_polygons<CGAL::Polygon, CGAL::SimplePolygon>)
        .def("intersection", &intersect_polygons<CGAL::Polygon, CGAL::Polygon>);

    py::class_<CGAL::MultiPolygon, std::unique_ptr<CGAL::MultiPolygon>>(m, "MultiPolygon")
        .def(py::init(
            [] (const CGAL::Polygon& polygon) {
                CGAL::MultiPolygon self;
                self.push_back(CGAL::Polygon(polygon));
                return self;
            }))
        .def(py::init(
            [] (const CGAL::SimplePolygon& polygon) {
                CGAL::MultiPolygon self;
                self.push_back(CGAL::Polygon(polygon));
                return self;
            }))
        .def("num_parts", &CGAL::MultiPolygon::size)
        .def("parts",
            [] (const CGAL::MultiPolygon& self) {
                py::list result;
                for (auto p = self.begin(); p != self.end(); ++p)
                    result.append(*p);
                return result;
            })
        .def("join",
            [] (const CGAL::MultiPolygon a, CGAL::SimplePolygon b) {
                return join_multipolygons(a, CGAL::MultiPolygon({static_cast<CGAL::Polygon>(b)}));
            })
        .def("join",
            [] (const CGAL::MultiPolygon a, CGAL::Polygon b) {
                return join_multipolygons(a, CGAL::MultiPolygon({b}));
            })
        .def("join",
            [] (const CGAL::MultiPolygon a, CGAL::MultiPolygon b) {
                return join_multipolygons(a, b);
            })
        .def("__getitem__", [] (const CGAL::MultiPolygon& self, int idx) {
                    return self.at(idx);
                    }, py::return_value_policy::reference_internal);

    bind_rasterdata<float>(m, "raster_data_float");
    bind_rasterdata<double>(m, "raster_data_double");

    bind_tin_from_raster<rasputin::RasterData<float>>(m);
    bind_tin_from_raster<rasputin::RasterData<double>>(m);
    bind_tin_from_raster<rasputin::RasterData<float>, CGAL::SimplePolygon>(m);

    bind_make_mesh<std::vector<rasputin::RasterData<float>>, CGAL::SimplePolygon>(m);
    bind_make_mesh<std::vector<rasputin::RasterData<double>>, CGAL::SimplePolygon>(m);
    bind_make_mesh<rasputin::RasterData<float>, CGAL::SimplePolygon>(m);
    bind_make_mesh<rasputin::RasterData<double>, CGAL::SimplePolygon>(m);

    // Polygons with hole not yet supported
    // bind_tin_from_raster<rasputin::RasterData<float>, CGAL::Polygon>(m);
    bind_tin_from_raster<rasputin::RasterData<double>, CGAL::SimplePolygon>(m);
    // bind_tin_from_raster<rasputin::RasterData<double>, CGAL::Polygon>(m);
    bind_tin_from_raster<std::vector<rasputin::RasterData<float>>, CGAL::SimplePolygon>(m);
    // bind_tin_from_raster<std::vector<rasputin::RasterData<float>>, CGAL::Polygon>(m);
    bind_tin_from_raster<std::vector<rasputin::RasterData<double>>, CGAL::SimplePolygon>(m);
    // bind_tin_from_raster<std::vector<rasputin::RasterData<double>>, CGAL::Polygon>(m);

    py::class_<rasputin::Mesh, std::unique_ptr<rasputin::Mesh>>(m, "Mesh")
        .def("lindstrom_turk_by_ratio",
            [] (const rasputin::Mesh& self, double ratio) {
                return self.coarsen(SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
                                    SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                    SMS::LindstromTurk_cost<CGAL::Mesh>());
            }, py::return_value_policy::take_ownership)
        .def("lindstrom_turk_by_size",
            [] (const rasputin::Mesh& self, int max_size) {
                return self.coarsen(SMS::Count_stop_predicate<CGAL::Mesh>(max_size),
                                    SMS::LindstromTurk_placement<CGAL::Mesh>(),
                                    SMS::LindstromTurk_cost<CGAL::Mesh>());
            }, py::return_value_policy::take_ownership)
        .def("copy", &rasputin::Mesh::copy, py::return_value_policy::take_ownership)

        .def_property_readonly("num_vertices", &rasputin::Mesh::num_vertices)
        .def_property_readonly("num_edges", &rasputin::Mesh::num_edges)
        .def_property_readonly("num_faces", &rasputin::Mesh::num_faces)

        .def_property_readonly("points", &rasputin::Mesh::get_points, py::return_value_policy::reference_internal)
        .def_property_readonly("faces", &rasputin::Mesh::get_faces, py::return_value_policy::reference_internal);

    py::class_<std::vector<rasputin::RasterData<float>>, std::unique_ptr<std::vector<rasputin::RasterData<float>>>> (m, "raster_list")
        .def(py::init( [] () {std::vector<rasputin::RasterData<float>> self; return self;}))
        .def("add_raster",
            [] (std::vector<rasputin::RasterData<float>>& self, rasputin::RasterData<float> raster_data) {
                self.push_back(raster_data);
            }, py::keep_alive<1,2>())
        .def("__getitem__",
            [] (std::vector<rasputin::RasterData<float>>& self, int index) {
                return self.at(index);
            }, py::return_value_policy::reference_internal);

    py::bind_vector<std::vector<int>>(m, "int_vector");
    py::bind_vector<std::vector<std::vector<int>>>(m, "shadow_vector");

    m.def("compute_shadow", &rasputin::compute_shadow, "Compute shadows for given sun ray direction.")
     .def("compute_shadows", &rasputin::compute_shadows, "Compute shadows for a series of times and ray directions.")
     .def("surface_normals", &rasputin::surface_normals,
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
