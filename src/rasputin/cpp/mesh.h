#pragma once

#include <concepts>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "types.h"
#include "raster_data.h"

namespace py = pybind11;

namespace rasputin {
namespace traits {
template<typename Mesh, typename VertexIndex, typename face_descriptor>
struct MeshInit {
    static_assert("Trait MeshInit not implemented!");
};

template<typename Mesh, typename RasterData>
struct Triangulation {
    static_assert("Trait Trianguation not implemented!");
};

template<typename PredRatio>
concept PredRatioType =
    std::same_as<PredRatio, double> || (std::unsigned_integral<PredRatio> && !std::same_as<PredRatio, char>);

template<
    typename Mesh,
    typename VertexIndex,
    typename face_descriptor,
    typename Constraints,
    typename PointList,
    typename Polygon,
    typename SimplePolygon,
    typename MultiPolygon
>
class MeshBase {
    public:
    using Mesh_t = Mesh;
    using VertexIndex_t = VertexIndex;
    using face_descriptor_t = face_descriptor;
    using PointList_t = PointList;
    using Polygon_t = Polygon;
    using SimplePolygon_t = SimplePolygon;
    using MultiPolygon_t = MultiPolygon;
    using Constraints_t = Constraints;

    const Mesh mesh;
    const std::string proj4_str;
    point3_vector points;
    face_vector faces;

    MeshBase(Mesh mesh, const std::string proj4_str): mesh(mesh), proj4_str(proj4_str) {
        MeshInit<Mesh_t, VertexIndex_t, face_descriptor_t>::set_points_faces(mesh, points, faces);
    }

    template<PredRatioType PredRatio>
    Mesh coarsen(PredRatio ratio) const;

    size_t num_edges() const;
    size_t num_vertices() const;
    size_t num_faces() const;

    Mesh_t construct_sub_mesh(const std::vector<int>& face_indices) const {
        std::map<int, int> remap;
        point3_vector new_points;
        face_vector new_faces;
        int counter = 0;
        for (auto face_idx: face_indices) {
            std::array<int, 3> new_face;
            int i = 0;
            for (auto idx: faces[face_idx]) {
                if (remap.count(idx) == 0) {
                    remap[idx] = counter++;
                    new_points.emplace_back(points[idx]);
                }
                new_face[i++] = remap[idx];
            }
            new_faces.emplace_back(new_face);
        }

        VertexIndexMap index_map;
        FaceDescrMap face_map;

        return MeshInit<Mesh_t, VertexIndex_t, face_descriptor_t>::construct_mesh(
            new_points, new_faces, index_map, face_map
        );
    }

    const point3_vector& get_points() const {
        return points;
    }

    const face_vector& get_faces() const {
        return faces;
    }
};

template<typename Mesh>
concept MeshInitType =
    requires {
        typename Mesh::Mesh_t;
        typename Mesh::VertexIndex_t;
        typename Mesh::face_descriptor_t;
        typename Mesh::VertexIndexMap;
        typename Mesh::FaceDescrMap;
    }
    && requires {
        { Mesh::VertexIndexMap } -> std::same_as<std::map<int, typename Mesh::VertexIndex_t>>;
        { Mesh::FaceDescrMap } -> std::same_as<std::map<typename Mesh::face_descriptor_t, int>>;
    };

template<typename Mesh>
concept MeshType =
    requires(Mesh m) {
        { m.copy() } -> std::same_as<Mesh>;
        { m.coarsen_ratio_name[0] } -> std::same_as<const char&>;
        { m.coarsen_ratio_description[0] } -> std::same_as<const char&>;
        { m.coarsen_size_name[0] } -> std::same_as<const char&>;
        { m.coarsen_size_description[0] } -> std::same_as<const char&>;
    }
    && requires {
        typename Mesh::Mesh_t;
        typename Mesh::VertexIndex_t;
        typename Mesh::face_descriptor_t;
        typename Mesh::Constraints_t;
        typename Mesh::PointList_t;
        typename Mesh::Polygon_t;
        typename Mesh::SimplePolygon_t;
        typename Mesh::MultiPolygon_t;
    }
    && requires(Mesh m, const std::vector<int>& face_indices) {
        { m.extract_sub_mesh(face_indices) } -> std::same_as<Mesh>;
    }
    && std::derived_from<
        Mesh,
        MeshBase<
            typename Mesh::Mesh_t,
            typename Mesh::VertexIndex_t,
            typename Mesh::face_descriptor_t,
            typename Mesh::Constraints_t,
            typename Mesh::PointList_t,
            typename Mesh::Polygon_t,
            typename Mesh::SimplePolygon_t,
            typename Mesh::MultiPolygon_t
        >
    >;

} // namespace traits

template<typename Mesh, typename VertexIndex, typename face_descriptor>
struct MeshInit {
    using M = traits::MeshInit<Mesh, VertexIndex, face_descriptor>;
    using Mesh_t = Mesh;
    using VertexIndex_t = VertexIndex;
    using face_descriptor_t = face_descriptor;
    using VertexIndexMap = std::map<int, VertexIndex_t>;
    using FaceDescrMap = std::map<face_descriptor_t, int>;

    static void set_points_faces(Mesh_t const& mesh, point3_vector& points, face_vector& faces) {
        M::set_points_faces(mesh);
    }

    static Mesh_t construct_mesh (
        const point3_vector &pts, const face_vector &faces, VertexIndexMap& index_map, FaceDescrMap& face_map
    ) {
        return M::construct_mesh(pts, faces, index_map, face_map);
    }
};

template<traits::MeshType Mesh, traits::RasterType RasterData>
struct Triangulation {
    using RasterData_t = RasterData;
    using Mesh_t = Mesh;
    using VertexIndex_t = typename Mesh::VertexIndex_t;
    using face_descriptor_t = typename Mesh::face_descriptor_t;
    using Constraints_t = typename Mesh::Constraints_t;
    using T = traits::Triangulation<Mesh, RasterData>;
    using VertexIndexMap = std::map<int, VertexIndex_t>;
    using FaceDescrMap = std::map<face_descriptor_t, int>;

    template<PolygonType<RasterData_t> P>
    static Constraints_t interpolate_boundary_points(const RasterData_t& raster, const P& boundary_polygon) {
        return T::template interpolate_boundary_points<RasterData_t, P>(raster, boundary_polygon);
    }

    template<PolygonType<RasterData_t> Pgn>
    static Mesh make_mesh(
        const typename RasterData_t::PointList_t &pts,
        const Pgn& inclusion_polygon,
        const Constraints_t &constraints,
        const std::string proj4_str
    ) {
        return T::template make_mesh<Pgn>(pts, inclusion_polygon, constraints, proj4_str);
    }

    static Mesh make_mesh(const typename RasterData_t::PointList_t &pts,  const std::string proj4_str) {
        return T::template make_mesh(pts, proj4_str);
    }

    template<PolygonType<RasterData_t> Pgn>
    static Mesh mesh_from_raster(
        const std::vector<RasterData_t>& raster_list, const Pgn& boundary_polygon, const std::string proj4_str
    ) {
        return T::template mesh_from_raster<Pgn>(raster_list, boundary_polygon, proj4_str);
    }

    static Mesh mesh_from_raster(const std::vector<RasterData_t>& raster_list, const std::string proj4_str) {
        return T::mesh_from_raster(raster_list, proj4_str);
    }

    template<PolygonType<RasterData_t> Pgn>
    static Mesh mesh_from_raster(const RasterData_t& raster, const Pgn& boundary_polygon, const std::string proj4_str) {
        return T::template mesh_from_raster<Pgn>(raster, proj4_str);
    }

    static Mesh mesh_from_raster(const RasterData_t& raster, const std::string proj4_str) {
        return T::mesh_from_raster(raster, proj4_str);
    }
};

template<typename T>
concept TriangulationType =
    requires {
        typename T::RasterData_t;
        typename T::Mesh_t;
        typename T::VertexIndex_t;
        typename T::face_descriptor_t;
        typename T::Constraints_t;
    }
    && traits::MeshType<typename T::Mesh_t>;

template<TriangulationType T, PolygonType<typename T::RasterData_t> P>
void bind_make_mesh(py::module &m) {
    using R = typename T::RasterData_t;
    using Mesh_t = typename T::Mesh_t::Mesh_t;
    using M = MeshInit<Mesh_t, typename T::VertexIndex_t, typename T::face_descriptor_t>;

    m.def("make_mesh",
        [] (const R& raster_data, const P polygon, const std::string proj4_str) {
            return T::mesh_from_raster(raster_data, polygon, proj4_str);
        }, py::return_value_policy::take_ownership)
    .def("make_mesh",
        [] (const R& raster_data, const std::string proj4_str) {
            return T::mesh_from_raster(raster_data, proj4_str);
        }, py::return_value_policy::take_ownership)
    .def("make_mesh",
        [] (const std::vector<R>& raster_data, const P polygon, const std::string proj4_str) {
            return T::mesh_from_raster(raster_data, polygon, proj4_str);
        }, py::return_value_policy::take_ownership)
    .def("make_mesh",
        [] (const std::vector<R>& raster_data, const std::string proj4_str) {
            return T::mesh_from_raster(raster_data, proj4_str);
        }, py::return_value_policy::take_ownership)
    .def("construct_mesh",
        [] (const rasputin::point3_vector& points, const rasputin::face_vector & faces, const std::string proj4_str) {
            typename T::VertexIndexMap index_map;
            typename T::FaceDescrMap face_map;
            return typename T::Mesh_t(M::construct_mesh(points, faces, index_map, face_map), proj4_str);
        }, py::return_value_policy::take_ownership);
}

template<TriangulationType T>
void bind_make_mesh(py::module& m) {
    using R = typename T::RasterData_t;
    rasputin::bind_make_mesh<T, typename R::SimplePolygon_t>(m);
    rasputin::bind_make_mesh<T, typename R::Polygon_t>(m);
}

template<traits::MeshType Mesh, traits::RasterType ...RasterData>
void bind_mesh(py::module& m) {
    py::class_<Mesh, std::unique_ptr<Mesh>>(m, "Mesh")
        .def(Mesh::coarsen_ratio_name,
            [] (const Mesh& self, double ratio) { return self.template coarsen<double>(ratio); },
            py::return_value_policy::take_ownership,
            Mesh::coarsen_ratio_description
        )
        .def(Mesh::coarsen_size_name,
            [] (const Mesh& self, std::size_t max_size) { return self.template coarsen<std::size_t>(max_size); },
            py::return_value_policy::take_ownership,
            Mesh::coarsen_size_description
        )
        .def("copy", &Mesh::copy, py::return_value_policy::take_ownership)
        .def("extract_sub_mesh", &Mesh::extract_sub_mesh, py::return_value_policy::take_ownership)
        .def_property_readonly("num_vertices", &Mesh::num_vertices)
        .def_property_readonly("num_edges", &Mesh::num_edges)
        .def_property_readonly("num_faces", &Mesh::num_faces)
        .def_property_readonly("points", &Mesh::get_points, py::return_value_policy::reference_internal)
        .def_property_readonly("faces", &Mesh::get_faces, py::return_value_policy::reference_internal);

    (
    bind_make_mesh<traits::Triangulation<Mesh, RasterData>>(m) ,...);
}
} // namespace rasputin
