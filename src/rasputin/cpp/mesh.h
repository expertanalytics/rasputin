#pragma once
#include <pybind11/pybind11.h>

#include "types.h"
#include "polygon.h"

namespace rasputin {
namespace py = pybind11;

namespace traits {
template<typename Mesh>
struct MeshBase {
    private:
    point3_vector points;
    face_vector faces;

    void set_points_faces();

    public:
    const Mesh mesh;
    const std::string proj4_str;

    MeshBase(Mesh mesh, const std::string proj4_str) : mesh(mesh), proj4_str(proj4_str) {
        this->set_points_faces();
    }

    MeshBase copy() const {
        // deep copy
        return MeshBase(Mesh(this->mesh), this->proj4_str);
    }

    const point3_vector& get_points() const {
        return this->points;
    }

    const face_vector& get_faces() const {
        return this->faces;
    }

    virtual MeshBase extract_sub_mesh(const std::vector<int> &face_indices) const {
        std::map<int, int> remap;
        point3_vector new_points;
        face_vector new_faces;
        int counter = 0;
        for (auto face_idx: face_indices) {
            std::array<int, 3> new_face;
            int i = 0;
            for (auto idx: this->faces[face_idx]) {
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
        return MeshBase(this->construct_mesh(new_points, new_faces, index_map, face_map), proj4_str);
    }

    static Mesh construct_mesh(const point3_vector&, const face_vector&, VertexIndexMap&, FaceDescrMap&);

    size_t num_edges() const;
    size_t num_vertices() const;
    size_t num_faces() const;

    template<typename S, typename P, typename C>
    Mesh coarsen(const S& stop, const P& placement, const C& cost) const;
};

template<typename Mesh, typename M>
concept DerivedFromMesh = std::derived_from<Mesh, traits::MeshBase<M>>;

template<
    template<typename> class RasterData,
    typename PointList,
    typename Mesh,
    typename Constraints
>
struct Triangulation {
    static_assert("Triangulation trait not implemented");
};
} // namespace traits

template<
    template<typename> class RasterData,
    typename PointList,
    typename Mesh,
    typename M,
    typename Constraints
>
requires std::derived_from<Mesh, traits::MeshBase<M>>
struct Triangulation {
    using D = traits::Triangulation<RasterData, PointList, Mesh, Constraints>;

    template<typename T, typename P>
    static Constraints interpolate_boundary_points(const RasterData<T>& raster, const P& boundary_polygon) {
        D::template interpolate_boundary_points<T, P>(raster, boundary_polygon);
    }

    template<typename Pgn>
    static Mesh make_mesh(
        const PointList &pts, const Pgn& inclusion_polygon, const Constraints &constraints, const std::string proj4_str
    ) {
        return D::template make_mesh<Pgn>(pts, inclusion_polygon, constraints, proj4_str);
    }

    static Mesh make_mesh(const PointList &pts,  const std::string proj4_str) {
        return D::make_mesh(pts, proj4_str);
    }

    template<typename T, typename Pgn>
    static Mesh mesh_from_raster(
        const std::vector<RasterData<T>>& raster_list, const Pgn& boundary_polygon, const std::string proj4_str
    ) {
        return D::template mesh_from_raster<T, Pgn>(raster_list, boundary_polygon, proj4_str);
    }

    template<typename T>
    static Mesh mesh_from_raster(const std::vector<RasterData<T>>& raster_list, const std::string proj4_str) {
        return D::template mesh_from_raster<T>(raster_list, proj4_str);
    }

    template<typename T, typename Pgn>
    static Mesh mesh_from_raster(
        const RasterData<T>& raster, const Pgn& boundary_polygon, const std::string proj4_str
    ) {
        return D::template mesh_from_raster<T, Pgn>(raster, boundary_polygon, proj4_str);
    }

    template<typename T>
    static Mesh mesh_from_raster(const RasterData<T>& raster, const std::string proj4_str) {
        return D::template mesh_from_raster<T>(raster, proj4_str);
    }
};

template<
    template<typename> class R,
    typename T,
    typename P,
    typename PointList,
    typename Mesh,
    typename M,
    typename Constraints
>
static void bind_make_mesh(py::module &m) {
    m.def("make_mesh",
        [] (const R<T>& raster_data, const P polygon, const std::string proj4_str) {
            return Triangulation<R, PointList, Mesh, M, Constraints>::template mesh_from_raster<T, P>(raster_data, polygon, proj4_str);
        }, py::return_value_policy::take_ownership)
    .def("make_mesh",
        [] (const R<T>& raster_data, const std::string proj4_str) {
            return Triangulation<R, PointList, Mesh, M, Constraints>::template mesh_from_raster<T>(raster_data, proj4_str);
        }, py::return_value_policy::take_ownership);
    m.def("make_mesh",
        [] (const std::vector<R<T>>& raster_data, const P polygon, const std::string proj4_str) {
            return Triangulation<R, PointList, Mesh, M, Constraints>::template mesh_from_raster<T, P>(raster_data, polygon, proj4_str);
        }, py::return_value_policy::take_ownership)
    .def("make_mesh",
        [] (const std::vector<R<T>>& raster_data, const std::string proj4_str) {
            return Triangulation<R, PointList, Mesh, M, Constraints>::template mesh_from_raster<T>(raster_data, proj4_str);
        }, py::return_value_policy::take_ownership);
}
} // namespace rasputin
