#pragma once
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include "types.h"
// #include "raster_data.h"
#include "cgal_raster_data.h"
#include "polygon.h"

namespace rasputin {
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

    Mesh extract_sub_mesh(const std::vector<int> &face_indices) const {
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
        return Mesh(this->construct_mesh(new_points, new_faces, index_map, face_map), proj4_str);
    }

    Mesh construct_mesh(const point3_vector&, const face_vector&, VertexIndexMap&, FaceDescrMap&);

    size_t num_edges() const;
    size_t num_vertices() const;
    size_t num_faces() const;

    template<typename S, typename P, typename C>
    Mesh coarsen(const S& stop, const P& placement, const C& cost) const;
};

template<template<typename> class RasterData, typename Constraints>
struct Delaunay { };
} // namespace traits

namespace rasputin {
template<template<typename> class RasterData, typename Constraints>
struct Delaunay {
    template<typename T, typename P>
    static Constraints interpolate_boundary_points(const RasterData<T>& raster, const P& boundary_polygon) {
        traits::Delaunay<RasterData, Constraints>::interpolate_boundary_points(raster, boundary_polygon);
    }
};


template<typename Pgn>
Mesh make_mesh(const CGAL::PointList &pts,
               const Pgn& inclusion_polygon,
               const CGAL::DelaunayConstraints &constraints,
               const std::string proj4_str) {

    CGAL::ConstrainedDelaunay dtin;
    for (const auto p: pts)
        dtin.insert(p);

    for (auto point_sequence: constraints)
        dtin.insert_constraint(point_sequence.begin(), point_sequence.end(), false);

    CGAL::Mesh mesh;
    CGAL::PointVertexMap pvm;

    auto index = [&pvm, &mesh] (CGAL::Point v) -> CGAL::VertexIndex {
        // Find index of point in mesh, adding it if needed
        auto iter = pvm.find(v);
        if (iter == pvm.end())
            iter = pvm.emplace(std::make_pair(v, mesh.add_vertex(v))).first;
        return iter->second;
    };

    for (auto f = dtin.finite_faces_begin(); f != dtin.finite_faces_end(); ++f) {
        CGAL::Point u = f->vertex(0)->point();
        CGAL::Point v = f->vertex(1)->point();
        CGAL::Point w = f->vertex(2)->point();

        // Add face if midpoint is contained
        CGAL::Point2 face_midpoint(u.x()/3 + v.x()/3 + w.x()/3,
                                   u.y()/3 + v.y()/3 + w.y()/3);
        if (CGAL::point_inside_polygon(face_midpoint, inclusion_polygon)) {
            mesh.add_face(index(u), index(v), index(w));
        }
    }
    pvm.clear();
    return Mesh(mesh, proj4_str);
};


inline Mesh make_mesh(const CGAL::PointList &pts,  const std::string proj4_str) {

    CGAL::ConstrainedDelaunay dtin;
    for (const auto p: pts)
        dtin.insert(p);

    CGAL::PointVertexMap pvm;
    CGAL::Mesh mesh;

    for (auto v = dtin.finite_vertices_begin(); v != dtin.finite_vertices_end(); ++v) {
        CGAL::Point2 pt(v->point().x(), v->point().y());
        pvm.emplace(std::make_pair(v->point(), mesh.add_vertex(v->point())));
    }

    for (auto f = dtin.finite_faces_begin(); f != dtin.finite_faces_end(); ++f) {
        CGAL::Point u = f->vertex(0)->point();
        CGAL::Point v = f->vertex(1)->point();
        CGAL::Point w = f->vertex(2)->point();
        mesh.add_face(pvm[u], pvm[v], pvm[w]);
    }
    pvm.clear();

    return Mesh(mesh, proj4_str);
};

template<typename T, typename Pgn>
Mesh mesh_from_raster(const std::vector<rasputin::RasterData<T>>& raster_list,
                      const Pgn& boundary_polygon,
                      const std::string proj4_str) {
    CGAL::PointList raster_points;
    CGAL::DelaunayConstraints boundary_points;
    for (auto raster : raster_list) {
        CGAL::PointList new_points = raster.raster_points();
        raster_points.insert(raster_points.end(),
                             std::make_move_iterator(new_points.begin()),
                             std::make_move_iterator(new_points.end()));
        CGAL::DelaunayConstraints new_constraints = interpolate_boundary_points(raster, boundary_polygon);
        boundary_points.insert(boundary_points.end(),
                               std::make_move_iterator(new_constraints.begin()),
                               std::make_move_iterator(new_constraints.end()));
    }
    return make_mesh(raster_points, boundary_polygon, boundary_points, proj4_str);
}


template<typename T>
Mesh mesh_from_raster(const std::vector<rasputin::RasterData<T>>& raster_list, const std::string proj4_str) {
    CGAL::PointList raster_points;
    for (auto raster : raster_list) {
        CGAL::PointList new_points = raster.raster_points();
        raster_points.insert(raster_points.end(),
                             std::make_move_iterator(new_points.begin()),
                             std::make_move_iterator(new_points.end()));
    }
    return make_mesh(raster_points, proj4_str);
}


template<typename T, typename Pgn>
Mesh mesh_from_raster(const rasputin::RasterData<T>& raster,
                      const Pgn& boundary_polygon,
                      const std::string proj4_str) {
    CGAL::PointList raster_points = raster.raster_points();
    CGAL::DelaunayConstraints boundary_points = interpolate_boundary_points(raster, boundary_polygon);

    return make_mesh(raster_points, boundary_polygon, boundary_points, proj4_str);
}


template<typename T>
Mesh mesh_from_raster(const rasputin::RasterData<T>& raster, const std::string proj4_str) {
    return make_mesh(raster.raster_points(), proj4_str);
}

} // namespace rasputin
