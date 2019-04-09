//
// Created by Ola Skavhaug on 08/10/2018.
//

#ifndef RASPUTIN_TRIANGULATE_DEM_H_H
#define RASPUTIN_TRIANGULATE_DEM_H_H

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Triangulation_face_base_2.h>

#include <armadillo>
#include <cmath>
#include <fstream>
#include <map>
#include <tuple>
#include <numeric>

namespace CGAL {
using K = Exact_predicates_inexact_constructions_kernel;
using Gt = Projection_traits_xy_3<K>;
using Delaunay = Delaunay_triangulation_2<Gt>;
using Point = K::Point_3;
using Vector = K::Vector_3;
using PointList = std::map<Point, int>;
using Mesh = Surface_mesh<Point>;
using VertexIndex = Mesh::Vertex_index;
using PointVertexMap = std::map<Point, VertexIndex>;
using Ray = K::Ray_3;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using Ray_intersection = boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
}

namespace rasputin {
using Point = std::array<double, 3>;
using PointList = std::vector<Point>;
using VectorList = PointList;
using ScalarList = std::vector<double>;
using Vector = Point;
using Face = std::array<int, 3>;
using FaceList = std::vector<Face>;
using VertexIndexMap = std::map<int, CGAL::VertexIndex>;
using FaceDescrMap = std::map<CGAL::face_descriptor, int>;

CGAL::Mesh construct_mesh(const PointList &pts,
                          const FaceList &faces,
                          VertexIndexMap& index_map,
                          FaceDescrMap& face_map){
    CGAL::Mesh mesh;
    index_map.clear();
    face_map.clear();
    size_t i = 0;
    size_t j = 0;
    for (auto p: pts)
        index_map[i++] = mesh.add_vertex(CGAL::Point(p[0], p[1], p[2]));
    for (auto f: faces)
        face_map[mesh.add_face(index_map[f[0]], index_map[f[1]], index_map[f[2]])] = j++;
    return mesh;
};


template<typename S, typename P, typename C>
std::tuple<PointList, FaceList> make_tin(const PointList &pts, const S &stop,
                                         const P &placement, const C &cost) {

    CGAL::Delaunay dtin;
    for (const auto p: pts)
        dtin.insert(CGAL::Point(p[0], p[1], p[2]));

    CGAL::PointVertexMap pvm;
    CGAL::Mesh mesh;
    for (auto v = dtin.finite_vertices_begin(); v != dtin.finite_vertices_end(); ++v)
        pvm.emplace(std::make_pair(v->point(), mesh.add_vertex(v->point())));

    for (auto f = dtin.finite_faces_begin(); f != dtin.finite_faces_end(); ++f)
        mesh.add_face(pvm[f->vertex(0)->point()],
                      pvm[f->vertex(1)->point()],
                      pvm[f->vertex(2)->point()]);
    pvm.clear();
    CGAL::Surface_mesh_simplification::edge_collapse(mesh,
                                                     stop,
                                                     CGAL::parameters::get_cost(cost)
                                                                      .get_placement(placement)
    );
    std::map<CGAL::VertexIndex, int> reindex;
    PointList o_points;
    o_points.reserve(mesh.num_vertices());
    int n = 0;
    for (auto v: mesh.vertices()) {
        const auto pt = mesh.point(v);
        o_points.emplace_back(Vector{pt.x(), pt.y(), pt.z()});
        reindex.emplace(v, n++);
    }
    FaceList faces;
    faces.reserve(mesh.num_faces());
    for (auto f: mesh.faces()) {
        std::array<int, 3> fl;
        size_t idx = 0;
        for (auto v: mesh.vertices_around_face(mesh.halfedge(f)))
            fl[idx++] = reindex[v];
        faces.emplace_back(Face{fl[0], fl[1], fl[2]});
    }
    return std::make_pair(std::move(o_points), std::move(faces));
};

std::vector<int> compute_shadow(const PointList &pts,
                                const FaceList &faces,
                                const Vector &sun_direction) {
    std::vector<int> shade;
    VertexIndexMap index_map;
    FaceDescrMap face_map;
    auto mesh = construct_mesh(pts, faces, index_map, face_map);
    const CGAL::Vector sun_vec(sun_direction[0], sun_direction[1], sun_direction[2]);
    CGAL::Tree tree(CGAL::faces(mesh).first, CGAL::faces(mesh).second, mesh);

    for (auto fd: CGAL::faces(mesh)) {
        auto hd = halfedge(fd, mesh);
        auto p = CGAL::centroid(mesh.point(source(hd, mesh)),
                                mesh.point(target(hd, mesh)),
                                mesh.point(target(next(hd, mesh), mesh)));
        auto v = CGAL::Polygon_mesh_processing::compute_face_normal(fd, mesh);
        if ( v[0]*sun_vec[0] + v[1]*sun_vec[1] + v[2]*sun_vec[2] > 0.0 )
            shade.emplace_back(face_map[fd]);
        else {
            CGAL::Ray sun_ray(p, -sun_vec);
            auto intersection = tree.first_intersection(sun_ray,
                                                        [fd] (const CGAL::face_descriptor &t) { return (t == fd); });
            if (intersection)
                shade.emplace_back(face_map[fd]);
        }
    }
    return shade;
};

std::vector<std::vector<int>> compute_shadows(const PointList &pts,
                                              const FaceList &faces,
                                              const std::vector<std::pair<int, Vector>> & sun_rays) {
    std::vector<std::vector<int>>  result;
    result.reserve(sun_rays.size());
    CGAL::Mesh mesh;
    std::map<size_t, CGAL::VertexIndex> index_map;
    std::map<CGAL::face_descriptor, size_t> face_map;
    size_t i = 0;
    size_t j = 0;
    for (auto p: pts)
        index_map[i++] = mesh.add_vertex(CGAL::Point(p[0], p[1], p[2]));
    for (auto f: faces)
        face_map[mesh.add_face(index_map[f[0]], index_map[f[1]], index_map[f[2]])] = j++;
    CGAL::Tree tree(CGAL::faces(mesh).first, CGAL::faces(mesh).second, mesh);
    for (auto item: sun_rays) {
        std::vector<int> shade;
        shade.reserve(mesh.num_faces());
        auto utc_time = item.first;
        const CGAL::Vector sun_ray(item.second[0],
                                   item.second[1],
                                   item.second[2]);

        for (auto fd: CGAL::faces(mesh)) {
            auto hd = halfedge(fd, mesh);
            auto p = CGAL::centroid(mesh.point(source(hd, mesh)),
                                    mesh.point(target(hd, mesh)),
                                    mesh.point(target(next(hd, mesh), mesh)));
            auto v = CGAL::Polygon_mesh_processing::compute_face_normal(fd, mesh);
            if ( v[0]*sun_ray[0] + v[1]*sun_ray[1] + v[2]*sun_ray[2] >= 0.0 )
                shade.emplace_back(face_map[fd]);
            else {
                CGAL::Ray ray_towards_sun(p, -sun_ray);
                auto intersection = tree.first_intersection(ray_towards_sun,
                                                            [fd] (const CGAL::face_descriptor &t) {
                                                                return (t == fd);
                                                            });
                if (intersection)
                    shade.emplace_back(face_map[fd]);
            }
        }
        result.emplace_back(std::move(shade));
    }
    return result;
};

VectorList orient_tin(const PointList &pts, FaceList &faces) {
    VectorList result;
    result.reserve(faces.size());
    for (auto& face: faces) {
        // Compute ccw normal
        const auto p0 = pts[face[0]];
        const auto p1 = pts[face[1]];
        const auto p2 = pts[face[2]];
        const arma::vec::fixed<3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        const arma::vec::fixed<3> n = arma::cross(v0, v1);
        double c = arma::norm(n);

        // Reverse triangle orientation if it is negatively oriented relative to xy plane
        if (n[2] < 0.0) {
          c *= -1.0;
          std::reverse(face.begin(), face.end());
        }

        // Store normalised and correctly oriented normal vector
        result.push_back(Vector{n[0]/c, n[1]/c, n[2]/c});
    }
    return result;
};

double compute_slope(const Point & normal) {
    return std::atan2(pow(pow(normal[0], 2) + pow(normal[1], 2), 0.5), normal[2]);
}

ScalarList compute_slopes(const VectorList &normals) {
    ScalarList result;
    result.reserve(normals.size());

    for (const auto &n : normals)
        result.push_back(compute_slope(n));

    return result;
};

double compute_aspect(const Point &normal) {
    return std::atan2(normal[0], normal[1]);
}

ScalarList compute_aspects(const VectorList &normals) {
    ScalarList result;
    result.reserve(normals.size());
    for (const auto &n: normals)
        result.emplace_back(compute_aspect(n));
    return result;
};

Vector normal(const Point p0, const Point p1, const Point p2) {
        const arma::vec::fixed<3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        arma::vec::fixed<3> n = arma::cross(v0, v1);
        n /= arma::norm(n);
        return n[2] >= 0.0 ? Vector{n[0], n[1], n[2]} : Vector{-n[0], -n[1], -n[2]};
}

VectorList surface_normals(const PointList &pts, const FaceList &faces) {
    VectorList result;
    result.reserve(faces.size());
    for (const auto face: faces)
        result.emplace_back(normal(pts[face[0]], pts[face[1]], pts[face[2]]));
    return result;
};

template <typename T> void iadd(T& v, const T& o) {
    v[0] += o[0];
    v[1] += o[1];
    v[2] += o[2];
}

VectorList point_normals(const PointList &pts, const FaceList &faces) {
    VectorList result(pts.size(), {0.0, 0.0, 0.0});
    for (auto face: faces) {
        const auto p0 = pts[face[0]];
        const auto p1 = pts[face[1]];
        const auto p2 = pts[face[2]];
        const arma::vec::fixed<3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        arma::vec::fixed<3> n = arma::cross(v0, v1);
        n /= arma::norm(n);
        const auto v = (n[2] >= 0.0) ? Vector{n[0], n[1], n[2]} : Vector{-n[0], -n[1], -n[2]};
        iadd(result[face[0]], v);
        iadd(result[face[1]], v);
        iadd(result[face[2]], v);
    }
    for (auto p: result) {
        const double norm = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        p[0] /= norm;
        p[1] /= norm;
        p[2] /= norm;
    }
    return result;
}

template <typename CB>
std::tuple<FaceList, FaceList> partition(const PointList &pts, const FaceList &faces, CB criterion) {
    FaceList part1, part2;
    for (auto face: faces) {
        if (criterion(pts[face[0]], pts[face[1]], pts[face[2]]))
            part1.emplace_back(face);
        else
            part2.emplace_back(face);
    }
    return std::make_pair(std::move(part1), std::move(part2));
}

std::tuple<FaceList, FaceList> extract_lakes(const PointList &pts, const FaceList &faces) {
   return partition(pts, faces, [] (const Point &p0, const Point &p1, const Point &p2){
       return compute_slope(normal(p0, p1, p2)) < 1.0e-2;
   });
}

std::tuple<FaceList, FaceList> extract_avalanche_expositions(const PointList &pts,
        const FaceList &faces,
        const std::vector<std::array<double, 2>> &exposed_intervals,
        const std::vector<std::array<double, 2>> &height_intervals){
    return partition(pts, faces, [exposed_intervals, height_intervals](const Point &p0, const Point &p1, const Point &p2){
        const auto max_height = std::max(p0[2], std::max(p1[2], p2[2]));
        const auto min_height = std::min(p0[2], std::min(p1[2], p2[2]));
        bool inside = false;
        for (auto height_interval: height_intervals) {
            if ((max_height <= height_interval[1] && max_height >= height_interval[0]) ||
                (min_height <= height_interval[1] && min_height >= height_interval[0]))
                inside = true;
        }
        if (not inside)
            return false;
        const auto cell_normal = normal(p0, p1, p2);
        const auto cell_slope = compute_slope(cell_normal);
        if (cell_slope < 30/180*M_PI)
            return false;
        const auto aspect = compute_aspect(cell_normal);
        for (auto exposition: exposed_intervals) {
            if ((exposition[0] < aspect) && (aspect < exposition[1]))
                return true;
            else if ((exposition[0] > exposition[1]) && ((exposition[0] < aspect) || (aspect < exposition[1])))
                return true;
        }
        return false;
    });
}
}

#endif //RASPUTIN_TRIANGULATE_DEM_H_H
