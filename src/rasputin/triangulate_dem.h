//
// Created by Ola Skavhaug on 08/10/2018.
//

#ifndef RASPUTIN_TRIANGULATE_DEM_H_H
#define RASPUTIN_TRIANGULATE_DEM_H_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>


#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

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
    using Point = std::tuple<double, double, double>;
    using Vector = Point;
    using PointList = std::vector<Point>;
    using Face = std::tuple<int, int, int>;
    using FaceList = std::vector<Face>;

    template<typename S, typename P, typename C>
    std::tuple<PointList, FaceList> make_tin(const PointList &pts,
                                             const S& stop,
                                             const P& placement,
                                             const C& cost) {

        CGAL::Delaunay dtin;
        for (const auto p: pts)
            dtin.insert(CGAL::Point(std::get<0>(p),
                                    std::get<1>(p),
                                    std::get<2>(p)));

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
            o_points.emplace_back(std::make_tuple(pt.x(), pt.y(), pt.z()));
            reindex.emplace(v, n++);
        }
        FaceList faces;
        faces.reserve(mesh.num_faces());
        for (auto f: mesh.faces()) {
            std::array<int, 3> fl;
            size_t idx = 0;
            for (auto v: mesh.vertices_around_face(mesh.halfedge(f)))
                fl[idx++] = reindex[v];
            faces.emplace_back(std::make_tuple(fl[0], fl[1], fl[2]));
        }
        return std::make_pair(std::move(o_points), std::move(faces));
    };

    std::vector<int> compute_shadow(const PointList &pts,
                                    const FaceList &faces,
                                    const Vector &sun_direction) {
        std::vector<int> shade;
        CGAL::Mesh mesh;
        std::map<int, CGAL::VertexIndex> index_map;
        std::map<CGAL::face_descriptor, int> face_map;
        size_t i = 0;
        size_t j = 0;
        const CGAL::Vector sun_vec(std::get<0>(sun_direction),
                                   std::get<1>(sun_direction),
                                   std::get<2>(sun_direction));
        for (auto p: pts)
            index_map[i++] = mesh.add_vertex(CGAL::Point(std::get<0>(p), std::get<1>(p), std::get<2>(p)));
        for (auto f: faces)
            face_map[mesh.add_face(index_map[std::get<0>(f)], index_map[std::get<1>(f)], index_map[std::get<2>(f)])] = j++;

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
        std::map<int, CGAL::VertexIndex> index_map;
        std::map<CGAL::face_descriptor, int> face_map;
        size_t i = 0;
        size_t j = 0;
        for (auto p: pts)
            index_map[i++] = mesh.add_vertex(CGAL::Point(std::get<0>(p), std::get<1>(p), std::get<2>(p)));
        for (auto f: faces)
            face_map[mesh.add_face(index_map[std::get<0>(f)], index_map[std::get<1>(f)], index_map[std::get<2>(f)])] = j++;
        CGAL::Tree tree(CGAL::faces(mesh).first, CGAL::faces(mesh).second, mesh);
        for (auto item: sun_rays) {
            std::vector<int> shade;
            shade.reserve(mesh.num_faces());
            auto utc_time = item.first;
            const CGAL::Vector sun_ray(std::get<0>(item.second),
                                       std::get<1>(item.second),
                                       std::get<2>(item.second));

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
    }
}

#endif //RASPUTIN_TRIANGULATE_DEM_H_H