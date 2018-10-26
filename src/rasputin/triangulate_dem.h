//
// Created by Ola Skavhaug on 08/10/2018.
//

#ifndef RASPUTIN_TINIFY_H_H
#define RASPUTIN_TINIFY_H_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <fstream>
#include <map>
#include <tuple>

namespace CGAL {
    using K = Exact_predicates_inexact_constructions_kernel;
    using Gt = Projection_traits_xy_3<K>;
    using Delaunay = Delaunay_triangulation_2<Gt>;
    using Point = K::Point_3;
    using PointList = std::map<Point, int>;
    using Mesh = Surface_mesh<Point>;
    using VertexIndex = Mesh::Vertex_index;
    using PointVertexMap = std::map<Point, VertexIndex>;
}

namespace rasputin {
    using Point = std::tuple<double, double, double>;
    using PointList = std::vector<Point>;
    using Face = std::tuple<int, int, int>;
    using FaceList = std::vector<Face>;

    template<typename S, typename P, typename C>
    void make_tin(const PointList &pts,
                  const S& stop,
                  const P& placement,
                  const C& cost,
                  PointList &o_points,
                  FaceList &faces) {

        o_points.clear();
        faces.clear();

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
        int n = 0;
        for (auto v: mesh.vertices()) {
            const auto pt = mesh.point(v);
            o_points.emplace_back(std::make_tuple(pt.x(), pt.y(), pt.z()));
            reindex.emplace(v, n++);
        }
        for (auto f: mesh.faces()) {
            std::vector<int> fl;
            for (auto v: mesh.vertices_around_face(mesh.halfedge(f)))
                fl.emplace_back(reindex[v]);
            faces.emplace_back(std::make_tuple(fl[0], fl[1], fl[2]));
        }
    }
}
#endif //RASPUTIN_TINIFY_H_H