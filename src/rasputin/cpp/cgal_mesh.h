#pragma once

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_placement.h>

#include "cgal_polygon.h"
#include "cgal_raster_data.h"
#include "mesh.h"

namespace SMS = CGAL::Surface_mesh_simplification;

namespace rasputin {
using VertexIndex = CGAL::VertexIndex;
using face_descriptor = CGAL::face_descriptor;
using Constraints = CGAL::DelaunayConstraints;

using Base = traits::MeshBase<
    CGAL::Mesh, VertexIndex, face_descriptor, Constraints, PointList, Polygon, SimplePolygon, MultiPolygon
>;

struct Mesh: Base {
    Mesh(CGAL::Mesh mesh, const std::string proj4_str): Base(mesh, proj4_str) { }

    static constexpr char coarsen_ratio_name[] = "lindstrom_turk_by_ratio";
    static constexpr char coarsen_size_name[] = "lindstrom_turk_by_size";
    static constexpr char coarsen_ratio_description[] =
        "Simplify the mesh.\n"
        "\nThe LindstromTurk cost and placement strategy is used, "
        "and simplification process stops when the number of undirected edges drops below the size threshold.";
    static constexpr char coarsen_size_description[] =
        "Simplify the mesh.\n"
        "\nThe LindstromTurk cost and placement strategy is used, "
        "and simplification process stops when the number of undirected edges drops below the size threshold.";

    Mesh copy() const {
        // Call CGAL::Mesh copy constructor to do a deep copy
        return Mesh(CGAL::Mesh(this->mesh), this->proj4_str);
    }

    template<typename Stop>
    requires std::same_as<Stop, SMS::Count_ratio_stop_predicate<CGAL::Mesh>>
    Mesh _coarsen(
        const Stop& stop,
        const SMS::LindstromTurk_placement<CGAL::Mesh>& placement,
        const SMS::LindstromTurk_cost<CGAL::Mesh>& cost
    ) const {
        CGAL::Mesh new_cgal_mesh = CGAL::Mesh(this->mesh);
        CGAL::Surface_mesh_simplification::edge_collapse(
            new_cgal_mesh, stop, CGAL::parameters::get_cost(cost).get_placement(placement)
        );
        return Mesh(new_cgal_mesh, this->proj4_str);
    }

    template<traits::PredRatioType PredRatio>
    Mesh coarsen(PredRatio ratio) const {
        return this->_coarsen<SMS::Count_ratio_stop_predicate<CGAL::Mesh>>(
            SMS::Count_ratio_stop_predicate<CGAL::Mesh>(ratio),
            SMS::LindstromTurk_placement<CGAL::Mesh>(),
            SMS::LindstromTurk_cost<CGAL::Mesh>()
        );
    }

    Mesh extract_sub_mesh(const std::vector<int> &face_indices) const {
        return Mesh(this->construct_sub_mesh(face_indices), this->proj4_str);
    }

    size_t num_edges() const {
        return this->mesh.number_of_edges();
    }
    size_t num_vertices() const {
        return this->mesh.number_of_vertices();
    }
    size_t num_faces() const {
        return this->mesh.number_of_faces();
    }
};

namespace traits {
template<>
struct MeshInit<CGAL::Mesh, VertexIndex, face_descriptor> {
    using Mesh_t = CGAL::Mesh;
    using VertexIndex_t = VertexIndex;
    using face_descriptor_t = face_descriptor;
    using VertexIndexMap = std::map<int, VertexIndex>;
    using FaceDescrMap = std::map<face_descriptor, int>;

    static void set_points_faces(Mesh_t& mesh, point3_vector& points, face_vector& faces) {
        points.clear();
        faces.clear();

        points.reserve(mesh.num_vertices());
        faces.reserve(mesh.num_faces());

        int n = 0;
        std::map<CGAL::VertexIndex, int> reindex;
        for (auto f: mesh.faces()) {
            std::array<int, 3> fl;
            size_t idx = 0;
            for (auto v: mesh.vertices_around_face(mesh.halfedge(f))) {
                if (reindex.count(v) == 0) {
                    reindex.emplace(v, n++);
                    const auto pt = mesh.point(v);
                    points.emplace_back(point3{pt.x(), pt.y(), pt.z()});
                }
                fl[idx++] = reindex[v];
            }
            faces.emplace_back(face{fl[0], fl[1], fl[2]});
        }
    }

    static Mesh_t construct_mesh(
        const point3_vector &pts, const face_vector &faces, VertexIndexMap& index_map, FaceDescrMap& face_map
    ) {
        Mesh_t mesh;
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
};

template<traits::RasterType RasterData>
struct Triangulation<Mesh, RasterData> {
    using PolygonUtils_t = PolygonUtils<typename RasterData::Polygon_t, typename RasterData::SimplePolygon_t>;
    using T = traits::Triangulation<Mesh, RasterData>;
    using RasterData_t = RasterData;
    using Mesh_t = Mesh;
    using VertexIndex_t = typename Mesh_t::VertexIndex_t;
    using face_descriptor_t = typename Mesh_t::face_descriptor_t;
    using Constraints_t = typename Mesh_t::Constraints_t;

    using VertexIndexMap = std::map<int, VertexIndex>;
    using FaceDescrMap = std::map<face_descriptor, int>;

    template<rasputin::PolygonType<RasterData> P>
    static Constraints interpolate_boundary_points(const RasterData& raster, const P& boundary_polygon) {
        // First we need to determine intersection points between the raster domain and the polygon
        CGAL::MultiPolygon intersection_polygon = raster.compute_intersection(boundary_polygon);

        // Iterate over edges of the intersection polygon and interpolate points
        // TODO: Handle holes. Only needed if the intersecting polygon has holes.
        CGAL::DelaunayConstraints interpolated_points;
        for (const auto& part : PolygonUtils_t::extract_boundaries(intersection_polygon)) {
            for (auto e = part.edges_begin(); e != part.edges_end(); ++e) {
                CGAL::Point2 first_vertex = e->vertex(0);
                CGAL::Point2 second_vertex = e->vertex(1);

                // We can skip edges that are aligend with the raster boundary
                // TODO: Consider checking using CGAL and exact arithmentic
                bool edge_is_aligned = not raster.contains((first_vertex.x() + second_vertex.x())/2,
                                                           (first_vertex.y() + second_vertex.y())/2);
                if (edge_is_aligned) {
                    continue;
                }

                // Sample with the approximately the same resolution as the raster data along the boundary edges
                double edge_len_x = second_vertex.x() - first_vertex.x();
                double edge_len_y = second_vertex.y() - first_vertex.y();
                std::size_t num_subedges = static_cast<int>(std::max<double>(std::fabs(edge_len_x/raster.delta_x),
                                                                             std::fabs(edge_len_y/raster.delta_y)));
                num_subedges = std::max<int>(1, num_subedges);

                double edge_dx = edge_len_x / num_subedges; // signed distance
                double edge_dy = edge_len_y / num_subedges;

                // Iterate over edge samples
                std::vector<CGAL::Point> interpolated_points_on_edge;
                for (std::size_t k=0; k < num_subedges + 1; ++k) {
                    double x = first_vertex.x() + k * edge_dx;
                    double y = first_vertex.y() + k * edge_dy;
                    double z = raster.get_interpolated_value_at_point(x, y);
                    interpolated_points_on_edge.emplace_back(x, y, z);
                }
                interpolated_points.push_back(std::move(interpolated_points_on_edge));
            }
        }
        return interpolated_points;
    }

    template<rasputin::PolygonType<RasterData> Pgn>
    static Mesh_t make_mesh(
        const typename RasterData::PointList_t &pts,
        const Pgn& inclusion_polygon,
        const Constraints &constraints,
        const std::string proj4_str
    ) {
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
            if (PolygonUtils_t::point_inside_polygon(face_midpoint, inclusion_polygon)) {
                mesh.add_face(index(u), index(v), index(w));
            }
        }
        pvm.clear();
        return Mesh_t(mesh, proj4_str);
    }

    static Mesh_t make_mesh(const typename RasterData::PointList_t &pts,  const std::string proj4_str) {
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

        return Mesh_t(mesh, proj4_str);
    }

    template<rasputin::PolygonType<RasterData> Pgn>
    static Mesh_t mesh_from_raster(
        const std::vector<RasterData>& raster_list, const Pgn& boundary_polygon, const std::string proj4_str
    ) {
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

    static Mesh_t mesh_from_raster(const std::vector<RasterData>& raster_list, const std::string proj4_str) {
        CGAL::PointList raster_points;
        for (auto raster : raster_list) {
            CGAL::PointList new_points = raster.raster_points();
            raster_points.insert(raster_points.end(),
                                 std::make_move_iterator(new_points.begin()),
                                 std::make_move_iterator(new_points.end()));
        }
        return make_mesh(raster_points, proj4_str);
    }

    template<rasputin::PolygonType<RasterData> Pgn>
    static Mesh_t mesh_from_raster(const RasterData& raster, const Pgn& boundary_polygon, const std::string proj4_str) {
        CGAL::PointList raster_points = raster.raster_points();
        CGAL::DelaunayConstraints boundary_points = interpolate_boundary_points(raster, boundary_polygon);

        return make_mesh(raster_points, boundary_polygon, boundary_points, proj4_str);
    }

    static Mesh_t mesh_from_raster(const RasterData& raster, const std::string proj4_str) {
        return make_mesh(raster.raster_points(), proj4_str);
    }
};
} // namespace traits
} // namespace rasputin
