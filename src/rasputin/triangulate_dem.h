//
// Created by Ola Skavhaug on 08/10/2018.
//

#pragma once

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <boost/geometry.hpp>
#include <boost/geometry/srs/epsg.hpp>
#include <boost/geometry/srs/projection.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/srs/transformation.hpp>

#include <sstream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <map>
#include <tuple>
#include <numeric>
#include <pybind11/numpy.h>
#include <cstdint>
#include "solar_position.h"

namespace CGAL {
using K = Exact_predicates_inexact_constructions_kernel;
using Gt = Projection_traits_xy_3<K>;
using Delaunay = Delaunay_triangulation_2<Gt>;
using ConstrainedDelaunay = Constrained_Delaunay_triangulation_2<Gt, CGAL::Default, CGAL::Exact_predicates_tag>;

using Point = Gt::Point_2;
using Point3 = K::Point_3;
using Point2 = K::Point_2;

using Vector = K::Vector_3;
using PointList = std::vector<Point>;
using Mesh = Surface_mesh<Point>;
using VertexIndex = Mesh::Vertex_index;
using FaceIndex = Mesh::Face_index;
using PointVertexMap = std::map<Point, VertexIndex>;
using Ray = K::Ray_3;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using Ray_intersection = boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

using PointSequence = std::vector<Point>;
using DelaunayConstraints = std::vector<PointSequence>;

using SimplePolygon = Polygon_2<K>;
using Polygon = Polygon_with_holes_2<K>;
using MultiPolygon = std::vector<Polygon>;

bool point_inside_polygon(const Point2 &x, const SimplePolygon &polygon) {
        return polygon.has_on_bounded_side(x);
}

bool point_inside_polygon(const Point2 &x, const Polygon &polygon) {
    if (not polygon.outer_boundary().has_on_bounded_side(x))
        return false;
    if (not polygon.has_holes())
        return true;
    for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
        if (it->has_on_bounded_side(x))
            return false;
    return true;
}

bool point_inside_polygon(const Point2 &x, const MultiPolygon &polygon) {
    for (const auto& part: polygon)
        if (point_inside_polygon(x, part))
            return true;
    return false;
}

std::vector<SimplePolygon> extract_boundaries(const SimplePolygon &polygon) {
    std::vector<SimplePolygon> ret;
    ret.emplace_back(polygon);

    return ret;
}
std::vector<SimplePolygon> extract_boundaries(const Polygon &polygon) {
    std::vector<SimplePolygon> ret;
    ret.emplace_back(polygon.outer_boundary());
    for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
        ret.emplace_back(*it);

    return ret;
}

std::vector<SimplePolygon> extract_boundaries(const MultiPolygon &polygon) {
    std::vector<SimplePolygon> ret;
    for (const auto& part: polygon) {
        ret.emplace_back(part.outer_boundary());
        for (auto it = part.holes_begin(); it != part.holes_end(); ++it)
            ret.emplace_back(*it);
    }

    return ret;
}
}

namespace boost::geometry::traits {

    template<> struct tag<CGAL::Point2>
    { typedef point_tag type; };

    template<> struct coordinate_type<CGAL::Point2>
    { typedef double type; };

    template<> struct coordinate_system<CGAL::Point2>
    {typedef cs::cartesian type; };

    template<> struct dimension<CGAL::Point2> : boost::mpl::int_<2> {};

    template <>
    struct access<CGAL::Point2, 0>
    {
        static double get(CGAL::Point2 const& p) { return p.x(); }
        static void set(CGAL::Point2& p, double const& value) { p = CGAL::Point2{value, p.y()}; }
    };

    template <>
    struct access<CGAL::Point2, 1>
    {
        static double get(CGAL::Point2 const& p) { return p.y(); }
        static void set(CGAL::Point2& p, double const& value) { p = CGAL::Point2{p.x(), value}; }
    };

    template<> struct tag<CGAL::Point3>
    { typedef point_tag type; };

    template<> struct coordinate_type<CGAL::Point3>
    { typedef double type; };

    template<> struct coordinate_system<CGAL::Point3>
    {typedef cs::cartesian type; };

    template<> struct dimension<CGAL::Point3> : boost::mpl::int_<3> {};

    template <>
    struct access<CGAL::Point3, 0>
    {
        static double get(CGAL::Point3 const& p) { return p.x(); }
        static void set(CGAL::Point3& p, double const& value) { p = CGAL::Point3{value, p.y(), p.z()}; }
    };

    template <>
    struct access<CGAL::Point3, 1>
    {
        static double get(CGAL::Point3 const& p) { return p.y(); }
        static void set(CGAL::Point3& p, double const& value) { p = CGAL::Point3{p.x(), value, p.z()}; }
    };

    template <>
    struct access<CGAL::Point3, 2>
    {
        static double get(CGAL::Point3 const& p) { return p.z(); }
        static void set(CGAL::Point3& p, double const& value) { p = CGAL::Point3{p.x(), p.y(), value}; }
    };
}


namespace rasputin {
using point3 = std::array<double, 3>;
using point3_vector = std::vector<std::array<double, 3>>;
using point2 = std::array<double, 2>;
using point2_vector= std::vector<point2>;
using face = std::array<int, 3>;
using index = std::array<unsigned int, 2>;
using face_vector = std::vector<face>;
using index_vector = std::vector<index>;
using double_vector = std::vector<double>;
using uint8_vector = std::vector<std::uint8_t>;

// Clean up below
using VertexIndexMap = std::map<int, CGAL::VertexIndex>;
using FaceDescrMap = std::map<CGAL::face_descriptor, int>;

CGAL::Mesh construct_mesh(const point3_vector &pts,
                          const face_vector &faces,
                          VertexIndexMap& index_map,
                          FaceDescrMap& face_map);
struct Mesh {
    const CGAL::Mesh cgal_mesh;
    const unsigned int epsg_id;
    Mesh(CGAL::Mesh cgal_mesh, const unsigned int epsg_id) : cgal_mesh(cgal_mesh), epsg_id(epsg_id) {set_points_faces();}

    template<typename S, typename P, typename C>
    Mesh coarsen(const S& stop, const P& placement, const C& cost) const {
        CGAL::Mesh new_cgal_mesh = CGAL::Mesh(this->cgal_mesh);
        CGAL::Surface_mesh_simplification::edge_collapse(new_cgal_mesh,
                                                         stop,
                                                         CGAL::parameters::get_cost(cost)
                                                                          .get_placement(placement));
        return Mesh(new_cgal_mesh, epsg_id);
    }

    Mesh copy() const {
        // Call CGAL::Mesh copy constructor to do a deep copy
        return Mesh(CGAL::Mesh(cgal_mesh), epsg_id);
    }

    const point3_vector& get_points() const {
        return points;
    }

    const face_vector& get_faces() const {
        return faces;
    }

    size_t num_edges() const {return cgal_mesh.number_of_edges();}
    size_t num_vertices() const {return cgal_mesh.number_of_vertices();}
    size_t num_faces() const {return cgal_mesh.number_of_faces();}

    Mesh extract_sub_mesh(const std::vector<int> &face_indices) const {
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
        return Mesh(construct_mesh(new_points, new_faces, index_map, face_map), epsg_id);
    }

    private:
    point3_vector points;
    face_vector faces;

    void set_points_faces() {
        points.clear();
        faces.clear();

        points.reserve(cgal_mesh.num_vertices());
        faces.reserve(cgal_mesh.num_faces());

        int n = 0;
        std::map<CGAL::VertexIndex, int> reindex;
        for (auto f: cgal_mesh.faces()) {
            std::array<int, 3> fl;
            size_t idx = 0;
            for (auto v: cgal_mesh.vertices_around_face(cgal_mesh.halfedge(f))) {
                if (reindex.count(v) == 0) {
                    reindex.emplace(v, n++);
                    const auto pt = cgal_mesh.point(v);
                    points.emplace_back(point3{pt.x(), pt.y(), pt.z()});
                }
                fl[idx++] = reindex[v];
            }
            faces.emplace_back(face{fl[0], fl[1], fl[2]});
        }
    }
};

CGAL::Mesh construct_mesh(const point3_vector &pts,
                          const face_vector &faces,
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


template<typename FT>
struct RasterData {
    RasterData(double x_min, double y_max, double delta_x, double delta_y,
               std::size_t num_points_x, std::size_t num_points_y,
               FT* data)
        : x_min(x_min), delta_x(delta_x), num_points_x(num_points_x),
          y_max(y_max), delta_y(delta_y), num_points_y(num_points_y),
          data(data) {}

    double x_min;
    double delta_x;
    std::size_t num_points_x;

    double y_max;
    double delta_y;
    std::size_t num_points_y;

    FT* data;

    double get_x_max() const {return x_min + (num_points_x - 1)*delta_x; }
    double get_y_min() const {return y_max - (num_points_y - 1)*delta_y; }

    CGAL::PointList raster_points() const {
        CGAL::PointList points;
        points.reserve(num_points_x * num_points_y);

        for (std::size_t i = 0; i < num_points_x; ++i)
            for (std::size_t j = 0; j < num_points_y; ++j)
                points.emplace_back(x_min + i*delta_x, y_max - j*delta_y, data[i*num_points_y + j]);
        return points;
    }

    // For every point inside the raster rectangle we identify indices (i, j) of the upper-left vertex of the cell containing the point
    std::pair<int, int> get_indices(double x, double y) const {
        int i = std::min<int>(std::max<int>(static_cast<int>((x-x_min) /delta_x), 0), num_points_x - 1);
        int j = std::min<int>(std::max<int>(static_cast<int>((y_max-y) /delta_y), 0), num_points_y - 1);

        return std::make_pair(i,j);
    }

    // Interpolate data using using a bilinear interpolation rule on each cell
    FT get_interpolated_value_at_point(double x, double y) const {
        // Determine indices of the cell containing (x, y)
        auto [i, j] = get_indices(x, y);

        // Determine the cell corners
        double x_0 = x_min + i * delta_x,
               y_0 = y_max - j * delta_y,
               x_1 = x_min + (i+1) * delta_x,
               y_1 = y_max - (j+1) * delta_y;

        // Using bilinear interpolation on the celll
        double h = data[(i + 0)*num_points_y + j + 0] * (x_1 - x)/delta_x * (y - y_1)/delta_y
                 + data[(i + 1)*num_points_y + j + 0] * (x - x_0)/delta_x * (y - y_1)/delta_y
                 + data[(i + 0)*num_points_y + j + 1] * (x_1 - x)/delta_x * (y_0 - y)/delta_y
                 + data[(i + 1)*num_points_y + j + 1] * (x - x_0)/delta_x * (y_0 - y)/delta_y;

        return h;
    }

    // Boundary of the raster domain as a CGAL polygon
    CGAL::SimplePolygon exterior() const {
        CGAL::SimplePolygon rectangle;
        rectangle.push_back(CGAL::Point2(x_min, get_y_min()));
        rectangle.push_back(CGAL::Point2(get_x_max(), get_y_min()));
        rectangle.push_back(CGAL::Point2(get_x_max(), y_max));
        rectangle.push_back(CGAL::Point2(x_min, y_max));

        return rectangle;
    }

    // Compute intersection of raster rectangle with polygon
    template<typename P>
    CGAL::MultiPolygon compute_intersection(const P& polygon) const {
        CGAL::SimplePolygon rectangle = exterior();

        CGAL::MultiPolygon intersection_polygon;
        CGAL::intersection(rectangle, polygon, std::back_inserter(intersection_polygon));

        return intersection_polygon;
    }

    // Determine if a point (x, y) is is strictly inside the raster domain
    bool contains(double x, double y) const {
        double eps = pow(pow(delta_x, 2) + pow(delta_y, 2), 0.5) * 1e-10;
        return ((x > x_min + eps) and (x < get_x_max() - eps)
                and (y > get_y_min() + eps) and (y < y_max - eps));
    }
};

template<typename T, typename P>
CGAL::DelaunayConstraints interpolate_boundary_points(const RasterData<T>& raster,
                                                      const P& boundary_polygon) {
    // First we need to determine intersection points between the raster domain and the polygon
    CGAL::MultiPolygon intersection_polygon = raster.compute_intersection(boundary_polygon);

    // Iterate over edges of the intersection polygon and interpolate points
    // TODO: Handle holes. Only needed if the intersecting polygon has holes.
    CGAL::DelaunayConstraints interpolated_points;
    for (const auto& part : CGAL::extract_boundaries(intersection_polygon)) {
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


template<typename Pgn>
Mesh make_mesh(const CGAL::PointList &pts,
               const Pgn& inclusion_polygon,
               const CGAL::DelaunayConstraints &constraints,
               const unsigned int epsg_id) {

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
    return Mesh(mesh, epsg_id);
};


Mesh make_mesh(const CGAL::PointList &pts,  const unsigned int epsg_id) {

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

    return Mesh(mesh, epsg_id);
};


template<typename T, typename Pgn>
Mesh mesh_from_raster(const std::vector<RasterData<T>>& raster_list,
                      const Pgn& boundary_polygon,
                      const unsigned int epsg_id) {
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
    return make_mesh(raster_points, boundary_polygon, boundary_points, epsg_id);
}


template<typename T>
Mesh mesh_from_raster(const std::vector<RasterData<T>>& raster_list, const unsigned int epsg_id) {
    CGAL::PointList raster_points;
    for (auto raster : raster_list) {
        CGAL::PointList new_points = raster.raster_points();
        raster_points.insert(raster_points.end(),
                             std::make_move_iterator(new_points.begin()),
                             std::make_move_iterator(new_points.end()));
    }
    return make_mesh(raster_points, epsg_id);
}


template<typename T, typename Pgn>
Mesh mesh_from_raster(const RasterData<T>& raster,
                      const Pgn& boundary_polygon,
                      const unsigned int epsg_id) {
    CGAL::PointList raster_points = raster.raster_points();
    CGAL::DelaunayConstraints boundary_points = interpolate_boundary_points(raster, boundary_polygon);

    return make_mesh(raster_points, boundary_polygon, boundary_points, epsg_id);
}


template<typename T>
Mesh mesh_from_raster(const RasterData<T>& raster, const unsigned int epsg_id) {
    return make_mesh(raster.raster_points(), epsg_id);
}

CGAL::Point3 centroid(const Mesh& mesh, const CGAL::face_descriptor &face) {
    CGAL::Point3 c{0, 0, 0};
    for (auto v: mesh.cgal_mesh.vertices_around_face(mesh.cgal_mesh.halfedge(face))) {
        const auto pt = mesh.cgal_mesh.point(v);
        c = CGAL::Point3{c.x() + pt.x(), c.y() + pt.y(), c.z() + pt.z()};
    }
    return CGAL::Point3{c.x()/3.0, c.y()/3.0, c.z()/3.0};
}

std::vector<int> compute_shadow(const Mesh & mesh,
                                const point3 &sun_direction) {
    std::vector<int> shade;
    auto cgal_mesh = mesh.cgal_mesh;
    const CGAL::Vector sun_vec(-sun_direction[0], -sun_direction[1], -sun_direction[2]);

    int i = 0;
    for (auto fd: cgal_mesh.faces()) {
        auto v = CGAL::Polygon_mesh_processing::compute_face_normal(fd, cgal_mesh);
        if ( v[0]*sun_vec[0] + v[1]*sun_vec[1] + v[2]*sun_vec[2] > 0.0 )
            shade.emplace_back(i);
        else {
            CGAL::Tree tree(CGAL::faces(cgal_mesh).first, CGAL::faces(cgal_mesh).second, cgal_mesh);
            const auto c = centroid(mesh, fd);
                CGAL::Ray sun_ray(c, -sun_vec);
                auto intersection = tree.first_intersection(sun_ray,
                                                            [fd] (const CGAL::face_descriptor &t) { return (t == fd); });
                if (intersection)
                    shade.emplace_back(i);
        }
        ++i;
    }
    return shade;
};

bool is_shaded(const CGAL::Tree &tree,
        const CGAL::face_descriptor &fd,
        const CGAL::Vector &face_normal,
        const CGAL::Point3 &face_center,
        const double azimuth,
        const double elevation) {
    const arma::vec::fixed<3> sd = arma::normalise(arma::vec::fixed<3>{sin(azimuth*M_PI/180.0),
                                                                       cos(azimuth*M_PI/180.0),
                                                                       tan(elevation*M_PI/180.0)});
    const CGAL::Vector sun_vec(-sd[0], -sd[1], -sd[2]);
    if ( face_normal[0]*sun_vec[0] + face_normal[1]*sun_vec[1] + face_normal[2]*sun_vec[2] > 0.0 )
        return true;
    CGAL::Ray sun_ray(face_center, -sun_vec);
    if (tree.first_intersection(sun_ray, [fd] (const CGAL::face_descriptor &t) { return (t == fd); }))
        return true;
    return false;
}

auto shade(const Mesh &mesh,
           const std::chrono::system_clock::time_point tp) {
    std::vector<bool> shade_vec;
    shade_vec.reserve(mesh.num_faces());
    namespace bg = boost::geometry;
    using point_car = bg::model::point<double, 2, bg::cs::cartesian>;
    using point_geo = bg::model::point<double, 2, bg::cs::geographic<bg::degree>>;
    bg::srs::transformation<> tr{
        bg::srs::epsg(mesh.epsg_id),
        bg::srs::epsg(4326)
    };
    CGAL::Tree tree(CGAL::faces(mesh.cgal_mesh).first, CGAL::faces(mesh.cgal_mesh).second, mesh.cgal_mesh);
    for (auto fd: mesh.cgal_mesh.faces()) {
        const auto c = centroid(mesh, fd);
        const point_car x_car{c.x(), c.y()};
        point_geo x_geo;
        tr.forward(x_car, x_geo);
        const auto [azimuth, elevation] = solar_position::time_point_solar_position(
                tp,
                bg::get<1>(x_geo),
                bg::get<0>(x_geo),
                c.z(),
                rasputin::solar_position::collectors::azimuth_and_elevation(),
                rasputin::solar_position::delta_t_calculator::coarse_timestamp_calc()
        );
        auto face_normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd, mesh.cgal_mesh);
        shade_vec.emplace_back(is_shaded(tree, fd, face_normal, c, azimuth, elevation));
    }
    return shade_vec;
}

std::vector<int> compute_shadow(const Mesh &mesh,
                                const double azimuth,
                                const double elevation) {
    // Topocentric azimuth and elevation
    const arma::vec::fixed<3> sd = arma::normalise(arma::vec::fixed<3>{sin(azimuth*M_PI/180.0), 
                                                                       cos(azimuth*M_PI/180.0),
                                                                       tan(elevation*M_PI/180.0)});
    return compute_shadow(mesh, point3{sd[0], sd[1], sd[2]});
};

std::vector<std::vector<int>> compute_shadows(const Mesh &mesh,
                                              const std::vector<std::pair<int, point3>> & sun_rays) {
    // TODO: Rewrite totally!
    std::vector<std::vector<int>>  result;
    return result;
    /*
    result.reserve(sun_rays.size());
    CGAL::Mesh cgal_mesh = mesh.cgal_mesh;
    std::map<size_t, CGAL::VertexIndex> index_map;
    std::map<CGAL::face_descriptor, size_t> face_map;
    size_t i = 0;
    size_t j = 0;
    CGAL::Tree tree(CGAL::faces(cgal_mesh).first, CGAL::faces(cgal_mesh).second, cgal_mesh);
    for (auto item: sun_rays) {
        std::vector<int> shade;
        shade.reserve(cgal_mesh.num_faces());
        auto utc_time = item.first;
        const CGAL::Vector sun_ray(item.second[0],
                                   item.second[1],
                                   item.second[2]);

        for (auto fd: CGAL::faces(cgal_mesh)) {
            auto hd = halfedge(fd, cgal_mesh);
            auto p = CGAL::centroid(cgal_mesh.point(source(hd, cgal_mesh)),
                                    cgal_mesh.point(target(hd, cgal_mesh)),
                                    cgal_mesh.point(target(next(hd, cgal_mesh), cgal_mesh)));
            auto v = CGAL::Polygon_mesh_processing::compute_face_normal(fd, cgal_mesh);
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
    */
};

point3_vector orient_tin(const point3_vector &pts, face_vector &faces) {
    point3_vector result;
    result.reserve(faces.size());
    for (auto& face: faces) {
        // Compute ccw normal
        const auto p0 = pts[face[0]];
        const auto p1 = pts[face[1]];
        const auto p2 = pts[face[2]];
        const arma::vec::fixed<3> v0{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        const arma::vec::fixed<3> n = arma::cross(v0, v1);
        double c = arma::norm(n);

        // Reverse triangle orientation if it is negatively oriented relative to xy plane
        if (n[2] < 0.0) {
          c *= -1.0;
          std::reverse(face.begin(), face.end());
        }

        // Store normalised and correctly oriented normal vector
        result.push_back(point3{n[0]/c, n[1]/c, n[2]/c});
    }
    return result;
};

double compute_slope(const point3 &normal) {
    return std::atan2(pow(pow(normal[0], 2) + pow(normal[1], 2), 0.5), normal[2]);
}

double_vector compute_slopes(const point3_vector &normals) {
    double_vector result;
    result.reserve(normals.size());

    for (const auto &n : normals)
        result.push_back(compute_slope(n));

    return result;
};

double compute_aspect(const point3 &normal) {
    return std::atan2(normal[0], normal[1]);
}

double_vector compute_aspects(const point3_vector &normals) {
    double_vector result;
    result.reserve(normals.size());
    for (const auto &n: normals)
        result.emplace_back(compute_aspect(n));
    return result;
};

point3 normal(const point3 &p0, const point3 &p1, const point3 &p2) {
        const arma::vec::fixed<3> v0{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        arma::vec::fixed<3> n = arma::cross(v0, v1);
        n /= arma::norm(n);
        return n[2] >= 0.0 ? point3{n[0], n[1], n[2]} : point3{-n[0], -n[1], -n[2]};
}

point3_vector surface_normals(const point3_vector &pts, const face_vector &faces) {
    point3_vector result;
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

point3_vector point_normals(const point3_vector &pts, const face_vector &faces) {
    point3_vector result(pts.size(), {0.0, 0.0, 0.0});
    for (auto face: faces) {
        const auto p0 = pts[face[0]];
        const auto p1 = pts[face[1]];
        const auto p2 = pts[face[2]];
        const arma::vec::fixed<3> v0{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        arma::vec::fixed<3> n = arma::cross(v0, v1);
        n /= arma::norm(n);
        const auto v = (n[2] >= 0.0) ? point3{n[0], n[1], n[2]} : point3{-n[0], -n[1], -n[2]};
        iadd(result[face[0]], v);
        iadd(result[face[1]], v);
        iadd(result[face[2]], v);
    }
    for (int i = 0; i < result.size(); ++i) {
        point3 &p = result[i];
        const double norm = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        if (norm > 1.0e-16) {
            p[0] /= norm;
            p[1] /= norm;
            p[2] /= norm;
        }
    }
    return result;
}

template <typename CB>
std::tuple<face_vector, face_vector> partition(const point3_vector &pts,
                                               const face_vector &faces,
                                               CB criterion) {
    face_vector part1, part2;
    for (auto face: faces) {
        if (criterion(pts[face[0]], pts[face[1]], pts[face[2]]))
            part1.emplace_back(face);
        else
            part2.emplace_back(face);
    }
    return std::make_pair(std::move(part1), std::move(part2));
}

std::tuple<face_vector, face_vector> extract_lakes(const point3_vector&pts, const face_vector &faces) {
   return partition(pts, faces, [] (const point3 &p0, const point3 &p1, const point3 &p2){
       return compute_slope(normal(p0, p1, p2)) < 1.0e-2;
   });
}

std::tuple<face_vector, face_vector> extract_avalanche_expositions(const point3_vector &pts,
        const face_vector &faces,
        const point2_vector &exposed_intervals,
        const point2_vector &height_intervals){
    return partition(pts, faces, [exposed_intervals, height_intervals](const point3 &p0, const point3 &p1, const point3 &p2){
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
        if (cell_slope < 30./180.*M_PI)
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

point3_vector cell_centers(const point3_vector& points, const face_vector & faces) {
    point3_vector result;
    result.reserve(faces.size());
    for (auto f: faces) {
        auto p0 = points[f[0]];
        auto p1 = points[f[1]];
        auto p2 = points[f[2]];
        auto x = (p0[0] + p1[0] + p2[0])/3.0;
        auto y = (p0[1] + p1[1] + p2[1])/3.0;
        auto z = (p0[2] + p1[2] + p2[2])/3.0;
        result.emplace_back(point3{x, y, z});
    }
    return result;
}

index_vector coordinates_to_indices(double x0,
        double y1,
        double dx,
        double dy,
        unsigned int M,
        unsigned int N,
        point2_vector pts) {

    index_vector indices;
    indices.reserve(pts.size());
    for (auto pt: pts)
        indices.emplace_back(std::array<unsigned int, 2>{(unsigned int)((pt[0]-x0)/dx), (unsigned int)((y1 - pt[1])/dy)});
    return indices;
}

template <typename T>
std::vector<T> extract_buffer_values(const index_vector& indices, pybind11::array_t<T>& array) {
    std::vector<T> result;
    result.reserve(indices.size());
    auto buffer = array.request();
    unsigned long M = (unsigned long)buffer.shape[0];
    unsigned long N = (unsigned long)buffer.shape[1];
    T* ptr = (T *)buffer.ptr;
    for (auto idx: indices)
        result.emplace_back(ptr[idx[0]*N + idx[1]]);  // TODO: Implement range check?
    return result;
}

std::tuple<point3_vector, face_vector> consolidate(const point3_vector &points, const face_vector &faces){
    face_vector new_faces;
    point3_vector new_points;
    new_faces.reserve(faces.size());
    std::map<int, int> point_map;
    int n = 0;
    for (const auto _face: faces) {
        face new_face;
        int i = 0;
        for (const auto f: _face) {
            if (not point_map.count(f)) {
                point_map.insert(std::make_pair(f, n++));
                new_points.emplace_back(points[f]);
            }
            new_face[i++] = point_map[f];
        }
        new_faces.emplace_back(new_face);
    }
    return std::make_pair(std::move(new_points), std::move(new_faces));
}
}
