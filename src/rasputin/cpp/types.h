#pragma once

#include <vector>
#include <map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>

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
}
