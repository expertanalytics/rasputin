#pragma once

#include <concepts>
#include <tuple>
#include <cmath>
#include <vector>
#include <algorithm>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/srs/epsg.hpp>
#include <boost/geometry/srs/projection.hpp>
#include <boost/geometry/srs/transformation.hpp>

namespace rasputin {
namespace bg = boost::geometry;

template<typename Polygon>
constexpr bool CLOCKWISE = std::same_as<
    std::integral_constant<std::size_t, bg::point_order<Polygon>::value>,
    std::integral_constant<std::size_t, bg::clockwise>
>;

template<std::size_t Dim>
using DimConst = std::integral_constant<std::size_t, Dim>;

template <typename Point, typename FT>
requires
    std::same_as<DimConst<bg::dimension<Point>::value>, DimConst<3>>
    && std::floating_point<FT>
    && std::same_as<typename bg::coordinate_type<Point>::type, FT>
struct RasterData {
    public:
        using MultiPoint = bg::model::multi_point<Point>;
        using LineString = bg::model::linestring<Point>;
        using Box = bg::model::box<Point>;
        using Polygon = bg::model::polygon<Point>;

    private:

        FT x_min;
        FT x_max;
        FT delta_x;
        std::size_t num_points_x;

        FT y_min;
        FT y_max;
        FT delta_y;
        std::size_t num_points_y;

        FT* data;
        MultiPoint points;

    public:
        RasterData(
            FT x_min, FT y_max, FT delta_x, FT delta_y, std::size_t num_points_x, std::size_t num_points_y, FT* data
        ) :
            x_min(x_min),
            x_max(x_min + (num_points_x - 1) * delta_x),
            y_min(y_max - (num_points_y - 1) * delta_y),
            y_max(y_max),
            delta_x(delta_x),
            delta_y(delta_y),
            num_points_x(num_points_x),
            num_points_y(num_points_y),
            data(data)
        {
            for (std::size_t i = 0; i < num_points_x; ++i) {
                for (std::size_t j = 0; j < num_points_y; ++j) {
                    bg::append(
                        this->points,
                        Point(x_min + i * delta_x, y_max - j * delta_y, data[j * num_points_x + i])
                    );
                }
            }
        }

        const MultiPoint& get_points() const {
            return points;
        }

        const std::size_t& get_num_x() const {
            return this->num_points_x;
        }

        const std::size_t& get_num_y() const {
            return this->num_points_y;
        }

        const FT& get_xmin() const {
            return this->x_min;
        }

        const FT& get_ymin() const {
            return this->y_min;
        }

        const FT& get_xmax() const {
            return this->x_max;
        }

        const FT& get_ymax() const {
            return this->y_max;
        }

        // For every point inside the raster rectangle we identify indices (i, j) of
        // the upper-left vertex of the cell containing the point
        std::pair<int, int> get_indices(FT x, FT y) const {
            int i = std::min<int>(std::max<int>(static_cast<int>((y_max - y) / delta_y), 0), num_points_y - 1);
            int j = std::min<int>(std::max<int>(static_cast<int>((x - x_min) / delta_x), 0), num_points_x - 1);

            return std::make_pair(i, j);
        }

        // Interpolate data using using a bilinear interpolation rule on each cell
        FT get_interpolated_value_at_point(FT x, FT y) const {
            // Determine indices of the cell containing (x, y)
            auto [i, j] = get_indices(x, y);

            // Determine the cell corners
            //     (x0, y0) -- upper left
            //     (x1, y1) -- lower right
            const FT x_0 = x_min + j * delta_x, y_0 = y_max - i * delta_y;
            const FT x_1 = x_min + (j + 1) * delta_x, y_1 = y_max - (i + 1) * delta_y;

            // Using bilinear interpolation on the cell
            const FT h = data[(i + 0) * num_points_x + j + 0] * (x_1 - x) / delta_x * (y - y_1) / delta_y // (x0, y0)
                       + data[(i + 0) * num_points_x + j + 1] * (x - x_0) / delta_x * (y - y_1) / delta_y // (x1, y0)
                       + data[(i + 1) * num_points_x + j + 0] * (x_1 - x) / delta_x * (y_0 - y) / delta_y // (x0, y1)
                       + data[(i + 1) * num_points_x + j + 1] * (x - x_0) / delta_x * (y_0 - y) / delta_y; // (x1, y1)

            return h;
        }

        // Boundary of the raster domain as a boost geometry polygon
        Polygon exterior() const {
            Polygon rectangle{{{x_min, y_min}, {x_max, y_min}, {x_max, y_max}, {x_min, y_max}}};
            return rectangle;
        }

        // Compute intersection of raster rectangle with polygon
        template <typename P>
        requires std::same_as<P, Polygon>
        Polygon compute_intersection(const P &polygon) const {
            Box rectangle = exterior();

            Polygon intersection_polygon;
            bg::intersection(rectangle, polygon, intersection_polygon);

            return intersection_polygon;
        }

        // Determine if a point (x, y) is strictly inside the raster domain
        bool contains(FT x, FT y) const {
            FT eps = pow(pow(delta_x, 2) + pow(delta_y, 2), 0.5) * 1e-10;
            return ((x > x_min + eps) and (x < x_max - eps) and (y > y_min + eps) and (y < y_max - eps));
        }
};

template <
    typename Point,
    typename FT,
    typename MultiPoint=bg::model::multi_point<Point>,
    typename Segment=bg::model::segment<Point>,
    typename LineString=bg::model::linestring<Point>,
    typename MultiLineString=bg::model::multi_linestring<Point>,
    typename Box=bg::model::box<Point>,
    typename Polygon=bg::model::polygon<Point>,
    typename MultiPolygon=bg::model::multi_polygon<Polygon>
>
requires std::same_as<DimConst<bg::dimension<Point>::value>, DimConst<3>> && std::floating_point<FT>
class Triangulation {
public:
    using PointPtr = const Point*;
    using Face = std::tuple<PointPtr, PointPtr, PointPtr>;
    using Faces = std::vector<Face>;
    using Edge = std::pair<PointPtr, PointPtr>;

    struct Mesh {
        private:
            const std::string proj4_str;

        public:
            const Faces triangles;

            explicit Mesh(Faces& triangles, const std::string proj4_str) :
                proj4_str(proj4_str),
                triangles(triangles)
            {
            }

            size_t num_faces() const {
                return this->triangles.size();
            }
    };

    inline static bool _is_inside_ccirc(const PointPtr& point, const Face& triangle) {
        const auto [a, b, c] = triangle;
        const FT px = bg::get<0>(point);
        const FT py = bg::get<1>(point);
        const FT pxa = bg::get<0>(a) - px;
        const FT pya = bg::get<1>(a) - py;
        const FT pxb = bg::get<0>(b) - px;
        const FT pyb = bg::get<1>(b) - py;
        const FT pxc = bg::get<0>(c) - px;
        const FT pyc = bg::get<1>(c) - py;

        // |pxa  pya  a2  pxa*pxa + pya*pya|
        // |pxb  pyb  b2  pxb*pxb + pyb*pyb|
        // |pxc  pyc  c2  pxc*pxc + pyc*pyc|
        const FT det =
            (pxa*pxa + pya*pya)*(pxb*pyc - pxc*pyb)
            - (pxb*pxb + pyb*pyb)*(pxa*pyc - pxc*pya)
            + (pxc*pxc + pyc*pyc)*(pxa*pyb - pxb*pya);

        if constexpr(CLOCKWISE<Polygon>) {
            return det < 0 ? true : false;
        } else {
            return det > 0 ? true : false;
        }
    }

    inline static Face _sort_points(const PointPtr& a, const PointPtr& b, const PointPtr& c) {
        using tri_type = std::pair<PointPtr, FT>;

        const FT ax = bg::get<0>(*a);
        const FT ay = bg::get<1>(*a);
        const FT bx = bg::get<0>(*b);
        const FT by = bg::get<1>(*b);
        const FT cx = bg::get<0>(*c);
        const FT cy = bg::get<1>(*c);

        FT ccw = (bx - ax)*(cy - ay) - (cx - ax)*(by - ay);

        if constexpr(CLOCKWISE<Polygon>) {
            return ccw < 0 ? Face{a, c, b} : Face{a, b, c};
        } else {
            return ccw > 0 ? Face{a, c, b} : Face{a, b, c};
        }
    }

    template <std::size_t DN = 20>
    static Mesh delaunay(const RasterData<Point, FT>& raster, const std::string proj4_str) {
        auto add_new_triangle = [](Faces& triangles, const Edge& edge, const PointPtr& point) {
            /* add new triangles formed by 2 points in edge and point. Sort points to be clockwise with respect to
               centroid */
            Face triangle = _sort_points(edge.first, edge.second, point);

            triangles.push_back(triangle);
        };

        // add super triangle
        const FT& lower_x = raster.get_xmin();
        const FT& lower_y = raster.get_ymin();
        const FT& upper_x = raster.get_xmax();
        const FT& upper_y = raster.get_ymax();

        const FT diff_max = std::max(upper_x - lower_x, upper_y - lower_y);
        const FT mid_x = 0.5 * (upper_x + lower_x);
        const FT mid_y = 0.5 * (upper_y + lower_y);

        Faces triangles;
        const Point sta = {mid_x - DN*diff_max, mid_y - diff_max};
        const Point stb = {mid_x, mid_y + DN*diff_max};
        const Point stc = {mid_x + DN*diff_max, mid_y - diff_max};
        const Face super_triangle = std::make_tuple(&sta, &stb, &stc);
        triangles.push_back(super_triangle);

        const MultiPoint& points = raster.get_points();

        add_new_triangle(triangles, {std::get<0>(super_triangle), std::get<1>(super_triangle)}, &(points[0]));
        add_new_triangle(triangles, {std::get<1>(super_triangle), std::get<2>(super_triangle)}, &(points[0]));
        add_new_triangle(triangles, {std::get<2>(super_triangle), std::get<0>(super_triangle)}, &(points[0]));

        for (typename MultiPoint::const_iterator point=points.begin() + 1; point != points.end(); ++point) {
            auto start_bad = std::remove_if(
                triangles.begin(),
                triangles.end(),
                [&point](const Face& triangle) -> bool {
                    return _is_inside_ccirc(&*point, triangle);
                }
            );

            std::vector<Edge> new_edges{};
            for (typename Faces::iterator triangle1=start_bad; triangle1 != triangles.end(); ++triangle1) {
                const auto [a1, b1, c1] = *triangle1;
                for (const Edge e1 : {std::make_pair(a1, b1), std::make_pair(b1, c1), std::make_pair(c1, a1)}) {
                    bool edge_good = true;
                    for (
                        typename Faces::iterator triangle2=triangle1 + 1;
                        triangle2 != triangles.end() && edge_good;
                        ++triangle2
                    ) {
                        const auto [a2, b2, c2] = *triangle2;
                        for (const Edge e2 : {
                            std::make_pair(a2, b2), std::make_pair(b2, c2), std::make_pair(c2, a2)
                        }) {
                            if (
                                ((e1.first == e2.first) && (e1.second == e2.second))
                                || ((e1.second == e2.first) && (e1.first == e2.second))
                            ) {
                                edge_good = false;
                                break;
                            }
                        }
                    }

                    if (edge_good) {
                        new_edges.push_back(e1);
                    }
                }
            }

            triangles.erase(start_bad, triangles.end());

            for (const Edge& edge : new_edges) {
                add_new_triangle(triangles, edge, &*point);
            }
        }

        // remove all triangles connected to super triangle
        std::erase_if(
            triangles,
            [&super_triangle](const Face& triangle) -> bool {
                const auto [a, b, c] = triangle;
                const auto [s1, s2, s3] = super_triangle;

                if (
                    (a == s1)
                    || (a == s2)
                    || (a == s3)
                    || (b == s1)
                    || (b == s2)
                    || (b == s3)
                    || (c == s1)
                    || (c == s2)
                    || (c == s3)
                ) {
                    return true;
                }

                return false;
            }
        );

        return Mesh(triangles, proj4_str);
    }

    template<typename P>
    requires std::same_as<P, Polygon>
    Mesh DelaunayConstrained(const RasterData<Point, FT> &raster, P &bound);

    template<typename P>
    requires std::same_as<P, Segment>
    Mesh DelaunayConstrained(const RasterData<Point, FT> &raster, P &bound);

    template<typename P>
    requires std::same_as<P, LineString>
    Mesh DelaunayConstrained(const RasterData<Point, FT> &raster, P &bound);

    template<typename P>
    requires std::same_as<P, MultiLineString>
    Mesh DelaunayConstrained(const RasterData<Point, FT> &raster, P &bound);

  // template<typename FT, typename P>
  // requires std::same_as<P, Polygon> && std::floating_point<FT>
  // std::vector<Points> interpolate_boundary_points(const RasterData<Point,
  // FT>& raster, const P& boundary_polygon) { MultiPolygon intersection_polygon
  // = raster.compute_intersection(boundary_polygon);
  //
  // std::vector<Points> interpolated_points;
  //
  // return interpolated_points;
  // }
};
}; // namespace rasputin
