#pragma once

#include <CGAL/Boolean_set_operations_2.h>

#include "types.h"

namespace rasputin {
template<typename FT>
struct RasterData {
    RasterData(
        double x_min,
        double y_max,
        double delta_x,
        double delta_y,
        std::size_t num_points_x,
        std::size_t num_points_y,
        FT* data
    ) :
        x_min(x_min),
        delta_x(delta_x),
        num_points_x(num_points_x),
        y_max(y_max),
        delta_y(delta_y),
        num_points_y(num_points_y),
        data(data) {
    }

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

        for (std::size_t i = 0; i < num_points_y; ++i)
            for (std::size_t j = 0; j < num_points_x; ++j)
                points.emplace_back(x_min + j*delta_x, y_max - i*delta_y, data[i*num_points_x + j]);
        return points;
    }

    // For every point inside the raster rectangle we identify indices (i, j) of the upper-left vertex of the cell containing the point
    std::pair<int, int> get_indices(double x, double y) const {
        int i = std::min<int>(std::max<int>(static_cast<int>((y_max-y) /delta_y), 0), num_points_y - 1);
        int j = std::min<int>(std::max<int>(static_cast<int>((x-x_min) /delta_x), 0), num_points_x - 1);

        return std::make_pair(i,j);
    }

    // Interpolate data using using a bilinear interpolation rule on each cell
    FT get_interpolated_value_at_point(double x, double y) const {
        // Determine indices of the cell containing (x, y)
        auto [i, j] = get_indices(x, y);

        // Determine the cell corners
        //     (x0, y0) -- upper left
        //     (x1, y1) -- lower right
        double x_0 = x_min + j * delta_x,
               y_0 = y_max - i * delta_y,
               x_1 = x_min + (j+1) * delta_x,
               y_1 = y_max - (i+1) * delta_y;

        // Using bilinear interpolation on the celll
        double h = data[(i + 0)*num_points_x + j + 0] * (x_1 - x)/delta_x * (y - y_1)/delta_y   // (x0, y0)
                 + data[(i + 0)*num_points_x + j + 1] * (x - x_0)/delta_x * (y - y_1)/delta_y   // (x1, y0)
                 + data[(i + 1)*num_points_x + j + 0] * (x_1 - x)/delta_x * (y_0 - y)/delta_y   // (x0, y1)
                 + data[(i + 1)*num_points_x + j + 1] * (x - x_0)/delta_x * (y_0 - y)/delta_y;  // (x1, y1)

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
}
