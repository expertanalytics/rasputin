#pragma once

#include <CGAL/Boolean_set_operations_2.h>

#include "raster_data.h"

namespace rasputin {
using PointList = CGAL::PointList;
using Polygon = CGAL::Polygon;
using SimplePolygon = CGAL::SimplePolygon;
using MultiPolygon = CGAL::MultiPolygon;

template<typename FT>
requires std::floating_point<FT>
struct RasterData: public traits::RasterBase<FT, PointList, PolygonUtils<Polygon, SimplePolygon>> {
    RasterData(
        double x_min,
        double y_max,
        double delta_x,
        double delta_y,
        std::size_t num_points_x,
        std::size_t num_points_y,
        FT* data
    ):
        traits::RasterBase<FT, PointList, PolygonUtils<Polygon, SimplePolygon>>(
            x_min, y_max, delta_x, delta_y, num_points_x, num_points_y, data
        )
    {
    }

    PointList raster_points() const {
        PointList points;
        points.reserve(this->num_points_x * this->num_points_y);

        for (std::size_t i = 0; i < this->num_points_y; ++i)
            for (std::size_t j = 0; j < this->num_points_x; ++j)
                points.emplace_back(
                    this->x_min + j*this->delta_x, this->y_max - i*this->delta_y, this->data[i*this->num_points_x + j]
                );
        return points;
    }

    // Boundary of the raster domain as a CGAL polygon
    SimplePolygon exterior() const {
        SimplePolygon rectangle;
        rectangle.push_back(CGAL::Point2(this->x_min, this->get_y_min()));
        rectangle.push_back(CGAL::Point2(this->get_x_max(), this->get_y_min()));
        rectangle.push_back(CGAL::Point2(this->get_x_max(), this->y_max));
        rectangle.push_back(CGAL::Point2(this->x_min, this->y_max));

        return rectangle;
    }

    // Compute intersection of raster rectangle with polygon
    template<traits::PolygonType<Polygon, SimplePolygon> P>
    MultiPolygon compute_intersection(const P& polygon) const {
        SimplePolygon rectangle = exterior();

        MultiPolygon intersection_polygon;
        intersection(rectangle, polygon, std::back_inserter(intersection_polygon));

        return intersection_polygon;
    }
};
}
