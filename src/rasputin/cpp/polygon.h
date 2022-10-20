#pragma once

#include "types.h"

namespace CGAL {
inline bool point_inside_polygon(const Point2 &x, const SimplePolygon &polygon) {
    return polygon.has_on_bounded_side(x);
}

inline bool point_inside_polygon(const Point2 &x, const Polygon &polygon) {
    if (not polygon.outer_boundary().has_on_bounded_side(x))
        return false;
    if (not polygon.has_holes())
        return true;
    for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
        if (it->has_on_bounded_side(x))
            return false;
    return true;
}

inline bool point_inside_polygon(const Point2 &x, const MultiPolygon &polygon) {
    for (const auto& part: polygon)
        if (point_inside_polygon(x, part))
            return true;
    return false;
}

inline std::vector<SimplePolygon> extract_boundaries(const SimplePolygon &polygon) {
    std::vector<SimplePolygon> ret;
    ret.emplace_back(polygon);

    return ret;
}

inline std::vector<SimplePolygon> extract_boundaries(const Polygon &polygon) {
    std::vector<SimplePolygon> ret;
    ret.emplace_back(polygon.outer_boundary());
    for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
        ret.emplace_back(*it);

    return ret;
}

inline std::vector<SimplePolygon> extract_boundaries(const MultiPolygon &polygon) {
    std::vector<SimplePolygon> ret;
    for (const auto& part: polygon) {
        ret.emplace_back(part.outer_boundary());
        for (auto it = part.holes_begin(); it != part.holes_end(); ++it)
            ret.emplace_back(*it);
    }

    return ret;
}
}
