#pragma once

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "types.h"
#include "polygon.h"

namespace rasputin::traits {
using Point2 = CGAL::Point2;
using Polygon = CGAL::Polygon;
using SimplePolygon = CGAL::SimplePolygon;

template<>
struct PolygonUtils<Polygon, SimplePolygon> {
    using Polygon_t = Polygon;
    using SimplePolygon_t = SimplePolygon;
    using MultiPolygon_t = std::vector<Polygon>;

    static bool point_inside_polygon(const Point2 &x, const SimplePolygon &polygon) {
        return polygon.has_on_bounded_side(x);
    }

    static bool point_inside_polygon(const Point2 &x, const Polygon &polygon) {
        if (not polygon.outer_boundary().has_on_bounded_side(x))
            return false;
        if (not polygon.has_holes())
            return true;
        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
            if (it->has_on_bounded_side(x))
                return false;
        return true;
    }

    static std::vector<SimplePolygon> extract_boundaries(const Polygon &polygon) {
        std::vector<SimplePolygon> ret;
        ret.emplace_back(polygon.outer_boundary());
        for (auto it = polygon.holes_begin(); it != polygon.holes_end(); ++it)
            ret.emplace_back(*it);

        return ret;
    }

    static std::vector<SimplePolygon> extract_boundaries(const MultiPolygon_t &polygon) {
        std::vector<SimplePolygon> ret;
        for (const auto& part: polygon) {
            ret.emplace_back(part.outer_boundary());
            for (auto it = part.holes_begin(); it != part.holes_end(); ++it)
                ret.emplace_back(*it);
        }

        return ret;
    }

    static SimplePolygon exterior(const Polygon& polygon) {
        return polygon.outer_boundary();
    }

    static py::list holes(const Polygon& polygon) {
        py::list result;
        for (auto h = polygon.holes_begin(); h != polygon.holes_end(); ++h)
            result.append(*h);
        return result;
    }

    static py::array_t<double> polygon_array(const SimplePolygon& polygon) {
        py::array_t<double> result(polygon.size() * 2);
        result.resize(py::array::ShapeContainer({static_cast<long int>(polygon.size()), 2}), true);
        auto info = result.request();
        double* data = static_cast<double*>(info.ptr);
        for (auto v = polygon.vertices_begin(); v != polygon.vertices_end(); ++v) {
            data[0] = v->x();
            data[1] = v->y();
            data += 2;
        }
        return result;
    }

    template<typename P0, typename P1>
    static MultiPolygon_t difference_polygons(const P0& polygon0, const P1& polygon1) {
        MultiPolygon_t difference_result;
        CGAL::difference(polygon0, polygon1, std::back_inserter(difference_result));

        return difference_result;
    }

    template<typename P0, typename P1>
    static MultiPolygon_t intersect_polygons(const P0& polygon0, const P1& polygon1) {
        MultiPolygon_t intersection_result;
        CGAL::intersection(polygon0, polygon1, std::back_inserter(intersection_result));

        return intersection_result;
    }

    template<typename P0, typename P1>
    static MultiPolygon_t join_polygons(const P0& polygon0, const P1& polygon1) {
        Polygon joined;
        MultiPolygon_t join_result;
        if (CGAL::join(polygon0, polygon1, joined))
            join_result.push_back(std::move(joined));
        else {
            join_result.push_back(Polygon(polygon0));
            join_result.push_back(Polygon(polygon1));
        }
        return join_result;
    }

    static MultiPolygon_t join_multipolygons(const MultiPolygon_t& polygon0, const MultiPolygon_t& polygon1) {
        MultiPolygon_t join_result;
        CGAL::join(polygon0.begin(), polygon0.end(), polygon1.begin(), polygon1.end(), std::back_inserter(join_result));

        return join_result;
    }

    static SimplePolygon polygon_from_numpy(py::array_t<double>& buf) {
            auto info = buf.request();
            if (info.ndim < 1 or info.ndim > 2)
                throw py::type_error("Can only convert one- and two-dimensional arrays.");

            // Make sure total size can be preserved
            auto size = info.shape[0];
            if (info.ndim == 1 and size % 2)
                throw py::type_error("Size of one-dimensional array must be divisible by 2.");

            if (info.ndim == 2 and info.shape[1] != 2)
                throw py::type_error("Second dimension does not have size equal to 2.");

            if (info.ndim == 2)
                size *= 2;

            SimplePolygon exterior;

            // Copy values
            double* p = static_cast<double*>(info.ptr);
            double* end = p + size;
            for (; p < end; p += 2)
                exterior.push_back(CGAL::Point2(*p, *(p+1)));
            return exterior;
    }

    static std::size_t num_vertices(const SimplePolygon& polygon) {
        return polygon.size();
    }
};
} // namespace rasputin::traits
