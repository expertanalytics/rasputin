//
// Created by Ola Skavhaug on 08/10/2018.
//

#pragma once

#include <sstream>
#include <cmath>
#include <fstream>
#include <map>
#include <tuple>
#include <numeric>
#include <cstdint>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_2_algorithms.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include "types.h"

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
