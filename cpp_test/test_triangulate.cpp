#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch_test_macros.hpp>

#include <vector>
#include <random>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#include <triangulate_dem.h>

namespace bg = boost::geometry;
using Point = bg::model::point<double, 3, bg::cs::cartesian>;
namespace boost::geometry::traits {
    using shared_type = std::shared_ptr<const Point>;
    BOOST_GEOMETRY_DETAIL_SPECIALIZE_POINT_TRAITS(shared_type, 3, double, cs::cartesian)
    template<>
    struct access<shared_type, 0> {
        static inline double get(shared_type const& p) { return p->get<0>(); }
    };
    template<>
    struct access<shared_type, 1> {
        static inline double get(shared_type const& p) { return p->get<1>(); }
    };
    template<>
    struct access<shared_type, 2> {
        static inline double get(shared_type const& p) { return p->get<2>(); }
    };
}

namespace rasputin::test_utils {
    namespace bg = boost::geometry;
    using Point = bg::model::point<double, 3, bg::cs::cartesian>;
    using Triangulation = rasputin::Triangulation<Point, double>;

    void check_points_exists(const std::vector<Triangulation::Face>& triangles, const std::vector<Point>& points) {
        for (const Point& point : points) {
            bool found = false;
            for (const Triangulation::Face& triangle : triangles) {
                const auto& [a, b, c] = triangle;
                if (&point == a || &point == b || &point == c) {
                    found = true;
                    break;
                }
            }

            INFO("Point (" << bg::get<0>(point) << " " << bg::get<1>(point) << " " << bg::get<2>(point) << ") not found");
            CHECK(found);
        }
    }

    void check_points_inside_ccirc(
        const std::vector<Triangulation::Face>& triangles, const std::vector<Point>& points
    ) {
        for (
            Triangulation::Faces::const_iterator triangle1=triangles.begin();
            triangle1 != triangles.end();
            ++triangle1
        ) {
            auto [t1a, t1b, t1c] = *triangle1;
            for (
                Triangulation::Faces::const_iterator triangle2=triangle1 + 1;
                triangle2 != triangles.end();
                ++triangle2
            ) {
                bool ta_inside = Triangulation::_is_inside_ccirc(t1a, *triangle2);
                bool tb_inside = Triangulation::_is_inside_ccirc(t1b, *triangle2);
                bool tc_inside = Triangulation::_is_inside_ccirc(t1c, *triangle2);
                INFO("Point (" << bg::get<0>(t1a) << " " << bg::get<1>(t1a) << " " << bg::get<2>(t1a) << ") inside");
                CHECK(!ta_inside);
                INFO("Point (" << bg::get<0>(t1a) << " " << bg::get<1>(t1a) << " " << bg::get<2>(t1a) << ") inside");
                CHECK(!tb_inside);
                INFO("Point (" << bg::get<0>(t1c) << " " << bg::get<1>(t1c) << " " << bg::get<2>(t1c) << ") inside");
                CHECK(!tc_inside);
            }
        }
    }

    void check_triangles_unique(const std::vector<Triangulation::Face>& triangles) {
        for (
            Triangulation::Faces::const_iterator triangle1=triangles.begin();
            triangle1 != triangles.end();
            ++triangle1
        ) {
            const auto& [a1, b1, c1] = *triangle1;
            for (
                Triangulation::Faces::const_iterator triangle2=triangle1 + 1;
                triangle2 != triangles.end();
                ++triangle2
            ) {
                const auto& [a2, b2, c2] = *triangle2;
                bool same =
                    (&a1 == &a2 && &b1 == &b2 && &c1 == &c2)
                    || (&a1 == &c2 && &b1 == &a2 && &c1 == &b2)
                    || (&a1 == &b2 && &b1 == &c2 && &c1 == &a2);
                INFO(
                    "Duplicate Triangle [("
                    << bg::get<0>(a1) << " " << bg::get<1>(a1) << " " << bg::get<2>(a1) << ") "
                    << bg::get<0>(b1) << " " << bg::get<1>(b1) << " " << bg::get<2>(b1) << ") "
                    << bg::get<0>(c1) << " " << bg::get<1>(c1) << " " << bg::get<2>(c1) << ")]"
                );
                CHECK(!same);
            }
        }
    }

    TEST_CASE("[Test sort polar]", "[Triangulation]") {
        std::random_device rd;
        std::default_random_engine engine(rd());
        std::uniform_real_distribution distr(0.0, 1.0);

        const Point a = {0.0, 0.0, 0.0};
        const Point b = {distr(engine), distr(engine), 0.0};
        const Point c = {distr(engine), 0.0, 0.0};

        Triangulation::Face triangle = Triangulation::_sort_points(&a, &b, &c);
        REQUIRE(std::get<0>(triangle) == &a);
        REQUIRE(std::get<1>(triangle) == &c);
        REQUIRE(std::get<2>(triangle) == &b);

        triangle = Triangulation::_sort_points(&a, &c, &b);
        REQUIRE(std::get<0>(triangle) == &a);
        REQUIRE(std::get<1>(triangle) == &c);
        REQUIRE(std::get<2>(triangle) == &b);

        triangle = Triangulation::_sort_points(&c, &a, &b);
        REQUIRE(std::get<0>(triangle) == &c);
        REQUIRE(std::get<1>(triangle) == &b);
        REQUIRE(std::get<2>(triangle) == &a);

        triangle = Triangulation::_sort_points(&c, &b, &a);
        REQUIRE(std::get<0>(triangle) == &c);
        REQUIRE(std::get<1>(triangle) == &b);
        REQUIRE(std::get<2>(triangle) == &a);
    }

    TEST_CASE("[Test inside circumcircle]", "[Triangulation]") {
        const Point a = {0.0, 0.0, 0.0};
        const Point b = {0.5, 1.0, 0.0};
        const Point c = {1.0, 0.0, 0.0};
        Triangulation::Face triangle = {&a, &b, &c};

        Point inside_point = {0.25, 0.25, 0.0};
        bool is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == true);

        inside_point = {100.25, 100.25, 0.0};
        is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == false);

        inside_point = {0.85, 0.25, 0.0};
        is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == true);

        inside_point = {0.15, 0.25, 0.0};
        is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == true);

        inside_point = {0.95, 0.25, 0.0};
        is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == true);

        inside_point = {0.05, 0.25, 0.0};
        is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == true);

        inside_point = {1.15, 0.25, 0.0};
        is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == false);

        inside_point = {-0.85, 0.25, 0.0};
        is_inside = Triangulation::_is_inside_ccirc(&inside_point, triangle);
        REQUIRE(is_inside == false);
    }

    TEST_CASE("[Test unconstrained Delaunay triangulation]", "[Triangulation]") {
        SECTION("square") {
            std::vector<double> raw_data = {0, 0, 0, 0, 0, 0, 0, 0};
            RasterData<Point, double> data = {0.0, 1.0, 1.0, 1.0, 2, 2, raw_data.data()};

            auto mesh = Triangulation::delaunay(data, "");

            check_points_exists(mesh.triangles, data.get_points());
            check_points_inside_ccirc(mesh.triangles, data.get_points());
            check_triangles_unique(mesh.triangles);

            REQUIRE(mesh.num_faces() == 2);
        }

        SECTION("rectangle 1") {
            std::vector<double> raw_data = {0, 0, 0, 0, 0, 0, 0, 0};
            RasterData<Point, double> data = {0.0, 1.0, 2.0, 1.0, 2, 2, raw_data.data()};

            auto mesh = Triangulation::delaunay(data, "");

            check_points_exists(mesh.triangles, data.get_points());
            check_points_inside_ccirc(mesh.triangles, data.get_points());
            check_triangles_unique(mesh.triangles);

            REQUIRE(mesh.num_faces() == 2);
        }

        SECTION("rectangle 2") {
            std::vector<double> raw_data = {0, 0, 0, 0, 0, 0, 0, 0};
            RasterData<Point, double> data = {0.0, 2.0, 1.0, 2.0, 2, 2, raw_data.data()};

            auto mesh = Triangulation::delaunay(data, "");

            check_points_exists(mesh.triangles, data.get_points());
            check_points_inside_ccirc(mesh.triangles, data.get_points());
            check_triangles_unique(mesh.triangles);

            REQUIRE(mesh.num_faces() == 2);
        }

        SECTION("trapezoid 1") {
            std::vector<double> raw_data = {0, 0, 0, 0, 0, 0, 0, 0};
            RasterData<Point, double> data = {0.0, 3.0, 1.0, 3.0, 2, 2, raw_data.data()};

            auto mesh = Triangulation::delaunay(data, "");

            check_points_exists(mesh.triangles, data.get_points());
            check_points_inside_ccirc(mesh.triangles, data.get_points());
            check_triangles_unique(mesh.triangles);

            REQUIRE(mesh.num_faces() == 2);
        }
    }

    // TEST_CASE("[Test random unconstrained Delaunay triangulation]", "[Triangulation]") {
    //     std::random_device rd;
    //     std::default_random_engine engine(rd());
    //     std::uniform_real_distribution distr(0.0, 1.0);
    //
    //     std::vector<double> raw_data = {0, 0, 0, 0, 0, 0, 0, 0};
    //     RasterData<Point, double> data = {0.0, 1.0, 1.0, 1.0, 2, 2, raw_data.data()};
    //
    //     auto mesh = Triangulation::delaunay(data, "");
    //
    //     REQUIRE(mesh.num_faces() == 2);
    //
    //     check_points_exists(mesh.triangles, data.points);
    //     check_points_inside_ccirc(mesh.triangles, data.points);
    // }
}
