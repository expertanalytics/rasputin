#include <boost/geometry.hpp>
#include <boost/geometry/srs/epsg.hpp>
#include <boost/geometry/srs/projection.hpp>
#include <boost/geometry/srs/transformation.hpp>

#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <armadillo>

#include "types.h"
#include "mesh.h"
#include "solar_position.h"

namespace rasputin {
inline CGAL::Point3 centroid(const Mesh& mesh, const CGAL::face_descriptor &face) {
    CGAL::Point3 c{0, 0, 0};
    for (auto v: mesh.cgal_mesh.vertices_around_face(mesh.cgal_mesh.halfedge(face))) {
        const auto pt = mesh.cgal_mesh.point(v);
        c = CGAL::Point3{c.x() + pt.x(), c.y() + pt.y(), c.z() + pt.z()};
    }
    return CGAL::Point3{c.x()/3.0, c.y()/3.0, c.z()/3.0};
}

inline std::vector<int> compute_shadow(const Mesh & mesh, const point3 &sun_direction) {
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

inline bool is_shaded(
    const CGAL::Tree &tree,
    const CGAL::face_descriptor &fd,
    const CGAL::Vector &face_normal,
    const CGAL::Point3 &face_center,
    const double azimuth,
    const double elevation
) {
    if (elevation < 0.0)
        return true;
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

inline std::vector<bool> shade(const Mesh &mesh, const std::chrono::system_clock::time_point tp) {
    std::vector<bool> shade_vec;
    shade_vec.reserve(mesh.num_faces());
    namespace bg = boost::geometry;
    using point_car = bg::model::point<double, 2, bg::cs::cartesian>;
    using point_geo = bg::model::point<double, 2, bg::cs::geographic<bg::degree>>;
    bg::srs::transformation<> tr{
        bg::srs::proj4(mesh.proj4_str),
        bg::srs::epsg(4326)
    };
    const CGAL::Tree tree(CGAL::faces(mesh.cgal_mesh).first, CGAL::faces(mesh.cgal_mesh).second, mesh.cgal_mesh);
    for (auto fd: mesh.cgal_mesh.faces()) {
        const auto c = centroid(mesh, fd);
        const point_car x_car{c.x(), c.y()};
        point_geo x_geo;
        tr.forward(x_car, x_geo);
        const auto lat = bg::get<1>(x_geo);
        const auto lon = bg::get<0>(x_geo);
        const auto [azimuth, elevation] = solar_position::time_point_solar_position(
                tp,
                lat,
                lon,
                c.z(),
                rasputin::solar_position::collectors::azimuth_and_elevation(),
                rasputin::solar_position::delta_t_calculator::coarse_timestamp_calc()
        );
        auto face_normal = CGAL::Polygon_mesh_processing::compute_face_normal(fd, mesh.cgal_mesh);
        shade_vec.emplace_back(is_shaded(tree, fd, face_normal, c, azimuth, elevation));
    }
    return shade_vec;
}

inline std::vector<int> compute_shadow(const Mesh &mesh, const double azimuth, const double elevation) {
    // Topocentric azimuth and elevation
    const arma::vec::fixed<3> sd = arma::normalise(arma::vec::fixed<3>{sin(azimuth*M_PI/180.0),
                                                                       cos(azimuth*M_PI/180.0),
                                                                       tan(elevation*M_PI/180.0)});
    return compute_shadow(mesh, point3{sd[0], sd[1], sd[2]});
};

inline std::vector<std::vector<int>> compute_shadows(
    const Mesh &mesh, const std::vector<std::pair<int, point3>> & sun_rays
) {
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
}