#include <chrono>

#include <armadillo>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "types.h"
#include "mesh.h"
#include "solar_position.h"

namespace py = pybind11;

namespace rasputin {
namespace traits {
template<typename Mesh>
struct Shade {
    static_assert("Shade trait not implemented!");
};
} // namespace traits

template<traits::MeshType Mesh>
struct Shade {
    using S = traits::Shade<Mesh>;
    using Mesh_t = Mesh;
    using face_descriptor_t = typename Mesh::face_descriptor_t;

    template<typename Point3>
    static Point3 centroid(const Mesh_t& mesh, const face_descriptor_t& face) {
        return S::centroid(mesh, face);
    }

    static std::vector<int> compute_shadow(const Mesh_t & mesh, const point3 &sun_direction) {
        return S::compute_shadow(mesh, sun_direction);
    }

    template<typename Tree, typename Vector, typename Point3>
    static bool is_shaded(
        const Tree& tree,
        const face_descriptor_t& fd,
        const Vector& face_normal,
        const Point3& face_center,
        const double azimuth,
        const double elevation
    ) {
        // return S::template is_shaded<Tree, Vector, Point3>(tree, fd, face_normal, face_center, azimuth, elevation);
        return S::is_shaded(tree, fd, face_normal, face_center, azimuth, elevation);
    }

    static std::vector<int> compute_shadow(const Mesh_t &mesh, const double azimuth, const double elevation) {
        const arma::vec::fixed<3> sd = arma::normalise(
            arma::vec::fixed<3>{sin(azimuth*M_PI/180.0), cos(azimuth*M_PI/180.0), tan(elevation*M_PI/180.0)}
        );
        return compute_shadow(mesh, point3{sd[0], sd[1], sd[2]});
    }

    static std::vector<std::vector<int>> compute_shadows(
        const Mesh_t& mesh, const std::vector<std::pair<int, point3>>& sun_rays
    ) {
        // TODO: Rewrite totally!
        std::vector<std::vector<int>>  result;
        return result;
        /*
        result.reserve(sun_rays.size());
        CGAL::Mesh cgal_mesh = mesh.mesh;
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
    }
};

template<typename S>
concept ShadeType =
    requires {
        typename S::Mesh_t;
        typename S::face_descriptor_t;
    }
    &&
    traits::MeshType<typename S::Mesh_t>;

template<ShadeType Shade>
void bind_shades(py::module& m) {
    using namespace std::chrono;
    using Mesh = typename Shade::Mesh_t;
    m.def("compute_shadow",
        (std::vector<int> (*)(const Mesh&, const rasputin::point3 &))&Shade::compute_shadow,
        "Compute shadows for given topocentric sun position.")
     .def("compute_shadow",
        (std::vector<int> (*)(const Mesh&, const double, const double))&Shade::compute_shadow,
        "Compute shadows for given azimuth and elevation.")
     // .def("compute_shadows",
     //    (std::vector<std::vector<int>> (*)(const Mesh&, const std::vector<std::pair<int, rasputin::point3>>&))&Shade::compute_shadows,
     //    "Compute shadows for a series of times and ray directions.")
     .def("shade", [] (const Mesh& mesh, const double timestamp) {
         const auto secs = seconds(int(std::round(timestamp)));
         const auto millisecs = milliseconds(int(round(1000*fmod(timestamp, 1))));
         const auto tp = sys_days{1970y / January / 1} + secs + millisecs;
         const auto shade_vec = Shade::shade(mesh, tp);
         py::array_t<bool> result(shade_vec.size());
         auto info = result.request();
         auto *data = static_cast<bool*>(info.ptr);
         std::copy(std::begin(shade_vec), std::end(shade_vec), data);
         return result;
     });
}

template<traits::MeshType Mesh>
void bind_shades(py::module& m) {
    bind_shades<traits::Shade<Mesh>>(m);
}
} // namespace rasputin
