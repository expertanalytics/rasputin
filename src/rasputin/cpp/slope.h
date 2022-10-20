#include <armadillo>

#include <pybind11/numpy.h>

#include "types.h"

namespace rasputin {
inline point3_vector orient_tin(const point3_vector &pts, face_vector &faces) {
    point3_vector result;
    result.reserve(faces.size());
    for (auto& face: faces) {
        // Compute ccw normal
        const auto p0 = pts[face[0]];
        const auto p1 = pts[face[1]];
        const auto p2 = pts[face[2]];
        const arma::vec::fixed<3> v0{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        const arma::vec::fixed<3> n = arma::cross(v0, v1);
        double c = arma::norm(n);

        // Reverse triangle orientation if it is negatively oriented relative to xy plane
        if (n[2] < 0.0) {
          c *= -1.0;
          std::reverse(face.begin(), face.end());
        }

        // Store normalised and correctly oriented normal vector
        result.push_back(point3{n[0]/c, n[1]/c, n[2]/c});
    }
    return result;
};

inline double compute_slope(const point3 &normal) {
    return std::atan2(pow(pow(normal[0], 2) + pow(normal[1], 2), 0.5), normal[2]);
}

inline double_vector compute_slopes(const point3_vector &normals) {
    double_vector result;
    result.reserve(normals.size());

    for (const auto &n : normals)
        result.push_back(compute_slope(n));

    return result;
};

inline double compute_aspect(const point3 &normal) {
    return std::atan2(normal[0], normal[1]);
}

inline double_vector compute_aspects(const point3_vector &normals) {
    double_vector result;
    result.reserve(normals.size());
    for (const auto &n: normals)
        result.emplace_back(compute_aspect(n));
    return result;
};

inline point3 normal(const point3 &p0, const point3 &p1, const point3 &p2) {
    const arma::vec::fixed<3> v0{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    const arma::vec::fixed<3> v1{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    arma::vec::fixed<3> n = arma::cross(v0, v1);
    n /= arma::norm(n);
    return n[2] >= 0.0 ? point3{n[0], n[1], n[2]} : point3{-n[0], -n[1], -n[2]};
}

inline point3_vector surface_normals(const point3_vector &pts, const face_vector &faces) {
    point3_vector result;
    result.reserve(faces.size());
    for (const auto face: faces)
        result.emplace_back(normal(pts[face[0]], pts[face[1]], pts[face[2]]));
    return result;
};

template <typename T> void iadd(T& v, const T& o) {
    v[0] += o[0];
    v[1] += o[1];
    v[2] += o[2];
}

inline point3_vector point_normals(const point3_vector &pts, const face_vector &faces) {
    point3_vector result(pts.size(), {0.0, 0.0, 0.0});
    for (auto face: faces) {
        const auto p0 = pts[face[0]];
        const auto p1 = pts[face[1]];
        const auto p2 = pts[face[2]];
        const arma::vec::fixed<3> v0{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
        const arma::vec::fixed<3> v1{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        arma::vec::fixed<3> n = arma::cross(v0, v1);
        n /= arma::norm(n);
        const auto v = (n[2] >= 0.0) ? point3{n[0], n[1], n[2]} : point3{-n[0], -n[1], -n[2]};
        iadd(result[face[0]], v);
        iadd(result[face[1]], v);
        iadd(result[face[2]], v);
    }
    for (int i = 0; i < result.size(); ++i) {
        point3 &p = result[i];
        const double norm = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        if (norm > 1.0e-16) {
            p[0] /= norm;
            p[1] /= norm;
            p[2] /= norm;
        }
    }
    return result;
}

template <typename CB>
std::tuple<face_vector, face_vector> partition(const point3_vector &pts,
                                               const face_vector &faces,
                                               CB criterion) {
    face_vector part1, part2;
    for (auto face: faces) {
        if (criterion(pts[face[0]], pts[face[1]], pts[face[2]]))
            part1.emplace_back(face);
        else
            part2.emplace_back(face);
    }
    return std::make_pair(std::move(part1), std::move(part2));
}

inline std::tuple<face_vector, face_vector> extract_lakes(const point3_vector&pts, const face_vector &faces) {
   return partition(pts, faces, [] (const point3 &p0, const point3 &p1, const point3 &p2){
       return compute_slope(normal(p0, p1, p2)) < 1.0e-2;
   });
}

inline std::tuple<face_vector, face_vector> extract_avalanche_expositions(
    const point3_vector &pts,
    const face_vector &faces,
    const point2_vector &exposed_intervals,
    const point2_vector &height_intervals
) {
    return partition(pts, faces, [exposed_intervals, height_intervals](const point3 &p0, const point3 &p1, const point3 &p2){
        const auto max_height = std::max(p0[2], std::max(p1[2], p2[2]));
        const auto min_height = std::min(p0[2], std::min(p1[2], p2[2]));
        bool inside = false;
        for (auto height_interval: height_intervals) {
            if ((max_height <= height_interval[1] && max_height >= height_interval[0]) ||
                (min_height <= height_interval[1] && min_height >= height_interval[0]))
                inside = true;
        }
        if (not inside)
            return false;
        const auto cell_normal = normal(p0, p1, p2);
        const auto cell_slope = compute_slope(cell_normal);
        if (cell_slope < 30./180.*M_PI)
            return false;
        const auto aspect = compute_aspect(cell_normal);
        for (auto exposition: exposed_intervals) {
            if ((exposition[0] < aspect) && (aspect < exposition[1]))
                return true;
            else if ((exposition[0] > exposition[1]) && ((exposition[0] < aspect) || (aspect < exposition[1])))
                return true;
        }
        return false;
    });
}

inline point3_vector cell_centers(const point3_vector& points, const face_vector & faces) {
    point3_vector result;
    result.reserve(faces.size());
    for (auto f: faces) {
        auto p0 = points[f[0]];
        auto p1 = points[f[1]];
        auto p2 = points[f[2]];
        auto x = (p0[0] + p1[0] + p2[0])/3.0;
        auto y = (p0[1] + p1[1] + p2[1])/3.0;
        auto z = (p0[2] + p1[2] + p2[2])/3.0;
        result.emplace_back(point3{x, y, z});
    }
    return result;
}

inline index_vector coordinates_to_indices(
    double x0,
    double y1,
    double dx,
    double dy,
    unsigned int M,
    unsigned int N,
    point2_vector pts
) {

    index_vector indices;
    indices.reserve(pts.size());
    for (auto pt: pts)
        indices.emplace_back(std::array<unsigned int, 2>{(unsigned int)((pt[0]-x0)/dx), (unsigned int)((y1 - pt[1])/dy)});
    return indices;
}

template <typename T>
std::vector<T> extract_buffer_values(const index_vector& indices, pybind11::array_t<T>& array) {
    std::vector<T> result;
    result.reserve(indices.size());
    auto buffer = array.request();
    unsigned long M = (unsigned long)buffer.shape[0];
    unsigned long N = (unsigned long)buffer.shape[1];
    T* ptr = (T *)buffer.ptr;
    for (auto idx: indices)
        result.emplace_back(ptr[idx[0]*N + idx[1]]);  // TODO: Implement range check?
    return result;
}

inline std::tuple<point3_vector, face_vector> consolidate(const point3_vector &points, const face_vector &faces) {
    face_vector new_faces;
    point3_vector new_points;
    new_faces.reserve(faces.size());
    std::map<int, int> point_map;
    int n = 0;
    for (const auto _face: faces) {
        face new_face;
        int i = 0;
        for (const auto f: _face) {
            if (not point_map.count(f)) {
                point_map.insert(std::make_pair(f, n++));
                new_points.emplace_back(points[f]);
            }
            new_face[i++] = point_map[f];
        }
        new_faces.emplace_back(new_face);
    }
    return std::make_pair(std::move(new_points), std::move(new_faces));
}
}
