#include "mesh.h"

namespace rasputin {
struct Mesh: public traits::MeshBase<CGAL::Mesh> {
    Mesh(CGAL::Mesh cgal_mesh, const std::string proj4_str) : traits::MeshBase<CGAL::Mesh>(cgal_mesh, proj4_str) { }

    template<typename S, typename P, typename C>
    Mesh coarsen(const S& stop, const P& placement, const C& cost) const {
        CGAL::Mesh new_cgal_mesh = CGAL::Mesh(this->mesh);
        CGAL::Surface_mesh_simplification::edge_collapse(new_cgal_mesh,
                                                         stop,
                                                         CGAL::parameters::get_cost(cost)
                                                                          .get_placement(placement));
        return Mesh(new_cgal_mesh, proj4_str);
    }

    size_t num_edges() const {
        return mesh.number_of_edges();
    }
    size_t num_vertices() const {
        return mesh.number_of_vertices();
    }
    size_t num_faces() const {
        return mesh.number_of_faces();
    }

    CGAL::Mesh construct_mesh(
        const point3_vector &pts, const face_vector &faces, VertexIndexMap& index_map, FaceDescrMap& face_map
    ) {
        CGAL::Mesh mesh;
        index_map.clear();
        face_map.clear();
        size_t i = 0;
        size_t j = 0;
        for (auto p: pts)
            index_map[i++] = mesh.add_vertex(CGAL::Point(p[0], p[1], p[2]));
        for (auto f: faces)
            face_map[mesh.add_face(index_map[f[0]], index_map[f[1]], index_map[f[2]])] = j++;
        return mesh;
    };

    private:
    point3_vector points;
    face_vector faces;

    void set_points_faces() {
        points.clear();
        faces.clear();

        points.reserve(this->mesh.num_vertices());
        faces.reserve(this->num_faces());

        int n = 0;
        std::map<CGAL::VertexIndex, int> reindex;
        for (auto f: this->mesh.faces()) {
            std::array<int, 3> fl;
            size_t idx = 0;
            for (auto v: this->mesh.vertices_around_face(this->mesh.halfedge(f))) {
                if (reindex.count(v) == 0) {
                    reindex.emplace(v, n++);
                    const auto pt = this->mesh.point(v);
                    points.emplace_back(point3{pt.x(), pt.y(), pt.z()});
                }
                fl[idx++] = reindex[v];
            }
            faces.emplace_back(face{fl[0], fl[1], fl[2]});
        }
    }
};

namespace traits {
template<>
struct Delaunay<RasterData, CGAL::DelaunayConstraints> {
    template<typename T, typename P>
    CGAL::DelaunayConstraints interpolate_boundary_points(const RasterData<T>& raster, const P& boundary_polygon) {
        // First we need to determine intersection points between the raster domain and the polygon
        CGAL::MultiPolygon intersection_polygon = raster.compute_intersection(boundary_polygon);

        // Iterate over edges of the intersection polygon and interpolate points
        // TODO: Handle holes. Only needed if the intersecting polygon has holes.
        CGAL::DelaunayConstraints interpolated_points;
        for (const auto& part : CGAL::extract_boundaries(intersection_polygon)) {
            for (auto e = part.edges_begin(); e != part.edges_end(); ++e) {
                CGAL::Point2 first_vertex = e->vertex(0);
                CGAL::Point2 second_vertex = e->vertex(1);

                // We can skip edges that are aligend with the raster boundary
                // TODO: Consider checking using CGAL and exact arithmentic
                bool edge_is_aligned = not raster.contains((first_vertex.x() + second_vertex.x())/2,
                                                           (first_vertex.y() + second_vertex.y())/2);
                if (edge_is_aligned) {
                    continue;
                }

                // Sample with the approximately the same resolution as the raster data along the boundary edges
                double edge_len_x = second_vertex.x() - first_vertex.x();
                double edge_len_y = second_vertex.y() - first_vertex.y();
                std::size_t num_subedges = static_cast<int>(std::max<double>(std::fabs(edge_len_x/raster.delta_x),
                                                                             std::fabs(edge_len_y/raster.delta_y)));
                num_subedges = std::max<int>(1, num_subedges);

                double edge_dx = edge_len_x / num_subedges; // signed distance
                double edge_dy = edge_len_y / num_subedges;

                // Iterate over edge samples
                std::vector<CGAL::Point> interpolated_points_on_edge;
                for (std::size_t k=0; k < num_subedges + 1; ++k) {
                    double x = first_vertex.x() + k * edge_dx;
                    double y = first_vertex.y() + k * edge_dy;
                    double z = raster.get_interpolated_value_at_point(x, y);
                    interpolated_points_on_edge.emplace_back(x, y, z);
                }
                interpolated_points.push_back(std::move(interpolated_points_on_edge));
            }
        }
        return interpolated_points;
    }
};
} // namespace traits

} // namespace rasputin
