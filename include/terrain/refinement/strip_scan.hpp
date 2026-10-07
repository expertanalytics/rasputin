#pragma once

// The edge strip's machinery inside the refinement loop
// (docs/increments/15f-edge-strip.md, D4, L1-L5, L10, L12, L13, L16): the
// per-triangle scan result shared by every point set, the constrained sub-edge
// map, the strip scan by owned sub-edge, the cut with consumption, the guarded
// insertion, L16's coincidence radius (coincidence_radius, which L12 and L14
// share), and L12's test for a point a hair off a constrained edge.
// refine_points.hpp's point_loop is its one user; constraint_points.hpp (the
// generator, which 23b reuses without the loop) does not depend on it.

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <tuple>
#include <unordered_map>

namespace terrain::refinement::detail {

enum class PointSet : std::uint8_t { Source, Strip, Dem };

// One triangle's scan: the winner across sets (for the split), and the
// non-strip set's own largest error and every set's uncovered count (L4).
struct PointScan {
    double max_error = 0.0;                  // the non-strip set's own; 0 when void
    std::optional<mesh::MeshVertex> point;  // the winner, or a void carve point
    double z = 0.0;                          // its z
    double error = 0.0;                      // its error
    NodeLocation where = NodeLocation::Inside;
    bool is_void = false;       // the winner is a carve point (source: the triangle is void)
    std::size_t uncovered = 0;  // points of every set in void triangles or on void sub-edges
    PointSet set = PointSet::Source;
    std::size_t strip_index = 0;  // a strip winner's flag index
    double s = 0.0;               // a strip winner's parameter
    std::size_t on_frozen = 0;    // stored points on a frozen edge (23b, N6), each in one triangle
    double frozen_error = 0.0;    // their largest error against the edge

    // Offered in set order: the first void candidate wins, else the strictly larger error.
    void offer(const PointScan& c) {
        if ((point && is_void) || !(c.is_void || !point || c.error > error))
            return;
        const auto own = std::tuple{max_error, uncovered, on_frozen, frozen_error};
        *this = c;
        std::tie(max_error, uncovered, on_frozen, frozen_error) = own;
    }
};

// A constrained sub-edge of strip edge k, ends a and b in the order of s
// (s_a < s_b, a nearer P0), keyed by its unordered vertex pair (D4 step 1, L5).
struct SubEdge {
    std::size_t k;
    std::uint32_t a, b;
    double s_a, s_b;
};
using SubEdges = std::unordered_map<std::uint64_t, SubEdge>;

[[nodiscard]] inline std::uint64_t edge_key(std::uint32_t u, std::uint32_t v) {
    return std::uint64_t{std::min(u, v)} << 32 | std::max(u, v);
}

// The mesh value along a sub-edge at s: z_a + sigma (z_b - z_a), and an end's
// own z at its own s (step 6, L5). The scan and step 6 share it, so they agree.
[[nodiscard]] inline double along(const SubEdge& e, std::span<const double> zt, double s) {
    if (s == e.s_a || s == e.s_b)
        return zt[s == e.s_a ? e.a : e.b];
    return zt[e.a] + (s - e.s_a) / (e.s_b - e.s_a) * (zt[e.b] - zt[e.a]);
}

// The strip's candidate in t, offered to r: over the sub-edges t owns, the
// unrefused points strictly inside, the worst by strictly larger error (lower
// edge index, then lower s, on ties); on a void sub-edge the point nearest an
// invalid end, its points counted in `uncovered` (D4 step 2).
inline void scan_strip(const ConstraintCheckPoints& strip, std::span<const std::size_t> offset,
                       std::span<const char> refused, const SubEdges& subs, const mesh::LatticeMesh& m,
                       std::span<const double> zt, std::uint32_t t, PointScan& r) {
    const auto& tri = m.triangles()[t];
    PointScan best;
    double nearest = std::numeric_limits<double>::infinity();
    for (unsigned e = 0; e < 3; ++e) {
        const std::uint32_t u = tri[e], v = tri[(e + 1) % 3];
        if (!m.is_constrained(t, e) || (u > v && m.neighbours(t)[e] != mesh::kNoNeighbour))
            continue;
        const auto it = subs.find(edge_key(u, v));
        if (it == subs.end())
            continue;
        const SubEdge& se = it->second;
        const auto pts = strip.on_edge(se.k);
        const bool void_a = std::isnan(zt[se.a]), void_b = std::isnan(zt[se.b]);
        auto j = static_cast<std::size_t>(
            std::partition_point(pts.begin(), pts.end(), [&](const ConstraintPoint& p) { return p.s <= se.s_a; })
            - pts.begin());
        for (; j < pts.size() && pts[j].s < se.s_b; ++j) {
            const ConstraintPoint& p = pts[j];
            if (refused[offset[se.k] + j] != 0)
                continue;
            PointScan c;
            c.set = PointSet::Strip;
            c.point = p.at;
            c.z = p.z;
            c.where = static_cast<NodeLocation>(e + 1);
            c.strip_index = offset[se.k] + j;
            c.s = p.s;
            if (void_a || void_b) {
                ++r.uncovered;
                const double inf = std::numeric_limits<double>::infinity();
                const double d = std::min(void_a ? p.s - se.s_a : inf, void_b ? se.s_b - p.s : inf);
                c.is_void = true;
                if (d < nearest) {
                    nearest = d;
                    best = c;
                }
            } else if (c.error = std::abs(p.z - along(se, zt, p.s)); !best.is_void && c.error > best.error) {
                best = c;
            }
        }
    }
    if (best.point)
        r.offer(best);
}

// D4 step 4: a split at vertex q of the constrained edge (u, v) replaces its
// sub-edge record by two. s_q is a strip point's own s; for any other point,
// that of the strip point at exactly q's position (consumed, L5), else q's
// projection onto the sub-edge mapped into [s_a, s_b].
inline void cut(SubEdges& subs, const ConstraintCheckPoints& strip, const mesh::LatticeMesh& m,
                std::uint32_t u, std::uint32_t v, std::uint32_t q, std::optional<double> s) {
    const auto it = subs.find(edge_key(u, v));
    if (it == subs.end())
        return;
    const SubEdge se = it->second;
    const mesh::MeshVertex p = m.vertices()[q], a = m.vertices()[se.a], b = m.vertices()[se.b];
    for (const ConstraintPoint& c : strip.on_edge(se.k))
        if (!s && c.s > se.s_a && c.s < se.s_b && c.at == p)
            s = c.s;
    if (!s) {
        const double dc = b.col - a.col, dr = b.row - a.row;
        const double sigma =
            std::clamp(((p.col - a.col) * dc + (p.row - a.row) * dr) / (dc * dc + dr * dr), 0.0, 1.0);
        s = se.s_a + sigma * (se.s_b - se.s_a);
    }
    subs.erase(it);
    subs[edge_key(se.a, q)] = SubEdge{se.k, se.a, q, se.s_a, *s};
    subs[edge_key(q, se.b)] = SubEdge{se.k, q, se.b, *s, se.s_b};
}

// foot_fits, and the split's new edges q-c (and q-d) locally Delaunay. A point
// a hair off its sub-edge near an end can lie outside t's circle, and
// legalise_around tests only the edges opposite q.
inline bool strip_fits(const mesh::LatticeMesh& m, std::uint32_t t, unsigned e, mesh::MeshVertex f,
                       const mesh::LatticeFrame& fr) {
    using K = pred::DefaultKernel;
    if (!mesh::detail::foot_fits(m, t, e, f))
        return false;
    // Whether edge (p, q) of the counter-clockwise (p, q, r) must flip, s across it (lawson.hpp's must_flip).
    const auto flips = [&](mesh::MeshVertex p, mesh::MeshVertex q, mesh::MeshVertex r, mesh::MeshVertex s) {
        const Point2 a = fr.at(p), b = fr.at(q), c = fr.at(r), d = fr.at(s);
        if (K::orient2d(a, b, c) == pred::Orientation::CounterClockwise)
            return K::incircle(a, b, c, d) == pred::Incircle::Inside;
        return K::orient2d(b, a, d) == pred::Orientation::CounterClockwise
            && K::incircle(b, a, d, c) == pred::Incircle::Inside;
    };
    const mesh::MeshVertex a = m.corner(t, e), b = m.corner(t, (e + 1) % 3), c = m.corner(t, (e + 2) % 3);
    if (flips(f, c, a, b))
        return false;
    const std::uint32_t u = m.neighbours(t)[e];
    if (u == mesh::kNoNeighbour)
        return true;
    unsigned k = 0;
    while (m.neighbours(u)[k] != t)
        ++k;
    return !flips(f, m.corner(u, (k + 2) % 3), b, a);
}

// L16: the coincidence radius of L12 and L14, max(1e-10, 64 ulp(M)) lattice
// units, M = max(cols, rows) - 1 the largest lattice coordinate.
[[nodiscard]] inline double coincidence_radius(const raster::RasterGeometry& g) {
    const auto m = static_cast<double>(std::max(g.cols(), g.rows()) - 1);
    return std::max(1e-10, 64.0 * (std::nextafter(m, std::numeric_limits<double>::infinity()) - m));
}

// L12: the constrained edge of t whose line p lies within `radius` lattice
// units of, (col, row) Euclidean, with p's projection strictly inside it; the
// nearest, ties to the lower edge index. Never a frozen edge (23b, N7): the
// scans already skip a point this near one, but by a copy of this expression
// that the compiler may contract differently, so the edge is excluded here too.
inline std::optional<unsigned> near_constraint(const mesh::LatticeMesh& m, std::uint32_t t, mesh::MeshVertex p,
                                               double radius) {
    std::optional<unsigned> best;
    double best_d = radius;
    for (unsigned e = 0; e < 3; ++e) {
        const mesh::MeshVertex a = m.corner(t, e), b = m.corner(t, (e + 1) % 3);
        const double dc = b.col - a.col, dr = b.row - a.row, len2 = dc * dc + dr * dr;
        const double sigma = ((p.col - a.col) * dc + (p.row - a.row) * dr) / len2;
        const double d = std::abs(dc * (p.row - a.row) - dr * (p.col - a.col)) / std::sqrt(len2);
        if (m.is_constrained(t, e) && !m.is_frozen(t, e) && sigma > 0.0 && sigma < 1.0 && (d < best_d || (!best && d == best_d))) {
            best = e;
            best_d = d;
        }
    }
    return best;
}

}  // namespace terrain::refinement::detail
