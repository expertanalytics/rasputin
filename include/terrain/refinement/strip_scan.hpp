#pragma once

// The edge strip's machinery inside the refinement loop: the
// per-triangle scan result shared by every point set, the constrained sub-edge
// map, the strip scan by owned sub-edge and the cut with consumption. The
// guarded insertion (strip_fits), the coincidence radius and the test for a
// point a hair off a constrained edge (near_constraint) are in refine.hpp,
// which uses them too. refine_points.hpp's point_loop is its one user;
// constraint_points.hpp (the generator, which the seam pass reuses without
// the loop) does not depend on it.

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
#include <type_traits>
#include <unordered_map>

namespace terrain::refinement::detail {

enum class PointSet : std::uint8_t { Source, Strip, Dem };

// One triangle's scan: the winner across sets (for the split), and the
// non-strip set's own largest error and every set's uncovered count.
// With a slope a second one is the slope's winner; SlopeRank holds
// the winner's error.
struct PointScan : SlopeRank {
    double max_error = 0.0;                  // the non-strip set's own; 0 when void
    std::optional<mesh::MeshVertex> point;  // the winner, or a void carve point
    double z = 0.0;                          // its z
    NodeLocation where = NodeLocation::Inside;
    bool is_void = false;       // the winner is a carve point (source: the triangle is void)
    std::size_t uncovered = 0;  // points of every set in void triangles or on void sub-edges
    PointSet set = PointSet::Source;
    std::size_t strip_index = 0;  // a strip winner's flag index
    double s = 0.0;               // a strip winner's parameter
    std::size_t on_frozen = 0;    // stored points on a frozen edge, each in one triangle
    double frozen_error = 0.0;    // their largest error against the edge

    // Offered in set order: the first void candidate wins, else the strictly larger error.
    void offer(const PointScan& c) {
        if ((point && is_void) || !(c.is_void || !point || c.error > error))
            return;
        take(c);
    }
    // c's winner, keeping this triangle's own figures.
    void take(const PointScan& c) {
        const auto own = std::tuple{max_error, uncovered, on_frozen, frozen_error};
        *this = c;
        std::tie(max_error, uncovered, on_frozen, frozen_error) = own;
    }
    // The slope's winner across sets: `over` is the OR, and the strictly
    // larger (ratio, error) wins, so ties go to the earlier set.
    void offer_steep(const PointScan& c) {
        const bool o = over || c.over;
        if (c.ratio > ratio || (c.ratio == ratio && c.error > error))
            *this = c;
        over = o;
    }
};

// A constrained sub-edge of strip edge k, ends a and b in the order of s
// (s_a < s_b, a nearer P0), keyed by its unordered vertex pair.
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
// own z at its own s. The scan and the final measure of every strip point
// share it, so they agree.
[[nodiscard]] inline double along(const SubEdge& e, std::span<const double> zt, double s) {
    if (s == e.s_a || s == e.s_b)
        return zt[s == e.s_a ? e.a : e.b];
    return zt[e.a] + (s - e.s_a) / (e.s_b - e.s_a) * (zt[e.b] - zt[e.a]);
}

// The strip's candidate in t, offered to r: over the sub-edges t owns, the
// unrefused points strictly inside, the worst by strictly larger error (lower
// edge index, then lower s, on ties); on a void sub-edge the point nearest an
// invalid end, its points counted in `uncovered`. With a slope,
// every valid point is also offered to `steep`, held to its cell's class.
template <class S = NoSlope>
void scan_strip(const ConstraintCheckPoints& strip, std::span<const std::size_t> offset,
                std::span<const char> refused, const SubEdges& subs, const mesh::LatticeMesh& m,
                std::span<const double> zt, std::uint32_t t, PointScan& r, const S* slope = nullptr,
                PointScan* steep = nullptr) {
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
            } else {
                c.error = std::abs(p.z - along(se, zt, p.s));
                if (!best.is_void && c.error > best.error)
                    best = c;
                if constexpr (!std::is_same_v<S, NoSlope>) {
                    const std::uint8_t k = slope->cell_class(p.at);
                    std::tie(c.over, c.ratio) = std::pair{c.error > slope->allowed(k), c.error * slope->weight(k)};
                    steep->offer_steep(c);
                }
            }
        }
    }
    if (best.point)
        r.offer(best);
}

// A split at vertex q of the constrained edge (u, v) replaces its
// sub-edge record by two. s_q is a strip point's own s; for any other point,
// that of the strip point at exactly q's position (consumed), else q's
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

}  // namespace terrain::refinement::detail
