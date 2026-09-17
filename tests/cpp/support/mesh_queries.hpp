#pragma once

// Test-only queries over an IndexedMesh2: the adjacency the mesh deliberately
// does NOT carry, and the exact visibility predicate the Delaunay property
// needs.
//
// IndexedMesh2 promises no neighbour array, no half-edge and no vertex-to-
// triangle map -- adjacency is `mesh`'s vocabulary and putting it in core/
// would be putting the ternary tree's vocabulary there. So a suite that needs
// to ask "which triangles share this edge" rebuilds it here, in O(T log T),
// which is also the honest way to test a flat handoff: the test derives the
// topology independently instead of trusting a field.
//
// std::map, not std::unordered_map, and for a different reason than
// constrained_edges.hpp gives: the iteration order of a failing property has to
// be reproducible across runs and standard libraries, or a shrunk counter-
// example does not reproduce.
//
// Everything numeric here is a kernel call or a comparison between coordinates
// the mesh already holds. No constructed intersection point, no tolerance, no
// division -- the same rule point_in_ring keeps, and for the same reason: an
// oracle that rounds cannot judge a predicate that does not.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <utility>
#include <vector>

namespace terrain::test {

// Where an undirected edge sits: triangle index and the edge slot e within it,
// so a caller can ask the mask about the same edge it found.
struct EdgeUse {
    std::size_t triangle{};
    std::size_t slot{};
};

[[nodiscard]] inline std::uint64_t undirected_key(std::uint32_t a, std::uint32_t b) noexcept {
    const std::uint32_t lo = std::min(a, b);
    const std::uint32_t hi = std::max(a, b);
    return (static_cast<std::uint64_t>(lo) << 32) | hi;
}

// Every undirected edge of the mesh, mapped to the triangle slots that use it.
// A well-formed CDT gives every edge either one use (a boundary edge) or two.
[[nodiscard]] inline std::map<std::uint64_t, std::vector<EdgeUse>> edge_uses(
    const IndexedMesh2& mesh) {
    std::map<std::uint64_t, std::vector<EdgeUse>> uses;
    const std::span<const TriangleIndices> tris = mesh.triangles();
    for (std::size_t t = 0; t < tris.size(); ++t) {
        for (std::size_t e = 0; e < 3; ++e) {
            const std::uint64_t k = undirected_key(tris[t][e], tris[t][(e + 1) % 3]);
            uses[k].push_back(EdgeUse{t, e});
        }
    }
    return uses;
}

// The vertex of triangle t opposite edge slot e -- the apex the incircle test
// asks about.
[[nodiscard]] inline std::uint32_t apex(const IndexedMesh2& mesh, const EdgeUse& u) noexcept {
    return mesh.triangles()[u.triangle][(u.slot + 2) % 3];
}

// Does the OPEN segment s meet the closed segment c anywhere other than at a
// shared endpoint of s? That is what "blocked" means for visibility: a proper
// crossing blocks, and so does a constraint endpoint sitting on the segment's
// interior, which is the T-junction case.
//
// Exact under an exact kernel: four orientation calls and, for the touching
// arms, on_segment -- no parameter t and no constructed point.
template <pred::GeometryKernel K>
[[nodiscard]] bool blocks_visibility(const Segment2& s, const Segment2& c) {
    using pred::Orientation;

    const Orientation o1 = K::orient2d(s.a, s.b, c.a);
    const Orientation o2 = K::orient2d(s.a, s.b, c.b);
    const Orientation o3 = K::orient2d(c.a, c.b, s.a);
    const Orientation o4 = K::orient2d(c.a, c.b, s.b);

    // A proper crossing: each segment strictly separates the other's endpoints.
    if (o1 != o2 && o3 != o4 && o1 != Orientation::Collinear && o2 != Orientation::Collinear &&
        o3 != Orientation::Collinear && o4 != Orientation::Collinear) {
        return true;
    }
    // A constraint endpoint strictly inside the open segment s.
    const auto strictly_inside = [&](const Point2& p) {
        return p != s.a && p != s.b && on_segment<K>(s, p);
    };
    return strictly_inside(c.a) || strictly_inside(c.b);
}

}  // namespace terrain::test
