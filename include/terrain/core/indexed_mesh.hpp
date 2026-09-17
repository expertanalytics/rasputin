#pragma once

// The flat indexed mesh: what a CDT returns and what `mesh` will later build
// its ternary tree from.
//
// It lives in core/, not in cdt/ and not in mesh/. The CDT *returns* one and
// `mesh` *consumes* one, so siting it under either module would make that
// module's vocabulary a dependency of the other; sited here, both depend
// downward and neither depends on the other. Same argument, same shape, as
// Pslg's.
//
// Structure-of-arrays: triangles_ and constrained_edges_ are parallel and
// index-aligned, which is what lets a refinement pass stream the mask without
// touching coordinates.
//
// THE CONSTRUCTOR IS PUBLIC, and that is a deliberate break from Pslg. Pslg's
// is private because it has exactly one legitimate producer and holding one is
// a proof that validation ran. IndexedMesh2 has plural producers on purpose --
// DetriaBackend today, the seam suite's FakeCdtBackend, mesh's flatten step
// later, an eventual ingress handing Python an existing TIN. The decisive one
// is the second: a fake backend that cannot construct a mesh cannot prove the
// backend seam is real.
//
// What a constructed mesh guarantees:
//   1. triangles().size() == constrained_edges().size().
//   2. Every index in triangles() is < vertices().size().
//   3. Every mask is < 8.
//   4. vertices() BEGINS WITH the producing PSLG's vertex buffer, element-wise
//      and in order, so index k means the same point coming out as going in --
//      the same guarantee, and the same reason, as Pslg guarantee 9.
//   5. Immutable after construction; a const IndexedMesh2 is safe for
//      concurrent read from any number of threads, which refinement requires.
//
// 1-3 are DEBUG ASSERTS, not release checks and not a status channel: the mesh
// is produced in bulk by code inside this repo, so a violation is a programmer
// error, and the asan+ubsan Debug job runs every assert on every fixture. The
// O(3T) scan is not worth paying in release on the one path where throughput
// matters.
//
// What it explicitly does NOT guarantee: that every vertex is referenced by a
// triangle (a vertex outside the outer ring or inside a hole is referenced by
// no in-domain triangle), that the mesh is Delaunay, conforming or non-empty
// (those are claims about a particular producer), or ANY ADJACENCY. There is no
// neighbour array, no half-edge and no vertex-to-triangle map: `mesh` owns
// topology and this is the flat handoff.
//
// Like Pslg, IndexedMesh2 stores NO SPAN INTO ITSELF. Every view is computed on
// demand, which is what makes the defaulted copy correct.

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace terrain {

using TriangleIndices = std::array<std::uint32_t, 3>;

// Bit e of a triangle's mask is set iff the edge (v[e], v[(e + 1) % 3]) is a
// constraint edge of the input. Bit 0 is v0->v1, bit 1 is v1->v2, bit 2 is
// v2->v0 -- and NOT the "edge e is opposite vertex e" convention CGAL uses,
// which is a rotation of this one. See terrain/cdt/constrained_edges.hpp.
inline constexpr std::uint8_t kEdgeMask01 = 1u << 0;
inline constexpr std::uint8_t kEdgeMask12 = 1u << 1;
inline constexpr std::uint8_t kEdgeMask20 = 1u << 2;

class IndexedMesh2 {
public:
    IndexedMesh2() = default;

    IndexedMesh2(std::vector<Point2> vertices, std::vector<TriangleIndices> triangles,
                 std::vector<std::uint8_t> constrained_edges)
        : vertices_{std::move(vertices)},
          triangles_{std::move(triangles)},
          constrained_edges_{std::move(constrained_edges)} {
        assert(triangles_.size() == constrained_edges_.size());
#ifndef NDEBUG
        for (const TriangleIndices& t : triangles_) {
            for (const std::uint32_t i : t) {
                assert(static_cast<std::size_t>(i) < vertices_.size());
            }
        }
        for (const std::uint8_t mask : constrained_edges_) {
            assert(mask < 8u);
        }
#endif
    }

    [[nodiscard]] std::span<const Point2> vertices() const noexcept { return vertices_; }
    [[nodiscard]] std::span<const TriangleIndices> triangles() const noexcept {
        return triangles_;
    }
    [[nodiscard]] std::span<const std::uint8_t> constrained_edges() const noexcept {
        return constrained_edges_;
    }

    [[nodiscard]] std::size_t triangle_count() const noexcept { return triangles_.size(); }

    // Asks about TRIANGLES and not about vertices, and the distinction is the
    // whole "forgot addOutline" failure mode: a point set with no outline
    // triangulates successfully and yields zero interior triangles over a full
    // vertex array.
    [[nodiscard]] bool empty() const noexcept { return triangles_.empty(); }

    // The segment from v[e] to v[(e + 1) % 3] -- the same edge bit e of the
    // mask is about, so the accessor and the convention cannot drift apart.
    // Precondition: t < triangle_count(), e < 3.
    [[nodiscard]] Segment2 edge(std::size_t t, std::size_t e) const noexcept {
        assert(t < triangles_.size());
        assert(e < 3);
        const TriangleIndices& tri = triangles_[t];
        return Segment2{vertices_[tri[e]], vertices_[tri[(e + 1) % 3]]};
    }

    // Precondition: t < triangle_count(), e < 3.
    [[nodiscard]] bool is_constrained(std::size_t t, std::size_t e) const noexcept {
        assert(t < constrained_edges_.size());
        assert(e < 3);
        return (constrained_edges_[t] & static_cast<std::uint8_t>(1u << e)) != 0u;
    }

private:
    std::vector<Point2> vertices_;
    std::vector<TriangleIndices> triangles_;
    std::vector<std::uint8_t> constrained_edges_;
};

}  // namespace terrain
