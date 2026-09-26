#pragma once

// Lawson legalisation of a LatticeMesh (docs/increments/14b-delaunay-insertion.md,
// R1, R3 to R5).
//
// An edge is flipped only when the neighbour's apex is strictly inside the
// triangle's circle (`Incircle::Inside`), never on a cocircular tie, and never
// when it is constrained or has no neighbour. Each such flip strictly lowers
// the lifted triangulation, so the loops end; that needs the sign to be exact
// for fixed points, which `FilteredKernel<DetriaExact>` gives.
//
// The circle test runs in a LatticeFrame, (col * dx, -(row * dy)) on the
// fractional coordinates (docs/increments/16-domain-polygon.md, R2): world
// coordinates without the translation, so the Delaunay property holds in the
// world when dx != dy and the coordinates stay small for the filter.
//
// Both legalisers use an explicit stack popped last-in first-out, so the
// result depends only on the mesh and the seeds. They return the flip count
// and report every slot a flip writes through on_write (repeats allowed).
//
// Depends on core, predicates and lattice_mesh.hpp; knows no raster.

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/predicates/kernel.hpp>

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace terrain::mesh {

struct LatticeFrame {
    double dx = 1.0;
    double dy = 1.0;

    [[nodiscard]] Point2 at(MeshVertex v) const noexcept {
        return Point2{v.col * dx, -(v.row * dy)};
    }
};

namespace detail {

// Whether t's edge e must flip: interior, unconstrained, and the apex across
// it strictly inside t's circle in the frame.
template <pred::GeometryKernel K>
[[nodiscard]] bool must_flip(const LatticeMesh& m, std::uint32_t t, unsigned e,
                             const LatticeFrame& f) {
    const auto u = m.neighbours(t)[e];
    if (u == kNoNeighbour || m.is_constrained(t, e))
        return false;
    const auto& tri = m.triangles()[t];
    unsigned j = 0;
    while (m.triangles()[u][j] != tri[(e + 1) % 3])
        ++j;
    const auto v = m.vertices();
    return K::incircle(f.at(v[tri[0]]), f.at(v[tri[1]]), f.at(v[tri[2]]),
                       f.at(v[m.triangles()[u][(j + 2) % 3]]))
        == pred::Incircle::Inside;
}

}  // namespace detail

// Legalise around the vertex q after it was inserted. For each seed slot that
// contains q, the edge opposite q is tested; a flip leaves q at index 0 of both
// slots it writes, and both go back on the stack.
template <pred::GeometryKernel K, class OnWrite>
std::size_t legalise_around(LatticeMesh& m, std::uint32_t q, std::span<const std::uint32_t> seeds,
                            const LatticeFrame& f, OnWrite&& on_write) {
    std::vector<std::uint32_t> stack(seeds.begin(), seeds.end());
    std::size_t flips = 0;
    while (!stack.empty()) {
        const auto t = stack.back();
        stack.pop_back();
        const auto& tri = m.triangles()[t];
        unsigned i = 0;
        while (i < 3 && tri[i] != q)
            ++i;
        if (i == 3 || !detail::must_flip<K>(m, t, (i + 1) % 3, f))
            continue;
        const auto u = m.neighbours(t)[(i + 1) % 3];
        m.flip(t, (i + 1) % 3);
        ++flips;
        on_write(t);
        on_write(u);
        stack.push_back(t);
        stack.push_back(u);
    }
    return flips;
}

// Legalise the whole mesh: every interior edge once on the stack, and each
// flip pushes the four outer sides of its quad. A stale entry (its slot was
// rewritten since) names some current edge and is simply tested.
template <pred::GeometryKernel K, class OnWrite>
std::size_t legalise_all(LatticeMesh& m, const LatticeFrame& f, OnWrite&& on_write) {
    std::vector<std::pair<std::uint32_t, unsigned>> stack;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto u = m.neighbours(t)[k]; u != kNoNeighbour && u > t)
                stack.emplace_back(t, k);
    std::size_t flips = 0;
    while (!stack.empty()) {
        const auto [t, e] = stack.back();
        stack.pop_back();
        if (!detail::must_flip<K>(m, t, e, f))
            continue;
        const auto u = m.neighbours(t)[e];
        m.flip(t, e);  // t = (c, a, d), u = (c, d, b)
        ++flips;
        on_write(t);
        on_write(u);
        stack.insert(stack.end(), {{t, 0u}, {t, 1u}, {u, 1u}, {u, 2u}});
    }
    return flips;
}

}  // namespace terrain::mesh
