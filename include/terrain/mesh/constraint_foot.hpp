#pragma once

// The foot of a point on a nearby constraint edge
// (docs/increments/20c-soft-quality.md, R1): the one geometric helper the
// quality start, refinement and the final check ask before they insert a
// point near a constraint. Each path answers the status in its own way.
//
// Search order, fixed so the result is a function of the mesh and p only:
// t's edges in edge order, then, for each of t's unconstrained edges in edge
// order, the neighbour's other two edges. Never across a constrained edge
// (frozen or not), and never onto a frozen edge (23b, N5). "Closer than
// delta" is strict, measured in the world frame f to the closed segment; an
// edge p lies exactly on is skipped. The first edge found decides.
//
// Depends on core, predicates, lattice_mesh.hpp and lawson.hpp (for
// LatticeFrame); knows no raster and reads no height.

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>

namespace terrain::mesh {

enum class FootStatus : std::uint8_t { None, Hit, NearEnd, NotCounterClockwise };

struct FootSearch {
    FootStatus status = FootStatus::None;
    std::uint32_t owner = 0;  // the triangle whose edge it is (t or a neighbour)
    unsigned edge = 0;
    MeshVertex at{};  // the foot; meaningful for Hit only
};

namespace detail {

// Every child of splitting t's edge e (and its neighbour's) at f is strictly
// counter-clockwise (20b R2 step 4).
[[nodiscard]] inline bool foot_fits(const LatticeMesh& m, std::uint32_t t, unsigned e, MeshVertex f) {
    const MeshVertex a = m.corner(t, e), b = m.corner(t, (e + 1) % 3), c = m.corner(t, (e + 2) % 3);
    if (orient_sign(a, f, c) <= 0 || orient_sign(f, b, c) <= 0)
        return false;
    const std::uint32_t u = m.neighbours(t)[e];
    if (u == kNoNeighbour)
        return true;
    unsigned k = 0;
    while (m.neighbours(u)[k] != t)
        ++k;
    const MeshVertex d = m.corner(u, (k + 2) % 3);
    return orient_sign(b, f, d) > 0 && orient_sign(f, a, d) > 0;
}

// The verdict of edge e of triangle o, or nothing when it is not a candidate.
[[nodiscard]] inline std::optional<FootSearch> foot_on(const LatticeMesh& m, std::uint32_t o, unsigned e,
                                                       MeshVertex p, double delta, const LatticeFrame& f) {
    const MeshVertex a = m.corner(o, e), b = m.corner(o, (e + 1) % 3);
    if (!m.is_constrained(o, e) || m.is_frozen(o, e) || orient_sign(a, b, p) == 0)
        return std::nullopt;
    const double ux = (b.col - a.col) * f.dx, uy = (b.row - a.row) * f.dy;
    const double px = (p.col - a.col) * f.dx, py = (p.row - a.row) * f.dy;
    const double s = std::clamp((px * ux + py * uy) / (ux * ux + uy * uy), 0.0, 1.0);
    if (std::hypot(px - s * ux, py - s * uy) >= delta)
        return std::nullopt;
    const double len = std::hypot(ux, uy);
    if (s * len < delta || (1.0 - s) * len < delta)
        return FootSearch{FootStatus::NearEnd, o, e, {}};
    const MeshVertex at{a.col + s * (b.col - a.col), a.row + s * (b.row - a.row)};
    return FootSearch{foot_fits(m, o, e, at) ? FootStatus::Hit : FootStatus::NotCounterClockwise, o, e, at};
}

}  // namespace detail

[[nodiscard]] inline FootSearch constraint_foot(const LatticeMesh& m, std::uint32_t t, MeshVertex p, double delta,
                                                const LatticeFrame& f) {
    for (unsigned e = 0; e < 3; ++e)
        if (const auto r = detail::foot_on(m, t, e, p, delta, f))
            return *r;
    for (unsigned e = 0; e < 3; ++e) {
        const std::uint32_t u = m.neighbours(t)[e];
        if (m.is_constrained(t, e) || u == kNoNeighbour)
            continue;
        unsigned k = 0;
        while (m.neighbours(u)[k] != t)
            ++k;
        for (const unsigned j : {(k + 1) % 3, (k + 2) % 3})
            if (const auto r = detail::foot_on(m, u, j, p, delta, f))
                return *r;
    }
    return {};
}

}  // namespace terrain::mesh
