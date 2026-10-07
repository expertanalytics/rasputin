#pragma once

// Test-only support for increment 20c, PR 20c-2
// (docs/increments/20c-soft-quality.md, R7, R8 and "Tests @tester writes red
// first", 20c-2): hand-built LatticeMesh fixtures, a recorder that keeps the
// mesh before and after every insertion the quality start makes, and the
// oracles T-P1 and T-P2 compare against. No Catch2 include; the checks that
// assert live in the test files.
//
// What the oracles read, and why it is not the producer's records:
//   - an insertion is seen as two snapshots of the mesh, taken from inside
//     the validity callable improve is given (it is asked before every
//     insertion, about the node, and about the foot it inserts instead), and
//     once after improve returns. The triangles legalise_around wrote are the
//     difference of the two triangle sets, as sets of corner positions;
//   - the angles are this file's own: each corner from its two edge vectors,
//     atan2(|cross|, dot), in degrees, in the frame the producer uses
//     (col * dx, -(row * dy)).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <span>
#include <utility>
#include <vector>

namespace quality_gain {

using terrain::TriangleIndices;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::MeshVertex;
using terrain::mesh::orient_sign;

using Pair = std::pair<std::uint32_t, std::uint32_t>;

// A hand-built mesh in (col, row): triangles counter-clockwise on (col, -row),
// constraint edges as undirected (min, max) vertex pairs with their masks, and
// the node rectangle improve is given.
struct Fixture {
    std::vector<MeshVertex> vertices;
    std::vector<TriangleIndices> triangles;
    std::map<Pair, std::uint32_t> constraints;
    std::size_t rows = 0;
    std::size_t cols = 0;
};

[[nodiscard]] inline std::optional<LatticeMesh> build(const Fixture& f) {
    std::vector<std::uint8_t> bits(f.triangles.size(), 0);
    std::vector<std::array<std::uint32_t, 3>> masks(f.triangles.size(), {0, 0, 0});
    for (std::size_t t = 0; t < f.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto it = f.constraints.find(std::minmax(f.triangles[t][k], f.triangles[t][(k + 1) % 3]));
                it != f.constraints.end()) {
                bits[t] |= static_cast<std::uint8_t>(1u << k);
                masks[t][k] = it->second;
            }
    return LatticeMesh::build(f.vertices, f.triangles, std::move(bits), std::move(masks));
}

// LS1's fixture, shared with the refine scenes: a long constraint line e from
// E1 (2, 10) to E2 (30, 10.5), mask 2, in a quadrilateral outline E1, T
// (16, 0), E2, D (16, 22), mask 1, with P (10, 8), Q (20, 8) and R (13, 7)
// above e, inserted by Lawson (the arrays as that printed them). The worst
// triangle is P Q R, obtuse at R: its circumcentre (15, 18) is a node beyond
// e, so the walk to it stops at e. Its foot on e is not e's midpoint
// (16, 10.25).
[[nodiscard]] inline Fixture line_beyond() {
    Fixture f;
    f.vertices = {{2, 10}, {30, 10.5}, {16, 0}, {16, 22}, {10, 8}, {20, 8}, {13, 7}};
    f.triangles = {{5, 0, 1}, {1, 0, 3}, {1, 2, 5}, {2, 0, 4}, {2, 4, 6}, {5, 4, 0}, {4, 5, 6}, {5, 2, 6}};
    f.constraints = {{{0u, 1u}, 2u}, {{0u, 2u}, 1u}, {{1u, 2u}, 1u}, {{1u, 3u}, 1u}, {{0u, 3u}, 1u}};
    f.rows = 25;
    f.cols = 33;
    return f;
}

// A triangle by its corner positions, rotated so the smallest (col, row)
// leads; the orientation is kept, so a clockwise triple never equals a
// counter-clockwise one.
using Corner = std::pair<double, double>;
using Tri = std::array<Corner, 3>;

[[nodiscard]] inline Tri tri(MeshVertex a, MeshVertex b, MeshVertex c) {
    const std::array<Corner, 3> x{{{a.col, a.row}, {b.col, b.row}, {c.col, c.row}}};
    const auto i = static_cast<std::size_t>(std::min_element(x.begin(), x.end()) - x.begin());
    return {x[i], x[(i + 1) % 3], x[(i + 2) % 3]};
}

[[nodiscard]] inline Tri tri_of(const LatticeMesh& m, std::uint32_t t) {
    return tri(m.corner(t, 0), m.corner(t, 1), m.corner(t, 2));
}

[[nodiscard]] inline std::set<Tri> triangle_set(const LatticeMesh& m) {
    std::set<Tri> s;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) s.insert(tri_of(m, t));
    return s;
}

// What an insertion did: the triangles before that are gone, and the ones
// after that are new.
struct Diff {
    std::set<Tri> removed;
    std::set<Tri> created;
};

[[nodiscard]] inline Diff diff(const LatticeMesh& before, const LatticeMesh& after) {
    const auto b = triangle_set(before), a = triangle_set(after);
    Diff d;
    std::set_difference(b.begin(), b.end(), a.begin(), a.end(), std::inserter(d.removed, d.removed.end()));
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(d.created, d.created.end()));
    return d;
}

// The smallest interior angle, degrees, in frame f (see the header).
[[nodiscard]] inline double min_angle_deg(const Tri& t, const LatticeFrame& f) {
    double best = 180.0;
    for (std::size_t k = 0; k < 3; ++k) {
        const auto at = [&](std::size_t i) { return f.at(MeshVertex{t[i % 3].first, t[i % 3].second}); };
        const auto o = at(k), u = at(k + 1), v = at(k + 2);
        const double ux = u.x - o.x, uy = u.y - o.y, vx = v.x - o.x, vy = v.y - o.y;
        best = std::min(best, std::atan2(std::abs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 / std::numbers::pi);
    }
    return best;
}

// The worst angle of a set of triangles, capped at theta (R7's old and new).
[[nodiscard]] inline double worst_capped(const std::set<Tri>& ts, const LatticeFrame& f, double theta) {
    double w = theta;
    for (const auto& t : ts) w = std::min(w, min_angle_deg(t, f));
    return w;
}

// a is b moved by one vector, corner for corner up to rotation: every edge
// vector is the same double pair, so any angle formula built from edge
// vectors gives both the same bits.
[[nodiscard]] inline bool translate_of(const Tri& a, const Tri& b) {
    for (std::size_t r = 0; r < 3; ++r) {
        bool same = true;
        for (std::size_t k = 1; k < 3 && same; ++k)
            same = a[k].first - a[0].first == b[(k + r) % 3].first - b[r].first
                && a[k].second - a[0].second == b[(k + r) % 3].second - b[r].second;
        if (same) return true;
    }
    return false;
}

// Where p goes into m: `on` 0..2 is t's edge (split_edge), 3 strictly inside t
// (split_inside).
struct Location {
    std::uint32_t t = 0;
    unsigned on = 3;
};

// p by the exact orientation: strictly inside a triangle, or on exactly one
// of its edges. Empty on a vertex or outside the mesh.
[[nodiscard]] inline std::optional<Location> locate(const LatticeMesh& m, MeshVertex p) {
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        int zeros = 0;
        unsigned on = 3;
        bool inside = true;
        for (unsigned k = 0; k < 3; ++k) {
            const int s = orient_sign(m.corner(t, k), m.corner(t, (k + 1) % 3), p);
            inside = inside && s >= 0;
            if (s == 0) {
                ++zeros;
                on = k;
            }
        }
        if (inside) return zeros > 1 ? std::nullopt : std::optional<Location>{Location{t, zeros == 0 ? 3u : on}};
    }
    return std::nullopt;
}

// The slot and edge of `before` holding the undirected edge (a, b).
[[nodiscard]] inline std::optional<Location> edge_of(const LatticeMesh& m, std::uint32_t a, std::uint32_t b) {
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (std::minmax(m.triangles()[t][k], m.triangles()[t][(k + 1) % 3]) == std::minmax(a, b))
                return Location{t, k};
    return std::nullopt;
}

// improve's insertion, as include/terrain/mesh/quality.hpp has it: split_inside, or split_edge with
// the neighbour across among the seeds, then legalise_around from the slots
// the split wrote. Returns the new vertex's index.
template <terrain::pred::GeometryKernel K>
std::uint32_t replay(LatticeMesh& m, Location at, MeshVertex p, const LatticeFrame& f) {
    const auto before = static_cast<std::uint32_t>(m.triangle_count());
    std::vector<std::uint32_t> seeds{at.t, before, before + 1};
    std::uint32_t q = 0;
    if (at.on == 3) {
        q = m.split_inside(at.t, p);
    } else {
        const auto u = m.neighbours(at.t)[at.on];
        q = m.split_edge(at.t, at.on, p);
        if (u == kNoNeighbour)
            seeds.pop_back();
        else
            seeds.push_back(u);
    }
    terrain::mesh::legalise_around<K>(m, q, std::span<const std::uint32_t>{seeds}, f, [](std::uint32_t) {});
    return q;
}

// One insertion as the recorder saw it.
struct Step {
    LatticeMesh before;
    LatticeMesh after;
    MeshVertex p;  // the vertex it added
    Location at;   // where p went in `before`
};

// Keeps a copy of the mesh whenever its vertex count has changed since the
// last copy. Call `snap` from the validity callable and once after improve.
class Recorder {
public:
    explicit Recorder(const LatticeMesh& m) : m_{&m}, snaps_{m} {}

    void snap() {
        if (m_->vertices().size() != snaps_.back().vertices().size()) snaps_.push_back(*m_);
    }

    [[nodiscard]] std::size_t snapshots() const noexcept { return snaps_.size(); }

    // Every insertion in order. `at` is read from the two meshes: when the
    // new vertex has exactly two constrained edges after, to a and b, and
    // (a, b) is an edge before, the insertion split that edge (a foot is a
    // rounded projection, so it need not lie on the line exactly); otherwise
    // it is the exact location of p in `before`. Empty `at` (a vertex count
    // that grew by more than one, or a p that cannot be located) is
    // reported by `ok`.
    [[nodiscard]] std::vector<Step> steps(bool* ok) const {
        std::vector<Step> out;
        *ok = true;
        for (std::size_t i = 1; i < snaps_.size(); ++i) {
            const LatticeMesh& b = snaps_[i - 1];
            const LatticeMesh& a = snaps_[i];
            if (a.vertices().size() != b.vertices().size() + 1) {
                *ok = false;
                continue;
            }
            const auto q = static_cast<std::uint32_t>(b.vertices().size());
            std::set<std::uint32_t> ends;
            for (std::uint32_t t = 0; t < a.triangle_count(); ++t)
                for (unsigned k = 0; k < 3; ++k) {
                    const auto x = a.triangles()[t][k], y = a.triangles()[t][(k + 1) % 3];
                    if (a.is_constrained(t, k) && (x == q || y == q)) ends.insert(x == q ? y : x);
                }
            std::optional<Location> at;
            if (ends.size() == 2) at = edge_of(b, *ends.begin(), *ends.rbegin());
            if (!at) at = locate(b, a.vertices().back());
            if (!at) {
                *ok = false;
                continue;
            }
            out.push_back(Step{b, a, a.vertices().back(), *at});
        }
        return out;
    }

private:
    const LatticeMesh* m_;
    std::vector<LatticeMesh> snaps_;
};

}  // namespace quality_gain
