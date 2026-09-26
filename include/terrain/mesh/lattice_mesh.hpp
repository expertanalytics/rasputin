#pragma once

// A triangle mesh over the DEM's lattice, with neighbour links
// (docs/increments/14-adaptive-refinement.md, R3 and R4). Every vertex that
// refinement inserts is a DEM node; a start vertex may lie anywhere in the
// node rectangle (docs/increments/16-domain-polygon.md, R2).
//
// A flat triangle array, not a tree: R3's edge split changes two triangles at
// once, and only adjacency can find the second. A split reuses the parent's
// slot for its first child and appends the rest, so the array stays dense and
// is the output as it stands -- there is no flatten step.
//
// Frame. Vertices are fractional (col, row) lattice coordinates, integers for
// a node. Orientation is taken in the world's handedness, x = col and y = -row
// (rows grow downward), so a triangle that is counter-clockwise in world
// coordinates is counter-clockwise here too. `orient` on nodes is exact
// integer arithmetic (uint32 coordinates, so a product of two differences fits
// int64); `orient_sign` on any vertices is DefaultKernel on (col, -row), exact
// on its inputs and of the same sign for nodes.
//
// Edge k of a triangle runs from vertex k to vertex k+1, as in IndexedMesh2.
// Per edge a triangle stores the neighbour across it (or kNoNeighbour), a
// constrained bit and a 32-bit property mask. A split keeps the parent's bit
// and mask on both halves of a split edge, and gives new interior edges
// neither: constraints gain Steiner points but never move.
//
// Depends on core and predicates. It knows no raster; (col, row) are just
// numbers.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace terrain::mesh {

struct LatticeVertex {
    std::uint32_t row{};
    std::uint32_t col{};

    friend constexpr bool operator==(const LatticeVertex&, const LatticeVertex&) = default;
};

// A mesh vertex in fractional lattice coordinates: col = (x - x_min) / dx and
// row = (y_max - y) / dy, integers exactly for a node. A LatticeVertex converts
// to it exactly. The conversion back is checked: exact for a node, and a throw
// for an off-node vertex, never a truncation.
struct MeshVertex {
    double col{};
    double row{};

    constexpr MeshVertex() = default;
    constexpr MeshVertex(double c, double r) noexcept : col{c}, row{r} {}
    constexpr MeshVertex(LatticeVertex v) noexcept  // NOLINT(google-explicit-constructor)
        : col{static_cast<double>(v.col)}, row{static_cast<double>(v.row)} {}

    [[nodiscard]] bool is_node() const noexcept {
        return col == std::floor(col) && row == std::floor(row);
    }
    operator LatticeVertex() const {  // NOLINT(google-explicit-constructor)
        if (!is_node())
            throw std::domain_error("MeshVertex: an off-node vertex is not a LatticeVertex");
        return LatticeVertex{static_cast<std::uint32_t>(row), static_cast<std::uint32_t>(col)};
    }
    [[nodiscard]] Point2 frame() const noexcept { return Point2{col, -row}; }

    friend constexpr bool operator==(const MeshVertex&, const MeshVertex&) = default;
    friend constexpr bool operator==(const MeshVertex& a, const LatticeVertex& b) noexcept {
        return a == MeshVertex{b};
    }
};

inline constexpr std::uint32_t kNoNeighbour = std::numeric_limits<std::uint32_t>::max();

// Twice the signed area of (a, b, c), positive when counter-clockwise in the
// world frame (x = col, y = -row).
[[nodiscard]] constexpr std::int64_t orient(LatticeVertex a, LatticeVertex b,
                                            LatticeVertex c) noexcept {
    const auto dr1 = std::int64_t{b.row} - a.row, dc1 = std::int64_t{b.col} - a.col;
    const auto dr2 = std::int64_t{c.row} - a.row, dc2 = std::int64_t{c.col} - a.col;
    return dr1 * dc2 - dc1 * dr2;
}

// The sign of orient on any vertices, from the exact kernel on (col, -row).
[[nodiscard]] inline int orient_sign(MeshVertex a, MeshVertex b, MeshVertex c) {
    return static_cast<int>(pred::DefaultKernel::orient2d(a.frame(), b.frame(), c.frame()));
}

class LatticeMesh {
public:
    // Refuses (nullopt) mismatched array lengths, an index out of range, a
    // triangle whose orientation is not positive, and a directed edge used
    // twice (an overlap or a flipped triangle). Adjacency is derived here.
    // The LatticeVertex overload is a template so a braced vertex list always
    // picks the MeshVertex one.
    template <std::same_as<LatticeVertex> V>
    [[nodiscard]] static std::optional<LatticeMesh> build(
        std::vector<V> vertices, std::vector<TriangleIndices> triangles,
        std::vector<std::uint8_t> constrained, std::vector<std::array<std::uint32_t, 3>> masks) {
        return build(std::vector<MeshVertex>(vertices.begin(), vertices.end()),
                     std::move(triangles), std::move(constrained), std::move(masks));
    }
    [[nodiscard]] static std::optional<LatticeMesh> build(
        std::vector<MeshVertex> vertices, std::vector<TriangleIndices> triangles,
        std::vector<std::uint8_t> constrained, std::vector<std::array<std::uint32_t, 3>> masks) {
        const std::size_t n = triangles.size();
        if (constrained.size() != n || masks.size() != n || n >= kNoNeighbour)
            return std::nullopt;
        LatticeMesh m;
        m.vertices_ = std::move(vertices);
        m.triangles_ = std::move(triangles);
        m.constrained_ = std::move(constrained);
        m.masks_ = std::move(masks);
        m.neighbours_.assign(n, {kNoNeighbour, kNoNeighbour, kNoNeighbour});

        std::unordered_map<std::uint64_t, std::uint32_t> directed;  // (from, to) -> triangle
        auto key = [](std::uint32_t a, std::uint32_t b) { return std::uint64_t{a} << 32 | b; };
        for (std::uint32_t t = 0; t < n; ++t) {
            const auto& tri = m.triangles_[t];
            for (const auto v : tri)
                if (v >= m.vertices_.size())
                    return std::nullopt;
            if (orient_sign(m.vertices_[tri[0]], m.vertices_[tri[1]], m.vertices_[tri[2]]) <= 0)
                return std::nullopt;
            for (unsigned k = 0; k < 3; ++k)
                if (!directed.emplace(key(tri[k], tri[(k + 1) % 3]), t).second)
                    return std::nullopt;
        }
        for (std::uint32_t t = 0; t < n; ++t) {
            const auto& tri = m.triangles_[t];
            for (unsigned k = 0; k < 3; ++k)
                if (const auto it = directed.find(key(tri[(k + 1) % 3], tri[k]));
                    it != directed.end())
                    m.neighbours_[t][k] = it->second;
        }
        return m;
    }

    [[nodiscard]] std::span<const MeshVertex> vertices() const noexcept { return vertices_; }
    [[nodiscard]] std::span<const TriangleIndices> triangles() const noexcept {
        return triangles_;
    }
    [[nodiscard]] std::size_t triangle_count() const noexcept { return triangles_.size(); }
    [[nodiscard]] std::array<std::uint32_t, 3> neighbours(std::size_t t) const noexcept {
        return neighbours_[t];
    }
    [[nodiscard]] bool is_constrained(std::size_t t, unsigned e) const noexcept {
        return (constrained_[t] >> e & 1u) != 0;
    }
    [[nodiscard]] std::uint32_t mask(std::size_t t, unsigned e) const noexcept {
        return masks_[t][e];
    }
    [[nodiscard]] MeshVertex corner(std::size_t t, unsigned k) const noexcept {
        return vertices_[triangles_[t][k]];
    }

    // Fan t into three around p, which must lie strictly inside t. The
    // children are (v0, v1, p) in t's slot, then (v1, v2, p) and (v2, v0, p)
    // appended. t's edges are unchanged, so no neighbour sees a new vertex.
    // Returns p's vertex index.
    std::uint32_t split_inside(std::uint32_t t, LatticeVertex p) {
        const auto q = add_vertex(p);
        const auto [v0, v1, v2] = triangles_[t];
        const std::array<Side, 3> s{side(t, 0), side(t, 1), side(t, 2)};
        const auto b = next_slot(), c = b + 1;
        put(t, {v0, v1, q}, {s[0], spoke(b), spoke(c)});
        put(b, {v1, v2, q}, {s[1], spoke(c), spoke(t)});
        put(c, {v2, v0, q}, {s[2], spoke(t), spoke(b)});
        repoint(s[1].neighbour, t, b);
        repoint(s[2].neighbour, t, c);
        return q;
    }

    // Split t's edge e at p, which must lie strictly inside that edge. If a
    // neighbour u shares the edge it is split at p too (2 -> 4), otherwise t
    // alone (1 -> 2). With the edge as (a, b) and c opposite: (a, p, c) takes
    // t's slot and (p, b, c) is appended; u's (b, a, d) becomes (b, p, d) in
    // u's slot and (p, a, d) appended. Returns p's vertex index.
    std::uint32_t split_edge(std::uint32_t t, unsigned e, LatticeVertex p) {
        const auto q = add_vertex(p);
        const auto a = triangles_[t][e], b = triangles_[t][(e + 1) % 3],
                   c = triangles_[t][(e + 2) % 3];
        const Side half = side(t, e), bc = side(t, (e + 1) % 3), ca = side(t, (e + 2) % 3);
        const auto u = half.neighbour;
        const auto t2 = next_slot();
        const bool shared = u != kNoNeighbour;
        const auto u2 = shared ? t2 + 1 : kNoNeighbour;

        put(t, {a, q, c}, {Side{u2, half.constrained, half.mask}, spoke(t2), ca});
        put(t2, {q, b, c}, {Side{shared ? u : kNoNeighbour, half.constrained, half.mask}, bc,
                            spoke(t)});
        repoint(bc.neighbour, t, t2);
        if (!shared)
            return q;

        unsigned f = 0;
        while (neighbours_[u][f] != t)
            ++f;
        const auto d = triangles_[u][(f + 2) % 3];
        const Side ad = side(u, (f + 1) % 3), db = side(u, (f + 2) % 3);
        put(u, {b, q, d}, {Side{t2, half.constrained, half.mask}, spoke(u2), db});
        put(u2, {q, a, d}, {Side{t, half.constrained, half.mask}, ad, spoke(u)});
        repoint(ad.neighbour, u, u2);
        return q;
    }

    // Flip t's edge e, which must have a neighbour u and a strictly convex
    // quad (increment 14b, R4). With t = (a, b, c), e the edge (a, b) and
    // u = (b, a, d): t's slot becomes (c, a, d) and u's (c, d, b). c is at
    // index 0 of both, so an inserted c is opposite edge 1 in both. The four
    // outer sides keep their neighbour, bit and mask; the new diagonal c-d
    // has neither bit nor mask.
    void flip(std::uint32_t t, unsigned e) {
        const auto a = triangles_[t][e], b = triangles_[t][(e + 1) % 3],
                   c = triangles_[t][(e + 2) % 3];
        const auto u = neighbours_[t][e];
        unsigned f = 0;
        while (triangles_[u][f] != b)
            ++f;
        const auto d = triangles_[u][(f + 2) % 3];
        assert(orient_sign(vertices_[c], vertices_[a], vertices_[d]) > 0
               && orient_sign(vertices_[c], vertices_[d], vertices_[b]) > 0);
        const Side bc = side(t, (e + 1) % 3), ca = side(t, (e + 2) % 3),
                   ad = side(u, (f + 1) % 3), db = side(u, (f + 2) % 3);
        put(t, {c, a, d}, {ca, ad, spoke(u)});
        put(u, {c, d, b}, {spoke(t), db, bc});
        repoint(ad.neighbour, u, t);
        repoint(bc.neighbour, t, u);
    }

    // The constrained edges, each once (from the lower-indexed triangle when
    // both sides are present), as vertex-index pairs, with their masks.
    [[nodiscard]] std::pair<std::vector<std::array<std::uint32_t, 2>>, std::vector<std::uint32_t>>
    constraint_edges() const {
        std::pair<std::vector<std::array<std::uint32_t, 2>>, std::vector<std::uint32_t>> out;
        for (std::uint32_t t = 0; t < triangles_.size(); ++t)
            for (unsigned k = 0; k < 3; ++k) {
                const auto n = neighbours_[t][k];
                if (is_constrained(t, k) && (n == kNoNeighbour || n > t)) {
                    out.first.push_back({triangles_[t][k], triangles_[t][(k + 1) % 3]});
                    out.second.push_back(masks_[t][k]);
                }
            }
        return out;
    }

private:
    struct Side {
        std::uint32_t neighbour;
        bool constrained;
        std::uint32_t mask;
    };

    LatticeMesh() = default;

    [[nodiscard]] Side side(std::uint32_t t, unsigned k) const noexcept {
        return Side{neighbours_[t][k], is_constrained(t, k), masks_[t][k]};
    }
    [[nodiscard]] static constexpr Side spoke(std::uint32_t neighbour) noexcept {
        return Side{neighbour, false, 0};
    }
    [[nodiscard]] std::uint32_t next_slot() const noexcept {
        return static_cast<std::uint32_t>(triangles_.size());
    }
    std::uint32_t add_vertex(LatticeVertex p) {
        vertices_.push_back(p);
        return static_cast<std::uint32_t>(vertices_.size() - 1);
    }

    // Write a triangle into slot t, appending when t is one past the end.
    // Slots are always filled in order, so t is never further out than that.
    void put(std::uint32_t t, TriangleIndices tri, std::array<Side, 3> s) {
        if (t == triangles_.size()) {
            triangles_.emplace_back();
            neighbours_.emplace_back();
            constrained_.emplace_back();
            masks_.emplace_back();
        }
        triangles_[t] = tri;
        constrained_[t] = 0;
        for (unsigned k = 0; k < 3; ++k) {
            neighbours_[t][k] = s[k].neighbour;
            constrained_[t] |= static_cast<std::uint8_t>(s[k].constrained ? 1u << k : 0u);
            masks_[t][k] = s[k].mask;
        }
    }

    // Tell n that its neighbour `from` is now `to`.
    void repoint(std::uint32_t n, std::uint32_t from, std::uint32_t to) noexcept {
        if (n == kNoNeighbour)
            return;
        for (auto& x : neighbours_[n])
            if (x == from)
                x = to;
    }

    std::vector<MeshVertex> vertices_;
    std::vector<TriangleIndices> triangles_;
    std::vector<std::array<std::uint32_t, 3>> neighbours_;
    std::vector<std::uint8_t> constrained_;
    std::vector<std::array<std::uint32_t, 3>> masks_;
};

}  // namespace terrain::mesh
