// Increment 14 (docs/increments/14-adaptive-refinement.md, R3 and R4):
// LatticeMesh and its three splits. INVARIANT-CRITICAL: conformity is decided
// here and nowhere else, so @reviewer mutation-tests this suite.
//
// Interface: include/terrain/mesh/lattice_mesh.hpp.
//
// Edge k runs from vertex k to vertex k+1, as in IndexedMesh2. The split
// functions return the new vertex's index. build() refuses (nullopt) a
// triangle whose orientation is not positive in the WORLD's handedness
// (x = col, y = -row) -- R4 says "positive integer orientation" without a
// frame; this suite picks the one that keeps a CDT's CCW output CCW.
//
// Child ORDER is pinned only where R4 pins it: the parent's slot holds a child.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/mesh/lattice_mesh.hpp>

#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using refinement_fixtures::in_closed;
using refinement_fixtures::on_open_segment;
using refinement_fixtures::orient;
using refinement_fixtures::RC;

namespace {

RC rc(LatticeVertex v) { return RC{v.row, v.col}; }
LatticeVertex lv(std::uint32_t row, std::uint32_t col) { return LatticeVertex{row, col}; }

std::array<RC, 3> corners(const LatticeMesh& m, std::size_t t) {
    const auto& tri = m.triangles()[t];
    const auto v = m.vertices();
    return {rc(v[tri[0]]), rc(v[tri[1]]), rc(v[tri[2]])};
}

std::int64_t area2_sum(const LatticeMesh& m) {
    std::int64_t s = 0;
    for (std::size_t t = 0; t < m.triangle_count(); ++t) {
        const auto c = corners(m, t);
        s += orient(c[0], c[1], c[2]);
    }
    return s;
}

// Recomputes adjacency from the triangles alone and compares; checks
// orientation, edge manifoldness and the absence of hanging vertices.
void check_topology(const LatticeMesh& m) {
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> directed;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        const auto c = corners(m, t);
        CAPTURE(t);
        REQUIRE(orient(c[0], c[1], c[2]) > 0);
        for (unsigned k = 0; k < 3; ++k) {
            const auto inserted = directed.emplace(std::pair{tri[k], tri[(k + 1) % 3]}, t).second;
            REQUIRE(inserted);  // a directed edge in two triangles: overlap or flip
        }
    }
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        for (unsigned k = 0; k < 3; ++k) {
            CAPTURE(t, k);
            const auto it = directed.find({tri[(k + 1) % 3], tri[k]});
            const std::uint32_t expected = it == directed.end() ? kNoNeighbour : it->second;
            REQUIRE(m.neighbours(t)[k] == expected);
        }
    }
    const auto v = m.vertices();
    for (const auto& [edge, t] : directed)
        for (std::size_t i = 0; i < v.size(); ++i) {
            CAPTURE(edge.first, edge.second, i);
            REQUIRE_FALSE(on_open_segment(rc(v[edge.first]), rc(v[edge.second]), rc(v[i])));
        }
}

// The (constrained, mask) of the undirected edge (a, b), read from whichever
// triangle has it. Interior edges must agree from both sides.
std::optional<std::pair<bool, std::uint32_t>> edge_info(const LatticeMesh& m, RC a, RC b) {
    std::optional<std::pair<bool, std::uint32_t>> found;
    for (std::size_t t = 0; t < m.triangle_count(); ++t) {
        const auto c = corners(m, t);
        for (unsigned k = 0; k < 3; ++k) {
            const RC p = c[k], q = c[(k + 1) % 3];
            if ((p == a && q == b) || (p == b && q == a)) {
                const std::pair info{m.is_constrained(t, k), m.mask(t, k)};
                if (found) REQUIRE(*found == info);
                found = info;
            }
        }
    }
    return found;
}

// The 3 x 3-node square split along its (0,0)-(2,2) diagonal.
//   t0 = (0,0) (2,0) (2,2): left side mask 8, bottom mask 4, diagonal free
//   t1 = (0,0) (2,2) (0,2): diagonal free, right mask 2, top mask 1
LatticeMesh square() {
    auto m = LatticeMesh::build({lv(0, 0), lv(2, 0), lv(2, 2), lv(0, 2)},
                                {TriangleIndices{0, 1, 2}, TriangleIndices{0, 2, 3}},
                                {0b011, 0b110},
                                {std::array<std::uint32_t, 3>{8, 4, 0}, {0, 2, 1}});
    REQUIRE(m.has_value());
    return *m;
}

// One triangle whose only strictly interior node is (2, 1) and whose edges
// hold two nodes each (Pick: A = 4.5, B = 9, I = 1).
LatticeMesh wedge() {
    auto m = LatticeMesh::build({lv(0, 0), lv(3, 0), lv(3, 3)}, {TriangleIndices{0, 1, 2}},
                                {0b111}, {std::array<std::uint32_t, 3>{8, 4, 16}});
    REQUIRE(m.has_value());
    return *m;
}

}  // namespace

TEST_CASE("build derives adjacency across the shared edge", "[refinement][lattice]") {
    const auto m = square();
    REQUIRE(m.neighbours(0) == std::array<std::uint32_t, 3>{kNoNeighbour, kNoNeighbour, 1});
    REQUIRE(m.neighbours(1) == std::array<std::uint32_t, 3>{0, kNoNeighbour, kNoNeighbour});
    REQUIRE(m.is_constrained(0, 0));
    REQUIRE_FALSE(m.is_constrained(0, 2));
    REQUIRE(m.mask(1, 2) == 1);
    check_topology(m);
}

TEST_CASE("build refuses a clockwise or a zero-area triangle", "[refinement][lattice]") {
    SECTION("clockwise") {
        REQUIRE_FALSE(LatticeMesh::build({lv(0, 0), lv(2, 0), lv(2, 2)},
                                         {TriangleIndices{0, 2, 1}}, {0},
                                         {std::array<std::uint32_t, 3>{}})
                          .has_value());
    }
    SECTION("collinear") {
        REQUIRE_FALSE(LatticeMesh::build({lv(0, 0), lv(1, 1), lv(2, 2)},
                                         {TriangleIndices{0, 1, 2}}, {0},
                                         {std::array<std::uint32_t, 3>{}})
                          .has_value());
    }
}

TEST_CASE("split_inside fans one triangle into three around the point", "[refinement][lattice]") {
    auto m = wedge();
    const auto before = area2_sum(m);
    const auto p = lv(2, 1);
    const auto idx = m.split_inside(0, p);

    REQUIRE(idx == 3);
    REQUIRE(m.vertices()[idx] == p);
    REQUIRE(m.triangle_count() == 3);
    for (std::size_t t = 0; t < 3; ++t) {
        const auto& tri = m.triangles()[t];
        REQUIRE(std::count(tri.begin(), tri.end(), idx) == 1);
    }
    REQUIRE(area2_sum(m) == before);
    check_topology(m);

    // The parent's edges keep their bit and mask; the three spokes are free.
    REQUIRE(edge_info(m, {0, 0}, {3, 0}) == std::pair{true, std::uint32_t{8}});
    REQUIRE(edge_info(m, {3, 0}, {3, 3}) == std::pair{true, std::uint32_t{4}});
    REQUIRE(edge_info(m, {3, 3}, {0, 0}) == std::pair{true, std::uint32_t{16}});
    for (const RC corner : {RC{0, 0}, RC{3, 0}, RC{3, 3}})
        REQUIRE(edge_info(m, corner, rc(p)) == std::pair{false, std::uint32_t{0}});
}

TEST_CASE("split_edge on a shared edge splits both sides: two become four",
          "[refinement][lattice]") {
    // From either side, the same four triangles.
    const auto [t, e] = GENERATE(std::pair{0u, 2u}, std::pair{1u, 0u});
    CAPTURE(t, e);
    auto m = square();
    const auto before = area2_sum(m);
    const auto idx = m.split_edge(t, e, lv(1, 1));

    REQUIRE(m.triangle_count() == 4);
    REQUIRE(m.vertices()[idx] == lv(1, 1));
    for (std::size_t k = 0; k < 4; ++k) {
        const auto& tri = m.triangles()[k];
        REQUIRE(std::count(tri.begin(), tri.end(), idx) == 1);
    }
    REQUIRE(area2_sum(m) == before);
    check_topology(m);  // includes: (1,1) hangs on no edge
    REQUIRE_FALSE(edge_info(m, {0, 0}, {2, 2}).has_value());  // the old edge is gone
    REQUIRE(edge_info(m, {0, 0}, {1, 1}) == std::pair{false, std::uint32_t{0}});
    REQUIRE(edge_info(m, {2, 0}, {2, 2}) == std::pair{true, std::uint32_t{4}});
    REQUIRE(edge_info(m, {0, 2}, {0, 0}) == std::pair{true, std::uint32_t{1}});
}

TEST_CASE("split_edge on an interior constrained edge: both halves inherit on both sides",
          "[refinement][lattice]") {
    // square() with its diagonal constrained under mask 16 on both sides. The
    // only fixture where a constrained edge has a triangle on each side, so
    // the only one that sees u's and u2's inherited (bit, mask).
    const auto [t, e] = GENERATE(std::pair{0u, 2u}, std::pair{1u, 0u});
    CAPTURE(t, e);
    auto built = LatticeMesh::build({lv(0, 0), lv(2, 0), lv(2, 2), lv(0, 2)},
                                    {TriangleIndices{0, 1, 2}, TriangleIndices{0, 2, 3}},
                                    {0b111, 0b111},
                                    {std::array<std::uint32_t, 3>{8, 4, 16}, {16, 2, 1}});
    REQUIRE(built.has_value());
    auto m = *built;
    m.split_edge(t, e, lv(1, 1));

    REQUIRE(m.triangle_count() == 4);
    check_topology(m);
    const std::pair constrained16{true, std::uint32_t{16}};
    // Read each half from every triangle carrying it: two halves, two sides.
    std::size_t sides_seen = 0;
    for (std::size_t tri = 0; tri < m.triangle_count(); ++tri) {
        const auto c = corners(m, tri);
        for (unsigned k = 0; k < 3; ++k) {
            const RC a = c[k], b = c[(k + 1) % 3];
            const bool on_diagonal = (a == RC{1, 1} || b == RC{1, 1}) &&
                                     (a == RC{0, 0} || a == RC{2, 2} || b == RC{0, 0} ||
                                      b == RC{2, 2});
            if (!on_diagonal) continue;
            CAPTURE(tri, k);
            REQUIRE(std::pair{m.is_constrained(tri, k), m.mask(tri, k)} == constrained16);
            ++sides_seen;
        }
    }
    REQUIRE(sides_seen == 4);
    REQUIRE(edge_info(m, {0, 0}, {1, 1}) == constrained16);
    REQUIRE(edge_info(m, {1, 1}, {2, 2}) == constrained16);
    REQUIRE(edge_info(m, {2, 0}, {1, 1}) == std::pair{false, std::uint32_t{0}});
    REQUIRE(edge_info(m, {0, 2}, {1, 1}) == std::pair{false, std::uint32_t{0}});
}

TEST_CASE("split_edge on the boundary splits one triangle; both halves inherit",
          "[refinement][lattice]") {
    auto m = square();
    const auto before = area2_sum(m);
    m.split_edge(1, 2, lv(0, 1));  // t1's top edge, (0,2) -> (0,0), mask 1

    REQUIRE(m.triangle_count() == 3);
    REQUIRE(area2_sum(m) == before);
    check_topology(m);
    REQUIRE_FALSE(edge_info(m, {0, 0}, {0, 2}).has_value());
    REQUIRE(edge_info(m, {0, 0}, {0, 1}) == std::pair{true, std::uint32_t{1}});
    REQUIRE(edge_info(m, {0, 1}, {0, 2}) == std::pair{true, std::uint32_t{1}});
    REQUIRE(edge_info(m, {0, 1}, {2, 2}) == std::pair{false, std::uint32_t{0}});
    // t0 was not touched.
    REQUIRE(corners(m, 0) == std::array<RC, 3>{RC{0, 0}, RC{2, 0}, RC{2, 2}});
}

TEST_CASE("a long deterministic sequence of splits stays conforming", "[refinement][lattice]") {
    // 9 x 9 nodes at stride 4: eight triangles. Split, in turn, the triangle at
    // a rotating index at the first node (row-major) of its closed set that is
    // not a vertex -- inside or on an edge, whichever it is -- and check the
    // whole mesh after every split.
    std::vector<LatticeVertex> vs;
    for (std::uint32_t r = 0; r <= 8; r += 4)
        for (std::uint32_t c = 0; c <= 8; c += 4) vs.push_back(lv(r, c));
    std::vector<TriangleIndices> tris;
    std::vector<std::uint8_t> con;
    std::vector<std::array<std::uint32_t, 3>> masks;
    for (std::uint32_t i = 0; i < 2; ++i)
        for (std::uint32_t j = 0; j < 2; ++j) {
            const auto tl = i * 3 + j, tr = tl + 1, bl = tl + 3, br = tl + 4;
            tris.push_back({tl, bl, br});
            tris.push_back({tl, br, tr});
            // Boundary edges: left (tl-bl) when j == 0, bottom (bl-br) when i == 1,
            // right (br-tr) when j == 1, top (tr-tl) when i == 0.
            con.push_back(static_cast<std::uint8_t>((j == 0 ? 1 : 0) | (i == 1 ? 2 : 0)));
            masks.push_back({j == 0 ? 8u : 0u, i == 1 ? 4u : 0u, 0u});
            con.push_back(static_cast<std::uint8_t>((j == 1 ? 2 : 0) | (i == 0 ? 4 : 0)));
            masks.push_back({0u, j == 1 ? 2u : 0u, i == 0 ? 1u : 0u});
        }
    auto built = LatticeMesh::build(vs, tris, con, masks);
    REQUIRE(built.has_value());
    auto m = *built;
    const auto before = area2_sum(m);

    int splits = 0;
    for (std::uint32_t step = 0; step < 400 && splits < 40; ++step) {
        const auto t = static_cast<std::uint32_t>((step * 7u) % m.triangle_count());
        const auto c = corners(m, t);
        std::optional<RC> pick;
        for (std::int64_t r = 0; r <= 8 && !pick; ++r)
            for (std::int64_t col = 0; col <= 8 && !pick; ++col) {
                const RC q{r, col};
                if (q != c[0] && q != c[1] && q != c[2] && in_closed(c[0], c[1], c[2], q)) pick = q;
            }
        if (!pick) continue;
        const auto p = lv(static_cast<std::uint32_t>(pick->row), static_cast<std::uint32_t>(pick->col));
        unsigned edge = 3;
        for (unsigned k = 0; k < 3; ++k)
            if (orient(c[k], c[(k + 1) % 3], *pick) == 0) edge = k;
        if (edge == 3) m.split_inside(t, p);
        else m.split_edge(t, edge, p);
        ++splits;
        CAPTURE(step, t, edge);
        REQUIRE(area2_sum(m) == before);
        check_topology(m);
    }
    REQUIRE(splits == 40);

    // Every boundary edge is constrained with its side's mask; nothing else is.
    for (std::size_t t = 0; t < m.triangle_count(); ++t) {
        const auto c = corners(m, t);
        for (unsigned k = 0; k < 3; ++k) {
            const RC a = c[k], b = c[(k + 1) % 3];
            std::uint32_t side = 0;
            if (a.row == 0 && b.row == 0) side = 1;
            else if (a.col == 8 && b.col == 8) side = 2;
            else if (a.row == 8 && b.row == 8) side = 4;
            else if (a.col == 0 && b.col == 0) side = 8;
            CAPTURE(t, k);
            REQUIRE(m.is_constrained(t, k) == (side != 0));
            REQUIRE(m.mask(t, k) == side);
            REQUIRE((m.neighbours(t)[k] == kNoNeighbour) == (side != 0));
        }
    }
}
