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
#include <cmath>
#include <cstdint>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::MeshVertex;
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
    return {rc(v[tri[0]].as_node().value()), rc(v[tri[1]].as_node().value()), rc(v[tri[2]].as_node().value())};
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
            REQUIRE_FALSE(on_open_segment(rc(v[edge.first].as_node().value()), rc(v[edge.second].as_node().value()), rc(v[i].as_node().value())));
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

// Increment 16 review ruling: the checked back-conversion MeshVertex ->
// LatticeVertex is a named query, not an implicit throwing conversion, so a
// node/off-node mix-up is a compile error or an empty optional, never a throw
// found at run time. The widening LatticeVertex -> MeshVertex stays implicit.
static_assert(std::is_convertible_v<LatticeVertex, MeshVertex>);
static_assert(!std::is_convertible_v<MeshVertex, LatticeVertex>);
static_assert(std::is_same_v<decltype(std::declval<const MeshVertex&>().as_node()),
                             std::optional<LatticeVertex>>);

TEST_CASE("as_node is exact for a node and empty for an off-node vertex", "[refinement][lattice]") {
    REQUIRE(MeshVertex{7.0, 3.0}.as_node() == LatticeVertex{3, 7});
    REQUIRE(MeshVertex{lv(4'000'000'000u, 0)}.as_node() == lv(4'000'000'000u, 0));
    REQUIRE(MeshVertex{0.0, 0.0}.as_node() == lv(0, 0));
    REQUIRE_FALSE(MeshVertex{7.5, 3.0}.as_node().has_value());
    REQUIRE_FALSE(MeshVertex{7.0, 3.0 + 1e-9}.as_node().has_value());
    REQUIRE_FALSE(MeshVertex{0.625, 0.625}.as_node().has_value());
}

// ---------------------------------------------------------------------------
// Increment 15f-4 (docs/increments/15f-edge-strip.md, "Settled after 15f-3's
// acceptance", A4, B1 to B3): build's refusals and its neighbour table, pinned
// against the unordered_map build of 15f-3 before A2 replaces it with a
// vertex-bucketed table. These are guards: they pass before the change and
// must still pass after it.
// ---------------------------------------------------------------------------

namespace {

using Masks = std::vector<std::array<std::uint32_t, 3>>;

// (a, b, c) or (a, c, b), whichever is counter-clockwise in build's frame, so
// a fixture never depends on hand-ordering corners. REQUIREs a non-degenerate
// triple: a zero-area fixture would make every refusal below vacuous.
TriangleIndices ccw(const std::vector<MeshVertex>& v, std::uint32_t a, std::uint32_t b,
                    std::uint32_t c) {
    const int s = terrain::mesh::orient_sign(v[a], v[b], v[c]);
    REQUIRE(s != 0);
    return s > 0 ? TriangleIndices{a, b, c} : TriangleIndices{a, c, b};
}

std::optional<LatticeMesh> build_plain(const std::vector<MeshVertex>& v,
                                       const std::vector<TriangleIndices>& tris) {
    return LatticeMesh::build(v, tris, std::vector<std::uint8_t>(tris.size(), 0),
                              Masks(tris.size(), {0, 0, 0}));
}

// The oracle: for every ordered pair of triangles and every pair of edges, u
// is t's neighbour across edge k exactly when u holds that edge reversed.
// Quadratic, independent of build's lookup structure.
std::vector<std::array<std::uint32_t, 3>> brute_neighbours(
    const std::vector<TriangleIndices>& tris) {
    std::vector<std::array<std::uint32_t, 3>> out(tris.size(),
                                                  {kNoNeighbour, kNoNeighbour, kNoNeighbour});
    for (std::uint32_t t = 0; t < tris.size(); ++t)
        for (std::uint32_t u = 0; u < tris.size(); ++u)
            for (unsigned k = 0; k < 3; ++k)
                for (unsigned j = 0; j < 3; ++j)
                    if (u != t && tris[t][k] == tris[u][(j + 1) % 3]
                        && tris[t][(k + 1) % 3] == tris[u][j]) {
                        REQUIRE(out[t][k] == kNoNeighbour);  // the fixture is a manifold
                        out[t][k] = u;
                    }
    return out;
}

// Builds with per-triangle bits and masks drawn from rng, and checks the
// neighbour table against the oracle and every slot's bits and masks against
// the input (build moves them; a re-indexing table must not reorder them).
void check_against_oracle(const std::vector<MeshVertex>& v,
                          const std::vector<TriangleIndices>& tris, std::mt19937& rng) {
    std::vector<std::uint8_t> bits;
    Masks masks;
    for (std::size_t t = 0; t < tris.size(); ++t) {
        bits.push_back(static_cast<std::uint8_t>(rng() % 8));
        masks.push_back({static_cast<std::uint32_t>(rng()), static_cast<std::uint32_t>(rng()),
                         static_cast<std::uint32_t>(rng())});
    }
    const auto m = LatticeMesh::build(v, tris, bits, masks);
    REQUIRE(m.has_value());
    REQUIRE(m->triangle_count() == tris.size());
    const auto expected = brute_neighbours(tris);
    std::size_t boundary = 0;
    for (std::size_t t = 0; t < tris.size(); ++t) {
        CAPTURE(t);
        REQUIRE(m->triangles()[t] == tris[t]);
        REQUIRE(m->neighbours(t) == expected[t]);
        for (unsigned k = 0; k < 3; ++k) {
            REQUIRE(m->is_constrained(t, k) == ((bits[t] >> k & 1u) != 0));
            REQUIRE(m->mask(t, k) == masks[t][k]);
            boundary += expected[t][k] == kNoNeighbour ? 1 : 0;
        }
    }
    REQUIRE(boundary > 0);  // the oracle saw both kinds of edge
    REQUIRE(boundary < 3 * tris.size());
}

// Fisher-Yates over mt19937's raw output: std::shuffle and the standard
// distributions are implementation-defined, this is the same on every
// standard library.
template <typename T>
void shuffle(std::vector<T>& xs, std::mt19937& rng) {
    for (std::size_t i = xs.size(); i > 1; --i)
        std::swap(xs[i - 1], xs[rng() % i]);
}

}  // namespace

TEST_CASE("B1: build refuses a directed edge used twice", "[refinement][lattice][15f-4]") {
    // (col, row): a-b is the shared edge, c and e above it, d below it.
    const std::vector<MeshVertex> v{{0, 2}, {2, 2}, {1, 0}, {1, 4}, {1, 1}};
    const auto abc = ccw(v, 0, 1, 2), abd = ccw(v, 0, 1, 3), abe = ccw(v, 0, 1, 4);
    REQUIRE(build_plain(v, {abc, abd}).has_value());  // control: one edge, two sides

    SECTION("two coincident triangles, the same slots") {
        REQUIRE_FALSE(build_plain(v, {abc, abc}).has_value());
    }
    SECTION("two coincident triangles, rotated") {
        const TriangleIndices rotated{abc[1], abc[2], abc[0]};
        REQUIRE_FALSE(build_plain(v, {abc, rotated}).has_value());
    }
    SECTION("coincident triangles that are not adjacent in the array") {
        REQUIRE_FALSE(build_plain(v, {abc, abd, abc}).has_value());
    }
    SECTION("three triangles on one edge") {
        REQUIRE_FALSE(build_plain(v, {abc, abd, abe}).has_value());
        REQUIRE_FALSE(build_plain(v, {abe, abd, abc}).has_value());
    }
    SECTION("two triangles overlapping on one side of an edge") {
        REQUIRE_FALSE(build_plain(v, {abc, abe}).has_value());
    }
}

TEST_CASE("B1: build refuses a degenerate triangle from a duplicate vertex",
          "[refinement][lattice][15f-4]") {
    // Two indices at one position: a zero-area triangle, refused like 14's
    // collinear case. A repeated index is the same thing.
    const std::vector<MeshVertex> v{{0, 2}, {2, 2}, {1, 0}, {0, 2}};
    REQUIRE_FALSE(build_plain(v, {TriangleIndices{0, 3, 1}}).has_value());
    REQUIRE_FALSE(build_plain(v, {TriangleIndices{0, 0, 1}}).has_value());
    REQUIRE(build_plain(v, {ccw(v, 3, 1, 2)}).has_value());  // the duplicate alone is fine
}

TEST_CASE("B2: build refuses a vertex index out of range", "[refinement][lattice][15f-4]") {
    const std::vector<MeshVertex> v{{0, 2}, {2, 2}, {1, 0}, {1, 4}};
    const auto abc = ccw(v, 0, 1, 2), abd = ccw(v, 0, 1, 3);
    REQUIRE(build_plain(v, {abc, abd}).has_value());  // control

    const auto n = static_cast<std::uint32_t>(v.size());
    const auto index = GENERATE_COPY(n, n + 1, kNoNeighbour);
    const auto slot = GENERATE(0u, 1u, 2u);
    CAPTURE(index, slot);
    auto bad = abd;
    bad[slot] = index;
    REQUIRE_FALSE(build_plain(v, {bad}).has_value());
    REQUIRE_FALSE(build_plain(v, {abc, bad}).has_value());  // after a valid triangle
    REQUIRE_FALSE(build_plain({}, {abc}).has_value());     // no vertices at all
}

TEST_CASE("B2: build refuses mismatched array lengths", "[refinement][lattice][15f-4]") {
    const std::vector<MeshVertex> v{{0, 2}, {2, 2}, {1, 0}, {1, 4}};
    const std::vector<TriangleIndices> tris{ccw(v, 0, 1, 2), ccw(v, 0, 1, 3)};
    const std::vector<std::uint8_t> bits(2, 0);
    const Masks masks(2, {0, 0, 0});
    REQUIRE(LatticeMesh::build(v, tris, bits, masks).has_value());  // control

    SECTION("constrained short or long") {
        REQUIRE_FALSE(LatticeMesh::build(v, tris, {0}, masks).has_value());
        REQUIRE_FALSE(LatticeMesh::build(v, tris, {0, 0, 0}, masks).has_value());
        REQUIRE_FALSE(LatticeMesh::build(v, tris, {}, masks).has_value());
    }
    SECTION("masks short or long") {
        REQUIRE_FALSE(LatticeMesh::build(v, tris, bits, Masks(1, {0, 0, 0})).has_value());
        REQUIRE_FALSE(LatticeMesh::build(v, tris, bits, Masks(3, {0, 0, 0})).has_value());
        REQUIRE_FALSE(LatticeMesh::build(v, tris, bits, Masks{}).has_value());
    }
    SECTION("no triangles, but bits or masks") {
        REQUIRE_FALSE(LatticeMesh::build(v, {}, {0}, Masks{}).has_value());
        REQUIRE_FALSE(LatticeMesh::build(v, {}, {}, Masks(1, {0, 0, 0})).has_value());
    }
}

TEST_CASE("B2: build accepts an empty mesh and unreferenced vertices",
          "[refinement][lattice][15f-4]") {
    const auto empty = LatticeMesh::build(std::vector<MeshVertex>{}, {}, {}, Masks{});
    REQUIRE(empty.has_value());
    REQUIRE(empty->triangle_count() == 0);

    // One triangle among unused vertices (to_lattice passes every start
    // vertex, referenced or not): all three edges are boundary.
    const std::vector<MeshVertex> v{{5, 5}, {0, 2}, {2, 2}, {9, 9}, {1, 0}};
    const auto one = build_plain(v, {ccw(v, 1, 2, 4)});
    REQUIRE(one.has_value());
    REQUIRE(one->vertices().size() == v.size());
    REQUIRE(one->neighbours(0) == std::array<std::uint32_t, 3>{kNoNeighbour, kNoNeighbour,
                                                                kNoNeighbour});
}

TEST_CASE("B3: build's neighbours equal the brute-force oracle on a shuffled grid",
          "[refinement][lattice][15f-4]") {
    // A cols x rows node grid, each cell cut along a random diagonal, then
    // triangle order and each triangle's first corner shuffled. Long straight
    // grid lines give long collinear vertex runs along the boundary.
    const auto seed = GENERATE(1u, 2u, 3u);
    std::mt19937 rng{seed};
    constexpr std::uint32_t kCols = 31, kRows = 23;
    std::vector<MeshVertex> v;
    for (std::uint32_t r = 0; r < kRows; ++r)
        for (std::uint32_t c = 0; c < kCols; ++c)
            v.emplace_back(static_cast<double>(c), static_cast<double>(r));
    std::vector<TriangleIndices> tris;
    for (std::uint32_t r = 0; r + 1 < kRows; ++r)
        for (std::uint32_t c = 0; c + 1 < kCols; ++c) {
            const auto tl = r * kCols + c, tr = tl + 1, bl = tl + kCols, br = bl + 1;
            if (rng() % 2 == 0) {
                tris.push_back(ccw(v, tl, bl, br));
                tris.push_back(ccw(v, tl, br, tr));
            } else {
                tris.push_back(ccw(v, tl, bl, tr));
                tris.push_back(ccw(v, tr, bl, br));
            }
        }
    shuffle(tris, rng);
    for (auto& t : tris)
        std::rotate(t.begin(), t.begin() + rng() % 3, t.end());
    CAPTURE(seed);
    check_against_oracle(v, tris, rng);
}

TEST_CASE("B3: build's neighbours equal the brute-force oracle on a high-degree fan",
          "[refinement][lattice][15f-4]") {
    // A centre joined to n cocircular rim points: degree n >= 1000 at the
    // centre. Closed (the rim wraps around) and open (one wedge missing, so
    // the centre has two boundary edges). Also at a lattice offset of 1e6
    // with a 1e-3 radius, where the rim's spacing is about 1e-10 of the
    // coordinates' magnitude.
    const auto n = GENERATE(1000u, 1531u);
    const bool closed = GENERATE(true, false);
    const auto [offset, radius] =
        GENERATE(std::pair{5000.0, 4000.0}, std::pair{1.0e6, 1.0e-3});
    CAPTURE(n, closed, offset, radius);
    constexpr double kTau = 6.283185307179586;
    std::vector<MeshVertex> v{{offset, offset}};
    for (std::uint32_t i = 0; i < n; ++i) {
        const double a = kTau * i / n;
        v.emplace_back(offset + radius * std::cos(a), offset + radius * std::sin(a));
    }
    std::vector<TriangleIndices> tris;
    for (std::uint32_t i = 0; i + (closed ? 0 : 1) < n; ++i)
        tris.push_back(ccw(v, 0, 1 + i, 1 + (i + 1) % n));
    std::mt19937 rng{n};
    shuffle(tris, rng);
    check_against_oracle(v, tris, rng);

    // Every spoke is interior except, on the open fan, the two at the gap.
    const auto m = build_plain(v, tris);
    REQUIRE(m.has_value());
    std::size_t centre_boundary = 0;
    for (std::size_t t = 0; t < tris.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if ((tris[t][k] == 0 || tris[t][(k + 1) % 3] == 0)
                && m->neighbours(t)[k] == kNoNeighbour)
                ++centre_boundary;
    REQUIRE(centre_boundary == (closed ? 0u : 2u));
}
