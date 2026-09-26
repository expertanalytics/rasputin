// Increment 14b (docs/increments/14b-delaunay-insertion.md, R1, R3, R4 and R5):
// LatticeMesh::flip and the two Lawson legalisers. INVARIANT-CRITICAL: the
// Delaunay property and conformity after a flip are decided here, so
// @reviewer mutation-tests this suite. Tests L1 to L4 of the design, plus L5,
// legalise_around after a real split.
//
// Interface: include/terrain/mesh/lawson.hpp and LatticeMesh::flip. Names the
// design leaves open, chosen here:
//
//   struct LatticeFrame { double dx; double dy; };   // point = (col*dx, -(row*dy))
//   void LatticeMesh::flip(std::uint32_t t, unsigned e);
//   template <pred::GeometryKernel K, class OnWrite>
//   std::size_t legalise_around(LatticeMesh&, std::uint32_t q,
//                               std::span<const std::uint32_t> seeds,
//                               const LatticeFrame&, OnWrite&& on_write);
//   template <pred::GeometryKernel K, class OnWrite>
//   std::size_t legalise_all(LatticeMesh&, const LatticeFrame&, OnWrite&& on_write);
//
// Both legalisers return the number of flips and call on_write(slot) for
// every slot a flip writes (repeats allowed).
//
// flip's vertex order is pinned, following split_edge's style: with t =
// (a, b, c), e the edge (a, b), and u = (b, a, d) across it, t's slot becomes
// (c, a, d) and u's slot becomes (c, d, b). c is at index 0 of both, so when c
// is the inserted vertex q, "the edge opposite q" is edge 1 in both (R1).
//
// The Delaunay oracle below is this file's own: it rebuilds the frame points
// from (row, col) and asks DefaultKernel::incircle, the producer's predicate
// (R3) but not its records.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <random>
#include <set>
#include <span>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::MeshVertex;
using terrain::mesh::orient_sign;
using terrain::mesh::legalise_all;
using terrain::mesh::legalise_around;
using terrain::pred::DefaultKernel;
using terrain::pred::Incircle;
using terrain::pred::Orientation;
using refinement_fixtures::on_open_segment;
using refinement_fixtures::orient;
using refinement_fixtures::RC;

namespace {

using Masks = std::array<std::uint32_t, 3>;

LatticeVertex lv(std::uint32_t row, std::uint32_t col) { return LatticeVertex{row, col}; }
RC rc(LatticeVertex v) { return RC{v.row, v.col}; }

// Takes a MeshVertex so off-node vertices work too; a LatticeVertex converts
// to it exactly.
Point2 at(const LatticeFrame& f, MeshVertex v) { return Point2{v.col * f.dx, -(v.row * f.dy)}; }

// Records every slot on_write reports.
struct Writes {
    std::set<std::uint32_t> slots;
    auto sink() {
        return [this](std::uint32_t s) { slots.insert(s); };
    }
};

// Adjacency recomputed from the triangles alone; orientation, manifold edges,
// no hanging vertex. Independent of LatticeMesh's own bookkeeping.
void check_topology(const LatticeMesh& m) {
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> directed;
    const auto v = m.vertices();
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        CAPTURE(t);
        REQUIRE(orient(rc(v[tri[0]].as_node().value()), rc(v[tri[1]].as_node().value()), rc(v[tri[2]].as_node().value())) > 0);
        for (unsigned k = 0; k < 3; ++k) {
            const bool fresh = directed.emplace(std::pair{tri[k], tri[(k + 1) % 3]}, t).second;
            REQUIRE(fresh);
        }
    }
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        for (unsigned k = 0; k < 3; ++k) {
            CAPTURE(t, k);
            const auto it = directed.find({tri[(k + 1) % 3], tri[k]});
            REQUIRE(m.neighbours(t)[k] == (it == directed.end() ? kNoNeighbour : it->second));
            if (it != directed.end()) {
                // Both sides of an interior edge agree on its bit and mask.
                const auto& o = m.triangles()[it->second];
                unsigned j = 0;
                while (o[j] != tri[(k + 1) % 3]) ++j;
                REQUIRE(m.is_constrained(t, k) == m.is_constrained(it->second, j));
                REQUIRE(m.mask(t, k) == m.mask(it->second, j));
            }
        }
    }
    for (const auto& [edge, t] : directed)
        for (std::size_t i = 0; i < v.size(); ++i)
            REQUIRE_FALSE(on_open_segment(rc(v[edge.first].as_node().value()), rc(v[edge.second].as_node().value()), rc(v[i].as_node().value())));
}

// Every unconstrained interior edge: the neighbour's apex is not strictly
// inside the triangle's circle, in the given frame.
std::size_t delaunay_violations(const LatticeMesh& m, const LatticeFrame& f) {
    std::size_t bad = 0;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto u = m.neighbours(t)[k];
            if (u == kNoNeighbour || m.is_constrained(t, k)) continue;
            const auto& tu = m.triangles()[u];
            const auto& tt = m.triangles()[t];
            std::uint32_t apex = tu[0];
            for (const auto x : tu)
                if (x != tt[k] && x != tt[(k + 1) % 3]) apex = x;
            const auto v = m.vertices();
            if (DefaultKernel::incircle(at(f, v[tt[0]]), at(f, v[tt[1]]), at(f, v[tt[2]]), at(f, v[apex]))
                == Incircle::Inside)
                ++bad;
        }
    return bad;
}

std::vector<TriangleIndices> snapshot(const LatticeMesh& m) {
    return {m.triangles().begin(), m.triangles().end()};
}

// Every slot whose triangle changed was reported through on_write. This is
// what refine's `touched` is built from (R2), so a slot missed here is a stale
// scan result there.
void check_writes_cover_changes(const std::vector<TriangleIndices>& before, const LatticeMesh& m,
                                const Writes& w) {
    for (std::uint32_t t = 0; t < before.size(); ++t) {
        CAPTURE(t);
        if (m.triangles()[t] != before[t]) REQUIRE(w.slots.count(t) == 1);
    }
    for (const auto s : w.slots) REQUIRE(s < m.triangle_count());
}

// The diamond A (2,0), B (2,4), C (1,2), D (3,2). Diagonal A-B is not
// Delaunay under dx = dy: D is inside the circle of (A, B, C) (centre
// (2, -3.5) in the frame, radius 2.5; D at distance 0.5).
//   t0 = (A, B, C), t1 = (B, A, D): the quad, diagonal as edge 0 of both.
// With `ring`, four outer triangles close each quad side, one per side, so
// the flip's repoints are observable:
//   t2 = (C, B, E) across t0's B-C, t3 = (A, C, F) across t0's C-A,
//   t4 = (D, A, G) across t1's A-D, t5 = (B, D, H) across t1's D-B.
// Quad side masks 0x10 (B-C, constrained), 0x20 (C-A), 0x40 (A-D,
// constrained), 0x80 (D-B).
enum : std::uint32_t { A, B, C, D, E, F, G, H };

LatticeMesh diamond(bool ring, bool diagonal_constrained = false) {
    const std::uint8_t diag = diagonal_constrained ? 0b001 : 0;
    std::vector<LatticeVertex> v{lv(2, 0), lv(2, 4), lv(1, 2), lv(3, 2),
                                 lv(0, 4), lv(0, 0), lv(4, 0), lv(4, 4)};
    std::vector<TriangleIndices> t{{A, B, C}, {B, A, D}};
    std::vector<std::uint8_t> bits{static_cast<std::uint8_t>(diag | 0b010),
                                   static_cast<std::uint8_t>(diag | 0b010)};
    std::vector<Masks> masks{{diagonal_constrained ? 0x1u : 0u, 0x10, 0x20},
                             {diagonal_constrained ? 0x1u : 0u, 0x40, 0x80}};
    if (ring) {
        t.insert(t.end(), {{C, B, E}, {A, C, F}, {D, A, G}, {B, D, H}});
        bits.insert(bits.end(), {0b001, 0b000, 0b001, 0b000});
        masks.insert(masks.end(), {Masks{0x10, 0, 0}, {0x20, 0, 0}, {0x40, 0, 0}, {0x80, 0, 0}});
    } else {
        v.resize(4);
    }
    auto m = LatticeMesh::build(std::move(v), std::move(t), std::move(bits), std::move(masks));
    REQUIRE(m.has_value());
    return *m;
}

const LatticeFrame kSquare{1.0, 1.0};
// GENERATE splits on the commas inside a braced initialiser, so frames are
// generated by index.
const std::array<LatticeFrame, 3> kFrames{LatticeFrame{1.0, 1.0}, LatticeFrame{1.0, 3.0},
                                          LatticeFrame{7.0, 2.0}};

}  // namespace

// ---------------------------------------------------------------- L1

TEST_CASE("L1: flip rewrites both slots and carries the four outer sides", "[lawson][flip]") {
    auto m = diamond(true);
    m.flip(0, 0);

    REQUIRE(m.triangle_count() == 6);
    REQUIRE(m.triangles()[0] == TriangleIndices{C, A, D});
    REQUIRE(m.triangles()[1] == TriangleIndices{C, D, B});

    // All six links of the two flipped slots, and the four outer triangles'
    // links back: t2 and t4 moved slot, t3 and t5 did not.
    REQUIRE(m.neighbours(0) == std::array<std::uint32_t, 3>{3, 4, 1});
    REQUIRE(m.neighbours(1) == std::array<std::uint32_t, 3>{0, 5, 2});
    REQUIRE(m.neighbours(2)[0] == 1);
    REQUIRE(m.neighbours(3)[0] == 0);
    REQUIRE(m.neighbours(4)[0] == 0);
    REQUIRE(m.neighbours(5)[0] == 1);

    // Bits and masks moved with their edges; the new diagonal has neither.
    REQUIRE_FALSE(m.is_constrained(0, 0));
    REQUIRE(m.mask(0, 0) == 0x20);  // C-A
    REQUIRE(m.is_constrained(0, 1));
    REQUIRE(m.mask(0, 1) == 0x40);  // A-D
    REQUIRE_FALSE(m.is_constrained(0, 2));
    REQUIRE(m.mask(0, 2) == 0);     // D-C, the diagonal
    REQUIRE_FALSE(m.is_constrained(1, 0));
    REQUIRE(m.mask(1, 0) == 0);     // C-D, the diagonal
    REQUIRE_FALSE(m.is_constrained(1, 1));
    REQUIRE(m.mask(1, 1) == 0x80);  // D-B
    REQUIRE(m.is_constrained(1, 2));
    REQUIRE(m.mask(1, 2) == 0x10);  // B-C

    check_topology(m);
}

TEST_CASE("L1: flip works from either side and twice restores the diagonal", "[lawson][flip]") {
    auto m = diamond(true);
    m.flip(1, 0);  // from u's side: u = (B, A, D), across is t0 = (A, B, C)
    // t = (B, A, D) with c = D: slot 1 is (D, B, C), slot 0 is (D, C, A).
    REQUIRE(m.triangles()[1] == TriangleIndices{D, B, C});
    REQUIRE(m.triangles()[0] == TriangleIndices{D, C, A});
    check_topology(m);
    REQUIRE(m.constraint_edges().first.size() == diamond(true).constraint_edges().first.size());
}

// ---------------------------------------------------------------- L2

TEST_CASE("L2: a cocircular square is never flipped whichever diagonal it has", "[lawson][ties]") {
    const bool other = GENERATE(false, true);
    const LatticeFrame f = kFrames[GENERATE(0u, 1u, 2u)];
    CAPTURE(other, f.dx, f.dy);
    // (0,0) (1,0) (1,1) (0,1); a rectangle in any frame, so always cocircular.
    auto m = other ? LatticeMesh::build({lv(0, 0), lv(1, 0), lv(1, 1), lv(0, 1)},
                                        {TriangleIndices{0, 1, 2}, TriangleIndices{0, 2, 3}},
                                        {0b011, 0b110}, {Masks{1, 2, 0}, {0, 3, 4}})
                   : LatticeMesh::build({lv(0, 0), lv(1, 0), lv(1, 1), lv(0, 1)},
                                        {TriangleIndices{0, 1, 3}, TriangleIndices{1, 2, 3}},
                                        {0b101, 0b011}, {Masks{1, 0, 4}, {2, 3, 0}});
    REQUIRE(m.has_value());
    const auto before = snapshot(*m);
    Writes w;
    REQUIRE(legalise_all<DefaultKernel>(*m, f, w.sink()) == 0);
    REQUIRE(snapshot(*m) == before);
    REQUIRE(w.slots.empty());
    const std::array<std::uint32_t, 2> seeds{0, 1};
    for (std::uint32_t q = 0; q < 4; ++q)
        REQUIRE(legalise_around<DefaultKernel>(*m, q, std::span<const std::uint32_t>{seeds}, f,
                                               w.sink())
                == 0);
    REQUIRE(snapshot(*m) == before);
}

namespace {

// A 17 x 17 lattice triangulated strip by strip, each strip a seeded random
// monotone merge of its two rows. Valid, CCW, and full of needles, so it is
// far from Delaunay. Built without flip, so the fixture cannot inherit a flip
// bug.
LatticeMesh scrambled_lattice(std::uint32_t n, std::uint32_t seed) {
    std::vector<LatticeVertex> v;
    for (std::uint32_t r = 0; r < n; ++r)
        for (std::uint32_t c = 0; c < n; ++c) v.push_back(lv(r, c));
    std::mt19937 gen{seed};
    std::vector<TriangleIndices> t;
    for (std::uint32_t r = 0; r + 1 < n; ++r) {
        std::uint32_t j = 0, k = 0;  // top column, bottom column
        const auto top = [&](std::uint32_t c) { return r * n + c; };
        const auto bot = [&](std::uint32_t c) { return (r + 1) * n + c; };
        while (j + 1 < n || k + 1 < n) {
            const bool advance_top = k + 1 == n || (j + 1 < n && (gen() & 1u) != 0);
            if (advance_top) {
                t.push_back({top(j), bot(k), top(j + 1)});
                ++j;
            } else {
                t.push_back({top(j), bot(k), bot(k + 1)});
                ++k;
            }
        }
    }
    const auto count = t.size();
    auto m = LatticeMesh::build(std::move(v), std::move(t), std::vector<std::uint8_t>(count, 0),
                                std::vector<Masks>(count, Masks{0, 0, 0}));
    REQUIRE(m.has_value());
    return *m;
}

}  // namespace

TEST_CASE("L2: a scrambled 17 x 17 lattice legalises to Delaunay in bounded flips", "[lawson][ties]") {
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    const LatticeFrame f = kFrames[GENERATE(0u, 1u)];
    CAPTURE(seed, f.dx, f.dy);
    const std::uint32_t n = 17;
    auto m = scrambled_lattice(n, seed);
    REQUIRE(delaunay_violations(m, f) > 0);  // the fixture has work in it
    const auto before = snapshot(m);
    Writes w;
    const std::size_t flips = legalise_all<DefaultKernel>(m, f, w.sink());
    CAPTURE(flips);
    // Lawson: an edge flipped out never returns, so flips <= point pairs. A
    // flip on Cocircular breaks this bound by cycling (R3, R4).
    REQUIRE(flips > 0);
    REQUIRE(flips <= std::size_t{n * n} * (n * n - 1) / 2);
    REQUIRE(m.triangle_count() == before.size());
    REQUIRE(m.vertices().size() == std::size_t{n} * n);
    REQUIRE(delaunay_violations(m, f) == 0);
    check_topology(m);
    check_writes_cover_changes(before, m, w);
    REQUIRE(w.slots.size() <= 2 * flips);

    // Idempotent: a Delaunay mesh has nothing left to flip, ties included.
    const auto done = snapshot(m);
    REQUIRE(legalise_all<DefaultKernel>(m, f, Writes{}.sink()) == 0);
    REQUIRE(snapshot(m) == done);
}

// ---------------------------------------------------------------- L3

TEST_CASE("L3: a constrained diagonal is never flipped", "[lawson][constraints]") {
    const bool constrained = GENERATE(false, true);
    const bool around = GENERATE(false, true);
    CAPTURE(constrained, around);
    auto m = diamond(false, constrained);
    Writes w;
    const std::array<std::uint32_t, 1> seeds{0};
    const std::size_t flips =
        around ? legalise_around<DefaultKernel>(m, C, std::span<const std::uint32_t>{seeds}, kSquare,
                                                w.sink())
               : legalise_all<DefaultKernel>(m, kSquare, w.sink());
    if (constrained) {
        REQUIRE(flips == 0);
        REQUIRE(m.triangles()[0] == TriangleIndices{A, B, C});
        REQUIRE(m.triangles()[1] == TriangleIndices{B, A, D});
        REQUIRE(m.is_constrained(0, 0));
        REQUIRE(m.mask(0, 0) == 0x1);
        REQUIRE(w.slots.empty());
    } else {
        REQUIRE(flips == 1);
        if (around) {  // seeded from slot 0 across its edge 0, so flip(0, 0)
            REQUIRE(m.triangles()[0] == TriangleIndices{C, A, D});
            REQUIRE(m.triangles()[1] == TriangleIndices{C, D, B});
        }
        REQUIRE(m.triangle_count() == 2);
        for (const auto& tri : m.triangles())  // both halves hold the new diagonal C-D
            REQUIRE((std::count(tri.begin(), tri.end(), C) == 1 && std::count(tri.begin(), tri.end(), D) == 1));
        REQUIRE(w.slots == std::set<std::uint32_t>{0, 1});
        REQUIRE(delaunay_violations(m, kSquare) == 0);
    }
    check_topology(m);
}

TEST_CASE("L3: a boundary edge is never a candidate", "[lawson][constraints]") {
    // One triangle: three boundary edges, nothing to flip, whatever q is.
    auto m = LatticeMesh::build({lv(0, 0), lv(5, 0), lv(1, 9)}, {TriangleIndices{0, 1, 2}}, {0},
                                {Masks{0, 0, 0}});
    REQUIRE(m.has_value());
    const std::array<std::uint32_t, 1> seeds{0};
    for (std::uint32_t q = 0; q < 3; ++q)
        REQUIRE(legalise_around<DefaultKernel>(*m, q, std::span<const std::uint32_t>{seeds}, kSquare,
                                               Writes{}.sink())
                == 0);
    REQUIRE(legalise_all<DefaultKernel>(*m, kSquare, Writes{}.sink()) == 0);
}

// ---------------------------------------------------------------- L4

TEST_CASE("L4: the circle test runs in the scaled frame", "[lawson][frame]") {
    // Kite A (1,0), B (1,4), C (0,2), D (2,2), diagonal C-D. In (col, -row)
    // the half-diagonals are 2 across and 1 down, so C-D is Delaunay (strict).
    // With dy = 3 the down half is 3 > 2 and C-D must flip to A-B. With
    // dx = 3 it is 6 across and C-D stays: a dx/dy swap cannot pass.
    enum : std::uint32_t { a, b, c, d };
    const auto build = [] {
        auto m = LatticeMesh::build({lv(1, 0), lv(1, 4), lv(0, 2), lv(2, 2)},
                                    {TriangleIndices{d, c, a}, TriangleIndices{c, d, b}}, {0, 0},
                                    {Masks{0, 0, 0}, {0, 0, 0}});
        REQUIRE(m.has_value());
        return *m;
    };
    SECTION("raw lattice frame: nothing to flip") {
        auto m = build();
        REQUIRE(legalise_all<DefaultKernel>(m, LatticeFrame{1.0, 1.0}, Writes{}.sink()) == 0);
    }
    SECTION("dx = 3: nothing to flip") {
        auto m = build();
        REQUIRE(legalise_all<DefaultKernel>(m, LatticeFrame{3.0, 1.0}, Writes{}.sink()) == 0);
    }
    SECTION("dy = 3: C-D flips to A-B") {
        const bool around = GENERATE(false, true);
        CAPTURE(around);
        auto m = build();
        const LatticeFrame f{1.0, 3.0};
        const std::array<std::uint32_t, 1> seeds{0};
        const std::size_t flips =
            around ? legalise_around<DefaultKernel>(m, a, std::span<const std::uint32_t>{seeds}, f,
                                                    Writes{}.sink())
                   : legalise_all<DefaultKernel>(m, f, Writes{}.sink());
        REQUIRE(flips == 1);
        REQUIRE(delaunay_violations(m, f) == 0);
        check_topology(m);
        bool has_ab = false;
        for (const auto& tri : m.triangles())
            for (unsigned k = 0; k < 3; ++k)
                has_ab = has_ab || (tri[k] == a && tri[(k + 1) % 3] == b)
                      || (tri[k] == b && tri[(k + 1) % 3] == a);
        REQUIRE(has_ab);
    }
}

// ---------------------------------------------------------------- L4b

TEST_CASE("L4b: an edge illegal only from the side that is collinear in the frame still flips",
          "[lawson][frame]") {
    // Guards c23583b's must_flip directly. Orientation is exact on (col, -row);
    // the circle test runs in the LatticeFrame (col * dx, -(row * dy)), which
    // rounds. C sits one ulp of row off the diagonal A-B, so t0 = (A, B, C) is
    // a counter-clockwise sliver in the mesh but collinear in the frame at
    // dx = 3, dy = 1.5, and the kernel answers Cocircular there. From t1 =
    // (B, A, D) the edge is plainly illegal: C lies on the chord A-B, strictly
    // inside the circle through B, A, D. Slot 0 holds the collinear side, so
    // legalise_all tests A-B only from t0, and legalise_around from C tests it
    // only from t0 as well. The pre-fix one-sided test flips nothing in both.
    //   A (col 0, row 0), B (col 1, row 1), D (col 0, row 1),
    //   C (col 0.34, row = the double just below 0.34).
    enum : std::uint32_t { a, b, c, d };
    const LatticeFrame f{3.0, 1.5};
    const std::vector<MeshVertex> v{lv(0, 0), lv(1, 1), MeshVertex{0.34, std::nextafter(0.34, 0.0)},
                                    lv(1, 0)};
    // Preconditions: without both, the case tests nothing.
    REQUIRE(orient_sign(v[a], v[b], v[c]) > 0);
    REQUIRE(DefaultKernel::orient2d(at(f, v[a]), at(f, v[b]), at(f, v[c])) == Orientation::Collinear);
    REQUIRE(DefaultKernel::orient2d(at(f, v[b]), at(f, v[a]), at(f, v[d])) == Orientation::CounterClockwise);
    REQUIRE(DefaultKernel::incircle(at(f, v[b]), at(f, v[a]), at(f, v[d]), at(f, v[c])) == Incircle::Inside);

    const bool around = GENERATE(false, true);
    CAPTURE(around);
    auto m = LatticeMesh::build(v, {TriangleIndices{a, b, c}, TriangleIndices{b, a, d}}, {0, 0},
                                {Masks{0, 0, 0}, {0, 0, 0}});
    REQUIRE(m.has_value());
    REQUIRE(delaunay_violations(*m, f) == 1);  // seen from t1 only
    const std::array<std::uint32_t, 1> seeds{0};
    const std::size_t flips =
        around ? legalise_around<DefaultKernel>(*m, c, std::span<const std::uint32_t>{seeds}, f, Writes{}.sink())
               : legalise_all<DefaultKernel>(*m, f, Writes{}.sink());
    REQUIRE(flips == 1);
    REQUIRE(delaunay_violations(*m, f) == 0);
    // The diagonal is now C-D, and both triangles are counter-clockwise in the mesh.
    bool has_cd = false;
    for (const auto& tri : m->triangles()) {
        REQUIRE(orient_sign(v[tri[0]], v[tri[1]], v[tri[2]]) > 0);
        for (unsigned k = 0; k < 3; ++k)
            has_cd = has_cd || (tri[k] == c && tri[(k + 1) % 3] == d) || (tri[k] == d && tri[(k + 1) % 3] == c);
    }
    REQUIRE(has_cd);
}

// ---------------------------------------------------------------- L5

TEST_CASE("L5: legalise_around after a split restores Delaunay", "[lawson][around]") {
    // The square (0,0) (3,0) (3,3) (0,3) at dx = 10, dy = 5, fanned at (2,1).
    // In the frame (0,3) is strictly inside the circle of (3,3) (0,0) (2,1):
    // centre (27.5, 17.5), radius^2 = 1062.5; (30, 0) is at 312.5. So at
    // least one flip is required.
    auto m = LatticeMesh::build({lv(0, 0), lv(3, 0), lv(3, 3), lv(0, 3)},
                                {TriangleIndices{0, 1, 2}, TriangleIndices{0, 2, 3}}, {0b011, 0b110},
                                {Masks{8, 4, 0}, {0, 2, 1}});
    REQUIRE(m.has_value());
    const LatticeFrame f{10.0, 5.0};
    REQUIRE(delaunay_violations(*m, f) == 0);
    const auto q = m->split_inside(0, lv(2, 1));
    const std::array<std::uint32_t, 3> seeds{0, 2, 3};
    const auto before = snapshot(*m);
    Writes w;
    const std::size_t flips =
        legalise_around<DefaultKernel>(*m, q, std::span<const std::uint32_t>{seeds}, f, w.sink());
    REQUIRE(flips >= 1);
    REQUIRE(m->triangle_count() == 4);
    REQUIRE(delaunay_violations(*m, f) == 0);
    check_topology(*m);
    check_writes_cover_changes(before, *m, w);
    // The perimeter is still constrained with its masks.
    REQUIRE(m->constraint_edges().first.size() == 4);
}

TEST_CASE("L5: legalise_around after an edge split on a constraint keeps both halves constrained",
          "[lawson][around]") {
    // A kite whose constrained diagonal (1,0)-(1,6) is split at (1,3); the
    // spokes to (0,3) and (4,3) are free, and nothing may flip the halves.
    auto m = LatticeMesh::build({lv(1, 0), lv(1, 6), lv(0, 3), lv(4, 3)},
                                {TriangleIndices{0, 1, 2}, TriangleIndices{1, 0, 3}}, {0b001, 0b001},
                                {Masks{0x5, 0, 0}, {0x5, 0, 0}});
    REQUIRE(m.has_value());
    const auto q = m->split_edge(0, 0, lv(1, 3));
    const std::array<std::uint32_t, 4> seeds{0, 1, 2, 3};
    legalise_around<DefaultKernel>(*m, q, std::span<const std::uint32_t>{seeds}, kSquare,
                                   Writes{}.sink());
    check_topology(*m);
    const auto [edges, masks] = m->constraint_edges();
    REQUIRE(edges.size() == 2);
    for (const auto mask : masks) REQUIRE(mask == 0x5);
    for (const auto& e : edges) REQUIRE((e[0] == q || e[1] == q));
    REQUIRE(delaunay_violations(*m, kSquare) == 0);
}
