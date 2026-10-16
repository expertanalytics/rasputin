// Increment 20c, PR 20c-1 (docs/increments/20c-soft-quality.md, R1 and "Tests
// @tester writes red first", 20c-1, CF1): the geometric helper every insertion
// path asks before it inserts a point near a constraint. INVARIANT-CRITICAL
// (the design names test_constraint_foot): @reviewer mutation-tests it.
//
// Interface, as R1 names it (include/terrain/mesh/constraint_foot.hpp):
//
//   enum class FootStatus : std::uint8_t { None, Hit, NearEnd, NotCounterClockwise };
//   struct FootSearch { FootStatus status; std::uint32_t owner; unsigned edge; MeshVertex at; };
//   FootSearch constraint_foot(const LatticeMesh&, std::uint32_t t, MeshVertex p,
//                              double delta, const LatticeFrame&);
//
// PINNED HERE, where R1 leaves it open (listed for @architect in the handback):
//   - the namespace is terrain::mesh, as the header's directory says;
//   - "closer than delta" is strict (a distance equal to delta is not near),
//     and the distance is to the closed segment in the world frame f
//     (col * dx, row * dy), as 20b's foot_of measured it;
//   - the foot `at` is the orthogonal projection of p on the edge's line in
//     that world frame, returned in (col, row);
//   - owner and edge are asserted only for a Hit; `at` only for a Hit, as R1
//     says.
//
// Search order (R1): t's edges in edge order, then, for each of t's
// unconstrained edges in edge order, the neighbour's other two edges; the
// first constrained, non-frozen edge closer than delta decides. The cases
// below never put two candidate edges in one neighbour, so the order inside a
// neighbour is not pinned.
//
// Fixtures are hand-built LatticeMesh objects in (col, row); the cases written
// in world (x, y) go through V(x, y) = (x, 20 - y), so "above" in the comments
// is the world's. Distances and projections are recomputed here from the
// stored coordinates with this file's own formulas, never read from the
// helper. Every Hit/None table states its scale: delta = 0.5 world units in a
// frame of dx = dy = 1 (cells), the largest coordinate 20 cells.
//
// Mutants each case is meant to kill (design, "Mutants to kill", CF1):
//   M-end   the end check removed       NearEnd cases, own and neighbour, and
//                                       the foot a hair (1e-9 cells) from a vertex
//   M-cross the neighbour search crosses a constrained edge
//                                       "None across a constrained edge" (on it
//                                       exactly, and frozen)
//   M-mid   the midpoint instead of the foot
//                                       every `at` check, the tilted world edge
//   and, beyond the design's list: a strict/non-strict delta flip (the point
//   exactly delta from the edge; 0.4 and 0.6 alone do not tell the two apart),
//   a row/col swap in the world distance (the dx = 10, dy = 5 case), the
//   frozen test dropped, and "nearest edge" instead of "first in order".

#include <catch2/catch_test_macros.hpp>

#include <terrain/mesh/constraint_foot.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <utility>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::constraint_foot;
using terrain::mesh::foot_reachable;
using terrain::mesh::FootSearch;
using terrain::mesh::FootStatus;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::MeshVertex;
using terrain::mesh::orient_sign;

namespace {

constexpr double kDelta = 0.5;
constexpr double kAt = 1e-12;  // cells: "the foot is the orthogonal projection (to 1e-12 cells)"

using Pair = std::pair<std::uint32_t, std::uint32_t>;

MeshVertex V(double x, double y) { return MeshVertex{x, 20.0 - y}; }

// A mesh from vertices, triangles and {undirected vertex pair -> mask} for
// the constrained edges.
LatticeMesh build(const std::vector<MeshVertex>& v, const std::vector<TriangleIndices>& tris,
                  const std::map<Pair, std::uint32_t>& constrained) {
    std::vector<std::uint8_t> bits(tris.size(), 0);
    std::vector<std::array<std::uint32_t, 3>> masks(tris.size(), {0, 0, 0});
    for (std::size_t t = 0; t < tris.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto it = constrained.find(std::minmax(tris[t][k], tris[t][(k + 1) % 3]));
                it != constrained.end()) {
                bits[t] |= static_cast<std::uint8_t>(1u << k);
                masks[t][k] = it->second;
            }
    auto m = LatticeMesh::build(v, tris, std::move(bits), std::move(masks));
    REQUIRE(m.has_value());
    return std::move(*m);
}

// World-frame distance from p to the closed segment a-b, and the orthogonal
// projection on its line, in (col, row).
struct Projection {
    double distance;
    double along;  // world distance from a to the projection
    MeshVertex at;
};

Projection project(MeshVertex a, MeshVertex b, MeshVertex p, const LatticeFrame& f) {
    const double ux = (b.col - a.col) * f.dx, uy = (b.row - a.row) * f.dy;
    const double px = (p.col - a.col) * f.dx, py = (p.row - a.row) * f.dy;
    const double len = std::hypot(ux, uy);
    const double s = (px * ux + py * uy) / (len * len);
    return Projection{std::abs(px * uy - py * ux) / len, s * len,
                      MeshVertex{a.col + s * (b.col - a.col), a.row + s * (b.row - a.row)}};
}

void require_hit(const FootSearch& r, std::uint32_t owner, unsigned edge, MeshVertex at) {
    REQUIRE(r.status == FootStatus::Hit);
    REQUIRE(r.owner == owner);
    REQUIRE(r.edge == edge);
    CAPTURE(r.at.col, r.at.row, at.col, at.row);
    REQUIRE(std::abs(r.at.col - at.col) <= kAt);
    REQUIRE(std::abs(r.at.row - at.row) <= kAt);
}

// ------------------------------------------------------------------ fixture A
//
// World: a (0, 0), b (20, 0), c (10, 10); t = (a, b, c) in slot 0 above a-b,
// and u = (b, a, d) in slot 1 below it, a sliver with d (10, -0.01). t's
// edges: 0 a-b, 1 b-c, 2 c-a. u's: 0 b-a, 1 a-d, 2 d-b. A point (5, h) in t
// is h from a-b and about h + 0.005 from a-d, its projection on either about
// 5 from a: far from both ends.
struct FixtureA {
    std::vector<MeshVertex> v{V(0, 0), V(20, 0), V(10, 10), V(10, -0.01)};
    std::vector<TriangleIndices> tris{{0, 1, 2}, {1, 0, 3}};
};

const LatticeFrame kUnit{1.0, 1.0};

}  // namespace

// ------------------------------------------------------------------ own edge

TEST_CASE("CF1: a point 0.01 and 0.4 cells from a constrained edge of its triangle is a Hit there; 0.6 is None",
          "[mesh][constraint_foot][CF1]") {
    const FixtureA fx;
    const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 1u}});  // a-b only
    for (const double h : {0.01, 0.4}) {
        CAPTURE(h);
        const MeshVertex p = V(5, h);
        const auto pr = project(fx.v[0], fx.v[1], p, kUnit);
        REQUIRE(std::abs(pr.distance - h) <= 1e-12);  // the fixture is what it claims
        require_hit(constraint_foot(m, 0, p, kDelta, kUnit), 0, 0, pr.at);
    }
    REQUIRE(constraint_foot(m, 0, V(5, 0.6), kDelta, kUnit).status == FootStatus::None);
}

TEST_CASE("CF1: a point exactly delta from a constrained edge of its triangle is None (closer than delta is strict)",
          "[mesh][constraint_foot][CF1]") {
    // V(5, 0.5) is (5, 19.5): every coordinate and the 0.5 cells to a-b are
    // exact in doubles, so the helper's distance is delta itself, not a
    // rounding of it. Kills ">= delta" made "> delta".
    const FixtureA fx;
    const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 1u}});  // a-b only
    const MeshVertex p = V(5, 0.5);
    const auto pr = project(fx.v[0], fx.v[1], p, kUnit);
    REQUIRE(pr.distance == kDelta);  // the fixture is what it claims: exactly delta
    REQUIRE(pr.along == 5.0);        // far from both ends: a near point would be a Hit, not NearEnd
    REQUIRE(constraint_foot(m, 0, p, kDelta, kUnit).status == FootStatus::None);
}

TEST_CASE("CF1: the foot is the orthogonal projection, not the edge's midpoint", "[mesh][constraint_foot][CF1]") {
    const FixtureA fx;
    const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 1u}});
    const MeshVertex p = V(3.25, 0.1);  // projection (3.25, 0); the midpoint is (10, 0)
    const auto r = constraint_foot(m, 0, p, kDelta, kUnit);
    require_hit(r, 0, 0, V(3.25, 0));
    REQUIRE(std::abs(r.at.col - 10.0) > 1.0);  // M-mid
}

// ------------------------------------------------------------------ the neighbour

TEST_CASE("CF1: a point 0.01 and 0.4 cells from a neighbour's constrained edge is a Hit owned by the neighbour; 0.6 is None",
          "[mesh][constraint_foot][CF1]") {
    const FixtureA fx;
    const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 3u}, 2u}});  // a-d only; a-b is free
    REQUIRE(m.neighbours(0)[0] == 1u);
    for (const double want : {0.01, 0.4, 0.6}) {
        CAPTURE(want);
        // (5, h) is h + 0.005 from a-d, up to the slope's cosine.
        const MeshVertex p = V(5, want - 0.005);
        const auto pr = project(fx.v[0], fx.v[3], p, kUnit);
        REQUIRE(std::abs(pr.distance - want) <= 1e-5);
        REQUIRE(pr.along > 4.9);
        const auto r = constraint_foot(m, 0, p, kDelta, kUnit);
        if (want < kDelta)
            require_hit(r, 1, 1, pr.at);
        else
            REQUIRE(r.status == FootStatus::None);
    }
}

TEST_CASE("CF1: t's own constrained edges come first, in edge order, however near a later one is",
          "[mesh][constraint_foot][CF1]") {
    // t = (a, b, c') with c' (10, 0.8): a sliver over a-b. u as in fixture A.
    std::vector<MeshVertex> v{V(0, 0), V(20, 0), V(10, 0.8), V(10, -0.01)};
    const std::vector<TriangleIndices> tris{{0, 1, 2}, {1, 0, 3}};

    SECTION("own edge c'-a (edge 2) before the neighbour's a-d, which is nearer") {
        const LatticeMesh m = build(v, tris, {{{0u, 2u}, 4u}, {{0u, 3u}, 2u}});
        const MeshVertex p = V(5, 0.1);
        const auto own = project(v[2], v[0], p, kUnit), nbr = project(v[0], v[3], p, kUnit);
        REQUIRE(own.distance > nbr.distance);  // the fixture is what it claims: 0.299 against 0.105
        REQUIRE(own.distance < kDelta);
        require_hit(constraint_foot(m, 0, p, kDelta, kUnit), 0, 2, own.at);
    }
    SECTION("own edge a-b (edge 0) before own edge c'-a (edge 2), which is nearer") {
        const LatticeMesh m = build(v, tris, {{{0u, 1u}, 1u}, {{0u, 2u}, 4u}});
        const MeshVertex p = V(5, 0.35);
        const auto ab = project(v[0], v[1], p, kUnit), ca = project(v[2], v[0], p, kUnit);
        REQUIRE(ab.distance > ca.distance);  // 0.35 against 0.05
        require_hit(constraint_foot(m, 0, p, kDelta, kUnit), 0, 0, ab.at);
    }
}

// ------------------------------------------------------------------ NearEnd

TEST_CASE("CF1: a foot within delta of either end of its edge is NearEnd, on t's edge and on a neighbour's",
          "[mesh][constraint_foot][CF1]") {
    const FixtureA fx;
    SECTION("own edge a-b") {
        const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 1u}});
        REQUIRE(constraint_foot(m, 0, V(0.3, 0.1), kDelta, kUnit).status == FootStatus::NearEnd);   // 0.3 from a
        REQUIRE(constraint_foot(m, 0, V(19.8, 0.1), kDelta, kUnit).status == FootStatus::NearEnd);  // 0.2 from b
        // A foot a hair (1e-9 cells) from a, p inside t (below c-a, y < x).
        REQUIRE(constraint_foot(m, 0, V(1e-9, 5e-10), kDelta, kUnit).status == FootStatus::NearEnd);
        require_hit(constraint_foot(m, 0, V(0.6, 0.1), kDelta, kUnit), 0, 0, V(0.6, 0));            // the control
    }
    SECTION("the neighbour's a-d") {
        const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 3u}, 2u}});
        const MeshVertex near_a = V(0.3, 0.001);
        REQUIRE(project(fx.v[0], fx.v[3], near_a, kUnit).distance < kDelta);
        REQUIRE(constraint_foot(m, 0, near_a, kDelta, kUnit).status == FootStatus::NearEnd);
    }
}

// ------------------------------------------------------------------ None

TEST_CASE("CF1: a point exactly on the constrained edge is None", "[mesh][constraint_foot][CF1]") {
    const FixtureA fx;
    const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 1u}});
    const MeshVertex p = V(5, 0);
    REQUIRE(orient_sign(fx.v[0], fx.v[1], p) == 0);
    REQUIRE(constraint_foot(m, 0, p, kDelta, kUnit).status == FootStatus::None);
}

TEST_CASE("CF1: no search across a constrained edge", "[mesh][constraint_foot][CF1]") {
    // p lies exactly on a-b, so a-b itself is skipped ("p not exactly on it");
    // u's a-d is 0.1 from p. Across a free a-b that is a Hit in u (the
    // control); across a constrained a-b it must be None (M-cross).
    const FixtureA fx;
    const MeshVertex p = V(5, 0);
    REQUIRE(project(fx.v[0], fx.v[3], p, kUnit).distance < kDelta);
    SECTION("a-b free: the control") {
        const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 3u}, 2u}});
        require_hit(constraint_foot(m, 0, p, kDelta, kUnit), 1, 1, project(fx.v[0], fx.v[3], p, kUnit).at);
    }
    SECTION("a-b constrained") {
        const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 1u}, {{0u, 3u}, 2u}});
        REQUIRE(constraint_foot(m, 0, p, kDelta, kUnit).status == FootStatus::None);
    }
    SECTION("a-b constrained and frozen, p 0.1 off it") {
        // a-b is near but frozen, so skipped; it is still a constraint, so the
        // search does not cross it to u's a-d (0.105 away).
        LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 4u}, {{0u, 3u}, 2u}});
        m.set_frozen_mask(4u);
        REQUIRE(constraint_foot(m, 0, V(5, 0.1), kDelta, kUnit).status == FootStatus::None);
    }
}

TEST_CASE("CF1: a frozen edge is never a foot, on t or on a neighbour", "[mesh][constraint_foot][CF1]") {
    const FixtureA fx;
    SECTION("own edge") {
        LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 4u}});
        require_hit(constraint_foot(m, 0, V(5, 0.1), kDelta, kUnit), 0, 0, V(5, 0));  // unfrozen: the control
        m.set_frozen_mask(4u);
        REQUIRE(constraint_foot(m, 0, V(5, 0.1), kDelta, kUnit).status == FootStatus::None);
    }
    SECTION("the neighbour's edge") {
        LatticeMesh m = build(fx.v, fx.tris, {{{0u, 3u}, 4u}});
        m.set_frozen_mask(4u);
        REQUIRE(constraint_foot(m, 0, V(5, 0.1), kDelta, kUnit).status == FootStatus::None);
    }
}

// ------------------------------------------------------------------ foot_reachable

// foot_reachable(m, t) is the contract the speed fix 69f37d1c rests on: where
// it is false, constraint_foot finds nothing at t. Kills "ignores the
// triangle's own edges", "ignores the neighbours" and "always true".
TEST_CASE("CF1: foot_reachable is true where t's own edge or a neighbour's is constrained and live, and false "
          "where constraint_foot can find nothing",
          "[mesh][constraint_foot][CF1][reachable]") {
    const FixtureA fx;
    const MeshVertex p = V(5, 0.1);  // 0.1 from a-b, 0.105 from a-d
    SECTION("only t's own a-b: true") {
        const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 1u}});
        REQUIRE(foot_reachable(m, 0));
        REQUIRE(constraint_foot(m, 0, p, kDelta, kUnit).status == FootStatus::Hit);
    }
    SECTION("only the neighbour's a-d, across t's free a-b: true") {
        const LatticeMesh m = build(fx.v, fx.tris, {{{0u, 3u}, 2u}});
        for (unsigned k = 0; k < 3; ++k) REQUIRE_FALSE(m.is_constrained(0, k));
        REQUIRE(foot_reachable(m, 0));
        REQUIRE(constraint_foot(m, 0, p, kDelta, kUnit).status == FootStatus::Hit);
    }
    SECTION("the only live edge, a-d, lies beyond t's constrained (frozen) a-b: false") {
        LatticeMesh m = build(fx.v, fx.tris, {{{0u, 1u}, 4u}, {{0u, 3u}, 2u}});
        m.set_frozen_mask(4u);
        REQUIRE_FALSE(m.is_frozen(1, 1));  // a-d is live
        REQUIRE_FALSE(foot_reachable(m, 0));
        REQUIRE(constraint_foot(m, 0, p, kDelta, kUnit).status == FootStatus::None);
    }
    SECTION("the only constrained edge is frozen: false, on t and on the neighbour") {
        for (const Pair& e : {Pair{0u, 1u}, Pair{0u, 3u}}) {
            CAPTURE(e.first, e.second);
            LatticeMesh m = build(fx.v, fx.tris, {{e, 4u}});
            m.set_frozen_mask(4u);
            REQUIRE_FALSE(foot_reachable(m, 0));
            REQUIRE(constraint_foot(m, 0, p, kDelta, kUnit).status == FootStatus::None);
        }
    }
}

// ------------------------------------------------------------------ world distance

TEST_CASE("CF1: distance and foot are measured in the world frame, dx = 10 and dy = 5",
          "[mesh][constraint_foot][CF1][frame]") {
    // One triangle (col, row): V0 (0, 0), V1 (0, 20), V2 (20, 20); edges 0
    // V0-V1 (a column line), 1 V1-V2 (a row line), 2 V2-V0 (the diagonal).
    // delta = min(dx, dy) / 2 = 2.5 m.
    const LatticeFrame f{10.0, 5.0};
    const double delta = 2.5;
    const std::vector<MeshVertex> v{{0, 0}, {0, 20}, {20, 20}};
    SECTION("0.3 columns is 3 m, None; 0.3 rows is 1.5 m, a Hit") {
        const LatticeMesh m = build(v, {{0, 1, 2}}, {{{0u, 1u}, 1u}, {{1u, 2u}, 2u}});
        REQUIRE(constraint_foot(m, 0, MeshVertex{0.3, 10}, delta, f).status == FootStatus::None);
        require_hit(constraint_foot(m, 0, MeshVertex{10, 19.7}, delta, f), 0, 1, MeshVertex{10, 20});
    }
    SECTION("on the diagonal the foot is the world projection, not the lattice one") {
        const LatticeMesh m = build(v, {{0, 1, 2}}, {{{0u, 2u}, 4u}});
        const MeshVertex p{10, 10.2};
        const auto world = project(v[2], v[0], p, f);
        const auto lattice = project(v[2], v[0], p, kUnit);
        REQUIRE(world.distance < delta);
        REQUIRE(std::abs(world.at.col - lattice.at.col) > 0.01);  // the two differ: (10.04, ...) against (10.1, ...)
        require_hit(constraint_foot(m, 0, p, delta, f), 0, 2, world.at);
    }
}

// ------------------------------------------------------------------ NotCounterClockwise

TEST_CASE("CF1: a foot whose split would leave a child not strictly counter-clockwise is NotCounterClockwise",
          "[mesh][constraint_foot][CF1]") {
    // a-b is constrained and shared by two slivers, t = (a, b, c) above and
    // s = (b, a, d) below, whose apexes c and d are the doubles nearest the
    // line a-b on either side within 64 ulps of p's foot (found with exact
    // rational arithmetic; no double in that window lies on the line). So
    // EVERY double the foot can round to splits one of the four children flat
    // or folded: checked below for the whole window. p (0.2 from a-b, its foot
    // 4.45 from either end) lies in w1 = (a, c, e) and reaches a-b through t,
    // a neighbour of w1 across the free edge a-c.
    const MeshVertex a{1.25, 6.5}, b{9.1, 2.3}, c{5.1750000000000504, 4.399999999999973},
        d{5.174999999999949, 4.400000000000027}, e{3.76, 1.755};
    const MeshVertex p{5.080649239959596, 4.22365393659115};
    const LatticeMesh m = build({a, b, c, d, e}, {{0, 1, 2}, {1, 0, 3}, {0, 2, 4}, {2, 1, 4}}, {{{0u, 1u}, 1u}});
    const std::uint32_t w1 = 2;
    for (unsigned k = 0; k < 3; ++k) REQUIRE(orient_sign(m.corner(w1, k), m.corner(w1, (k + 1) % 3), p) > 0);
    const auto pr = project(a, b, p, kUnit);
    REQUIRE(pr.distance < kDelta);
    REQUIRE(pr.along > 4.0);

    // The fixture is what it claims: every double within 64 ulps of the
    // projection, in each coordinate, folds a child.
    auto ulps = [](double x, int n) {
        for (; n > 0; --n) x = std::nextafter(x, std::numeric_limits<double>::infinity());
        for (; n < 0; ++n) x = std::nextafter(x, -std::numeric_limits<double>::infinity());
        return x;
    };
    std::size_t fits = 0;
    for (int i = -64; i <= 64; ++i)
        for (int j = -64; j <= 64; ++j) {
            const MeshVertex f{ulps(pr.at.col, i), ulps(pr.at.row, j)};
            REQUIRE(orient_sign(a, b, f) != 0);
            if (orient_sign(a, f, c) > 0 && orient_sign(f, b, c) > 0 && orient_sign(b, f, d) > 0
                && orient_sign(f, a, d) > 0)
                ++fits;
        }
    REQUIRE(fits == 0);

    REQUIRE(constraint_foot(m, w1, p, kDelta, kUnit).status == FootStatus::NotCounterClockwise);
}
