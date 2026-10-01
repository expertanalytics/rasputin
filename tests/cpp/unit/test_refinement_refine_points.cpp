// Increment 15c-1 (docs/increments/15c-geographic-dem.md, D5 and "Tests for
// @tester", RP1, RP2, RP5 to RP8): the final check, `refine_points`, on
// hand-built start meshes and a few placed check points. RP3 and RP4 (the J2
// oracle and determinism) are in property/prop_refinement_refine_points.cpp.
//
// Interface, as D5 names it, with what this suite CHOOSES where D5 is silent:
//
//   PointRefineOptions{tolerance, threads}
//   PointRefineOutcome derives publicly from RefineOutcome (D5: "a RefineOutcome
//        ... plus coincident and coincident_max_error"), so status, vertices,
//        z, valid, triangles, edges, masks, rounds, inserted, max_error,
//        uncovered and carved are read as refine's are
//   refine_points(store, start, z, valid, edges, masks, options)
//   LatticeMesh::split_inside(t, MeshVertex)    D5 step 3: takes an off-node
//        vertex (a LatticeVertex converts to it exactly)
//   an inserted vertex is output at (x_min + col h, y_max - row h), its z the
//        point's own, valid
//
// Positions are dyadic fractions of a cell on an integral frame and heights
// are dyadic too, so the stored position is the given one and the planes are
// exact: equality is asserted where the arithmetic is exact.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/refine_points.hpp>

#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <map>
#include <random>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::LatticeMesh;
using terrain::mesh::MeshVertex;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::raster::RasterGeometry;
using terrain::refinement::CheckPoints;
using terrain::refinement::PointRefineOptions;
using terrain::refinement::PointRefineOutcome;
using terrain::refinement::RefineOutcome;
using terrain::refinement::RefineStatus;
using terrain::refinement::refine_points;

static_assert(std::derived_from<PointRefineOutcome, RefineOutcome>,
              "this suite reads PointRefineOutcome as a RefineOutcome plus two fields");

namespace {

constexpr double kH = 30.0, kX = 1000.0, kY = 2000.0;

RasterGeometry square(std::size_t n) { return RasterGeometry{kX, kY, kH, kH, n, n}; }

Point2 world(double col, double row) { return Point2{kX + col * kH, kY - row * kH}; }

struct Frac {
    double col;
    double row;
};
Frac frac(Point2 p) { return Frac{(p.x - kX) / kH, (kY - p.y) / kH}; }

// A start mesh as numbers, the way phase 1 hands it over.
struct Start {
    IndexedMesh2 mesh;
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    std::vector<std::array<std::uint32_t, 2>> edges;
    std::vector<std::uint32_t> masks;
};

// refinement_fixtures::grid_mesh over `n` nodes a side, z from f(col, row).
template <class F>
Start grid_start(std::size_t n, std::size_t stride, F f) {
    auto s = refinement_fixtures::grid_mesh(square(n), stride);
    Start out{std::move(s.mesh), {}, {}, std::move(s.edges), std::move(s.masks)};
    for (const auto& rc : s.lattice) {
        out.z.push_back(f(static_cast<double>(rc.col), static_cast<double>(rc.row)));
        out.valid.push_back(1);
    }
    return out;
}

// A start mesh from (col, row) vertices and triangles; every listed edge is a
// constraint with its mask. z from f(col, row).
template <class F>
Start hand_start(const std::vector<std::array<double, 2>>& cr, std::vector<TriangleIndices> tris,
                 std::vector<std::pair<std::array<std::uint32_t, 2>, std::uint32_t>> constraints, F f) {
    Start s;
    std::vector<Point2> xy;
    for (const auto& [c, r] : cr) {
        xy.push_back(world(c, r));
        s.z.push_back(f(c, r));
        s.valid.push_back(1);
    }
    std::vector<std::uint8_t> bits(tris.size(), 0);
    s.mesh = IndexedMesh2{std::move(xy), std::move(tris), std::move(bits)};
    for (const auto& [e, m] : constraints) {
        s.edges.push_back(e);
        s.masks.push_back(m);
    }
    return s;
}

struct Pt {
    double col;
    double row;
    float z;
};

CheckPoints store(std::size_t n, const std::vector<Pt>& pts) {
    CheckPoints cp{square(n)};
    std::vector<Point2> xy;
    std::vector<float> z;
    for (const auto& p : pts) {
        xy.push_back(world(p.col, p.row));
        z.push_back(p.z);
    }
    cp.add(std::span<const Point2>{xy}, std::span<const float>{z});
    cp.freeze();
    return cp;
}

PointRefineOutcome run(const CheckPoints& cp, const Start& s, double tol, unsigned threads = 1) {
    PointRefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    return refine_points(cp, s.mesh, std::span<const double>{s.z}, std::span<const std::uint8_t>{s.valid},
                         std::span<const std::array<std::uint32_t, 2>>{s.edges},
                         std::span<const std::uint32_t>{s.masks}, o);
}

// Output constraint edges as {(lo, hi) -> mask}.
std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> constraint_map(const RefineOutcome& out) {
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> m;
    for (std::size_t i = 0; i < out.edges.size(); ++i)
        m[std::minmax(out.edges[i][0], out.edges[i][1])] = out.masks[i];
    return m;
}

bool has_edge(const RefineOutcome& out, std::uint32_t a, std::uint32_t b) {
    for (const auto& t : out.triangles)
        for (unsigned k = 0; k < 3; ++k)
            if (std::minmax(t[k], t[(k + 1) % 3]) == std::minmax(a, b)) return true;
    return false;
}

Point2 lattice_frame(Point2 world_point) {
    const Frac f = frac(world_point);
    return Point2{f.col, -f.row};
}

// p in the closed triangle t of out, by the exact kernel on (col, -row).
bool in_closed(const RefineOutcome& out, const TriangleIndices& t, Point2 p) {
    const Point2 a = lattice_frame(out.vertices[t[0]]), b = lattice_frame(out.vertices[t[1]]),
                 c = lattice_frame(out.vertices[t[2]]), q = lattice_frame(p);
    return DefaultKernel::orient2d(a, b, q) != Orientation::Clockwise
        && DefaultKernel::orient2d(b, c, q) != Orientation::Clockwise
        && DefaultKernel::orient2d(c, a, q) != Orientation::Clockwise;
}

double flat10(double, double) { return 10.0; }

}  // namespace

// --------------------------------------------------------------------------
// split_inside in lattice_mesh.hpp, with an off-node vertex
// --------------------------------------------------------------------------

TEST_CASE("split_inside takes an off-node MeshVertex", "[lattice_mesh][split_inside]") {
    // (0,0), (0,4), (4,4) in (col, row): counter-clockwise in (col, -row).
    auto m = LatticeMesh::build(std::vector<MeshVertex>{{0.0, 0.0}, {0.0, 4.0}, {4.0, 4.0}}, {{0, 1, 2}},
                                {0}, {{0, 0, 0}});
    REQUIRE(m.has_value());
    const MeshVertex p{1.25, 2.5};
    const auto q = m->split_inside(0, p);
    REQUIRE(q == 3);
    REQUIRE(m->vertices()[3] == p);
    REQUIRE(m->triangle_count() == 3);
    REQUIRE(m->triangles()[0] == TriangleIndices{0, 1, 3});
    REQUIRE(m->triangles()[1] == TriangleIndices{1, 2, 3});
    REQUIRE(m->triangles()[2] == TriangleIndices{2, 0, 3});
    for (std::size_t t = 0; t < 3; ++t)
        REQUIRE(terrain::mesh::orient_sign(m->corner(t, 0), m->corner(t, 1), m->corner(t, 2)) > 0);
}

// --------------------------------------------------------------------------
// RP1, a plane
// --------------------------------------------------------------------------

TEST_CASE("RP1: check points on the start mesh's planes insert nothing", "[refine_points][RP1]") {
    auto plane = [](double c, double r) { return 10.0 + 0.5 * c - 0.25 * r; };
    const Start s = grid_start(17, 4, plane);
    std::mt19937 gen{151};
    std::vector<Pt> pts;
    for (int i = 0; i < 300; ++i) {
        const double c = static_cast<double>(gen() % (16u * 1024u)) / 1024.0;
        const double r = static_cast<double>(gen() % (16u * 1024u)) / 1024.0;
        pts.push_back(Pt{c, r, static_cast<float>(plane(c, r))});
    }
    const CheckPoints cp = store(17, pts);
    const auto out = run(cp, s, 1e-9);
    REQUIRE(out.ok());
    REQUIRE(out.message.empty());
    REQUIRE(out.inserted == 0);
    REQUIRE(out.rounds == 1);
    REQUIRE(out.max_error <= 1e-9);
    REQUIRE(out.vertices == std::vector<Point2>(s.mesh.vertices().begin(), s.mesh.vertices().end()));
    REQUIRE(out.z == s.z);
    REQUIRE(out.triangles.size() == s.mesh.triangle_count());
    REQUIRE(out.uncovered == 0);
}

// --------------------------------------------------------------------------
// RP2, one bump
// --------------------------------------------------------------------------

TEST_CASE("RP2: above tolerance is inserted at its own position with its own z; at or below is not",
          "[refine_points][RP2]") {
    const Start s = grid_start(17, 4, flat10);
    const std::size_t n0 = s.mesh.vertices().size();
    // Three triangles far apart, none of the points on an edge or a diagonal.
    const Pt above{2.25, 1.5, 12.0f}, below{13.5, 14.25, 10.5f}, at_tol{2.5, 13.75, 11.0f};
    const CheckPoints cp = store(17, {above, below, at_tol});
    const auto out = run(cp, s, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 1);
    REQUIRE(out.rounds == 2);
    REQUIRE(out.vertices.size() == n0 + 1);
    REQUIRE(out.vertices.back() == world(above.col, above.row));
    REQUIRE(out.z.back() == 12.0);
    REQUIRE(out.valid.back() == 1);
    REQUIRE(out.triangles.size() == s.mesh.triangle_count() + 2);
    // `>` as needs_split: the point exactly at tolerance stays, and is the max.
    REQUIRE(out.max_error == 1.0);
}

TEST_CASE("RP2: one point exactly at tolerance is not inserted", "[refine_points][RP2]") {
    const Start s = grid_start(17, 4, flat10);
    const CheckPoints cp = store(17, {{2.5, 13.75, 11.0f}});
    const auto out = run(cp, s, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 0);
    REQUIRE(out.max_error == 1.0);
}

TEST_CASE("RP2: equal errors in one triangle go to the first point in store order",
          "[refine_points][RP2][tiebreak]") {
    // D5: the worst point wins by strictly larger error, so a tie goes to the
    // first in store order, whatever the thread count. Both points are in the
    // start triangle (0,0) (0,4) (4,4), off its edges, |z - 10| = 2 for both.
    // Store order is (cell row, ...): `first` (cell row 2) before `second`
    // (cell row 3), though `second` is added first and has the smaller z.
    const Start s = grid_start(17, 4, flat10);
    const std::size_t n0 = s.mesh.vertices().size();
    const Pt first{1.25, 2.5, 12.0f}, second{0.75, 3.5, 8.0f};
    for (const unsigned threads : {1u, 4u}) {
        CAPTURE(threads);
        const auto out = run(store(17, {second, first}), s, 1.0, threads);
        REQUIRE(out.ok());
        REQUIRE(out.inserted >= 1);
        REQUIRE(out.vertices[n0] == world(first.col, first.row));
        REQUIRE(out.z[n0] == 12.0);
    }
}

TEST_CASE("refine_points refuses a store that was never frozen", "[refine_points][refusal]") {
    const Start s = grid_start(9, 4, flat10);
    CheckPoints cp{square(9)};
    const std::vector<Point2> xy{world(1.5, 1.5)};
    const std::vector<float> z{12.0f};
    cp.add(std::span<const Point2>{xy}, std::span<const float>{z});
    REQUIRE_THROWS_MATCHES(run(cp, s, 1.0), std::logic_error,
                           Catch::Matchers::MessageMatches(Catch::Matchers::ContainsSubstring("not frozen")));
}

// --------------------------------------------------------------------------
// RP5, constraints and edges
// --------------------------------------------------------------------------

namespace {
// A square (0,0) (0,8) (8,8) (8,0) in (col, row), its diagonal (0,0)-(8,8)
// constrained with mask 16 and its sides with mask 1.
Start diagonal_square() {
    return hand_start({{0, 0}, {0, 8}, {8, 8}, {8, 0}}, {{0, 1, 2}, {0, 2, 3}},
                      {{{0, 1}, 1}, {{1, 2}, 1}, {{2, 3}, 1}, {{3, 0}, 1}, {{0, 2}, 16}},
                      [](double, double) { return 0.0; });
}
}  // namespace

TEST_CASE("RP5: a point exactly on a constrained edge splits it; both halves keep bit and mask",
          "[refine_points][RP5]") {
    const Start s = diagonal_square();
    const CheckPoints cp = store(9, {{2.5, 2.5, 50.0f}});
    const auto out = run(cp, s, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 1);
    REQUIRE(out.vertices.size() == 5);
    REQUIRE(out.vertices[4] == world(2.5, 2.5));
    REQUIRE(out.triangles.size() == 4);  // 2 -> 4
    const auto cm = constraint_map(out);
    REQUIRE(cm.size() == 6);
    REQUIRE(cm.at({0, 4}) == 16);
    REQUIRE(cm.at({2, 4}) == 16);
    REQUIRE_FALSE(cm.contains({0, 2}));
    REQUIRE(cm.at({0, 1}) == 1);
    REQUIRE(cm.at({1, 2}) == 1);
    REQUIRE(cm.at({2, 3}) == 1);
    REQUIRE(cm.at({0, 3}) == 1);
    REQUIRE(out.max_error == 0.0);
}

TEST_CASE("RP5: a point on the domain's boundary splits 1 -> 2", "[refine_points][RP5]") {
    const Start s = diagonal_square();
    const CheckPoints cp = store(9, {{0.0, 3.5, 50.0f}});  // on side (0,0)-(0,8)
    const auto out = run(cp, s, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 1);
    REQUIRE(out.vertices[4] == world(0.0, 3.5));
    REQUIRE(out.triangles.size() == 3);
    const auto cm = constraint_map(out);
    REQUIRE(cm.at({0, 4}) == 1);
    REQUIRE(cm.at({1, 4}) == 1);
    REQUIRE_FALSE(cm.contains({0, 1}));
    REQUIRE(cm.at({0, 2}) == 16);
}

TEST_CASE("RP5: a constrained edge is never flipped, even where Delaunay wants it",
          "[refine_points][RP5]") {
    // A kite: L (0,4), B (4,5), R (8,4), T (4,3); L-R constrained (mask 32),
    // its sides mask 1. The point (4, 3.5) goes into (L, R, T), and B is then
    // strictly inside the circle of (L, R, p): an unconstrained L-R would flip.
    const Start s = hand_start({{0, 4}, {4, 5}, {8, 4}, {4, 3}}, {{0, 1, 2}, {0, 2, 3}},
                               {{{0, 1}, 1}, {{1, 2}, 1}, {{2, 3}, 1}, {{3, 0}, 1}, {{0, 2}, 32}},
                               [](double, double) { return 0.0; });
    const CheckPoints cp = store(9, {{4.0, 3.5, 40.0f}});
    const auto out = run(cp, s, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 1);
    REQUIRE(constraint_map(out).at({0, 2}) == 32);
    REQUIRE(has_edge(out, 0, 2));
    REQUIRE_FALSE(has_edge(out, 1, 4));  // B to p: the flip a free edge would take
    REQUIRE(out.triangles.size() == 4);
}

// --------------------------------------------------------------------------
// RP6, void
// --------------------------------------------------------------------------

TEST_CASE("RP6: a NaN start vertex is carved around; uncovered 0 at the end", "[refine_points][RP6]") {
    Start s = grid_start(9, 4, [](double, double) { return 0.0; });
    // Vertex 4 is the centre node (4, 4): NoData, as refine outputs it.
    REQUIRE(s.mesh.vertices()[4] == world(4.0, 4.0));
    s.valid[4] = 0;
    std::mt19937 gen{606};
    std::vector<Pt> pts;
    for (int i = 0; i < 120; ++i) {
        const double c = static_cast<double>(1 + gen() % (8u * 1024u - 1)) / 1024.0;
        const double r = static_cast<double>(1 + gen() % (8u * 1024u - 1)) / 1024.0;
        pts.push_back(Pt{c, r, 0.0f});
    }
    const CheckPoints cp = store(9, pts);
    const auto out = run(cp, s, 1e6);  // only carving can insert
    REQUIRE(out.ok());
    REQUIRE(out.uncovered == 0);
    REQUIRE(out.carved > 0);
    REQUIRE(out.inserted == out.carved);
    REQUIRE(out.valid[4] == 0);
    for (std::size_t i = s.mesh.vertices().size(); i < out.vertices.size(); ++i) {
        REQUIRE(out.valid[i] == 1);
        REQUIRE(out.z[i] == 0.0);
    }
    // Independently: no check point is left in a closed triangle with a NoData
    // corner, other than as that triangle's corner.
    for (const auto& p : pts) {
        const Point2 w = world(p.col, p.row);
        for (const auto& t : out.triangles) {
            if (out.valid[t[0]] && out.valid[t[1]] && out.valid[t[2]]) continue;
            if (w == out.vertices[t[0]] || w == out.vertices[t[1]] || w == out.vertices[t[2]]) continue;
            CAPTURE(p.col, p.row);
            REQUIRE_FALSE(in_closed(out, t, w));
        }
    }
}

// --------------------------------------------------------------------------
// RP7, coincident
// --------------------------------------------------------------------------

TEST_CASE("RP7: a point at a start vertex is counted, never inserted", "[refine_points][RP7]") {
    const Start s = grid_start(9, 4, flat10);
    REQUIRE(s.mesh.vertices()[4] == world(4.0, 4.0));
    const CheckPoints cp = store(9, {{4.0, 4.0, 13.25f}, {1.5, 6.25, 10.0f}, {6.75, 2.5, 10.0f}});
    const auto out = run(cp, s, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 0);
    REQUIRE(out.coincident == 1);
    REQUIRE(out.coincident_max_error == 3.25);
    REQUIRE(out.max_error == 0.0);
    REQUIRE(out.z[4] == 10.0);  // the start vertex keeps its own z
}

TEST_CASE("RP7: no coincident point, both counters 0", "[refine_points][RP7]") {
    const Start s = grid_start(9, 4, flat10);
    const CheckPoints cp = store(9, {{1.5, 6.25, 10.0f}});
    const auto out = run(cp, s, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.coincident == 0);
    REQUIRE(out.coincident_max_error == 0.0);
}

// --------------------------------------------------------------------------
// RP8, outside
// --------------------------------------------------------------------------

namespace {
// One triangle (4,4) (4,12) (12,12) in a 17-node grid; z 0.
Start lone_triangle() {
    return hand_start({{4, 4}, {4, 12}, {12, 12}}, {{0, 1, 2}},
                      {{{0, 1}, 1}, {{1, 2}, 1}, {{2, 0}, 1}}, [](double, double) { return 0.0; });
}
const std::vector<Pt> kOutside{{12.0, 4.5, 1000.0f}, {1.5, 1.5, -1000.0f}, {15.5, 2.25, 1000.0f},
                               {8.0, 7.75, 1000.0f}};  // the last a hair above the hypotenuse
}  // namespace

TEST_CASE("RP8: points outside every triangle are never inserted and never raise max_error",
          "[refine_points][RP8]") {
    const auto out = run(store(17, kOutside), lone_triangle(), 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 0);
    REQUIRE(out.max_error == 0.0);
    REQUIRE(out.triangles.size() == 1);
}

TEST_CASE("RP8: the control, one point inside beside them is inserted", "[refine_points][RP8]") {
    auto pts = kOutside;
    pts.push_back({6.5, 9.25, 1000.0f});
    const auto out = run(store(17, pts), lone_triangle(), 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 1);
    REQUIRE(out.vertices.back() == world(6.5, 9.25));
    REQUIRE(out.max_error == 0.0);
}

// --------------------------------------------------------------------------
// Refusals, as refine's
// --------------------------------------------------------------------------

TEST_CASE("refine_points refuses a tolerance that is negative or not finite", "[refine_points][refusal]") {
    const Start s = grid_start(9, 4, flat10);
    const CheckPoints cp = store(9, {{1.5, 1.5, 10.0f}});
    for (const double tol : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity()}) {
        CAPTURE(tol);
        const auto out = run(cp, s, tol);
        REQUIRE_FALSE(out.ok());
        REQUIRE(out.status == RefineStatus::InvalidTolerance);
        REQUIRE_FALSE(out.message.empty());
        REQUIRE(out.vertices.empty());
    }
}

TEST_CASE("refine_points refuses a start vertex outside the store's grid", "[refine_points][refusal]") {
    const Start s = hand_start({{0, 0}, {0, 8}, {9.5, 8}}, {{0, 1, 2}}, {},
                               [](double, double) { return 0.0; });
    const auto out = run(store(9, {{1.5, 6.5, 10.0f}}), s, 1.0);
    REQUIRE(out.status == RefineStatus::OutsideGrid);
    REQUIRE_FALSE(out.message.empty());
}
