// Increment 20c, PR 20c-1 (docs/increments/20c-soft-quality.md, R2, R5 and
// "Tests @tester writes red first", 20c-1, CF2): the quality start puts a node
// that lies near a constraint onto the constraint, at its foot. Part of the
// invariant-critical test_constraint_foot suite; its own target so a compile
// failure here leaves CF1, CF3 and CF4 built and run.
//
// Interface, as R2 and R5 name it:
//
//   QualityOptions::constraint_feet   bool, false by default; set by member assignment
//   QualityOutcome::feet              feet inserted instead of the node
//   QualityOutcome::skipped_near_line NearEnd, NotCounterClockwise, or a foot the
//                                     validity callable refuses
//   improve<K>(m, f, o, valid)        `valid` is called with a MeshVertex (R2.3)
//   RefineOutcome::quality_feet       the quality start's feet, through refine
//   RefineOutcome::quality_skipped    still one total, skipped_near_line in it (R5)
//
// PINNED HERE, where the design leaves it open (listed for @architect):
//   - QualityOutcome::inserted counts DEM nodes only, as its comment says, and
//     a foot is counted in `feet` alone; so the output holds n0 + inserted +
//     feet vertices, and through refine n0 + quality_inserted + quality_feet +
//     inserted (refinement's own feet stay inside `inserted`, as 20b has them);
//   - the exact counts on the two hand fixtures (one foot, no node), recorded
//     by running today's improve on the mesh a foot leaves (see the handback).
//
// Oracles are this file's own: topology and Delaunay on the LatticeMesh by the
// exact kernel in the frame, and the constraint lines (Q3's check, to 1e-9
// cells because a foot is a rounded projection). Through refine they are
// support/constraint_foot_oracles.hpp's, the QA rules' section D pair.
//
// Mutants these cases are meant to kill: the foot branch skipped (node goes in;
// "the output has the foot and not the node"), the neighbour search missing in
// the quality path (the neighbour fixture), a refused foot inserting the node
// instead of skipping, skipped_near_line left out of quality_skipped, and the
// switch ignored (feet with constraint_feet false).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/mesh/quality.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/refinement/refine.hpp>

#include "constraint_foot_oracles.hpp"
#include "feet_fixtures.hpp"
#include "quality_fixtures.hpp"
#include "refine_digest.hpp"
#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <set>
#include <span>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::improve;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::legalise_all;
using terrain::mesh::MeshVertex;
using terrain::mesh::orient_sign;
using terrain::mesh::QualityOptions;
using terrain::mesh::QualityOutcome;
using terrain::pred::DefaultKernel;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;

namespace cfo = constraint_foot_oracles;

namespace {

constexpr double kOnLine = 1e-9;  // cells, for a 21 x 21 grid

using Pair = std::pair<std::uint32_t, std::uint32_t>;

struct Fixture {
    std::vector<MeshVertex> vertices;
    std::vector<TriangleIndices> triangles;
    std::map<Pair, std::uint32_t> constraints;  // undirected pair -> mask
    std::size_t rows = 0;
    std::size_t cols = 0;
};

LatticeMesh build(const Fixture& f) {
    std::vector<std::uint8_t> bits(f.triangles.size(), 0);
    std::vector<std::array<std::uint32_t, 3>> masks(f.triangles.size(), {0, 0, 0});
    for (std::size_t t = 0; t < f.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto it = f.constraints.find(std::minmax(f.triangles[t][k], f.triangles[t][(k + 1) % 3]));
                it != f.constraints.end()) {
                bits[t] |= static_cast<std::uint8_t>(1u << k);
                masks[t][k] = it->second;
            }
    auto m = LatticeMesh::build(f.vertices, f.triangles, std::move(bits), std::move(masks));
    REQUIRE(m.has_value());
    return std::move(*m);
}

// The own-edge fixture: Q3's right triangle with its hypotenuse moved 0.01
// rows up. A (0, 0.99), C (2, 6.99), B (20, 0.99): right-angled at C, so its
// circumcentre is the hypotenuse's midpoint (10, 0.99) and the snapped node is
// N = (10, 1), inside the triangle, 0.01 cells from the constrained A-B. The
// angle at B is 18.4 degrees: bad at 25. Masks A-B 1, A-C 2, C-B 4.
Fixture own_edge() {
    Fixture f;
    f.vertices = {{0, 0.99}, {2, 6.99}, {20, 0.99}};
    f.triangles = {{0, 1, 2}};
    f.constraints = {{{0u, 2u}, 1u}, {{0u, 1u}, 2u}, {{1u, 2u}, 4u}};
    f.rows = 8;
    f.cols = 21;
    return f;
}

// The neighbour fixture, in (col, row) with the world's "up" as smaller row.
// a (8, 10) and b (12, 10) on row 10, c on the circle of diameter a-b at 40
// degrees from b's end: t = (a, b, c) is right-angled at c, 20 degrees at a,
// its circumcentre exactly N = (10, 10), on the FREE edge a-b. Below a-b, u =
// (b, a, d) with d 6 cells from a at 10 degrees: u's edge a-d is constrained
// and passes 0.347 cells from N, N's foot 1.97 cells from a. legalise_all
// leaves both (no flip), and u, the worse, is offered first and blocked. So
// the node t's slot names is reached only through t's neighbour.
Fixture neighbour() {
    const double deg = std::numbers::pi / 180.0;
    Fixture f;
    f.vertices = {{8, 10}, {12, 10}, {10 + 2 * std::cos(40 * deg), 10 - 2 * std::sin(40 * deg)},
                  {8 + 6 * std::cos(10 * deg), 10 + 6 * std::sin(10 * deg)}};
    f.triangles = {{0, 1, 2}, {1, 0, 3}};
    f.constraints = {{{1u, 2u}, 1u}, {{0u, 2u}, 2u}, {{0u, 3u}, 4u}, {{1u, 3u}, 8u}};
    f.rows = f.cols = 21;
    return f;
}

QualityOutcome run(LatticeMesh& m, const Fixture& f, std::optional<bool> feet) {
    const LatticeFrame frame{1.0, 1.0};
    legalise_all<DefaultKernel>(m, frame, [](std::uint32_t) {});
    QualityOptions o{25.0, f.rows, f.cols};
    if (feet) o.constraint_feet = *feet;
    return improve<DefaultKernel>(m, frame, o);
}

// The orthogonal projection of p on a-b, frame (1, 1).
MeshVertex foot_of(MeshVertex a, MeshVertex b, MeshVertex p) {
    const double ux = b.col - a.col, uy = b.row - a.row;
    const double s = ((p.col - a.col) * ux + (p.row - a.row) * uy) / (ux * ux + uy * uy);
    return MeshVertex{a.col + s * ux, a.row + s * uy};
}

bool has_vertex(const LatticeMesh& m, MeshVertex p, double eps = 1e-12) {
    return std::any_of(m.vertices().begin(), m.vertices().end(), [&](MeshVertex v) {
        return std::abs(v.col - p.col) <= eps && std::abs(v.row - p.row) <= eps;
    });
}

bool has_triangle(const LatticeMesh& m, const TriangleIndices& t) {
    for (const auto& x : m.triangles())
        for (unsigned r = 0; r < 3; ++r)
            if (x[r] == t[0] && x[(r + 1) % 3] == t[1] && x[(r + 2) % 3] == t[2]) return true;
    return false;
}

// Topology from the triangles alone: positive orientation, each directed edge
// once, neighbour links and constraint bits and masks agreeing on both sides.
void check_topology(const LatticeMesh& m) {
    std::map<Pair, std::uint32_t> directed;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        CAPTURE(t);
        REQUIRE(orient_sign(m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)) > 0);
        for (unsigned k = 0; k < 3; ++k)
            REQUIRE(directed.emplace(Pair{m.triangles()[t][k], m.triangles()[t][(k + 1) % 3]}, t).second);
    }
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto& tri = m.triangles()[t];
            const auto it = directed.find({tri[(k + 1) % 3], tri[k]});
            CAPTURE(t, k);
            REQUIRE(m.neighbours(t)[k] == (it == directed.end() ? kNoNeighbour : it->second));
            if (it == directed.end()) continue;
            unsigned j = 0;
            while (m.triangles()[it->second][j] != tri[(k + 1) % 3]) ++j;
            REQUIRE(m.is_constrained(t, k) == m.is_constrained(it->second, j));
            REQUIRE(m.mask(t, k) == m.mask(it->second, j));
        }
}

std::size_t delaunay_violations(const LatticeMesh& m, const LatticeFrame& f) {
    std::size_t bad = 0;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto u = m.neighbours(t)[k];
            if (u == kNoNeighbour || m.is_constrained(t, k)) continue;
            unsigned j = 0;
            while (m.neighbours(u)[j] != t) ++j;
            if (DefaultKernel::incircle(f.at(m.corner(t, 0)), f.at(m.corner(t, 1)), f.at(m.corner(t, 2)),
                                        f.at(m.corner(u, (j + 2) % 3)))
                == terrain::pred::Incircle::Inside)
                ++bad;
        }
    return bad;
}

// Q3's check, to kOnLine: every constraint edge on exactly one input segment
// with its mask, and each segment's pieces chaining end to end.
void check_lines(const LatticeMesh& m, const Fixture& f) {
    const auto v = m.vertices();
    const auto [edges, masks] = m.constraint_edges();
    std::map<Pair, std::vector<std::pair<double, double>>> pieces;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        CAPTURE(i);
        std::size_t owners = 0;
        for (const auto& [seg, mask] : f.constraints) {
            const MeshVertex a = f.vertices[seg.first], b = f.vertices[seg.second];
            const cfo::Frac fa{a.col, a.row}, fb{b.col, b.row};
            const auto [tp, dp] = cfo::param_dist(fa, fb, cfo::Frac{v[edges[i][0]].col, v[edges[i][0]].row});
            const auto [tq, dq] = cfo::param_dist(fa, fb, cfo::Frac{v[edges[i][1]].col, v[edges[i][1]].row});
            if (dp > kOnLine || dq > kOnLine || std::min(tp, tq) < -kOnLine || std::max(tp, tq) > 1 + kOnLine) continue;
            ++owners;
            REQUIRE(masks[i] == mask);  // bits and masks on both halves
            pieces[seg].push_back(std::minmax(tp, tq));
        }
        REQUIRE(owners == 1);
    }
    REQUIRE(pieces.size() == f.constraints.size());
    for (auto& [seg, ps] : pieces) {
        std::sort(ps.begin(), ps.end());
        REQUIRE(std::abs(ps.front().first) <= kOnLine);
        REQUIRE(std::abs(ps.back().second - 1.0) <= kOnLine);
        for (std::size_t j = 1; j < ps.size(); ++j) REQUIRE(std::abs(ps[j].first - ps[j - 1].second) <= kOnLine);
    }
}

void check_mesh(const LatticeMesh& m, const Fixture& f) {
    check_topology(m);
    REQUIRE(delaunay_violations(m, LatticeFrame{1.0, 1.0}) == 0);
    check_lines(m, f);
}

}  // namespace

// ------------------------------------------------------------------ the switch

TEST_CASE("CF2: QualityOptions leaves feet off by default", "[mesh][quality][constraint_foot][CF2]") {
    REQUIRE_FALSE(QualityOptions{}.constraint_feet);
}

// ------------------------------------------------------------------ own edge

TEST_CASE("CF2: a node 0.01 cells from its triangle's constrained edge goes in as its foot",
          "[mesh][quality][constraint_foot][CF2]") {
    const Fixture fx = own_edge();
    const MeshVertex node{10, 1};
    const MeshVertex foot = foot_of(fx.vertices[0], fx.vertices[2], node);
    REQUIRE(std::abs(foot.row - 0.99) <= 1e-15);  // the fixture is what it claims

    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, fx, true);
    REQUIRE(q.feet == 1);
    REQUIRE(q.inserted == 0);
    REQUIRE(q.skipped_near_line == 0);
    REQUIRE(m.vertices().size() == fx.vertices.size() + q.inserted + q.feet);
    REQUIRE(has_vertex(m, foot));
    REQUIRE_FALSE(has_vertex(m, node));
    REQUIRE_FALSE(has_triangle(m, {0, 1, 2}));  // the bad triangle is gone
    check_mesh(m, fx);
}

TEST_CASE("CF2: with feet off the same fixture gets the node, as today", "[mesh][quality][constraint_foot][CF2]") {
    const Fixture fx = own_edge();
    const auto off = GENERATE(std::optional<bool>{}, std::optional<bool>{false});
    CAPTURE(off.has_value());
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, fx, off);
    REQUIRE(q.inserted == 1);
    REQUIRE(q.feet == 0);
    REQUIRE(q.skipped_near_line == 0);
    REQUIRE(has_vertex(m, MeshVertex{10, 1}, 0.0));
}

TEST_CASE("CF2: a foot the validity callable refuses is a skip, and neither foot nor node goes in",
          "[mesh][quality][constraint_foot][CF2]") {
    // R2.1 and R2.3: the callable takes a MeshVertex; this one accepts nodes
    // only, so the node N is valid and its foot is not.
    const Fixture fx = own_edge();
    LatticeMesh m = build(fx);
    const LatticeFrame frame{1.0, 1.0};
    legalise_all<DefaultKernel>(m, frame, [](std::uint32_t) {});
    QualityOptions o{25.0, fx.rows, fx.cols};
    o.constraint_feet = true;
    const QualityOutcome q =
        improve<DefaultKernel>(m, frame, o, [](const MeshVertex& v) { return v.is_node(); });
    REQUIRE(q.skipped_near_line == 1);
    REQUIRE(q.feet == 0);
    REQUIRE(q.inserted == 0);
    REQUIRE(m.vertices().size() == fx.vertices.size());
}

// ------------------------------------------------------------------ the neighbour

TEST_CASE("CF2: a node on a free edge, near the neighbour's constrained edge, goes in as its foot there",
          "[mesh][quality][constraint_foot][CF2]") {
    const Fixture fx = neighbour();
    const MeshVertex node{10, 10};
    const MeshVertex foot = foot_of(fx.vertices[0], fx.vertices[3], node);
    {  // the fixture is what it claims: N on a-b, 0.347 from a-d, 1.97 along it
        REQUIRE(orient_sign(fx.vertices[0], fx.vertices[1], node) == 0);
        const auto [s, d] = cfo::param_dist({8, 10}, {fx.vertices[3].col, fx.vertices[3].row}, {10, 10});
        REQUIRE(d < 0.35);
        REQUIRE(d > 0.34);
        REQUIRE(s * 6.0 > 1.9);
    }
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, fx, true);
    REQUIRE(q.feet == 1);
    REQUIRE(q.inserted == 0);
    REQUIRE(m.vertices().size() == fx.vertices.size() + q.inserted + q.feet);
    REQUIRE(has_vertex(m, foot));
    REQUIRE_FALSE(has_vertex(m, node));
    REQUIRE_FALSE(has_triangle(m, {0, 1, 2}));
    check_mesh(m, fx);
}

TEST_CASE("CF2: with feet off the neighbour fixture gets the node on the free edge, as today",
          "[mesh][quality][constraint_foot][CF2]") {
    const Fixture fx = neighbour();
    LatticeMesh m = build(fx);
    const QualityOutcome q = run(m, fx, false);
    REQUIRE(q.inserted == 1);
    REQUIRE(q.feet == 0);
    REQUIRE(has_vertex(m, MeshVertex{10, 10}, 0.0));
}

// ------------------------------------------------------------------ through refine

namespace {

// The own-edge fixture in world: 21 x 8 nodes, x_min 0, y_max 7, dx = dy = 1,
// so world (col, 7 - row). A flat DEM at 3 m: refinement adds nothing.
RasterGeometry own_geometry() { return RasterGeometry{0.0, 7.0, 1.0, 1.0, 21, 8}; }

quality_fixtures::Start own_start() {
    const auto w = [](double c, double r) { return Point2{c, 7.0 - r}; };
    quality_fixtures::Start s;
    s.mesh = terrain::IndexedMesh2{{w(0, 0.99), w(2, 6.99), w(20, 0.99)}, {{0, 1, 2}}, std::vector<std::uint8_t>(1, 0)};
    s.edges = {{0, 2}, {0, 1}, {1, 2}};
    s.masks = {1, 2, 4};
    return s;
}

Raster<float> flat(const RasterGeometry& g, std::optional<std::pair<std::size_t, std::size_t>> hole = std::nullopt) {
    std::vector<float> z(g.rows() * g.cols(), 3.0f);
    if (hole) z[hole->first * g.cols() + hole->second] = std::numeric_limits<float>::quiet_NaN();
    return Raster<float>{g, std::move(z)};
}

auto run_refine(const Raster<float>& dem, const quality_fixtures::Start& s, double tol, bool feet,
                double min_angle = 25.0, unsigned threads = 1) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.min_angle_deg = min_angle;
    o.constraint_feet = feet;
    return refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks}, o);
}

template <typename Outcome>
void section_d(const Raster<float>& dem, const quality_fixtures::Start& s, const Outcome& out, double tol) {
    REQUIRE(out.ok());
    const auto n = cfo::node_findings(dem, out, tol);
    CHECK(n.not_ccw == 0);
    CHECK(n.over == 0);
    CHECK(out.max_error <= tol);
    CHECK(cfo::delaunay_violations(dem.geometry(), out) == 0);
    const std::vector<Point2> start(s.mesh.vertices().begin(), s.mesh.vertices().end());
    const auto l = cfo::line_findings(dem.geometry(), start, s.edges, s.masks, out);
    CHECK(l.unplaced == 0);
    CHECK(l.wrong_mask == 0);
    CHECK(l.broken_chain == 0);
    // Every vertex is a start vertex, a node or a foot, counted once each (pinned).
    CHECK(out.vertices.size() == s.mesh.vertices().size() + out.quality_inserted + out.quality_feet + out.inserted);
}

}  // namespace

TEST_CASE("CF2 through refine: the quality start's foot is counted in quality_feet, with a bilinear z",
          "[refinement][quality][constraint_foot][CF2]") {
    const auto g = own_geometry();
    const auto dem = flat(g);
    const auto s = own_start();
    const auto out = run_refine(dem, s, 1.0, true);
    section_d(dem, s, out, 1.0);
    REQUIRE(out.quality_feet == 1);
    REQUIRE(out.quality_inserted == 0);
    const auto f = cfo::vertex_at(g, out, cfo::Frac{10.0, 0.99}, 1e-9);
    REQUIRE(f.has_value());
    REQUIRE(out.valid[*f] == 1);
    // 1e-9 m on a 3 m flat surface: the output's lattice bilinear against the
    // world-point bilinear, which may differ in the last bits.
    REQUIRE(std::abs(out.z[*f] - terrain::raster::bilinear(dem, out.vertices[*f]).value()) <= 1e-9);
    REQUIRE_FALSE(cfo::vertex_at(g, out, cfo::Frac{10.0, 1.0}, 0.0).has_value());

    const auto off = run_refine(dem, s, 1.0, false);  // R2.6: refine's switch is the pass's
    REQUIRE(off.quality_feet == 0);
    REQUIRE(off.quality_inserted == 1);
    REQUIRE(cfo::vertex_at(g, off, cfo::Frac{10.0, 1.0}, 0.0).has_value());
}

TEST_CASE("CF2 through refine: a foot on NoData is a skip in quality_skipped, and the node does not go in",
          "[refinement][quality][constraint_foot][CF2]") {
    // Node (row 0, col 10) is NoData: a corner of the foot's cell (rows 0-1),
    // not of the node's own (the node N = (10, 1) is valid), and outside the
    // triangle (rows >= 0.99), so refinement has no void to carve.
    const auto g = own_geometry();
    const auto dem = flat(g, std::pair<std::size_t, std::size_t>{0, 10});
    const auto s = own_start();
    const auto out = run_refine(dem, s, 1.0, true);
    section_d(dem, s, out, 1.0);
    REQUIRE(out.quality_feet == 0);
    REQUIRE(out.quality_inserted == 0);
    REQUIRE(out.quality_skipped >= 1);  // R5: one total, the near-line skip in it
    REQUIRE_FALSE(cfo::vertex_at(g, out, cfo::Frac{10.0, 1.0}, 0.0).has_value());
    REQUIRE_FALSE(cfo::vertex_at(g, out, cfo::Frac{10.0, 0.99}, 1e-9).has_value());
}

// ------------------------------------------------------------------ 14b T3 and T6, rule on

TEST_CASE("CF2 T3: with the quality start and feet on, the section D oracles hold on the needle and the ring",
          "[refinement][quality][constraint_foot][CF2][T3]") {
    const double tol = GENERATE(0.0, 0.5, 3.0);
    const bool needle = GENERATE(true, false);
    CAPTURE(tol, needle);
    if (needle) {
        const auto dem = feet_fixtures::needle_dem(6.0, 1.0);
        const auto s = feet_fixtures::needle_start(dem.geometry(), 1.0 / 7.77);
        section_d(dem, s, run_refine(dem, s, tol, true), tol);
    } else {
        const std::size_t n = 33;
        const Raster<float> dem{refinement_fixtures::geometry(n, n), refinement_fixtures::rough_dem(n, n, 11)};
        const auto s = quality_fixtures::fan(dem.geometry(), quality_fixtures::q7_ring());
        section_d(dem, s, run_refine(dem, s, tol, true), tol);
    }
}

TEST_CASE("CF2 T6: with the quality start and feet on, the output is bit-identical for 1, 2, 7 and all threads",
          "[refinement][quality][constraint_foot][CF2][T6]") {
    const auto dem = feet_fixtures::needle_dem(6.0, 1.0);
    const auto s = feet_fixtures::needle_start(dem.geometry(), 1.0 / 7.77);
    const auto ref = run_refine(dem, s, 0.1, true, 25.0, 1);
    REQUIRE(ref.ok());
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    const auto other = run_refine(dem, s, 0.1, true, 25.0, threads);
    REQUIRE(refine_digest::digest(other) == refine_digest::digest(ref));
    REQUIRE(other.quality_feet == ref.quality_feet);
    REQUIRE(other.quality_inserted == ref.quality_inserted);
    REQUIRE(other.feet == ref.feet);
}
