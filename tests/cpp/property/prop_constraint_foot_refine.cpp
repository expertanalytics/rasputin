// Increment 20c, PR 20c-1 (docs/increments/20c-soft-quality.md, R3 and "Tests
// @tester writes red first", 20c-1, CF3): refinement's foot rule (20b R1)
// widened to the holding triangle's neighbours. Part of the invariant-critical
// test_constraint_foot suite; its own target. It names nothing new, so on the
// red commit it builds and fails by assertion.
//
// What R3 says, and how each part is observed from refine's output alone:
//   - a worst node whose holding triangle has no constrained edge, but whose
//     neighbour has one within eps, puts its foot on the neighbour's edge:
//     the foot is an output vertex, and out.feet counts it;
//   - the holding triangle stays active (the round's skipped list), so its
//     node still goes in when its error stays above the tolerance: the
//     fixture's holding triangle is untouched by the foot (no flip reaches
//     it; checked below by replaying the split on the start mesh), so only
//     that rule brings the node in, and the tolerance oracle fails without it
//     (the mutant "holding triangle dropped from the active set", CF3);
//   - with that neighbour touched this round, the node waits a round: the
//     deferral fixture holds the section D oracles and is bit-identical over
//     threads. (Which round the foot lands in is not visible from the output,
//     so the wait itself is not asserted.)
//
// The QA rules' section D oracles come from support/constraint_foot_oracles.hpp.
//
// Off switch: a features start (an interior constraint line) with the
// quality start on and feet off gives master's output. The topology digest
// was RECORDED FROM c074f900 (master ed125121's production code; the branch
// changed only docs) with refine_digest::topology_digest, before any
// production change. No commit may update it to agree with new code.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/scan.hpp>

#include "constraint_foot_fixtures.hpp"
#include "constraint_foot_oracles.hpp"
#include "quality_fixtures.hpp"
#include "refine_digest.hpp"
#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <variant>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::TriangleIndices;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;
using quality_fixtures::Start;

namespace cfo = constraint_foot_oracles;
namespace cff = constraint_foot_fixtures;

namespace {

// World nodes at integer (x, y): x_min -5, y_max 32, dx = dy = 1, 31 columns
// and 39 rows, so node (row r, col c) is (c - 5, 32 - r).
RasterGeometry geometry() { return RasterGeometry{-5.0, 32.0, 1.0, 1.0, 31, 39}; }

// The constraint P-Q runs along y = -0.3 from x = 0 to x = 20. Above it, the
// sliver u = (P, Q, V) with V (4, 0.1), then t = (P, V, W) and t2 = (V, Q, W)
// with W (2, 29.7) far up, so t's angle at W is small and its circle bulges
// little past its chord P-V. The node N (2, 0) lies in t, 0.3 above P-Q; t has
// no constrained edge, its neighbour u has P-Q. Only P-Q is a constraint.
Start start() {
    Start s;
    s.mesh = IndexedMesh2{{{0, -0.3}, {20, -0.3}, {4, 0.1}, {2, 29.7}},
                          {{0, 1, 2}, {0, 2, 3}, {2, 1, 3}},
                          std::vector<std::uint8_t>(3, 0)};
    s.edges = {{0, 1}};
    s.masks = {1};
    return s;
}

constexpr std::size_t kNRow = 32, kNCol = 7;  // N = (2, 0)
const cfo::Frac kN{7.0, 32.0};
const cfo::Frac kFoot{7.0, 32.3};  // N's foot on P-Q: (2, -0.3) in world

// Zero, with a 1 m bump at N, and optionally at (6, 0), which lies in u.
Raster<float> dem(bool second_bump = false) {
    const auto g = geometry();
    std::vector<float> z(g.rows() * g.cols(), 0.0f);
    z[kNRow * g.cols() + kNCol] = 1.0f;
    if (second_bump) z[kNRow * g.cols() + 11] = 1.0f;
    return Raster<float>{g, std::move(z)};
}

constexpr double kTol = 0.5;

auto run(const Raster<float>& d, const Start& s, double tol, bool feet, unsigned threads = 1, double min_angle = 0.0) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.min_angle_deg = min_angle;
    o.constraint_feet = feet;
    return refine(d, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks}, o);
}

template <typename Outcome>
void section_d(const Raster<float>& d, const Start& s, const Outcome& out, double tol) {
    REQUIRE(out.ok());
    const auto n = cfo::node_findings(d, out, tol);
    CHECK(n.not_ccw == 0);
    CHECK(n.over == 0);
    CHECK(out.max_error <= tol);
    CHECK(cfo::delaunay_violations(d.geometry(), out) == 0);
    const std::vector<Point2> sv(s.mesh.vertices().begin(), s.mesh.vertices().end());
    const auto l = cfo::line_findings(d.geometry(), sv, s.edges, s.masks, out);
    CHECK(l.unplaced == 0);
    CHECK(l.wrong_mask == 0);
    CHECK(l.broken_chain == 0);
}

}  // namespace

// ------------------------------------------------------------------ the fixture

TEST_CASE("CF3: the fixture is what it claims", "[refinement][constraint_foot][CF3]") {
    // Replayed with refine's own pieces on the start mesh, as refine's first
    // round sees it: the holding triangle of the worst node, its neighbour
    // with the constraint, eps, and that the foot's split leaves the holding
    // triangle's slot unwritten with its error still above the tolerance.
    namespace rd = terrain::refinement::detail;
    const auto d = dem();
    const auto s = start();
    auto built = rd::to_lattice(d.geometry(), s.mesh, s.edges, s.masks);
    REQUIRE(std::holds_alternative<terrain::mesh::LatticeMesh>(built));
    auto& m = std::get<terrain::mesh::LatticeMesh>(built);
    const auto frame = terrain::mesh::lattice_frame(1.0, 1.0, d.geometry().rows(), d.geometry().cols());
    REQUIRE(terrain::mesh::legalise_all<terrain::pred::DefaultKernel>(m, frame, [](std::uint32_t) {}) == 0);

    const auto r = terrain::refinement::scan(d, m, 1);
    REQUIRE(r.node.has_value());
    REQUIRE(r.node->row == kNRow);
    REQUIRE(r.node->col == kNCol);
    REQUIRE(r.max_error > kTol);
    for (unsigned k = 0; k < 3; ++k) REQUIRE_FALSE(m.is_constrained(1, k));
    REQUIRE(m.neighbours(1)[0] == 0u);  // across P-V
    REQUIRE(m.is_constrained(0, 0));    // u's P-Q
    REQUIRE(rd::foot_epsilon(d, *r.node, kTol) > 0.3);  // N is 0.3 from P-Q

    const auto before = m.triangles()[1];
    const auto q = m.split_edge(0, 0, terrain::mesh::MeshVertex{kFoot.col, kFoot.row});
    const std::array<std::uint32_t, 2> seeds{0, static_cast<std::uint32_t>(m.triangle_count() - 1)};
    std::vector<std::uint32_t> written;
    terrain::mesh::legalise_around<terrain::pred::DefaultKernel>(m, q, std::span<const std::uint32_t>{seeds}, frame,
                                                                [&](std::uint32_t w) { written.push_back(w); });
    REQUIRE(m.triangles()[1] == before);
    REQUIRE(std::find(written.begin(), written.end(), 1u) == written.end());
    REQUIRE(terrain::refinement::scan(d, m, 1).max_error > kTol);
}

// ------------------------------------------------------------------ CF3

TEST_CASE("CF3: the foot goes on the neighbour's edge, and the node still goes in from the rescanned holding triangle",
          "[refinement][constraint_foot][CF3]") {
    const auto d = dem();
    const auto s = start();
    const auto out = run(d, s, kTol, true);
    REQUIRE(out.ok());
    REQUIRE(cfo::vertex_at(d.geometry(), out, kFoot).has_value());  // the foot, on P-Q
    REQUIRE(out.feet >= 1);
    REQUIRE(out.feet_refused == 0);
    REQUIRE(cfo::vertex_at(d.geometry(), out, kN, 0.0).has_value());  // N, after its foot (20b R5's fallback)
    section_d(d, s, out, kTol);
}

TEST_CASE("CF3: with feet off the node goes in and nothing lands on P-Q", "[refinement][constraint_foot][CF3]") {
    const auto d = dem();
    const auto s = start();
    const auto out = run(d, s, kTol, false);
    REQUIRE(out.feet == 0);
    REQUIRE_FALSE(cfo::vertex_at(d.geometry(), out, kFoot).has_value());
    REQUIRE(cfo::vertex_at(d.geometry(), out, kN, 0.0).has_value());
    section_d(d, s, out, kTol);
}

TEST_CASE("CF3: with the neighbour touched earlier in the round the oracles hold, for any thread count",
          "[refinement][constraint_foot][CF3]") {
    // A second bump at (6, 0) inside u: u (slot 0) is split first in round 1,
    // so when t (slot 1) asks for its foot, the owner was touched this round
    // (R3: the node waits).
    const auto d = dem(true);
    const auto s = start();
    const auto ref = run(d, s, kTol, true);
    section_d(d, s, ref, kTol);
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    REQUIRE(refine_digest::digest(run(d, s, kTol, true, threads)) == refine_digest::digest(ref));
}

// ------------------------------------------------------------------ T3 and the off switch

TEST_CASE("CF3 T3: a features start on rough ground, feet on, holds the section D oracles",
          "[refinement][constraint_foot][CF3][T3]") {
    const double tol = GENERATE(0.0, 1.0, 5.0);
    const double min_angle = GENERATE(0.0, 25.0);
    CAPTURE(tol, min_angle);
    const Raster<float> d{cff::geometry(), refinement_fixtures::rough_dem(cff::kN, cff::kN, 7)};
    const auto s = cff::line_start(4.3);
    section_d(d, s, run(d, s, tol, true, 1, min_angle), tol);
}

TEST_CASE("CF3 off switch: a features start with the quality start on and feet off is master's output",
          "[refinement][constraint_foot][off]") {
    const Raster<float> d{cff::geometry(), refinement_fixtures::rough_dem(cff::kN, cff::kN, 7)};
    const auto s = cff::line_start(4.3);
    const auto off = run(d, s, 1.0, false, 1, 25.0);
    REQUIRE(off.ok());
    REQUIRE(refine_digest::topology_digest(off) == 0x1253792f2761a97dull);  // RECORDED, see the header
    REQUIRE(off.feet == 0);
    REQUIRE(off.quality_inserted > 0);  // the quality start ran
}
