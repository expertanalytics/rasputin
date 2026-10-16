// Increment 20c, PR 20c-1 (docs/increments/20c-soft-quality.md, R4, R5, R6 and
// "Tests @tester writes red first", 20c-1, CF4 and the off switch): the final
// check puts a stored point (or, on the projected path, a DEM node) that lies
// near a constraint onto it, at its foot, and inserts the point itself later
// if its error stays above the tolerance. Part of the invariant-critical
// test_constraint_foot suite; its own target.
//
// Interface, as R4.5, R4.6 and R5 name it:
//
//   PointRefineOptions::constraint_feet   bool, false by default; set by member assignment
//   PointRefineOutcome::feet              feet inserted (RefineOutcome's field)
//   PointRefineOutcome::feet_fallback     footed points later inserted as themselves
//   PointRefineOutcome::feet_refused      a foot found and not inserted (R5)
//
// PINNED HERE, where the design leaves it open (listed for @architect):
//   - feet are counted in `inserted`, as refine's are (20b), so the output has
//     n0 + inserted vertices on every path; feet_fallback <= feet;
//   - R4.3's z on the reprojected path is asserted to 1e-9 m against the
//     strip's own values, and to 1e-9 m against vertex_z on the projected one;
//   - the stored points here come from an exact-position store
//     (support/constraint_foot_fixtures.hpp, ExactStore), which refine_points
//     accepts as it accepts any Store with geometry(), frozen() and for_each_in.
//
// Oracles: on the reprojected path the tolerance oracle is J2 at every stored
// point (support/j2_oracle.hpp), as ES9 has it; on the projected path the DEM
// node oracle; Delaunay and the constraint lines on both
// (support/constraint_foot_oracles.hpp). The Delaunay oracle excuses an apex
// inside by no more than 1e-7 * min(dx, dy) (20c-1's green-step ruling 2) and
// every case prints the depth of each quad it excuses. Feet are recognised from the output:
// an inserted vertex that is neither a stored point, a strip point nor a node.
//
// Mutants these cases are meant to kill (design, "Mutants to kill"):
//   the foot's z from the point          CF4 z: 62, not the point's 80
//   footed-once removed                  the two-lines case (no second foot on A1-B1);
//                                        the wedge's bound does not reach it
//   the fallback removed                 CF4 J2 at tolerance 0 and 3
//   and beyond the list: L12 not first (the point within r(g) goes in at its
//   own position), a point on an edge footed, the switch ignored.
//
// Off switch: refine_points (with a strip) and refine_strip with
// constraint_feet = false give master's output. The topology digests were
// RECORDED FROM c074f900 (master ed125121's production code; the branch
// changed only docs) with refine_digest::topology_digest and the default
// options, before any production change. No commit may update them to agree
// with new code.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/mesh/constraint_foot.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>

#include "constraint_foot_fixtures.hpp"
#include "constraint_foot_oracles.hpp"
#include "feet_fixtures.hpp"
#include "j2_oracle.hpp"
#include "quality_fixtures.hpp"
#include "refine_digest.hpp"
#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <variant>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::MeshVertex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::constraint_check_points;
using terrain::refinement::ConstraintCheckPoints;
using terrain::refinement::PointRefineOptions;
using terrain::refinement::PointRefineOutcome;
using terrain::refinement::refine;
using terrain::refinement::refine_points;
using terrain::refinement::refine_strip;
using terrain::refinement::RefineOptions;
using quality_fixtures::Start;

namespace cfo = constraint_foot_oracles;
namespace cff = constraint_foot_fixtures;

namespace {

// What a final check starts from: a mesh in world coordinates with z, valid,
// and its constraint edges.
struct Begin {
    IndexedMesh2 mesh;
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    std::vector<std::array<std::uint32_t, 2>> edges;
    std::vector<std::uint32_t> masks;
};

Begin with_z(const Start& s, std::vector<double> z) {
    return Begin{s.mesh, std::move(z), std::vector<std::uint8_t>(s.mesh.vertices().size(), 1), s.edges, s.masks};
}

// Phase 1, refine with feet off, as the CLI hands its output on.
Begin phase1(const Raster<float>& dem, const Start& s, double tol) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = 1;
    auto r = refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                    std::span<const std::uint32_t>{s.masks}, o);
    REQUIRE(r.ok());
    const std::size_t nt = r.triangles.size();
    return Begin{IndexedMesh2{std::move(r.vertices), std::move(r.triangles), std::vector<std::uint8_t>(nt, 0)},
                 std::move(r.z), std::move(r.valid), std::move(r.edges), std::move(r.masks)};
}

ConstraintCheckPoints strip_of(const Raster<float>& dem, const Begin& b) {
    return constraint_check_points(dem, b.mesh.vertices(), std::span<const std::array<std::uint32_t, 2>>{b.edges});
}

PointRefineOptions options(double tol, std::optional<bool> feet, unsigned threads = 1) {
    PointRefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    if (feet) o.constraint_feet = *feet;
    return o;
}

template <class Store>
PointRefineOutcome run_points(const Store& store, const ConstraintCheckPoints* strip, const Begin& b,
                              const PointRefineOptions& o) {
    return refine_points(store, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                         std::span<const std::array<std::uint32_t, 2>>{b.edges},
                         std::span<const std::uint32_t>{b.masks}, o, strip);
}

PointRefineOutcome run_strip(const Raster<float>& dem, const ConstraintCheckPoints& strip, const Begin& b,
                             const PointRefineOptions& o) {
    return refine_strip(dem, strip, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                        std::span<const std::array<std::uint32_t, 2>>{b.edges},
                        std::span<const std::uint32_t>{b.masks}, o);
}

// J2 at every stored point that is not a start vertex, in (col, -row).
std::size_t j2_over(const RasterGeometry& g, const cff::ExactStore& store, const Begin& b,
                    const PointRefineOutcome& out, double tol) {
    std::vector<Point2> v, p;
    std::vector<float> pz;
    std::vector<std::uint8_t> skip;
    for (const auto& w : out.vertices) {
        const auto f = cfo::frac(g, w);
        v.push_back(Point2{f.col, -f.row});
    }
    std::set<std::pair<double, double>> start;
    for (const auto& w : b.mesh.vertices()) {
        const auto f = cfo::frac(g, w);
        start.insert({f.col, f.row});
    }
    for (const auto& [q, z] : store.points()) {
        p.push_back(Point2{q.col, -q.row});
        pz.push_back(z);
        skip.push_back(start.count({q.col, q.row}) != 0 ? 1 : 0);
    }
    const auto f = j2_oracle::violations(v, out.z, out.valid, out.triangles, p, pz, skip, g.cols(), g.rows(), tol);
    CHECK(f.not_ccw == 0);
    return f.over;
}

void shape(const RasterGeometry& g, const Begin& b, const PointRefineOutcome& out) {
    REQUIRE(out.ok());
    CHECK(out.vertices.size() == b.mesh.vertices().size() + out.inserted);
    CHECK(out.feet_fallback <= out.feet);
    std::vector<double> excused;
    CHECK(cfo::delaunay_violations(g, out, &excused) == 0);
    for (const double depth : excused) WARN("Delaunay oracle excused a quad, apex inside by " << depth << " m");
    const std::vector<Point2> sv(b.mesh.vertices().begin(), b.mesh.vertices().end());
    const auto l = cfo::line_findings(g, sv, b.edges, b.masks, out);
    CHECK(l.unplaced == 0);
    CHECK(l.wrong_mask == 0);
    CHECK(l.broken_chain == 0);
}

// The strip's DEM for line_start(4.25): 60 m, but 62 m on columns 7 and 8, so
// along A-B the strip reads 62 at cols 7, 7.5 and 8 and 61 at 6.5 and 8.5,
// within 3 m of the start's 60 m everywhere.
Raster<float> step_dem() {
    const auto g = cff::geometry();
    std::vector<float> z(g.rows() * g.cols(), 60.0f);
    for (std::size_t r = 0; r < g.rows(); ++r) z[r * g.cols() + 7] = z[r * g.cols() + 8] = 62.0f;
    return Raster<float>{g, std::move(z)};
}

}  // namespace

// ------------------------------------------------------------------ the switch

TEST_CASE("CF4: PointRefineOptions leaves feet off by default", "[refine_points][constraint_foot][CF4]") {
    REQUIRE_FALSE(PointRefineOptions{}.constraint_feet);
}

// ------------------------------------------------------------------ z, R4.3

TEST_CASE("CF4: a stored point 0.02 cells from a constrained edge goes in as its foot, z from the strip, then itself",
          "[refine_points][constraint_foot][CF4]") {
    // A-B along row 4.25; the point (7.25, 4.27) is 0.02 cells below it with
    // z 80 against the start's 60: error 20, over tolerance 3. Its foot
    // (7.25, 4.25) lies between the strip points at cols 7 and 7.5, both 62:
    // R4.3 gives the foot 62. The point's own z would give 80 and the line's
    // ends 60, so the three are told apart. After the foot the point is still
    // 18 off: it goes in as itself (the fallback).
    const auto g = cff::geometry();
    const auto dem = step_dem();
    const Begin b = with_z(cff::line_start(4.25), std::vector<double>(6, 60.0));
    const auto strip = strip_of(dem, b);
    cff::ExactStore store{g};
    store.add(MeshVertex{7.25, 4.27}, 80.0f);
    const auto out = run_points(store, &strip, b, options(3.0, true));
    shape(g, b, out);
    REQUIRE(out.feet == 1);
    REQUIRE(out.feet_fallback == 1);
    REQUIRE(out.feet_refused == 0);
    const auto foot = cfo::vertex_at(g, out, cfo::Frac{7.25, 4.25}, 1e-12);
    REQUIRE(foot.has_value());
    CAPTURE(out.z[*foot]);
    REQUIRE(std::abs(out.z[*foot] - 62.0) <= 1e-9);  // not 80 (the point's), not 60 (the ends')
    const auto self = cfo::vertex_at(g, out, cfo::Frac{7.25, 4.27}, 1e-12);
    REQUIRE(self.has_value());
    REQUIRE(out.z[*self] == 80.0);
    REQUIRE(j2_over(g, store, b, out, 3.0) == 0);
}

TEST_CASE("CF4: with no strip the foot's z is linear between the edge's two ends",
          "[refine_points][constraint_foot][CF4]") {
    // Ends A 60 m and B 90 m: at 6.75 of 15 cells along, 73.5 m.
    const auto g = cff::geometry();
    const Begin b = with_z(cff::line_start(4.25), {60.0, 60.0, 60.0, 60.0, 60.0, 90.0});
    cff::ExactStore store{g};
    store.add(MeshVertex{7.25, 4.27}, 80.0f);
    const auto out = run_points(store, nullptr, b, options(1.0, true));
    shape(g, b, out);
    REQUIRE(out.feet == 1);
    const auto foot = cfo::vertex_at(g, out, cfo::Frac{7.25, 4.25}, 1e-12);
    REQUIRE(foot.has_value());
    CAPTURE(out.z[*foot]);
    REQUIRE(std::abs(out.z[*foot] - 73.5) <= 1e-9);
    REQUIRE(j2_over(g, store, b, out, 1.0) == 0);
}

TEST_CASE("CF4: on the projected path a footed DEM node's foot takes vertex_z",
          "[refine_strip][constraint_foot][CF4]") {
    // 20b's needle: nodes on column 8 sit 0.0035 to 0.18 cells (3.5 cm to
    // 1.8 m) from the tilted side, inside delta_p = 2.5 m. Phase 1 with feet
    // off leaves them; the strip run at tolerance 0 rescans every slot it
    // writes along the side.
    const auto dem = feet_fixtures::needle_dem(6.0, 1.0);
    const auto& g = dem.geometry();
    const Begin b = phase1(dem, feet_fixtures::needle_start(g), 0.5);
    const auto strip = strip_of(dem, b);
    const auto out = run_strip(dem, strip, b, options(0.0, true));
    shape(g, b, out);
    REQUIRE(out.feet > 0);
    // The node oracle at phase 1's 0.5 m: the strip run holds 0 m only in the
    // slots it writes (15f L3), phase 1's 0.5 m everywhere else.
    const auto n = cfo::node_findings(dem, out, 0.5);
    CHECK(n.over == 0);
    CHECK(n.not_ccw == 0);
    // Every foot: off-node, not a start vertex, not a strip point; z = vertex_z.
    std::set<std::pair<double, double>> strip_at;
    for (std::size_t k = 0; k < strip.edge_count(); ++k)
        for (const auto& p : strip.on_edge(k)) strip_at.insert({p.at.col, p.at.row});
    std::size_t feet = 0;
    for (std::size_t i = b.mesh.vertices().size(); i < out.vertices.size(); ++i) {
        const auto f = cfo::frac(g, out.vertices[i]);
        if (f.col == std::round(f.col) && f.row == std::round(f.row)) continue;  // a node
        bool on_strip = false;
        for (const auto& [c, r] : strip_at) on_strip |= std::abs(c - f.col) <= 1e-9 && std::abs(r - f.row) <= 1e-9;
        if (on_strip) continue;
        ++feet;
        const auto z = terrain::raster::bilinear(dem, out.vertices[i]);
        REQUIRE(z.has_value());
        CAPTURE(i, out.z[i], *z);
        REQUIRE(std::abs(out.z[i] - *z) <= 1e-9 * std::max(1.0, std::abs(*z)));
    }
    REQUIRE(feet == out.feet);
}

// ------------------------------------------------------------------ R4.0

TEST_CASE("CF4: with a strip, a point within r(g) of the edge goes in by L12 at its own position and z",
          "[refine_points][constraint_foot][CF4]") {
    // r(g) = 1e-10 cells here (L16's floor). 5e-11 below A-B: L12 runs first
    // and splits the edge at the point itself; R4 never sees it.
    const auto g = cff::geometry();
    const auto dem = step_dem();
    const Begin b = with_z(cff::line_start(4.25), std::vector<double>(6, 60.0));
    const auto strip = strip_of(dem, b);
    const MeshVertex p{7.25, 4.25 + 5e-11};
    REQUIRE(p.row != 4.25);
    cff::ExactStore store{g};
    store.add(p, 80.0f);
    const auto out = run_points(store, &strip, b, options(3.0, true));
    shape(g, b, out);
    REQUIRE(out.feet == 0);
    const auto self = cfo::vertex_at(g, out, cfo::Frac{p.col, p.row}, 1e-13);
    REQUIRE(self.has_value());
    REQUIRE(out.z[*self] == 80.0);
    REQUIRE_FALSE(cfo::vertex_at(g, out, cfo::Frac{7.25, 4.25}, 1e-13).has_value());
}

TEST_CASE("CF4: a point exactly on a constrained edge is not footed", "[refine_points][constraint_foot][CF4]") {
    const auto g = cff::geometry();
    const auto dem = step_dem();
    const Begin b = with_z(cff::line_start(4.25), std::vector<double>(6, 60.0));
    const auto strip = strip_of(dem, b);
    const bool with_strip = GENERATE(false, true);
    CAPTURE(with_strip);
    cff::ExactStore store{g};
    store.add(MeshVertex{7.25, 4.25}, 80.0f);
    const auto out = run_points(store, with_strip ? &strip : nullptr, b, options(3.0, true));
    shape(g, b, out);
    REQUIRE(out.feet == 0);
    const auto self = cfo::vertex_at(g, out, cfo::Frac{7.25, 4.25}, 0.0);
    REQUIRE(self.has_value());
    REQUIRE(out.z[*self] == 80.0);
}

// ------------------------------------------------------------------ tolerance 0 and the bound

namespace {

// A wedge: the square [0, 16]^2 fanned from V (8, 8), with two constraint
// lines from V to the right side, V-E1 along row 8 and V-E2 at 30 degrees
// above it (world). Points in the wedge 0.6 to 1.8 cells from V are within
// delta_p of one line or both, their feet at least 0.5 from V: a point footed
// on one line can find the other next round.
struct Wedge {
    Begin begin;
    cff::ExactStore store;
};

Wedge wedge(std::uint32_t seed) {
    const auto g = cff::geometry();
    const double t30 = std::tan(std::numbers::pi / 6.0);
    // Ring, counter-clockwise in world: BL, BR, E1, E2, TR, TL; V last.
    const std::vector<std::array<double, 2>> cr{{0, 16}, {16, 16}, {16, 8}, {16, 8 - 8 * t30}, {16, 0}, {0, 0}, {8, 8}};
    std::vector<Point2> xy;
    for (const auto& [c, r] : cr) xy.push_back(cff::world(c, r));
    std::vector<TriangleIndices> tris;
    for (std::uint32_t i = 0; i < 6; ++i) tris.push_back({6, i, (i + 1) % 6});
    Start s;
    s.mesh = IndexedMesh2{std::move(xy), std::move(tris), std::vector<std::uint8_t>(6, 0)};
    s.edges = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {0, 5}, {2, 6}, {3, 6}};
    s.masks = {1, 1, 1, 1, 1, 1, 2, 4};
    Wedge w{with_z(s, std::vector<double>(7, 0.0)), cff::ExactStore{g}};
    std::mt19937 gen{seed};
    for (int i = 0; i < 40; ++i) {
        const double r = 0.6 + 1.2 * static_cast<double>(gen() % 1000u) / 1000.0;
        const double phi = (5.0 + 20.0 * static_cast<double>(gen() % 1000u) / 1000.0) * std::numbers::pi / 180.0;
        w.store.add(MeshVertex{8.0 + r * std::cos(phi), 8.0 - r * std::sin(phi)},
                    static_cast<float>(1 + gen() % 900u) / 100.0f);  // 0.01 to 9 m, never on the 0 m start
    }
    return w;
}

}  // namespace

TEST_CASE("CF4: at tolerance 0 every footed point goes in after its foot, the run ends, and no point owns two feet",
          "[refine_points][constraint_foot][CF4]") {
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    CAPTURE(seed);
    const auto g = cff::geometry();
    const Wedge w = wedge(seed);
    const auto out = run_points(w.store, nullptr, w.begin, options(0.0, true));
    shape(g, w.begin, out);
    REQUIRE(j2_over(g, w.store, w.begin, out, 0.0) == 0);  // the fallback: every point within 0 m

    // Feet from the output: inserted vertices that are no stored point (to
    // 1e-12 cells: the world round trip (col, 16 - row) and back rounds).
    // Each is the projection of a stored point on an input line; count per
    // point.
    const auto is_stored = [&](cfo::Frac f) {
        return std::any_of(w.store.points().begin(), w.store.points().end(), [&](const auto& pz) {
            return std::abs(pz.first.col - f.col) <= 1e-12 && std::abs(pz.first.row - f.row) <= 1e-12;
        });
    };
    const std::vector<Point2> sv(w.begin.mesh.vertices().begin(), w.begin.mesh.vertices().end());
    std::vector<std::size_t> owned(w.store.size(), 0);
    std::size_t feet = 0;
    for (std::size_t i = sv.size(); i < out.vertices.size(); ++i) {
        const auto f = cfo::frac(g, out.vertices[i]);
        if (is_stored(f)) continue;
        ++feet;
        for (std::size_t k = 0; k < w.store.size(); ++k) {
            const MeshVertex p = w.store.points()[k].first;
            for (const auto& e : w.begin.edges) {
                const auto a = cfo::frac(g, sv[e[0]]), bb = cfo::frac(g, sv[e[1]]);
                const auto [s, d] = cfo::param_dist(a, bb, cfo::Frac{p.col, p.row});
                if (d >= 0.5 || s <= 0.0 || s >= 1.0) continue;
                const cfo::Frac proj{a.col + s * (bb.col - a.col), a.row + s * (bb.row - a.row)};
                if (std::abs(proj.col - f.col) <= 1e-9 && std::abs(proj.row - f.row) <= 1e-9) ++owned[k];
            }
        }
    }
    REQUIRE(feet == out.feet);
    REQUIRE(out.feet > 0);
    std::size_t footed = 0;
    for (std::size_t k = 0; k < owned.size(); ++k) {
        CAPTURE(k);
        REQUIRE(owned[k] <= 1);  // footed once (R4.2)
        footed += owned[k];
    }
    REQUIRE(out.feet <= footed);  // feet <= footed points
}

// ------------------------------------------------------------------ footed once, on a second line

namespace {

// @architect's fixture for the footed-once guard ("Rulings on the survivors",
// (b)): two parallel constraint lines A1 (0, 0) - B1 (10, 0) and A2 (0, 0.6) -
// B2 (10, 0.6) in world (x, y), the strip between them as (A1, B1, B2) and
// (A1, B2, A2), its two short sides constrained too (mask 1), z 0 m. In
// (col, row) = (x, 16 - y) the lines are rows 16 and 15.4; delta_p = 0.5
// cells. The stored point P (5, 0.35), 5 m, lies in (A1, B2, A2), 0.25 cells
// from A2-B2 and 0.35 from A1-B1. Its foot F (5, 0.6) on A2-B2 makes Lawson
// flip A1-B2, after which P lies in (A1, B1, F), whose A1-B1 gives a second
// Hit at (5, 0): the footed line is cut, but the other line is not.
Start parallel_lines() {
    Start s;
    s.mesh = IndexedMesh2{{cff::world(0, 16), cff::world(10, 16), cff::world(10, 16 - 0.6), cff::world(0, 16 - 0.6)},
                          {{0, 1, 2}, {0, 2, 3}},
                          std::vector<std::uint8_t>(2, 0)};
    s.edges = {{0, 1}, {3, 2}, {0, 3}, {1, 2}};
    s.masks = {2, 4, 1, 1};
    return s;
}

const MeshVertex kTwoLinesP{5.0, 16.0 - 0.35};

// The triangle of m that holds p strictly inside.
std::optional<std::uint32_t> holding(const terrain::mesh::LatticeMesh& m, MeshVertex p) {
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        bool in = true;
        for (unsigned k = 0; k < 3; ++k) in = in && terrain::mesh::orient_sign(m.corner(t, k), m.corner(t, (k + 1) % 3), p) > 0;
        if (in) return t;
    }
    return std::nullopt;
}

}  // namespace

TEST_CASE("CF4: the two-lines fixture is what it claims: after its first foot the point has a Hit on the other line",
          "[refine_points][constraint_foot][CF4]") {
    // Replayed with refine_points' own pieces: to_lattice, legalise_all, the
    // search, the foot's split and legalise_around (seeds as refine_points
    // gives them for an outline edge: the owner and the new triangle).
    namespace mesh = terrain::mesh;
    const auto g = cff::geometry();
    const Start s = parallel_lines();
    auto built = terrain::refinement::detail::to_lattice(g, s.mesh, s.edges, s.masks);
    REQUIRE(std::holds_alternative<mesh::LatticeMesh>(built));
    auto& m = std::get<mesh::LatticeMesh>(built);
    const auto frame = mesh::lattice_frame(1.0, 1.0, g.rows(), g.cols());
    REQUIRE(mesh::legalise_all<terrain::pred::DefaultKernel>(m, frame, [](std::uint32_t) {}) == 0);  // a rectangle: a tie, no flip

    const auto t = holding(m, kTwoLinesP);
    REQUIRE(t.has_value());
    const auto first = mesh::constraint_foot(m, *t, kTwoLinesP, 0.5, frame);
    REQUIRE(first.status == mesh::FootStatus::Hit);
    REQUIRE(first.owner == *t);
    REQUIRE(m.corner(first.owner, first.edge).row == m.corner(first.owner, (first.edge + 1) % 3).row);
    REQUIRE(m.corner(first.owner, first.edge).row != 16.0);  // on A2-B2
    REQUIRE(std::abs(first.at.col - 5.0) <= 1e-12);

    const std::uint32_t before = static_cast<std::uint32_t>(m.triangle_count());
    const auto q = m.split_edge(first.owner, first.edge, first.at);
    const std::array<std::uint32_t, 2> seeds{first.owner, before};
    mesh::FlipStack stack;
    REQUIRE(mesh::legalise_around<terrain::pred::DefaultKernel>(m, q, std::span<const std::uint32_t>{seeds}, frame,
                                                                stack, [](std::uint32_t) {}) >= 1);  // A1-B2 flipped

    const auto t2 = holding(m, kTwoLinesP);
    REQUIRE(t2.has_value());
    const auto second = mesh::constraint_foot(m, *t2, kTwoLinesP, 0.5, frame);
    REQUIRE(second.status == mesh::FootStatus::Hit);  // what the guard is for
    REQUIRE(m.corner(second.owner, second.edge).row == 16.0);  // on A1-B1
    REQUIRE(m.corner(second.owner, (second.edge + 1) % 3).row == 16.0);
    REQUIRE(std::abs(second.at.col - 5.0) <= 1e-12);
}

TEST_CASE("CF4: a point footed on one line is not footed again on a second line; it goes in as itself",
          "[refine_points][constraint_foot][CF4]") {
    // Kills the footed-once guard dropped (the final check's `&& !was_footed`):
    // a second foot at (5, 0), feet 2 and the point owning two feet.
    const auto g = cff::geometry();
    const Begin b = with_z(parallel_lines(), std::vector<double>(4, 0.0));
    cff::ExactStore store{g};
    store.add(kTwoLinesP, 5.0f);
    const auto out = run_points(store, nullptr, b, options(0.0, true));
    shape(g, b, out);
    CHECK(out.feet_refused == 0);
    REQUIRE(cfo::vertex_at(g, out, cfo::Frac{5.0, 16.0 - 0.6}, 1e-12).has_value());  // F, on A2-B2
    REQUIRE_FALSE(cfo::vertex_at(g, out, cfo::Frac{5.0, 16.0}, 1e-9).has_value());   // no foot on A1-B1
    REQUIRE(out.feet == 1);
    REQUIRE(out.feet_fallback == 1);
    const auto self = cfo::vertex_at(g, out, cfo::Frac{kTwoLinesP.col, kTwoLinesP.row}, 1e-12);
    REQUIRE(self.has_value());
    REQUIRE(out.z[*self] == 5.0);
    REQUIRE(j2_over(g, store, b, out, 0.0) == 0);
}

// ------------------------------------------------------------------ R4.2, a foot on a neighbour's edge

namespace {

// CF3's neighbour fixture (prop_constraint_foot_refine.cpp), for the final
// check: nodes at integer world (x, y), x_min -5, y_max 32, dx = dy = 1, 31
// columns and 39 rows, so (col, row) = (x + 5, 32 - y). The only constraint
// P (0, -0.3) - Q (20, -0.3); above it the sliver u = (P, Q, V) with V (4, 0.1),
// then t = (P, V, W) and (V, Q, W) with W (2, 29.7) far up. The stored point
// N (2, 0), 1 m on 0 m ground, lies in t, which has no constrained edge; P-Q,
// u's edge across t's free P-V, is 0.3 cells from N. The foot's split is on
// the outline, and no flip after it reaches t.
RasterGeometry neighbour_geometry() { return RasterGeometry{-5.0, 32.0, 1.0, 1.0, 31, 39}; }

Start neighbour_start() {
    Start s;
    s.mesh = IndexedMesh2{{{0, -0.3}, {20, -0.3}, {4, 0.1}, {2, 29.7}},
                          {{0, 1, 2}, {0, 2, 3}, {2, 1, 3}},
                          std::vector<std::uint8_t>(3, 0)};
    s.edges = {{0, 1}};
    s.masks = {1};
    return s;
}

const MeshVertex kNeighbourN{7.0, 32.0};
const cfo::Frac kNeighbourFoot{7.0, 32.3};

}  // namespace

TEST_CASE("CF4: the neighbour fixture is what it claims: the foot is on the neighbour's edge and no flip touches t",
          "[refine_points][constraint_foot][CF4]") {
    namespace mesh = terrain::mesh;
    const auto g = neighbour_geometry();
    const Start s = neighbour_start();
    auto built = terrain::refinement::detail::to_lattice(g, s.mesh, s.edges, s.masks);
    REQUIRE(std::holds_alternative<mesh::LatticeMesh>(built));
    auto& m = std::get<mesh::LatticeMesh>(built);
    const auto frame = mesh::lattice_frame(1.0, 1.0, g.rows(), g.cols());
    REQUIRE(mesh::legalise_all<terrain::pred::DefaultKernel>(m, frame, [](std::uint32_t) {}) == 0);

    const auto t = holding(m, kNeighbourN);
    REQUIRE(t.has_value());
    for (unsigned k = 0; k < 3; ++k) REQUIRE_FALSE(m.is_constrained(*t, k));
    const auto foot = mesh::constraint_foot(m, *t, kNeighbourN, 0.5, frame);
    REQUIRE(foot.status == mesh::FootStatus::Hit);
    REQUIRE(foot.owner != *t);  // the neighbour's edge
    REQUIRE(std::abs(foot.at.col - kNeighbourFoot.col) <= 1e-12);
    REQUIRE(std::abs(foot.at.row - kNeighbourFoot.row) <= 1e-12);
    REQUIRE(m.neighbours(foot.owner)[foot.edge] == mesh::kNoNeighbour);  // the outline: seeds owner and new

    const auto held = m.triangles()[*t];
    const std::uint32_t before = static_cast<std::uint32_t>(m.triangle_count());
    const auto q = m.split_edge(foot.owner, foot.edge, foot.at);
    const std::array<std::uint32_t, 2> seeds{foot.owner, before};
    mesh::FlipStack stack;
    std::vector<std::uint32_t> written;
    mesh::legalise_around<terrain::pred::DefaultKernel>(m, q, std::span<const std::uint32_t>{seeds}, frame, stack,
                                                        [&](std::uint32_t w) { written.push_back(w); });
    REQUIRE(m.triangles()[*t] == held);  // t unchanged ...
    REQUIRE(std::find(written.begin(), written.end(), *t) == written.end());  // ... and not touched
    REQUIRE(holding(m, kNeighbourN) == t);
}

TEST_CASE("CF4: the foot goes on the neighbour's edge, and the point still goes in from the holding triangle (R4.2)",
          "[refine_points][constraint_foot][CF4]") {
    // Kills the final check's holding triangle dropped from the active set
    // when its foot goes on a neighbour's edge (refine_points.hpp, R4.2's
    // `if (owner != t) skipped.push_back(t)`): t is untouched, so without it N
    // is never scanned again and stays 1 m off.
    const auto g = neighbour_geometry();
    const Begin b = with_z(neighbour_start(), std::vector<double>(4, 0.0));
    cff::ExactStore store{g};
    store.add(kNeighbourN, 1.0f);
    const auto out = run_points(store, nullptr, b, options(0.0, true));
    shape(g, b, out);
    CHECK(out.feet_refused == 0);
    REQUIRE(cfo::vertex_at(g, out, kNeighbourFoot, 1e-12).has_value());  // the foot, on P-Q
    REQUIRE(out.feet == 1);
    const auto self = cfo::vertex_at(g, out, cfo::Frac{kNeighbourN.col, kNeighbourN.row}, 0.0);
    REQUIRE(self.has_value());  // N, after its foot
    REQUIRE(out.z[*self] == 1.0);
    REQUIRE(out.feet_fallback == 1);
    REQUIRE(j2_over(g, store, b, out, 0.0) == 0);  // the tolerance oracle over every stored point
}

// ------------------------------------------------------------------ 14b T3 and T6, rule on

TEST_CASE("CF4 T3 T6: refine_points with a strip and feet on keeps J2 and is bit-identical over threads",
          "[refine_points][constraint_foot][CF4][T6]") {
    const auto g = cff::geometry();
    const Raster<float> dem{g, refinement_fixtures::rough_dem(cff::kN, cff::kN, 7)};
    const Begin b = phase1(dem, cff::line_start(4.3), 5.0);
    const auto strip = strip_of(dem, b);
    const cff::ExactStore store = cff::every_node(dem);
    const auto ref = run_points(store, &strip, b, options(0.5, true));
    shape(g, b, ref);
    REQUIRE(j2_over(g, store, b, ref, 0.5) == 0);
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    const auto other = run_points(store, &strip, b, options(0.5, true, threads));
    REQUIRE(refine_digest::digest(other) == refine_digest::digest(ref));
    REQUIRE(other.feet == ref.feet);
    REQUIRE(other.feet_fallback == ref.feet_fallback);
    REQUIRE(other.feet_refused == ref.feet_refused);
}

TEST_CASE("CF4 T3 T6: refine_strip with feet on keeps the node oracle and is bit-identical over threads",
          "[refine_strip][constraint_foot][CF4][T6]") {
    // Phase 1 at the strip run's own tolerance: refine_strip rescans DEM
    // nodes only in the slots it writes (15f L3), so E2 holds over the whole
    // mesh only when phase 1 met the same tolerance.
    const auto g = cff::geometry();
    const Raster<float> dem{g, refinement_fixtures::rough_dem(cff::kN, cff::kN, 7)};
    const Begin b = phase1(dem, cff::line_start(4.3), 0.5);
    const auto strip = strip_of(dem, b);
    const auto ref = run_strip(dem, strip, b, options(0.5, true));
    shape(g, b, ref);
    REQUIRE(cfo::node_findings(dem, ref, 0.5).over == 0);
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    const auto other = run_strip(dem, strip, b, options(0.5, true, threads));
    REQUIRE(refine_digest::digest(other) == refine_digest::digest(ref));
    REQUIRE(other.feet == ref.feet);
    REQUIRE(other.feet_fallback == ref.feet_fallback);
}

// ------------------------------------------------------------------ the off switch

TEST_CASE("CF4 off switch: constraint_feet false is master's final check on both paths",
          "[refine_points][refine_strip][constraint_foot][off]") {
    const auto g = cff::geometry();
    const Raster<float> dem{g, refinement_fixtures::rough_dem(cff::kN, cff::kN, 7)};
    const Begin b = phase1(dem, cff::line_start(4.3), 5.0);
    const auto strip = strip_of(dem, b);
    const cff::ExactStore store = cff::every_node(dem);
    const auto points = run_points(store, &strip, b, options(0.5, false));
    REQUIRE(points.ok());
    REQUIRE(refine_digest::topology_digest(points) == 0xa460db59ad797984ull);  // RECORDED, see the header
    REQUIRE(points.feet == 0);
    REQUIRE(points.feet_fallback == 0);
    const auto strip_run = run_strip(dem, strip, b, options(0.5, false));
    REQUIRE(strip_run.ok());
    REQUIRE(refine_digest::topology_digest(strip_run) == 0x93e787a9ab2b77fdull);  // RECORDED, see the header
    REQUIRE(strip_run.feet == 0);
}
