// Increment 23b (docs/increments/23-basin-scale.md, "The seam protocol" step 4,
// K1, K2, the degeneracy policy, and "Tests @tester can write red", FE1 to
// FE5): refinement with frozen edges, through refine, refine_points and
// refine_strip. FE2 to FE5 are invariant-critical (mutation runs at green).
//
// Interface, as the design names it and as "Settled after 23b's red step"
// rules what it left open (N1-N19; here N2, N4 to N7 and N16):
//
//   RefineOptions::frozen_mask        std::uint32_t, default 0; set by member assignment
//   PointRefineOptions::frozen_mask   std::uint32_t, default 0; the same, for
//                                     refine_points and refine_strip
//   PointRefineOutcome::on_frozen            stored check points on a frozen
//                                            edge (exactly; with a strip, within
//                                            r(g) too, N7), each counted once
//   PointRefineOutcome::on_frozen_max_error  their largest |z - the frozen edge's
//                                            linear z there|, ends valid; the
//                                            edge is never split, so "there" is
//                                            between the start edge's two ends
//   RefineOutcome::quality_skipped    includes QualityOutcome::skipped_frozen
//                                     ("every reason summed")
//   Feet (FE3): a frozen edge is never a foot's edge; the node goes in itself
//     and is not counted in feet_refused (no foot was refused; none was taken)
//
// Every check is an oracle on the output (support/frozen_oracle.hpp,
// support/strip_oracle.hpp), shown able to fail in
// unit/test_refinement_frozen_oracle.cpp:
//   K2       frozen_findings: every frozen start edge is an output constraint
//            edge between the same two vertices with the same mask, its ends
//            unmoved, and no output vertex on or within 1e-9 cells of it;
//   §3D      the tolerance oracle (node_findings_off_frozen: every valid DEM
//            node not exactly on a frozen edge, every closed triangle with three
//            valid vertices, plane recomputed from the output), and the
//            constrained-Delaunay oracle in the producer's frame
//            (delaunay_violations, frozen edges being constraints);
//   K1       a frozen mask that meets no edge gives master's output bit for bit
//            (FE1; the existing suites' golden digests are the other half).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>

#include "feet_fixtures.hpp"
#include "frozen_oracle.hpp"
#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using frozen_oracle::kFeature;
using frozen_oracle::kNoEdgeHasThis;
using frozen_oracle::kOutline;
using frozen_oracle::kSeam;
using frozen_oracle::Lat;
using frozen_oracle::Side;
using strip_oracle::Edges;
using strip_oracle::Mesh;
using strip_oracle::Start;
using strip_oracle::Terrain;
using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::CheckPoints;
using terrain::refinement::constraint_check_points;
using terrain::refinement::ConstraintCheckPoints;
using terrain::refinement::PointRefineOptions;
using terrain::refinement::PointRefineOutcome;
using terrain::refinement::refine;
using terrain::refinement::refine_points;
using terrain::refinement::refine_strip;
using terrain::refinement::RefineOptions;
using terrain::refinement::RefineOutcome;

namespace {

using EdgeSpan = std::span<const std::array<std::uint32_t, 2>>;
using MaskSpan = std::span<const std::uint32_t>;

struct Knobs {
    double tol = 0.5;
    std::uint32_t frozen = 0;
    bool feet = false;
    double angle = 0.0;
    unsigned threads = 1;
};

RefineOutcome run(const Raster<float>& dem, const Start& s, const Knobs& k) {
    RefineOptions o;
    o.tolerance = k.tol;
    o.threads = k.threads;
    o.min_angle_deg = k.angle;
    o.constraint_feet = k.feet;
    o.frozen_mask = k.frozen;
    return refine(dem, s.mesh, EdgeSpan{s.edges}, MaskSpan{s.masks}, o);
}

// refine with frozen_mask never named: master's call, byte for byte.
RefineOutcome run_master(const Raster<float>& dem, const Start& s, const Knobs& k) {
    RefineOptions o;
    o.tolerance = k.tol;
    o.threads = k.threads;
    o.min_angle_deg = k.angle;
    o.constraint_feet = k.feet;
    return refine(dem, s.mesh, EdgeSpan{s.edges}, MaskSpan{s.masks}, o);
}

PointRefineOptions point_options(double tol, std::uint32_t frozen, unsigned threads = 1) {
    PointRefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.frozen_mask = frozen;
    return o;
}

// The start of a refine_points / refine_strip run: a mesh with its z.
struct Begin {
    IndexedMesh2 mesh;
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    Edges edges;
    std::vector<std::uint32_t> masks;
};

Begin begin_from(const Raster<float>& dem, const Start& s) {
    auto [z, valid] = strip_oracle::start_z(dem, s.mesh.vertices());
    return Begin{s.mesh, std::move(z), std::move(valid), s.edges, s.masks};
}

Begin begin_from(const RefineOutcome& out) {
    return Begin{IndexedMesh2{out.vertices, out.triangles, std::vector<std::uint8_t>(out.triangles.size(), 0)},
                 out.z, out.valid, out.edges, out.masks};
}

PointRefineOutcome run_points(const CheckPoints& store, const Begin& b, const PointRefineOptions& o,
                              const ConstraintCheckPoints* strip = nullptr) {
    return refine_points(store, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                         EdgeSpan{b.edges}, MaskSpan{b.masks}, o, strip);
}

PointRefineOutcome run_strip(const Raster<float>& dem, const ConstraintCheckPoints& strip, const Begin& b,
                             const PointRefineOptions& o) {
    return refine_strip(dem, strip, b.mesh, std::span<const double>{b.z}, std::span<const std::uint8_t>{b.valid},
                        EdgeSpan{b.edges}, MaskSpan{b.masks}, o);
}

// The strip a caller builds under 23b: the non-frozen constraint edges only
// ("the edge strip makes no check points on frozen edges"; 15f: "a start
// constraint edge with no strip edge is allowed: that is how 23b leaves frozen
// edges out").
ConstraintCheckPoints strip_of(const Raster<float>& dem, const Begin& b, std::uint32_t frozen) {
    Edges open;
    for (std::size_t k = 0; k < b.edges.size(); ++k)
        if ((b.masks[k] & frozen) == 0) open.push_back(b.edges[k]);
    return constraint_check_points(dem, b.mesh.vertices(), EdgeSpan{open});
}

CheckPoints store_of(const RasterGeometry& g, const std::vector<std::pair<Lat, double>>& pts) {
    CheckPoints cp{g};
    std::vector<Point2> xy;
    std::vector<float> z;
    for (const auto& [p, h] : pts) {
        xy.push_back(strip_oracle::world(g, p.col, p.row));
        z.push_back(static_cast<float>(h));
    }
    cp.add(std::span<const Point2>{xy}, std::span<const float>{z});
    cp.freeze();
    return cp;
}

// Random dyadic positions inside the jittered outline, z the DEM's there.
std::vector<std::pair<Lat, double>> random_points(const Raster<float>& dem, std::uint32_t seed, std::size_t n) {
    std::mt19937 gen{seed * 31u + 7u};
    std::vector<std::pair<Lat, double>> out;
    while (out.size() < n) {
        const Lat p{4.0 + static_cast<double>(gen() % 1537u) / 64.0, 4.0 + static_cast<double>(gen() % 1537u) / 64.0};
        if (const auto z = strip_oracle::height(dem, p)) out.push_back({p, *z});
    }
    return out;
}

// K2, the §3D oracles and the DEM-node guarantee off the frozen edges, on one
// output against the start the run was given.
void check(const Raster<float>& dem, const Start& start, std::uint32_t frozen, const Mesh& out, double tol) {
    const RasterGeometry& g = dem.geometry();
    const auto& sv = start.mesh.vertices();
    const auto f = frozen_oracle::frozen_findings(g, sv, start.edges, start.masks, frozen, out);
    CAPTURE(f.frozen_edges, f.missing, f.wrong_mask, f.moved, f.on, f.near);
    REQUIRE(f.frozen_edges > 0);
    REQUIRE(frozen_oracle::clean(f));
    const auto segs = frozen_oracle::frozen_segments(g, sv, start.edges, start.masks, frozen);
    const auto nf = frozen_oracle::node_findings_off_frozen(dem, out, tol, segs);
    CAPTURE(nf.over, nf.worst, nf.not_ccw);
    REQUIRE(nf.not_ccw == 0);
    REQUIRE(nf.over == 0);
    REQUIRE(strip_oracle::delaunay_violations(g, out) == 0);
}

template <class A, class B>
void require_same_mesh(const A& a, const B& b) {
    REQUIRE(a.vertices.size() == b.vertices.size());
    for (std::size_t i = 0; i < a.vertices.size(); ++i) {
        REQUIRE(a.vertices[i].x == b.vertices[i].x);
        REQUIRE(a.vertices[i].y == b.vertices[i].y);
    }
    REQUIRE(a.z == b.z);
    REQUIRE(a.valid == b.valid);
    REQUIRE(a.triangles == b.triangles);
    REQUIRE(a.edges == b.edges);
    REQUIRE(a.masks == b.masks);
}

void require_same(const RefineOutcome& a, const RefineOutcome& b) {
    REQUIRE(a.status == b.status);
    require_same_mesh(a, b);
    REQUIRE(a.rounds == b.rounds);
    REQUIRE(a.inserted == b.inserted);
    REQUIRE(a.flips == b.flips);
    REQUIRE(a.max_error == b.max_error);
    REQUIRE(a.uncovered == b.uncovered);
    REQUIRE(a.carved == b.carved);
    REQUIRE(a.quality_inserted == b.quality_inserted);
    REQUIRE(a.quality_skipped == b.quality_skipped);
    REQUIRE(a.feet == b.feet);
    REQUIRE(a.feet_refused == b.feet_refused);
}

void require_same(const PointRefineOutcome& a, const PointRefineOutcome& b) {
    require_same(static_cast<const RefineOutcome&>(a), static_cast<const RefineOutcome&>(b));
    REQUIRE(a.coincident == b.coincident);
    REQUIRE(a.coincident_max_error == b.coincident_max_error);
    REQUIRE(a.strip_points == b.strip_points);
    REQUIRE(a.strip_inserted == b.strip_inserted);
    REQUIRE(a.strip_max_error == b.strip_max_error);
    REQUIRE(a.strip_refused == b.strip_refused);
    REQUIRE(a.strip_refused_max_error == b.strip_refused_max_error);
    REQUIRE(a.nodes_inserted == b.nodes_inserted);
    REQUIRE(a.on_frozen == 0);
    REQUIRE(b.on_frozen == 0);
}

bool has_vertex(const RasterGeometry& g, const Mesh& m, Lat p) {
    return std::any_of(m.vertices.begin(), m.vertices.end(), [&](Point2 v) { return strip_oracle::lat(g, v) == p; });
}

// ------------------------------------------------------------------ fixtures

const RasterGeometry kSquare = strip_oracle::exact_geometry(9, 9);

// Smooth ground over a 9 x 9 lattice, dyadic values, plus named spikes.
Raster<float> ground(const RasterGeometry& g, std::vector<std::pair<std::pair<std::size_t, std::size_t>, float>> spikes,
                     std::optional<float> nodata = std::nullopt) {
    std::vector<float> v(g.rows() * g.cols());
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c)
            v[r * g.cols() + c] = static_cast<float>(
                std::round((20.0 + 3.0 * std::sin(0.9 * static_cast<double>(c)) + 2.0 * std::cos(0.7 * static_cast<double>(r))) * 64.0)
                / 64.0);
    for (const auto& [rc, z] : spikes) v[rc.first * g.cols() + rc.second] = z;
    if (nodata) return Raster<float>{g, std::move(v), *nodata};
    return Raster<float>{g, std::move(v)};
}

// The 9 x 9 square cut by the grid-line seam col 4, row 0 to row 8.
Start grid_seam(std::uint32_t seam_mask = kSeam) {
    return frozen_oracle::seam_domain(kSquare, 8, 8, 4, 4, {}, Side::Both, seam_mask);
}

}  // namespace

// ------------------------------------------------------------------------ FE1

TEST_CASE("FE1: refine with a frozen mask that meets no edge is bit-identical to master", "[refinement][frozen][k1]") {
    // K1: "with frozen_mask 0, every mesh is bit-identical to master's". A
    // mask meeting no edge is the empty frozen set; 0 is master's default.
    const std::uint32_t seed = GENERATE(1u, 2u);
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth, Terrain::SmoothWithHole);
    const double tol = GENERATE(0.0, 0.5, 5.0);
    const bool feet = GENERATE(false, true);
    const double angle = GENERATE(0.0, 25.0);
    CAPTURE(seed, static_cast<int>(kind), tol, feet, angle);
    const auto g = strip_oracle::exact_geometry();
    const auto dem = strip_oracle::terrain_dem(g, kind, seed);
    const auto start = strip_oracle::jittered(g, seed);
    const Knobs k{tol, 0, feet, angle, 1};
    const RefineOutcome master = run_master(dem, start, k);
    REQUIRE(master.ok());
    require_same(master, run(dem, start, k));
    require_same(master, run(dem, start, Knobs{tol, kNoEdgeHasThis, feet, angle, 1}));
}

TEST_CASE("FE1: refine_points and refine_strip with an empty frozen set are bit-identical to master",
          "[refinement][frozen][k1]") {
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    const double tol = GENERATE(0.0, 0.5);
    CAPTURE(seed, tol);
    const auto g = strip_oracle::exact_geometry();
    const auto dem = strip_oracle::terrain_dem(g, Terrain::Smooth, seed);
    const auto start = strip_oracle::jittered(g, seed);
    const RefineOutcome phase1 = run_master(dem, start, Knobs{2.0});
    REQUIRE(phase1.ok());
    const Begin b = begin_from(phase1);
    const auto strip = strip_of(dem, b, 0);
    const CheckPoints store = store_of(g, random_points(dem, seed, 200));

    PointRefineOptions master;  // frozen_mask never named
    master.tolerance = tol;
    master.threads = 1;
    for (const ConstraintCheckPoints* s : {static_cast<const ConstraintCheckPoints*>(nullptr), &strip}) {
        CAPTURE(s != nullptr);
        const auto m = run_points(store, b, master, s);
        REQUIRE(m.ok());
        require_same(m, run_points(store, b, point_options(tol, kNoEdgeHasThis), s));
    }
    const auto m = run_strip(dem, strip, b, master);
    REQUIRE(m.ok());
    require_same(m, run_strip(dem, strip, b, point_options(tol, kNoEdgeHasThis)));
}

// ------------------------------------------------------------------------ FE2

TEST_CASE("FE2: a node exactly on a frozen grid-line seam is never inserted", "[refinement][frozen][fe2]") {
    // A 50 m spike on node (row 4, col 4), on the seam. Unfrozen it goes in
    // on the seam (the control, so the test can fail); frozen it never does,
    // and every other node ends within tolerance (the oracle skips only the
    // seam's own nodes, which are the seam pass's).
    const double tol = GENERATE(0.0, 0.5, 2.0);
    const bool feet = GENERATE(false, true);
    const double angle = GENERATE(0.0, 25.0);
    CAPTURE(tol, feet, angle);
    const auto dem = ground(kSquare, {{{4, 4}, 70.0f}});
    const auto start = grid_seam();

    const auto control = strip_oracle::mesh_of(run(dem, start, Knobs{tol, 0, feet, angle}));
    REQUIRE(has_vertex(kSquare, control, Lat{4, 4}));

    const RefineOutcome out = run(dem, start, Knobs{tol, kSeam, feet, angle});
    REQUIRE(out.ok());
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE_FALSE(has_vertex(kSquare, m, Lat{4, 4}));
    check(dem, start, kSeam, m, tol);
    REQUIRE(out.max_error <= tol);  // the seam's nodes are not in it
}

TEST_CASE("FE2: freezing one edge does not stop work beside it", "[refinement][frozen][fe2]") {
    // The spike one column off the seam goes in, frozen or not.
    const auto dem = ground(kSquare, {{{4, 3}, 70.0f}});
    const auto start = grid_seam();
    const RefineOutcome out = run(dem, start, Knobs{0.5, kSeam});
    REQUIRE(out.ok());
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE(has_vertex(kSquare, m, Lat{3, 4}));
    check(dem, start, kSeam, m, 0.5);
}

TEST_CASE("FE2: a seam along a feature edge is frozen by the seam's bit and keeps both bits",
          "[refinement][frozen][fe2]") {
    // Degeneracy policy: "one edge with both bits; the feature keeps its bits,
    // and the seam's rules apply".
    const auto dem = ground(kSquare, {{{4, 4}, 70.0f}});
    const auto start = grid_seam(kSeam | kFeature);
    const RefineOutcome out = run(dem, start, Knobs{0.5, kSeam});
    REQUIRE(out.ok());
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE_FALSE(has_vertex(kSquare, m, Lat{4, 4}));
    check(dem, start, kSeam, m, 0.5);  // wrong_mask would catch a dropped feature bit
}

TEST_CASE("FE2: a general seam through DEM nodes is never split at them", "[refinement][frozen][fe2]") {
    // tester.md §3A, seams through nodes: T (2, 0) to B (6, 8) passes exactly
    // through nodes (3, 2), (4, 4) and (5, 6) (col, row); each carries a spike.
    const auto dem = ground(kSquare, {{{2, 3}, 60.0f}, {{4, 4}, -40.0f}, {{6, 5}, 55.0f}});
    const auto start = frozen_oracle::seam_domain(kSquare, 8, 8, 2, 6, {}, Side::Both);
    const double tol = GENERATE(0.0, 1.0);
    const bool feet = GENERATE(false, true);
    CAPTURE(tol, feet);
    const auto control = strip_oracle::mesh_of(run(dem, start, Knobs{tol, 0, feet}));
    REQUIRE(frozen_oracle::frozen_findings(kSquare, start.mesh.vertices(), start.edges, start.masks, kSeam, control).on
            >= 1);
    const RefineOutcome out = run(dem, start, Knobs{tol, kSeam, feet});
    REQUIRE(out.ok());
    check(dem, start, kSeam, strip_oracle::mesh_of(out), tol);
}

TEST_CASE("FE2: a void triangle beside a seam never carves on it", "[refinement][frozen][fe2]") {
    // The seam's top end (row 0, col 4) is NoData: every triangle around it is
    // void and carves the valid node nearest it, which is the seam's own node
    // (row 1, col 4) unless the seam is frozen. The stopping rule still
    // empties every void triangle of valid nodes off the seam.
    const auto dem = ground(kSquare, {{{0, 4}, -9999.0f}}, -9999.0f);
    const auto start = grid_seam();
    const auto control = strip_oracle::mesh_of(run(dem, start, Knobs{0.5, 0}));
    REQUIRE(has_vertex(kSquare, control, Lat{4, 1}));
    const RefineOutcome out = run(dem, start, Knobs{0.5, kSeam});
    REQUIRE(out.ok());
    REQUIRE(out.uncovered == 0);
    check(dem, start, kSeam, strip_oracle::mesh_of(out), 0.5);
}

TEST_CASE("FE2: the jittered sweep keeps K2, the tolerance and Delaunay on every path", "[refinement][frozen][fe2]") {
    // tester.md §3A on the 15f fixture: chain 2 runs ALONG grid line 16,
    // chain 8 passes THROUGH six nodes per edge, chain 1 is the jittered
    // outline (a seam along the outline). §3D on every path: feet on/off,
    // quality on/off, tolerance 0 included, NoData included.
    const std::uint32_t frozen = GENERATE(2u, 8u, 2u | 8u, 1u | 2u | 4u | 8u);
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth, Terrain::SmoothWithHole);
    const double tol = GENERATE(0.0, 0.5, 5.0);
    const bool feet = GENERATE(false, true);
    const double angle = GENERATE(0.0, 25.0);
    const std::uint32_t seed = GENERATE(1u, 2u);
    CAPTURE(frozen, static_cast<int>(kind), tol, feet, angle, seed);
    const auto g = strip_oracle::exact_geometry();
    const auto dem = strip_oracle::terrain_dem(g, kind, seed);
    const auto start = strip_oracle::jittered(g, seed);
    const RefineOutcome out = run(dem, start, Knobs{tol, frozen, feet, angle});
    REQUIRE(out.ok());
    REQUIRE(out.uncovered == 0);
    REQUIRE(out.max_error <= tol);
    check(dem, start, frozen, strip_oracle::mesh_of(out), tol);
}

TEST_CASE("FE2: the output with frozen edges does not depend on the thread count", "[refinement][frozen][fe2]") {
    const std::uint32_t seed = GENERATE(1u, 2u);
    CAPTURE(seed);
    const auto g = strip_oracle::exact_geometry();
    const auto dem = strip_oracle::terrain_dem(g, Terrain::Rough, seed);
    const auto start = strip_oracle::jittered(g, seed);
    const RefineOutcome one = run(dem, start, Knobs{0.5, 2u | 8u, true, 25.0, 1});
    REQUIRE(one.ok());
    for (const unsigned threads : {2u, 8u}) {
        CAPTURE(threads);
        require_same(one, run(dem, start, Knobs{0.5, 2u | 8u, true, 25.0, threads}));
    }
}

TEST_CASE("FE2: a seam at the far edge of a 16,385-node lattice at UTM scale", "[refinement][frozen][fe2][scale]") {
    // tester.md §3A, extreme scales: world coordinates near (5e5, 7e6), cells
    // 10 m by 5 m, a lattice wider than 8,191 nodes so r(g) is above its
    // floor (L16), the seam on column 16,380 next to the node rectangle's
    // right edge.
    const RasterGeometry g{500000.0, 7000000.0, 10.0, 5.0, 16385, 9};
    std::vector<float> v(g.rows() * g.cols());
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c)
            v[r * g.cols() + c] = static_cast<float>(300.0 + 0.25 * static_cast<double>(c % 7) + 0.5 * static_cast<double>(r));
    v[4 * g.cols() + 16380] = 400.0f;  // on the seam
    v[4 * g.cols() + 16378] = 380.0f;  // beside it
    const Raster<float> dem{g, std::move(v)};
    const auto start = frozen_oracle::seam_domain(g, 8, 8, 16380, 16380, {}, Side::Both, kSeam, 16376);
    const RefineOutcome out = run(dem, start, Knobs{0.5, kSeam, true});
    REQUIRE(out.ok());
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE_FALSE(has_vertex(g, m, Lat{16380, 4}));
    REQUIRE(has_vertex(g, m, Lat{16378, 4}));
    check(dem, start, kSeam, m, 0.5);
}

// ------------------------------------------------------------------------ FE3

TEST_CASE("FE3: a node within eps of a frozen edge goes in itself, not as a foot", "[refinement][frozen][fe3]") {
    // 20b's needle: node (row 16, col 8) is 0.0035 cells from the tilted side
    // A-B (mask 1). With feet on and A-B unfrozen, a foot lands on A-B (the
    // control); with A-B frozen, nothing lands on it and the node itself is a
    // vertex. The other sides (masks 2, 4, 8) may still take feet.
    const double tol = GENERATE(0.1, 0.5);
    const double steep = GENERATE(0.0, 1.0);
    CAPTURE(tol, steep);
    const auto dem = feet_fixtures::needle_dem(6.0, steep);
    const RasterGeometry& g = dem.geometry();
    const auto q = feet_fixtures::needle_start(g);
    const Start start{q.mesh, q.edges, q.masks};

    const auto control = strip_oracle::mesh_of(run(dem, start, Knobs{tol, 0, true}));
    const auto c = frozen_oracle::frozen_findings(g, start.mesh.vertices(), start.edges, start.masks, 1u, control);
    REQUIRE(c.on + c.near >= 1);

    const RefineOutcome out = run(dem, start, Knobs{tol, 1u, true});
    REQUIRE(out.ok());
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE(has_vertex(g, m, Lat{static_cast<double>(feet_fixtures::kNeedleCol),
                                 static_cast<double>(feet_fixtures::kNeedleRow)}));
    check(dem, start, 1u, m, tol);
}

TEST_CASE("FE3: a frozen edge in the triangle does not stop a foot on its other constraint",
          "[refinement][frozen][fe3]") {
    // N5: a frozen edge is skipped as an unconstrained one is; another
    // non-frozen constrained edge of the same triangle that qualifies still
    // takes its foot. One start triangle (A, B, E): A-B is 20b's needle side
    // (mask 1), 0.0035 cells from node (row 16, col 8), the bump's peak and the
    // triangle's worst node in the first round; B-E is frozen (mask 4); E-A is
    // mask 1. The foot lands on A-B in that first round, in a triangle that
    // holds the frozen edge.
    const double tol = GENERATE(0.1, 0.5);
    const double steep = GENERATE(0.0, 1.0);
    CAPTURE(tol, steep);
    const auto dem = feet_fixtures::needle_dem(6.0, steep);
    const RasterGeometry& g = dem.geometry();
    const auto ring = feet_fixtures::needle_ring();
    std::vector<Point2> xy{strip_oracle::world(g, ring[0][0], ring[0][1]),
                           strip_oracle::world(g, ring[1][0], ring[1][1]), strip_oracle::world(g, 14.3, 16.2)};
    Start start;
    start.mesh = IndexedMesh2{std::move(xy), {{0, 1, 2}}, {0}};
    start.edges = {{0, 1}, {1, 2}, {0, 2}};
    start.masks = {1u, kSeam, 1u};

    const RefineOutcome out = run(dem, start, Knobs{tol, kSeam, true});
    REQUIRE(out.ok());
    REQUIRE(out.feet >= 1);
    const Mesh m = strip_oracle::mesh_of(out);
    const auto ab = frozen_oracle::frozen_findings(g, start.mesh.vertices(), start.edges, start.masks, 1u, m);
    REQUIRE(ab.on + ab.near >= 1);  // a vertex on A-B: the foot
    check(dem, start, kSeam, m, tol);
}

// ------------------------------------------------------------------------ FE4

TEST_CASE("FE4: refine's quality pass never splits a frozen edge and sums skipped_frozen into quality_skipped",
          "[refinement][frozen][fe4]") {
    // unit/test_mesh_quality.cpp's Q3, mirrored across its hypotenuse
    // (0, 7)-(20, 7): two right triangles whose circumcentres are node (10, 7)
    // on the shared edge, a seam. Flat ground, so refine itself inserts
    // nothing; only the quality pass could split the seam.
    const RasterGeometry g{0.0, 0.0, 1.0, 1.0, 21, 14};
    const Raster<float> dem{g, std::vector<float>(21 * 14, 5.0f)};
    std::vector<Point2> xy;
    for (const auto& [c, r] : std::vector<std::array<double, 2>>{{0, 7}, {2, 13}, {20, 7}, {18, 1}})
        xy.push_back(strip_oracle::world(g, c, r));
    Start start;
    start.mesh = IndexedMesh2{std::move(xy), {{0, 1, 2}, {2, 3, 0}}, {0, 0}};
    start.edges = {{0, 1}, {1, 2}, {2, 3}, {0, 3}, {0, 2}};
    start.masks = {kOutline, kOutline, kOutline, kOutline, kSeam};

    const RefineOutcome control = run(dem, start, Knobs{1.0, 0, false, 25.0});
    REQUIRE(control.quality_inserted >= 1);
    REQUIRE(frozen_oracle::frozen_findings(g, start.mesh.vertices(), start.edges, start.masks, kSeam,
                                           strip_oracle::mesh_of(control))
                .on
            >= 1);

    const RefineOutcome out = run(dem, start, Knobs{1.0, kSeam, false, 25.0});
    REQUIRE(out.ok());
    REQUIRE(out.quality_inserted == 0);
    REQUIRE(out.quality_skipped >= 2);  // both triangles' skips, summed in
    check(dem, start, kSeam, strip_oracle::mesh_of(out), 1.0);
}

// ------------------------------------------------------------------------ FE5

TEST_CASE("FE5: refine_points never inserts a check point on a frozen edge and reports it", "[refinement][frozen][fe5]") {
    // The square [0, 8]^2 cut by the seam col 4, z 0 at every start vertex.
    // Three check points, all 100 m off: (4, 2.5) exactly on the seam, (2, 5)
    // inside the left piece, (6, 8) on the bottom outline (mask 1, not
    // frozen). Unfrozen all three go in (the control); frozen the seam's is
    // skipped and reported with its error, the other two go in.
    const auto start = grid_seam();
    Begin b{start.mesh, std::vector<double>(start.mesh.vertices().size(), 0.0),
            std::vector<std::uint8_t>(start.mesh.vertices().size(), 1), start.edges, start.masks};
    const CheckPoints store = store_of(kSquare, {{Lat{4, 2.5}, 100.0}, {Lat{2, 5}, 100.0}, {Lat{6, 8}, 100.0}});

    const auto control = run_points(store, b, point_options(1.0, 0));
    REQUIRE(control.ok());
    REQUIRE(control.on_frozen == 0);
    REQUIRE(has_vertex(kSquare, strip_oracle::mesh_of(control), Lat{4, 2.5}));

    const auto out = run_points(store, b, point_options(1.0, kSeam));
    REQUIRE(out.ok());
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE_FALSE(has_vertex(kSquare, m, Lat{4, 2.5}));
    REQUIRE(has_vertex(kSquare, m, Lat{2, 5}));
    REQUIRE(has_vertex(kSquare, m, Lat{6, 8}));
    REQUIRE(out.on_frozen == 1);  // once, though two triangles hold it
    REQUIRE(out.on_frozen_max_error == 100.0);
    const auto f = frozen_oracle::frozen_findings(kSquare, start.mesh.vertices(), start.edges, start.masks, kSeam, m);
    REQUIRE(frozen_oracle::clean(f));
    REQUIRE(strip_oracle::delaunay_violations(kSquare, m) == 0);
}

TEST_CASE("FE5: on_frozen counts every point on a frozen chain with its error against the chain",
          "[refinement][frozen][fe5]") {
    // The jittered fixture with chains 2 (along row 16) and 8 (the diagonal
    // through nodes) frozen. Random points plus points planted EXACTLY on the
    // two chains, between their vertices. The count and the largest error are
    // recomputed here: each planted point against the linear z of the start
    // edge holding it (frozen edges are never split).
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    const double tol = GENERATE(0.0, 0.5);
    CAPTURE(seed, tol);
    const auto g = strip_oracle::exact_geometry();
    const auto dem = strip_oracle::terrain_dem(g, Terrain::Smooth, seed);
    const auto start = strip_oracle::jittered(g, seed);
    const Begin b = begin_from(dem, start);
    const std::uint32_t frozen = 2u | 8u;

    // Planted first, so a random point at a planted position is dropped here
    // rather than by the store; every z rounded to the store's float.
    std::map<Lat, double> unique;
    std::size_t planted = 0;
    std::mt19937 gen{seed};
    std::vector<std::array<Lat, 2>> segs;
    for (std::size_t k = 0; k < start.edges.size(); ++k) {
        if ((start.masks[k] & frozen) == 0) continue;
        const Lat a = strip_oracle::lat(g, start.mesh.vertices()[start.edges[k][0]]);
        const Lat c = strip_oracle::lat(g, start.mesh.vertices()[start.edges[k][1]]);
        segs.push_back({a, c});
        for (const double t : {0.25, 0.5, 0.75}) {
            const Lat p{a.col + t * (c.col - a.col), a.row + t * (c.row - a.row)};
            if (!frozen_oracle::on_open(a, c, p)) continue;  // only points exactly on the edge
            planted += unique.emplace(p, 40.0 + static_cast<double>(gen() % 64u)).second ? 1 : 0;
        }
    }
    REQUIRE(planted >= 12);
    for (const auto& [p, z] : random_points(dem, seed, 150)) unique.emplace(p, static_cast<double>(static_cast<float>(z)));

    // The oracle's count and error, over every stored point (a random one can
    // land exactly on row 16 or the diagonal too): exact incidence with the
    // open start edge, linear z along it by the projection parameter.
    std::vector<std::pair<Lat, double>> pts(unique.begin(), unique.end());
    std::size_t expect = 0;
    double worst = 0.0;
    for (const auto& [p, z] : pts)
        for (std::size_t k = 0, j = 0; k < start.edges.size(); ++k) {
            if ((start.masks[k] & frozen) == 0) continue;
            const auto& seg = segs[j++];
            if (!frozen_oracle::on_open(seg[0], seg[1], p)) continue;
            ++expect;
            const double t = strip_oracle::param_dist(seg[0], seg[1], p).first;
            const double za = b.z[start.edges[k][0]], zc = b.z[start.edges[k][1]];
            worst = std::max(worst, std::abs(z - (za + t * (zc - za))));
        }
    REQUIRE(expect >= planted);
    const CheckPoints store = store_of(g, pts);
    const auto out = run_points(store, b, point_options(tol, frozen));
    REQUIRE(out.ok());
    REQUIRE(out.on_frozen == expect);
    REQUIRE(std::abs(out.on_frozen_max_error - worst) <= 1e-9 * std::max(1.0, worst));
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE(frozen_oracle::clean(
        frozen_oracle::frozen_findings(g, start.mesh.vertices(), start.edges, start.masks, frozen, m)));
    REQUIRE(strip_oracle::delaunay_violations(g, m) == 0);
}

TEST_CASE("FE5: refine_strip and refine_points with a strip keep the frozen edges and the oracles",
          "[refinement][frozen][fe5]") {
    // The 15f loop after a frozen refine, as a piece runs it: the strip built
    // on the non-frozen constraints only; the DEM rescan (refine_strip) and
    // L12's on-edge insertion must not reach a frozen edge either.
    const std::uint32_t frozen = GENERATE(2u, 8u, 2u | 8u);
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth, Terrain::SmoothWithHole);
    const double tol = GENERATE(0.0, 0.5);
    const std::uint32_t seed = GENERATE(1u, 2u);
    CAPTURE(frozen, static_cast<int>(kind), tol, seed);
    const auto g = strip_oracle::exact_geometry();
    const auto dem = strip_oracle::terrain_dem(g, kind, seed);
    const auto start = strip_oracle::jittered(g, seed);
    const RefineOutcome phase1 = run(dem, start, Knobs{tol, frozen, true});
    REQUIRE(phase1.ok());
    const Begin b = begin_from(phase1);
    const auto strip = strip_of(dem, b, frozen);

    const auto out = run_strip(dem, strip, b, point_options(tol, frozen));
    REQUIRE(out.ok());
    check(dem, start, frozen, strip_oracle::mesh_of(out), tol);

    const CheckPoints store = store_of(g, random_points(dem, seed, 200));
    const auto pts = run_points(store, b, point_options(tol, frozen), &strip);
    REQUIRE(pts.ok());
    const Mesh m = strip_oracle::mesh_of(pts);
    REQUIRE(frozen_oracle::clean(
        frozen_oracle::frozen_findings(g, start.mesh.vertices(), start.edges, start.masks, frozen, m)));
    REQUIRE(strip_oracle::delaunay_violations(g, m) == 0);
}

TEST_CASE("FE5: with a strip, a check point within r(g) of a frozen edge is treated as on it",
          "[refinement][frozen][fe5][pinned]") {
    // Ruled by N7 (pinned here at the red step): L12 puts a point
    // within the coincidence radius r(g) of a constrained edge onto that edge
    // when a strip is given. On a frozen edge it must not, and inserting it by
    // split_inside instead leaves a vertex 1e-11 cells from the seam, a sliver
    // no flip can remove (F1, L14). So it is skipped and counted in on_frozen,
    // as a point exactly on the edge is. Without a strip (15c's path, radius
    // 0) nothing changes, so this case gives a strip.
    const auto start = grid_seam();
    const auto dem = ground(kSquare, {});
    Begin b = begin_from(dem, start);
    const auto strip = strip_of(dem, b, kSeam);
    const CheckPoints store = store_of(kSquare, {{Lat{4 + 1e-11, 2.5}, 100.0}});
    const auto out = run_points(store, b, point_options(0.5, kSeam), &strip);
    REQUIRE(out.ok());
    REQUIRE(out.on_frozen == 1);
    // N6: the error at its projection on the seam, sigma = 2.5 / 8 from (4, 0).
    const auto za = strip_oracle::height(dem, Lat{4, 0}), zb = strip_oracle::height(dem, Lat{4, 8});
    REQUIRE(std::abs(out.on_frozen_max_error - std::abs(100.0 - (*za + 2.5 / 8.0 * (*zb - *za)))) <= 1e-6);
    const Mesh m = strip_oracle::mesh_of(out);
    REQUIRE(frozen_oracle::clean(
        frozen_oracle::frozen_findings(kSquare, start.mesh.vertices(), start.edges, start.masks, kSeam, m)));
}

TEST_CASE("FE5: without a strip the frozen test is exact: a point 1e-11 off the seam goes in",
          "[refinement][frozen][fe5]") {
    // N7: "Without a strip, the frozen test is exact (radius 0), as L12 and L14
    // are; K1 and 15c's path are unchanged." The same point as the case above,
    // no strip: it is inside the right piece, not on the seam, so 15c inserts
    // it by split_inside and on_frozen stays 0.
    const auto start = grid_seam();
    const auto dem = ground(kSquare, {});
    const Begin b = begin_from(dem, start);
    const CheckPoints store = store_of(kSquare, {{Lat{4 + 1e-11, 2.5}, 100.0}});
    const auto out = run_points(store, b, point_options(0.5, kSeam));
    REQUIRE(out.ok());
    REQUIRE(out.on_frozen == 0);
    REQUIRE(out.inserted == 1);
}

// ------------------------------------------------------------------------ N16

namespace {

// A strip over EVERY constraint edge of the start, the frozen seam included:
// what a caller that forgot to leave frozen edges out would build.
struct WholeStrip {
    Start start = grid_seam();
    Raster<float> dem = ground(kSquare, {});
    Begin b = begin_from(dem, start);
    ConstraintCheckPoints strip = strip_of(dem, b, 0);
};

}  // namespace

TEST_CASE("N16: refine_strip refuses a strip edge that is frozen", "[refinement][frozen][n16]") {
    // 23-basin-scale.md N16: std::logic_error, its text starting with the entry
    // point's name, as 15f L2's other programming errors are. The same strip
    // with frozen_mask 0 runs (the control).
    const WholeStrip w;
    REQUIRE(run_strip(w.dem, w.strip, w.b, point_options(0.5, 0)).ok());
    REQUIRE_THROWS_MATCHES(run_strip(w.dem, w.strip, w.b, point_options(0.5, kSeam)), std::logic_error,
                           Catch::Matchers::MessageMatches(Catch::Matchers::StartsWith("refine_strip: ")));
}

TEST_CASE("N16: refine_points with a strip refuses a strip edge that is frozen", "[refinement][frozen][n16]") {
    const WholeStrip w;
    const CheckPoints store = store_of(kSquare, {{Lat{2, 5}, 30.0}});
    REQUIRE(run_points(store, w.b, point_options(0.5, 0), &w.strip).ok());
    REQUIRE_THROWS_MATCHES(run_points(store, w.b, point_options(0.5, kSeam), &w.strip), std::logic_error,
                           Catch::Matchers::MessageMatches(Catch::Matchers::StartsWith("refine_points: ")));
}

TEST_CASE("N16: the frozen-edge refusal comes after 15f's not-a-constraint-edge check",
          "[refinement][frozen][n16]") {
    // N16: "checked after L2's check (3) ... and before (4)". A strip with a
    // frozen edge AND an edge that is no constraint of the start is refused
    // for the second, under either entry point.
    const WholeStrip w;
    std::set<std::pair<std::uint32_t, std::uint32_t>> constrained;
    for (const auto& e : w.start.edges) constrained.insert(std::minmax(e[0], e[1]));
    Edges edges = w.start.edges;
    const auto& tri = w.start.mesh.triangles()[0];
    for (unsigned k = 0; k < 3 && edges.size() == w.start.edges.size(); ++k)
        if (!constrained.contains(std::minmax(tri[k], tri[(k + 1) % 3]))) edges.push_back({tri[k], tri[(k + 1) % 3]});
    REQUIRE(edges.size() == w.start.edges.size() + 1);
    const auto strip = constraint_check_points(w.dem, w.b.mesh.vertices(), EdgeSpan{edges});
    const auto not_constraint = Catch::Matchers::ContainsSubstring("is not a constraint edge of the start mesh");
    REQUIRE_THROWS_MATCHES(run_strip(w.dem, strip, w.b, point_options(0.5, kSeam)), std::logic_error,
                           Catch::Matchers::MessageMatches(not_constraint));
    const CheckPoints store = store_of(kSquare, {{Lat{2, 5}, 30.0}});
    REQUIRE_THROWS_MATCHES(run_points(store, w.b, point_options(0.5, kSeam), &strip), std::logic_error,
                           Catch::Matchers::MessageMatches(not_constraint));
}

TEST_CASE("FE5: with a strip, on_frozen survives a strip point winning the triangle", "[refinement][frozen][fe5]") {
    // One stored point exactly on the seam, a strip on every other constraint
    // edge, and a tolerance nothing exceeds, so the run is one scan. The
    // triangle that counts the seam point also owns an outline sub-edge whose
    // strip points have nonzero error, and the worst of them is offered over
    // the (skipped) stored point. The count must survive that offer.
    const auto start = grid_seam();
    const auto dem = ground(kSquare, {});
    const Begin b = begin_from(dem, start);
    const auto strip = strip_of(dem, b, kSeam);
    REQUIRE(strip.size() > 0);
    const CheckPoints store = store_of(kSquare, {{Lat{4, 2.5}, 100.0}});
    const auto out = run_points(store, b, point_options(1e9, kSeam), &strip);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 0);
    REQUIRE(out.strip_max_error > 0.0);  // the strip had something to offer
    REQUIRE(out.on_frozen == 1);
}
