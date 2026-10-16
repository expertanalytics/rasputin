// Increment 20c, PR 20c-2 (docs/increments/20c-soft-quality.md, R7, R8 and
// "Tests @tester writes red first", 20c-2): the soft criterion and the line
// split through refine. T-P3 (gain -1 is 20c-1's output), determinism as
// 14b's T6, the QA rules' section D oracles with the gain on, and the
// counters refine reports. Not mutation-critical by the design; the
// invariant-critical suite is tests/cpp/unit/test_quality_gain.cpp.
//
// Interface, as R7 and R8 name it ("RefineOptions and the binding carry it"):
//
//   RefineOptions::min_gain_deg   passed to QualityOptions::min_gain_deg
//
// PINNED HERE, where the design leaves it open (listed for @architect):
//   - RefineOptions::min_gain_deg defaults to a negative value (off), as
//     QualityOptions' does; the CLI passes 0 (R7);
//   - RefineOutcome::quality_no_gain (QualityOutcome::skipped_no_gain) and
//     RefineOutcome::quality_line_splits (QualityOutcome::line_splits);
//     quality_skipped keeps its one total and holds quality_no_gain (R7,
//     "start_quality_points_skipped keeps its total");
//   - a line split is a vertex of its own kind: in none of quality_inserted,
//     quality_feet, inserted, so the output holds start + quality_inserted +
//     quality_feet + quality_line_splits + inserted vertices;
//   - the exact counts on CF2's own-edge start (one refusal, nothing else
//     skipped), from running a prototype of R7 (see the handback).
//
// T-P3's digests were RECORDED FROM 41bda81a (master with 20c-1 merged, no
// 20c-2 production change) with refine_digest::topology_digest, on
// tests/cpp/support/quality_gain_fixtures.hpp's scenes and options, by a scratch
// program that included that header. No commit may update them to agree
// with new code.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>

#include "constraint_foot_oracles.hpp"
#include "quality_fixtures.hpp"
#include "quality_gain_fixtures.hpp"
#include "refine_digest.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

using terrain::Point2;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::RefineOptions;
using terrain::refinement::RefineOutcome;

namespace cfo = constraint_foot_oracles;
namespace qgf = quality_gain_fixtures;

namespace {

// 20c-1's topology digests of qgf::scene(i) under qgf::options, RECORDED (see
// the header). Integers only, so they hold across compilers and libm.
constexpr std::array<std::uint64_t, qgf::kScenes> kRecorded{
    0x1253792f2761a97dull,  // the features start (equal to 20c-1's own feet-off record: no foot there)
    0x228cb3c7622c28c1ull,  // the Q7 ring
    0x876c079d4439001aull,  // the line fixture
};
static_assert(qgf::kScenes == 3, "the GENERATE lists below name every scene");

RefineOptions with_gain(const qgf::Scene& s, double gain, bool feet = true, unsigned threads = 1) {
    RefineOptions o = qgf::options(s, threads);
    o.constraint_feet = feet;
    o.min_gain_deg = gain;
    return o;
}

// The QA rules' section D pair and the lines, from the output alone
// (tests/cpp/support/constraint_foot_oracles.hpp), and the vertex identity (pinned).
void section_d(const Raster<float>& dem, const quality_fixtures::Start& s, const RefineOutcome& out, double tol) {
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
    CHECK(out.vertices.size()
          == s.mesh.vertices().size() + out.quality_inserted + out.quality_feet + out.quality_line_splits
                 + out.inserted);
    CHECK(out.quality_skipped >= out.quality_no_gain);
}

}  // namespace

TEST_CASE("R7: RefineOptions leaves the gain test off by default", "[refinement][quality][gain][R7]") {
    REQUIRE(RefineOptions{}.min_gain_deg < 0.0);  // pinned
}

// ------------------------------------------------------------------ T-P3

TEST_CASE("T-P3: gain -1 gives 20c-1's output, recorded before the change", "[refinement][quality][gain][T-P3]") {
    const auto i = GENERATE(std::size_t{0}, std::size_t{1}, std::size_t{2});
    CAPTURE(i);
    const auto s = qgf::scene(i);
    const auto off = qgf::run(s, with_gain(s, -1.0));
    REQUIRE(off.ok());
    REQUIRE(refine_digest::topology_digest(off) == kRecorded[i]);  // RECORDED, see the header
    REQUIRE(off.quality_no_gain == 0);
    REQUIRE(off.quality_line_splits == 0);
    // The default is the same switch, bit for bit, doubles included.
    REQUIRE(refine_digest::digest(qgf::run(s, qgf::options(s))) == refine_digest::digest(off));
}

TEST_CASE("T-P3: gain 0 changes the scenes that have work for it", "[refinement][quality][gain][T-P3]") {
    // Not vacuous: on the features start R7 refuses, on the line fixture R7
    // refuses and R8 splits, and both outputs differ from gain -1's.
    const auto features = qgf::scene(0);
    const auto f0 = qgf::run(features, with_gain(features, 0.0));
    REQUIRE(f0.ok());
    CHECK(f0.quality_no_gain > 0);
    CHECK(refine_digest::topology_digest(f0) != kRecorded[0]);

    const auto line = qgf::scene(2);
    const auto l0 = qgf::run(line, with_gain(line, 0.0));
    REQUIRE(l0.ok());
    CHECK(l0.quality_no_gain > 0);
    CHECK(l0.quality_line_splits > 0);
    CHECK(refine_digest::topology_digest(l0) != kRecorded[2]);
}

// ------------------------------------------------------------------ section D, the gain on

TEST_CASE("Section D: with the gain on, the oracles hold on every scene, feet on and off, at any tolerance",
          "[refinement][quality][gain][T3]") {
    const auto i = GENERATE(std::size_t{0}, std::size_t{1}, std::size_t{2});
    const double gain = GENERATE(0.0, 2.0);
    const bool feet = GENERATE(true, false);
    const double tol = GENERATE(0.0, 0.5, 3.0);
    CAPTURE(i, gain, feet, tol);
    const auto s = qgf::scene(i);
    RefineOptions o = with_gain(s, gain, feet);
    o.tolerance = tol;
    section_d(s.dem, s.start, qgf::run(s, o), tol);
}

// ------------------------------------------------------------------ the counters

namespace {

// 20c-1's CF2 own-edge start through refine (tests/cpp/unit/test_constraint_foot_quality.cpp):
// A (0, 0.99), C (2, 6.99), B (20, 0.99) on a 21 x 8 grid with 1 m cells, a
// flat DEM at 3 m. Its node N = (10, 1) lies 0.01 cells inside A-B, so with
// feet off the node's triangles hold a sliver: R7 refuses it.
RasterGeometry own_geometry() { return RasterGeometry{0.0, 7.0, 1.0, 1.0, 21, 8}; }

quality_fixtures::Start own_start() {
    const auto w = [](double c, double r) { return Point2{c, 7.0 - r}; };
    quality_fixtures::Start s;
    s.mesh = terrain::IndexedMesh2{{w(0, 0.99), w(2, 6.99), w(20, 0.99)}, {{0, 1, 2}}, std::vector<std::uint8_t>(1, 0)};
    s.edges = {{0, 2}, {0, 1}, {1, 2}};
    s.masks = {1, 2, 4};
    return s;
}

RefineOutcome run_own(double gain) {
    const Raster<float> dem{own_geometry(), std::vector<float>(21 * 8, 3.0f)};
    const auto s = own_start();
    RefineOptions o;
    o.tolerance = 1.0;
    o.threads = 1;
    o.min_angle_deg = 25.0;
    o.constraint_feet = false;
    o.min_gain_deg = gain;
    return terrain::refinement::refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                                       std::span<const std::uint32_t>{s.masks}, o);
}

}  // namespace

TEST_CASE("R7 through refine: a refusal is counted in quality_no_gain and in the quality_skipped total",
          "[refinement][quality][gain][R7]") {
    const auto hard = run_own(-1.0);
    REQUIRE(hard.ok());
    REQUIRE(hard.quality_inserted == 1);  // the premise: the hard rule takes N
    REQUIRE(hard.quality_no_gain == 0);
    const auto soft = run_own(0.0);
    REQUIRE(soft.ok());
    REQUIRE(soft.quality_no_gain == 1);  // pinned, from the prototype run
    REQUIRE(soft.quality_inserted == 0);
    REQUIRE(soft.quality_line_splits == 0);
    REQUIRE(soft.quality_skipped == 1);  // the one total holds the refusal
    REQUIRE(soft.vertices.size() == 3 + soft.inserted);
}

// ------------------------------------------------------------------ T6, the gain on

TEST_CASE("T6: with the gain on, the output is bit-identical for 1, 2, 7 and all threads",
          "[refinement][quality][gain][T6]") {
    const auto s = qgf::scene(2);  // the scene where R7 refuses and R8 splits
    const auto ref = qgf::run(s, with_gain(s, 0.0, true, 1));
    REQUIRE(ref.ok());
    REQUIRE(ref.quality_line_splits > 0);
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    const auto other = qgf::run(s, with_gain(s, 0.0, true, threads));
    REQUIRE(refine_digest::digest(other) == refine_digest::digest(ref));
    REQUIRE(other.quality_no_gain == ref.quality_no_gain);
    REQUIRE(other.quality_line_splits == ref.quality_line_splits);
    REQUIRE(other.quality_inserted == ref.quality_inserted);
}
