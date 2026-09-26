// Increment 17 (docs/increments/17-mesh-stats.md, R6 and C2): refine's phase
// split and carving counter on RefineOutcome. Not invariant-critical, no
// mutation round: nothing here changes a mesh.
//
// The design names this suite "test_refinement_refine, extended"; no such
// binary existed (the refine tests live in prop_refinement_refine), so this is
// that name as a new, small unit suite. It is in the TSan job's lists because
// refine starts threads through for_each_chunk.
//
// Timings are never compared with each other or across thread counts (R6: not
// part of the determinism guarantee). They are checked only as non-negative
// and bounded by a steady_clock interval around the call. Ola chose C2 (a):
// `carved` counts splits of void triangles.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>

#include "refinement_fixtures.hpp"

#include <array>
#include <chrono>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <vector>

using terrain::raster::Raster;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;
using terrain::refinement::RefineOutcome;
using namespace refinement_fixtures;

namespace {

struct Timed {
    RefineOutcome out;
    double wall_seconds;
};

Timed timed_run(const Raster<float>& dem, const StartMesh& s, double tol, unsigned threads) {
    const auto t0 = std::chrono::steady_clock::now();
    auto out = refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                      std::span<const std::uint32_t>{s.masks},
                      RefineOptions{.tolerance = tol, .threads = threads});
    const auto t1 = std::chrono::steady_clock::now();
    return {std::move(out), std::chrono::duration<double>(t1 - t0).count()};
}

void check_seconds(const Timed& t) {
    const RefineOutcome& out = t.out;
    REQUIRE(out.legalise_seconds >= 0.0);
    REQUIRE(out.scan_seconds >= 0.0);
    REQUIRE(out.split_seconds >= 0.0);
    REQUIRE(out.legalise_seconds + out.scan_seconds + out.split_seconds <= t.wall_seconds);
}

Raster<float> rough(std::size_t n) { return Raster<float>{geometry(n, n), rough_dem(n, n, 7)}; }

Raster<float> nodata_corner(std::size_t n, bool sentinel) {
    auto z = rough_dem(n, n, 3);
    const float hole = sentinel ? -32767.0f : std::numeric_limits<float>::quiet_NaN();
    for (std::size_t r = 0; r <= 4; ++r)
        for (std::size_t c = 0; c <= 6; ++c) z[r * n + c] = hole;
    return Raster<float>{geometry(n, n), std::move(z),
                         sentinel ? std::optional<float>{-32767.0f} : std::nullopt};
}

}  // namespace

TEST_CASE("R6: the three phase seconds are non-negative and within the call's wall time",
          "[refinement][refine][stats]") {
    const unsigned threads = GENERATE(1u, 4u);
    CAPTURE(threads);
    const auto dem = rough(33);
    const auto t = timed_run(dem, grid_mesh(dem.geometry(), 8), 1.0, threads);
    REQUIRE(t.out.ok());
    REQUIRE(t.out.rounds > 1);  // the scan and split intervals were actually entered
    check_seconds(t);
}

TEST_CASE("R6: a run that needs no split still reports sane seconds", "[refinement][refine][stats]") {
    const std::size_t n = 9;
    const Raster<float> dem{geometry(n, n), std::vector<float>(n * n, 5.0f)};
    const auto t = timed_run(dem, grid_mesh(dem.geometry(), 4), 0.0, 1);
    REQUIRE(t.out.ok());
    REQUIRE(t.out.inserted == 0);
    check_seconds(t);
    REQUIRE(t.out.carved == 0);
}

TEST_CASE("R6: a refusal reports zero seconds and zero carved", "[refinement][refine][stats]") {
    const auto dem = rough(9);
    const auto t = timed_run(dem, grid_mesh(dem.geometry(), 4), -1.0, 1);  // negative tolerance
    REQUIRE_FALSE(t.out.ok());
    REQUIRE(t.out.legalise_seconds == 0.0);
    REQUIRE(t.out.scan_seconds == 0.0);
    REQUIRE(t.out.split_seconds == 0.0);
    REQUIRE(t.out.carved == 0);
}

TEST_CASE("C2: carved is 0 on a DEM without NoData", "[refinement][refine][stats]") {
    const auto dem = rough(33);
    const auto t = timed_run(dem, grid_mesh(dem.geometry(), 8), 1.0, 1);
    REQUIRE(t.out.ok());
    REQUIRE(t.out.inserted > 0);
    REQUIRE(t.out.carved == 0);
}

TEST_CASE("C2: carved is positive with a NoData corner, and never above inserted",
          "[refinement][refine][stats]") {
    const bool sentinel = GENERATE(false, true);
    CAPTURE(sentinel);
    const auto dem = nodata_corner(17, sentinel);
    const auto t = timed_run(dem, grid_mesh(dem.geometry(), 8), 2.0, 1);
    REQUIRE(t.out.ok());
    REQUIRE(t.out.carved > 0);
    REQUIRE(t.out.carved <= t.out.inserted);
    check_seconds(t);
}

TEST_CASE("C2: carved is the same for 1 and 4 threads", "[refinement][refine][stats]") {
    // Counters are part of the output, unlike the seconds (R6).
    const auto dem = nodata_corner(17, false);
    const auto start = grid_mesh(dem.geometry(), 8);
    const auto one = timed_run(dem, start, 2.0, 1);
    const auto four = timed_run(dem, start, 2.0, 4);
    REQUIRE(one.out.carved == four.out.carved);
}
