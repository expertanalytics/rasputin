// Increment 34 (docs/increments/34-slope-tolerance.md, sections 3, 4.3, 4.5
// and 9): terrain::refinement::SlopeRamp and SlopeTolerance on their own.
//
//   test 3      the ramp: t at 0, START, the middle, END and 89.5; the step;
//               allowed() never increasing over classes 0-180; entries
//               181-255 equal F; weight(c) = 1 / allowed(c).
//   test 6 (b)  cell_class on V1n's cell (31, 30), which has the NoData node
//               as a corner: 86, the largest VALID corner, not 255 (M9); 0
//               for a cell with no valid corner.
//   test 9      make's refusals, one sentence each naming the bound ("near",
//               "far", "start", "end"), and the accepted edges.
//   histogram   entry 255 counts NoData; the entries sum to the node count;
//               row() is steepness()'s row.
//
// Interface used, as section 4.3 writes it: SlopeRamp{near, far, start_deg,
// end_deg}.at(s); SlopeTolerance::make(dem, ramp, threads, why) ->
// std::optional; geometry(), row(r), cell_class(mesh::MeshVertex{col, row}),
// allowed(c), weight(c), histogram(), near(), far(). Starts threads only in
// make (one thread here). Of section 9's invariant-critical suite (tests 2,
// 4, 5, 6 and 8), this file holds test 6 (b)'s cell_class case.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/steepness.hpp>
#include <terrain/refinement/slope_tolerance.hpp>

#include "slope_oracle.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

using namespace slope_oracle;
using terrain::mesh::MeshVertex;
using terrain::refinement::SlopeRamp;
using terrain::refinement::SlopeTolerance;

namespace {

SlopeTolerance built(const Raster<float>& dem, SlopeRamp r) {
    std::string why;
    auto s = SlopeTolerance::make(dem, r, 1, why);
    INFO(why);
    REQUIRE(s.has_value());
    return std::move(*s);
}

// The refusal's sentence, which must name `word`.
std::string refused(SlopeRamp r) {
    std::string why;
    const auto s = SlopeTolerance::make(v1(), r, 1, why);
    CHECK_FALSE(s.has_value());
    CHECK_FALSE(why.empty());
    return why;
}

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInf = std::numeric_limits<double>::infinity();

}  // namespace

// ------------------------------------------------------------------- test 3

TEST_CASE("34 test 3: the ramp at 0, START, the middle, END and 89.5 degrees", "[slope][ramp][test3]") {
    const SlopeRamp r{2.0, 10.0, 25.0, 35.0};
    CHECK(r.at(0.0) == 10.0);
    CHECK(r.at(25.0) == 10.0);
    CHECK(r.at(30.0) == 6.0);
    CHECK(r.at(35.0) == 2.0);
    CHECK(r.at(89.5) == 2.0);
    const Ramp oracle{2.0, 10.0, 25.0, 35.0};
    for (int c = 0; c <= 180; ++c) {
        CAPTURE(c);
        CHECK(std::abs(r.at(c / 2.0) - oracle.at(c / 2.0)) <= 1e-12 * 10.0);
    }
}

TEST_CASE("34 test 3: a step (START = END) gives N from END up, F below", "[slope][ramp][test3]") {
    const SlopeRamp step{2.0, 10.0, 30.0, 30.0};
    CHECK(step.at(0.0) == 10.0);
    CHECK(step.at(29.5) == 10.0);
    CHECK(step.at(30.0) == 2.0);
    CHECK(step.at(30.5) == 2.0);
    const SlopeRamp every{2.0, 10.0, 0.0, 0.0};  // START = END = 0: every node held to N (G6)
    CHECK(every.at(0.0) == 2.0);
    CHECK(every.at(45.0) == 2.0);
}

TEST_CASE("34 test 3: the tables: allowed(c) is the ramp at c / 2, never increasing, F from 181 to 255",
          "[slope][ramp][test3]") {
    const auto ramp = GENERATE(SlopeRamp{2.0, 10.0, 25.0, 35.0}, SlopeRamp{2.0, 10.0, 30.0, 30.0},
                               SlopeRamp{0.5, 0.5, 10.0, 60.0}, SlopeRamp{1.0, 7.0, 0.0, 0.0});
    CAPTURE(ramp.near, ramp.far, ramp.start_deg, ramp.end_deg);
    const auto s = built(v1(), ramp);
    CHECK(s.near() == ramp.near);
    CHECK(s.far() == ramp.far);
    double before = kInf;
    for (int c = 0; c <= 180; ++c) {
        CAPTURE(c);
        const auto k = static_cast<std::uint8_t>(c);
        CHECK(s.allowed(k) == ramp.at(c / 2.0));
        CHECK(s.allowed(k) <= before);
        CHECK(s.weight(k) == 1.0 / s.allowed(k));
        before = s.allowed(k);
    }
    for (int c = 181; c <= 255; ++c) {
        CAPTURE(c);
        const auto k = static_cast<std::uint8_t>(c);
        CHECK(s.allowed(k) == ramp.far);
        CHECK(s.weight(k) == 1.0 / ramp.far);
    }
}

// ------------------------------------------------------------- test 6 (b)

TEST_CASE("34 test 6 (b): cell_class beside NoData is the largest valid corner, not 255 (M9)",
          "[slope][cell_class][test6]") {
    // V1n's four cells around the NoData node (row 32, column 30): largest
    // valid corner classes 83, 86, 83 and 81 (README item 4b).
    const auto s = built(v1n(), SlopeRamp{2.0, 10.0, 30.0, 30.0});
    CHECK(static_cast<int>(s.cell_class(MeshVertex{30.5, 31.5})) == 86);
    CHECK(static_cast<int>(s.cell_class(MeshVertex{29.25, 31.75})) == 83);
    CHECK(static_cast<int>(s.cell_class(MeshVertex{29.5, 32.5})) == 83);
    CHECK(static_cast<int>(s.cell_class(MeshVertex{30.125, 32.875})) == 81);
    // A node is filed in the cell it opens, (floor(row), floor(col)), and is
    // still judged by that cell's corners (section 4.3).
    CHECK(static_cast<int>(s.cell_class(MeshVertex{30.0, 31.0})) == 86);
    // The same against the oracle's classes, every cell of V1n.
    const auto cls = classes(v1n());
    const auto& g = s.geometry();
    for (std::size_t r = 0; r + 1 < g.rows(); ++r)
        for (std::size_t c = 0; c + 1 < g.cols(); ++c) {
            const double col = static_cast<double>(c) + 0.5, row = static_cast<double>(r) + 0.5;
            const auto want = cell_class(g, cls, col, row);
            const auto got = s.cell_class(MeshVertex{col, row});
            if (got != want) {
                CAPTURE(r, c, static_cast<int>(want), static_cast<int>(got));
                CHECK(got >= want);  // the producer's class may sit one above at a boundary
                CHECK(got <= want + 1);
            }
        }
}

TEST_CASE("34 test 6 (b): cell_class clamps to the last cell, and is 0 with no valid corner",
          "[slope][cell_class][test6]") {
    std::vector<float> z(5 * 5);
    for (std::size_t i = 0; i < z.size(); ++i) z[i] = static_cast<float>(4 * (i % 5));  // 21.8 degrees in x
    for (const std::size_t i : {0u, 1u, 5u, 6u}) z[i] = kNoData;  // cell (0, 0) has no valid corner
    const Raster<float> dem{geometry(5, 5, 10.0, 10.0), std::move(z), kNoData};
    const auto s = built(dem, SlopeRamp{2.0, 10.0, 30.0, 30.0});
    CHECK(static_cast<int>(s.cell_class(MeshVertex{0.5, 0.5})) == 0);
    const auto cls = classes(dem);
    // The far corner node and a point on the far edges file in the last cell, (3, 3).
    CHECK(s.cell_class(MeshVertex{4.0, 4.0}) == cell_class(dem.geometry(), cls, 3.5, 3.5));
    CHECK(s.cell_class(MeshVertex{4.0, 3.25}) == cell_class(dem.geometry(), cls, 3.5, 3.25));
    CHECK(static_cast<int>(cell_class(dem.geometry(), cls, 3.5, 3.5)) == 44);  // atan(0.4) = 21.8 degrees
}

// ------------------------------------------------------------------- test 9

TEST_CASE("34 test 9: make refuses a bad ramp with one sentence naming the bound", "[slope][make][test9]") {
    CHECK(refused(SlopeRamp{0.0, 10.0, 25.0, 35.0}).find("near") != std::string::npos);    // N > 0
    CHECK(refused(SlopeRamp{-1.0, 10.0, 25.0, 35.0}).find("near") != std::string::npos);
    CHECK(refused(SlopeRamp{12.0, 10.0, 25.0, 35.0}).find("near") != std::string::npos);   // N <= F
    CHECK(refused(SlopeRamp{kNaN, 10.0, 25.0, 35.0}).find("near") != std::string::npos);
    CHECK(refused(SlopeRamp{2.0, kNaN, 25.0, 35.0}).find("far") != std::string::npos);
    CHECK(refused(SlopeRamp{2.0, kInf, 25.0, 35.0}).find("far") != std::string::npos);
    CHECK(refused(SlopeRamp{2.0, 10.0, -1.0, 35.0}).find("start") != std::string::npos);   // START >= 0
    CHECK(refused(SlopeRamp{2.0, 10.0, 36.0, 35.0}).find("start") != std::string::npos);   // START <= END
    CHECK(refused(SlopeRamp{2.0, 10.0, kNaN, 35.0}).find("start") != std::string::npos);
    CHECK(refused(SlopeRamp{2.0, 10.0, 25.0, 90.0}).find("end") != std::string::npos);     // END < 90
    CHECK(refused(SlopeRamp{2.0, 10.0, 25.0, kInf}).find("end") != std::string::npos);
}

TEST_CASE("34 test 9: make accepts the edges of the bounds", "[slope][make][test9]") {
    for (const SlopeRamp r : {SlopeRamp{10.0, 10.0, 25.0, 35.0}, SlopeRamp{2.0, 10.0, 0.0, 0.0},
                              SlopeRamp{2.0, 10.0, 89.5, 89.5}, SlopeRamp{1e-9, 10.0, 0.0, 89.999}}) {
        CAPTURE(r.near, r.far, r.start_deg, r.end_deg);
        std::string why;
        CHECK(SlopeTolerance::make(v1(), r, 1, why).has_value());
        CHECK(why.empty());
    }
}

// ---------------------------------------------------------------- histogram

TEST_CASE("34: histogram() counts every node once, NoData in entry 255, as row() holds them",
          "[slope][histogram]") {
    const bool hole = GENERATE(false, true);
    const auto dem = hole ? v1n() : v1();
    const auto s = built(dem, SlopeRamp{2.0, 10.0, 30.0, 30.0});
    const auto h = s.histogram();
    std::size_t total = 0;
    for (const std::size_t n : h) total += n;
    CHECK(total == 65u * 65u);
    CHECK(h[255] == (hole ? 1u : 0u));
    for (std::size_t c = 181; c < 255; ++c) CHECK(h[c] == 0u);
    const auto all = terrain::raster::steepness(dem, 1);
    std::array<std::size_t, 256> mine{};
    for (std::size_t r = 0; r < 65; ++r) {
        const auto row = s.row(r);
        REQUIRE(row.size() == 65u);
        for (std::size_t c = 0; c < 65; ++c) {
            CHECK(row[c] == all[r * 65 + c]);
            ++mine[row[c]];
        }
    }
    CHECK(mine == h);
    // Section 9's figure: with the step at 30, 1 496 of V1's 4 225 nodes are
    // held to 2 m (1 495 of V1n's 4 224 valid ones).
    std::size_t steep = 0;
    for (std::size_t c = 60; c <= 180; ++c) steep += h[c];
    CHECK(steep == (hole ? 1495u : 1496u));
}
