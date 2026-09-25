// Increment 14 (docs/increments/14-adaptive-refinement.md, R2 and R6): the
// per-triangle scan. INVARIANT-CRITICAL: the tolerance guarantee is decided
// here, so @reviewer mutation-tests this suite.
//
// Interface: include/terrain/refinement/scan.hpp (LatticeMesh from
// include/terrain/mesh/lattice_mesh.hpp). Edge k runs from vertex k to vertex k+1.
//
// Every DEM below is float and every expected error is a dyadic rational, so
// the equalities are exact: z = 3*col - 2*row + 7 plus a deviation in {0.25,
// 0.5, 1}. The base triangle is (0,0) (8,0) (8,8) in (row, col); its closed
// node set is { 0 <= col <= row <= 8 }, 45 nodes.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/scan.hpp>

#include "refinement_fixtures.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::raster::Raster;
using terrain::refinement::NodeLocation;
using terrain::refinement::scan;

namespace {

constexpr std::size_t N = 9;
constexpr float kSentinel = -32767.0f;
const float kNaN = std::numeric_limits<float>::quiet_NaN();

struct Dem {
    std::vector<float> z;
    std::optional<float> nodata;
    Dem() : z(N * N) {
        for (std::size_t r = 0; r < N; ++r)
            for (std::size_t c = 0; c < N; ++c)
                z[r * N + c] = static_cast<float>(3.0 * static_cast<double>(c)
                                                  - 2.0 * static_cast<double>(r) + 7.0);
    }
    Dem& bump(std::size_t r, std::size_t c, float dz) {
        z[r * N + c] += dz;
        return *this;
    }
    Dem& set(std::size_t r, std::size_t c, float v) {
        z[r * N + c] = v;
        return *this;
    }
    [[nodiscard]] Raster<float> raster() const {
        return Raster<float>{refinement_fixtures::geometry(N, N), z, nodata};
    }
};

// The base triangle, its vertex order rotated by `rot`: triangle edge j runs
// from V[(rot + j) % 3] to V[(rot + j + 1) % 3].
LatticeMesh base(unsigned rot = 0) {
    const std::array<LatticeVertex, 3> v{LatticeVertex{0, 0}, LatticeVertex{8, 0}, LatticeVertex{8, 8}};
    auto m = LatticeMesh::build({v[0], v[1], v[2]},
                                {TriangleIndices{rot % 3, (rot + 1) % 3, (rot + 2) % 3}}, {0},
                                {std::array<std::uint32_t, 3>{}});
    REQUIRE(m.has_value());
    return *m;
}

NodeLocation edge(unsigned j) {
    constexpr std::array<NodeLocation, 3> e{NodeLocation::Edge0, NodeLocation::Edge1,
                                            NodeLocation::Edge2};
    return e[j];
}

}  // namespace

TEST_CASE("scan: a plane DEM has error exactly zero", "[refinement][scan]") {
    const auto dem = Dem{}.raster();
    REQUIRE(scan(dem, base(), 0).max_error == 0.0);
    REQUIRE_FALSE(scan(dem, base(), 0).is_void);
}

TEST_CASE("scan: the error is |z - plane| at a strictly interior node", "[refinement][scan]") {
    const auto dem = Dem{}.bump(5, 2, 0.25f).raster();
    const auto r = scan(dem, base(), 0);
    REQUIRE(r.max_error == 0.25);
    REQUIRE(r.node == LatticeVertex{5, 2});
    REQUIRE(r.where == NodeLocation::Inside);
}

TEST_CASE("scan: a node on an edge is in the set, and says which edge", "[refinement][scan]") {
    // Original edge k of V: 0 is (0,0)-(8,0), 1 is (8,0)-(8,8), 2 is (8,8)-(0,0).
    const unsigned k = GENERATE(0u, 1u, 2u);
    const unsigned rot = GENERATE(0u, 1u, 2u);
    CAPTURE(k, rot);
    constexpr std::array<std::array<std::size_t, 2>, 3> on{{{4, 0}, {8, 3}, {5, 5}}};
    const auto dem = Dem{}.bump(on[k][0], on[k][1], -0.5f).raster();
    const auto r = scan(dem, base(rot), 0);
    REQUIRE(r.max_error == 0.5);
    REQUIRE(r.node == LatticeVertex{static_cast<std::uint32_t>(on[k][0]),
                                    static_cast<std::uint32_t>(on[k][1])});
    REQUIRE(r.where == edge((k + 3 - rot) % 3));
}

TEST_CASE("scan: the three vertices are not in the set", "[refinement][scan]") {
    // A triangle whose closed node set is exactly its vertices, on a DEM far
    // off any plane: error 0 and no node.
    Dem d;
    d.set(0, 0, 0.0f).set(1, 0, 100.0f).set(1, 1, -50.0f);
    const auto dem = d.raster();
    auto m = LatticeMesh::build({LatticeVertex{0, 0}, LatticeVertex{1, 0}, LatticeVertex{1, 1}},
                                {TriangleIndices{0, 1, 2}}, {0}, {std::array<std::uint32_t, 3>{}});
    REQUIRE(m.has_value());
    const auto r = scan(dem, *m, 0);
    REQUIRE(r.max_error == 0.0);
    REQUIRE_FALSE(r.node.has_value());
}

TEST_CASE("scan: a node in the bounding box but outside the triangle is not in the set",
          "[refinement][scan]") {
    const auto dem = Dem{}.bump(2, 5, 1000.0f).bump(7, 8, 1000.0f).raster();  // col > row
    REQUIRE(scan(dem, base(), 0).max_error == 0.0);
}

TEST_CASE("scan: NoData nodes are skipped, NaN and sentinel alike", "[refinement][scan]") {
    SECTION("NaN") {
        const auto dem = Dem{}.set(5, 2, kNaN).bump(6, 3, 0.5f).raster();
        const auto r = scan(dem, base(), 0);
        REQUIRE(r.max_error == 0.5);
        REQUIRE(r.node == LatticeVertex{6, 3});
    }
    SECTION("sentinel") {
        Dem d;
        d.set(5, 2, kSentinel).bump(6, 3, 0.5f).nodata = kSentinel;
        const auto r = scan(d.raster(), base(), 0);
        REQUIRE(r.max_error == 0.5);
        REQUIRE(r.node == LatticeVertex{6, 3});
        REQUIRE_FALSE(r.is_void);
    }
    SECTION("every non-vertex node NoData: error 0, no node") {
        Dem d;
        for (std::size_t r = 0; r < N; ++r)
            for (std::size_t c = 0; c <= r; ++c)
                if (!((r == 0 && c == 0) || (r == 8 && c == 0) || (r == 8 && c == 8)))
                    d.set(r, c, kNaN);
        const auto r = scan(d.raster(), base(), 0);
        REQUIRE(r.max_error == 0.0);
        REQUIRE_FALSE(r.node.has_value());
    }
}

TEST_CASE("scan: ties go to the smallest (row, col)", "[refinement][scan]") {
    const unsigned rot = GENERATE(0u, 1u, 2u);
    CAPTURE(rot);
    SECTION("row decides before column") {
        const auto dem = Dem{}.bump(6, 1, 1.0f).bump(5, 2, 1.0f).raster();
        REQUIRE(scan(dem, base(rot), 0).node == LatticeVertex{5, 2});
    }
    SECTION("then column") {
        const auto dem = Dem{}.bump(6, 4, 1.0f).bump(6, 1, 1.0f).raster();
        REQUIRE(scan(dem, base(rot), 0).node == LatticeVertex{6, 1});
    }
    SECTION("the sign of the deviation does not matter") {
        const auto dem = Dem{}.bump(7, 3, -1.0f).bump(4, 1, 1.0f).raster();
        const auto r = scan(dem, base(rot), 0);
        REQUIRE(r.max_error == 1.0);
        REQUIRE(r.node == LatticeVertex{4, 1});
    }
    SECTION("an edge node against an interior node") {
        const auto dem = Dem{}.bump(6, 3, 1.0f).bump(4, 0, 1.0f).raster();
        const auto r = scan(dem, base(rot), 0);
        REQUIRE(r.node == LatticeVertex{4, 0});
        REQUIRE(r.where == edge((3 - rot) % 3));
    }
    SECTION("a strictly larger error beats a smaller index") {
        const auto dem = Dem{}.bump(4, 1, 0.5f).bump(7, 6, 1.0f).raster();
        REQUIRE(scan(dem, base(rot), 0).node == LatticeVertex{7, 6});
    }
}

TEST_CASE("scan: a void triangle is carved at the valid node nearest a NoData vertex",
          "[refinement][scan]") {
    SECTION("one NoData vertex") {
        const auto dem = Dem{}.set(0, 0, kNaN).raster();
        const auto r = scan(dem, base(), 0);
        REQUIRE(r.is_void);
        REQUIRE(r.node == LatticeVertex{1, 0});  // squared distance 1 from (0,0)
        REQUIRE(r.where == NodeLocation::Edge0);
        REQUIRE(r.uncovered == 42);  // 45 closed nodes minus 3 vertices, all valid
    }
    SECTION("two NoData vertices, equally near candidates: smallest (row, col)") {
        // (1,0) is 1 from (0,0); (8,7) is 1 from (8,8).
        const auto dem = Dem{}.set(0, 0, kNaN).set(8, 8, kNaN).raster();
        const auto r = scan(dem, base(), 0);
        REQUIRE(r.is_void);
        REQUIRE(r.node == LatticeVertex{1, 0});
    }
    SECTION("the nearest node is NoData itself: the next nearest valid one") {
        // (1,0) is gone; at squared distance 2 only (1,1) remains, on edge 2.
        const auto dem = Dem{}.set(0, 0, kNaN).set(1, 0, kNaN).raster();
        const auto r = scan(dem, base(), 0);
        REQUIRE(r.node == LatticeVertex{1, 1});
        REQUIRE(r.where == NodeLocation::Edge2);
        REQUIRE(r.uncovered == 41);
    }
    SECTION("a void triangle holding no valid node is done") {
        auto m = LatticeMesh::build({LatticeVertex{0, 0}, LatticeVertex{1, 0}, LatticeVertex{1, 1}},
                                    {TriangleIndices{0, 1, 2}}, {0},
                                    {std::array<std::uint32_t, 3>{}});
        REQUIRE(m.has_value());
        const auto r = scan(Dem{}.set(1, 0, kNaN).raster(), *m, 0);
        REQUIRE(r.is_void);
        REQUIRE_FALSE(r.node.has_value());
        REQUIRE(r.uncovered == 0);
    }
}
