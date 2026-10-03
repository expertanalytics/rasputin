// Increment 23b (docs/increments/23-basin-scale.md, "The seam protocol" step 4,
// first bullet, and FE2): the scan does not count a node lying exactly on a
// frozen edge of the triangle; only triangles that have a frozen edge take
// that path. Part of FE2, which is invariant-critical.
//
// The scan reads the frozen mask from the mesh (LatticeMesh::set_frozen_mask,
// pinned in unit/test_mesh_frozen.cpp); its signature is unchanged. What this
// suite PINS where the design is silent:
//   - "does not count" means what L14's coincidence radius means for a skipped
//     node: never the argmax, never a void triangle's carve point, and not
//     counted in `uncovered` (a void triangle beside a seam would otherwise
//     carve on the seam, a sixth way onto it);
//   - a node on a constrained edge that is NOT frozen is scanned as today.
//
// Expected values are computed here by brute force with integer orientation,
// not read from the scan.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/scan.hpp>

#include "refinement_fixtures.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <utility>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::MeshVertex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::NodeLocation;
using terrain::refinement::scan;
using terrain::refinement::ScanResult;

namespace {

constexpr std::size_t kN = 9;
constexpr float kNoData = -9999.0f;

RasterGeometry square() { return RasterGeometry{0.0, 0.0, 1.0, 1.0, kN, kN}; }

Raster<float> dem(std::map<std::pair<int, int>, float> at, float base = 0.0f) {  // (row, col) -> z
    std::vector<float> z(kN * kN, base);
    for (const auto& [rc, v] : at) z[static_cast<std::size_t>(rc.first) * kN + static_cast<std::size_t>(rc.second)] = v;
    return Raster<float>{square(), std::move(z), kNoData};
}

// Two triangles on the seam col 4, from `top` to `bottom` (col, row):
//   L = (top, (0.25 or 0, 8 or 7.75), bottom), R = (top, bottom, (8, 8)).
// The seam is L's edge 2 (bottom -> top) and R's edge 0 (top -> bottom), mask
// 4; every other side is outline, mask 1.
LatticeMesh seam_pair(MeshVertex top, MeshVertex bottom, MeshVertex left) {
    const std::vector<MeshVertex> v{top, left, bottom, MeshVertex{8, 8}};
    const std::vector<TriangleIndices> t{{0, 1, 2}, {0, 2, 3}};
    const std::vector<std::uint8_t> bits{0b111, 0b111};
    const std::vector<std::array<std::uint32_t, 3>> masks{{1, 1, 4}, {4, 1, 1}};
    auto m = LatticeMesh::build(v, t, bits, masks);
    REQUIRE(m.has_value());
    return std::move(*m);
}

LatticeMesh on_node_pair() { return seam_pair({4, 0}, {4, 8}, {0, 8}); }
LatticeMesh off_node_pair() { return seam_pair({4, 0.5}, {4, 7.5}, {0.25, 7.75}); }

constexpr std::uint32_t kL = 0, kR = 1;

}  // namespace

TEST_CASE("FS1: a node exactly on a frozen edge is never the argmax", "[refinement][scan][frozen]") {
    // Node (row 4, col 4) is on the seam, 50 above flat ground. Unfrozen it is
    // L's argmax on edge 2 (the control); frozen, L's set holds nothing over 0.
    const bool on_node = GENERATE(true, false);
    CAPTURE(on_node);
    LatticeMesh m = on_node ? on_node_pair() : off_node_pair();
    const auto d = dem({{{4, 4}, 50.0f}});

    const ScanResult before = scan(d, m, kL);
    REQUIRE(before.node == LatticeVertex{4, 4});
    REQUIRE(before.where == NodeLocation::Edge2);
    REQUIRE(before.max_error == 50.0);

    m.set_frozen_mask(4u);
    for (const auto t : {kL, kR}) {
        CAPTURE(t);
        const ScanResult r = scan(d, m, t);
        REQUIRE_FALSE(r.node.has_value());
        REQUIRE(r.max_error == 0.0);
        REQUIRE_FALSE(r.is_void);
    }
}

TEST_CASE("FS1: a node off the frozen edge still wins when it is worst", "[refinement][scan][frozen]") {
    // The seam node is 50 high, an interior node of L 30: the interior one is
    // L's argmax once the seam node is not counted.
    const bool on_node = GENERATE(true, false);
    CAPTURE(on_node);
    LatticeMesh m = on_node ? on_node_pair() : off_node_pair();
    m.set_frozen_mask(4u);
    const ScanResult r = scan(dem({{{4, 4}, 50.0f}, {{4, 3}, 30.0f}}), m, kL);
    REQUIRE(r.node == LatticeVertex{4, 3});
    REQUIRE(r.where == NodeLocation::Inside);
    REQUIRE(r.max_error == 30.0);
}

TEST_CASE("FS1: a node on a constrained edge that is not frozen is scanned as today", "[refinement][scan][frozen]") {
    // Node (row 8, col 2) lies on L's outline edge (0, 8)-(4, 8), mask 1;
    // only the seam (mask 4) is frozen.
    LatticeMesh m = on_node_pair();
    m.set_frozen_mask(4u);
    const ScanResult r = scan(dem({{{8, 2}, 20.0f}, {{4, 4}, 50.0f}}), m, kL);
    REQUIRE(r.node == LatticeVertex{8, 2});
    REQUIRE(r.where == NodeLocation::Edge1);
    REQUIRE(r.max_error == 20.0);
}

TEST_CASE("FS2: a void triangle never carves on its frozen edge", "[refinement][scan][frozen]") {
    // The seam's top end (row 0, col 4) is NoData, so L is void. The valid
    // node nearest it in L is (1, 4), ON the seam: the carve point unfrozen,
    // never frozen. Frozen, the carve point is the nearest valid node of L off
    // the seam, and the seam's nodes are not in `uncovered`.
    LatticeMesh m = on_node_pair();
    const auto d = dem({{{0, 4}, kNoData}}, 1.0f);
    const ScanResult before = scan(d, m, kL);
    REQUIRE(before.is_void);
    REQUIRE(before.node == LatticeVertex{1, 4});

    // L's node set, by brute force: closed triangle (4,0) (0,8) (4,8) in
    // (col, row), its corners excluded; split into seam (col 4) and the rest.
    using refinement_fixtures::RC;
    const RC a{0, 4}, b{8, 0}, c{8, 4};
    std::size_t on_seam = 0, off_seam = 0;
    RC nearest{-1, -1};
    std::int64_t best = 1 << 30;
    for (std::int64_t r = 0; r < static_cast<std::int64_t>(kN); ++r)
        for (std::int64_t col = 0; col < static_cast<std::int64_t>(kN); ++col) {
            const RC p{r, col};
            if (!refinement_fixtures::in_closed(a, b, c, p) || p == a || p == b || p == c) continue;
            if (col == 4) {
                ++on_seam;
                continue;
            }
            ++off_seam;
            const std::int64_t d2 = (r - a.row) * (r - a.row) + (col - a.col) * (col - a.col);
            if (d2 < best) {  // row-major walk: ties keep the smallest (row, col)
                best = d2;
                nearest = p;
            }
        }
    REQUIRE(before.uncovered == on_seam + off_seam);

    m.set_frozen_mask(4u);
    const ScanResult r = scan(d, m, kL);
    REQUIRE(r.is_void);
    REQUIRE(r.node == LatticeVertex{static_cast<std::uint32_t>(nearest.row), static_cast<std::uint32_t>(nearest.col)});
    REQUIRE(r.node->col != 4u);
    REQUIRE(r.uncovered == off_seam);
}

TEST_CASE("FS3: a frozen mask that meets no edge leaves every scan result as it was", "[refinement][scan][frozen]") {
    // K1 at the scan: only triangles with a frozen edge take the new path.
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    const bool on_node = GENERATE(true, false);
    CAPTURE(seed, on_node);
    const Raster<float> d{square(), refinement_fixtures::rough_dem(kN, kN, seed)};
    LatticeMesh m = on_node ? on_node_pair() : off_node_pair();
    const std::array<ScanResult, 2> before{scan(d, m, kL), scan(d, m, kR)};
    m.set_frozen_mask(1u << 31);
    for (const auto t : {kL, kR}) {
        const ScanResult r = scan(d, m, t);
        REQUIRE(r.node == before[t].node);
        REQUIRE(r.where == before[t].where);
        REQUIRE(r.max_error == before[t].max_error);
        REQUIRE(r.is_void == before[t].is_void);
        REQUIRE(r.uncovered == before[t].uncovered);
    }
}
