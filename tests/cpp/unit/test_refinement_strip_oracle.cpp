// Increment 15f-2 (docs/increments/15f-edge-strip.md, "Tests for @tester"):
// the oracles of the edge-strip suite (support/strip_oracle.hpp), shown to
// fail. Everything here runs against refine() and code that exists before
// 15f-2, so it builds and passes while prop_refinement_edge_strip.cpp, which
// calls refine_strip, does not compile.
//
//   - ES1, the defect: refine() alone leaves the edge strip's crossings over
//     the tolerance on the bump fixture (Surprise 3 in miniature), while its
//     own guarantee at DEM nodes holds. The repair half is in the red suite.
//   - ES2's control: the strip oracle reports nothing on a planar DEM, and
//     reports every point once the constraint vertices' z are shifted by
//     twice the tolerance.
//   - The DEM-node oracle (tester.md section 3D's tolerance oracle, E2) and
//     the constrained-Delaunay oracle, each planted.
//   - The oracle's own check points, by hand, on two edges.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>

#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <set>
#include <span>
#include <vector>

using namespace strip_oracle;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;

namespace {

Mesh refined(const Raster<float>& dem, const Start& s, double tol) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = 1;
    const auto r = refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                          std::span<const std::uint32_t>{s.masks}, o);
    REQUIRE(r.ok());
    return mesh_of(r);
}

}  // namespace

TEST_CASE("ES1 (the defect): refine alone leaves crossings of the bump side over the tolerance",
          "[edge_strip][oracle][ES1]") {
    const auto dem = bump_dem();
    const Start s = bump_start();
    const double tol = 1.0;
    const Mesh out = refined(dem, s, tol);
    REQUIRE(out.vertices.size() == s.mesh.vertices().size());  // nothing inside needed a node

    const auto pts = ruled_points(dem, s.mesh.vertices(), s.edges);
    const auto f = strip_findings(dem.geometry(), pts, out, tol);
    CHECK(f.unfiled == 0);
    CHECK(f.on_void == 0);
    REQUIRE(f.over > 0);
    CHECK(f.worst == 5.0);  // (0 + 10) / 2 at the crossings of columns 3 and 4, against 0
    const std::set<Lat> over(f.over_at.begin(), f.over_at.end());
    CHECK(over.contains(Lat{3.0, 3.5}));
    CHECK(over.contains(Lat{4.0, 3.5}));

    // refine's own guarantee holds: the DEM-node oracle finds nothing.
    const auto n = node_findings(dem, out, tol);
    CHECK(n.not_ccw == 0);
    CHECK(n.over == 0);
    CHECK(delaunay_violations(dem.geometry(), out) == 0);
}

TEST_CASE("ES2's control: the strip oracle is silent on a plane and fires when the constraints are shifted",
          "[edge_strip][oracle][ES2]") {
    const auto g = exact_geometry();
    const std::uint32_t seed = GENERATE(1u, 2u, 3u);
    CAPTURE(seed);
    const auto dem = terrain_dem(g, Terrain::Plane, seed);
    const Start s = jittered(g, seed);
    const double tol = 0.25;
    const Mesh out = refined(dem, s, tol);
    const auto pts = ruled_points(dem, out.vertices, out.edges);
    REQUIRE(pts.size() > 100);

    const auto clean = strip_findings(g, pts, out, tol);
    CHECK(clean.unfiled == 0);
    CHECK(clean.on_void == 0);
    CHECK(clean.over == 0);
    CHECK(clean.worst < 1e-9);

    Mesh planted = out;
    std::set<std::uint32_t> on_constraint;
    for (const auto& e : out.edges) on_constraint.insert({e[0], e[1]});
    for (const auto v : on_constraint) planted.z[v] += 2 * tol;
    const auto shifted = strip_findings(g, pts, planted, tol);
    CHECK(shifted.over == pts.size());
}

TEST_CASE("The DEM-node oracle is silent on refine's output and fires when z is shifted",
          "[edge_strip][oracle][ES3]") {
    const auto g = exact_geometry();
    const auto kind = GENERATE(Terrain::Rough, Terrain::Smooth, Terrain::SmoothWithHole);
    const double tol = GENERATE(0.0, 0.5, 2.0);
    CAPTURE(static_cast<int>(kind), tol);
    const auto dem = terrain_dem(g, kind, 1);
    const Start s = jittered(g, 1);
    const Mesh out = refined(dem, s, tol);
    const auto n = node_findings(dem, out, tol);
    CHECK(n.not_ccw == 0);
    CHECK(n.over == 0);
    CHECK(delaunay_violations(g, out) == 0);

    // At tolerance 0 every valid node is a vertex and the oracle skips
    // vertices, so a shift has nothing to show there.
    if (tol > 0.0) {
        Mesh planted = out;
        for (auto& z : planted.z) z += 2 * tol;
        CHECK(node_findings(dem, planted, tol).over > 0);
    }
}

TEST_CASE("The constrained-Delaunay oracle fires on the long diagonal unless it is a constraint",
          "[edge_strip][oracle]") {
    // A diamond in (col, row): a (0, 2), b (3, 0), c (6, 2), d (3, 4). In the
    // producer's frame (2 col, -row) it is 12 wide and 4 tall, so the
    // diagonal a-c is the non-Delaunay one: b lies inside the circle of a, d, c.
    const auto g = exact_geometry(9, 9);
    Mesh m;
    m.vertices = {world(g, 0, 2), world(g, 3, 0), world(g, 6, 2), world(g, 3, 4)};
    m.z = {0, 0, 0, 0};
    m.valid = {1, 1, 1, 1};
    m.triangles = {{0, 3, 2}, {0, 2, 1}};
    CHECK(node_findings(Raster<float>{g, std::vector<float>(81, 0.0f)}, m, 0.0).not_ccw == 0);
    CHECK(delaunay_violations(g, m) == 2);
    m.edges = {{0, 2}};
    m.masks = {1};
    CHECK(delaunay_violations(g, m) == 0);
}

TEST_CASE("The oracle's ruled points, by hand", "[edge_strip][oracle]") {
    const auto dem = bump_dem();
    const auto g = dem.geometry();
    SECTION("the bump's bottom side: seven column crossings and eight midpoints") {
        const Start s = bump_start();
        // From the lower index, v0 at col 7.5, to v3 at col 0.5; read by column here.
        auto pts = ruled_points(dem, s.mesh.vertices(), Edges{{0, 3}});
        std::sort(pts.begin(), pts.end(), [](const OraclePoint& x, const OraclePoint& y) { return x.at < y.at; });
        std::vector<Lat> want;
        double prev = 0.5;
        for (int k = 1; k <= 7; ++k) {
            want.push_back(Lat{(prev + k) / 2, 3.5});
            want.push_back(Lat{static_cast<double>(k), 3.5});
            prev = k;
        }
        want.push_back(Lat{(7.0 + 7.5) / 2, 3.5});
        REQUIRE(pts.size() == want.size());
        for (std::size_t i = 0; i < pts.size(); ++i) {
            CAPTURE(i);
            CHECK(pts[i].at == want[i]);
        }
        CHECK(pts[5].z == 5.0);   // the crossing of column 3
        CHECK(pts[4].z == 2.5);   // the midpoint at col 2.5: a quarter of node (4, 3)
    }
    SECTION("a diagonal through three nodes: each node once, and four midpoints") {
        const std::vector<terrain::Point2> v{world(g, 0.5, 0.5), world(g, 3.5, 3.5)};
        const auto pts = ruled_points(dem, v, Edges{{0, 1}});
        REQUIRE(pts.size() == 7);
        const std::vector<Lat> want{{0.75, 0.75}, {1, 1}, {1.5, 1.5}, {2, 2}, {2.5, 2.5}, {3, 3}, {3.25, 3.25}};
        for (std::size_t i = 0; i < pts.size(); ++i) {
            CAPTURE(i);
            CHECK(pts[i].at.col == want[i].col);
            CHECK(std::abs(pts[i].at.row - want[i].row) <= 1e-15);
        }
    }
}
