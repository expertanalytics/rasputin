// Increment 16 (docs/increments/16-domain-polygon.md, R2, "Tests for
// @tester"): the scan on triangles with off-node vertices. INVARIANT-CRITICAL:
// the tolerance guarantee near the domain boundary is decided here, so
// @reviewer mutation-tests this suite. S4 (node-only triangles bit-identical
// to 14b) is test_refinement_scan.cpp, unchanged: it must still pass.
//
// Interface assumed (the design's "The vertex type"; names it leaves open are
// chosen here and stated in the handback):
//   struct terrain::mesh::MeshVertex { double col; double row; bool is_node() const; };
//     built here as MeshVertex{col, row}, which works for an aggregate and for
//     a (double, double) constructor alike;
//   LatticeMesh::build(std::vector<MeshVertex>, triangles, constrained, masks);
//   scan(dem, mesh, t) keeps its signature; ScanResult::node stays a
//     LatticeVertex (only nodes are ever inserted), and an off-node vertex's z
//     is raster::bilinear at its position, derived from the DEM and the
//     vertex alone.
//
// THE ORACLE IS THE PRODUCER'S RELATION (design S3): membership is three
// DefaultKernel orientations on (col, -row), "on an edge" is Collinear; the
// plane uses double barycentric weights; an off-node vertex's z is bilinear
// in the fractional frame. It shares no code with scan.hpp.
//
// Fixture DEM: 9 x 9 (17 x 17 for the random cases), geometry of
// refinement_fixtures (x_min 500 000, dx 10, dy 5). Every hand-placed
// coordinate is dyadic, so the world point RasterGeometry gives back is exact
// and bilinear on the plane z = 3 col - 2 row + 7 is the plane exactly.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/scan.hpp>

#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <vector>

using Catch::Matchers::WithinAbs;
using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::MeshVertex;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::raster::Raster;
using terrain::refinement::NodeLocation;
using terrain::refinement::scan;

namespace {

const float kNaN = std::numeric_limits<float>::quiet_NaN();

MeshVertex mv(double col, double row) { return MeshVertex{col, row}; }

std::vector<float> plane(std::size_t n) {
    std::vector<float> z(n * n);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            z[r * n + c] = static_cast<float>(3.0 * static_cast<double>(c) - 2.0 * static_cast<double>(r) + 7.0);
    return z;
}

Raster<float> dem_of(std::size_t n, std::vector<float> z) {
    return Raster<float>{refinement_fixtures::geometry(n, n), std::move(z)};
}

LatticeMesh one(MeshVertex a, MeshVertex b, MeshVertex c) {
    auto m = LatticeMesh::build({a, b, c}, {TriangleIndices{0, 1, 2}}, {0}, {std::array<std::uint32_t, 3>{}});
    REQUIRE(m.has_value());
    return *m;
}

Point2 frame(double col, double row) { return Point2{col, -row}; }

// Bilinear in the fractional frame, written from R0's text. nullopt on NoData.
std::optional<double> z_at(const Raster<float>& dem, double col, double row) {
    const auto n_r = dem.geometry().rows(), n_c = dem.geometry().cols();
    const auto r0 = std::min(static_cast<std::size_t>(std::floor(row)), n_r - 2);
    const auto c0 = std::min(static_cast<std::size_t>(std::floor(col)), n_c - 2);
    const double ty = row - static_cast<double>(r0), tx = col - static_cast<double>(c0);
    if (tx == 0.0 && ty == 0.0) {  // a node: its own value (Z1)
        if (dem.is_nodata({r0, c0})) return std::nullopt;
        return static_cast<double>(dem.value_at({r0, c0}));
    }
    double z = 0.0;
    for (std::size_t dr = 0; dr < 2; ++dr)
        for (std::size_t dc = 0; dc < 2; ++dc) {
            if (dem.is_nodata({r0 + dr, c0 + dc})) return std::nullopt;
            z += static_cast<double>(dem.value_at({r0 + dr, c0 + dc})) * (dr ? ty : 1.0 - ty)
               * (dc ? tx : 1.0 - tx);
        }
    return z;
}

struct Expected {
    double max_error = 0.0;
    std::optional<LatticeVertex> node;
    NodeLocation where = NodeLocation::Inside;
    std::vector<std::array<double, 3>> errors;  // (row, col, error) of every member
};

// The brute-force oracle over every node of the grid. Three valid vertices only.
Expected oracle(const Raster<float>& dem, std::array<std::array<double, 2>, 3> v) {  // (col, row)
    std::array<Point2, 3> f{};
    std::array<double, 3> z{};
    for (unsigned k = 0; k < 3; ++k) {
        f[k] = frame(v[k][0], v[k][1]);
        const auto zk = z_at(dem, v[k][0], v[k][1]);
        REQUIRE(zk.has_value());
        z[k] = *zk;
    }
    auto cross = [](Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); };
    const double two_a = cross(f[0], f[1], f[2]);
    Expected e;
    for (std::uint32_t r = 0; r < dem.geometry().rows(); ++r)
        for (std::uint32_t c = 0; c < dem.geometry().cols(); ++c) {
            const Point2 p = frame(c, r);
            std::array<Orientation, 3> o{};
            bool out = false;
            for (unsigned k = 0; k < 3; ++k) {
                o[k] = DefaultKernel::orient2d(f[k], f[(k + 1) % 3], p);
                out = out || o[k] == Orientation::Clockwise;
            }
            const bool vertex = p == f[0] || p == f[1] || p == f[2];
            if (out || vertex || dem.is_nodata({r, c})) continue;
            const double plane_z = (cross(p, f[1], f[2]) * z[0] + cross(f[0], p, f[2]) * z[1]
                                    + cross(f[0], f[1], p) * z[2]) / two_a;
            const double err = std::abs(static_cast<double>(dem.value_at({r, c})) - plane_z);
            e.errors.push_back({static_cast<double>(r), static_cast<double>(c), err});
            if (err > e.max_error) {
                e.max_error = err;
                e.node = LatticeVertex{r, c};
                e.where = o[0] == Orientation::Collinear   ? NodeLocation::Edge0
                          : o[1] == Orientation::Collinear ? NodeLocation::Edge1
                          : o[2] == Orientation::Collinear ? NodeLocation::Edge2
                                                           : NodeLocation::Inside;
            }
        }
    return e;
}

// Two triangles sharing the edge A-B on the line col = row + shift:
//   left  = (A, C, B), the shared edge is its edge 2;
//   right = (A, B, D), the shared edge is its edge 0.
// Every vertex is off-node and no vertex's cell touches (3,3).
struct Pair {
    std::array<double, 2> a, b, c{0.25, 6.75}, d{6.75, 0.25};
    explicit Pair(double shift) : a{0.625 + shift, 0.625}, b{6.5 + shift, 6.5} {}
    [[nodiscard]] LatticeMesh left() const { return one(mv(a[0], a[1]), mv(c[0], c[1]), mv(b[0], b[1])); }
    [[nodiscard]] LatticeMesh right() const { return one(mv(a[0], a[1]), mv(b[0], b[1]), mv(d[0], d[1])); }
};

}  // namespace

TEST_CASE("scan off-node: S1 a node exactly on an edge between two off-node vertices is in both sets",
          "[refinement][scan][offnode]") {
    auto z = plane(9);
    z[3 * 9 + 3] += 1.0f;
    const auto dem = dem_of(9, std::move(z));
    const Pair p{0.0};
    const auto l = scan(dem, p.left(), 0);
    const auto r = scan(dem, p.right(), 0);
    REQUIRE_FALSE(l.is_void);
    REQUIRE(l.node == LatticeVertex{3, 3});
    REQUIRE(l.where == NodeLocation::Edge2);
    REQUIRE_THAT(l.max_error, WithinAbs(1.0, 1e-9));
    REQUIRE(r.node == LatticeVertex{3, 3});
    REQUIRE(r.where == NodeLocation::Edge0);
    REQUIRE_THAT(r.max_error, WithinAbs(1.0, 1e-9));
}

TEST_CASE("scan off-node: S2 a node a nanometre off such an edge is in exactly one set",
          "[refinement][scan][offnode]") {
    // shift 2^-30 cell: the line is col = row + 2^-30, so (3,3) lies strictly
    // on the left triangle's side. Truncating the vertices to integers (the
    // int64 path on a mixed triangle) would put it on the edge of both.
    auto z = plane(9);
    z[3 * 9 + 3] += 1.0f;
    const auto dem = dem_of(9, std::move(z));
    const Pair p{std::ldexp(1.0, -30)};
    const auto l = scan(dem, p.left(), 0);
    const auto r = scan(dem, p.right(), 0);
    REQUIRE(l.node == LatticeVertex{3, 3});
    REQUIRE(l.where == NodeLocation::Inside);
    REQUIRE_THAT(l.max_error, WithinAbs(1.0, 1e-9));
    REQUIRE(r.node != std::optional<LatticeVertex>{LatticeVertex{3, 3}});
    REQUIRE(r.max_error < 1e-6);
}

TEST_CASE("scan off-node: the node an off-node vertex rounds to is still in the set",
          "[refinement][scan][offnode]") {
    // A = (0.625, 0.625) rounds to (1,1), which lies on edge A-B. An off-node
    // vertex compared to nodes by rounded (row, col) would drop it.
    auto z = plane(9);
    z[1 * 9 + 1] += 4.0f;
    const auto dem = dem_of(9, std::move(z));
    const Pair p{0.0};
    const auto l = scan(dem, p.left(), 0);
    const auto want = oracle(dem, {p.a, p.c, p.b});
    REQUIRE(want.node == LatticeVertex{1, 1});
    REQUIRE(l.node == LatticeVertex{1, 1});
    REQUIRE(l.where == NodeLocation::Edge2);
    REQUIRE_THAT(l.max_error, WithinAbs(want.max_error, 1e-9));
}

TEST_CASE("scan off-node: the box reaches the last node below a fractional maximum",
          "[refinement][scan][offnode]") {
    // Vertices at 0.5 and 7.5: nodes 1..7 on both axes can be members. A box
    // with floor on both ends, or ceil on both, loses row/col 0 or 7 errors.
    auto z = plane(9);
    z[7 * 9 + 1] += 2.0f;  // (row 7, col 1): near the lower-left corner
    const auto dem = dem_of(9, std::move(z));
    const auto m = one(mv(0.5, 0.5), mv(0.5, 7.5), mv(7.5, 7.5));
    const auto s = scan(dem, m, 0);
    const auto want = oracle(dem, {{{0.5, 0.5}, {0.5, 7.5}, {7.5, 7.5}}});
    REQUIRE(want.node == LatticeVertex{7, 1});
    REQUIRE(s.node == LatticeVertex{7, 1});
    REQUIRE_THAT(s.max_error, WithinAbs(want.max_error, 1e-9));
}

TEST_CASE("scan off-node: an off-node vertex carries bilinear z not the rounded node's value",
          "[refinement][scan][offnode]") {
    // A plane DEM with one node bumped under an off-node vertex's cell: the
    // vertex's bilinear z moves by the bump times its weight, so the error at
    // every member depends on which z the vertex was given.
    auto z = plane(9);
    z[1 * 9 + 1] += 8.0f;
    const auto dem = dem_of(9, std::move(z));
    const std::array<std::array<double, 2>, 3> v{{{0.75, 1.25}, {0.5, 7.5}, {7.5, 7.5}}};
    const auto s = scan(dem, one(mv(v[0][0], v[0][1]), mv(v[1][0], v[1][1]), mv(v[2][0], v[2][1])), 0);
    const auto want = oracle(dem, v);
    REQUIRE(want.node.has_value());
    REQUIRE(s.node == want.node);
    REQUIRE(s.where == want.where);
    REQUIRE_THAT(s.max_error, WithinAbs(want.max_error, 1e-9));
}

TEST_CASE("scan off-node: S3 the error agrees with the brute-force oracle on random mixed triangles",
          "[refinement][scan][offnode]") {
    const std::uint32_t seed = GENERATE(1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u);
    const bool with_node_vertex = GENERATE(false, true);
    CAPTURE(seed, with_node_vertex);
    constexpr std::size_t n = 17;
    const auto dem = dem_of(n, refinement_fixtures::rough_dem(n, n, seed));
    std::mt19937 gen{seed * 7919u};
    auto coord = [&] { return static_cast<double>(gen() % 160000u) / 10000.0; };  // [0, 16)
    int checked = 0;
    for (int trial = 0; trial < 40; ++trial) {
        std::array<std::array<double, 2>, 3> v{{{coord(), coord()}, {coord(), coord()}, {coord(), coord()}}};
        if (with_node_vertex) v[0] = {std::round(v[0][0]), std::round(v[0][1])};
        const auto o = DefaultKernel::orient2d(frame(v[0][0], v[0][1]), frame(v[1][0], v[1][1]),
                                               frame(v[2][0], v[2][1]));
        if (o == Orientation::Collinear) continue;
        if (o == Orientation::Clockwise) std::swap(v[1], v[2]);
        CAPTURE(trial, v[0][0], v[0][1], v[1][0], v[1][1], v[2][0], v[2][1]);
        const auto s = scan(dem, one(mv(v[0][0], v[0][1]), mv(v[1][0], v[1][1]), mv(v[2][0], v[2][1])), 0);
        const auto want = oracle(dem, v);
        REQUIRE_FALSE(s.is_void);
        REQUIRE_THAT(s.max_error, WithinAbs(want.max_error, 1e-9));
        REQUIRE(s.node.has_value() == want.node.has_value());
        if (s.node) {
            // The argmax may differ from the oracle's only on a tie within rounding.
            const auto it = std::find_if(want.errors.begin(), want.errors.end(), [&](const auto& e) {
                return e[0] == s.node->row && e[1] == s.node->col;
            });
            REQUIRE(it != want.errors.end());
            REQUIRE_THAT((*it)[2], WithinAbs(want.max_error, 1e-9));
        }
        ++checked;
    }
    REQUIRE(checked >= 30);
}

TEST_CASE("scan off-node: an off-node vertex in a cell with a NoData corner makes a void triangle",
          "[refinement][scan][offnode][nodata]") {
    // A = (0.625, 0.625) sits in the cell whose corner (0,0) is NaN, so
    // bilinear refuses and A is a NoData vertex. The carve point is the valid
    // node nearest A in the fractional frame: (1,1), on edge A-B.
    auto z = plane(9);
    z[0] = kNaN;
    const auto dem = dem_of(9, std::move(z));
    const Pair p{0.0};
    const auto l = scan(dem, p.left(), 0);
    REQUIRE(l.is_void);
    REQUIRE(l.node == LatticeVertex{1, 1});
    REQUIRE(l.where == NodeLocation::Edge2);
    REQUIRE(l.uncovered > 0);
}

TEST_CASE("scan off-node: a NoData corner with zero weight still voids an off-node vertex",
          "[refinement][scan][offnode][nodata]") {
    // A = (0.625, 0.0) is off-node but on row 0: the cell's lower corners have
    // weight 0. bilinear refuses on any NoData corner (increment 12), so A is
    // still a NoData vertex (design R2, "NoData near an off-node vertex").
    auto z = plane(9);
    z[1 * 9 + 1] = kNaN;
    const auto dem = dem_of(9, std::move(z));
    const auto s = scan(dem, one(mv(0.625, 0.0), mv(0.5, 7.5), mv(7.5, 7.5)), 0);
    REQUIRE(s.is_void);
}
