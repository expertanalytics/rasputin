// Increment 15f-1 (docs/increments/15f-edge-strip.md, D2 step 1 and D5): the
// lattice position of one world point, extracted from detail::to_lattice as
//
//   detail::lattice_position(const raster::RasterGeometry& g, Point2 p) -> mesh::MeshVertex
//
// "a node exactly when g.node of the rounded position gives p back bit for
// bit, otherwise the clamped fractional position. to_lattice then calls it."
// The extraction changes no behaviour, so the suite pins two things: the rule
// itself, on hand-built points, and agreement bit for bit with what
// to_lattice puts in the lattice mesh for the same points. The second is what
// lets the generator claim its ends are the loop's vertices.
//
// The precondition (p inside the node rectangle) is to_lattice's refusal and
// is not exercised here; constraint_check_points' refusals are in
// test_refinement_constraint_points.cpp.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/refine.hpp>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <variant>
#include <vector>

using terrain::Point2;
using terrain::mesh::MeshVertex;
using terrain::raster::RasterGeometry;
using terrain::refinement::detail::lattice_position;

namespace {

// UTM-scale, non-square cells, so a row/column swap cannot pass.
constexpr double kX0 = 500000.0, kY0 = 7000000.0, kDx = 10.0, kDy = 5.0;
constexpr std::size_t kCols = 9, kRows = 7;

RasterGeometry grid() { return RasterGeometry{kX0, kY0, kDx, kDy, kCols, kRows}; }

// Exact for the dyadic (col, row) used below.
Point2 world(double col, double row) { return Point2{kX0 + col * kDx, kY0 - row * kDy}; }

}  // namespace

TEST_CASE("LP1: a node, bit for bit, is the integer lattice vertex", "[lattice_position]") {
    const auto g = grid();
    for (std::size_t r = 0; r < kRows; ++r)
        for (std::size_t c = 0; c < kCols; ++c) {
            const MeshVertex v = lattice_position(g, g.node({r, c}));
            CHECK(v.is_node());
            CHECK(v == MeshVertex{static_cast<double>(c), static_cast<double>(r)});
        }
}

TEST_CASE("LP2: an off-node point is its fractional (col, row)", "[lattice_position]") {
    const auto g = grid();
    CHECK(lattice_position(g, world(0.625, 2.25)) == MeshVertex{0.625, 2.25});
    CHECK(lattice_position(g, world(3.0, 0.5)) == MeshVertex{3.0, 0.5});  // on a column line
    CHECK(lattice_position(g, world(4.5, 6.0)) == MeshVertex{4.5, 6.0});  // on the last row line
    CHECK(lattice_position(g, world(8.0, 3.75)) == MeshVertex{8.0, 3.75});  // on the far side
}

TEST_CASE("LP3: a hair off a node is not the node", "[lattice_position]") {
    const auto g = grid();
    const Point2 node = g.node({3, 4});
    const Point2 east{std::nextafter(node.x, std::numeric_limits<double>::infinity()), node.y};
    const Point2 south{node.x, std::nextafter(node.y, -std::numeric_limits<double>::infinity())};
    const MeshVertex ve = lattice_position(g, east);
    const MeshVertex vs = lattice_position(g, south);
    CHECK_FALSE(ve.is_node());
    CHECK_FALSE(vs.is_node());
    CHECK(ve.col > 4.0);
    CHECK(ve.row == 3.0);
    CHECK(vs.row > 3.0);
    CHECK(vs.col == 4.0);
}

TEST_CASE("LP4: agrees bit for bit with the vertices to_lattice builds", "[lattice_position]") {
    const auto g = grid();
    // Four corners of the node rectangle (nodes) and points that are off-node
    // in several ways: dyadic, non-dyadic, a hair off a node, on grid lines.
    const Point2 hair{std::nextafter(g.node({2, 2}).x, 0.0), g.node({2, 2}).y};
    const std::vector<Point2> xy{
        g.node({0, 0}), g.node({kRows - 1, 0}), g.node({kRows - 1, kCols - 1}),
        g.node({0, kCols - 1}), world(1.3, 1.7), world(2.0, 4.1), world(6.7, 3.0),
        hair, world(5.5, 5.25),
    };
    // Two CCW triangles over the four corners; the other vertices are
    // unreferenced. Only the vertices matter here, but to_lattice builds a mesh.
    const std::vector<terrain::TriangleIndices> tris{{0, 1, 2}, {0, 2, 3}};
    terrain::IndexedMesh2 start{xy, tris, std::vector<std::uint8_t>(tris.size(), 0)};
    const std::vector<std::array<std::uint32_t, 2>> edges;
    const std::vector<std::uint32_t> masks;
    auto built = terrain::refinement::detail::to_lattice(g, start, edges, masks);
    REQUIRE(std::holds_alternative<terrain::mesh::LatticeMesh>(built));
    const auto& m = std::get<terrain::mesh::LatticeMesh>(built);
    REQUIRE(m.vertices().size() == xy.size());
    for (std::size_t i = 0; i < xy.size(); ++i) {
        INFO("vertex " << i);
        const MeshVertex v = lattice_position(g, xy[i]);
        CHECK(std::bit_cast<std::uint64_t>(v.col) == std::bit_cast<std::uint64_t>(m.vertices()[i].col));
        CHECK(std::bit_cast<std::uint64_t>(v.row) == std::bit_cast<std::uint64_t>(m.vertices()[i].row));
    }
}

// ---------------------------------------------------------------------------
// Increment 15f-4 (docs/increments/15f-edge-strip.md, "Settled after 15f-3's
// acceptance", A2 and A4, B4): to_lattice's constraint lookup, pinned against
// today's std::map before A2 replaces it with a sorted vector. A guard: it
// passes before the change and must still pass after it.
// ---------------------------------------------------------------------------

namespace {

// The (constrained, mask) of the undirected edge (a, b) on every side that
// holds it, in triangle order.
std::vector<std::pair<bool, std::uint32_t>> sides(const terrain::mesh::LatticeMesh& m,
                                                  std::uint32_t a, std::uint32_t b) {
    std::vector<std::pair<bool, std::uint32_t>> out;
    for (std::size_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto p = m.triangles()[t][k], q = m.triangles()[t][(k + 1) % 3];
            if ((p == a && q == b) || (p == b && q == a))
                out.emplace_back(m.is_constrained(t, k), m.mask(t, k));
        }
    return out;
}

using Side = std::pair<bool, std::uint32_t>;
using Edges = std::vector<std::array<std::uint32_t, 2>>;
using MaskList = std::vector<std::uint32_t>;

// Four nodes of grid(): 0 top-left, 1 bottom-left, 2 bottom-right, 3
// top-right; triangles (0, 1, 2) and (0, 2, 3) share the diagonal 0-2.
terrain::mesh::LatticeMesh square(const Edges& edges, const MaskList& masks) {
    const auto g = grid();
    const std::vector<Point2> xy{g.node({0, 0}), g.node({4, 0}), g.node({4, 6}), g.node({0, 6})};
    const std::vector<terrain::TriangleIndices> tris{{0, 1, 2}, {0, 2, 3}};
    // The start's own constraint bits are not what to_lattice reads: set them
    // all, so a lookup that fell back to them would show.
    const terrain::IndexedMesh2 start{xy, tris, std::vector<std::uint8_t>(tris.size(), 0b111)};
    auto built = terrain::refinement::detail::to_lattice(g, start, edges, masks);
    REQUIRE(std::holds_alternative<terrain::mesh::LatticeMesh>(built));
    return std::get<terrain::mesh::LatticeMesh>(std::move(built));
}

}  // namespace

TEST_CASE("B4: to_lattice gives an edge listed twice the later mask", "[lattice_position][15f-4]") {
    SECTION("the same direction") {
        const auto m = square({{0, 1}, {0, 1}}, {5, 9});
        CHECK(sides(m, 0, 1) == std::vector<Side>{{true, 9}});
    }
    SECTION("reversed") {
        const auto m = square({{1, 0}, {0, 1}}, {5, 9});
        CHECK(sides(m, 0, 1) == std::vector<Side>{{true, 9}});
    }
    SECTION("reversed, the other way round, with other edges between") {
        const auto m = square({{0, 1}, {2, 3}, {1, 0}}, {9, 7, 5});
        CHECK(sides(m, 0, 1) == std::vector<Side>{{true, 5}});
        CHECK(sides(m, 2, 3) == std::vector<Side>{{true, 7}});
    }
    SECTION("an interior edge: both sides get the later mask") {
        const auto m = square({{2, 0}, {0, 2}, {2, 0}}, {1, 2, 4});
        CHECK(sides(m, 0, 2) == std::vector<Side>{{true, 4}, {true, 4}});
    }
    SECTION("a later mask of 0 still constrains the edge") {
        const auto m = square({{0, 1}, {1, 0}}, {5, 0});
        CHECK(sides(m, 0, 1) == std::vector<Side>{{true, 0}});
    }
}

TEST_CASE("B4: to_lattice ignores edges past masks.size()", "[lattice_position][15f-4]") {
    SECTION("an edge with no mask is not constrained") {
        const auto m = square({{0, 1}, {2, 3}, {3, 0}}, {5, 6});
        CHECK(sides(m, 0, 1) == std::vector<Side>{{true, 5}});
        CHECK(sides(m, 2, 3) == std::vector<Side>{{true, 6}});
        CHECK(sides(m, 3, 0) == std::vector<Side>{{false, 0}});
    }
    SECTION("a repeat past masks.size() does not override the earlier mask") {
        const auto m = square({{0, 1}, {1, 0}}, {5});
        CHECK(sides(m, 0, 1) == std::vector<Side>{{true, 5}});
    }
    SECTION("no masks: nothing is constrained") {
        const auto m = square({{0, 1}, {0, 2}}, {});
        for (const auto& [a, b] : Edges{{0, 1}, {1, 2}, {2, 0}, {2, 3}, {3, 0}})
            for (const auto& s : sides(m, a, b)) CHECK(s == Side{false, 0});
    }
    SECTION("masks past edges.size() are ignored, and so is a non-edge") {
        const auto m = square({{1, 3}, {1, 2}}, {5, 6, 7, 8});
        CHECK(sides(m, 1, 2) == std::vector<Side>{{true, 6}});
        CHECK(sides(m, 1, 3).empty());
        for (const auto& [a, b] : Edges{{0, 1}, {2, 0}, {2, 3}, {3, 0}})
            for (const auto& s : sides(m, a, b)) CHECK(s == Side{false, 0});
    }
}
