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
