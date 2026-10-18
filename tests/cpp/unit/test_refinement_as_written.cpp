// The as-written swap pass of refine_points, as three named functions in
// terrain::refinement::detail, tested without running the refinement loop:
//
// - written_points: each vertex as the loop writes it, (x, -y), and whether
//   that written point is an exact translate of its grid point;
// - written_kinds: per triangle, Faithful (three faithful corners), Candidate
//   (counter-clockwise as written) or Flat (flat or turned over as written);
// - legalise_as_written: the pass itself. It leaves the mesh as given when no
//   edge between candidates is illegal as written.
//
// The fixtures are fed through detail::to_lattice and lattice_frame as
// refine_points feeds its start mesh. Every expected value below was measured
// on the loop's code before the split, by a scratch copy of its body run on
// the same fixture; the quad's values match what
// prop_refinement_written_delaunay.cpp asserts through refine_points.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/refine_points.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <variant>
#include <vector>

namespace {

using terrain::IndexedMesh2;
using terrain::pred::Orientation;
using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::LatticeMesh;
using terrain::pred::DefaultKernel;
using terrain::raster::RasterGeometry;
using terrain::refinement::detail::AsWrittenPass;
using terrain::refinement::detail::legalise_as_written;
using terrain::refinement::detail::to_lattice;
using terrain::refinement::detail::written_kinds;
using terrain::refinement::detail::written_points;
using terrain::refinement::detail::WrittenKind;
using terrain::refinement::detail::WrittenPoints;

using Edges = std::vector<std::array<std::uint32_t, 2>>;

// A start mesh on the grid as refine_points builds it: every edge listed is
// constrained with mask 1, nothing frozen.
struct Fixture {
    RasterGeometry g;
    IndexedMesh2 start;
    LatticeMesh m;
};

Fixture fixture(const RasterGeometry& g, std::vector<Point2> xy, std::vector<TriangleIndices> tris, const Edges& edges) {
    const std::size_t nt = tris.size();
    IndexedMesh2 start{std::move(xy), std::move(tris), std::vector<std::uint8_t>(nt, 0)};
    const std::vector<std::uint32_t> masks(edges.size(), 1);
    auto built = to_lattice(g, start, std::span<const std::array<std::uint32_t, 2>>{edges},
                            std::span<const std::uint32_t>{masks});
    REQUIRE(std::holds_alternative<LatticeMesh>(built));
    auto m = std::get<LatticeMesh>(std::move(built));
    m.set_frozen_mask(0);
    return Fixture{g, std::move(start), std::move(m)};
}

terrain::mesh::LatticeFrame frame_of(const RasterGeometry& g) {
    return terrain::mesh::lattice_frame(g.delta_x(), g.delta_y(), g.rows(), g.cols());
}

std::vector<TriangleIndices> triangles_of(const LatticeMesh& m) {
    return {m.triangles().begin(), m.triangles().end()};
}

// The quad of prop_refinement_written_delaunay.cpp's T3: A, B, C on the line
// x = y as written (exact orientation 0), counter-clockwise on the grid frame
// (about +9.25e-20, which that test checks); D makes (A, C, D)
// counter-clockwise. All five edges are constrained.
Fixture collinear_quad() {
    const RasterGeometry g{0.0, 30.0, 3.0, 3.0, 11, 11};
    const Point2 a{0.5, 0.5}, b{0.51, 0.51}, c{0.52, 0.52}, d{0.4, 0.6};
    REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::Collinear);
    REQUIRE(DefaultKernel::orient2d(a, c, d) == Orientation::CounterClockwise);
    return fixture(g, {a, b, c, d}, {{0, 1, 2}, {0, 2, 3}}, {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {2, 3}});
}

// One grid cell's four corner nodes (x_min 0, y_max 30, 3 m spacing), split
// by an unconstrained diagonal: every written point is its grid point exactly.
Fixture node_square() {
    const RasterGeometry g{0.0, 30.0, 3.0, 3.0, 11, 11};
    return fixture(g, {{0.0, 0.0}, {3.0, 0.0}, {3.0, 3.0}, {0.0, 3.0}}, {{0, 1, 2}, {0, 2, 3}},
                   {{0, 1}, {1, 2}, {2, 3}, {3, 0}});
}

}  // namespace

TEST_CASE("written_points writes each start vertex as (x, -y), and an off-node point is not faithful",
          "[refine_points][as_written]") {
    const auto f = collinear_quad();
    const WrittenPoints w = written_points(f.m, f.start, f.g, frame_of(f.g));
    REQUIRE(w.at.size() == 4);
    REQUIRE(w.faithful.size() == 4);
    for (std::size_t i = 0; i < 4; ++i) {
        CAPTURE(i);
        CHECK(w.at[i].col == f.start.vertices()[i].x);
        CHECK(w.at[i].row == -f.start.vertices()[i].y);
        CHECK(w.faithful[i] == 0);
    }
}

TEST_CASE("written_points marks grid nodes faithful", "[refine_points][as_written]") {
    const auto f = node_square();
    const WrittenPoints w = written_points(f.m, f.start, f.g, frame_of(f.g));
    CHECK(w.faithful == std::vector<char>(4, 1));
}

TEST_CASE("written_kinds: the collinear triangle is Flat, its counter-clockwise neighbour a Candidate",
          "[refine_points][as_written]") {
    const auto f = collinear_quad();
    const auto kinds = written_kinds(f.m, written_points(f.m, f.start, f.g, frame_of(f.g)));
    CHECK(kinds == std::vector<WrittenKind>{WrittenKind::Flat, WrittenKind::Candidate});
}

TEST_CASE("written_kinds: triangles with three faithful corners are Faithful", "[refine_points][as_written]") {
    const auto f = node_square();
    const auto kinds = written_kinds(f.m, written_points(f.m, f.start, f.g, frame_of(f.g)));
    CHECK(kinds == std::vector<WrittenKind>(2, WrittenKind::Faithful));
}

TEST_CASE("legalise_as_written counts the flat triangle and leaves the mesh as given",
          "[refine_points][as_written]") {
    auto f = collinear_quad();
    const auto before = triangles_of(f.m);
    std::vector<char> touched(f.m.triangle_count(), 0);
    const AsWrittenPass pass =
        legalise_as_written(f.m, written_points(f.m, f.start, f.g, frame_of(f.g)), 0, touched);
    CHECK(pass.flat == 1);
    CHECK_FALSE(pass.changed);
    CHECK(pass.flips == 0);
    CHECK_FALSE(pass.left_illegal.has_value());
    CHECK(triangles_of(f.m) == before);
    CHECK(touched == std::vector<char>(before.size(), 0));
}

TEST_CASE("legalise_as_written on grid nodes finds nothing flat and changes nothing",
          "[refine_points][as_written]") {
    auto f = node_square();
    const auto before = triangles_of(f.m);
    std::vector<char> touched(f.m.triangle_count(), 0);
    const AsWrittenPass pass =
        legalise_as_written(f.m, written_points(f.m, f.start, f.g, frame_of(f.g)), 0, touched);
    CHECK(pass.flat == 0);
    CHECK_FALSE(pass.changed);
    CHECK(pass.flips == 0);
    CHECK_FALSE(pass.left_illegal.has_value());
    CHECK(triangles_of(f.m) == before);
    CHECK(touched == std::vector<char>(before.size(), 0));
}
