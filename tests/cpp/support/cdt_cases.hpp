#pragma once

// Test-only scaffolding for the CDT suites: the named constraint sets each row
// of the degeneracy table in docs/increments/04-cdt.md is measured against, and
// the disjoint rectangle-with-holes generator the property suite draws from.
//
// WHY THIS HEADER NAMES NO KERNEL. Every fixture returns a PslgBuilder, not a
// Pslg, and the caller names the kernel through build_fixture<K>. That is not
// ceremony: test_cdt_backend_seam.cpp links terrain_headers and Catch2 and
// NOTHING ELSE -- that link line is the whole proof that the seam is real --
// so it must build its Pslg under the header-only FastKernel. A fixture header
// that included default_kernel.hpp would drag the compiled terrain_predicates
// target into the one suite whose value is that it needs no compiled target.
// ring.hpp and segment.hpp decline to name a default kernel for the same
// reason.
//
// WHY THE FIXTURES CANNOT REUSE valid_chain_specs FROM pslg_cases.hpp. That
// generator scatters breaklines and holes over one disc independently, so a
// generated breakline routinely crosses a generated hole ring. For increment 3
// that was irrelevant -- a Pslg promises nothing about disjointness -- but a
// crossing is a CdtStatus::NotNoded failure with no mesh, so a CDT property
// over those families would assert almost nothing. The generator below places
// every feature in its own grid cell with a margin, which is what makes
// "the input triangulates" a precondition rather than a coin toss. The
// hand-written fixtures still use pslg_cases' spans and rings.
//
// Nothing here asserts, matching pslg_cases.hpp: a helper that a REQUIRE lives
// inside reports failures at the helper's line rather than the test's.
// build_fixture throws instead, because a fixture the PSLG validator rejects is
// a bug in THIS header and should not read as a failure of the code under test.
//
// Every generator takes an explicit std::mt19937_64: a test seeds once and the
// whole sequence is reproducible.

#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/predicates/kernel.hpp>

#include <pslg_cases.hpp>
#include <ring_cases.hpp>

#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace terrain::test {

// ---------------------------------------------------------------------------
// Building a fixture
// ---------------------------------------------------------------------------

template <pred::GeometryKernel K>
[[nodiscard]] Pslg build_fixture(PslgBuilder b) {
    PslgBuildResult r = std::move(b).build<K>();
    if (!r.ok()) {
        throw std::logic_error{"cdt_cases fixture is not a valid Pslg:\n" + render(r)};
    }
    return std::move(*r.pslg);
}

// ---------------------------------------------------------------------------
// Rectangles, the shape every fixture below is assembled from
// ---------------------------------------------------------------------------

// Counterclockwise, no stored closure. Integer corners, so no kernel can
// disagree about the winding and no fixture fails for a second reason.
[[nodiscard]] inline std::vector<Point2> ccw_rect(double x0, double y0, double x1, double y1) {
    return {Point2{x0, y0}, Point2{x1, y0}, Point2{x1, y1}, Point2{x0, y1}};
}

// The same rectangle wound clockwise: a legal Hole.
[[nodiscard]] inline std::vector<Point2> cw_rect(double x0, double y0, double x1, double y1) {
    return flipped(points(ccw_rect(x0, y0, x1, y1)));
}

// ---------------------------------------------------------------------------
// The degeneracy table, one fixture per row
// ---------------------------------------------------------------------------
//
// The expected status and, for the Ok rows, the exact interior triangle count
// live in the tests rather than here. A fixture that carried its own expected
// answer would let a reader check the two against each other instead of against
// the header, which is the one comparison that means anything.

// The minimal valid input: one Outer triangle.
[[nodiscard]] inline PslgBuilder single_triangle_domain() {
    PslgBuilder b;
    b.add_chain(points(std::vector<Point2>{Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{0.0, 4.0}}),
                ChainRole::Outer);
    return b;
}

// A square, nothing else. Two triangles by Euler: 2*4 - 4 - 2.
[[nodiscard]] inline PslgBuilder square_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    return b;
}

// A square and one disjoint square hole. 2*8 - 8 - 2 + 2 = 8.
[[nodiscard]] inline PslgBuilder disjoint_hole_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(cw_rect(3.0, 3.0, 6.0, 6.0)), ChainRole::Hole);
    return b;
}

// The goal line: an outer ring, two disjoint holes and two breaklines, none of
// them touching. This is the fixture the whole increment exists to produce a
// mesh for, and it is also the multi-ring multi-hole ASan fixture -- with a
// single ring a dangling span usually still points at live memory and passes.
[[nodiscard]] inline PslgBuilder polygon_with_holes_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 20.0, 20.0)), ChainRole::Outer);
    b.add_chain(points(cw_rect(3.0, 3.0, 7.0, 7.0)), ChainRole::Hole);
    b.add_chain(points(cw_rect(13.0, 13.0, 17.0, 17.0)), ChainRole::Hole);
    b.add_chain(points(std::vector<Point2>{Point2{2.0, 10.0}, Point2{10.0, 11.0},
                                           Point2{18.0, 10.0}}),
                ChainRole::Breakline, /*properties=*/kRiver);
    b.add_chain(points(std::vector<Point2>{Point2{9.0, 15.0}, Point2{11.0, 12.0}}),
                ChainRole::Breakline);
    return b;
}

// A hole meeting its outer ring at ONE SHARED INDEX -- the index-taking
// add_chain is the encoding that exists for exactly this. Measured: 5 interior
// triangles, and the rings are no longer disjoint, so Euler predicts 6 and is
// wrong. See the comment on the test.
[[nodiscard]] inline PslgBuilder corner_touching_hole_domain() {
    PslgBuilder b;
    const std::vector<Point2> outer = ccw_rect(0.0, 0.0, 10.0, 10.0);
    b.add_chain(points(outer), ChainRole::Outer);
    const std::uint32_t first = b.append_vertices(
        points(std::vector<Point2>{Point2{1.0, 3.0}, Point2{3.0, 1.0}}));
    const std::vector<std::uint32_t> hole = {0u, first, first + 1u};
    b.add_chain(indices(hole), ChainRole::Hole);
    return b;
}

// The same touch spelled with a SECOND COINCIDENT VERTEX instead of a shared
// index. How you spell the touch is what decides whether it works.
[[nodiscard]] inline PslgBuilder coincident_touch_hole_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(std::vector<Point2>{Point2{0.0, 0.0}, Point2{1.0, 3.0}, Point2{3.0, 1.0}}),
                ChainRole::Hole);
    return b;
}

// A hole sharing a whole EDGE with its outer ring: a notch in the boundary
// spelled as a hole.
[[nodiscard]] inline PslgBuilder hole_sharing_edge_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    const std::uint32_t apex = b.append_vertices(points(std::vector<Point2>{Point2{5.0, 3.0}}));
    const std::vector<std::uint32_t> hole = {0u, apex, 1u};
    b.add_chain(indices(hole), ChainRole::Hole);
    return b;
}

// A vertex lying exactly ON a constraint edge it is not part of: the T-junction
// the noder owns. Unreferenced, which a valid Pslg permits.
[[nodiscard]] inline PslgBuilder foreign_vertex_on_constraint_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.append_vertices(points(std::vector<Point2>{Point2{5.0, 0.0}}));
    return b;
}

// Collinear consecutive vertices WITHIN a ring -- accepted, and easily confused
// with the row above. Redundant ring vertices are fine; a foreign vertex on a
// constraint is not.
[[nodiscard]] inline PslgBuilder collinear_ring_vertices_domain() {
    PslgBuilder b;
    b.add_chain(points(std::vector<Point2>{Point2{0.0, 0.0}, Point2{2.0, 0.0}, Point2{4.0, 0.0},
                                           Point2{4.0, 4.0}, Point2{0.0, 4.0}}),
                ChainRole::Outer);
    return b;
}

// A repeated consecutive index in a ring: {0, 0, 1, 2, 3}. A valid Pslg
// (increment 2 degeneracy 4 accepts it) that detria refuses.
[[nodiscard]] inline PslgBuilder repeated_index_ring_domain() {
    PslgBuilder b;
    b.append_vertices(points(ccw_rect(0.0, 0.0, 10.0, 10.0)));
    const std::vector<std::uint32_t> ring = {0u, 0u, 1u, 2u, 3u};
    b.add_chain(indices(ring), ChainRole::Outer);
    return b;
}

// Two vertices with identical coordinates, neither referenced by any chain.
// detria's duplicate scan runs over the whole point array, so this fails --
// and this is the most likely failure on real data.
[[nodiscard]] inline PslgBuilder duplicate_unreferenced_vertex_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.append_vertices(points(std::vector<Point2>{Point2{1.0, 1.0}, Point2{1.0, 1.0}}));
    return b;
}

// Two crossing breaklines: the un-noded input increment 5 owns.
[[nodiscard]] inline PslgBuilder crossing_breaklines_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(std::vector<Point2>{Point2{2.0, 2.0}, Point2{8.0, 8.0}}),
                ChainRole::Breakline);
    b.add_chain(points(std::vector<Point2>{Point2{2.0, 8.0}, Point2{8.0, 2.0}}),
                ChainRole::Breakline);
    return b;
}

// A hole lying outside every outline. The nesting check increment 3 deferred.
[[nodiscard]] inline PslgBuilder hole_outside_outline_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 4.0, 4.0)), ChainRole::Outer);
    b.add_chain(points(cw_rect(10.0, 10.0, 14.0, 14.0)), ChainRole::Hole);
    return b;
}

// A hole that CONTAINS the outer ring. A different detria arm from the row
// above and a different enumerator, mapping to the same status -- which is the
// point of grouping the mapping by what the caller should do.
[[nodiscard]] inline PslgBuilder hole_containing_outline_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 4.0, 4.0)), ChainRole::Outer);
    b.add_chain(points(cw_rect(-5.0, -5.0, 15.0, 15.0)), ChainRole::Hole);
    return b;
}

// An island in a lake: an outline inside a hole inside an outline. detria
// computes nesting itself, which is why increment 3 stores no parent map.
[[nodiscard]] inline PslgBuilder island_in_lake_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 12.0, 12.0)), ChainRole::Outer);
    b.add_chain(points(cw_rect(2.0, 2.0, 10.0, 10.0)), ChainRole::Hole);
    b.add_chain(points(std::vector<Point2>{Point2{5.0, 5.0}, Point2{7.0, 5.0}, Point2{6.0, 7.0}}),
                ChainRole::Outer);
    return b;
}

// A vertex outside the outer ring: kept in the vertex array for index
// identity, referenced by no interior triangle.
[[nodiscard]] inline PslgBuilder outside_vertex_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 4.0, 4.0)), ChainRole::Outer);
    b.append_vertices(points(std::vector<Point2>{Point2{9.0, 9.0}}));
    return b;
}

// A sliver outer ring: 4 units long and 1e-13 tall. There is no area threshold
// anywhere in this project, so this is an ordinary one-triangle domain.
[[nodiscard]] inline PslgBuilder sliver_domain() {
    PslgBuilder b;
    b.add_chain(
        points(std::vector<Point2>{Point2{0.0, 0.0}, Point2{4.0, 0.0}, Point2{4.0, 1e-13}}),
        ChainRole::Outer);
    return b;
}

// A closed-loop breakline inside a square -- a contour, or a ring road. Fed
// edge by edge it must leave the interior INTACT; auto-detected as a hole it
// would carve a void out of the domain. That difference is the test.
//
// THE LOOP CLOSES ON A SHARED INDEX, not on a repeated point, and that is the
// same lesson the corner-touching hole teaches: pslg_cases' closed_polyline()
// spells the same shape with a coincident fifth vertex, which is a perfectly
// valid Pslg and a DuplicatePointsFound failure at the backend. Measured, not
// assumed -- it is what this fixture did before it was corrected.
[[nodiscard]] inline PslgBuilder closed_loop_breakline_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    const std::uint32_t first = b.append_vertices(
        points(std::vector<Point2>{Point2{3.0, 3.0}, Point2{7.0, 3.0}, Point2{7.0, 7.0},
                                   Point2{3.0, 7.0}}));
    const std::vector<std::uint32_t> loop = {first, first + 1u, first + 2u, first + 3u, first};
    b.add_chain(indices(loop), ChainRole::Breakline);
    return b;
}

// The same ring road spelled with a coincident closing VERTEX instead of a
// shared index. A valid Pslg, and a hard backend failure -- the row of the
// degeneracy table that says duplicate coordinates are not pre-checked by us.
[[nodiscard]] inline PslgBuilder coincident_closed_breakline_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(std::vector<Point2>{Point2{3.0, 3.0}, Point2{7.0, 3.0}, Point2{7.0, 7.0},
                                           Point2{3.0, 7.0}, Point2{3.0, 3.0}}),
                ChainRole::Breakline);
    return b;
}

// A breakline lying entirely OUTSIDE the outer ring, and a second one entirely
// INSIDE a hole. Both are legal in a Pslg -- it promises no disjointness and no
// nesting -- both triangulate Ok, and neither appears in any interior triangle.
// This is the input for which "every constraint edge survives into the mesh" is
// false, which is why the property says so and this fixture pins it.
[[nodiscard]] inline PslgBuilder out_of_domain_breaklines_domain() {
    PslgBuilder b;
    b.add_chain(points(ccw_rect(0.0, 0.0, 10.0, 10.0)), ChainRole::Outer);
    b.add_chain(points(cw_rect(3.0, 3.0, 7.0, 7.0)), ChainRole::Hole);
    b.add_chain(points(std::vector<Point2>{Point2{20.0, 20.0}, Point2{25.0, 25.0}}),
                ChainRole::Breakline);
    b.add_chain(points(std::vector<Point2>{Point2{4.0, 4.0}, Point2{6.0, 6.0}}),
                ChainRole::Breakline);
    return b;
}

// An L-shaped outer ring, for the auto-close characterisation. Non-convex on
// purpose: its closing edge (0,6)->(0,0) is a wall of the domain, so a version
// of detria that stopped closing an open polyline could not produce this answer
// by accident. 2*6 - 6 - 2 = 4 triangles.
[[nodiscard]] inline PslgBuilder l_shaped_domain() {
    PslgBuilder b;
    b.add_chain(points(std::vector<Point2>{Point2{0.0, 0.0}, Point2{6.0, 0.0}, Point2{6.0, 2.0},
                                           Point2{2.0, 2.0}, Point2{2.0, 6.0}, Point2{0.0, 6.0}}),
                ChainRole::Outer);
    return b;
}

// ---------------------------------------------------------------------------
// The generated family
// ---------------------------------------------------------------------------

// The structure a fixture DECLARES, which is where the Euler check takes its
// numbers from. Reading them back off the mesh would make the check a tautology
// and reading n off vertices() would make it wrong -- a vertex inside a hole or
// outside the outer ring is in the array and in no triangle.
//
// rings_disjoint is not decoration either: Euler's formula holds for a domain
// that is a disk with h holes, and a hole touching its outer ring at a shared
// vertex is not one. The corner-touching fixture predicts 6 against a measured
// 5, so it carries an explicit count instead of coming through here.
struct DomainShape {
    std::size_t referenced_vertices{};  // n
    std::size_t boundary_vertices{};    // b, on any ring, outer or hole
    std::size_t holes{};                // h
};

[[nodiscard]] constexpr std::size_t euler_triangles(const DomainShape& s) noexcept {
    return 2 * s.referenced_vertices - s.boundary_vertices - 2 + 2 * s.holes;
}

struct GeneratedDomain {
    PslgBuilder builder;
    DomainShape shape;
};

// A rectangle-with-holes family laid out on a grid of 10x10 cells, one feature
// per cell, every feature inset 2 units from its cell. Nothing touches anything
// -- which is the precondition for a mesh existing at all before increment 5's
// noder -- and every vertex is referenced and strictly inside the domain, which
// is what lets DomainShape be filled in from the construction.
//
// `resolution` subdivides each side of the outer rectangle, so a generated
// outer ring carries long collinear runs. That is a degeneracy the table above
// says is accepted, and putting it in the generator is what keeps it accepted
// under every later edit rather than only in its one named fixture.
[[nodiscard]] inline GeneratedDomain generated_domain(std::mt19937_64& rng) {
    constexpr double cell = 10.0;
    std::uniform_int_distribution<std::size_t> dim{1, 3};
    std::uniform_int_distribution<std::size_t> res{1, 3};
    std::uniform_int_distribution<int> feature{0, 2};  // 0 none, 1 hole, 2 breakline

    const std::size_t cols = dim(rng);
    const std::size_t rows = dim(rng);
    const std::size_t resolution = res(rng);
    const double w = cell * static_cast<double>(cols);
    const double h = cell * static_cast<double>(rows);

    // The outer ring, counterclockwise, with `resolution` vertices per side.
    std::vector<Point2> outer;
    const auto side = [&](Point2 from, Point2 to) {
        for (std::size_t k = 0; k < resolution; ++k) {
            const double t = static_cast<double>(k) / static_cast<double>(resolution);
            outer.push_back(Point2{from.x + (to.x - from.x) * t, from.y + (to.y - from.y) * t});
        }
    };
    side(Point2{0.0, 0.0}, Point2{w, 0.0});
    side(Point2{w, 0.0}, Point2{w, h});
    side(Point2{w, h}, Point2{0.0, h});
    side(Point2{0.0, h}, Point2{0.0, 0.0});

    PslgBuilder b;
    b.add_chain(points(outer), ChainRole::Outer);

    DomainShape shape{outer.size(), outer.size(), 0};

    for (std::size_t j = 0; j < rows; ++j) {
        for (std::size_t i = 0; i < cols; ++i) {
            const double x = cell * static_cast<double>(i);
            const double y = cell * static_cast<double>(j);
            switch (feature(rng)) {
                case 1: {
                    b.add_chain(points(cw_rect(x + 2.0, y + 2.0, x + 8.0, y + 8.0)),
                                ChainRole::Hole);
                    shape.referenced_vertices += 4;
                    shape.boundary_vertices += 4;
                    ++shape.holes;
                    break;
                }
                case 2: {
                    // Three vertices, strictly interior to the cell and to the
                    // domain, so every breakline vertex is an interior point of
                    // the triangulation and contributes 2 triangles.
                    b.add_chain(points(std::vector<Point2>{Point2{x + 2.0, y + 3.0},
                                                           Point2{x + 5.0, y + 6.0},
                                                           Point2{x + 8.0, y + 3.0}}),
                                ChainRole::Breakline);
                    shape.referenced_vertices += 3;
                    break;
                }
                default: break;
            }
        }
    }
    return GeneratedDomain{std::move(b), shape};
}

}  // namespace terrain::test
