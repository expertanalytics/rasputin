// Unit tests for terrain/core/indexed_mesh.hpp: the flat mesh the CDT returns
// and `mesh` will later build its ternary tree from.
//
// NOT invariant-critical, and no mutation round is spent here. This file's
// failure mode is a typo in an accessor, not a topology decision -- the
// increment file says so and the split is deliberate: the spend goes to the
// backend, the mask and the property suite.
//
// The constructor is PUBLIC and that is a deliberate break from Pslg. Pslg's
// is private because it has exactly one legitimate producer and holding one is
// a proof that validation ran; IndexedMesh2 has plural producers on purpose --
// DetriaBackend, the seam suite's FakeCdtBackend, mesh's flatten step, and an
// eventual ingress handing Python an existing TIN. The decisive one is the
// second: a fake backend that cannot construct a mesh cannot prove the seam is
// real. So every fixture here hand-builds, which is the type's contract rather
// than a shortcut around it.
//
// Guarantees 1-3 (parallel arrays, in-range indices, masks < 8) are DEBUG
// ASSERTS in the constructor, not a status channel and not a release check, so
// no test here constructs a violating mesh: under the asan+ubsan Debug job that
// is a deliberate abort, and under Release it is undefined behaviour. What is
// testable is that a well-formed mesh reports what it was given.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>

#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::Segment2;
using terrain::TriangleIndices;
using terrain::kEdgeMask01;
using terrain::kEdgeMask12;
using terrain::kEdgeMask20;

namespace {

// The unit square as two triangles, sharing the diagonal (1, 3). Triangle 0
// carries exactly one constrained edge -- slot 0, v0->v1 -- and that single bit
// is what tells our convention apart from the opposite-vertex one. See the mask
// suite for the full argument.
[[nodiscard]] IndexedMesh2 two_triangle_mesh() {
    std::vector<Point2> vertices = {Point2{0.0, 0.0}, Point2{1.0, 0.0}, Point2{1.0, 1.0},
                                    Point2{0.0, 1.0}};
    std::vector<TriangleIndices> triangles = {TriangleIndices{0, 1, 3}, TriangleIndices{1, 2, 3}};
    std::vector<std::uint8_t> masks = {kEdgeMask01, kEdgeMask01 | kEdgeMask12};
    return IndexedMesh2{std::move(vertices), std::move(triangles), std::move(masks)};
}

}  // namespace

// ---------------------------------------------------------------------------
// The accessors
// ---------------------------------------------------------------------------

TEST_CASE("a default-constructed mesh is empty", "[cdt][mesh]") {
    const IndexedMesh2 mesh;

    CHECK(mesh.empty());
    CHECK(mesh.triangle_count() == 0);
    CHECK(mesh.vertices().empty());
    CHECK(mesh.triangles().empty());
    CHECK(mesh.constrained_edges().empty());
}

TEST_CASE("a mesh reports the buffers it was constructed from", "[cdt][mesh]") {
    const IndexedMesh2 mesh = two_triangle_mesh();

    REQUIRE(mesh.vertices().size() == 4);
    REQUIRE(mesh.triangle_count() == 2);
    CHECK_FALSE(mesh.empty());
    // Guarantee 1, as an observation rather than as the debug assert: the two
    // parallel arrays are index-aligned, which is what lets a refinement pass
    // stream the mask without touching coordinates.
    CHECK(mesh.triangles().size() == mesh.constrained_edges().size());
    CHECK(mesh.triangles()[0] == TriangleIndices{0, 1, 3});
    CHECK(mesh.vertices()[2] == Point2{1.0, 1.0});
}

// empty() is triangle_count() == 0 and NOT vertices().empty(). The distinction
// is the whole "forgot addOutline" failure mode: a point set with no outline
// triangulates successfully and yields zero interior triangles over a full
// vertex array.
TEST_CASE("empty() asks about triangles, not vertices", "[cdt][mesh]") {
    const IndexedMesh2 mesh{std::vector<Point2>{Point2{0.0, 0.0}, Point2{1.0, 0.0}}, {}, {}};

    CHECK(mesh.vertices().size() == 2);
    CHECK(mesh.empty());
    CHECK(mesh.triangle_count() == 0);
}

// ---------------------------------------------------------------------------
// edge() and is_constrained(), which must not drift apart
// ---------------------------------------------------------------------------

// Bit e is the edge (v[e], v[(e+1) % 3]) and edge(t, e) returns that same
// segment, directed the same way. Asserting the three slots by number is what
// keeps the accessor and the mask convention married; asserting them as a set
// would pass under any rotation.
TEST_CASE("edge(t, e) is the segment from v[e] to v[(e+1) % 3]", "[cdt][mesh]") {
    const IndexedMesh2 mesh = two_triangle_mesh();
    const std::span<const Point2> v = mesh.vertices();

    CHECK(mesh.edge(0, 0) == Segment2{v[0], v[1]});
    CHECK(mesh.edge(0, 1) == Segment2{v[1], v[3]});
    CHECK(mesh.edge(0, 2) == Segment2{v[3], v[0]});
}

TEST_CASE("is_constrained reads bit e of the triangle's mask", "[cdt][mesh]") {
    const IndexedMesh2 mesh = two_triangle_mesh();

    CHECK(mesh.is_constrained(0, 0));
    CHECK_FALSE(mesh.is_constrained(0, 1));
    CHECK_FALSE(mesh.is_constrained(0, 2));

    CHECK(mesh.is_constrained(1, 0));
    CHECK(mesh.is_constrained(1, 1));
    CHECK_FALSE(mesh.is_constrained(1, 2));
}

TEST_CASE("the mask constants are the three single bits", "[cdt][mesh]") {
    STATIC_REQUIRE(kEdgeMask01 == 1u);
    STATIC_REQUIRE(kEdgeMask12 == 2u);
    STATIC_REQUIRE(kEdgeMask20 == 4u);
    STATIC_REQUIRE((kEdgeMask01 | kEdgeMask12 | kEdgeMask20) == 7u);
}

// ---------------------------------------------------------------------------
// Value semantics
// ---------------------------------------------------------------------------

// Guarantee 5: immutable after construction, and a const IndexedMesh2 is safe
// for concurrent read from any number of threads -- which refinement requires.
// The type-level half of that is asserted here; there is no mutator to test and
// that absence is the point.
TEST_CASE("IndexedMesh2 is a copyable, movable value type", "[cdt][mesh]") {
    STATIC_REQUIRE(std::is_copy_constructible_v<IndexedMesh2>);
    STATIC_REQUIRE(std::is_move_constructible_v<IndexedMesh2>);
    STATIC_REQUIRE(std::is_nothrow_move_constructible_v<IndexedMesh2>);
    STATIC_REQUIRE(std::is_default_constructible_v<IndexedMesh2>);
}

// A copy owns its own buffers. IndexedMesh2 stores no span into itself, for the
// same reason Pslg stores none: a cached view member would make a copied mesh
// point at the original's storage, and that bug survives every test that never
// copies.
TEST_CASE("a copied mesh shares no storage with its source", "[cdt][mesh]") {
    const IndexedMesh2 original = two_triangle_mesh();
    const IndexedMesh2 copy = original;

    REQUIRE(copy.triangle_count() == original.triangle_count());
    CHECK(copy.vertices().data() != original.vertices().data());
    CHECK(copy.triangles().data() != original.triangles().data());
    CHECK(copy.constrained_edges().data() != original.constrained_edges().data());
    for (std::size_t t = 0; t < original.triangle_count(); ++t) {
        CHECK(copy.triangles()[t] == original.triangles()[t]);
        CHECK(copy.constrained_edges()[t] == original.constrained_edges()[t]);
    }
}

TEST_CASE("a moved-from mesh hands over its buffers", "[cdt][mesh]") {
    IndexedMesh2 source = two_triangle_mesh();
    const Point2* const vertex_storage = source.vertices().data();

    const IndexedMesh2 moved = std::move(source);

    CHECK(moved.vertices().data() == vertex_storage);
    CHECK(moved.triangle_count() == 2);
}
