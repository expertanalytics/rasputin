// Increment 20's fix (docs/increments/20-start-quality.md, "Fix: the
// start-quality pass skips NoData nodes", Q-V1): improve never inserts a node
// its validity callable rejects. A bad triangle whose snapped node is invalid
// stays as it is, and is counted in QualityOutcome::skipped_void.
//
// Interface, as the design gives it:
//
//   struct AllNodesValid { bool operator()(const LatticeVertex&) const; };  // true
//   template <class K, class Valid = AllNodesValid>
//   QualityOutcome improve(LatticeMesh&, const LatticeFrame&, const QualityOptions&,
//                          const Valid& valid = {});
//
// The callable is asked after the node is snapped and found inside the
// rectangle, before the visibility walk: the floor and outside skips come
// first, and a rejected node is never walked to (so it is a void skip, not a
// vertex or blocked one).
//
// A target of its own, not a section of test_mesh_quality.cpp: the
// field-by-field guard of the default call is there, and stays buildable
// whatever this file needs.
//
// Delaunay uses DefaultKernel's incircle in the frame the pass legalises in
// (LatticeFrame::at), the producer's predicate but not its records.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/mesh/quality.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <utility>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::AllNodesValid;
using terrain::mesh::improve;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::legalise_all;
using terrain::mesh::MeshVertex;
using terrain::mesh::QualityOptions;
using terrain::mesh::QualityOutcome;
using terrain::pred::DefaultKernel;
using terrain::pred::Incircle;

namespace {

constexpr double kTheta = 25.0;

struct Start {
    std::vector<MeshVertex> vertices;
    std::vector<TriangleIndices> triangles;
    std::vector<std::uint8_t> bits;  // every edge constrained unless set otherwise
    std::size_t rows = 0;
    std::size_t cols = 0;
};

LatticeMesh build(const Start& s) {
    const auto bits = s.bits.empty() ? std::vector<std::uint8_t>(s.triangles.size(), 0b111) : s.bits;
    auto m = LatticeMesh::build(s.vertices, s.triangles, bits,
                                std::vector<std::array<std::uint32_t, 3>>(s.triangles.size(), {1, 1, 1}));
    REQUIRE(m.has_value());
    return std::move(*m);
}

Start single(MeshVertex a, MeshVertex b, MeshVertex c, std::size_t rows, std::size_t cols) {
    return Start{{a, b, c}, {{0, 1, 2}}, {}, rows, cols};
}

// Q6's world-thin triangle (test_mesh_quality.cpp): dx 1, dy 3, 21.8° at its
// apex in the world. Its circumcentre snaps to node (col 15, row 13), inside
// it; the pass inserts that node and nothing else (the default-call guard
// pins it).
constexpr LatticeFrame kThinFrame{1.0, 3.0};
Start thin() { return single(MeshVertex{0, 0}, MeshVertex{15, 26}, MeshVertex{30, 0}, 27, 31); }
constexpr LatticeVertex kThinNode{13, 15};  // {row, col}

// The quarter-circle analogue of test_mesh_quality.cpp, a fan of slivers.
Start quarter_circle() {
    constexpr std::size_t kArc = 70;
    Start s;
    s.rows = s.cols = 64;
    const double c0 = 62.0, r0 = 62.0, radius = 60.0;
    s.vertices.push_back(MeshVertex{c0, r0});
    for (std::size_t i = 0; i <= kArc; ++i) {
        const double t = std::numbers::pi / 2.0 * (1.0 + static_cast<double>(i) / kArc);
        s.vertices.push_back(MeshVertex{c0 + radius * std::cos(t), r0 - radius * std::sin(t)});
    }
    s.vertices[1] = MeshVertex{c0, r0 - radius};
    s.vertices[kArc + 1] = MeshVertex{c0 - radius, r0};
    const auto n = static_cast<std::uint32_t>(s.vertices.size());
    for (std::uint32_t i = 1; i + 1 < n; ++i) {
        s.triangles.push_back({0, i, i + 1});
        // Edge 1 (the chord) always; edge 0 (0 -> i) only on the first, edge 2 only on the last.
        s.bits.push_back(static_cast<std::uint8_t>(0b010 | (i == 1 ? 0b001 : 0) | (i + 2 == n ? 0b100 : 0)));
    }
    return s;
}

QualityOptions options(const Start& s) {
    return QualityOptions{.min_angle_deg = kTheta, .rows = s.rows, .cols = s.cols};
}

std::size_t skips(const QualityOutcome& q) {
    return q.skipped_floor + q.skipped_outside + q.skipped_vertex + q.skipped_blocked + q.walk_bound_hits
         + q.skipped_void;
}

bool has_vertex(const LatticeMesh& m, MeshVertex p) {
    return std::find(m.vertices().begin(), m.vertices().end(), p) != m.vertices().end();
}

std::vector<TriangleIndices> triangles_of(const LatticeMesh& m) {
    return {m.triangles().begin(), m.triangles().end()};
}

std::size_t delaunay_violations(const LatticeMesh& m, const LatticeFrame& f) {
    const auto v = m.vertices();
    std::size_t bad = 0;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto u = m.neighbours(t)[k];
            if (u == kNoNeighbour || m.is_constrained(t, k))
                continue;
            const auto& tt = m.triangles()[t];
            std::uint32_t apex = 0;
            for (const auto x : m.triangles()[u])
                if (x != tt[k] && x != tt[(k + 1) % 3])
                    apex = x;
            if (DefaultKernel::incircle(f.at(v[tt[0]]), f.at(v[tt[1]]), f.at(v[tt[2]]), f.at(v[apex]))
                == Incircle::Inside)
                ++bad;
        }
    return bad;
}

}  // namespace

static_assert(AllNodesValid{}(LatticeVertex{0, 0}));

TEST_CASE("Q-V1: a bad triangle whose node is invalid is skipped and counted", "[mesh][quality][void]") {
    std::vector<LatticeVertex> asked;
    const auto valid = [&](const LatticeVertex& v) {
        asked.push_back(v);
        return !(v == kThinNode);
    };
    LatticeMesh m = build(thin());
    const QualityOutcome q = improve<DefaultKernel>(m, kThinFrame, options(thin()), valid);
    CHECK(q.inserted == 0);
    CHECK(q.skipped_void == 1);
    CHECK(skips(q) == 1);
    CHECK(asked == std::vector<LatticeVertex>{kThinNode});  // asked once, about the snapped node
    CHECK(m.vertices().size() == 3);
    CHECK(triangles_of(m) == std::vector<TriangleIndices>{{0, 1, 2}});  // left as it is
}

TEST_CASE("Q-V1: the control, every node valid, inserts that node", "[mesh][quality][void]") {
    const auto valid = [](const LatticeVertex&) { return true; };
    LatticeMesh m = build(thin());
    const QualityOutcome q = improve<DefaultKernel>(m, kThinFrame, options(thin()), valid);
    CHECK(q.inserted == 1);
    CHECK(q.skipped_void == 0);
    CHECK(has_vertex(m, MeshVertex{kThinNode}));
}

TEST_CASE("Q-V1: an invalid node elsewhere changes nothing", "[mesh][quality][void]") {
    // The callable is about the snapped node, not about some node near it.
    const LatticeVertex neighbour = GENERATE(LatticeVertex{12, 15}, LatticeVertex{13, 14}, LatticeVertex{14, 16});
    CAPTURE(neighbour.row, neighbour.col);
    const auto valid = [&](const LatticeVertex& v) { return !(v == neighbour); };
    LatticeMesh m = build(thin());
    const QualityOutcome q = improve<DefaultKernel>(m, kThinFrame, options(thin()), valid);
    CHECK(q.inserted == 1);
    CHECK(q.skipped_void == 0);
    CHECK(has_vertex(m, MeshVertex{kThinNode}));
}

TEST_CASE("Q-V1: the default call, AllNodesValid and an always-true callable agree", "[mesh][quality][void]") {
    const double dy = GENERATE(10.0, 5.0);
    CAPTURE(dy);
    const LatticeFrame frame{10.0, dy};
    const Start s = quarter_circle();
    const auto legalised = [&] {
        LatticeMesh m = build(s);
        legalise_all<DefaultKernel>(m, frame, [](std::uint32_t) {});
        return m;
    };
    LatticeMesh a = legalised(), b = legalised(), c = legalised();
    const QualityOutcome qa = improve<DefaultKernel>(a, frame, options(s));
    const QualityOutcome qb = improve<DefaultKernel>(b, frame, options(s), AllNodesValid{});
    const QualityOutcome qc = improve<DefaultKernel>(c, frame, options(s), [](const LatticeVertex&) { return true; });
    REQUIRE(qa.inserted > 0);  // the fixture exercises the pass
    for (const QualityOutcome* q : {&qb, &qc}) {
        CHECK(q->inserted == qa.inserted);
        CHECK(q->skipped_floor == qa.skipped_floor);
        CHECK(q->skipped_outside == qa.skipped_outside);
        CHECK(q->skipped_vertex == qa.skipped_vertex);
        CHECK(q->skipped_blocked == qa.skipped_blocked);
        CHECK(q->walk_bound_hits == qa.walk_bound_hits);
        CHECK(q->skipped_void == 0);
    }
    CHECK(qa.skipped_void == 0);
    for (const LatticeMesh* m : {&b, &c}) {
        CHECK(std::vector<MeshVertex>(m->vertices().begin(), m->vertices().end())
              == std::vector<MeshVertex>(a.vertices().begin(), a.vertices().end()));
        CHECK(triangles_of(*m) == triangles_of(a));
    }
}

TEST_CASE("Q-V1: the earlier skips come first; the void check precedes the walk", "[mesh][quality][void]") {
    const auto none = [](const LatticeVertex&) { return false; };
    SECTION("outside the rectangle: an outside skip, never asked") {
        // Q9's fixture: circumcentre at row -99.
        const Start s = single(MeshVertex{0, 0}, MeshVertex{20, 2}, MeshVertex{40, 0}, 3, 41);
        LatticeMesh m = build(s);
        const QualityOutcome q = improve<DefaultKernel>(m, LatticeFrame{10.0, 10.0}, options(s), none);
        CHECK(q.skipped_outside == 1);
        CHECK(q.skipped_void == 0);
        CHECK(skips(q) == 1);
    }
    SECTION("below the floor: a floor skip, never asked") {
        const Start s = single(MeshVertex{2.0, 5.0}, MeshVertex{2.3, 5.1}, MeshVertex{2.6, 5.0}, 11, 11);
        LatticeMesh m = build(s);
        const QualityOutcome q = improve<DefaultKernel>(m, LatticeFrame{10.0, 10.0}, options(s), none);
        CHECK(q.skipped_floor == 1);
        CHECK(q.skipped_void == 0);
        CHECK(skips(q) == 1);
    }
    SECTION("a node that is already a vertex: void, since the walk is not taken") {
        // test_mesh_quality.cpp's vertex-skip fixture, unlegalised: ABC is bad
        // and its circumcentre is exactly D = (10, 34), the far apex across AB.
        Start s{{MeshVertex{0, 10}, MeshVertex{20, 10}, MeshVertex{10, 8}, MeshVertex{10, 34}},
                {{0, 1, 2}, {0, 3, 1}},
                {0b110, 0b011},
                36,
                21};
        LatticeMesh m = build(s);
        const QualityOutcome q = improve<DefaultKernel>(m, LatticeFrame{1.0, 1.0}, options(s), none);
        CHECK(q.skipped_void == 1);
        CHECK(q.skipped_vertex == 0);
        CHECK(skips(q) == 1);
    }
    SECTION("a node on a constrained edge: void, not split") {
        // Q3's right triangle: the circumcentre is node (10, 1) on the hypotenuse.
        const Start s = single(MeshVertex{0, 1}, MeshVertex{2, 7}, MeshVertex{20, 1}, 8, 21);
        LatticeMesh m = build(s);
        const auto valid = [](const LatticeVertex& v) { return !(v == LatticeVertex{1, 10}); };
        const QualityOutcome q = improve<DefaultKernel>(m, LatticeFrame{1.0, 1.0}, options(s), valid);
        CHECK(q.inserted == 0);
        CHECK(q.skipped_void == 1);
        CHECK(triangles_of(m) == std::vector<TriangleIndices>{{0, 1, 2}});
    }
}

TEST_CASE("Q-V1: over a block of invalid nodes no vertex lands in it, and the mesh stays Delaunay",
          "[mesh][quality][void]") {
    // The quarter circle with every node in a band of columns invalid, as a
    // NoData stripe would be: the pass still works outside it, and every
    // vertex it adds is a valid node.
    const double dy = GENERATE(10.0, 5.0);
    CAPTURE(dy);
    const LatticeFrame frame{10.0, dy};
    const Start s = quarter_circle();
    const auto valid = [](const LatticeVertex& v) { return v.col < 20 || v.col > 40; };
    LatticeMesh m = build(s);
    legalise_all<DefaultKernel>(m, frame, [](std::uint32_t) {});
    const QualityOutcome q = improve<DefaultKernel>(m, frame, options(s), valid);
    CHECK(q.inserted > 0);
    CHECK(q.skipped_void > 0);
    CHECK(q.walk_bound_hits == 0);
    CHECK(m.vertices().size() == s.vertices.size() + q.inserted);
    for (std::size_t i = s.vertices.size(); i < m.vertices().size(); ++i) {
        const auto node = m.vertices()[i].as_node();
        REQUIRE(node.has_value());
        CAPTURE(node->row, node->col);
        CHECK(valid(*node));
    }
    CHECK(delaunay_violations(m, frame) == 0);
}
