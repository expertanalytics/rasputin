// Increment 23b (docs/increments/23-basin-scale.md, "The seam protocol" step 4,
// K2, and "Tests @tester can write red", FE4): the frozen mask on LatticeMesh,
// and the quality pass's skip on a frozen edge. FE4 is invariant-critical.
//
// Interface, as the PR table names it ("lattice_mesh.hpp: the frozen mask,
// is_frozen, the assertion in split_edge"; "quality.hpp: skip on frozen,
// skipped_frozen"), with what this suite PINS where the design is silent:
//
//   void          LatticeMesh::set_frozen_mask(std::uint32_t) noexcept
//   std::uint32_t LatticeMesh::frozen_mask() const noexcept      0 after build()
//   bool          LatticeMesh::is_frozen(std::size_t t, unsigned e) const noexcept
//                 (mask(t, e) & frozen_mask()) != 0: "an edge is frozen when
//                 its mask meets it"
//   std::size_t   QualityOutcome::skipped_frozen
//                 a bad triangle whose snapped node lies exactly on a frozen
//                 edge of the triangle the walk found: counted, nothing inserted
//
// The assertion in split_edge is not tested (FE6: no death-test harness, and
// Release defines NDEBUG).

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/mesh/quality.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::TriangleIndices;
using terrain::mesh::improve;
using terrain::mesh::lattice_frame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::legalise_all;
using terrain::mesh::MeshVertex;
using terrain::mesh::QualityOptions;
using terrain::mesh::QualityOutcome;
using terrain::pred::DefaultKernel;

namespace {

using Constraints = std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t>;

std::pair<std::uint32_t, std::uint32_t> undirected(std::uint32_t a, std::uint32_t b) { return std::minmax(a, b); }

LatticeMesh build(const std::vector<MeshVertex>& v, const std::vector<TriangleIndices>& tris, const Constraints& c) {
    std::vector<std::uint8_t> bits(tris.size(), 0);
    std::vector<std::array<std::uint32_t, 3>> masks(tris.size(), {0, 0, 0});
    for (std::size_t t = 0; t < tris.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto it = c.find(undirected(tris[t][k], tris[t][(k + 1) % 3])); it != c.end()) {
                bits[t] |= static_cast<std::uint8_t>(1u << k);
                masks[t][k] = it->second;
            }
    auto m = LatticeMesh::build(v, tris, bits, masks);
    REQUIRE(m.has_value());
    return std::move(*m);
}

// The square [0, 8]^2 cut by a vertical seam from (4, 0) to (4, 8):
// 0 (0,0), 1 (4,0), 2 (8,0), 3 (8,8), 4 (4,8), 5 (0,8). Outline mask 1, the
// seam 4 | `extra` (a seam along a feature carries both bits).
LatticeMesh seam_square(std::uint32_t extra = 0) {
    const std::vector<MeshVertex> v{{0, 0}, {4, 0}, {8, 0}, {8, 8}, {4, 8}, {0, 8}};
    const std::vector<TriangleIndices> t{{0, 5, 1}, {1, 5, 4}, {1, 4, 3}, {1, 3, 2}};
    const Constraints c{{undirected(0, 1), 1}, {undirected(1, 2), 1}, {undirected(2, 3), 1}, {undirected(3, 4), 1},
                        {undirected(4, 5), 1}, {undirected(5, 0), 1}, {undirected(1, 4), 4u | extra}};
    return build(v, t, c);
}

// Every (triangle, edge) whose two vertices are a and b, in either order.
std::vector<std::pair<std::uint32_t, unsigned>> sides_of(const LatticeMesh& m, MeshVertex a, MeshVertex b) {
    std::vector<std::pair<std::uint32_t, unsigned>> out;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const MeshVertex p = m.corner(t, k), q = m.corner(t, (k + 1) % 3);
            if ((p == a && q == b) || (p == b && q == a)) out.push_back({t, k});
        }
    return out;
}

std::size_t frozen_sides(const LatticeMesh& m) {
    std::size_t n = 0;
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) n += m.is_frozen(t, k) ? 1 : 0;
    return n;
}

// p strictly inside the open segment (a, b), exactly.
bool on_open(MeshVertex a, MeshVertex b, MeshVertex p) {
    if (p == a || p == b) return false;
    if (DefaultKernel::orient2d(a.frame(), b.frame(), p.frame()) != terrain::pred::Orientation::Collinear) return false;
    return std::min(a.col, b.col) <= p.col && p.col <= std::max(a.col, b.col) && std::min(a.row, b.row) <= p.row
        && p.row <= std::max(a.row, b.row);
}

std::size_t vertices_on(const LatticeMesh& m, MeshVertex a, MeshVertex b) {
    const auto v = m.vertices();
    return static_cast<std::size_t>(std::count_if(v.begin(), v.end(), [&](MeshVertex p) { return on_open(a, b, p); }));
}

QualityOutcome run_quality(LatticeMesh& m, std::size_t rows, std::size_t cols, double theta = 25.0) {
    const auto frame = lattice_frame(1.0, 1.0, rows, cols);
    legalise_all<DefaultKernel>(m, frame, [](std::uint32_t) {});
    return improve<DefaultKernel>(m, frame, QualityOptions{.min_angle_deg = theta, .rows = rows, .cols = cols});
}

}  // namespace

// ------------------------------------------------------------- the mask

TEST_CASE("FM1: a built mesh has frozen mask 0 and no frozen edge", "[mesh][frozen]") {
    const LatticeMesh m = seam_square();
    REQUIRE(m.frozen_mask() == 0u);
    REQUIRE(frozen_sides(m) == 0);
}

TEST_CASE("FM1: an edge is frozen exactly when its mask meets the frozen mask", "[mesh][frozen]") {
    LatticeMesh m = seam_square();
    m.set_frozen_mask(4u);
    REQUIRE(m.frozen_mask() == 4u);
    const auto seam = sides_of(m, MeshVertex{4, 0}, MeshVertex{4, 8});
    REQUIRE(seam.size() == 2);  // interior: both triangles see it
    for (const auto& [t, k] : seam) REQUIRE(m.is_frozen(t, k));
    REQUIRE(frozen_sides(m) == 2);  // and nothing else: the outline (1) and the diagonals (0) are not

    m.set_frozen_mask(1u);  // the outline instead
    for (const auto& [t, k] : seam) REQUIRE_FALSE(m.is_frozen(t, k));
    REQUIRE(frozen_sides(m) == 6);

    m.set_frozen_mask(0u);
    REQUIRE(frozen_sides(m) == 0);
}

TEST_CASE("FM1: a seam along a feature edge carries both bits and is frozen under the seam's", "[mesh][frozen]") {
    // Degeneracy policy, "A seam along a feature edge": one edge with both
    // bits; the seam's rules apply.
    LatticeMesh m = seam_square(2u);
    const auto seam = sides_of(m, MeshVertex{4, 0}, MeshVertex{4, 8});
    m.set_frozen_mask(4u);
    for (const auto& [t, k] : seam) {
        REQUIRE(m.mask(t, k) == 6u);
        REQUIRE(m.is_frozen(t, k));
    }
    m.set_frozen_mask(8u);  // a mask that meets neither bit
    for (const auto& [t, k] : seam) REQUIRE_FALSE(m.is_frozen(t, k));
}

TEST_CASE("FM2: the frozen edge stays frozen through splits and flips around it", "[mesh][frozen]") {
    // Work on both sides of the seam, never on it: a split inside each piece,
    // a split of an outline edge and a full legalisation. The seam keeps its
    // two vertices, its mask and its frozen state in whichever slots now hold
    // it, and no new edge is frozen (new interior edges carry no mask).
    LatticeMesh m = seam_square();
    m.set_frozen_mask(4u);
    const auto frame = lattice_frame(1.0, 1.0, 9, 9);
    auto find = [&](MeshVertex p) {
        for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
            bool in = true;
            for (unsigned k = 0; k < 3; ++k)
                in = in && DefaultKernel::orient2d(m.corner(t, k).frame(), m.corner(t, (k + 1) % 3).frame(), p.frame())
                               == terrain::pred::Orientation::CounterClockwise;
            if (in) return t;
        }
        FAIL("no triangle holds the point");
        return std::uint32_t{0};
    };
    m.split_inside(find(MeshVertex{1, 5}), MeshVertex{1, 5});
    m.split_inside(find(MeshVertex{6, 2}), MeshVertex{6, 2});
    // The outline's top edge (0,0)-(4,0) at (2, 0).
    const auto top = sides_of(m, MeshVertex{0, 0}, MeshVertex{4, 0});
    REQUIRE(top.size() == 1);
    m.split_edge(top[0].first, top[0].second, MeshVertex{2, 0});
    legalise_all<DefaultKernel>(m, frame, [](std::uint32_t) {});

    const auto seam = sides_of(m, MeshVertex{4, 0}, MeshVertex{4, 8});
    REQUIRE(seam.size() == 2);
    for (const auto& [t, k] : seam) {
        REQUIRE(m.is_constrained(t, k));
        REQUIRE(m.mask(t, k) == 4u);
        REQUIRE(m.is_frozen(t, k));
    }
    REQUIRE(frozen_sides(m) == 2);
    REQUIRE(m.frozen_mask() == 4u);
}

// ------------------------------------------------------------------------ FE4

TEST_CASE("FE4: the quality pass never splits a frozen boundary edge and counts the skip", "[mesh][quality][frozen]") {
    // test_mesh_quality's Q3: a right triangle's circumcentre is its
    // hypotenuse's midpoint, node (10, 1), exactly on that edge. Unfrozen the
    // pass splits the hypotenuse there (the control, so this can fail);
    // frozen it inserts nothing on it and counts skipped_frozen.
    const std::vector<MeshVertex> v{{0, 1}, {2, 7}, {20, 1}};
    const Constraints c{{undirected(0, 1), 1}, {undirected(1, 2), 1}, {undirected(2, 0), 4}};

    LatticeMesh control = build(v, {{0, 1, 2}}, c);
    const QualityOutcome q0 = run_quality(control, 8, 21);
    REQUIRE(q0.skipped_frozen == 0);
    REQUIRE(vertices_on(control, v[0], v[2]) >= 1);

    LatticeMesh m = build(v, {{0, 1, 2}}, c);
    m.set_frozen_mask(4u);
    const QualityOutcome q = run_quality(m, 8, 21);
    REQUIRE(q.skipped_frozen >= 1);
    REQUIRE(vertices_on(m, v[0], v[2]) == 0);
    const auto mv = m.vertices();
    REQUIRE(std::find(mv.begin(), mv.end(), MeshVertex{10, 1}) == mv.end());
    const auto hyp = sides_of(m, v[0], v[2]);
    REQUIRE(hyp.size() == 1);
    REQUIRE(m.is_frozen(hyp[0].first, hyp[0].second));
}

TEST_CASE("FE4: the quality pass never splits a frozen interior edge", "[mesh][quality][frozen]") {
    // Q3 mirrored across the hypotenuse (0, 7)-(20, 7): two right triangles,
    // both bad at 25 degrees, both with circumcentre node (10, 7) on the shared
    // edge, a seam. Unfrozen it is split 2 -> 4; frozen, neither triangle
    // splits it and each skip is counted.
    const std::vector<MeshVertex> v{{0, 7}, {2, 13}, {20, 7}, {18, 1}};
    const std::vector<TriangleIndices> t{{0, 1, 2}, {2, 3, 0}};
    const Constraints c{{undirected(0, 1), 1}, {undirected(1, 2), 1}, {undirected(2, 3), 1},
                        {undirected(3, 0), 1}, {undirected(0, 2), 4}};

    LatticeMesh control = build(v, t, c);
    run_quality(control, 14, 21);
    REQUIRE(vertices_on(control, v[0], v[2]) >= 1);

    LatticeMesh m = build(v, t, c);
    m.set_frozen_mask(4u);
    const QualityOutcome q = run_quality(m, 14, 21);
    REQUIRE(q.skipped_frozen >= 2);
    REQUIRE(vertices_on(m, v[0], v[2]) == 0);
    REQUIRE(sides_of(m, v[0], v[2]).size() == 2);
}

TEST_CASE("FE4: a frozen mask that meets no edge changes nothing in the quality pass", "[mesh][quality][frozen]") {
    // K1 at the pass: bit for bit the mesh mask 0 gives.
    const std::vector<MeshVertex> v{{0, 7}, {2, 13}, {20, 7}, {18, 1}};
    const std::vector<TriangleIndices> t{{0, 1, 2}, {2, 3, 0}};
    const Constraints c{{undirected(0, 1), 1}, {undirected(1, 2), 1}, {undirected(2, 3), 1},
                        {undirected(3, 0), 1}, {undirected(0, 2), 4}};
    LatticeMesh a = build(v, t, c), b = build(v, t, c);
    b.set_frozen_mask(1u << 31);
    const QualityOutcome qa = run_quality(a, 14, 21), qb = run_quality(b, 14, 21);
    REQUIRE(qa.inserted == qb.inserted);
    REQUIRE(qb.skipped_frozen == 0);
    REQUIRE(std::vector<MeshVertex>(a.vertices().begin(), a.vertices().end())
            == std::vector<MeshVertex>(b.vertices().begin(), b.vertices().end()));
    REQUIRE(std::vector<TriangleIndices>(a.triangles().begin(), a.triangles().end())
            == std::vector<TriangleIndices>(b.triangles().begin(), b.triangles().end()));
}
