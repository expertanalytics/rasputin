// Increment 21a (docs/increments/21-parallel-refine.md, end of section 3, and
// "Pinned by the red suite (21a)"): legalise_around with a caller-owned stack,
// so refine's serial phase stops allocating one per insertion. 21a is
// bit-identical, so the contract is only what is observable: on the same mesh
// and seeds, the same flips, the same on_write sequence and the same mesh as
// today's legalise_around. Allocation counts are not pinned.
//
// Interface pinned here (include/terrain/mesh/lawson.hpp):
//
//   using FlipStack = std::vector<std::uint32_t>;
//   template <pred::GeometryKernel K, class OnWrite>
//   std::size_t legalise_around(LatticeMesh&, std::uint32_t q,
//                               std::span<const std::uint32_t> seeds,
//                               const LatticeFrame&, FlipStack& stack,
//                               OnWrite&& on_write);
//
// The oracle is today's algorithm (lawson.hpp at 09b005b), copied below, not
// the five-argument overload, which 21a may turn into a wrapper of the new one.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include "refinement_fixtures.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

using terrain::TriangleIndices;
using terrain::mesh::FlipStack;
using terrain::mesh::kNoNeighbour;
using terrain::mesh::LatticeFrame;
using terrain::mesh::LatticeMesh;
using terrain::mesh::LatticeVertex;
using terrain::mesh::legalise_around;
using terrain::pred::DefaultKernel;
using refinement_fixtures::orient;
using refinement_fixtures::RC;

namespace {

// legalise_around as of 09b005b, with its own vector.
std::size_t reference_around(LatticeMesh& m, std::uint32_t q, std::span<const std::uint32_t> seeds,
                             const LatticeFrame& f, std::vector<std::uint32_t>& writes) {
    std::vector<std::uint32_t> stack(seeds.begin(), seeds.end());
    std::size_t flips = 0;
    while (!stack.empty()) {
        const auto t = stack.back();
        stack.pop_back();
        const auto& tri = m.triangles()[t];
        unsigned i = 0;
        while (i < 3 && tri[i] != q) ++i;
        if (i == 3 || !terrain::mesh::detail::must_flip<DefaultKernel>(m, t, (i + 1) % 3, f)) continue;
        const auto u = m.neighbours(t)[(i + 1) % 3];
        m.flip(t, (i + 1) % 3);
        ++flips;
        writes.push_back(t);
        writes.push_back(u);
        stack.push_back(t);
        stack.push_back(u);
    }
    return flips;
}

// A (n+1) x (n+1) grid of nodes `step` apart, two triangles per cell, the
// ring constrained. Diagonal direction alternates by cell when `mixed`, so
// the start already has Delaunay ties both ways.
LatticeMesh grid(std::uint32_t n, std::uint32_t step, bool mixed) {
    std::vector<LatticeVertex> v;
    for (std::uint32_t r = 0; r <= n; ++r)
        for (std::uint32_t c = 0; c <= n; ++c) v.push_back(LatticeVertex{r * step, c * step});
    const auto at = [n](std::uint32_t r, std::uint32_t c) { return r * (n + 1) + c; };
    std::vector<TriangleIndices> t;
    std::vector<std::uint8_t> bits;
    std::vector<std::array<std::uint32_t, 3>> masks;
    for (std::uint32_t r = 0; r < n; ++r)
        for (std::uint32_t c = 0; c < n; ++c) {
            const auto tl = at(r, c), tr = at(r, c + 1), bl = at(r + 1, c), br = at(r + 1, c + 1);
            if (!mixed || (r + c) % 2 == 0) {
                t.push_back({tl, bl, br});
                t.push_back({tl, br, tr});
            } else {
                t.push_back({tl, bl, tr});
                t.push_back({bl, br, tr});
            }
        }
    // Constrain every edge that lies on the ring (one triangle only has it).
    for (const auto& tri : t) {
        std::uint8_t b = 0;
        for (unsigned k = 0; k < 3; ++k) {
            const auto p = v[tri[k]], q = v[tri[(k + 1) % 3]];
            const bool on_ring = (p.row == q.row && (p.row == 0 || p.row == n * step)) ||
                                 (p.col == q.col && (p.col == 0 || p.col == n * step));
            if (on_ring) b = static_cast<std::uint8_t>(b | (1u << k));
        }
        bits.push_back(b);
        masks.push_back({(b & 1u) != 0 ? 1u : 0u, (b & 2u) != 0 ? 1u : 0u, (b & 4u) != 0 ? 1u : 0u});
    }
    auto m = LatticeMesh::build(std::move(v), std::move(t), std::move(bits), std::move(masks));
    REQUIRE(m.has_value());
    return *m;
}

RC rc(const LatticeMesh& m, std::uint32_t i) {
    const auto v = m.vertices()[i].as_node().value();
    return RC{v.row, v.col};
}

// The triangle holding node p, and the edge p lies on if any; nullopt when p
// is already a vertex. Integer orientation, independent of the mesh's own.
struct Hit {
    std::uint32_t t;
    std::optional<unsigned> edge;
};
std::optional<Hit> locate(const LatticeMesh& m, LatticeVertex p) {
    const RC x{p.row, p.col};
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t) {
        const auto& tri = m.triangles()[t];
        std::array<std::int64_t, 3> o{};
        for (unsigned k = 0; k < 3; ++k) o[k] = orient(rc(m, tri[k]), rc(m, tri[(k + 1) % 3]), x);
        if (o[0] < 0 || o[1] < 0 || o[2] < 0) continue;
        const int zeros = (o[0] == 0) + (o[1] == 0) + (o[2] == 0);
        if (zeros >= 2) return std::nullopt;
        if (zeros == 0) return Hit{t, std::nullopt};
        return Hit{t, static_cast<unsigned>(o[0] == 0 ? 0 : o[1] == 0 ? 1 : 2)};
    }
    return std::nullopt;
}

void require_same_mesh(const LatticeMesh& a, const LatticeMesh& b) {
    REQUIRE(a.triangle_count() == b.triangle_count());
    REQUIRE(a.vertices().size() == b.vertices().size());
    for (std::uint32_t t = 0; t < a.triangle_count(); ++t) {
        CAPTURE(t);
        REQUIRE(a.triangles()[t] == b.triangles()[t]);
        REQUIRE(a.neighbours(t) == b.neighbours(t));
        for (unsigned k = 0; k < 3; ++k) {
            REQUIRE(a.is_constrained(t, k) == b.is_constrained(t, k));
            REQUIRE(a.mask(t, k) == b.mask(t, k));
        }
    }
}

const std::array<LatticeFrame, 3> kFrames{LatticeFrame{1.0, 1.0}, LatticeFrame{10.0, 5.0},
                                          LatticeFrame{0.1, 0.3}};

}  // namespace

TEST_CASE("legalise_around with a caller-owned stack matches today's, insertion by insertion",
          "[lawson][around][stack]") {
    // Random nodes inserted as refine inserts them: split_inside for an
    // interior node, split_edge for one on an edge (the ring's constrained
    // edges included), with refine's seeds. One FlipStack serves every
    // insertion, as refine will reuse it. After each insertion the new
    // overload must agree with the reference on a copy: flip count, the
    // on_write sequence in order, and the whole mesh.
    const std::uint32_t seed = GENERATE(range(1u, 9u));
    const std::size_t frame = GENERATE(std::size_t{0}, 1, 2);
    const bool mixed = GENERATE(false, true);
    CAPTURE(seed, frame, mixed);
    const LatticeFrame f = kFrames[frame];
    auto m = grid(4, 6, mixed);  // nodes 0..24 both ways
    std::mt19937 rng{seed};
    FlipStack stack;
    std::size_t inserted = 0, flipped = 0, on_edges = 0;
    for (int attempt = 0; attempt < 300; ++attempt) {
        const LatticeVertex p{static_cast<std::uint32_t>(rng() % 25), static_cast<std::uint32_t>(rng() % 25)};
        const auto hit = locate(m, p);
        if (!hit) continue;
        const auto t = hit->t;
        const auto before = static_cast<std::uint32_t>(m.triangle_count());
        std::array<std::uint32_t, 4> seeds{t, before, before + 1, before + 1};
        std::size_t n_seeds = 3;
        std::uint32_t q = 0;
        if (!hit->edge) {
            q = m.split_inside(t, p);
        } else {
            const std::uint32_t u = m.neighbours(t)[*hit->edge];
            q = m.split_edge(t, *hit->edge, p);
            n_seeds = u != kNoNeighbour ? 4 : 2;
            if (n_seeds == 4) seeds[3] = u;
            ++on_edges;
        }
        const std::span<const std::uint32_t> s{seeds.data(), n_seeds};
        LatticeMesh expected = m;
        std::vector<std::uint32_t> expected_writes;
        const std::size_t expected_flips = reference_around(expected, q, s, f, expected_writes);

        std::vector<std::uint32_t> writes;
        const std::size_t flips = legalise_around<DefaultKernel>(
            m, q, s, f, stack, [&writes](std::uint32_t w) { writes.push_back(w); });
        CAPTURE(attempt, p.row, p.col);
        REQUIRE(flips == expected_flips);
        REQUIRE(writes == expected_writes);
        require_same_mesh(m, expected);
        ++inserted;
        flipped += flips;
    }
    // The run exercised what it claims to: many insertions, flips, and edge
    // splits (a run with none of these would compare nothing).
    REQUIRE(inserted >= 100);
    REQUIRE(flipped >= 20);
    REQUIRE(on_edges >= 5);
}

TEST_CASE("legalise_around with a caller-owned stack and no flip to make", "[lawson][around][stack]") {
    // Seeds that do not contain q and a q whose fan is already Delaunay: no
    // flip, no write, mesh unchanged, as today.
    auto m = grid(2, 4, false);
    const auto q = m.split_inside(0, LatticeVertex{3, 1});
    LatticeMesh expected = m;
    std::vector<std::uint32_t> expected_writes;
    const std::array<std::uint32_t, 5> seeds{1, 2, 3, 0, 1};
    const std::size_t expected_flips =
        reference_around(expected, q, std::span<const std::uint32_t>{seeds}, kFrames[0], expected_writes);
    FlipStack stack;
    std::vector<std::uint32_t> writes;
    REQUIRE(legalise_around<DefaultKernel>(m, q, std::span<const std::uint32_t>{seeds}, kFrames[0], stack,
                                           [&writes](std::uint32_t w) { writes.push_back(w); }) ==
            expected_flips);
    REQUIRE(writes == expected_writes);
    require_same_mesh(m, expected);
}
