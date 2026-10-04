// Increment 23b: the oracles of support/frozen_oracle.hpp can fail (FZ0).
//
// Each oracle is run on an output it must reject, built from refine() as it
// is today (no frozen mask: this target compiles and passes before 23b's code
// exists, so the oracles are known to work before the red suites lean on
// them), or from a planted copy of an output. The red suites are
// property/prop_refinement_frozen.cpp and property/prop_refinement_seam.cpp.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>

#include "frozen_oracle.hpp"
#include "strip_oracle.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

using frozen_oracle::kSeam;
using frozen_oracle::Lat;
using frozen_oracle::Side;
using strip_oracle::Mesh;
using terrain::Point2;
using terrain::raster::Raster;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;

namespace {

const auto kG = strip_oracle::exact_geometry(9, 9);

Raster<float> spike(int row, int col, float z = 50.0f) {
    std::vector<float> v(81, 0.0f);
    v[static_cast<std::size_t>(row) * 9 + static_cast<std::size_t>(col)] = z;
    return Raster<float>{kG, std::move(v)};
}

// The start, read as if it were an output: z by the oracle's own heights.
Mesh as_output(const Raster<float>& dem, const strip_oracle::Start& s) {
    const auto [z, valid] = strip_oracle::start_z(dem, s.mesh.vertices());
    return Mesh{{s.mesh.vertices().begin(), s.mesh.vertices().end()},
                z,
                valid,
                {s.mesh.triangles().begin(), s.mesh.triangles().end()},
                s.edges,
                s.masks};
}

std::size_t vertex_at(const Mesh& m, Lat p) {
    for (std::size_t i = 0; i < m.vertices.size(); ++i)
        if (strip_oracle::lat(kG, m.vertices[i]) == p) return i;
    FAIL("no vertex there");
    return 0;
}

}  // namespace

TEST_CASE("FZ0: frozen_findings sees a split, a nudge, a move and a changed mask", "[refinement][frozen][oracle]") {
    const auto start = frozen_oracle::seam_domain(kG, 8, 8, 4, 4, {}, Side::Both);
    const auto dem = spike(4, 4);
    const auto& sv = start.mesh.vertices();

    // Clean on the start itself.
    const Mesh still = as_output(dem, start);
    const auto f0 = frozen_oracle::frozen_findings(kG, sv, start.edges, start.masks, kSeam, still);
    REQUIRE(f0.frozen_edges == 1);
    REQUIRE(frozen_oracle::clean(f0));

    // Today's refine splits the seam at the spike: `on` and `missing`.
    RefineOptions o;
    o.tolerance = 0.5;
    o.threads = 1;
    const auto out = refine(dem, start.mesh, std::span<const std::array<std::uint32_t, 2>>{start.edges},
                            std::span<const std::uint32_t>{start.masks}, o);
    REQUIRE(out.ok());
    const Mesh split = strip_oracle::mesh_of(out);
    const auto f1 = frozen_oracle::frozen_findings(kG, sv, start.edges, start.masks, kSeam, split);
    REQUIRE(f1.on >= 1);
    REQUIRE(f1.missing == 1);

    // A vertex a hair off the seam: `near`, not `on`.
    Mesh nudged = split;
    nudged.vertices[vertex_at(split, Lat{4, 4})].x += 4e-10;  // 2e-10 cells at dx = 2
    const auto f2 = frozen_oracle::frozen_findings(kG, sv, start.edges, start.masks, kSeam, nudged);
    REQUIRE(f2.near >= 1);
    REQUIRE(f2.on == f1.on - 1);

    // An end moved, and a mask changed, on the unsplit start.
    Mesh moved = still;
    moved.vertices[vertex_at(still, Lat{4, 0})].y -= 1.0;
    REQUIRE(frozen_oracle::frozen_findings(kG, sv, start.edges, start.masks, kSeam, moved).moved == 1);
    Mesh remasked = still;
    for (auto& m : remasked.masks)
        if (m == kSeam) m = kSeam | 2u;
    REQUIRE(frozen_oracle::frozen_findings(kG, sv, start.edges, start.masks, kSeam, remasked).wrong_mask == 1);
}

TEST_CASE("FZ0: node_findings_off_frozen skips the seam's nodes and nothing else", "[refinement][frozen][oracle]") {
    const auto start = frozen_oracle::seam_domain(kG, 8, 8, 4, 4, {}, Side::Both);
    const auto frozen = frozen_oracle::frozen_segments(kG, start.mesh.vertices(), start.edges, start.masks, kSeam);
    REQUIRE(frozen.size() == 1);

    // The spike on the seam: the plain oracle reports it, the frozen one not.
    const auto on_seam = spike(4, 4);
    const Mesh m1 = as_output(on_seam, start);
    REQUIRE(strip_oracle::node_findings(on_seam, m1, 0.5).over >= 1);
    REQUIRE(frozen_oracle::node_findings_off_frozen(on_seam, m1, 0.5, frozen).over == 0);

    // A spike one column off the seam is reported by both.
    const auto off_seam = spike(4, 3);
    const Mesh m2 = as_output(off_seam, start);
    REQUIRE(frozen_oracle::node_findings_off_frozen(off_seam, m2, 0.5, frozen).over >= 1);
}

TEST_CASE("FZ0: seam_findings measures check points against the polyline", "[refinement][frozen][oracle]") {
    // The seam col 4 from row 0 to row 8, over a spike at (4, 4).
    const auto dem = spike(4, 4);
    const std::vector<Point2> ends{strip_oracle::world(kG, 4, 0), strip_oracle::world(kG, 4, 8)};
    const auto pts = strip_oracle::ruled_points(dem, ends, {{0, 1}});
    const Lat a{4, 0}, b{4, 8};

    // Through the ends alone: the spike and its midpoints are over.
    const frozen_oracle::Polyline ends_only{{0.0, 1.0}, {0.0, 0.0}, {true, true}};
    const auto f0 = frozen_oracle::seam_findings(a, b, ends_only, pts, 0.5);
    REQUIRE(f0.checked == pts.size());
    REQUIRE(f0.over >= 1);

    // Through every node of the seam: nothing is over, even at tolerance 0.
    frozen_oracle::Polyline nodes;
    for (int r = 0; r <= 8; ++r) {
        nodes.t.push_back(r / 8.0);
        nodes.z.push_back(r == 4 ? 50.0 : 0.0);
        nodes.valid.push_back(true);
    }
    REQUIRE(frozen_oracle::seam_findings(a, b, nodes, pts, 0.0).over == 0);

    // A piece with an invalid end is void, not over.
    frozen_oracle::Polyline void_end = ends_only;
    void_end.valid[0] = false;
    const auto f2 = frozen_oracle::seam_findings(a, b, void_end, pts, 0.5);
    REQUIRE(f2.on_void == pts.size());
    REQUIRE(f2.over == 0);
}

TEST_CASE("FZ0: seam_domain builds valid counter-clockwise pieces whose seams agree", "[refinement][frozen][oracle]") {
    // Left, right and whole, with two inner seam points: every triangle CCW in
    // the exact kernel, and the pieces' seam vertices at the same world points.
    const std::vector<Lat> inner{{4, 2.5}, {4, 6}};
    for (const auto side : {Side::Both, Side::Left, Side::Right}) {
        const auto s = frozen_oracle::seam_domain(kG, 8, 8, 4, 4, inner, side);
        const Mesh m = as_output(spike(4, 4), s);
        REQUIRE(strip_oracle::node_findings(spike(4, 4), m, 1e9).not_ccw == 0);
        const auto seq = frozen_oracle::seam_sequence(kG, m, Lat{4, 0}, Lat{4, 8});
        REQUIRE(seq.size() == 4);
        REQUIRE(seq[1].first == strip_oracle::world(kG, 4, 2.5));
        REQUIRE(seq[2].first == strip_oracle::world(kG, 4, 6));
    }
}
