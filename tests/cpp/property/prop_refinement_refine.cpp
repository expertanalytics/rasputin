// Increment 14 (docs/increments/14-adaptive-refinement.md, R4 to R8): the
// refinement loop. Tests T1 to T8. Property and integration tests, not
// mutation-tested (the design names test_mesh_lattice_split and
// test_refinement_scan as the invariant-critical pair).
//
// Interface: include/terrain/refinement/refine.hpp.
// A clockwise and a zero-area start triangle both give NotCounterClockwise.
//
// THE T3 ORACLE SHARES NO CODE WITH scan.hpp. It lives in this file and in
// refinement_fixtures.hpp: every DEM node against every output triangle, by
// integer orientation, with the plane from barycentric weights.
//
// "Valid DEM nodes not covered" (R6) is checked against this oracle: valid
// nodes that are neither a mesh vertex nor in the closed node set of any
// triangle with three valid vertices. R6 defines it as a per-void-triangle
// sum, which counts a node on an edge between two void triangles twice; see
// the handback.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/refine.hpp>

#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <thread>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::TriangleIndices;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;
using terrain::refinement::RefineStatus;
using namespace refinement_fixtures;

namespace {

const float kNaN = std::numeric_limits<float>::quiet_NaN();

using Edges = std::vector<std::array<std::uint32_t, 2>>;

template <typename Outcome>
std::vector<RC> lattice_of(const RasterGeometry& g, const Outcome& out) {
    std::vector<RC> v;
    for (const Point2& p : out.vertices) {
        const RC q = to_rc(g, p);
        const Point2 back = g.node({static_cast<std::size_t>(q.row), static_cast<std::size_t>(q.col)});
        REQUIRE(back.x == p.x);  // R4: every output vertex is exactly a node
        REQUIRE(back.y == p.y);
        v.push_back(q);
    }
    return v;
}

auto run(const Raster<float>& dem, const StartMesh& s, double tol, unsigned threads = 1) {
    return refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks},
                  RefineOptions{.tolerance = tol, .threads = threads});
}

std::int64_t area2(const std::vector<RC>& v, const std::vector<TriangleIndices>& t) {
    std::int64_t s = 0;
    for (const auto& tri : t) s += orient(v[tri[0]], v[tri[1]], v[tri[2]]);
    return s;
}

// A start mesh from lattice nodes, for the hand-built T2 fixtures. The whole
// perimeter is listed in `edges` with the given masks.
StartMesh hand_mesh(const RasterGeometry& g, std::vector<RC> nodes,
                    std::vector<TriangleIndices> tris, Edges edges, std::vector<std::uint32_t> masks) {
    std::vector<Point2> xy;
    for (const RC& n : nodes)
        xy.push_back(g.node({static_cast<std::size_t>(n.row), static_cast<std::size_t>(n.col)}));
    std::vector<std::uint8_t> constrained(tris.size(), 0);
    StartMesh s;
    s.mesh = IndexedMesh2{std::move(xy), std::move(tris), std::move(constrained)};
    s.edges = std::move(edges);
    s.masks = std::move(masks);
    s.lattice = std::move(nodes);
    return s;
}

// The 3 x 3-node square, diagonal (0,0)-(2,2); perimeter masks top 1, right 2,
// bottom 4, left 8.
StartMesh square(const RasterGeometry& g) {
    return hand_mesh(g, {{0, 0}, {2, 0}, {2, 2}, {0, 2}}, {{0, 1, 2}, {0, 2, 3}},
                     {{0, 3}, {2, 3}, {1, 2}, {0, 1}}, {1, 2, 4, 8});
}

double dem_at(const Raster<float>& dem, RC p) {
    return static_cast<double>(dem.value_at({static_cast<std::size_t>(p.row), static_cast<std::size_t>(p.col)}));
}

bool valid_at(const Raster<float>& dem, RC p) {
    return !dem.is_nodata({static_cast<std::size_t>(p.row), static_cast<std::size_t>(p.col)});
}

// T3, T4, T5 and R6's count, on one outcome. Independent of scan.hpp.
template <typename Outcome>
void check_properties(const Raster<float>& dem, const StartMesh& start, const Outcome& out, double tol) {
    const RasterGeometry& g = dem.geometry();
    const auto rows = static_cast<std::int64_t>(g.rows()), cols = static_cast<std::int64_t>(g.cols());
    const auto v = lattice_of(g, out);
    REQUIRE(out.z.size() == v.size());
    REQUIRE(out.valid.size() == v.size());
    REQUIRE(out.inserted == v.size() - start.lattice.size());
    double zmax = 0.0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        CAPTURE(i);
        REQUIRE(static_cast<bool>(out.valid[i]) == valid_at(dem, v[i]));
        if (valid_at(dem, v[i])) {
            REQUIRE(out.z[i] == dem_at(dem, v[i]));
            zmax = std::max(zmax, std::abs(out.z[i]));
        }
    }

    // T4: orientation, area, manifold edges, boundary edges on the tile edge,
    // no hanging vertex.
    const auto& T = out.triangles;
    for (const auto& tri : T) REQUIRE(orient(v[tri[0]], v[tri[1]], v[tri[2]]) > 0);
    REQUIRE(area2(v, T) == area2(start.lattice, std::vector<TriangleIndices>(
                                                    start.mesh.triangles().begin(),
                                                    start.mesh.triangles().end())));
    auto on_tile_edge = [&](RC a, RC b) {
        return (a.row == 0 && b.row == 0) || (a.row == rows - 1 && b.row == rows - 1)
            || (a.col == 0 && b.col == 0) || (a.col == cols - 1 && b.col == cols - 1);
    };
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
    for (const auto& tri : T)
        for (unsigned k = 0; k < 3; ++k)
            ++uses[{std::min(tri[k], tri[(k + 1) % 3]), std::max(tri[k], tri[(k + 1) % 3])}];
    for (const auto& [e, n] : uses) {
        CAPTURE(e.first, e.second);
        REQUIRE(n <= 2);
        REQUIRE((n == 1) == on_tile_edge(v[e.first], v[e.second]));
        for (std::size_t i = 0; i < v.size(); ++i) REQUIRE_FALSE(on_open_segment(v[e.first], v[e.second], v[i]));
    }

    // T5: constraint edges are exactly the boundary mesh edges, each once,
    // carrying the side's mask; per side, the pieces' lengths sum to the side.
    REQUIRE(out.edges.size() == out.masks.size());
    std::map<std::uint32_t, std::int64_t> side_length;
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> seen;
    for (std::size_t i = 0; i < out.edges.size(); ++i) {
        const auto [a, b] = out.edges[i];
        CAPTURE(i, a, b);
        REQUIRE(on_tile_edge(v[a], v[b]));
        REQUIRE(uses.count({std::min(a, b), std::max(a, b)}) == 1);
        REQUIRE(++seen[{std::min(a, b), std::max(a, b)}] == 1);
        const std::uint32_t side = v[a].row == 0 && v[b].row == 0               ? 1u
                                   : v[a].col == cols - 1 && v[b].col == cols - 1 ? 2u
                                   : v[a].row == rows - 1 && v[b].row == rows - 1 ? 4u
                                                                                  : 8u;
        REQUIRE(out.masks[i] == side);
        side_length[side] += std::abs(v[a].row - v[b].row) + std::abs(v[a].col - v[b].col);
    }
    REQUIRE(seen.size() == static_cast<std::size_t>(std::count_if(
                               uses.begin(), uses.end(), [](const auto& u) { return u.second == 1; })));
    REQUIRE(side_length[1] == cols - 1);
    REQUIRE(side_length[4] == cols - 1);
    REQUIRE(side_length[2] == rows - 1);
    REQUIRE(side_length[8] == rows - 1);

    // T3: brute force, every node against every triangle with three valid
    // vertices. Also collects R6's oracle.
    std::vector<char> covered(static_cast<std::size_t>(rows * cols), 0);
    for (const RC& p : v) covered[static_cast<std::size_t>(p.row * cols + p.col)] = 1;
    double worst = 0.0;
    for (const auto& tri : T) {
        const RC a = v[tri[0]], b = v[tri[1]], c = v[tri[2]];
        if (!(valid_at(dem, a) && valid_at(dem, b) && valid_at(dem, c))) continue;
        const double two_a = static_cast<double>(orient(a, b, c));
        for (std::int64_t r = 0; r < rows; ++r)
            for (std::int64_t col = 0; col < cols; ++col) {
                const RC p{r, col};
                if (!in_closed(a, b, c, p) || p == a || p == b || p == c || !valid_at(dem, p)) continue;
                covered[static_cast<std::size_t>(r * cols + col)] = 1;
                const double plane = (static_cast<double>(orient(b, c, p)) * dem_at(dem, a)
                                      + static_cast<double>(orient(c, a, p)) * dem_at(dem, b)
                                      + static_cast<double>(orient(a, b, p)) * dem_at(dem, c)) / two_a;
                const double err = std::abs(dem_at(dem, p) - plane);
                worst = std::max(worst, err);
                CAPTURE(p.row, p.col, err, tol);
                REQUIRE(err <= tol + 1e-9 * zmax);
            }
    }
    REQUIRE(out.max_error <= tol);
    REQUIRE(std::abs(out.max_error - worst) <= 1e-9 * std::max(1.0, zmax));

    std::size_t uncovered = 0;
    for (std::int64_t r = 0; r < rows; ++r)
        for (std::int64_t col = 0; col < cols; ++col)
            if (valid_at(dem, {r, col}) && !covered[static_cast<std::size_t>(r * cols + col)]) ++uncovered;
    REQUIRE(out.uncovered == uncovered);
}

Raster<float> plane_dem(std::size_t rows, std::size_t cols) {
    std::vector<float> z(rows * cols);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c)
            z[r * cols + c] = static_cast<float>(3.0 * static_cast<double>(c) - 2.0 * static_cast<double>(r) + 7.0);
    return Raster<float>{geometry(rows, cols), std::move(z)};
}

Raster<float> flat_with_peak(std::size_t n, RC peak) {
    std::vector<float> z(n * n, 0.0f);
    z[static_cast<std::size_t>(peak.row) * n + static_cast<std::size_t>(peak.col)] = 10.0f;
    return Raster<float>{geometry(n, n), std::move(z)};
}

}  // namespace

TEST_CASE("T1: a plane DEM needs no refinement", "[refinement][refine]") {
    const auto dem = plane_dem(9, 13);
    const auto start = grid_mesh(dem.geometry(), 4);
    const auto out = run(dem, start, 0.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 0);
    REQUIRE(out.rounds == 1);
    REQUIRE(out.max_error == 0.0);
    REQUIRE(std::vector<TriangleIndices>(start.mesh.triangles().begin(), start.mesh.triangles().end())
            == out.triangles);
    REQUIRE(out.vertices.size() == start.mesh.vertices().size());
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        REQUIRE(out.vertices[i].x == start.mesh.vertices()[i].x);
        REQUIRE(out.vertices[i].y == start.mesh.vertices()[i].y);
    }
    check_properties(dem, start, out, 0.0);
}

TEST_CASE("T2: a single peak refines around it and nowhere else", "[refinement][refine]") {
    SECTION("strictly inside a start triangle: one insertion, a 1 -> 3 fan") {
        // (2,1) is the only interior node of (0,0) (3,0) (3,3); the children
        // hold no node off the parent's edges, whose planes stay 0.
        const auto dem = flat_with_peak(4, {2, 1});
        const auto start = hand_mesh(dem.geometry(), {{0, 0}, {3, 0}, {3, 3}, {0, 3}},
                                     {{0, 1, 2}, {0, 2, 3}}, {{0, 3}, {2, 3}, {1, 2}, {0, 1}},
                                     {1, 2, 4, 8});
        const auto out = run(dem, start, 1.0);
        REQUIRE(out.ok());
        REQUIRE(out.inserted == 1);
        REQUIRE(out.triangles.size() == 4);
        REQUIRE(to_rc(dem.geometry(), out.vertices.back()) == RC{2, 1});
        REQUIRE(out.z.back() == 10.0);
        REQUIRE(out.max_error == 0.0);
        check_properties(dem, start, out, 1.0);
    }
    SECTION("on the interior start edge: one insertion, four triangles replace two") {
        const auto dem = flat_with_peak(3, {1, 1});
        const auto start = square(dem.geometry());
        const auto out = run(dem, start, 1.0);
        REQUIRE(out.ok());
        REQUIRE(out.inserted == 1);
        REQUIRE(out.triangles.size() == 4);
        REQUIRE(to_rc(dem.geometry(), out.vertices.back()) == RC{1, 1});
        check_properties(dem, start, out, 1.0);
    }
    SECTION("on the tile boundary: one insertion, both halves keep the parent's mask") {
        const auto dem = flat_with_peak(3, {0, 1});
        const auto start = square(dem.geometry());
        const auto out = run(dem, start, 1.0);
        REQUIRE(out.ok());
        REQUIRE(out.inserted == 1);
        REQUIRE(out.triangles.size() == 3);
        REQUIRE(out.edges.size() == 5);
        check_properties(dem, start, out, 1.0);  // includes: the top halves carry mask 1
    }
}

TEST_CASE("T3-T5: the tolerance oracle, conformity and constraints", "[refinement][refine]") {
    const bool rough = GENERATE(false, true);
    const double tol = GENERATE(0.0, 0.5, 5.0);
    const std::uint32_t seed = GENERATE(1u, 2u);
    CAPTURE(rough, tol, seed);
    const std::size_t rows = 17, cols = 15;
    const Raster<float> dem{geometry(rows, cols),
                            rough ? rough_dem(rows, cols, seed) : smooth_dem(rows, cols, seed)};
    const auto start = grid_mesh(dem.geometry(), 8);
    const auto out = run(dem, start, tol);
    REQUIRE(out.ok());
    check_properties(dem, start, out, tol);
    if (tol == 0.0 && rough) REQUIRE(out.vertices.size() >= rows * cols / 2);
}

TEST_CASE("T6: the output is bit-identical for 1, 2, 7 and all threads", "[refinement][refine]") {
    const std::size_t n = 33;
    const Raster<float> dem{geometry(n, n), rough_dem(n, n, 7)};
    const auto start = grid_mesh(dem.geometry(), 16);
    const auto ref = run(dem, start, 3.0, 1);
    REQUIRE(ref.ok());
    REQUIRE(ref.inserted > 0);
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    const auto out = run(dem, start, 3.0, threads);
    REQUIRE(out.ok());
    REQUIRE(out.vertices.size() == ref.vertices.size());
    for (std::size_t i = 0; i < ref.vertices.size(); ++i) {
        REQUIRE(out.vertices[i].x == ref.vertices[i].x);
        REQUIRE(out.vertices[i].y == ref.vertices[i].y);
    }
    REQUIRE(out.z == ref.z);
    REQUIRE(out.valid == ref.valid);
    REQUIRE(out.triangles == ref.triangles);
    REQUIRE(out.edges == ref.edges);
    REQUIRE(out.masks == ref.masks);
    REQUIRE(out.rounds == ref.rounds);
    REQUIRE(out.inserted == ref.inserted);
    REQUIRE(out.max_error == ref.max_error);
    REQUIRE(out.uncovered == ref.uncovered);
}

TEST_CASE("T7: NoData is carved, not refined, and what is lost is counted", "[refinement][refine]") {
    const std::size_t rows = 17, cols = 17;
    const bool sentinel = GENERATE(false, true);
    CAPTURE(sentinel);
    auto z = rough_dem(rows, cols, 3);
    const float hole = sentinel ? -32767.0f : kNaN;
    for (std::size_t r = 0; r <= 4; ++r)
        for (std::size_t c = 0; c <= 6; ++c) z[r * cols + c] = hole;
    const Raster<float> dem{geometry(rows, cols), std::move(z),
                            sentinel ? std::optional<float>{-32767.0f} : std::nullopt};
    const auto start = grid_mesh(dem.geometry(), 8);
    const auto out = run(dem, start, 2.0);
    REQUIRE(out.ok());
    check_properties(dem, start, out, 2.0);  // includes the uncovered-count oracle

    // Inserted vertices are always valid (R6).
    for (std::size_t i = start.lattice.size(); i < out.vertices.size(); ++i)
        REQUIRE(static_cast<bool>(out.valid[i]));
    // The carving went as far as R6 says: no void triangle holds a valid node
    // strictly inside it.
    const auto v = lattice_of(dem.geometry(), out);
    for (const auto& tri : out.triangles) {
        const RC a = v[tri[0]], b = v[tri[1]], c = v[tri[2]];
        if (valid_at(dem, a) && valid_at(dem, b) && valid_at(dem, c)) continue;
        for (std::int64_t r = 0; r < static_cast<std::int64_t>(rows); ++r)
            for (std::int64_t col = 0; col < static_cast<std::int64_t>(cols); ++col) {
                const RC p{r, col};
                const bool strictly = orient(a, b, p) > 0 && orient(b, c, p) > 0 && orient(c, a, p) > 0;
                if (strictly) REQUIRE_FALSE(valid_at(dem, p));
            }
    }
}

TEST_CASE("T7: an all-NoData DEM leaves no triangle with data", "[refinement][refine]") {
    const std::size_t n = 9;
    const Raster<float> dem{geometry(n, n), std::vector<float>(n * n, kNaN)};
    const auto start = grid_mesh(dem.geometry(), 4);
    const auto out = run(dem, start, 1.0);
    REQUIRE(out.ok());
    REQUIRE(out.inserted == 0);
    for (std::size_t i = 0; i < out.valid.size(); ++i) REQUIRE_FALSE(static_cast<bool>(out.valid[i]));
    REQUIRE(out.uncovered == 0);
}

TEST_CASE("T8: refusals come back as a status, not a crash", "[refinement][refine]") {
    const auto dem = plane_dem(9, 9);
    const auto good = grid_mesh(dem.geometry(), 4);

    SECTION("a tolerance that is negative, NaN or infinite") {
        const double tol = GENERATE(-1.0, std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::infinity());
        CAPTURE(tol);
        const auto out = run(dem, good, tol);
        REQUIRE(out.status == RefineStatus::InvalidTolerance);
        REQUIRE_FALSE(out.ok());
        REQUIRE_FALSE(out.message.empty());
    }
    SECTION("an off-lattice start vertex") {
        std::vector<Point2> xy(good.mesh.vertices().begin(), good.mesh.vertices().end());
        xy[4].x += 0.5 * dem.geometry().delta_x();  // the centre node, off by half a cell
        StartMesh bad = good;
        bad.mesh = IndexedMesh2{std::move(xy),
                                {good.mesh.triangles().begin(), good.mesh.triangles().end()},
                                {good.mesh.constrained_edges().begin(), good.mesh.constrained_edges().end()}};
        const auto out = run(dem, bad, 1.0);
        REQUIRE(out.status == RefineStatus::OffLattice);
        REQUIRE_FALSE(out.message.empty());
    }
    SECTION("a clockwise or a zero-area start triangle") {
        const bool degenerate = GENERATE(false, true);
        CAPTURE(degenerate);
        std::vector<TriangleIndices> tris(good.mesh.triangles().begin(), good.mesh.triangles().end());
        if (degenerate) tris[0] = {0, 1, 2};  // three nodes along the top row
        else std::swap(tris[0][1], tris[0][2]);
        StartMesh bad = good;
        bad.mesh = IndexedMesh2{{good.mesh.vertices().begin(), good.mesh.vertices().end()}, std::move(tris),
                                {good.mesh.constrained_edges().begin(), good.mesh.constrained_edges().end()}};
        const auto out = run(dem, bad, 1.0);
        REQUIRE(out.status == RefineStatus::NotCounterClockwise);
        REQUIRE_FALSE(out.message.empty());
    }
}
