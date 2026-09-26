// Increment 14 (docs/increments/14-adaptive-refinement.md, R4 to R8): the
// refinement loop. Tests T1 to T8. Property and integration tests, not
// mutation-tested (the design names test_mesh_lattice_split and
// test_refinement_scan as the invariant-critical pair).
//
// Increment 14b (docs/increments/14b-delaunay-insertion.md) re-runs them under
// Delaunay insertion. THE T3 ORACLE IS 14b's SECOND INVARIANT-CRITICAL SUITE:
// its mutant is a flipped slot not marked touched, whose stale scan result
// must surface as a tolerance or conformity failure here. check_properties
// now also asserts R10 (constrained Delaunay, in world coordinates, with
// DefaultKernel::incircle on the outcome's vertices; the geometry's
// translation is integral so no sign can move), and the T3 cases require
// flips > 0 on rough terrain so the oracle is known to see flipped slots.
// T2 is amended: its fixed counts assumed fans. T12 is new (R10 on a cone and
// on an island with a step coast). The outcome gains `flips` (design, Files).
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
#include <terrain/predicates/default_kernel.hpp>
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
#include <set>
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

// R10: every interior edge that is not a constraint edge has neither apex
// strictly inside the other triangle's circle. World coordinates, straight
// from the outcome.
template <typename Outcome>
void check_delaunay(const Outcome& out) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> constrained;
    for (const auto& e : out.edges) constrained.insert(std::minmax(e[0], e[1]));
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> sides;
    for (std::size_t t = 0; t < out.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            sides[std::minmax(out.triangles[t][k], out.triangles[t][(k + 1) % 3])].push_back(t);
    std::size_t bad = 0;
    for (const auto& [e, ts] : sides) {
        if (ts.size() != 2 || constrained.count(e) != 0) continue;
        for (unsigned s = 0; s < 2; ++s) {
            const auto& tri = out.triangles[ts[s]];
            const auto& other = out.triangles[ts[1 - s]];
            std::uint32_t apex = other[0];
            for (const auto x : other)
                if (x != e.first && x != e.second) apex = x;
            if (terrain::pred::DefaultKernel::incircle(out.vertices[tri[0]], out.vertices[tri[1]],
                                                       out.vertices[tri[2]], out.vertices[apex])
                == terrain::pred::Incircle::Inside) {
                UNSCOPED_INFO("edge " << e.first << "-" << e.second << " apex " << apex);
                ++bad;
            }
        }
    }
    REQUIRE(bad == 0);
}

// T3's brute force, restricted to each triangle's bounding box (the nodes
// outside it cannot be in the closed triangle). Returns the worst error and
// marks covered nodes.
template <typename Outcome>
double check_tolerance(const Raster<float>& dem, const Outcome& out, const std::vector<RC>& v,
                       double tol, double zmax, std::vector<char>& covered) {
    const auto cols = static_cast<std::int64_t>(dem.geometry().cols());
    double worst = 0.0;
    for (const auto& tri : out.triangles) {
        const RC a = v[tri[0]], b = v[tri[1]], c = v[tri[2]];
        if (!(valid_at(dem, a) && valid_at(dem, b) && valid_at(dem, c))) continue;
        const double two_a = static_cast<double>(orient(a, b, c));
        for (std::int64_t r = std::min({a.row, b.row, c.row}); r <= std::max({a.row, b.row, c.row}); ++r)
            for (std::int64_t col = std::min({a.col, b.col, c.col}); col <= std::max({a.col, b.col, c.col}); ++col) {
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
    return worst;
}

// T3, T4, T5, R6's count and R10, on one outcome. Independent of scan.hpp.
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
    const double worst = check_tolerance(dem, out, v, tol, zmax, covered);
    REQUIRE(out.max_error <= tol);
    REQUIRE(std::abs(out.max_error - worst) <= 1e-9 * std::max(1.0, zmax));

    std::size_t uncovered = 0;
    for (std::int64_t r = 0; r < rows; ++r)
        for (std::int64_t col = 0; col < cols; ++col)
            if (valid_at(dem, {r, col}) && !covered[static_cast<std::size_t>(r * cols + col)]) ++uncovered;
    REQUIRE(out.uncovered == uncovered);

    check_delaunay(out);  // R10
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
    REQUIRE(out.flips == 0);  // a grid of rectangles is cocircular: ties never flip (R3)
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

TEST_CASE("T2: a single peak is inserted first and the result is Delaunay", "[refinement][refine]") {
    // Amended for 14b. Under fans the counts below were children per split;
    // under Delaunay insertion the first split may flip, which exposes nodes
    // the fan kept on edges, so the count is not fixed. What stays fixed: the
    // peak is the first vertex inserted, and check_properties (tolerance,
    // conformity, constraints, R10) holds. Geometry: dx = 10, dy = 5.
    SECTION("strictly inside a start triangle: the fan must flip") {
        // After the 1 -> 3 fan at (2,1), (0,3) is strictly inside the circle
        // of (0,0) (3,3) (2,1) in world (centre offset (27.5, 17.5) from
        // (0,0), radius^2 1062.5; (0,3) at 312.5), so legalisation flips.
        const auto dem = flat_with_peak(4, {2, 1});
        const auto start = hand_mesh(dem.geometry(), {{0, 0}, {3, 0}, {3, 3}, {0, 3}},
                                     {{0, 1, 2}, {0, 2, 3}}, {{0, 3}, {2, 3}, {1, 2}, {0, 1}},
                                     {1, 2, 4, 8});
        const auto out = run(dem, start, 1.0);
        REQUIRE(out.ok());
        REQUIRE(out.inserted >= 1);
        REQUIRE(out.flips >= 1);
        REQUIRE(to_rc(dem.geometry(), out.vertices[4]) == RC{2, 1});
        REQUIRE(out.z[4] == 10.0);
        check_properties(dem, start, out, 1.0);
    }
    SECTION("on the interior start edge: one insertion, four triangles, no flip") {
        // Every edge opposite the new vertex is on the tile boundary.
        const auto dem = flat_with_peak(3, {1, 1});
        const auto start = square(dem.geometry());
        const auto out = run(dem, start, 1.0);
        REQUIRE(out.ok());
        REQUIRE(out.inserted == 1);
        REQUIRE(out.flips == 0);
        REQUIRE(out.triangles.size() == 4);
        REQUIRE(to_rc(dem.geometry(), out.vertices.back()) == RC{1, 1});
        check_properties(dem, start, out, 1.0);
    }
    SECTION("on the tile boundary: inserted first, both halves keep the parent's mask") {
        // After the boundary split at (0,1), (2,0) is strictly inside the
        // circle of (0,1) (0,0) (2,2), so the diagonal flips and (1,1) then
        // needs a vertex of its own.
        const auto dem = flat_with_peak(3, {0, 1});
        const auto start = square(dem.geometry());
        const auto out = run(dem, start, 1.0);
        REQUIRE(out.ok());
        REQUIRE(out.inserted >= 1);
        REQUIRE(out.flips >= 1);
        REQUIRE(to_rc(dem.geometry(), out.vertices[4]) == RC{0, 1});
        REQUIRE(out.edges.size() >= 5);
        check_properties(dem, start, out, 1.0);  // includes: the top halves carry mask 1
    }
}

TEST_CASE("T3-T5: the tolerance oracle, conformity and constraints", "[refinement][refine]") {
    const bool rough = GENERATE(false, true);
    const double tol = GENERATE(0.0, 0.5, 5.0);
    const std::uint32_t seed = GENERATE(1u, 2u, 3u, 4u);
    CAPTURE(rough, tol, seed);
    const std::size_t rows = 17, cols = 15;
    const Raster<float> dem{geometry(rows, cols),
                            rough ? rough_dem(rows, cols, seed) : smooth_dem(rows, cols, seed)};
    const auto start = grid_mesh(dem.geometry(), 8);
    const auto out = run(dem, start, tol);
    REQUIRE(out.ok());
    check_properties(dem, start, out, tol);
    if (tol == 0.0 && rough) REQUIRE(out.vertices.size() >= rows * cols / 2);
    // 14b: the oracle must see flipped slots, or a flipped slot left
    // untouched (stale scan result) could not show here.
    if (rough && tol < 5.0) REQUIRE(out.flips > 0);
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
    REQUIRE(out.flips == ref.flips);
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

// ---------------------------------------------------------------- T12 (14b)

namespace {

// Integer x_min, y_max, dx and dy, so the world coordinates are exact integers
// and check_delaunay's incircle answers the same question as refine's frame.
RasterGeometry integral_geometry(std::size_t n, double dx, double dy) {
    return RasterGeometry{1000.0, 2000.0, dx, dy, n, n};
}

// A cone, 100 m at the centre node, falling 1 m per cell in lattice distance.
Raster<float> cone(std::size_t n, double dx, double dy) {
    const double mid = static_cast<double>(n - 1) / 2.0;
    std::vector<float> z(n * n);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            z[r * n + c] = static_cast<float>(
                100.0 - std::hypot(static_cast<double>(r) - mid, static_cast<double>(c) - mid));
    return Raster<float>{integral_geometry(n, dx, dy), std::move(z)};
}

// Flat sea at exactly 0, and land that starts at 3 m on a circular coast and
// rises 0.5 m per cell inland: a step the tolerance cannot smooth over.
Raster<float> island(std::size_t n, double dx, double dy) {
    const double mid = static_cast<double>(n - 1) / 2.0, radius = 0.3 * static_cast<double>(n);
    std::vector<float> z(n * n, 0.0f);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c) {
            const double d = std::hypot(static_cast<double>(r) - mid, static_cast<double>(c) - mid);
            if (d <= radius) z[r * n + c] = static_cast<float>(3.0 + 0.5 * (radius - d));
        }
    return Raster<float>{integral_geometry(n, dx, dy), std::move(z)};
}

// R10 plus what is cheap at 129 x 129: orientation, area, the bounding-box
// tolerance oracle. (check_properties' hanging-vertex scan is O(edges x
// vertices) and is left to the small fixtures.)
template <typename Outcome>
void check_refined(const Raster<float>& dem, const StartMesh& start, const Outcome& out, double tol) {
    REQUIRE(out.ok());
    const auto v = lattice_of(dem.geometry(), out);
    for (const auto& tri : out.triangles) REQUIRE(orient(v[tri[0]], v[tri[1]], v[tri[2]]) > 0);
    REQUIRE(area2(v, out.triangles) == area2(start.lattice, std::vector<TriangleIndices>(
                                                               start.mesh.triangles().begin(),
                                                               start.mesh.triangles().end())));
    double zmax = 0.0;
    for (const double z : out.z) zmax = std::max(zmax, std::abs(z));
    std::vector<char> covered(dem.geometry().rows() * dem.geometry().cols(), 0);
    const double worst = check_tolerance(dem, out, v, tol, zmax, covered);
    REQUIRE(out.max_error <= tol);
    REQUIRE(std::abs(out.max_error - worst) <= 1e-9 * std::max(1.0, zmax));
    check_delaunay(out);
}

}  // namespace

TEST_CASE("T12: the refined mesh is constrained Delaunay on a cone and an island", "[refinement][refine][delaunay]") {
    const std::size_t n = 129;
    const bool on_island = GENERATE(false, true);
    const std::size_t stride = GENERATE(4u, 16u);
    const double tol = GENERATE(5.0, 1.0, 0.0);
    CAPTURE(on_island, stride, tol);
    const auto dem = on_island ? island(n, 2.0, 2.0) : cone(n, 2.0, 2.0);
    const auto start = grid_mesh(dem.geometry(), stride);
    // The start mesh's own error, measured by a refine that may insert nothing.
    const auto untouched = run(dem, start, 1e9);
    REQUIRE(untouched.inserted == 0);
    const double start_error = untouched.max_error;
    CAPTURE(start_error);

    const auto out = run(dem, start, tol);
    check_refined(dem, start, out, tol);
    // Refine inserts exactly when the start mesh is out of tolerance. At 5 m
    // every start mesh here is already within it (1.2 to 4.7 m), so that row
    // checks the no-op path.
    REQUIRE((out.inserted > 0) == (start_error > tol));
    if (out.inserted == 0) REQUIRE(out.flips == 0);
    // Tolerance 0 forces thousands of interior inserts on every combination,
    // so it must flip, and check_refined then catches a flipped slot that is
    // not rescanned. At 1 m a flip is not guaranteed: the cone at stride 4
    // inserts two nodes on start-grid edges and legitimately needs none.
    if (tol == 0.0) REQUIRE(out.flips > 0);
}

TEST_CASE("T12: constrained Delaunay in world coordinates when dx differs from dy", "[refinement][refine][delaunay]") {
    // dy = 3 dx: the lattice-frame Delaunay is not the world one (R3), so a
    // refine that legalised in (col, -row) fails check_delaunay here.
    const std::size_t n = 129;
    const bool on_island = GENERATE(false, true);
    CAPTURE(on_island);
    const auto dem = on_island ? island(n, 1.0, 3.0) : cone(n, 1.0, 3.0);
    const auto start = grid_mesh(dem.geometry(), 16);
    const auto out = run(dem, start, 1.0);
    check_refined(dem, start, out, 1.0);
    REQUIRE(out.flips > 0);
}

TEST_CASE("T12: an interior constraint edge through the cone's apex stays constrained in pieces", "[refinement][refine][delaunay]") {
    // grid_mesh at stride 16 on 129 nodes has 9 columns of vertices; column 4
    // is col 64, through the apex. Its eight vertical cell sides are added as
    // constraint edges with mask 16. The cone wants edges across that line;
    // R5 forbids flipping it, so at the end the line is still covered, piece
    // by piece, by constraint edges carrying mask 16.
    const std::size_t n = 129;
    const auto dem = cone(n, 2.0, 2.0);
    auto start = grid_mesh(dem.geometry(), 16);
    const std::uint32_t nc = 9, mid_col = 4;
    for (std::uint32_t i = 0; i + 1 < nc; ++i) {
        start.edges.push_back({i * nc + mid_col, (i + 1) * nc + mid_col});
        start.masks.push_back(16);
    }
    const double tol = GENERATE(1.0, 0.0);
    CAPTURE(tol);
    const auto out = run(dem, start, tol);
    check_refined(dem, start, out, tol);  // check_delaunay skips constraint edges

    const auto v = lattice_of(dem.geometry(), out);
    std::set<std::pair<std::uint32_t, std::uint32_t>> mesh_edges;
    for (const auto& tri : out.triangles)
        for (unsigned k = 0; k < 3; ++k) mesh_edges.insert(std::minmax(tri[k], tri[(k + 1) % 3]));
    std::int64_t covered_rows = 0;
    std::set<std::int64_t> starts;
    for (std::size_t i = 0; i < out.edges.size(); ++i) {
        if (out.masks[i] != 16) continue;
        const auto [a, b] = out.edges[i];
        CAPTURE(i, a, b);
        REQUIRE(v[a].col == 64);
        REQUIRE(v[b].col == 64);
        REQUIRE(mesh_edges.count(std::minmax(a, b)) == 1);
        REQUIRE(starts.insert(std::min(v[a].row, v[b].row)).second);  // pieces do not overlap
        covered_rows += std::abs(v[a].row - v[b].row);
    }
    REQUIRE(covered_rows == static_cast<std::int64_t>(n - 1));
    REQUIRE(out.flips > 0);
}
