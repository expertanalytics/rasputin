// Increment 20 (docs/increments/20-start-quality.md, R1, R9 to R11): refine with
// the start-quality pass on. Property and integration tests: 14b's T3 and T6
// re-run with quality on, Q7 (the off switch is bit-identical to increment
// 18), Q8 through refine, Q10 (NoData). Not the invariant-critical suite;
// that is unit/test_mesh_quality.cpp.
//
// Interface, as the design's R11 names it, with two choices made here:
//
//   RefineOptions::min_angle_deg      double, default 0 (off); set by member
//                                     assignment, so its declaration order is free
//   RefineOutcome::quality_inserted   std::size_t, DEM nodes the pass added
//   RefineOutcome::quality_skipped    std::size_t, the pass's skips, a total
//   RefineOutcome::quality_seconds    double, wall seconds, >= 0
//
// CHOSEN HERE: RefineOutcome::inserted stays DEM refinement's own count, so
// vertices = start + quality_inserted + inserted. And Q7's increment-18
// reference is the topology digest in support/refine_digest.hpp, recorded from
// increment 18's refine.hpp (edde638's tree) before any production change.
//
// The tolerance oracle is 16's off-node one (prop_refinement_refine.cpp,
// check_offnode), restated here: every DEM node against every output triangle
// with three valid vertices, membership by DefaultKernel::orient2d on
// (col, -row), the plane from barycentric weights; constrained Delaunay by
// DefaultKernel::incircle in the frame (col * dx, -(row * dy)).

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/mesh/quality.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/refinement/refine.hpp>

#include "quality_fixtures.hpp"
#include "refine_digest.hpp"
#include "refinement_fixtures.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <utility>
#include <vector>

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::pred::DefaultKernel;
using terrain::pred::Orientation;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::refine;
using terrain::refinement::RefineOptions;
using namespace quality_fixtures;

namespace {

constexpr double kTheta = 25.0;

// min_angle_deg absent: the increment-18 call, byte for byte.
template <typename Mesh>
auto run18(const Raster<float>& dem, const Mesh& s, double tol, unsigned threads = 1) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    return refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks}, o);
}

template <typename Mesh>
auto run(const Raster<float>& dem, const Mesh& s, double tol, double min_angle, unsigned threads = 1) {
    RefineOptions o;
    o.tolerance = tol;
    o.threads = threads;
    o.min_angle_deg = min_angle;
    return refine(dem, s.mesh, std::span<const std::array<std::uint32_t, 2>>{s.edges},
                  std::span<const std::uint32_t>{s.masks}, o);
}

struct Frac {
    double col;
    double row;
};

Frac frac(const RasterGeometry& g, Point2 p) {
    return Frac{(p.x - g.x_min()) / g.delta_x(), (g.y_max() - p.y) / g.delta_y()};
}

bool is_node(const RasterGeometry& g, Point2 p) {
    const Frac f = frac(g, p);
    const double c = std::round(f.col), r = std::round(f.row);
    if (!(c >= 0 && r >= 0 && c < static_cast<double>(g.cols()) && r < static_cast<double>(g.rows())))
        return false;
    return g.node({static_cast<std::size_t>(r), static_cast<std::size_t>(c)}) == p;
}

bool nodata(const Raster<float>& dem, std::int64_t r, std::int64_t c) {
    return dem.is_nodata({static_cast<std::size_t>(r), static_cast<std::size_t>(c)});
}

double at(const Raster<float>& dem, std::int64_t r, std::int64_t c) {
    return static_cast<double>(dem.value_at({static_cast<std::size_t>(r), static_cast<std::size_t>(c)}));
}

std::optional<double> expected_z(const Raster<float>& dem, Point2 p) {
    const RasterGeometry& g = dem.geometry();
    if (is_node(g, p)) {
        const Frac f = frac(g, p);
        const auto r = static_cast<std::int64_t>(std::round(f.row)), c = static_cast<std::int64_t>(std::round(f.col));
        if (nodata(dem, r, c)) return std::nullopt;
        return at(dem, r, c);
    }
    return terrain::raster::bilinear(dem, p);
}

// Start vertices as given, every other vertex a node (quality's and
// refinement's alike), z per R0, the tolerance oracle, constrained Delaunay
// in the frame, and the vertex count split between the two inserters.
template <typename Outcome>
void check(const Raster<float>& dem, const IndexedMesh2& start, const Outcome& out, double tol) {
    REQUIRE(out.ok());
    const RasterGeometry& g = dem.geometry();
    const std::size_t n0 = start.vertices().size();
    REQUIRE(out.vertices.size() == n0 + out.quality_inserted + out.inserted);
    REQUIRE(out.quality_seconds >= 0.0);
    std::vector<Frac> f;
    double zmax = 0.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        CAPTURE(i);
        const Point2 p = out.vertices[i];
        if (i < n0) {
            REQUIRE(p.x == start.vertices()[i].x);
            REQUIRE(p.y == start.vertices()[i].y);
        } else {
            REQUIRE(is_node(g, p));  // R3: a Steiner point is a DEM node
        }
        const auto z = expected_z(dem, p);
        REQUIRE(static_cast<bool>(out.valid[i]) == z.has_value());
        REQUIRE(out.z[i] == z.value_or(0.0));
        zmax = std::max(zmax, std::abs(out.z[i]));
        f.push_back(frac(g, p));
    }
    auto fp = [&](std::uint32_t i) { return Point2{f[i].col, -f[i].row}; };
    auto cross = [](Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); };
    double worst = 0.0;
    for (const auto& tri : out.triangles) {
        const Point2 a = fp(tri[0]), b = fp(tri[1]), c = fp(tri[2]);
        REQUIRE(DefaultKernel::orient2d(a, b, c) == Orientation::CounterClockwise);
        if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
        const double two_a = cross(a, b, c);
        const auto lo_c = static_cast<std::int64_t>(std::ceil(std::min({a.x, b.x, c.x})));
        const auto hi_c = static_cast<std::int64_t>(std::floor(std::max({a.x, b.x, c.x})));
        const auto lo_r = static_cast<std::int64_t>(std::ceil(-std::max({a.y, b.y, c.y})));
        const auto hi_r = static_cast<std::int64_t>(std::floor(-std::min({a.y, b.y, c.y})));
        for (std::int64_t r = std::max<std::int64_t>(lo_r, 0); r <= hi_r; ++r)
            for (std::int64_t col = std::max<std::int64_t>(lo_c, 0); col <= hi_c; ++col) {
                const Point2 p{static_cast<double>(col), -static_cast<double>(r)};
                if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                    continue;
                if (p == a || p == b || p == c || nodata(dem, r, col)) continue;
                const double plane = (cross(p, b, c) * out.z[tri[0]] + cross(a, p, c) * out.z[tri[1]]
                                      + cross(a, b, p) * out.z[tri[2]]) / two_a;
                const double err = std::abs(at(dem, r, col) - plane);
                worst = std::max(worst, err);
                CAPTURE(r, col, err, tol);
                REQUIRE(err <= tol + 1e-9 * std::max(1.0, zmax));
            }
    }
    REQUIRE(out.max_error <= tol);
    REQUIRE(std::abs(out.max_error - worst) <= 1e-9 * std::max(1.0, zmax));

    auto lf = [&](std::uint32_t i) { return Point2{f[i].col * g.delta_x(), -(f[i].row * g.delta_y())}; };
    std::set<std::pair<std::uint32_t, std::uint32_t>> constrained;
    for (const auto& e : out.edges) constrained.insert(std::minmax(e[0], e[1]));
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> sides;
    for (std::size_t t = 0; t < out.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            sides[std::minmax(out.triangles[t][k], out.triangles[t][(k + 1) % 3])].push_back(t);
    std::size_t bad = 0;
    for (const auto& [e, ts] : sides) {
        REQUIRE(ts.size() <= 2);
        if (ts.size() != 2 || constrained.count(e) != 0) continue;
        for (unsigned s = 0; s < 2; ++s) {
            const auto& tri = out.triangles[ts[s]];
            std::uint32_t apex = 0;
            for (const auto x : out.triangles[ts[1 - s]])
                if (x != e.first && x != e.second) apex = x;
            if (DefaultKernel::incircle(lf(tri[0]), lf(tri[1]), lf(tri[2]), lf(apex)) == terrain::pred::Incircle::Inside)
                ++bad;
        }
    }
    REQUIRE(bad == 0);
    for (const auto m : out.masks) REQUIRE(m == 1u);  // the ring's mask, on every piece
}

template <typename Outcome>
void require_identical(const Outcome& a, const Outcome& b) {
    REQUIRE(refine_digest::digest(a) == refine_digest::digest(b));
    REQUIRE(a.quality_inserted == b.quality_inserted);
    REQUIRE(a.quality_skipped == b.quality_skipped);
}

}  // namespace

// ------------------------------------------------- 14b T3 and T6, quality on

TEST_CASE("T3 with quality: off-node ring fans refine to tolerance and stay constrained Delaunay",
          "[refinement][quality]") {
    const bool rough = GENERATE(false, true);
    const double tol = GENERATE(0.0, 0.5, 5.0);
    const std::uint32_t seed = GENERATE(1u, 2u, 3u, 4u);
    CAPTURE(rough, tol, seed);
    const std::size_t n = 33;
    const Raster<float> dem{refinement_fixtures::geometry(n, n),
                            rough ? refinement_fixtures::rough_dem(n, n, seed)
                                  : refinement_fixtures::smooth_dem(n, n, seed)};
    // A dense ring of 40 to 55 off-node vertices: its fan is all slivers.
    const auto start = fan(dem.geometry(), circle_ring(16.37, 16.61, 14.3, 40 + 5 * seed, 0.01 * seed));
    const auto out = run(dem, start, tol, kTheta);
    check(dem, start.mesh, out, tol);
    REQUIRE(out.quality_inserted > 0);
}

TEST_CASE("T6 with quality: bit-identical for 1 2 7 and all threads", "[refinement][quality]") {
    const std::size_t n = 33;
    const Raster<float> dem{refinement_fixtures::geometry(n, n), refinement_fixtures::rough_dem(n, n, 11)};
    const auto start = fan(dem.geometry(), circle_ring(16.21, 15.87, 15.1, 48));
    const auto ref = run(dem, start, 3.0, kTheta, 1);
    REQUIRE(ref.ok());
    REQUIRE(ref.quality_inserted > 0);
    REQUIRE(ref.inserted > 0);
    const unsigned threads = GENERATE(2u, 7u, 0u);
    CAPTURE(threads);
    require_identical(run(dem, start, 3.0, kTheta, threads), ref);
}

// ------------------------------------------------------------------------ Q7

TEST_CASE("Q7: min_angle_deg 0 reproduces increment 18 on the T12 fixtures and a domain start",
          "[refinement][quality]") {
    // Topology digests recorded from increment 18's refine.hpp; see the header.
    struct Case {
        const char* name;
        std::uint64_t golden;
    };
    const auto c = GENERATE(Case{"cone", 0x2326c8c0a462ac4dull}, Case{"island", 0x39c68151cb4e1186ull},
                            Case{"ring", 0x2b9d05a87b859d47ull});
    CAPTURE(c.name);
    const std::string name = c.name;
    auto both = [&](const auto& dem, const auto& start, double tol) {
        const auto before = run18(dem, start, tol);
        const auto off = run(dem, start, tol, 0.0);
        REQUIRE(refine_digest::topology_digest(before) == c.golden);
        require_identical(off, before);  // every double too, inside this binary
        REQUIRE(off.quality_inserted == 0);
        REQUIRE(off.quality_skipped == 0);
    };
    if (name == "ring") {
        const std::size_t n = 33;
        const Raster<float> dem{refinement_fixtures::geometry(n, n), refinement_fixtures::rough_dem(n, n, 11)};
        both(dem, fan(dem.geometry(), q7_ring()), 3.0);
    } else {
        const auto dem = name == "island" ? island(129, 2.0, 2.0) : cone(129, 2.0, 2.0);
        both(dem, refinement_fixtures::grid_mesh(dem.geometry(), 16), 1.0);
    }
}

TEST_CASE("Q7: the default RefineOptions leave the pass off", "[refinement][quality]") {
    REQUIRE(RefineOptions{}.min_angle_deg == 0.0);
}

TEST_CASE("Q7: the pass on changes the domain start's output", "[refinement][quality]") {
    // The converse, so Q7 cannot pass with a pass that never runs.
    const std::size_t n = 33;
    const Raster<float> dem{refinement_fixtures::geometry(n, n), refinement_fixtures::rough_dem(n, n, 11)};
    const auto start = fan(dem.geometry(), q7_ring());
    const auto on = run(dem, start, 3.0, kTheta);
    REQUIRE(on.ok());
    REQUIRE(on.quality_inserted > 0);
    REQUIRE(refine_digest::topology_digest(on) != 0x2b9d05a87b859d47ull);
    check(dem, start.mesh, on, 3.0);
}

// ------------------------------------------------------------------------ Q8

TEST_CASE("Q8: a stride start with square cells gets no quality node through refine", "[refinement][quality]") {
    const std::size_t n = 129;
    const auto dem = cone(n, 2.0, 2.0);
    const auto start = refinement_fixtures::grid_mesh(dem.geometry(), 16);
    const auto on = run(dem, start, 1.0, kTheta);
    REQUIRE(on.ok());
    REQUIRE(on.quality_inserted == 0);
    REQUIRE(on.quality_skipped == 0);
    require_identical(on, run18(dem, start, 1.0));
}

// ----------------------------------------------------------------------- Q10

TEST_CASE("Q10: a domain over a NoData block terminates and its NoData vertices are invalid",
          "[refinement][quality]") {
    // R9: the pass reads no height, so it may insert a NoData node; that vertex
    // is invalid, its triangles void, and refine still ends. The trim that
    // drops them is Python's (tests/python/test_cli_start_quality.py).
    const std::size_t n = 33;
    const bool sentinel = GENERATE(false, true);
    CAPTURE(sentinel);
    auto z = refinement_fixtures::rough_dem(n, n, 3);
    const float hole = sentinel ? -32767.0f : std::numeric_limits<float>::quiet_NaN();
    for (std::size_t r = 10; r <= 22; ++r)
        for (std::size_t c = 10; c <= 22; ++c) z[r * n + c] = hole;
    const RasterGeometry g0 = refinement_fixtures::geometry(n, n);
    const Raster<float> dem{g0, std::move(z), sentinel ? std::optional<float>{-32767.0f} : std::nullopt};
    const auto start = fan(dem.geometry(), circle_ring(16.37, 16.61, 14.3, 48));
    // The ring is (nearly) cocircular about (16.37, 16.61), so the first
    // circumcentre snaps to node (16, 17), inside the hole.
    const auto out = run(dem, start, 1.0, kTheta);
    check(dem, start.mesh, out, 1.0);
    REQUIRE(out.quality_inserted > 0);
    REQUIRE(std::count(out.valid.begin(), out.valid.end(), std::uint8_t{0}) > 0);
}
