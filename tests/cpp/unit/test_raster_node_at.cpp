// Increment 27 (docs/increments/27-node-sampling.md, "Where it lives" and
// "Tests for @tester", S5): RasterGeometry::node_at, the node a point is, bit
// for bit, or nullopt.
//
// S5 is the invariant-critical case for "the tolerance path is unchanged":
// refine writes a start vertex's output z with raster::bilinear only when that
// vertex is not a node by its own test (refine.hpp, the output loop:
// `!given || (node && g.node(c) == p) ? vertex_z : bilinear`). If node_at
// engaged on any point refine sends to bilinear, the new early return would
// change refine's output. So node_at is checked against that predicate --
// the producer's relation, lattice_position followed by the output loop's
// `g.node(c) == p` -- not against lattice_position's is_node() alone, which
// is also true for an off-node point whose fractional (col, row) happen to
// round to integers.
//
// Its own target: it includes refine.hpp, which names DefaultKernel and runs
// threads, so it links the backend target and Threads; test_raster and
// test_raster_view stay free of both.

#include <catch2/catch_test_macros.hpp>

#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/refine.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <random>
#include <string_view>
#include <vector>

using terrain::Point2;
using terrain::raster::CellIndex;
using terrain::raster::RasterGeometry;

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

struct NamedGrid {
    const char* name;
    RasterGeometry g;
};

// Integer, dyadic and non-dyadic geometry. The non-dyadic one is where
// RasterGeometry::node may be contracted into an FMA and where cell_of can
// put a node into the cell before it ("Why this test").
std::vector<NamedGrid> grids() {
    return {
        {"integer", RasterGeometry{799750.0, 7950250.0, 10.0, 10.0, 41, 37}},
        {"dyadic", RasterGeometry{500000.5, 7900000.25, 0.5, 0.25, 41, 37}},
        {"non-dyadic", RasterGeometry{500000.3, 7900000.7, 0.7, 0.3, 41, 37}},
        // "A limit of the exact rule": origin and spacing 0.1, where NumPy's
        // unfused node expression and a fused RasterGeometry::node differ in
        // the last bit at about a third of the nodes.
        {"tenths", RasterGeometry{0.1, 100.1, 0.1, 0.1, 120, 120}},
    };
}

// What refine does with a start vertex at p, reduced to the one question that
// matters here: does its output z come from the node (vertex_z), i.e. not
// from raster::bilinear? Copied from refine's output loop, with
// lattice_position called as to_lattice calls it.
std::optional<CellIndex> refine_reads_as_node(const RasterGeometry& g, const Point2& p) {
    const terrain::mesh::MeshVertex v = terrain::refinement::detail::lattice_position(g, p);
    if (!v.is_node())
        return std::nullopt;
    const CellIndex c{static_cast<std::size_t>(v.row), static_cast<std::size_t>(v.col)};
    if (!(g.node(c) == p))
        return std::nullopt;
    return c;
}

// Node (r, c) by NumPy's expression, never fused: the product is forced
// through memory, so `x_min + product` is two roundings (subsample builds
// stride vertices this way, grid_domain.py:68-69).
Point2 unfused_node(const RasterGeometry& g, std::size_t r, std::size_t c) {
    volatile double cx = static_cast<double>(c) * g.delta_x();
    volatile double ry = static_cast<double>(r) * g.delta_y();
    return Point2{g.x_min() + cx, g.y_max() - ry};
}

// Finite points in the node rectangle (cell_of's test): every node, the
// unfused version of it, the four nextafter neighbours of each node,
// cell-side points on every node line, and seeded random points.
std::vector<Point2> probes(const RasterGeometry& g, unsigned seed) {
    std::mt19937 rng{seed};
    std::uniform_real_distribution<double> ux{g.x_min(), g.x_max()};
    std::uniform_real_distribution<double> uy{g.y_min(), g.y_max()};
    std::vector<Point2> pts;
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c) {
            const Point2 at = g.node(CellIndex{r, c});
            pts.push_back(at);
            pts.push_back(unfused_node(g, r, c));
            pts.push_back(Point2{std::nextafter(at.x, kInf), at.y});
            pts.push_back(Point2{std::nextafter(at.x, -kInf), at.y});
            pts.push_back(Point2{at.x, std::nextafter(at.y, kInf)});
            pts.push_back(Point2{at.x, std::nextafter(at.y, -kInf)});
            pts.push_back(Point2{at.x, uy(rng)});
            pts.push_back(Point2{ux(rng), at.y});
        }
    for (int i = 0; i < 5000; ++i)
        pts.push_back(Point2{ux(rng), uy(rng)});
    std::vector<Point2> inside;
    for (const Point2& p : pts)
        if (g.cell_of(p))
            inside.push_back(p);
    return inside;
}

} // namespace

TEST_CASE("node_at: every node is itself, at its own index", "[raster][geometry][increment27]") {
    for (const auto& [name, g] : grids()) {
        CAPTURE(name);
        for (std::size_t r = 0; r < g.rows(); ++r)
            for (std::size_t c = 0; c < g.cols(); ++c) {
                CAPTURE(r, c);
                const auto n = g.node_at(g.node(CellIndex{r, c}));
                REQUIRE(n.has_value());
                REQUIRE(*n == CellIndex{r, c});
            }
    }
}

TEST_CASE("node_at: one ulp from a node is not a node", "[raster][geometry][increment27]") {
    for (const auto& [name, g] : grids()) {
        CAPTURE(name);
        const Point2 at = g.node(CellIndex{3, 4});
        for (const Point2 p : {Point2{std::nextafter(at.x, kInf), at.y},
                               Point2{std::nextafter(at.x, -kInf), at.y},
                               Point2{at.x, std::nextafter(at.y, kInf)},
                               Point2{at.x, std::nextafter(at.y, -kInf)}}) {
            CAPTURE(p.x, p.y);
            REQUIRE_FALSE(g.node_at(p).has_value());
        }
    }
}

TEST_CASE("node_at: a point on a cell side is not a node", "[raster][geometry][increment27]") {
    for (const auto& [name, g] : grids()) {
        CAPTURE(name);
        const Point2 at = g.node(CellIndex{3, 4});
        REQUIRE_FALSE(g.node_at(Point2{at.x + 0.5 * g.delta_x(), at.y}).has_value());
        REQUIRE_FALSE(g.node_at(Point2{at.x, at.y - 0.5 * g.delta_y()}).has_value());
    }
}

TEST_CASE("node_at: NaN, infinities and points outside the node rectangle are nullopt",
          "[raster][geometry][edge][increment27]") {
    for (const auto& [name, g] : grids()) {
        CAPTURE(name);
        const Point2 nw = g.node(CellIndex{0, 0});
        const Point2 se = g.node(CellIndex{g.rows() - 1, g.cols() - 1});
        for (const Point2 p : {
                 Point2{kNaN, nw.y}, Point2{nw.x, kNaN}, Point2{kNaN, kNaN},
                 Point2{kInf, nw.y}, Point2{-kInf, nw.y}, Point2{nw.x, kInf}, Point2{nw.x, -kInf},
                 // Outside by one ulp and by one whole cell, on node lines, so
                 // a clamping implementation would round them onto a corner.
                 Point2{std::nextafter(nw.x, -kInf), nw.y}, Point2{nw.x, std::nextafter(nw.y, kInf)},
                 Point2{std::nextafter(se.x, kInf), se.y}, Point2{se.x, std::nextafter(se.y, -kInf)},
                 Point2{nw.x - g.delta_x(), nw.y}, Point2{nw.x, nw.y + g.delta_y()},
                 Point2{se.x + g.delta_x(), se.y}, Point2{se.x, se.y - g.delta_y()}}) {
            CAPTURE(p.x, p.y);
            REQUIRE_FALSE(g.node_at(p).has_value());
        }
    }
}

TEST_CASE("S5: node_at agrees with refine's node test on every probe",
          "[raster][geometry][refine][invariant][increment27]") {
    unsigned seed = 2700;
    for (const auto& [name, g] : grids()) {
        CAPTURE(name);
        const auto pts = probes(g, seed++);
        REQUIRE(pts.size() > g.size() * 5);  // the probe set was built, mostly inside
        std::size_t nodes = 0, integral_but_off_node = 0;
        for (const Point2& p : pts) {
            CAPTURE(p.x, p.y);
            const auto mine = g.node_at(p);
            const auto refines = refine_reads_as_node(g, p);
            // lattice_position's fractional (col, row) can both be integers at
            // a point that is not node(c) bit for bit. refine sends such a
            // point to bilinear, so node_at must say nullopt there too.
            if (!refines && terrain::refinement::detail::lattice_position(g, p).is_node())
                ++integral_but_off_node;
            REQUIRE(mine.has_value() == refines.has_value());
            if (mine) {
                REQUIRE(*mine == *refines);
                ++nodes;
            }
        }
        // Every node probe is a node to both; the probe could fail.
        REQUIRE(nodes >= g.size());
        // On the tenths grid the probe set holds points that separate refine's
        // test from lattice_position(...).is_node() alone, so an
        // is_node()-shaped node_at fails here. Measured at 120 x 120 with
        // Apple clang: 1300 such probes; 350 of them nextafter neighbours that
        // remain under -ffp-contract=off, so this does not rest on FMA.
        if (std::string_view{name} == "tenths")
            REQUIRE(integral_but_off_node > 0);
    }
}
