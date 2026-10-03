// Increment 15f-1 (docs/increments/15f-edge-strip.md, D2, D3 and "Tests for
// @tester", CC1 to CC6): the edge-strip generator
//
//   constraint_check_points(dem, std::span<const Point2> vertices,
//                           std::span<const std::array<std::uint32_t, 2>> edges)
//       -> ConstraintCheckPoints
//
// and its store. Per edge, from its lower vertex index P0 to its higher P1:
// every crossing with a column or row line strictly between the ends (a node
// crossing once, at the node exactly), and the midpoint between each two
// neighbouring points of the list with the ends added. Each point carries
// vertex_z there and its parameter s along P0 -> P1.
//
// What this suite CHOOSES where D2/D3 are silent, so @developer can see it:
//   - edge_count() equals the number of edges given, and each given edge
//     appears once as edge(k), lower index first. The suite never assumes
//     which k an edge lands on; it looks edges up by their pair.
//   - the refusals of D2 (an end outside the node rectangle, NaN or infinite
//     coordinates among them, and an index out of range) are
//     std::invalid_argument whose what() contains "constraint_check_points",
//     after refine's "refine: ..." house style.
//
// Frame. Lattice (col, row), with world x = x_min + col dx and
// y = y_max - row dy, so rows grow downward. Cells are 10 m by 5 m at UTM
// scale, so a row/column swap cannot pass. In CC1 to CC5 every end is a dyadic
// fraction of a cell and every DEM value an integer, so every position, s and
// z is exact and compared with ==. The property cases (random ends, extreme
// scales) compare z with a relative margin and positions by structure.
//
// The expected heights are written here from Q14's words (linear between the
// two nodes of a cell side; bilinear in the cell that holds a midpoint;
// NoData when any corner of that cell is NoData), never by calling vertex_z.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

using terrain::Point2;
using terrain::mesh::MeshVertex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::refinement::ConstraintCheckPoints;
using terrain::refinement::ConstraintPoint;
using terrain::refinement::constraint_check_points;

namespace {

using Edge = std::array<std::uint32_t, 2>;

constexpr double kX0 = 500000.0, kY0 = 7000000.0, kDx = 10.0, kDy = 5.0;
constexpr std::size_t kCols = 9, kRows = 7;

RasterGeometry geometry(std::size_t cols, std::size_t rows) {
    return RasterGeometry{kX0, kY0, kDx, kDy, cols, rows};
}

// Exact for the dyadic (col, row) of CC1 to CC5.
Point2 world(double col, double row) { return Point2{kX0 + col * kDx, kY0 - row * kDy}; }

// The DEM: integers, so linear and bilinear values at dyadic positions are
// exact; not linear in (row, col), so a crossing's height depends on which two
// nodes it is taken between.
double h(std::size_t row, std::size_t col) {
    return 1000.0 + 7.0 * static_cast<double>(col) + 13.0 * static_cast<double>(row)
         + 4.0 * static_cast<double>((5 * row + 3 * col) % 11);
}

struct NoData {
    std::size_t row;
    std::size_t col;
};

enum class Void { Sentinel, NaN };

constexpr float kSentinel = -9999.0f;

Raster<float> dem(std::size_t cols, std::size_t rows, std::vector<NoData> holes = {},
                  Void how = Void::Sentinel) {
    std::vector<float> data(cols * rows);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c)
            data[r * cols + c] = static_cast<float>(h(r, c));
    for (const auto& n : holes)
        data[n.row * cols + n.col] =
            how == Void::NaN ? std::numeric_limits<float>::quiet_NaN() : kSentinel;
    return how == Void::NaN ? Raster<float>{geometry(cols, rows), std::move(data)}
                            : Raster<float>{geometry(cols, rows), std::move(data), kSentinel};
}

// Q14: a crossing of column line K at fractional row r is linear between the
// nodes (floor r, K) and (floor r + 1, K); a crossing of row line R likewise.
double on_column(std::size_t k, double row) {
    const auto r0 = static_cast<std::size_t>(std::floor(row));
    const double f = row - static_cast<double>(r0);
    return h(r0, k) * (1.0 - f) + h(r0 + 1, k) * f;
}
double on_row(double col, std::size_t r) {
    const auto c0 = static_cast<std::size_t>(std::floor(col));
    const double f = col - static_cast<double>(c0);
    return h(r, c0) * (1.0 - f) + h(r, c0 + 1) * f;
}

// Bilinear in the cell holding (col, row), the last cell on the far sides; a
// node's own value at a node.
double bilinear(double col, double row, std::size_t cols = kCols, std::size_t rows = kRows) {
    if (col == std::floor(col) && row == std::floor(row))
        return h(static_cast<std::size_t>(row), static_cast<std::size_t>(col));
    const std::size_t r0 = std::min(static_cast<std::size_t>(row), rows - 2);
    const std::size_t c0 = std::min(static_cast<std::size_t>(col), cols - 2);
    const double ty = row - static_cast<double>(r0), tx = col - static_cast<double>(c0);
    return h(r0, c0) * (1 - tx) * (1 - ty) + h(r0, c0 + 1) * tx * (1 - ty)
         + h(r0 + 1, c0) * (1 - tx) * ty + h(r0 + 1, c0 + 1) * tx * ty;
}

struct Expected {
    double col;
    double row;
    double z;
    double s;
};

ConstraintCheckPoints generate(const Raster<float>& d, const std::vector<Point2>& vertices,
                               const std::vector<Edge>& edges) {
    return constraint_check_points(d, std::span<const Point2>{vertices},
                                   std::span<const Edge>{edges});
}

// The k under which `e` was filed, whatever order its pair was given in.
std::size_t index_of(const ConstraintCheckPoints& st, Edge e) {
    const Edge want{std::min(e[0], e[1]), std::max(e[0], e[1])};
    for (std::size_t k = 0; k < st.edge_count(); ++k)
        if (st.edge(k) == want)
            return k;
    FAIL("edge (" << want[0] << ", " << want[1] << ") is not in the store");
    return 0;
}

void check_exact(std::span<const ConstraintPoint> got, const std::vector<Expected>& want) {
    REQUIRE(got.size() == want.size());
    for (std::size_t i = 0; i < want.size(); ++i) {
        INFO("point " << i << ": expected (" << want[i].col << ", " << want[i].row << ") z "
                      << want[i].z << " s " << want[i].s << "; got (" << got[i].at.col << ", "
                      << got[i].at.row << ") z " << got[i].z << " s " << got[i].s);
        CHECK(got[i].at.col == want[i].col);
        CHECK(got[i].at.row == want[i].row);
        CHECK(got[i].z == want[i].z);
        CHECK(got[i].s == want[i].s);
    }
}

std::uint64_t bits(double v) { return std::bit_cast<std::uint64_t>(v); }

bool same_bits(const ConstraintPoint& a, const ConstraintPoint& b) {
    return bits(a.at.col) == bits(b.at.col) && bits(a.at.row) == bits(b.at.row)
        && bits(a.z) == bits(b.z) && bits(a.s) == bits(b.s);
}

// Integers strictly between a and b.
std::size_t strictly_between(double a, double b) {
    const double lo = std::min(a, b), hi = std::max(a, b);
    const double n = std::ceil(hi) - std::floor(lo) - 1.0;
    return n > 0.0 ? static_cast<std::size_t>(n) : 0;
}

bool integral(double v) { return v == std::floor(v); }

// D2's shape, on one edge of a store built on a DEM with no NoData, given the
// lattice ends of edge(k) in canonical order (P0 = the lower index):
//   - an odd count, midpoint and crossing alternating, midpoints at both ends;
//   - s strictly increasing inside (0, 1);
//   - a crossing lies on a grid line, inside the node rectangle, within 1e-9
//     cells of the segment, is neither end, and its s is its parameter;
//   - a midpoint is ((a.col + b.col)/2, (a.row + b.row)/2) of its neighbours,
//     with s (s_a + s_b)/2, the ends counting as neighbours at s = 0 and 1;
//   - z is the bilinear height there, to a relative 1e-12.
// Returns the number of crossings.
std::size_t check_structure(const ConstraintCheckPoints& st, std::size_t k, MeshVertex a,
                            MeshVertex b, std::size_t cols, std::size_t rows) {
    const auto pts = st.on_edge(k);
    INFO("edge " << k << " (" << st.edge(k)[0] << ", " << st.edge(k)[1] << ") from (" << a.col
                 << ", " << a.row << ") to (" << b.col << ", " << b.row << "), " << pts.size()
                 << " points");
    REQUIRE(pts.size() % 2 == 1);
    const double dc = b.col - a.col, dr = b.row - a.row;
    const double len2 = dc * dc + dr * dr;
    double prev_s = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const ConstraintPoint& p = pts[i];
        INFO("point " << i << " at (" << p.at.col << ", " << p.at.row << ") s " << p.s);
        REQUIRE(p.s > prev_s);
        REQUIRE(p.s < 1.0);
        prev_s = p.s;
        CHECK(p.z == Catch::Approx(bilinear(p.at.col, p.at.row, cols, rows)).epsilon(1e-12));
        if (i % 2 == 1) {
            CHECK((integral(p.at.col) || integral(p.at.row)));
            CHECK(p.at.col >= 0.0);
            CHECK(p.at.col <= static_cast<double>(cols - 1));
            CHECK(p.at.row >= 0.0);
            CHECK(p.at.row <= static_cast<double>(rows - 1));
            CHECK_FALSE(p.at == a);
            CHECK_FALSE(p.at == b);
            const double pc = p.at.col - a.col, pr = p.at.row - a.row;
            CHECK(std::abs(dc * pr - dr * pc) / std::sqrt(len2) <= 1e-9);
            CHECK(p.s == Catch::Approx((pc * dc + pr * dr) / len2).margin(1e-9));
        } else {
            const MeshVertex l = i == 0 ? a : pts[i - 1].at;
            const MeshVertex r = i + 1 == pts.size() ? b : pts[i + 1].at;
            const double sl = i == 0 ? 0.0 : pts[i - 1].s;
            const double sr = i + 1 == pts.size() ? 1.0 : pts[i + 1].s;
            CHECK(bits(p.at.col) == bits((l.col + r.col) / 2));
            CHECK(bits(p.at.row) == bits((l.row + r.row) / 2));
            CHECK(bits(p.s) == bits((sl + sr) / 2));
        }
    }
    return pts.size() / 2;
}

// The generator's ends, by the function D2 step 1 names (pinned on its own by
// test_refinement_lattice_position.cpp).
MeshVertex lattice(const RasterGeometry& g, Point2 p) {
    return terrain::refinement::detail::lattice_position(g, p);
}

// One edge, by hand: P0 = vertex 0 = (0.5, 0.5), P1 = vertex 1 = (4.5, 2.5).
// Column crossings at K = 1..4 (s = 1/8, 3/8, 5/8, 7/8), row crossings at
// R = 1, 2 (s = 1/4, 3/4); none at a node.
std::vector<Expected> cc1_expected() {
    return {
        {0.75, 0.625, bilinear(0.75, 0.625), 0.0625},
        {1.0, 0.75, on_column(1, 0.75), 0.125},
        {1.25, 0.875, bilinear(1.25, 0.875), 0.1875},
        {1.5, 1.0, on_row(1.5, 1), 0.25},
        {1.75, 1.125, bilinear(1.75, 1.125), 0.3125},
        {2.0, 1.25, on_column(2, 1.25), 0.375},
        {2.5, 1.5, bilinear(2.5, 1.5), 0.5},
        {3.0, 1.75, on_column(3, 1.75), 0.625},
        {3.25, 1.875, bilinear(3.25, 1.875), 0.6875},
        {3.5, 2.0, on_row(3.5, 2), 0.75},
        {3.75, 2.125, bilinear(3.75, 2.125), 0.8125},
        {4.0, 2.25, on_column(4, 2.25), 0.875},
        {4.25, 2.375, bilinear(4.25, 2.375), 0.9375},
    };
}

const std::vector<Point2> kCc1Vertices{world(0.5, 0.5), world(4.5, 2.5)};

// An edge along a grid line from fractional end to fractional end, with four
// node crossings at s = 1/8, 3/8, 5/8, 7/8: the nodes exactly, with their
// values, and the five side midpoints. `along` maps the parameter along the
// line (0.5 .. 4.5) to (col, row).
template <class Along>
std::vector<Expected> along_line_expected(Along along) {
    std::vector<Expected> out;
    auto add = [&](double u, double s) {
        const auto [c, r] = along(u);
        out.push_back({c, r, bilinear(c, r), s});
    };
    add(0.75, 0.0625);
    for (int n = 1; n <= 4; ++n) {
        add(static_cast<double>(n), (n - 0.5) / 4.0);
        if (n < 4)
            add(n + 0.5, n / 4.0);
    }
    add(4.25, 0.9375);
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// The store
// ---------------------------------------------------------------------------

TEST_CASE("ST1: a constraint point is 32 bytes", "[constraint_points][store]") {
    STATIC_REQUIRE(sizeof(ConstraintPoint) == 32);
}

TEST_CASE("ST2: no edges, no points; the store keeps the DEM's geometry",
          "[constraint_points][store]") {
    const auto d = dem(kCols, kRows);
    const auto st = generate(d, kCc1Vertices, {});
    CHECK(st.edge_count() == 0);
    CHECK(st.size() == 0);
    CHECK(st.no_data() == 0);
    CHECK(st.duplicates() == 0);
    const RasterGeometry& g = st.geometry();
    CHECK(g.x_min() == kX0);
    CHECK(g.y_max() == kY0);
    CHECK(g.delta_x() == kDx);
    CHECK(g.delta_y() == kDy);
    CHECK(g.cols() == kCols);
    CHECK(g.rows() == kRows);
}

TEST_CASE("ST3: every given edge is filed once, lower index first, and size() sums them",
          "[constraint_points][store]") {
    const auto d = dem(kCols, kRows);
    const std::vector<Point2> v{world(0.5, 0.5), world(4.5, 2.5), world(2.25, 5.5),
                                world(7.75, 0.25)};
    const std::vector<Edge> edges{{1, 0}, {1, 2}, {3, 2}, {0, 3}};
    const auto st = generate(d, v, edges);
    REQUIRE(st.edge_count() == edges.size());
    std::multiset<Edge> want, got;
    for (const auto& e : edges)
        want.insert({std::min(e[0], e[1]), std::max(e[0], e[1])});
    std::size_t total = 0;
    for (std::size_t k = 0; k < st.edge_count(); ++k) {
        CHECK(st.edge(k)[0] < st.edge(k)[1]);
        got.insert(st.edge(k));
        total += st.on_edge(k).size();
    }
    CHECK(got == want);
    CHECK(st.size() == total);
    CHECK(total > 0);
}

// ---------------------------------------------------------------------------
// CC1 to CC6
// ---------------------------------------------------------------------------

TEST_CASE("CC1: one off-node edge by hand: crossings, midpoints, z, s",
          "[constraint_points][CC1]") {
    const auto d = dem(kCols, kRows);
    // Given as (1, 0): the store files it as (0, 1) and measures s from vertex 0.
    const auto st = generate(d, kCc1Vertices, {{1, 0}});
    REQUIRE(st.edge_count() == 1);
    CHECK(st.edge(0) == Edge{0, 1});
    check_exact(st.on_edge(0), cc1_expected());
    CHECK(st.size() == 13);  // 2c + 1 for c = 6 crossings
    CHECK(st.duplicates() == 0);
    CHECK(st.no_data() == 0);
    for (const auto& p : st.on_edge(0)) {
        CHECK_FALSE(p.at == MeshVertex{0.5, 0.5});
        CHECK_FALSE(p.at == MeshVertex{4.5, 2.5});
    }
}

TEST_CASE("CC1: the lower index is P0 whichever end it is", "[constraint_points][CC1]") {
    // Vertex 0 is now the far end, so s runs the other way and the list reverses.
    const auto d = dem(kCols, kRows);
    const std::vector<Point2> v{kCc1Vertices[1], kCc1Vertices[0]};
    const auto st = generate(d, v, {{0, 1}});
    auto want = cc1_expected();
    std::reverse(want.begin(), want.end());
    for (auto& e : want)
        e.s = 1.0 - e.s;  // exact: every s is a multiple of 1/16
    check_exact(st.on_edge(0), want);
}

TEST_CASE("CC2: an edge along a column line gives the nodes and the side midpoints",
          "[constraint_points][CC2]") {
    const auto d = dem(kCols, kRows);
    SECTION("an interior column, 3") {
        const auto st = generate(d, {world(3, 0.5), world(3, 4.5)}, {{0, 1}});
        const auto want = along_line_expected([](double u) { return std::pair{3.0, u}; });
        check_exact(st.on_edge(0), want);
        CHECK(st.duplicates() == 0);
        // The side midpoints between nodes are the mean of their two nodes.
        CHECK(want[2].z == (h(1, 3) + h(2, 3)) / 2);
        for (std::size_t i = 1; i < st.on_edge(0).size(); i += 2)
            CHECK(st.on_edge(0)[i].at.is_node());
    }
    SECTION("the far column of the node rectangle, the domain outline's case") {
        const double far = static_cast<double>(kCols - 1);
        const auto st = generate(d, {world(far, 0.5), world(far, 4.5)}, {{0, 1}});
        check_exact(st.on_edge(0), along_line_expected([far](double u) {
                        return std::pair{far, u};
                    }));
        CHECK(st.on_edge(0)[0].z == on_column(kCols - 1, 0.75));
    }
}

TEST_CASE("CC2: an edge along a row line gives the nodes and the side midpoints",
          "[constraint_points][CC2]") {
    const auto d = dem(kCols, kRows);
    SECTION("an interior row, 2") {
        const auto st = generate(d, {world(0.5, 2), world(4.5, 2)}, {{0, 1}});
        check_exact(st.on_edge(0), along_line_expected([](double u) { return std::pair{u, 2.0}; }));
        CHECK(st.duplicates() == 0);
    }
    SECTION("the bottom row of the node rectangle") {
        const double last = static_cast<double>(kRows - 1);
        const auto st = generate(d, {world(0.5, last), world(4.5, last)}, {{0, 1}});
        check_exact(st.on_edge(0), along_line_expected([last](double u) {
                        return std::pair{u, last};
                    }));
        CHECK(st.on_edge(0)[0].z == on_row(0.75, kRows - 1));
    }
}

TEST_CASE("CC3: a diagonal through nodes: each node once, exactly; the merges counted",
          "[constraint_points][CC3]") {
    const auto d = dem(kCols, kRows);
    SECTION("(0.5, 0.5) to (3.5, 3.5) through (1, 1), (2, 2), (3, 3)") {
        const auto st = generate(d, {world(0.5, 0.5), world(3.5, 3.5)}, {{0, 1}});
        const auto pts = st.on_edge(0);
        REQUIRE(pts.size() == 7);
        CHECK(st.duplicates() == 3);
        const std::array<double, 7> at{0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.25};
        for (std::size_t i = 0; i < 7; ++i) {
            INFO("point " << i);
            CHECK(pts[i].at == MeshVertex{at[i], at[i]});
            CHECK(pts[i].z == bilinear(at[i], at[i]));
        }
        for (std::size_t n = 1; n <= 3; ++n) {
            CHECK(pts[2 * n - 1].at.is_node());
            CHECK(pts[2 * n - 1].z == h(n, n));
            CHECK(pts[2 * n - 1].s == Catch::Approx((n - 0.5) / 3.0).margin(1e-15));
        }
        check_structure(st, 0, {0.5, 0.5}, {3.5, 3.5}, kCols, kRows);
    }
    SECTION("(0.5, 0.5) to (3.5, 1.5), a shallow slope through the node (2, 1)") {
        const auto st = generate(d, {world(0.5, 0.5), world(3.5, 1.5)}, {{0, 1}});
        const auto pts = st.on_edge(0);
        REQUIRE(pts.size() == 7);  // column crossings K = 1, 2, 3; row R = 1 is K = 2
        CHECK(st.duplicates() == 1);
        CHECK(pts[3].at == MeshVertex{2.0, 1.0});
        CHECK(pts[3].z == h(1, 2));
        CHECK(std::count_if(pts.begin(), pts.end(), [](const ConstraintPoint& p) {
                  return p.at == MeshVertex{2.0, 1.0};
              }) == 1);
        check_structure(st, 0, {0.5, 0.5}, {3.5, 1.5}, kCols, kRows);
    }
    SECTION("a hair from the nodes: both crossings of each kept, none at a node") {
        // P1 moved one ulp of y (about 1.9e-10 rows): the line misses every
        // node, by far more than rounding, so nothing is snapped or merged.
        const Point2 p1 = world(3.5, 3.5);
        const std::vector<Point2> v{world(0.5, 0.5),
                                    {p1.x, std::nextafter(p1.y, -std::numeric_limits<double>::infinity())}};
        const auto st = generate(d, v, {{0, 1}});
        const auto pts = st.on_edge(0);
        CHECK(pts.size() == 13);
        CHECK(st.duplicates() == 0);
        for (const auto& p : pts)
            CHECK_FALSE(p.at.is_node());
        const auto g = d.geometry();
        check_structure(st, 0, lattice(g, v[0]), lattice(g, v[1]), kCols, kRows);
    }
}

TEST_CASE("CC4: the short, the parallel, and ends on grid lines", "[constraint_points][CC4]") {
    const auto d = dem(kCols, kRows);
    SECTION("an edge inside one cell has its one midpoint") {
        const auto st = generate(d, {world(2.25, 1.25), world(2.75, 1.5)}, {{0, 1}});
        check_exact(st.on_edge(0), {{2.5, 1.375, bilinear(2.5, 1.375), 0.5}});
        CHECK(st.duplicates() == 0);
    }
    SECTION("an edge shorter than a cell across one column line") {
        const auto st = generate(d, {world(0.75, 0.5), world(1.25, 0.5)}, {{0, 1}});
        check_exact(st.on_edge(0), {
                                       {0.875, 0.5, bilinear(0.875, 0.5), 0.25},
                                       {1.0, 0.5, on_column(1, 0.5), 0.5},
                                       {1.125, 0.5, bilinear(1.125, 0.5), 0.75},
                                   });
    }
    SECTION("a vertical edge between column lines, at col 2.25: row crossings only") {
        const auto st = generate(d, {world(2.25, 0.5), world(2.25, 4.5)}, {{0, 1}});
        auto want = along_line_expected([](double u) { return std::pair{2.25, u}; });
        for (std::size_t r = 1; r <= 4; ++r)
            CHECK(want[2 * r - 1].z == on_row(2.25, r));  // linear along the row side
        check_exact(st.on_edge(0), want);
        CHECK(st.duplicates() == 0);
    }
    SECTION("a horizontal edge between row lines, at row 3.75: column crossings only") {
        const auto st = generate(d, {world(0.5, 3.75), world(4.5, 3.75)}, {{0, 1}});
        auto want = along_line_expected([](double u) { return std::pair{u, 3.75}; });
        for (std::size_t c = 1; c <= 4; ++c)
            CHECK(want[2 * c - 1].z == on_column(c, 3.75));
        check_exact(st.on_edge(0), want);
    }
    SECTION("an end exactly on a column line has no point at that end") {
        const auto st = generate(d, {world(1.0, 0.25), world(3.0, 1.25)}, {{0, 1}});
        check_exact(st.on_edge(0), {
                                       {1.5, 0.5, bilinear(1.5, 0.5), 0.25},
                                       {2.0, 0.75, on_column(2, 0.75), 0.5},
                                       {2.25, 0.875, bilinear(2.25, 0.875), 0.625},
                                       {2.5, 1.0, on_row(2.5, 1), 0.75},
                                       {2.75, 1.125, bilinear(2.75, 1.125), 0.875},
                                   });
    }
    SECTION("an end exactly on a node has no point at that node") {
        const auto st = generate(d, {d.geometry().node({1, 1}), world(2.5, 2.5)}, {{0, 1}});
        const auto pts = st.on_edge(0);
        REQUIRE(pts.size() == 3);
        CHECK(pts[0].at == MeshVertex{1.5, 1.5});
        CHECK(pts[1].at == MeshVertex{2.0, 2.0});
        CHECK(pts[2].at == MeshVertex{2.25, 2.25});
        CHECK(st.duplicates() == 1);
        check_structure(st, 0, {1.0, 1.0}, {2.5, 2.5}, kCols, kRows);
    }
    SECTION("an edge between two nodes along a column line: interior nodes only") {
        const auto& g = d.geometry();
        const auto st = generate(d, {g.node({1, 4}), g.node({4, 4})}, {{0, 1}});
        const auto pts = st.on_edge(0);
        REQUIRE(pts.size() == 5);
        CHECK(pts[1].at == MeshVertex{4.0, 2.0});
        CHECK(pts[3].at == MeshVertex{4.0, 3.0});
        CHECK(pts[0].z == (h(1, 4) + h(2, 4)) / 2);
        check_structure(st, 0, {4.0, 1.0}, {4.0, 4.0}, kCols, kRows);
    }
}

TEST_CASE("CC5: points whose cell has a NoData corner are dropped and counted",
          "[constraint_points][CC5]") {
    const auto clean = generate(dem(kCols, kRows), kCc1Vertices, {{0, 1}});
    const auto clean_pts = clean.on_edge(0);
    REQUIRE(clean_pts.size() == 13);
    // NoData at node (row 3, col 4). It is a corner of the cells of the last
    // four points of CC1's edge: the row crossing at (3.5, 2) -- where that
    // corner has weight 0 -- the midpoint after it, the column crossing at
    // (4, 2.25) and the last midpoint. The first nine are untouched.
    for (const Void how : {Void::Sentinel, Void::NaN}) {
        INFO((how == Void::NaN ? "NaN" : "sentinel") << " NoData");
        const auto st = generate(dem(kCols, kRows, {{3, 4}}, how), kCc1Vertices, {{0, 1}});
        const auto pts = st.on_edge(0);
        CHECK(st.no_data() == 4);
        CHECK(st.duplicates() == 0);
        CHECK(st.size() == 9);
        REQUIRE(pts.size() == 9);
        for (std::size_t i = 0; i < 9; ++i) {
            INFO("point " << i);
            CHECK(same_bits(pts[i], clean_pts[i]));
        }
    }
}

TEST_CASE("CC5: a node crossing on a NoData node is dropped with the cells around it",
          "[constraint_points][CC5]") {
    // The diagonal of CC3 with node (2, 2) NoData: the node itself and the two
    // midpoints in cells that have it as a corner go; (1, 1), (3, 3) and the
    // end midpoints stay.
    const auto st =
        generate(dem(kCols, kRows, {{2, 2}}), {world(0.5, 0.5), world(3.5, 3.5)}, {{0, 1}});
    const auto pts = st.on_edge(0);
    CHECK(st.no_data() == 3);
    CHECK(st.duplicates() == 3);
    REQUIRE(pts.size() == 4);
    const std::array<double, 4> at{0.75, 1.0, 3.0, 3.25};
    for (std::size_t i = 0; i < 4; ++i) {
        INFO("point " << i);
        CHECK(pts[i].at == MeshVertex{at[i], at[i]});
        CHECK(pts[i].z == bilinear(at[i], at[i]));
    }
    CHECK(pts[0].s < pts[1].s);
    CHECK(pts[1].s < pts[2].s);
    CHECK(pts[2].s < pts[3].s);
}

namespace {

// splitmix64: a fixed, library-independent stream, so the fixture is the same
// on every standard library.
struct Rng {
    std::uint64_t state;
    std::uint64_t next() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double unit() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }
    std::uint32_t below(std::uint32_t n) { return static_cast<std::uint32_t>(next() % n); }
};

struct Fixture {
    std::size_t cols;
    std::size_t rows;
    std::vector<Point2> vertices;
    std::vector<Edge> edges;
};

// Off-node ends at non-dyadic positions, nodes, points on grid lines and on
// the rectangle's sides; edges at random, plus vertical and horizontal pairs,
// diagonals through nodes, and a fan sharing one vertex.
Fixture random_fixture(std::uint64_t seed) {
    Fixture f{40, 30, {}, {}};
    Rng rng{seed};
    const auto g = geometry(f.cols, f.rows);
    const double cmax = static_cast<double>(f.cols - 1), rmax = static_cast<double>(f.rows - 1);
    auto any = [&] { return Point2{kX0 + rng.unit() * cmax * kDx, kY0 - rng.unit() * rmax * kDy}; };
    for (int i = 0; i < 60; ++i)
        f.vertices.push_back(any());
    for (int i = 0; i < 10; ++i)
        f.vertices.push_back(g.node({rng.below(30), rng.below(40)}));
    for (int i = 0; i < 10; ++i) {  // on a column line, then on a row line
        const double c = static_cast<double>(rng.below(40));
        f.vertices.push_back({kX0 + c * kDx, kY0 - rng.unit() * rmax * kDy});
        const double r = static_cast<double>(rng.below(30));
        f.vertices.push_back({kX0 + rng.unit() * cmax * kDx, kY0 - r * kDy});
    }
    for (int i = 0; i < 4; ++i) {  // on the four sides of the node rectangle
        f.vertices.push_back(world(0.0, rng.unit() * rmax));
        f.vertices.push_back(world(cmax, rng.unit() * rmax));
        f.vertices.push_back(world(rng.unit() * cmax, 0.0));
        f.vertices.push_back(world(rng.unit() * cmax, rmax));
    }
    const auto n = static_cast<std::uint32_t>(f.vertices.size());
    std::set<Edge> seen;
    auto add = [&](std::uint32_t a, std::uint32_t b) {
        if (a != b && seen.insert({std::min(a, b), std::max(a, b)}).second)
            f.edges.push_back({a, b});
    };
    for (int i = 0; i < 150; ++i)
        add(rng.below(n), rng.below(n));
    // Exactly vertical and horizontal, at the same fractional column / row.
    for (int i = 0; i < 6; ++i) {
        const double c = rng.unit() * cmax, r = rng.unit() * rmax;
        const auto base = static_cast<std::uint32_t>(f.vertices.size());
        f.vertices.push_back({kX0 + c * kDx, kY0 - r * kDy});
        f.vertices.push_back({kX0 + c * kDx, kY0 - rng.unit() * rmax * kDy});
        f.vertices.push_back({kX0 + rng.unit() * cmax * kDx, kY0 - r * kDy});
        add(base, base + 1);
        add(base + 2, base);
    }
    // Diagonals through nodes, node end to node end.
    for (std::size_t i = 0; i < 4; ++i) {
        const auto base = static_cast<std::uint32_t>(f.vertices.size());
        f.vertices.push_back(g.node({2 + i, 3 + i}));
        f.vertices.push_back(g.node({2 + i + 7, 3 + i + 14}));
        add(base + 1, base);
    }
    // A fan from one node: the edges share only that vertex.
    const auto hub = static_cast<std::uint32_t>(f.vertices.size());
    f.vertices.push_back(g.node({15, 20}));
    for (std::uint32_t i = 0; i < 12; ++i)
        add(hub, rng.below(hub));
    return f;
}

}  // namespace

TEST_CASE("CC6: reversing every pair and permuting the edge list changes no point, bit for bit",
          "[constraint_points][CC6]") {
    for (const std::uint64_t seed : {1ull, 2ull, 0xC0FFEEull}) {
        INFO("seed " << seed);
        const auto f = random_fixture(seed);
        const auto d = dem(f.cols, f.rows);
        const auto a = generate(d, f.vertices, f.edges);

        auto edges = f.edges;
        for (auto& e : edges)
            std::swap(e[0], e[1]);
        Rng rng{seed ^ 0x5EEDull};
        for (std::size_t i = edges.size(); i > 1; --i)
            std::swap(edges[i - 1], edges[rng.below(static_cast<std::uint32_t>(i))]);
        const auto b = generate(d, f.vertices, edges);

        REQUIRE(a.edge_count() == b.edge_count());
        CHECK(a.size() == b.size());
        CHECK(a.duplicates() == b.duplicates());
        CHECK(a.no_data() == b.no_data());
        for (std::size_t k = 0; k < a.edge_count(); ++k) {
            const auto pa = a.on_edge(k);
            const auto pb = b.on_edge(index_of(b, a.edge(k)));
            INFO("edge (" << a.edge(k)[0] << ", " << a.edge(k)[1] << ")");
            REQUIRE(pa.size() == pb.size());
            for (std::size_t i = 0; i < pa.size(); ++i) {
                INFO("point " << i);
                CHECK(same_bits(pa[i], pb[i]));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// D2's shape over adversarial input
// ---------------------------------------------------------------------------

TEST_CASE("CC-P: every edge of a random fixture has D2's shape; crossings are all counted",
          "[constraint_points][property]") {
    for (const std::uint64_t seed : {1ull, 2ull, 0xC0FFEEull}) {
        INFO("seed " << seed);
        const auto f = random_fixture(seed);
        const auto d = dem(f.cols, f.rows);
        const auto& g = d.geometry();
        const auto st = generate(d, f.vertices, f.edges);
        CHECK(st.no_data() == 0);
        std::size_t crossings = 0, lines = 0;
        for (std::size_t k = 0; k < st.edge_count(); ++k) {
            const MeshVertex a = lattice(g, f.vertices[st.edge(k)[0]]);
            const MeshVertex b = lattice(g, f.vertices[st.edge(k)[1]]);
            crossings += check_structure(st, k, a, b, f.cols, f.rows);
            lines += strictly_between(a.col, b.col) + strictly_between(a.row, b.row);
        }
        // Every grid-line crossing is either a point or a counted merge.
        CHECK(crossings + st.duplicates() == lines);
    }
}

TEST_CASE("CC-P: edges sharing a vertex put no point on it, and no two of their points coincide",
          "[constraint_points][property]") {
    const auto d = dem(kCols, kRows);
    const auto& g = d.geometry();
    // A fan from the node (3, 4) and from the off-node vertex (2.5, 4.25), as
    // a shared vertex of a constrained outline and of interior polygons is.
    std::vector<Point2> v{g.node({3, 4}), world(2.5, 4.25), world(0.25, 0.75), world(7.5, 0.5),
                          world(8.0, 5.75), world(0.5, 6.0), world(4.0, 0.0), world(6.25, 3.0)};
    // Planar: these edges meet only at shared vertices, as a triangulation's
    // constraint edges do. {0, 7} runs along row line 3 from the node.
    const std::vector<Edge> edges{{0, 2}, {3, 0}, {0, 4}, {5, 0}, {0, 6}, {0, 7}, {1, 0},
                                  {1, 2}, {5, 1}};
    const auto st = generate(d, v, edges);
    std::set<std::pair<double, double>> positions;
    std::size_t total = 0;
    for (std::size_t k = 0; k < st.edge_count(); ++k)
        for (const auto& p : st.on_edge(k)) {
            CHECK_FALSE(p.at == MeshVertex{4.0, 3.0});
            CHECK_FALSE(p.at == MeshVertex{2.5, 4.25});
            positions.insert({p.at.col, p.at.row});
            ++total;
        }
    CHECK(positions.size() == total);
    for (std::size_t k = 0; k < st.edge_count(); ++k)
        check_structure(st, k, lattice(g, v[st.edge(k)[0]]), lattice(g, v[st.edge(k)[1]]), kCols,
                        kRows);
}

TEST_CASE("CC-X: extreme scales", "[constraint_points][extreme]") {
    SECTION("an edge across 99,998 columns: 199,999 points, every one in shape") {
        constexpr std::size_t cols = 100001, rows = 3;
        const auto d = dem(cols, rows);
        const std::vector<Point2> v{world(0.5, 0.5), world(99998.5, 1.5)};
        const auto st = generate(d, v, {{0, 1}});
        CHECK(st.size() == 199999);  // 99,998 column crossings and one row crossing
        CHECK(st.duplicates() == 0);
        CHECK(check_structure(st, 0, {0.5, 0.5}, {99998.5, 1.5}, cols, rows) == 99999);
    }
    SECTION("a sub-millimetre edge inside one cell has one midpoint") {
        const auto d = dem(kCols, kRows);
        const Point2 p = world(3.3, 2.7);
        const std::vector<Point2> v{p, {p.x + 1e-4, p.y - 1e-4}};
        const auto st = generate(d, v, {{0, 1}});
        REQUIRE(st.size() == 1);
        check_structure(st, 0, lattice(d.geometry(), v[0]), lattice(d.geometry(), v[1]), kCols,
                        kRows);
    }
    SECTION("a sub-millimetre edge across a column line has its crossing on the line") {
        const auto d = dem(kCols, kRows);
        const Point2 n = d.geometry().node({2, 3});
        const std::vector<Point2> v{{n.x - 5e-5, n.y - 1.25}, {n.x + 5e-5, n.y - 1.2501}};
        const auto st = generate(d, v, {{0, 1}});
        REQUIRE(st.size() == 3);
        CHECK(st.on_edge(0)[1].at.col == 3.0);
        check_structure(st, 0, lattice(d.geometry(), v[0]), lattice(d.geometry(), v[1]), kCols,
                        kRows);
    }
}

// ---------------------------------------------------------------------------
// Refusals (D2: "an end outside the node rectangle, or an index out of range")
// ---------------------------------------------------------------------------

// Both the type and the message: REQUIRE_THROWS_AS alone would pass on an
// invalid_argument thrown by something the generator calls.
#define REQUIRE_REFUSED(expr)                                                                   \
    do {                                                                                        \
        REQUIRE_THROWS_AS(expr, std::invalid_argument);                                         \
        REQUIRE_THROWS_WITH(expr, Catch::Matchers::ContainsSubstring("constraint_check_points")); \
    } while (false)

TEST_CASE("CC-R: refusals are std::invalid_argument naming the generator",
          "[constraint_points][refusal]") {
    const auto d = dem(kCols, kRows);
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    SECTION("an index out of range") {
        REQUIRE_REFUSED(generate(d, kCc1Vertices, {{0, 2}}));
    }
    SECTION("an end outside the node rectangle") {
        const std::vector<Point2> v{world(0.5, 0.5), world(8.5, 2.0)};
        REQUIRE_REFUSED(generate(d, v, {{0, 1}}));
    }
    SECTION("an end a hair beyond the rectangle's top side") {
        const Point2 top = d.geometry().node({0, 2});
        const std::vector<Point2> v{world(0.5, 0.5), {top.x, std::nextafter(top.y, inf)}};
        REQUIRE_REFUSED(generate(d, v, {{0, 1}}));
    }
    SECTION("NaN and infinite coordinates") {
        for (const Point2 bad : {Point2{nan, kY0 - 1.0}, Point2{kX0 + 1.0, nan},
                                 Point2{inf, kY0 - 1.0}, Point2{kX0 + 1.0, -inf}}) {
            const std::vector<Point2> v{world(0.5, 0.5), bad};
            REQUIRE_REFUSED(generate(d, v, {{1, 0}}));
        }
    }
}
