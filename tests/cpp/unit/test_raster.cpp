#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

using terrain::Point2;
using terrain::raster::CellIndex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::raster::bilinear;

namespace {

// Deliberately non-square: a transposed index is invisible on a square grid.
RasterGeometry grid_3x5() {
    return RasterGeometry{/*x_min=*/0.0, /*y_max=*/10.0, /*delta_x=*/1.0, /*delta_y=*/2.0,
                          /*cols=*/5, /*rows=*/3};
}

// z = a + b*x + c*y sampled on the grid; bilinear must reproduce a plane exactly.
Raster<double> plane_raster(const RasterGeometry& g, double a, double b, double c) {
    std::vector<double> v(g.size());
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c_ = 0; c_ < g.cols(); ++c_) {
            const Point2 p = g.node(CellIndex{r, c_});
            v[g.linear_index(CellIndex{r, c_})] = a + b * p.x + c * p.y;
        }
    return Raster<double>{g, std::move(v)};
}

// An independent RasterSource model: proves the concept is satisfiable by
// something other than the one class it was written around, and that
// sample.hpp binds to it unchanged. Increment 12 (R4) added row() to the
// concept, so the double carries one: every row is the same constant.
class ConstantRaster {
public:
    using value_type = double;
    ConstantRaster(RasterGeometry g, double v) : geometry_{g}, value_{v}, row_(g.cols(), v) {}
    [[nodiscard]] const RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] double value_at(const CellIndex&) const noexcept { return value_; }
    [[nodiscard]] bool is_nodata(const CellIndex&) const noexcept { return false; }
    [[nodiscard]] const std::optional<double>& nodata() const noexcept { return nodata_; }
    [[nodiscard]] std::span<const double> row(std::size_t) const noexcept { return row_; }

private:
    RasterGeometry geometry_;
    double value_;
    std::vector<double> row_;
    std::optional<double> nodata_;  // none; RasterSource requires nodata() (increment 18, C2 (b))
};

} // namespace

TEST_CASE("RasterGeometry rejects degenerate extents", "[raster][geometry][edge]") {
    REQUIRE_THROWS_AS(RasterGeometry(0.0, 0.0, 0.0, 1.0, 2, 2), std::invalid_argument);
    REQUIRE_THROWS_AS(RasterGeometry(0.0, 0.0, 1.0, -1.0, 2, 2), std::invalid_argument);
    REQUIRE_THROWS_AS(RasterGeometry(0.0, 0.0, 1.0, 1.0, 0, 2), std::invalid_argument);
    REQUIRE_THROWS_AS(RasterGeometry(0.0, 0.0, 1.0, 1.0, 2, 0), std::invalid_argument);
}

TEST_CASE("Grid registration: extent spans (n-1) cells, not n", "[raster][geometry]") {
    const auto g = grid_3x5();
    REQUIRE_THAT(g.x_max(), WithinAbs(4.0, 1e-12));   // 0 + (5-1)*1
    REQUIRE_THAT(g.y_min(), WithinAbs(6.0, 1e-12));   // 10 - (3-1)*2
}

TEST_CASE("North-up rows: row 0 is at y_max and y decreases with row", "[raster][geometry]") {
    const auto g = grid_3x5();
    REQUIRE_THAT(g.node(CellIndex{0, 0}).y, WithinAbs(10.0, 1e-12));
    REQUIRE_THAT(g.node(CellIndex{1, 0}).y, WithinAbs(8.0, 1e-12));
    REQUIRE_THAT(g.node(CellIndex{2, 0}).y, WithinAbs(6.0, 1e-12));
}

TEST_CASE("linear_index is row-major and never transposes", "[raster][geometry][regression]") {
    const auto g = grid_3x5();
    // Regression for the legacy defect: coordinates_to_indices emitted {col,row}
    // while extract_buffer_values indexed ptr[idx[0]*cols + idx[1]], i.e.
    // col*cols + row. On a 3x5 grid that both aliases and reads out of range.
    REQUIRE(g.linear_index(CellIndex{0, 0}) == 0);
    REQUIRE(g.linear_index(CellIndex{0, 4}) == 4);
    REQUIRE(g.linear_index(CellIndex{1, 0}) == 5);
    REQUIRE(g.linear_index(CellIndex{2, 4}) == 14);
    REQUIRE(g.linear_index(CellIndex{2, 4}) == g.size() - 1);

    // Every cell maps to a distinct, in-range slot.
    std::vector<bool> seen(g.size(), false);
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c) {
            const std::size_t k = g.linear_index(CellIndex{r, c});
            REQUIRE(k < g.size());
            REQUIRE_FALSE(seen[k]);
            seen[k] = true;
        }
}

TEST_CASE("node and cell_of round-trip over every cell", "[raster][geometry]") {
    const auto g = grid_3x5();
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c) {
            const auto back = g.cell_of(g.node(CellIndex{r, c}));
            REQUIRE(back.has_value());
            REQUIRE(*back == CellIndex{r, c});
        }
}

TEST_CASE("cell_of rejects out-of-domain, clamped_cell_of saturates", "[raster][geometry][edge]") {
    const auto g = grid_3x5();
    REQUIRE_FALSE(g.cell_of(Point2{-1.0, 8.0}).has_value());
    REQUIRE_FALSE(g.cell_of(Point2{99.0, 8.0}).has_value());
    REQUIRE_FALSE(g.cell_of(Point2{2.0, 99.0}).has_value());
    REQUIRE_FALSE(g.cell_of(Point2{2.0, -99.0}).has_value());

    // Legacy clamping behaviour, but callers must now ask for it by name.
    REQUIRE(g.clamped_cell_of(Point2{-1.0, 99.0}) == CellIndex{0, 0});
    REQUIRE(g.clamped_cell_of(Point2{99.0, -99.0}) == CellIndex{2, 4});
}

TEST_CASE("bilinear reproduces a plane exactly", "[raster][sample]") {
    const auto g = grid_3x5();
    const auto r = plane_raster(g, 3.0, 2.0, -5.0);
    for (double x = 0.0; x <= 4.0; x += 0.37)
        for (double y = 6.0; y <= 10.0; y += 0.53) {
            const auto z = bilinear(r, Point2{x, y});
            REQUIRE(z.has_value());
            REQUIRE_THAT(*z, WithinRel(3.0 + 2.0 * x - 5.0 * y, 1e-12));
        }
}

TEST_CASE("bilinear at the far corner does not read out of bounds", "[raster][sample][regression]") {
    // Regression for the legacy defect: get_indices clamped i to rows-1, then
    // the interpolant unconditionally read row i+1 -- past the buffer end.
    const auto g = grid_3x5();
    const auto r = plane_raster(g, 1.0, 1.0, 1.0);

    const auto corner = bilinear(r, Point2{g.x_max(), g.y_min()});
    REQUIRE(corner.has_value());
    REQUIRE_THAT(*corner, WithinRel(1.0 + g.x_max() + g.y_min(), 1e-12));

    REQUIRE(bilinear(r, Point2{g.x_max(), g.y_max()}).has_value());
    REQUIRE(bilinear(r, Point2{g.x_min(), g.y_min()}).has_value());
}

TEST_CASE("bilinear needs a 2x2 neighbourhood", "[raster][sample][edge]") {
    // A single row or column has no cell to interpolate over. The legacy read
    // out of bounds on the very first query here.
    const RasterGeometry one_row{0.0, 0.0, 1.0, 1.0, /*cols=*/4, /*rows=*/1};
    const Raster<double> r{one_row, std::vector<double>(4, 1.0)};
    REQUIRE_FALSE(bilinear(r, Point2{1.5, 0.0}).has_value());

    const RasterGeometry one_col{0.0, 0.0, 1.0, 1.0, /*cols=*/1, /*rows=*/4};
    const Raster<double> c{one_col, std::vector<double>(4, 1.0)};
    REQUIRE_FALSE(bilinear(c, Point2{0.0, -1.5}).has_value());
}

TEST_CASE("bilinear is bounded by its four corner values", "[raster][sample]") {
    const auto g = grid_3x5();
    std::vector<double> v{1, 9, 2, 8, 3, 7, 4, 6, 5, 0, 2, 4, 6, 8, 1};
    const Raster<double> r{g, std::move(v)};
    for (double x = 0.1; x < 4.0; x += 0.29)
        for (double y = 6.1; y < 10.0; y += 0.31) {
            const auto z = bilinear(r, Point2{x, y});
            REQUIRE(z.has_value());
            REQUIRE(*z >= 0.0);
            REQUIRE(*z <= 9.0);
        }
}

TEST_CASE("bilinear propagates NoData as nullopt, never as a silent value",
          "[raster][sample][edge]") {
    const auto g = grid_3x5();
    std::vector<double> v(g.size(), 1.0);
    v[g.linear_index(CellIndex{1, 2})] = -9999.0;
    const Raster<double> r{g, std::move(v), -9999.0};

    REQUIRE_FALSE(bilinear(r, Point2{2.5, 7.5}).has_value());  // touches the hole
    REQUIRE(bilinear(r, Point2{0.5, 9.5}).has_value());        // clear of it
}

TEST_CASE("bilinear survives UTM-scale coordinates with sub-metre spacing",
          "[raster][sample][edge]") {
    // UTM33 Norway with millimetre spacing: x_min + j*delta_x cancels badly if
    // the offset is accumulated naively.
    const RasterGeometry g{500000.0, 7900000.0, 0.001, 0.001, 8, 8};
    const auto r = plane_raster(g, 0.0, 1.0, 1.0);
    const Point2 p{500000.0035, 7899999.9965};
    const auto z = bilinear(r, p);
    REQUIRE(z.has_value());
    REQUIRE_THAT(*z, WithinRel(p.x + p.y, 1e-9));
}

TEST_CASE("non-finite query points are rejected, never answered plausibly",
          "[raster][geometry][sample][edge][regression]") {
    // Regression: the domain guard compared against NaN, and every comparison
    // against NaN is false, so the guard fell through to the saturating path
    // and returned a confident CellIndex{0,0}. bilinear then handed back an
    // engaged optional carrying NaN -- exactly the "plausible-looking edge
    // sample" this API exists to prevent.
    const auto g = grid_3x5();
    const auto r = plane_raster(g, 1.0, 1.0, 1.0);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    for (const Point2 p : {Point2{nan, nan}, Point2{nan, 8.0}, Point2{2.0, nan},
                           Point2{inf, 8.0}, Point2{2.0, inf},
                           Point2{-inf, 8.0}, Point2{2.0, -inf}}) {
        REQUIRE_FALSE(g.cell_of(p).has_value());
        REQUIRE_FALSE(g.bilinear_cell_of(p).has_value());
        REQUIRE_FALSE(bilinear(r, p).has_value());
    }
}

TEST_CASE("NaN in the raster data is treated as NoData", "[raster][sample][edge]") {
    // raster.hpp claims NaN "must never propagate silently into a mesh". That
    // invariant was asserted in a comment and covered by nothing.
    const auto g = grid_3x5();
    std::vector<double> v(g.size(), 1.0);
    v[g.linear_index(CellIndex{1, 2})] = std::numeric_limits<double>::quiet_NaN();
    const Raster<double> r{g, std::move(v)};

    REQUIRE(r.is_nodata(CellIndex{1, 2}));
    REQUIRE_FALSE(bilinear(r, Point2{2.5, 7.5}).has_value());  // touches the NaN
    REQUIRE(bilinear(r, Point2{0.5, 9.5}).has_value());        // clear of it
}

TEST_CASE("Raster rejects data that does not match its geometry", "[raster][edge]") {
    REQUIRE_THROWS_AS((Raster<double>{grid_3x5(), std::vector<double>(14, 0.0)}),
                      std::invalid_argument);
}

TEST_CASE("RasterSource is satisfiable by more than the class it was written around",
          "[raster][concept]") {
    STATIC_REQUIRE(terrain::raster::RasterSource<Raster<double>>);
    STATIC_REQUIRE(terrain::raster::RasterSource<Raster<float>>);
    STATIC_REQUIRE(terrain::raster::RasterSource<ConstantRaster>);

    // A second, independent model must interoperate with sample.hpp unchanged.
    const ConstantRaster c{grid_3x5(), 7.0};
    const auto z = bilinear(c, Point2{2.5, 7.5});
    REQUIRE(z.has_value());
    REQUIRE_THAT(*z, WithinAbs(7.0, 1e-12));
}

TEST_CASE("nodata is fixed at construction and readable", "[raster][edge]") {
    const auto g = grid_3x5();
    const Raster<double> plain{g, std::vector<double>(g.size(), 1.0)};
    REQUIRE_FALSE(plain.nodata().has_value());

    const Raster<double> sentinel{g, std::vector<double>(g.size(), 1.0), -9999.0};
    REQUIRE(sentinel.nodata().has_value());
    REQUIRE_THAT(*sentinel.nodata(), WithinAbs(-9999.0, 1e-12));
}

TEST_CASE("clamped_cell_of agrees with cell_of for interior points", "[raster][geometry]") {
    const auto g = grid_3x5();
    for (double x = 0.0; x <= 4.0; x += 0.5)
        for (double y = 6.0; y <= 10.0; y += 0.5) {
            const auto c = g.cell_of(Point2{x, y});
            REQUIRE(c.has_value());
            REQUIRE(*c == g.clamped_cell_of(Point2{x, y}));
        }
}

// ---------------------------------------------------------------------------
// Increment 27 (docs/increments/27-node-sampling.md, "The rule" and "Tests for
// @tester", S1-S4 and S6): a point that is a DEM node, bit for bit, reads that
// node and nothing else; every other point keeps the four-corner rule.
//
// These cases use only today's API (g.node and bilinear), so the target
// compiles before the change and S1 and S4 fail by assertion. S5, which needs
// RasterGeometry::node_at, is its own target (test_raster_node_at) so its
// compile failure does not hide these.
// ---------------------------------------------------------------------------

namespace {

constexpr double kSentinel27 = -32767.0;

// Non-square, unequal spacing, integer UTM-scale origin: the geometry real
// DEMs have. On it today's bilinear is already exact at a node it accepts, so
// S1's only failure before the change is the refusal it exists for.
RasterGeometry grid_6x7() {
    return RasterGeometry{/*x_min=*/500000.0, /*y_max=*/7900000.0, /*delta_x=*/10.0,
                          /*delta_y=*/5.0, /*cols=*/7, /*rows=*/6};
}

// Cell (r, c) holds 100 r + c + 0.25: distinct, and never the sentinel.
std::vector<double> distinct27(const RasterGeometry& g) {
    std::vector<double> v(g.size());
    for (std::size_t r = 0; r < g.rows(); ++r)
        for (std::size_t c = 0; c < g.cols(); ++c)
            v[g.linear_index(CellIndex{r, c})] = 100.0 * static_cast<double>(r)
                                               + static_cast<double>(c) + 0.25;
    return v;
}

enum class Void { Sentinel, NaN };

// The raster with `hole` set to NoData of the given kind. The sentinel kind
// declares the sentinel; the NaN kind declares none, so only NaN-ness voids.
Raster<double> with_hole(const RasterGeometry& g, CellIndex hole, Void kind) {
    auto v = distinct27(g);
    if (kind == Void::Sentinel) {
        v[g.linear_index(hole)] = kSentinel27;
        return Raster<double>{g, std::move(v), kSentinel27};
    }
    v[g.linear_index(hole)] = std::numeric_limits<double>::quiet_NaN();
    return Raster<double>{g, std::move(v)};
}

// An interior node, a node on the last row, one on the last column, and the
// four grid corners. On the last row and column the clamped bilinear cell is
// the one before, so there the NoData corner sits on the other side.
std::vector<CellIndex> s1_targets(const RasterGeometry& g) {
    const std::size_t lr = g.rows() - 1, lc = g.cols() - 1;
    return {{2, 3}, {lr, 3}, {2, lc}, {0, 0}, {0, lc}, {lr, 0}, {lr, lc}};
}

// The in-range neighbours of n among its eight.
std::vector<CellIndex> neighbours(const RasterGeometry& g, CellIndex n) {
    std::vector<CellIndex> out;
    for (int dr = -1; dr <= 1; ++dr)
        for (int dc = -1; dc <= 1; ++dc) {
            if (dr == 0 && dc == 0)
                continue;
            const auto r = static_cast<long long>(n.row) + dr;
            const auto c = static_cast<long long>(n.col) + dc;
            if (r < 0 || c < 0 || r >= static_cast<long long>(g.rows())
                || c >= static_cast<long long>(g.cols()))
                continue;
            out.push_back(CellIndex{static_cast<std::size_t>(r), static_cast<std::size_t>(c)});
        }
    return out;
}

// Today's four-corner expression, copied verbatim from sample.hpp before
// increment 27 (S3's reference). A reordered expression in sample.hpp changes
// the bits and fails S3; this copy must not be edited to follow it.
template <typename R>
std::optional<double> four_corner_reference(const R& raster, const Point2& p) {
    const RasterGeometry& g = raster.geometry();
    const auto cell = g.bilinear_cell_of(p);
    if (!cell)
        return std::nullopt;
    const CellIndex c00{cell->row, cell->col};
    const CellIndex c01{cell->row, cell->col + 1};
    const CellIndex c10{cell->row + 1, cell->col};
    const CellIndex c11{cell->row + 1, cell->col + 1};
    if (raster.is_nodata(c00) || raster.is_nodata(c01)
        || raster.is_nodata(c10) || raster.is_nodata(c11))
        return std::nullopt;
    const Point2 upper_left = g.node(c00);
    const double tx = (p.x - upper_left.x) / g.delta_x();
    const double ty = (upper_left.y - p.y) / g.delta_y();
    const double z00 = static_cast<double>(raster.value_at(c00));
    const double z01 = static_cast<double>(raster.value_at(c01));
    const double z10 = static_cast<double>(raster.value_at(c10));
    const double z11 = static_cast<double>(raster.value_at(c11));
    return z00 * (1.0 - tx) * (1.0 - ty)
         + z01 * tx * (1.0 - ty)
         + z10 * (1.0 - tx) * ty
         + z11 * tx * ty;
}

bool same_bits(double a, double b) {
    return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b);
}

} // namespace

TEST_CASE("S1: a node next to NoData is valid and reads the node value exactly",
          "[raster][sample][nodata][increment27]") {
    const auto g = grid_6x7();
    for (const Void kind : {Void::Sentinel, Void::NaN})
        for (const CellIndex n : s1_targets(g))
            for (const CellIndex hole : neighbours(g, n)) {
                CAPTURE(kind == Void::NaN, n.row, n.col, hole.row, hole.col);
                const auto r = with_hole(g, hole, kind);
                const auto z = bilinear(r, g.node(n));
                REQUIRE(z.has_value());
                REQUIRE(same_bits(*z, r.value_at(n)));
            }
}

TEST_CASE("S1: every neighbour NoData at once still leaves the node its own value",
          "[raster][sample][nodata][increment27]") {
    // The strongest form: all eight neighbours void, the node alone valid.
    const auto g = grid_6x7();
    const CellIndex n{2, 3};
    auto v = distinct27(g);
    for (const CellIndex h : neighbours(g, n))
        v[g.linear_index(h)] = kSentinel27;
    const Raster<double> r{g, std::move(v), kSentinel27};
    const auto z = bilinear(r, g.node(n));
    REQUIRE(z.has_value());
    REQUIRE(same_bits(*z, r.value_at(n)));
}

TEST_CASE("S2: a NoData node is refused, sentinel and NaN",
          "[raster][sample][nodata][increment27]") {
    const auto g = grid_6x7();
    for (const Void kind : {Void::Sentinel, Void::NaN})
        for (const CellIndex n : s1_targets(g)) {
            CAPTURE(kind == Void::NaN, n.row, n.col);
            const auto r = with_hole(g, n, kind);
            REQUIRE_FALSE(bilinear(r, g.node(n)).has_value());
        }
}

TEST_CASE("S3: an off-node point with a NoData corner is refused",
          "[raster][sample][nodata][increment27]") {
    const auto g = grid_6x7();
    const CellIndex n{2, 3};
    const Point2 at = g.node(n);
    for (const Void kind : {Void::Sentinel, Void::NaN}) {
        CAPTURE(kind == Void::NaN);

        SECTION("a cell-interior point, its far corner void") {
            const auto r = with_hole(g, CellIndex{3, 4}, kind);
            REQUIRE_FALSE(bilinear(r, Point2{at.x + 3.0, at.y - 1.25}).has_value());
        }
        SECTION("on a cell's top side, NoData across the side (cell sides are deferred)") {
            // y on the node line of row 2, x mid-cell: the row-3 corners weigh zero.
            const auto r = with_hole(g, CellIndex{3, 3}, kind);
            REQUIRE_FALSE(bilinear(r, Point2{at.x + 5.0, at.y}).has_value());
        }
        SECTION("on a cell's left side, NoData across the side (cell sides are deferred)") {
            // x on the node line of column 3, y mid-cell: the column-4 corners weigh zero.
            const auto r = with_hole(g, CellIndex{2, 4}, kind);
            REQUIRE_FALSE(bilinear(r, Point2{at.x, at.y - 2.5}).has_value());
        }
        SECTION("one ulp from a node toward its NoData neighbour, in x (exactness)") {
            const auto r = with_hole(g, CellIndex{2, 4}, kind);
            const Point2 p{std::nextafter(at.x, std::numeric_limits<double>::infinity()), at.y};
            REQUIRE(p.x != at.x);
            REQUIRE(g.bilinear_cell_of(p) == CellIndex{2, 3});  // the probe is in the cell meant
            REQUIRE_FALSE(bilinear(r, p).has_value());
        }
        SECTION("one ulp from a node toward its NoData neighbour, in y (exactness)") {
            const auto r = with_hole(g, CellIndex{3, 3}, kind);
            const Point2 p{at.x, std::nextafter(at.y, -std::numeric_limits<double>::infinity())};
            REQUIRE(p.y != at.y);
            REQUIRE(g.bilinear_cell_of(p) == CellIndex{2, 3});
            REQUIRE_FALSE(bilinear(r, p).has_value());
        }
    }
}

TEST_CASE("S3: off-node points on a NoData-free raster keep today's bits",
          "[raster][sample][increment27][invariant]") {
    // Non-dyadic geometry and values drawn from a seeded generator. Points:
    // random interior points, cell-side points (one coordinate on a node
    // line), and nextafter neighbours of nodes in all four directions -- none
    // of which is a node, so each must take the unchanged four-corner path.
    const RasterGeometry g{500000.3, 7900000.7, 0.7, 0.3, /*cols=*/23, /*rows=*/17};
    std::mt19937 rng{27};
    std::uniform_real_distribution<double> value{0.0, 1000.0};
    std::vector<double> v(g.size());
    for (double& x : v)
        x = value(rng);
    const Raster<double> r{g, std::move(v)};

    std::uniform_real_distribution<double> ux{g.x_min(), g.x_max()};
    std::uniform_real_distribution<double> uy{g.y_min(), g.y_max()};
    std::vector<Point2> pts;
    for (int i = 0; i < 2000; ++i)
        pts.push_back(Point2{ux(rng), uy(rng)});
    for (std::size_t row = 0; row < g.rows(); ++row)
        for (std::size_t col = 0; col < g.cols(); ++col) {
            const Point2 at = g.node(CellIndex{row, col});
            pts.push_back(Point2{at.x, uy(rng)});  // on a vertical node line
            pts.push_back(Point2{ux(rng), at.y});  // on a horizontal node line
            constexpr double inf = std::numeric_limits<double>::infinity();
            pts.push_back(Point2{std::nextafter(at.x, inf), at.y});
            pts.push_back(Point2{std::nextafter(at.x, -inf), at.y});
            pts.push_back(Point2{at.x, std::nextafter(at.y, inf)});
            pts.push_back(Point2{at.x, std::nextafter(at.y, -inf)});
        }

    std::size_t compared = 0;
    for (const Point2& p : pts) {
        CAPTURE(p.x, p.y);
        const auto expected = four_corner_reference(r, p);
        const auto z = bilinear(r, p);
        if (!expected) {  // outside the rectangle (a nextafter past the far edge)
            REQUIRE_FALSE(z.has_value());
            continue;
        }
        REQUIRE(z.has_value());
        REQUIRE(same_bits(*z, *expected));
        ++compared;
    }
    // Every random point and every side point is inside; only far-edge
    // nextafter probes may fall out, so most of the set was compared.
    REQUIRE(compared > pts.size() * 9 / 10);
}

TEST_CASE("S4: exact at every node of a non-dyadic grid, far edges included",
          "[raster][sample][increment27][invariant]") {
    // The measurement in "Why this test": on this geometry today's bilinear
    // returns a value a few ulps off the node value at most nodes (cell_of
    // can put a node into the cell before it). The green assertion holds on
    // any compiler: the points and node_at's comparison both come from
    // RasterGeometry::node.
    const RasterGeometry g{500000.3, 7900000.7, 0.7, 0.3, /*cols=*/200, /*rows=*/300};
    std::mt19937 rng{1};
    std::uniform_real_distribution<double> value{0.0, 1000.0};
    std::vector<double> v(g.size());
    for (double& x : v)
        x = value(rng);
    const Raster<double> r{g, std::move(v)};

    std::size_t refused = 0, inexact = 0;
    std::optional<CellIndex> first;
    for (std::size_t row = 0; row < g.rows(); ++row)
        for (std::size_t col = 0; col < g.cols(); ++col) {
            const CellIndex c{row, col};
            const auto z = bilinear(r, g.node(c));
            if (!z)
                ++refused;
            else if (!same_bits(*z, r.value_at(c)))
                ++inexact;
            else
                continue;
            if (!first)
                first = c;
        }
    const CellIndex first_bad = first.value_or(CellIndex{g.rows(), g.cols()});
    CAPTURE(refused, inexact, g.size(), first_bad.row, first_bad.col);
    REQUIRE(refused == 0);
    REQUIRE(inexact == 0);
}

TEST_CASE("S6: a raster too small for a cell still refuses its nodes",
          "[raster][sample][edge][increment27]") {
    // bilinear_cell_of's guard stays first, so a 1 x N or N x 1 raster has no
    // value anywhere, nodes included.
    const RasterGeometry one_row{500000.0, 7900000.0, 10.0, 5.0, /*cols=*/4, /*rows=*/1};
    const RasterGeometry one_col{500000.0, 7900000.0, 10.0, 5.0, /*cols=*/1, /*rows=*/4};
    const RasterGeometry one_node{500000.0, 7900000.0, 10.0, 5.0, /*cols=*/1, /*rows=*/1};
    for (const RasterGeometry& g : {one_row, one_col, one_node}) {
        const Raster<double> r{g, std::vector<double>(g.size(), 3.5)};
        for (std::size_t row = 0; row < g.rows(); ++row)
            for (std::size_t col = 0; col < g.cols(); ++col) {
                CAPTURE(g.rows(), g.cols(), row, col);
                REQUIRE_FALSE(bilinear(r, g.node(CellIndex{row, col})).has_value());
            }
    }
}

TEST_CASE("S6: outside and non-finite points still give nullopt",
          "[raster][sample][edge][increment27]") {
    const auto g = grid_6x7();
    const Raster<double> r{g, distinct27(g)};
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    const Point2 nw = g.node(CellIndex{0, 0});
    const Point2 se = g.node(CellIndex{g.rows() - 1, g.cols() - 1});
    for (const Point2 p : {
             Point2{std::nextafter(nw.x, -inf), nw.y},  // one ulp west of a corner node
             Point2{nw.x, std::nextafter(nw.y, inf)},   // one ulp north
             Point2{std::nextafter(se.x, inf), se.y},   // one ulp east
             Point2{se.x, std::nextafter(se.y, -inf)},  // one ulp south
             Point2{nw.x - 10.0, nw.y},                 // a whole cell out, on a node line
             Point2{nan, nw.y}, Point2{nw.x, nan}, Point2{inf, nw.y}, Point2{nw.x, -inf}}) {
        CAPTURE(p.x, p.y);
        REQUIRE_FALSE(bilinear(r, p).has_value());
    }
}
