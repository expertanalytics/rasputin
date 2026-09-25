#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>

#include <limits>
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
    [[nodiscard]] std::span<const double> row(std::size_t) const noexcept { return row_; }

private:
    RasterGeometry geometry_;
    double value_;
    std::vector<double> row_;
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

TEST_CASE("contains_strict excludes the boundary", "[raster][geometry]") {
    const auto g = grid_3x5();
    REQUIRE(g.contains_strict(Point2{2.0, 8.0}));
    // The strict test is what lets boundary draping skip edges lying on the
    // raster border instead of re-emitting them as duplicate constraints.
    REQUIRE_FALSE(g.contains_strict(Point2{0.0, 8.0}));
    REQUIRE_FALSE(g.contains_strict(Point2{4.0, 8.0}));
    REQUIRE_FALSE(g.contains_strict(Point2{2.0, 10.0}));
    REQUIRE_FALSE(g.contains_strict(Point2{2.0, 6.0}));
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
        // A non-finite point is not strictly inside anything.
        REQUIRE_FALSE(g.contains_strict(p));
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

TEST_CASE("boundary tolerance stays above rounding noise at UTM scale",
          "[raster][geometry][edge]") {
    // The legacy tolerance scaled only with cell size, which at UTM33
    // northings is about 1.5 ulp -- noise. It must scale with coordinate
    // magnitude too.
    const RasterGeometry g{500000.0, 7900000.0, 10.0, 10.0, 100, 100};

    const double one_ulp = std::nextafter(7900000.0, 1e9) - 7900000.0;
    REQUIRE(g.boundary_epsilon() > one_ulp * 100.0);   // well clear of noise
    REQUIRE(g.boundary_epsilon() < g.delta_x() * 1e-3); // still far below a cell

    // A point one ulp inside the border still reads as on the border.
    REQUIRE_FALSE(g.contains_strict(Point2{std::nextafter(g.x_min(), g.x_max()), 7899500.0}));
    REQUIRE_FALSE(g.contains_strict(Point2{500500.0, std::nextafter(g.y_max(), g.y_min())}));

    // A metre inside is unambiguously interior.
    REQUIRE(g.contains_strict(Point2{500001.0, 7899999.0}));
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
