// Increment 26b (docs/increments/26-land-cover.md, 1.1 "in which coordinates",
// 1.2 "The exact overlap", section 8 "Degeneracy policy" and section 10's
// "26b"): the exact overlap of triangles, given in the class grid's index
// space, with the grid's cells. Not named invariant-critical (section 9: 26b
// has no mutation round).
//
// Interface, as 1.2 fixes it, with the pins the handback lists:
//
//   #include <terrain/land_cover/overlap.hpp>
//   namespace terrain::land_cover {
//   struct ClassAreas {
//       std::vector<std::uint32_t> offsets;  // T + 1
//       std::vector<std::uint8_t>  classes;  // ascending within a row
//       std::vector<double>        areas;    // m^2, > 0
//   };
//   struct OutsideWindow : std::invalid_argument { ... };
//   ClassAreas class_areas(std::span<const Point2> corners,     // 3 per triangle
//                          const raster::RasterView<std::uint8_t>& classes,
//                          std::span<const double> row_area,    // one per row
//                          unsigned threads = 0);               // 0: all cores
//   }
//
// Index space (1.1): x is the column and y the row, so y grows downward and
// cell (r, c) is the unit square [c, c + 1] x [r, r + 1]; the window is
// [0, cols] x [0, rows], closed. The view's geometry is
// RasterGeometry{0, 0, 1, 1, cols, rows}, of which only rows(), cols() and
// row(r) are read; its nodata is not consulted (class 0 is a class, 1.1). A
// triangle may be given in either orientation and its areas are positive.
// OutsideWindow derives from std::invalid_argument; so does the refusal of a
// zero-area triangle (section 8: "none in a valid mesh; refused if met"), of a
// non-finite corner, of a corners span whose length is not a multiple of
// three, and of a row_area span whose length is not rows().
//
// ORACLES: every expected piece is a trapezoid or triangle worked out by hand
// in the test (the arithmetic is in the comments), the shoelace area of the
// whole triangle, and bitwise comparison between runs. The bound is 1e-12 in
// cell units (coordinates below 32 here; a piece's shoelace sum is taken in
// coordinates local to its cell, so its rounding is near 1e-16).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <terrain/core/point.hpp>
#include <terrain/land_cover/overlap.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/view.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using terrain::Point2;
using terrain::land_cover::class_areas;
using terrain::land_cover::ClassAreas;
using terrain::land_cover::OutsideWindow;
using terrain::raster::RasterGeometry;
using terrain::raster::RasterView;

namespace {

constexpr double kCell = 1e-12;  // in cell units; coordinates below 32 in this file

// A class grid in index space with one cell area per row.
struct Grid {
    std::size_t rows;
    std::size_t cols;
    std::vector<std::uint8_t> cells;
    std::vector<double> row_area;

    [[nodiscard]] RasterView<std::uint8_t> view() const {
        return RasterView<std::uint8_t>{RasterGeometry{0.0, 0.0, 1.0, 1.0, cols, rows},
                                        cells.data(), std::nullopt};
    }
    [[nodiscard]] std::uint8_t at(std::size_t r, std::size_t c) const { return cells[r * cols + c]; }
};

std::vector<double> ones(std::size_t rows) { return std::vector<double>(rows, 1.0); }

// Cell (r, c) holds 10 r + c + 1: every class distinct and none 0, so a piece
// in the wrong cell shows as the wrong class.
Grid distinct_grid(std::size_t rows, std::size_t cols, std::vector<double> row_area = {}) {
    Grid g{rows, cols, std::vector<std::uint8_t>(rows * cols), {}};
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c)
            g.cells[r * cols + c] = static_cast<std::uint8_t>(10 * r + c + 1);
    g.row_area = row_area.empty() ? ones(rows) : std::move(row_area);
    return g;
}

Grid uniform_grid(std::size_t rows, std::size_t cols, std::uint8_t cls,
                  std::vector<double> row_area = {}) {
    Grid g{rows, cols, std::vector<std::uint8_t>(rows * cols, cls), {}};
    g.row_area = row_area.empty() ? ones(rows) : std::move(row_area);
    return g;
}

// Sixteen classes including 0, mixed across the grid.
Grid mixed_grid(std::size_t rows, std::size_t cols) {
    Grid g{rows, cols, std::vector<std::uint8_t>(rows * cols), ones(rows)};
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c)
            g.cells[r * cols + c] = static_cast<std::uint8_t>((7 * r + 3 * c) % 16);
    return g;
}

ClassAreas areas_of(const Grid& g, const std::vector<Point2>& corners, unsigned threads = 1) {
    return class_areas(std::span<const Point2>{corners}, g.view(),
                       std::span<const double>{g.row_area}, threads);
}

using Row = std::map<std::uint8_t, double>;

Row row_of(const ClassAreas& out, std::size_t t) {
    Row row;
    for (std::uint32_t i = out.offsets.at(t); i < out.offsets.at(t + 1); ++i)
        row[out.classes.at(i)] = out.areas.at(i);
    return row;
}

double sum_of(const Row& row) {
    double s = 0.0;
    for (const auto& [cls, area] : row)
        s += area;
    return s;
}

double shoelace(const Point2& a, const Point2& b, const Point2& c) {
    return 0.5 * std::fabs(terrain::cross(b - a, c - a));
}

// The row structure of 1.2: T + 1 offsets from 0, classes strictly ascending
// within each row, every area positive.
void require_well_formed(const ClassAreas& out, std::size_t triangles) {
    REQUIRE(out.offsets.size() == triangles + 1);
    REQUIRE(out.offsets.front() == 0);
    REQUIRE(out.offsets.back() == out.classes.size());
    REQUIRE(out.classes.size() == out.areas.size());
    for (std::size_t t = 0; t < triangles; ++t) {
        REQUIRE(out.offsets[t] <= out.offsets[t + 1]);
        for (std::uint32_t i = out.offsets[t]; i < out.offsets[t + 1]; ++i) {
            REQUIRE(out.areas[i] > 0.0);
            if (i + 1 < out.offsets[t + 1])
                REQUIRE(out.classes[i] < out.classes[i + 1]);
        }
    }
}

void require_same_row(const Row& got, const Row& want, double tolerance = kCell) {
    REQUIRE(got.size() == want.size());
    for (const auto& [cls, area] : want) {
        REQUIRE(got.contains(cls));
        REQUIRE_THAT(got.at(cls), WithinAbs(area, tolerance));
    }
}

// A fixed linear congruential generator: the same triangles on every run and
// every platform, which is what a bitwise comparison between runs needs.
struct Lcg {
    std::uint64_t state;
    double next() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(state >> 11) * 0x1.0p-53;
    }
};

// n triangles inside [0, extent]^2 with positive area; every fifth has one
// vertex snapped to a grid corner and every seventh one on a grid line, so the
// degenerate positions of section 8 occur among them.
std::vector<Point2> random_triangles(std::size_t n, double extent, std::uint64_t seed) {
    Lcg rng{seed};
    std::vector<Point2> out;
    out.reserve(3 * n);
    while (out.size() < 3 * n) {
        Point2 p[3];
        for (auto& v : p)
            v = Point2{rng.next() * extent, rng.next() * extent};
        const std::size_t k = out.size() / 3;
        if (k % 5 == 0)
            p[0] = Point2{std::floor(p[0].x), std::floor(p[0].y)};
        if (k % 7 == 0)
            p[1].y = std::floor(p[1].y);
        if (shoelace(p[0], p[1], p[2]) < 1e-6)
            continue;
        out.insert(out.end(), p, p + 3);
    }
    return out;
}

bool same_bits(const ClassAreas& a, const ClassAreas& b) {
    return a.offsets == b.offsets && a.classes == b.classes &&
           a.areas.size() == b.areas.size() &&
           std::memcmp(a.areas.data(), b.areas.data(), a.areas.size() * sizeof(double)) == 0;
}

} // namespace

// ------------------------------------------------------------ hand cases (section 10)

TEST_CASE("1. a triangle wholly inside one cell: one entry, that class, its area times the row's",
          "[land_cover][overlap]") {
    const auto g = distinct_grid(3, 4, {887.5, 836.0, 800.0});
    const std::vector<Point2> t{{2.2, 1.3}, {2.9, 1.4}, {2.5, 1.9}};  // in cell (1, 2)
    const auto out = areas_of(g, t);
    require_well_formed(out, 1);
    REQUIRE(out.offsets == std::vector<std::uint32_t>{0, 1});
    REQUIRE(out.classes == std::vector<std::uint8_t>{g.at(1, 2)});
    REQUIRE_THAT(out.areas[0], WithinRel(shoelace(t[0], t[1], t[2]) * 836.0, 1e-12));
}

TEST_CASE("2. a right triangle with legs 2 on a 2 x 2 grid: 1, 0.5, 0.5, the fourth cell absent",
          "[land_cover][overlap]") {
    const auto g = distinct_grid(2, 2);
    // x + y <= 2: cell (0, 0) whole, (0, 1) and (1, 0) half, (1, 1) the one
    // point (1, 1), area zero, so no entry.
    const std::vector<Point2> t{{0.0, 0.0}, {2.0, 0.0}, {0.0, 2.0}};
    const auto out = areas_of(g, t);
    require_well_formed(out, 1);
    require_same_row(row_of(out, 0), {{g.at(0, 0), 1.0}, {g.at(0, 1), 0.5}, {g.at(1, 0), 0.5}});

    SECTION("one class everywhere: the pieces merge into one entry of area 2") {
        const auto one = areas_of(uniform_grid(2, 2, 33), t);
        require_well_formed(one, 1);
        require_same_row(row_of(one, 0), {{33, 2.0}});
    }
}

TEST_CASE("3. vertices on grid corners: no zero-area entry from the cells they touch",
          "[land_cover][overlap][degenerate]") {
    const auto g = distinct_grid(4, 4);
    // (1, 1) and (3, 1) are corners; (2, 3) too. Row 1 (y in [1, 2], s = y - 1):
    // left edge x = 1 + s/2, right x = 3 - s/2; cell (1, 1) holds width 1 - s/2,
    // area 3/4, cell (1, 2) the same. Row 2 (s = y - 2): x from 1.5 + s/2 to
    // 2.5 - s/2; cells (2, 1) and (2, 2) each width 0.5 - s/2, area 1/4. The
    // corner (1, 1) touches cells (0, 0), (0, 1) and (1, 0) in one point.
    const std::vector<Point2> t{{1.0, 1.0}, {3.0, 1.0}, {2.0, 3.0}};
    const auto out = areas_of(g, t);
    require_well_formed(out, 1);
    require_same_row(row_of(out, 0), {{g.at(1, 1), 0.75}, {g.at(1, 2), 0.75},
                                      {g.at(2, 1), 0.25}, {g.at(2, 2), 0.25}});
    REQUIRE_THAT(sum_of(row_of(out, 0)), WithinAbs(2.0, kCell));
}

TEST_CASE("4. an edge lying on a grid line: the row below it gets nothing",
          "[land_cover][overlap][degenerate]") {
    const auto g = distinct_grid(3, 3);
    // Base on y = 1 from x = 0.5 to 2.5, apex (1.5, 2.5). Row 1 (s = y - 1):
    // x from 0.5 + s/1.5 to 2.5 - s/1.5; cell (1, 0) holds width 0.5 - s/1.5
    // while s < 0.75, area 3/16; cell (1, 2) the same; cell (1, 1) the rest of
    // the trapezoid (2 + 2/3)/2 = 4/3, so 4/3 - 3/8 = 23/24. Row 2 (s = y - 2,
    // up to 0.5): width 2/3 - 2s/1.5, all in cell (2, 1), area 1/6. Total 3/2.
    const std::vector<Point2> t{{0.5, 1.0}, {2.5, 1.0}, {1.5, 2.5}};
    const auto out = areas_of(g, t);
    require_well_formed(out, 1);
    require_same_row(row_of(out, 0), {{g.at(1, 0), 3.0 / 16}, {g.at(1, 1), 23.0 / 24},
                                      {g.at(1, 2), 3.0 / 16}, {g.at(2, 1), 1.0 / 6}});
    REQUIRE_THAT(sum_of(row_of(out, 0)), WithinAbs(1.5, kCell));
}

TEST_CASE("5. a long sliver along a row: ten positive pieces that sum to the triangle's area",
          "[land_cover][overlap][degenerate]") {
    const auto g = distinct_grid(1, 10);
    const std::vector<Point2> t{{0.25, 0.5}, {9.75, 0.5}, {5.0, 0.5 + 1e-6}};
    const auto out = areas_of(g, t);
    require_well_formed(out, 1);
    const auto row = row_of(out, 0);
    REQUIRE(row.size() == 10);
    // The area is 4.75e-6; the pieces are built from crossings with x = 1 .. 9,
    // each placed to about 1e-16, so 1e-15 on the sum is loose by a decade.
    REQUIRE_THAT(sum_of(row), WithinAbs(shoelace(t[0], t[1], t[2]), 1e-15));
}

TEST_CASE("6. a triangle a nanometre across keeps its positive area",
          "[land_cover][overlap][degenerate][scale]") {
    const auto g = uniform_grid(2, 2, 3);
    const std::vector<Point2> t{{0.5, 0.5}, {0.5 + 1e-9, 0.5}, {0.5, 0.5 + 1e-9}};
    const auto out = areas_of(g, t);
    require_well_formed(out, 1);
    // The intended area is 0.5e-18, but 0.5 + 1e-9 rounds to a double whose
    // legs are 1e-9 - 2.8e-17 (relative 2.8e-8), so the triangle as rounded
    // has area 4.9999997e-19, relative 5.7e-8 below the intended one. Compare
    // against the shoelace area of the corners as rounded, which the clip must
    // reproduce to the rounding of its own sums (relative 1e-9 is generous).
    REQUIRE_THAT(row_of(out, 0).at(3), WithinRel(shoelace(t[0], t[1], t[2]), 1e-9));
    REQUIRE(row_of(out, 0).at(3) > 0.4e-18);
}

TEST_CASE("7. rows weigh differently: each piece is multiplied by its row's cell area",
          "[land_cover][overlap]") {
    const auto g = uniform_grid(3, 3, 7, {887.5, 836.0, 800.0});
    // x + y <= 3 above (0.5, 0.5): row 0 holds y in [0.5, 1], width 2 to 1.5,
    // area 7/8; row 1 width 1.5 to 0.5, area 1; row 2 y in [2, 2.5], width 0.5
    // to 0, area 1/8.
    const std::vector<Point2> t{{0.5, 0.5}, {2.5, 0.5}, {0.5, 2.5}};
    const auto out = areas_of(g, t);
    require_well_formed(out, 1);
    require_same_row(row_of(out, 0), {{7, 0.875 * 887.5 + 1.0 * 836.0 + 0.125 * 800.0}}, 1e-9);
}

// ------------------------------------------------------------ order and neighbours

TEST_CASE("8. the vertex order does not change the areas: a rotation, and the mirror orientation",
          "[land_cover][overlap][determinism]") {
    const auto g = mixed_grid(8, 8);
    const std::vector<Point2> ccw{{0.7, 0.3}, {6.2, 1.9}, {2.4, 5.6}};
    const std::vector<Point2> rotated{ccw[1], ccw[2], ccw[0]};
    const std::vector<Point2> mirrored{ccw[0], ccw[2], ccw[1]};
    const auto a = areas_of(g, ccw);
    require_well_formed(a, 1);
    REQUIRE_THAT(sum_of(row_of(a, 0)), WithinAbs(shoelace(ccw[0], ccw[1], ccw[2]), kCell));
    require_same_row(row_of(areas_of(g, rotated), 0), row_of(a, 0));
    require_same_row(row_of(areas_of(g, mirrored), 0), row_of(a, 0));
}

TEST_CASE("9. two triangles on a shared edge: in every cell their pieces sum to the cell's covered area",
          "[land_cover][overlap][determinism]") {
    // The square [0.3, 2.7]^2 cut by a diagonal. Covered per cell: the corners
    // 0.7^2, the edges 0.7, the middle 1. Every cell has its own class.
    const auto g = distinct_grid(3, 3);
    const Point2 a{0.3, 0.3}, b{2.7, 0.3}, c{2.7, 2.7}, d{0.3, 2.7};
    for (const auto& diagonal : {std::vector<Point2>{a, b, c, a, c, d}, std::vector<Point2>{a, b, d, b, c, d}}) {
        const auto out = areas_of(g, diagonal);
        require_well_formed(out, 2);
        Row total;
        for (std::size_t t = 0; t < 2; ++t)
            for (const auto& [cls, area] : row_of(out, t))
                total[cls] += area;
        Row want;
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t col = 0; col < 3; ++col) {
                const double w = col == 1 ? 1.0 : 0.7;
                const double h = r == 1 ? 1.0 : 0.7;
                want[g.at(r, col)] = w * h;
            }
        require_same_row(total, want);
    }
}

TEST_CASE("10. random triangles: every row well formed and summing to the shoelace area",
          "[land_cover][overlap][property]") {
    const auto g = mixed_grid(32, 32);
    const auto corners = random_triangles(2000, 32.0, 26);
    const auto out = areas_of(g, corners);
    require_well_formed(out, 2000);
    for (std::size_t t = 0; t < 2000; ++t) {
        const auto row = row_of(out, t);
        REQUIRE(row.size() <= 16);
        REQUIRE_THAT(sum_of(row), WithinAbs(shoelace(corners[3 * t], corners[3 * t + 1], corners[3 * t + 2]), kCell));
    }
}

// ------------------------------------------------------------ parallel and chunks

TEST_CASE("11. four threads, and all cores, give one thread's bits",
          "[land_cover][overlap][determinism][threads]") {
    const auto g = mixed_grid(32, 32);
    const auto corners = random_triangles(5000, 32.0, 2026);
    const auto serial = areas_of(g, corners, 1);
    require_well_formed(serial, 5000);
    REQUIRE(same_bits(serial, areas_of(g, corners, 4)));
    REQUIRE(same_bits(serial, areas_of(g, corners, 0)));
}

TEST_CASE("12. chunks of triangles concatenated equal one call, bit for bit",
          "[land_cover][overlap][determinism]") {
    const auto g = mixed_grid(32, 32);
    const auto corners = random_triangles(1000, 32.0, 7);
    const auto whole = areas_of(g, corners, 2);
    ClassAreas joined{{0}, {}, {}};
    for (std::size_t begin = 0; begin < corners.size(); begin += 3 * 7) {
        const std::size_t end = std::min(corners.size(), begin + 3 * 7);
        const std::vector<Point2> chunk(corners.begin() + static_cast<std::ptrdiff_t>(begin),
                                        corners.begin() + static_cast<std::ptrdiff_t>(end));
        const auto part = areas_of(g, chunk, 2);
        const std::uint32_t base = joined.offsets.back();
        for (std::size_t t = 1; t < part.offsets.size(); ++t)
            joined.offsets.push_back(base + part.offsets[t]);
        joined.classes.insert(joined.classes.end(), part.classes.begin(), part.classes.end());
        joined.areas.insert(joined.areas.end(), part.areas.begin(), part.areas.end());
    }
    REQUIRE(same_bits(whole, joined));
}

TEST_CASE("13. no triangles: one offset and nothing else", "[land_cover][overlap][empty]") {
    const auto out = areas_of(distinct_grid(2, 2), {});
    REQUIRE(out.offsets == std::vector<std::uint32_t>{0});
    REQUIRE(out.classes.empty());
    REQUIRE(out.areas.empty());
}

// ------------------------------------------------------------ the window and refusals

TEST_CASE("14. the window is [0, cols] x [0, rows], closed: on the edge is inside, past it is OutsideWindow",
          "[land_cover][overlap][refusal]") {
    const auto g = distinct_grid(4, 4);
    const auto right = areas_of(g, {{3.0, 0.5}, {4.0, 0.5}, {3.5, 1.5}});
    require_well_formed(right, 1);
    REQUIRE_THAT(sum_of(row_of(right, 0)), WithinAbs(0.5, kCell));
    const auto bottom = areas_of(g, {{0.5, 3.0}, {1.5, 4.0}, {0.5, 4.0}});
    require_well_formed(bottom, 1);
    REQUIRE_THAT(sum_of(row_of(bottom, 0)), WithinAbs(0.5, kCell));

    REQUIRE_THROWS_AS(areas_of(g, {{3.5, 0.5}, {4.5, 0.5}, {4.0, 1.5}}), OutsideWindow);
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, -0.1}, {1.5, 0.5}, {0.5, 0.9}}), OutsideWindow);
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, 3.5}, {1.5, 3.5}, {1.0, 4.1}}), OutsideWindow);
    // A whole batch is refused, not clamped: the refusal is an exception of
    // the std::invalid_argument family, as the raster headers throw.
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, 0.5}, {1.5, 0.5}, {1.0, 1.5}, {-1.0, 0.5}, {1.5, 0.5}, {1.0, 1.5}}),
                      std::invalid_argument);
}

TEST_CASE("15. a zero-area triangle is refused: collinear, and a repeated vertex",
          "[land_cover][overlap][refusal][degenerate]") {
    const auto g = distinct_grid(4, 4);
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, 0.5}, {1.5, 1.5}, {2.5, 2.5}}), std::invalid_argument);
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, 0.5}, {0.5, 0.5}, {1.5, 1.5}}), std::invalid_argument);
}

TEST_CASE("16. a non-finite corner is refused", "[land_cover][overlap][refusal][nonfinite]") {
    const auto g = distinct_grid(4, 4);
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, nan}, {1.5, 0.5}, {1.0, 1.5}}), std::invalid_argument);
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, 0.5}, {inf, 0.5}, {1.0, 1.5}}), std::invalid_argument);
}

TEST_CASE("17. spans of the wrong length are refused before any work",
          "[land_cover][overlap][refusal]") {
    const auto g = distinct_grid(4, 4);
    REQUIRE_THROWS_AS(areas_of(g, {{0.5, 0.5}, {1.5, 0.5}}), std::invalid_argument);
    Grid short_rows = g;
    short_rows.row_area.pop_back();
    REQUIRE_THROWS_AS(areas_of(short_rows, {{0.5, 0.5}, {1.5, 0.5}, {1.0, 1.5}}), std::invalid_argument);
}
