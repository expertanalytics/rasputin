// Increment 15c-1 (docs/increments/15c-geographic-dem.md, D4 and "Tests for
// @tester", CP1): the check-point store. Points arrive as world (x, y) in the
// target CRS with a float z, are filed by lattice cell, and after freeze() are
// read back by cell row and column range.
//
// Interface, as D4 names it, with what this suite CHOOSES where D4 is silent:
//
//   CheckPoints(raster::RasterGeometry g)      the target grid's frame, h by h
//   add(std::span<const Point2> xy, std::span<const float> z)
//        sizes must match (std::invalid_argument); after freeze, std::logic_error
//   freeze()
//   size(), duplicates(), outside(), geometry()
//   for_each_in(cell_row, c0, c1, f)           c0..c1 INCLUSIVE cell columns;
//        f(p, z) with p.col, p.row the reconstructed fractional lattice position
//        (the stored position, D4's "the check point is the stored position")
//
// Every position here is a dyadic fraction of a cell on an integral frame, so
// the stored position equals the given one bit for bit and the suite can
// compare exactly. One case checks the 2^-24-cell bound on a non-dyadic point.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/check_points.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <span>
#include <stdexcept>
#include <tuple>
#include <vector>

using terrain::Point2;
using terrain::raster::RasterGeometry;
using terrain::refinement::CheckPoints;

namespace {

constexpr std::size_t kCols = 9, kRows = 7;  // cells: 6 rows by 8 columns
constexpr double kH = 30.0;

RasterGeometry grid() { return RasterGeometry{1000.0, 2000.0, kH, kH, kCols, kRows}; }

Point2 world(double col, double row) { return Point2{1000.0 + col * kH, 2000.0 - row * kH}; }

struct Seen {
    std::size_t cell_row;  // the row queried
    double col;
    double row;
    float z;
    friend bool operator==(const Seen&, const Seen&) = default;
};

// Every stored point, by querying each cell row over every column. With
// inclusive or half-open columns alike, [0, kCols - 1] covers cells 0..kCols-2.
std::vector<Seen> everything(const CheckPoints& cp) {
    std::vector<Seen> out;
    for (std::size_t r = 0; r + 1 < kRows; ++r)
        cp.for_each_in(r, 0, kCols - 1, [&](const auto& p, auto z) {
            out.push_back(Seen{r, p.col, p.row, static_cast<float>(z)});
        });
    return out;
}

std::vector<Seen> cell(const CheckPoints& cp, std::size_t r, std::size_t c) {
    std::vector<Seen> out;
    cp.for_each_in(r, c, c, [&](const auto& p, auto z) {
        out.push_back(Seen{r, p.col, p.row, static_cast<float>(z)});
    });
    return out;
}

void add(CheckPoints& cp, const std::vector<Point2>& xy, const std::vector<float>& z) {
    cp.add(std::span<const Point2>{xy}, std::span<const float>{z});
}

// Seeded dyadic points: offsets k / 1024 of a cell, z in 1/64 steps.
struct Cloud {
    std::vector<Point2> xy;
    std::vector<float> z;
};

Cloud cloud(std::size_t n, std::uint32_t seed) {
    std::mt19937 gen{seed};
    Cloud c;
    for (std::size_t i = 0; i < n; ++i) {
        const double col = static_cast<double>(gen() % ((kCols - 1) * 1024u)) / 1024.0;
        const double row = static_cast<double>(gen() % ((kRows - 1) * 1024u)) / 1024.0;
        c.xy.push_back(world(col, row));
        c.z.push_back(static_cast<float>(gen() % 6400u) / 64.0f);
    }
    return c;
}

}  // namespace

TEST_CASE("CP1: the store keeps its geometry", "[check_points][CP1]") {
    const CheckPoints cp{grid()};
    REQUIRE(cp.geometry().x_min() == 1000.0);
    REQUIRE(cp.geometry().y_max() == 2000.0);
    REQUIRE(cp.geometry().delta_x() == kH);
    REQUIRE(cp.geometry().rows() == kRows);
    REQUIRE(cp.geometry().cols() == kCols);
}

TEST_CASE("CP1: points are filed by cell (floor(row), floor(col))", "[check_points][CP1]") {
    CheckPoints cp{grid()};
    add(cp, {world(2.5, 3.25), world(7.75, 0.5), world(0.125, 5.5), world(3.0, 4.0)},
        {1.0f, 2.0f, 3.0f, 4.0f});
    cp.freeze();
    REQUIRE(cp.size() == 4);
    REQUIRE(cp.outside() == 0);
    REQUIRE(cp.duplicates() == 0);

    REQUIRE(cell(cp, 3, 2) == std::vector<Seen>{{3, 2.5, 3.25, 1.0f}});
    REQUIRE(cell(cp, 0, 7) == std::vector<Seen>{{0, 7.75, 0.5, 2.0f}});
    REQUIRE(cell(cp, 5, 0) == std::vector<Seen>{{5, 0.125, 5.5, 3.0f}});
    // A point on a node is filed in the cell it is the upper-left corner of.
    REQUIRE(cell(cp, 4, 3) == std::vector<Seen>{{4, 3.0, 4.0, 4.0f}});
    REQUIRE(cell(cp, 3, 3).empty());
    REQUIRE(cell(cp, 2, 2).empty());
    REQUIRE(cell(cp, 3, 1).empty());
    // A column range finds exactly the cells in it.
    std::vector<Seen> row3;
    cp.for_each_in(3, 0, 1, [&](const auto&, auto) { row3.push_back({}); });
    REQUIRE(row3.empty());
    cp.for_each_in(3, 2, 7, [&](const auto&, auto) { row3.push_back({}); });
    REQUIRE(row3.size() == 1);
}

TEST_CASE("CP1: a point on the far edges is in the last cell", "[check_points][CP1][edge]") {
    CheckPoints cp{grid()};
    const double last_col = kCols - 1, last_row = kRows - 1;
    add(cp, {world(last_col, last_row), world(last_col, 2.5), world(4.5, last_row), world(0.0, 0.0)},
        {1.0f, 2.0f, 3.0f, 4.0f});
    cp.freeze();
    REQUIRE(cp.size() == 4);
    REQUIRE(cp.outside() == 0);
    REQUIRE(cell(cp, kRows - 2, kCols - 2) == std::vector<Seen>{{kRows - 2, last_col, last_row, 1.0f}});
    REQUIRE(cell(cp, 2, kCols - 2) == std::vector<Seen>{{2, last_col, 2.5, 2.0f}});
    REQUIRE(cell(cp, kRows - 2, 4) == std::vector<Seen>{{kRows - 2, 4.5, last_row, 3.0f}});
    REQUIRE(cell(cp, 0, 0) == std::vector<Seen>{{0, 0.0, 0.0, 4.0f}});
}

TEST_CASE("CP1: outside the node rectangle, or non-finite, is dropped and counted",
          "[check_points][CP1][nonfinite]") {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    constexpr double inf = std::numeric_limits<double>::infinity();
    const RasterGeometry g = grid();
    const double eps_x = g.x_max() * 1e-15;  // well past one ulp outside
    CheckPoints cp{g};
    add(cp,
        {
            Point2{g.x_min() - 1.0, 1950.0},           // west
            Point2{g.x_max() + 1.0, 1950.0},           // east
            Point2{1050.0, g.y_max() + 1.0},           // north
            Point2{1050.0, g.y_min() - 1.0},           // south
            Point2{g.x_max() + 4096.0 * eps_x, 1950.0},  // a hair east
            Point2{nan, 1950.0},
            Point2{1050.0, nan},
            Point2{inf, 1950.0},
            Point2{1050.0, -inf},
            world(1.5, 1.5),  // NaN z
            world(2.5, 2.5),  // infinite z
            world(3.5, 3.5),  // the one kept
        },
        {1, 1, 1, 1, 1, 1, 1, 1, 1, std::numeric_limits<float>::quiet_NaN(),
         std::numeric_limits<float>::infinity(), 7.0f});
    cp.freeze();
    REQUIRE(cp.outside() == 11);
    REQUIRE(cp.size() == 1);
    REQUIRE(cp.duplicates() == 0);
    REQUIRE(everything(cp) == std::vector<Seen>{{3, 3.5, 3.5, 7.0f}});
}

TEST_CASE("CP1: the order after freeze does not depend on the order of add",
          "[check_points][CP1][determinism]") {
    const Cloud c = cloud(400, 20261001);
    CheckPoints one{grid()};
    add(one, c.xy, c.z);
    one.freeze();

    // The same points, permuted, in three add calls of unequal size.
    std::vector<std::size_t> idx(c.xy.size());
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::shuffle(idx.begin(), idx.end(), std::mt19937{7});
    std::vector<Point2> xy;
    std::vector<float> z;
    for (const auto i : idx) {
        xy.push_back(c.xy[i]);
        z.push_back(c.z[i]);
    }
    CheckPoints three{grid()};
    const auto part = [&](std::size_t b, std::size_t e) {
        three.add(std::span<const Point2>{xy}.subspan(b, e - b), std::span<const float>{z}.subspan(b, e - b));
    };
    part(0, 13);
    part(13, 290);
    part(290, xy.size());
    three.freeze();

    REQUIRE(one.size() + one.duplicates() == c.xy.size());
    REQUIRE(three.size() == one.size());
    REQUIRE(three.duplicates() == one.duplicates());
    REQUIRE(everything(three) == everything(one));
}

TEST_CASE("CP1: freeze sorts by (cell row, cell col, row offset, col offset, z)",
          "[check_points][CP1]") {
    const Cloud c = cloud(400, 99);
    CheckPoints cp{grid()};
    add(cp, c.xy, c.z);
    cp.freeze();
    const auto seen = everything(cp);
    REQUIRE(seen.size() == cp.size());
    auto key = [](const Seen& s) {
        const double cr = std::floor(s.row), cc = std::floor(s.col);
        return std::tuple{cr, cc, s.row - cr, s.col - cc, s.z};
    };
    for (std::size_t i = 0; i < seen.size(); ++i) {
        CAPTURE(i);
        REQUIRE(static_cast<double>(seen[i].cell_row) == std::floor(seen[i].row));
        if (i > 0) REQUIRE(key(seen[i - 1]) < key(seen[i]));  // strictly: duplicates are gone
    }
}

TEST_CASE("CP1: two points at one position count one duplicate, whatever the order",
          "[check_points][CP1][duplicates]") {
    const Point2 p = world(4.25, 1.75);
    CheckPoints ab{grid()}, ba{grid()};
    add(ab, {p, world(1.5, 1.5), p}, {5.0f, 1.0f, 9.0f});
    add(ba, {p, world(1.5, 1.5), p}, {9.0f, 1.0f, 5.0f});
    ab.freeze();
    ba.freeze();
    REQUIRE(ab.size() == 2);
    REQUIRE(ab.duplicates() == 1);
    REQUIRE(ba.duplicates() == 1);
    REQUIRE(everything(ab) == everything(ba));  // the same one is kept

    CheckPoints three{grid()};
    add(three, {p, p, p}, {2.0f, 2.0f, 3.0f});
    three.freeze();
    REQUIRE(three.size() == 1);
    REQUIRE(three.duplicates() == 2);
}

TEST_CASE("CP1: add after freeze is refused, and the store is unchanged", "[check_points][CP1]") {
    CheckPoints cp{grid()};
    add(cp, {world(1.5, 1.5)}, {1.0f});
    cp.freeze();
    const std::vector<Point2> more{world(2.5, 2.5)};
    const std::vector<float> mz{2.0f};
    REQUIRE_THROWS_AS(cp.add(std::span<const Point2>{more}, std::span<const float>{mz}), std::logic_error);
    REQUIRE(cp.size() == 1);
    REQUIRE(everything(cp) == std::vector<Seen>{{1, 1.5, 1.5, 1.0f}});
}

TEST_CASE("CP1: add refuses xy and z of different lengths", "[check_points][CP1]") {
    CheckPoints cp{grid()};
    const std::vector<Point2> xy{world(1.5, 1.5), world(2.5, 2.5)};
    const std::vector<float> z{1.0f};
    REQUIRE_THROWS_AS(cp.add(std::span<const Point2>{xy}, std::span<const float>{z}),
                      std::invalid_argument);
}

TEST_CASE("CP1: an empty store freezes and yields nothing", "[check_points][CP1][empty]") {
    CheckPoints cp{grid()};
    cp.freeze();
    REQUIRE(cp.size() == 0);
    REQUIRE(cp.outside() == 0);
    REQUIRE(everything(cp).empty());
}

TEST_CASE("CP1: a non-dyadic point is stored within 2^-23 of a cell", "[check_points][CP1]") {
    CheckPoints cp{grid()};
    const double col = 2.1, row = 4.7;
    add(cp, {world(col, row)}, {1.0f});
    cp.freeze();
    const auto seen = cell(cp, 4, 2);
    REQUIRE(seen.size() == 1);
    REQUIRE(std::abs(seen[0].col - col) <= std::ldexp(1.0, -23));
    REQUIRE(std::abs(seen[0].row - row) <= std::ldexp(1.0, -23));
}
