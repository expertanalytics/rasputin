// Increment 15e, fix 4 (docs/increments/15e-memory-fixes.md): the check-point
// store kept in fixed chunks of kChunk points carved from large slabs, not one
// growing vector per cell row. This is the increment's invariant-critical suite.
//
// Interface this suite reads, beyond CP1's (test_refinement_check_points.cpp):
//
//   CheckPoints::kChunk                  points per chunk, 1024, a public
//                                        static constexpr std::size_t
//   reserved_points()                    chunks handed out x kChunk
//
// What it pins:
//   - the bound: reserved_points() <= N + R * kChunk for N points over R cell
//     rows, before and after freeze(), and freeze() takes no new chunk;
//   - a row longer than two chunks, added in reverse order with duplicates on
//     both sides of a chunk boundary, freezes to CP1's order, keeps CP1's
//     first-of-duplicates, and answers for_each_in column ranges that start and
//     end in different chunks.
//
// Every position is a dyadic fraction of a cell on an integral frame, so the
// stored position equals the given one bit for bit (as in CP1).
//
// This file went red at 9879805 because neither kChunk nor reserved_points()
// existed yet, so the target did not compile; its separate target kept every
// other suite building and running under ctest.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/check_points.hpp>

#include <algorithm>
#include <cstddef>
#include <span>
#include <tuple>
#include <vector>

using terrain::Point2;
using terrain::raster::RasterGeometry;
using terrain::refinement::CheckPoints;

namespace {

constexpr std::size_t kChunk = CheckPoints::kChunk;
static_assert(kChunk == 1024, "15e fixes the chunk at 1024 points (16 KiB)");

// Wide enough for four points per cell over more than three chunks of one row.
constexpr std::size_t kCols = 1025, kRows = 9;  // cells: 8 rows by 1024 columns
constexpr std::size_t kCellRows = kRows - 1;
constexpr std::size_t kPerCell = 4;  // column offsets 0, 1/4, 1/2, 3/4
constexpr double kH = 30.0;

RasterGeometry grid() { return RasterGeometry{1000.0, 2000.0, kH, kH, kCols, kRows}; }

Point2 world(double col, double row) { return Point2{1000.0 + col * kH, 2000.0 - row * kH}; }

struct Seen {
    double col;
    double row;
    float z;
    friend bool operator==(const Seen&, const Seen&) = default;
};

// The k-th point of a cell row in frozen order: cell column k / 4, column
// offset (k % 4) / 4, row offset 1/2.
Seen nth(std::size_t row, std::size_t k, float z) {
    return Seen{static_cast<double>(k / kPerCell) + static_cast<double>(k % kPerCell) / 4.0,
                static_cast<double>(row) + 0.5, z};
}

void add(CheckPoints& cp, const std::vector<Seen>& points) {
    std::vector<Point2> xy;
    std::vector<float> z;
    for (const auto& p : points) {
        xy.push_back(world(p.col, p.row));
        z.push_back(p.z);
    }
    cp.add(std::span<const Point2>{xy}, std::span<const float>{z});
}

std::vector<Seen> in(const CheckPoints& cp, std::size_t row, std::size_t c0, std::size_t c1) {
    std::vector<Seen> out;
    cp.for_each_in(row, c0, c1, [&](const auto& p, auto z) {
        out.push_back(Seen{p.col, p.row, static_cast<float>(z)});
    });
    return out;
}

std::size_t bound(std::size_t n) { return n + kCellRows * kChunk; }

}  // namespace

TEST_CASE("15e: reserved points stay within N + R * kChunk, before and after freeze",
          "[check_points][arena][bound]") {
    // One row of 3 * kChunk + 1 points (four chunks, the last holding one),
    // one point in each of six other rows, and one row left empty.
    constexpr std::size_t kLong = 3 * kChunk + 1;
    std::vector<Seen> long_row;
    for (std::size_t k = 0; k < kLong; ++k)
        long_row.push_back(nth(2, k, static_cast<float>(k % 64)));
    std::vector<Seen> singles;
    for (const std::size_t r : {0u, 1u, 3u, 4u, 5u, 7u})  // row 6 stays empty
        singles.push_back(nth(r, 5, 1.0f));
    const std::size_t n = kLong + singles.size();

    CheckPoints cp{grid()};
    add(cp, singles);
    add(cp, long_row);

    const std::size_t before = cp.reserved_points();
    CAPTURE(n, before, bound(n));
    REQUIRE(before % kChunk == 0);
    REQUIRE(before >= n);  // every point added is held somewhere
    REQUIRE(before <= bound(n));

    cp.freeze();
    const std::size_t after = cp.reserved_points();
    CAPTURE(after);
    REQUIRE(after <= bound(n));
    REQUIRE(after == before);  // gather/scatter back into the same chunks
    REQUIRE(cp.size() == n);
    REQUIRE(cp.duplicates() == 0);

    // The long row reads back whole, in order, across its four chunks.
    REQUIRE(in(cp, 2, 0, kCols - 1) == long_row);
    REQUIRE(in(cp, 6, 0, kCols - 1).empty());
    for (const auto& s : singles)
        REQUIRE(in(cp, static_cast<std::size_t>(s.row), 0, kCols - 1) == std::vector<Seen>{s});
}

TEST_CASE("15e: an empty store reserves nothing beyond the bound", "[check_points][arena][empty]") {
    CheckPoints cp{grid()};
    REQUIRE(cp.reserved_points() <= bound(0));
    cp.freeze();
    REQUIRE(cp.reserved_points() <= bound(0));
    REQUIRE(cp.size() == 0);
}

TEST_CASE("15e: exactly one chunk, and one chunk and one point, stay within the bound",
          "[check_points][arena][bound][edge]") {
    std::vector<Seen> points;
    for (std::size_t k = 0; k < kChunk; ++k)
        points.push_back(nth(0, k, 1.0f));
    for (std::size_t k = 0; k < kChunk + 1; ++k)
        points.push_back(nth(1, k, 1.0f));
    CheckPoints cp{grid()};
    add(cp, points);
    REQUIRE(cp.reserved_points() >= points.size());
    REQUIRE(cp.reserved_points() <= bound(points.size()));
    cp.freeze();
    REQUIRE(cp.reserved_points() <= bound(points.size()));
    REQUIRE(cp.size() == points.size());
    REQUIRE(in(cp, 0, 0, kCols - 1).size() == kChunk);
    REQUIRE(in(cp, 1, 0, kCols - 1).size() == kChunk + 1);
}

TEST_CASE("15e: a row across chunk boundaries, added in reverse with duplicates at the boundaries",
          "[check_points][arena][duplicates][determinism]") {
    constexpr std::size_t kRow = 3;
    constexpr std::size_t kLen = 2 * kChunk + 7;  // three chunks, the last holding seven

    // The frozen row: kLen distinct positions, z = 10 at each.
    std::vector<Seen> expected;
    for (std::size_t k = 0; k < kLen; ++k)
        expected.push_back(nth(kRow, k, 10.0f));

    // Duplicates of the points on both sides of each boundary, so that a
    // duplicate pair straddles index kChunk and index 2 * kChunk in the sorted,
    // not yet unique row. Two have a smaller z: the kept one is the first in
    // (col, row offset, col offset, z) order, so their z replaces 10.
    std::vector<Seen> added = expected;
    const std::vector<std::tuple<std::size_t, float>> dups{
        {kChunk - 1, 20.0f}, {kChunk, 5.0f}, {kChunk, 30.0f},
        {2 * kChunk - 1, 3.0f}, {2 * kChunk, 40.0f},
    };
    for (const auto& [k, z] : dups) {
        added.push_back(nth(kRow, k, z));
        if (z < expected[k].z)
            expected[k].z = z;
    }
    // Reverse of the frozen order: every point arrives after all that follow it.
    std::ranges::sort(added, [](const Seen& a, const Seen& b) {
        return std::tuple{a.col, a.z} > std::tuple{b.col, b.z};
    });

    CheckPoints cp{grid()};
    add(cp, added);
    cp.freeze();

    REQUIRE(cp.size() == kLen);
    REQUIRE(cp.duplicates() == dups.size());
    REQUIRE(cp.reserved_points() <= bound(added.size()));
    REQUIRE(in(cp, kRow, 0, kCols - 1) == expected);
    for (std::size_t r = 0; r < kCellRows; ++r)
        if (r != kRow)
            REQUIRE(in(cp, r, 0, kCols - 1).empty());

    // Column ranges whose ends fall in different chunks, on a boundary cell, at
    // the row's partial last chunk, and past its end. Cell c holds indices
    // 4c..4c+3, so cell 256 starts chunk 1 and cell 512 starts chunk 2.
    const std::vector<std::size_t> ends{0, 1, 200, 255, 256, 257, 300, 511, 512, 513, 514, 600, kCols - 2};
    for (const auto c0 : ends)
        for (const auto c1 : ends) {
            if (c1 < c0)
                continue;
            CAPTURE(c0, c1);
            std::vector<Seen> want;
            for (const auto& s : expected)
                if (s.col >= static_cast<double>(c0) && s.col < static_cast<double>(c1 + 1))
                    want.push_back(s);
            REQUIRE(in(cp, kRow, c0, c1) == want);
        }
}
