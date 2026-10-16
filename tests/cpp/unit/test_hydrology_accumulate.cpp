// Increment 29, PR 1 (docs/increments/29-nve-reference-catchments.md,
// "Accumulation, from the same flood (C++)" and "The red suites", PR 1 C++):
// how many nodes drain through each node, from the same flood `upstream` runs.
//
// THIS IS THE INVARIANT-CRITICAL SUITE OF INCREMENT 29: every sensitivity and
// every station result rests on it. Its oracle is exact and is `upstream`
// itself, called once per node:
//
//   accumulate(z).count[c] == upstream(z, {c}).nodes_in
//   bit 0 of reach[c]      == upstream(z, {c}).touches_edge
//   bit 1 of reach[c]      == upstream(z, {c}).touches_nodata
//   upstream(z, {flow_to[d] as a node}).mask[d] == 1 for every d that has one
//
// The interface (AccumulateOutcome's three arrays, their encodings, and the
// std::length_error refusal at 2^32 nodes) is stated once, in
// include/terrain/hydrology/accumulate.hpp; this suite does not restate it.

#include <catch2/catch_test_macros.hpp>

#include <terrain/hydrology/accumulate.hpp>
#include <terrain/hydrology/upstream.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/view.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

using terrain::hydrology::accumulate;
using terrain::hydrology::AccumulateOutcome;
using terrain::hydrology::upstream;
using terrain::hydrology::UpstreamOutcome;
using terrain::raster::CellIndex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::raster::RasterView;

namespace {

constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
constexpr float kSentinel = -9999.0F;
constexpr std::uint8_t kOutlet = 255;
constexpr std::uint8_t kWest = 3;      // dr = 0, dc = -1
constexpr std::uint8_t kNorthWest = 0; // dr = -1, dc = -1

using Mask = std::vector<std::uint8_t>;

RasterGeometry grid(std::size_t rows, std::size_t cols) {
    // UTM-scale origin and unequal spacing, as upstream's suite: the flood
    // depends on neither.
    return RasterGeometry{500000.0, 6600000.0, 10.0, 5.0, cols, rows};
}

struct Dem {
    std::size_t rows;
    std::size_t cols;
    std::vector<float> z;
    std::optional<float> nodata{};

    [[nodiscard]] std::size_t n() const { return rows * cols; }
    [[nodiscard]] std::size_t at(std::size_t r, std::size_t c) const { return r * cols + c; }
    [[nodiscard]] Raster<float> raster() const { return Raster<float>{grid(rows, cols), z, nodata}; }
    [[nodiscard]] bool invalid(std::size_t i) const {
        return std::isnan(z[i]) || (nodata && z[i] == *nodata);
    }
};

AccumulateOutcome run(const Dem& d) { return accumulate(d.raster()); }

UpstreamOutcome flood_from(const Dem& d, std::size_t node) {
    Mask seed(d.n(), 0);
    seed[node] = 1;
    return upstream(d.raster(), std::span<const std::uint8_t>{seed});
}

bool on_edge(const Dem& d, std::size_t i) {
    const std::size_t r = i / d.cols, c = i % d.cols;
    return r == 0 || c == 0 || r + 1 == d.rows || c + 1 == d.cols;
}

bool beside_nodata(const Dem& d, std::size_t i) {
    const auto r = static_cast<std::ptrdiff_t>(i / d.cols), c = static_cast<std::ptrdiff_t>(i % d.cols);
    for (std::ptrdiff_t dr = -1; dr <= 1; ++dr)
        for (std::ptrdiff_t dc = -1; dc <= 1; ++dc) {
            const std::ptrdiff_t rr = r + dr, cc = c + dc;
            if ((dr == 0 && dc == 0) || rr < 0 || cc < 0 || rr >= static_cast<std::ptrdiff_t>(d.rows)
                || cc >= static_cast<std::ptrdiff_t>(d.cols))
                continue;
            if (d.invalid(static_cast<std::size_t>(rr) * d.cols + static_cast<std::size_t>(cc)))
                return true;
        }
    return false;
}

// The node flow_to[i] names, or nullopt when the byte is not a valid in-window
// direction (4, out of range, or off the raster).
std::optional<std::size_t> target(const Dem& d, std::size_t i, std::uint8_t code) {
    if (code > 8 || code == 4)
        return std::nullopt;
    const auto dr = static_cast<std::ptrdiff_t>(code / 3) - 1;
    const auto dc = static_cast<std::ptrdiff_t>(code % 3) - 1;
    const auto r = static_cast<std::ptrdiff_t>(i / d.cols) + dr;
    const auto c = static_cast<std::ptrdiff_t>(i % d.cols) + dc;
    if (r < 0 || c < 0 || r >= static_cast<std::ptrdiff_t>(d.rows) || c >= static_cast<std::ptrdiff_t>(d.cols))
        return std::nullopt;
    return static_cast<std::size_t>(r) * d.cols + static_cast<std::size_t>(c);
}

void check_shapes(const Dem& d, const AccumulateOutcome& out) {
    REQUIRE(out.count.size() == d.n());
    REQUIRE(out.reach.size() == d.n());
    REQUIRE(out.flow_to.size() == d.n());
}

// flow_to's own rules, conservation, and the tree: everything that needs no
// call to upstream.
void check_tree(const Dem& d, const AccumulateOutcome& out) {
    check_shapes(d, out);
    std::vector<std::uint64_t> children_sum(d.n(), 0);
    std::uint64_t data = 0, outlets_sum = 0;
    for (std::size_t i = 0; i < d.n(); ++i) {
        INFO("node (" << i / d.cols << ", " << i % d.cols << ")");
        if (d.invalid(i)) {
            CHECK(out.count[i] == 0);
            CHECK(out.reach[i] == 0);
            CHECK(out.flow_to[i] == kOutlet);
            continue;
        }
        ++data;
        CHECK(out.count[i] >= 1);
        CHECK(out.reach[i] <= 3); // two bits, nothing else
        const bool outlet = on_edge(d, i) || beside_nodata(d, i);
        if (outlet) {
            CHECK(out.flow_to[i] == kOutlet);
            outlets_sum += out.count[i];
            continue;
        }
        const auto to = target(d, i, out.flow_to[i]);
        REQUIRE(to.has_value()); // a direction, never 4, never off the raster
        REQUIRE_FALSE(d.invalid(*to));
        children_sum[*to] += out.count[i];
    }
    // Conservation: the outlets' counts sum to the nodes with data.
    CHECK(outlets_sum == data);
    for (std::size_t i = 0; i < d.n(); ++i)
        if (!d.invalid(i)) {
            INFO("node (" << i / d.cols << ", " << i % d.cols << ")");
            CHECK(std::uint64_t{out.count[i]} == 1 + children_sum[i]);
        }
    // Following flow_to from any node reaches an outlet: no cycle.
    for (std::size_t i = 0; i < d.n(); ++i) {
        if (d.invalid(i))
            continue;
        std::size_t j = i, steps = 0;
        while (out.flow_to[j] != kOutlet && steps <= d.n()) {
            const auto to = target(d, j, out.flow_to[j]);
            REQUIRE(to.has_value());
            j = *to;
            ++steps;
        }
        INFO("from node (" << i / d.cols << ", " << i % d.cols << ")");
        CHECK(steps <= d.n());
    }
}

// The exact oracle: one upstream call per node with data.
void check_oracle(const Dem& d, const AccumulateOutcome& out) {
    check_shapes(d, out);
    for (std::size_t c = 0; c < d.n(); ++c) {
        if (d.invalid(c))
            continue;
        const UpstreamOutcome u = flood_from(d, c);
        INFO("node (" << c / d.cols << ", " << c % d.cols << ")");
        CHECK(std::size_t{out.count[c]} == u.nodes_in);
        CHECK(((out.reach[c] & 1U) != 0) == u.touches_edge);
        CHECK(((out.reach[c] & 2U) != 0) == u.touches_nodata);
        // Every node that names c as its flooder is in c's catchment.
        for (std::size_t e = 0; e < d.n(); ++e) {
            if (d.invalid(e) || out.flow_to[e] == kOutlet)
                continue;
            const auto to = target(d, e, out.flow_to[e]);
            if (to && *to == c) {
                INFO("drains in: (" << e / d.cols << ", " << e % d.cols << ")");
                CHECK(u.mask[e] == 1);
            }
        }
    }
}

Dem random_dem(int trial, std::mt19937& rng) {
    std::uniform_int_distribution<std::size_t> size(1, 14);
    std::uniform_int_distribution<int> level(0, 4); // few levels: pits, flats, equal heights
    std::bernoulli_distribution hole(0.05);
    Dem d{size(rng), size(rng), {}, std::nullopt};
    d.z.resize(d.n());
    const int kind = trial % 4; // 0 no NoData, 1 NaN holes, 2 sentinel holes, 3 sentinel border
    if (kind >= 2)
        d.nodata = kSentinel;
    for (std::size_t i = 0; i < d.n(); ++i) {
        const bool missing = (kind == 1 || kind == 2) && hole(rng);
        d.z[i] = missing ? (kind == 1 ? kNaN : kSentinel) : static_cast<float>(level(rng));
        if (kind == 3 && on_edge(d, i))
            d.z[i] = kSentinel;
    }
    return d;
}

// z = 10 * (steps to the nearest edge) + a distinct fraction: every interior
// node drains strictly towards the edge, one step at a time.
Dem towards_edge(std::size_t rows, std::size_t cols) {
    Dem d{rows, cols, std::vector<float>(rows * cols), std::nullopt};
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) {
            const std::size_t steps = std::min({r, c, rows - 1 - r, cols - 1 - c});
            // The fraction stays below 1, so steps decide every flooder.
            d.z[d.at(r, c)] = 10.0F * static_cast<float>(steps)
                              + static_cast<float>((r * 7 + c * 3) % 64) / 64.0F;
        }
    return d;
}

// A source that reports a geometry but holds no cells: the size check must
// refuse it before any per-node array is allocated or any cell is read.
struct HugeGeometry {
    using value_type = float;
    RasterGeometry g;
    std::optional<float> none{};
    [[nodiscard]] const RasterGeometry& geometry() const { return g; }
    [[nodiscard]] double value_at(CellIndex) const { throw std::logic_error("HugeGeometry: a cell was read"); }
    [[nodiscard]] bool is_nodata(CellIndex) const { throw std::logic_error("HugeGeometry: a cell was read"); }
    [[nodiscard]] const std::optional<float>& nodata() const { return none; }
    [[nodiscard]] std::span<const float> row(std::size_t) const {
        throw std::logic_error("HugeGeometry: a row was read");
    }
};
static_assert(terrain::raster::RasterSource<HugeGeometry>);

} // namespace

TEST_CASE("the oracle: every count, both bits and every flooder agree with upstream, node for node",
          "[hydrology][accumulate][oracle]") {
    std::mt19937 rng{29001};
    int nodata_cases = 0, deep = 0;
    for (int trial = 0; trial < 300; ++trial) {
        const Dem d = random_dem(trial, rng);
        INFO("trial " << trial << ", " << d.rows << " x " << d.cols);
        const AccumulateOutcome out = run(d);
        check_tree(d, out);
        check_oracle(d, out);
        nodata_cases += std::any_of(out.reach.begin(), out.reach.end(), [](auto b) { return (b & 2U) != 0; });
        deep += *std::max_element(out.count.begin(), out.count.end()) > 3 ? 1 : 0;
    }
    // The oracle must be able to disagree: many counts above one, and the
    // NoData bit set somewhere in many trials.
    CHECK(deep > 150);
    CHECK(nodata_cases > 50);
}

TEST_CASE("a V-valley: the outlet counts the whole valley, and the floor drains west node by node",
          "[hydrology][accumulate][valley]") {
    // upstream's V-valley (test_hydrology_upstream.cpp): 21 x 15, floor along
    // row 10 falling west at 0.37 a column, sides rising 1 a row to ridges at
    // rows 4 and 16; column 0 a wall at 1000 except the outlet (10, 0) at -1.
    constexpr std::size_t rows = 21, cols = 15, rv = 10, w = 6;
    Dem d{rows, cols, std::vector<float>(rows * cols), std::nullopt};
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) {
            const auto off = static_cast<float>(r > rv ? r - rv : rv - r);
            const float side = off <= static_cast<float>(w)
                                   ? off
                                   : static_cast<float>(w) - 2.9F * (off - static_cast<float>(w));
            d.z[d.at(r, c)] = c == 0 ? 1000.0F : 0.37F * static_cast<float>(c) + side;
        }
    d.z[d.at(rv, 0)] = -1.0F;

    const AccumulateOutcome out = run(d);
    check_tree(d, out);
    CHECK(out.count[d.at(rv, 0)] == 11 * 13 + 1);
    CHECK(out.flow_to[d.at(rv, 0)] == kOutlet);
    CHECK((out.reach[d.at(rv, 0)] & 1U) != 0);
    for (std::size_t c = 1; c + 1 < cols; ++c) {
        INFO("floor column " << c);
        CHECK(out.flow_to[d.at(rv, c)] == kWest);
        // Column cols - 1 is an edge outlet of its own, not upstream of the
        // floor, so the comparison stops one column short of it.
        if (c + 2 < cols)
            CHECK(out.count[d.at(rv, c)] > out.count[d.at(rv, c + 1)]);
    }
    check_oracle(d, out);
}

TEST_CASE("a single cell: count 1, an outlet on the edge; a NoData single cell: count 0, no bits",
          "[hydrology][accumulate][degenerate]") {
    const Dem one{1, 1, {7.0F}, std::nullopt};
    const AccumulateOutcome a = run(one);
    REQUIRE(a.count.size() == 1);
    CHECK(a.count[0] == 1);
    CHECK(a.flow_to[0] == kOutlet);
    CHECK(a.reach[0] == 1);

    const Dem missing{1, 1, {kSentinel}, kSentinel};
    const AccumulateOutcome b = run(missing);
    REQUIRE(b.count.size() == 1);
    CHECK(b.count[0] == 0);
    CHECK(b.flow_to[0] == kOutlet);
    CHECK(b.reach[0] == 0);
}

TEST_CASE("a flat plateau draining over one rim node: that node counts the plateau and itself",
          "[hydrology][accumulate][flat]") {
    // 7 x 7: the edge a wall at 100 except (0, 3) at 0; the 5 x 5 interior a
    // flat at 50, so all of it floods from (0, 3) before any wall node pops.
    Dem d{7, 7, std::vector<float>(49, 50.0F), std::nullopt};
    for (std::size_t i = 0; i < d.n(); ++i)
        if (on_edge(d, i))
            d.z[i] = 100.0F;
    d.z[d.at(0, 3)] = 0.0F;
    const AccumulateOutcome out = run(d);
    check_tree(d, out);
    CHECK(out.count[d.at(0, 3)] == 26);
    for (std::size_t i = 0; i < d.n(); ++i)
        if (on_edge(d, i) && i != d.at(0, 3))
            CHECK(out.count[i] == 1);
    check_oracle(d, out);
}

TEST_CASE("1 x n and n x 1 strips: every node an outlet of its own, bit 0 set",
          "[hydrology][accumulate][degenerate]") {
    for (const auto& [rows, cols] : std::array<std::array<std::size_t, 2>, 2>{{{1, 9}, {9, 1}}}) {
        Dem d{rows, cols, std::vector<float>(rows * cols), std::nullopt};
        for (std::size_t i = 0; i < d.n(); ++i)
            d.z[i] = static_cast<float>(i % 3); // ties and rises, deliberately
        const AccumulateOutcome out = run(d);
        INFO(rows << " x " << cols);
        check_tree(d, out);
        CHECK(std::all_of(out.count.begin(), out.count.end(), [](auto v) { return v == 1; }));
        CHECK(std::all_of(out.flow_to.begin(), out.flow_to.end(), [](auto v) { return v == kOutlet; }));
        CHECK(std::all_of(out.reach.begin(), out.reach.end(), [](auto v) { return v == 1; }));
    }
}

TEST_CASE("bit 0: clear on a node two cells from the edge, set on its neighbour one cell from it",
          "[hydrology][accumulate][flags]") {
    // Everything drains straight towards the nearest edge, so (2, 4)'s
    // catchment lies two or more cells in and (1, 4)'s reaches row 1.
    const Dem d = towards_edge(9, 9);
    const AccumulateOutcome out = run(d);
    check_tree(d, out);
    CHECK((out.reach[d.at(2, 4)] & 1U) == 0);
    CHECK((out.reach[d.at(1, 4)] & 1U) != 0);
    CHECK(out.count[d.at(2, 4)] >= 1);
    CHECK((out.reach[d.at(4, 4)] & 1U) == 0); // the centre
    check_oracle(d, out);
}

TEST_CASE("bit 1: set on an outlet beside NoData; a NoData node counts 0 and has no bits",
          "[hydrology][accumulate][flags]") {
    // Bit 1 away from the hole is the oracle's to judge (the random suite).
    Dem d = towards_edge(13, 13);
    d.nodata = kSentinel;
    d.z[d.at(6, 6)] = kSentinel; // a hole in the middle: its 8 neighbours are outlets
    const AccumulateOutcome out = run(d);
    check_tree(d, out);
    CHECK(out.flow_to[d.at(5, 5)] == kOutlet);
    CHECK((out.reach[d.at(5, 6)] & 2U) != 0); // an outlet beside NoData
    CHECK(out.count[d.at(6, 6)] == 0);
    CHECK(out.reach[d.at(6, 6)] == 0);
    check_oracle(d, out);
}

TEST_CASE("a two-node ridge of equal heights: the first-in rule decides each flooder, the same way twice",
          "[hydrology][accumulate][flow_to][tie]") {
    // 3 x 4: the edge outlets all at 0, the two interior nodes (1, 1) and
    // (1, 2) a ridge at 5. Outlets enter the queue in row-major order and pop
    // first in, first out at equal level, so (0, 0) pops first and floods
    // (1, 1), then (0, 1) floods (1, 2): both drain north-west.
    Dem d{3, 4, std::vector<float>(12, 0.0F), std::nullopt};
    d.z[d.at(1, 1)] = 5.0F;
    d.z[d.at(1, 2)] = 5.0F;
    const AccumulateOutcome first = run(d);
    const AccumulateOutcome second = run(d);
    check_tree(d, first);
    CHECK(first.flow_to[d.at(1, 1)] == kNorthWest);
    CHECK(first.flow_to[d.at(1, 2)] == kNorthWest);
    CHECK(first.count[d.at(0, 0)] == 2);
    CHECK(first.count[d.at(0, 1)] == 2);
    CHECK(first.flow_to == second.flow_to);
    CHECK(first.count == second.count);
    check_oracle(d, first);
}

TEST_CASE("determinism: the same input twice gives equal arrays", "[hydrology][accumulate][determinism]") {
    std::mt19937 rng{29002};
    for (int trial = 0; trial < 40; ++trial) {
        const Dem d = random_dem(trial, rng);
        const AccumulateOutcome a = run(d);
        const AccumulateOutcome b = run(d);
        INFO("trial " << trial);
        CHECK(a.count == b.count);
        CHECK(a.reach == b.reach);
        CHECK(a.flow_to == b.flow_to);
    }
}

TEST_CASE("an owning Raster and a RasterView over the same values give the same arrays",
          "[hydrology][accumulate][concept]") {
    std::mt19937 rng{29003};
    Dem d = random_dem(0, rng);
    while (d.n() < 20)
        d = random_dem(0, rng);
    std::vector<double> z64(d.z.begin(), d.z.end());
    const RasterView<double> view{grid(d.rows, d.cols), z64.data(), std::nullopt};
    const AccumulateOutcome from_view = accumulate(view);
    const AccumulateOutcome from_raster = run(d);
    CHECK(from_view.count == from_raster.count);
    CHECK(from_view.reach == from_raster.reach);
    CHECK(from_view.flow_to == from_raster.flow_to);
}

TEST_CASE("2^32 nodes or more is std::length_error, before a cell is read or an array allocated",
          "[hydrology][accumulate][refusal]") {
    constexpr std::size_t k2to16 = std::size_t{1} << 16;
    constexpr std::size_t k2to32 = std::size_t{1} << 32;
    const std::array<HugeGeometry, 4> huge{{
        {grid(k2to16, k2to16)},         // exactly 2^32
        {grid(1, k2to32)},              // one row of 2^32
        {grid(k2to32, 1)},              // one column of 2^32
        {grid(k2to16 + 1, k2to16)},     // just above
    }};
    for (const HugeGeometry& h : huge) {
        INFO(h.g.rows() << " x " << h.g.cols());
        CHECK_THROWS_AS(accumulate(h), std::length_error);
    }
}
