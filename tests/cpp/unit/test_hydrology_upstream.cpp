// Increment 22, PR 1 (docs/increments/22-auto-catchment.md, "Flow and
// membership: one flood" and "The red suites", PR 1 C++): the flood that
// labels every DEM node draining into a seed set.
//
// The design names this file tests/cpp/unit/hydrology_upstream.cpp; it carries
// the test_ prefix every other suite in tests/cpp/unit/ has.
//
// Interface assumed (the design fixes the header, the template on the
// RasterSource concept, and what is returned; the names are chosen here and
// stated in the handback):
//
//   #include <terrain/hydrology/upstream.hpp>
//   namespace terrain::hydrology {
//   struct UpstreamOutcome {
//       std::vector<std::uint8_t> mask; // row-major, geometry().size() bytes:
//                                       // 1 for a node in the catchment, 0 for
//                                       // every other (out, unreached, NoData)
//       std::size_t nodes_in;           // the number of 1s
//       std::size_t row_min, row_max, col_min, col_max; // the in-nodes' bounds,
//                                       // inclusive; meaningful when nodes_in > 0
//       bool touches_edge;              // see below
//       bool touches_nodata;
//   };
//   template <raster::RasterSource R>
//   UpstreamOutcome upstream(const R& z, std::span<const std::uint8_t> seed);
//   }
//
// `seed` is row-major, one byte per node, non-zero for a seed. A span whose
// size is not geometry().size() is std::invalid_argument (the binding's
// ValueError). A NoData seed is ignored: NoData is never in.
//
// THE FLAGS ARE PINNED ONLY WHERE THE DESIGN'S WORDS AND ANY REPAIR OF THEM
// AGREE. As written, every edge node and every NoData-adjacent node is an
// outlet and "outlets are labelled in if they are seeds", so "an in-node on
// the window's edge" can only be a seed; a catchment truncated by the window
// shows as in-nodes one node inside the edge, not on it. That is reported as a
// design gap. This suite pins the cases true under either reading: a seed on
// the edge sets touches_edge; a catchment with no in-node within one node of
// the edge (or within two of NoData) sets neither flag; an in-node on the edge
// (beside NoData) implies the flag. The Python window suite pins the
// behaviour that matters, truncation refused, whatever the flag means.
//
// THE ORACLE for the no-pit, no-tie family is lowest-neighbour descent,
// computed here and sharing nothing with the header: on
//   z = 10 * (steps to the nearest edge) + a distinct fraction in [0, 1)
// every interior node has a strictly lower neighbour, levels equal z, the
// flood pops in global z order, and so each node is flooded by its lowest
// neighbour. The flood's catchment must equal descent's node for node.

#include <catch2/catch_test_macros.hpp>

#include <terrain/hydrology/upstream.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/view.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

using terrain::hydrology::upstream;
using terrain::hydrology::UpstreamOutcome;
using terrain::raster::CellIndex;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using terrain::raster::RasterView;

namespace {

constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();
constexpr float kSentinel = -9999.0F;

using Mask = std::vector<std::uint8_t>;

RasterGeometry grid(std::size_t rows, std::size_t cols) {
    // UTM-scale origin, unequal spacing: nothing here depends on either, and a
    // flood that did would show it.
    return RasterGeometry{500000.0, 6600000.0, 10.0, 5.0, cols, rows};
}

struct Dem {
    std::size_t rows;
    std::size_t cols;
    std::vector<float> z;

    [[nodiscard]] std::size_t at(std::size_t r, std::size_t c) const { return r * cols + c; }
    [[nodiscard]] Raster<float> raster(std::optional<float> nodata = std::nullopt) const {
        return Raster<float>{grid(rows, cols), z, nodata};
    }
};

Mask none(const Dem& d) { return Mask(d.rows * d.cols, 0); }

UpstreamOutcome flood(const Dem& d, const Mask& seed, std::optional<float> nodata = std::nullopt) {
    return upstream(d.raster(nodata), std::span<const std::uint8_t>{seed});
}

bool invalid(const Dem& d, std::size_t i, std::optional<float> nodata) {
    const float v = d.z[i];
    return std::isnan(v) || (nodata && v == *nodata);
}

template <typename F>
void each_neighbour(const Dem& d, std::size_t r, std::size_t c, F&& f) {
    for (int dr = -1; dr <= 1; ++dr)
        for (int dc = -1; dc <= 1; ++dc) {
            if (dr == 0 && dc == 0)
                continue;
            const auto rr = static_cast<std::ptrdiff_t>(r) + dr;
            const auto cc = static_cast<std::ptrdiff_t>(c) + dc;
            if (rr < 0 || cc < 0 || rr >= static_cast<std::ptrdiff_t>(d.rows)
                || cc >= static_cast<std::ptrdiff_t>(d.cols))
                continue;
            f(static_cast<std::size_t>(rr), static_cast<std::size_t>(cc));
        }
}

bool on_edge(const Dem& d, std::size_t r, std::size_t c) {
    return r == 0 || c == 0 || r + 1 == d.rows || c + 1 == d.cols;
}

// ---- the no-pit, no-tie family and its descent oracle ----------------------

Dem no_pit_dem(std::size_t rows, std::size_t cols, std::mt19937& rng) {
    Dem d{rows, cols, std::vector<float>(rows * cols)};
    std::vector<std::size_t> order(rows * cols);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::shuffle(order.begin(), order.end(), rng);
    const auto n = static_cast<float>(rows * cols);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) {
            const std::size_t steps = std::min({r, c, rows - 1 - r, cols - 1 - c});
            // A distinct fraction per node, k / n with k a permutation: exact in
            // float32 for every grid here (n <= 900, steps <= 14).
            d.z[d.at(r, c)] = 10.0F * static_cast<float>(steps)
                              + static_cast<float>(order[d.at(r, c)]) / n;
        }
    return d;
}

// In iff the descent chain from the node (the node itself included) meets a
// seed. Edge nodes end their chain: they are outlets.
Mask descent_oracle(const Dem& d, const Mask& seed) {
    const std::size_t n = d.rows * d.cols;
    std::vector<std::size_t> next(n, n);
    for (std::size_t r = 0; r < d.rows; ++r)
        for (std::size_t c = 0; c < d.cols; ++c) {
            if (on_edge(d, r, c))
                continue;
            float best = std::numeric_limits<float>::infinity();
            each_neighbour(d, r, c, [&](std::size_t rr, std::size_t cc) {
                if (d.z[d.at(rr, cc)] < best) {
                    best = d.z[d.at(rr, cc)];
                    next[d.at(r, c)] = d.at(rr, cc);
                }
            });
            REQUIRE(best < d.z[d.at(r, c)]); // the family's premise: no pit
        }
    Mask in(n, 0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i; j != n; j = next[j])
            if (seed[j]) {
                in[i] = 1;
                break;
            }
    return in;
}

// ---- invariants any outcome must satisfy -----------------------------------

void check_invariants(const Dem& d, const Mask& seed, const UpstreamOutcome& out,
                      std::optional<float> nodata = std::nullopt) {
    const std::size_t n = d.rows * d.cols;
    REQUIRE(out.mask.size() == n);
    std::size_t count = 0;
    std::size_t r0 = d.rows, r1 = 0, c0 = d.cols, c1 = 0;
    bool edge_in = false, near_edge_in = false, beside_nodata_in = false, near_nodata_in = false;
    for (std::size_t r = 0; r < d.rows; ++r)
        for (std::size_t c = 0; c < d.cols; ++c) {
            const std::size_t i = d.at(r, c);
            REQUIRE((out.mask[i] == 0 || out.mask[i] == 1));
            if (invalid(d, i, nodata)) {
                REQUIRE(out.mask[i] == 0); // NoData is never in, seed or not
                continue;
            }
            if (seed[i])
                REQUIRE(out.mask[i] == 1); // every valid seed is in
            if (!out.mask[i])
                continue;
            ++count;
            r0 = std::min(r0, r);
            r1 = std::max(r1, r);
            c0 = std::min(c0, c);
            c1 = std::max(c1, c);
            edge_in = edge_in || on_edge(d, r, c);
            near_edge_in = near_edge_in || r <= 1 || c <= 1 || r + 2 >= d.rows || c + 2 >= d.cols;
            for (std::size_t rr = r > 2 ? r - 2 : 0; rr <= std::min(r + 2, d.rows - 1); ++rr)
                for (std::size_t cc = c > 2 ? c - 2 : 0; cc <= std::min(c + 2, d.cols - 1); ++cc)
                    if (invalid(d, d.at(rr, cc), nodata)) {
                        near_nodata_in = true;
                        const bool adjacent = (rr + 1 >= r && rr <= r + 1 && cc + 1 >= c && cc <= c + 1);
                        beside_nodata_in = beside_nodata_in || adjacent;
                    }
        }
    CHECK(out.nodes_in == count);
    if (count > 0) {
        CHECK(out.row_min == r0);
        CHECK(out.row_max == r1);
        CHECK(out.col_min == c0);
        CHECK(out.col_max == c1);
    }
    if (edge_in)
        CHECK(out.touches_edge);
    if (!near_edge_in)
        CHECK_FALSE(out.touches_edge);
    if (beside_nodata_in)
        CHECK(out.touches_nodata);
    if (!near_nodata_in)
        CHECK_FALSE(out.touches_nodata);

    // Every in-node has an 8-path through in-nodes to a seed.
    Mask reached(n, 0);
    std::deque<std::size_t> queue;
    for (std::size_t i = 0; i < n; ++i)
        if (seed[i] && out.mask[i]) {
            reached[i] = 1;
            queue.push_back(i);
        }
    while (!queue.empty()) {
        const std::size_t i = queue.front();
        queue.pop_front();
        each_neighbour(d, i / d.cols, i % d.cols, [&](std::size_t rr, std::size_t cc) {
            const std::size_t j = d.at(rr, cc);
            if (out.mask[j] && !reached[j]) {
                reached[j] = 1;
                queue.push_back(j);
            }
        });
    }
    CHECK(reached == out.mask);
}

Mask random_seeds(std::size_t n, double density, std::mt19937& rng) {
    std::bernoulli_distribution pick(density);
    Mask seed(n, 0);
    for (auto& s : seed)
        s = pick(rng) ? 1 : 0;
    return seed;
}

Mask blob(const Dem& d, std::mt19937& rng) {
    std::uniform_int_distribution<std::size_t> rr(0, d.rows - 1), cc(0, d.cols - 1);
    std::size_t a = rr(rng), b = rr(rng), e = cc(rng), f = cc(rng);
    Mask seed = none(d);
    for (std::size_t r = std::min(a, b); r <= std::max(a, b); ++r)
        for (std::size_t c = std::min(e, f); c <= std::max(e, f); ++c)
            seed[d.at(r, c)] = 1;
    return seed;
}

} // namespace

TEST_CASE("the flood equals lowest-neighbour descent on DEMs with no pit and no tie",
          "[hydrology][upstream][oracle]") {
    std::mt19937 rng{22001};
    std::uniform_int_distribution<std::size_t> size(3, 30);
    int cases = 0, nontrivial = 0;
    for (int trial = 0; trial < 300; ++trial) {
        const Dem d = no_pit_dem(size(rng), size(rng), rng);
        const std::size_t n = d.rows * d.cols;
        Mask seed = none(d);
        switch (trial % 3) {
        case 0: seed[std::uniform_int_distribution<std::size_t>(0, n - 1)(rng)] = 1; break;
        case 1: seed = blob(d, rng); break;
        default: seed = random_seeds(n, 0.05, rng); break;
        }
        const Mask expected = descent_oracle(d, seed);
        const UpstreamOutcome out = flood(d, seed);
        INFO("trial " << trial << ", " << d.rows << " x " << d.cols);
        REQUIRE(out.mask == expected);
        check_invariants(d, seed, out);
        ++cases;
        const auto seeds = static_cast<std::size_t>(std::count(seed.begin(), seed.end(), 1));
        nontrivial += out.nodes_in > seeds ? 1 : 0;
    }
    // The oracle must be able to disagree: most catchments are larger than
    // their seed set, so a flood that returned the seeds alone would fail.
    CHECK(cases == 300);
    CHECK(nontrivial > 150);
}

TEST_CASE("a closed depression is wholly in when it spills towards the seed, wholly out otherwise",
          "[hydrology][upstream][depression]") {
    // Every edge node is a high wall (100) except two outlets at 0: the seed
    // (7, 0) on the west edge, and a plain outlet (3, 14) on the east. The
    // interior is the lower of two funnels, one towards each.
    Dem d{15, 15, std::vector<float>(225)};
    for (std::size_t r = 0; r < 15; ++r)
        for (std::size_t c = 0; c < 15; ++c) {
            const float dr7 = std::abs(static_cast<float>(r) - 7.0F);
            const float dr3 = std::abs(static_cast<float>(r) - 3.0F);
            const float west = 10.0F + 2.0F * dr7 + static_cast<float>(c);
            const float east = 10.0F + 2.0F * dr3 + static_cast<float>(14 - c) + 0.5F;
            d.z[d.at(r, c)] = on_edge(d, r, c) ? 100.0F : std::min(west, east);
        }
    d.z[d.at(7, 0)] = 0.0F;
    d.z[d.at(3, 14)] = 0.0F;
    // Bowl A, 2 x 2 at rows 10-11, cols 4-5, in the west funnel (west 20-23,
    // east above 33): a pit well below its rim.
    for (std::size_t r : {10U, 11U})
        for (std::size_t c : {4U, 5U})
            d.z[d.at(r, c)] = 1.0F;
    // Bowl B, 2 x 3 at rows 2-3, cols 10-12, in the east funnel.
    for (std::size_t r : {2U, 3U})
        for (std::size_t c : {10U, 11U, 12U})
            d.z[d.at(r, c)] = 1.0F;

    Mask seed = none(d);
    seed[d.at(7, 0)] = 1;
    const UpstreamOutcome out = flood(d, seed);
    check_invariants(d, seed, out);
    for (std::size_t r : {10U, 11U})
        for (std::size_t c : {4U, 5U})
            CHECK(out.mask[d.at(r, c)] == 1);
    for (std::size_t r : {2U, 3U})
        for (std::size_t c : {10U, 11U, 12U})
            CHECK(out.mask[d.at(r, c)] == 0);
    CHECK(out.mask[d.at(3, 14)] == 0); // the other outlet
    CHECK(out.touches_edge);           // the seed is on the west edge
}

TEST_CASE("a flat lake with shelves: a shelf that spills into it is in, one that spills away is out",
          "[hydrology][upstream][flat]") {
    // 24 x 24. The lake is the flat block rows 10-13, cols 10-13, at 50. Around
    // it, by Chebyshev distance k from the block: a bowl rising to a rim at
    // k = 3 (70, 80, 90), then a slope falling to the edges (90 - 5 (k - 3),
    // 55 on the edge). The lake's outlet is a channel at 40 along row 11 from
    // the lake to the east edge, so the lake is not itself a depression.
    Dem d{24, 24, std::vector<float>(576)};
    const auto k_of = [](std::size_t r, std::size_t c) {
        const auto axis = [](std::size_t i) {
            return i < 10 ? 10 - i : (i > 13 ? i - 13 : std::size_t{0});
        };
        return std::max(axis(r), axis(c));
    };
    for (std::size_t r = 0; r < 24; ++r)
        for (std::size_t c = 0; c < 24; ++c) {
            const std::size_t k = k_of(r, c);
            d.z[d.at(r, c)] = k == 0   ? 50.0F
                              : k <= 3 ? 60.0F + 10.0F * static_cast<float>(k)
                                       : 90.0F - 5.0F * static_cast<float>(k - 3);
        }
    for (std::size_t c = 14; c < 24; ++c)
        d.z[d.at(11, c)] = 40.0F;
    // Shelf N: a flat at 60 outside the rim (k = 4, row 6, cols 8-15), closed
    // by the slope (80 at row 5) except for a channel at 55 through the rim
    // (col 12, rows 7-9) down to the lake: it spills into the lake.
    for (std::size_t c = 8; c <= 15; ++c)
        d.z[d.at(6, c)] = 60.0F;
    for (std::size_t r = 7; r <= 9; ++r)
        d.z[d.at(r, 12)] = 55.0F;
    // Shelf S: a flat at 60 outside the rim on the south (k = 5, row 18, cols
    // 8-15), whose lowest way out is outwards (k = 6, 75), not inwards (k = 4,
    // 85): it spills away.
    for (std::size_t c = 8; c <= 15; ++c)
        d.z[d.at(18, c)] = 60.0F;

    Mask seed = none(d);
    for (std::size_t r = 10; r <= 13; ++r)
        for (std::size_t c = 10; c <= 13; ++c)
            seed[d.at(r, c)] = 1;
    const UpstreamOutcome out = flood(d, seed);
    check_invariants(d, seed, out);

    for (std::size_t c = 8; c <= 15; ++c) {
        CHECK(out.mask[d.at(6, c)] == 1);
        CHECK(out.mask[d.at(18, c)] == 0);
    }
    for (std::size_t r = 7; r <= 9; ++r)
        CHECK(out.mask[d.at(r, 12)] == 1);
    // The outlet channel is downstream of the lake: out.
    for (std::size_t c = 14; c < 24; ++c)
        CHECK(out.mask[d.at(11, c)] == 0);
    // The bowl inside the rim drains to the lake, away from the channel's mouth.
    for (std::size_t r = 0; r < 24; ++r)
        for (std::size_t c = 0; c <= 12; ++c)
            if (k_of(r, c) <= 2)
                CHECK(out.mask[d.at(r, c)] == 1);
    // Far from every edge and with no NoData: no flag.
    CHECK_FALSE(out.touches_edge);
    CHECK_FALSE(out.touches_nodata);
}

TEST_CASE("a V-valley: the catchment is its side slopes up to the ridge lines, exactly",
          "[hydrology][upstream][valley]") {
    // 21 x 15. Valley along row 10, falling west at b = 0.37 per column, sides
    // rising at a = 1 per row to ridges at rows 4 and 16, outer slopes falling
    // at 2.9 per row to the north and south edges. Column 0 is a wall (1000)
    // except the seed (10, 0) at -1, the valley's outlet.
    constexpr std::size_t rows = 21, cols = 15, rv = 10, w = 6;
    Dem d{rows, cols, std::vector<float>(rows * cols)};
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) {
            const auto off = static_cast<float>(r > rv ? r - rv : rv - r);
            const float side = off <= static_cast<float>(w)
                                   ? off
                                   : static_cast<float>(w) - 2.9F * (off - static_cast<float>(w));
            d.z[d.at(r, c)] = c == 0 ? 1000.0F : 0.37F * static_cast<float>(c) + side;
        }
    d.z[d.at(rv, 0)] = -1.0F;

    Mask seed = none(d);
    seed[d.at(rv, 0)] = 1;
    const UpstreamOutcome out = flood(d, seed);
    check_invariants(d, seed, out);

    Mask expected = none(d);
    expected[d.at(rv, 0)] = 1;
    for (std::size_t r = rv - (w - 1); r <= rv + (w - 1); ++r)
        for (std::size_t c = 1; c + 1 < cols; ++c)
            expected[d.at(r, c)] = 1;
    CHECK(out.mask == expected);
    CHECK(out.nodes_in == 11 * 13 + 1);
    CHECK(out.row_min == 5);
    CHECK(out.row_max == 15);
    CHECK(out.col_min == 0);
    CHECK(out.col_max == 13);
    CHECK(out.touches_edge);
}

TEST_CASE("NoData is never in, a NoData seed is ignored, and a seed beside NoData sets the flag",
          "[hydrology][upstream][nodata]") {
    std::mt19937 rng{22002};
    Dem d = no_pit_dem(12, 12, rng);
    d.z[d.at(5, 5)] = kNaN;
    d.z[d.at(8, 3)] = kSentinel;

    SECTION("NaN and the sentinel, as seeds, are not in") {
        Mask seed = none(d);
        seed[d.at(5, 5)] = 1;
        seed[d.at(8, 3)] = 1;
        const UpstreamOutcome out = flood(d, seed, kSentinel);
        CHECK(out.nodes_in == 0);
        CHECK(std::all_of(out.mask.begin(), out.mask.end(), [](auto m) { return m == 0; }));
        CHECK_FALSE(out.touches_nodata);
        CHECK_FALSE(out.touches_edge);
    }
    SECTION("without a sentinel, -9999 is a value like any other") {
        Mask seed = none(d);
        seed[d.at(8, 3)] = 1;
        const UpstreamOutcome out = flood(d, seed);
        CHECK(out.mask[d.at(8, 3)] == 1);
    }
    SECTION("a seed beside NoData is in and sets touches_nodata") {
        Mask seed = none(d);
        seed[d.at(5, 6)] = 1;
        const UpstreamOutcome out = flood(d, seed, kSentinel);
        CHECK(out.mask[d.at(5, 6)] == 1);
        CHECK(out.mask[d.at(5, 5)] == 0);
        CHECK(out.touches_nodata);
        check_invariants(d, seed, out, kSentinel);
    }
    SECTION("every valid node seeded: every valid node in, no NoData in") {
        Mask seed(d.rows * d.cols, 1);
        const UpstreamOutcome out = flood(d, seed, kSentinel);
        CHECK(out.nodes_in == d.rows * d.cols - 2);
        check_invariants(d, seed, out, kSentinel);
    }
}

TEST_CASE("a seed on the window's edge is allowed and sets touches_edge; the bounds are the in-nodes'",
          "[hydrology][upstream][edge]") {
    std::mt19937 rng{22003};
    const Dem d = no_pit_dem(9, 11, rng);
    Mask seed = none(d);
    seed[d.at(0, 4)] = 1;
    const UpstreamOutcome out = flood(d, seed);
    CHECK(out.mask[d.at(0, 4)] == 1);
    CHECK(out.touches_edge);
    CHECK(out.row_min == 0);
    check_invariants(d, seed, out);
    // An edge node is an outlet: nothing drains through it, so on this family
    // the seed on the edge is its own whole catchment.
    CHECK(out.mask == descent_oracle(d, seed));
}

TEST_CASE("no seed: an empty catchment and no flag", "[hydrology][upstream][empty]") {
    std::mt19937 rng{22004};
    const Dem d = no_pit_dem(7, 6, rng);
    const UpstreamOutcome out = flood(d, none(d));
    CHECK(out.nodes_in == 0);
    CHECK(out.mask == none(d));
    CHECK_FALSE(out.touches_edge);
    CHECK_FALSE(out.touches_nodata);
}

TEST_CASE("a seed mask of the wrong size is std::invalid_argument", "[hydrology][upstream][refusal]") {
    std::mt19937 rng{22005};
    const Dem d = no_pit_dem(5, 6, rng);
    const Raster<float> raster = d.raster();
    const Mask short_by_one(29, 0), long_by_one(31, 0), empty;
    CHECK_THROWS_AS(upstream(raster, std::span<const std::uint8_t>{short_by_one}),
                    std::invalid_argument);
    CHECK_THROWS_AS(upstream(raster, std::span<const std::uint8_t>{long_by_one}),
                    std::invalid_argument);
    CHECK_THROWS_AS(upstream(raster, std::span<const std::uint8_t>{empty}), std::invalid_argument);
}

TEST_CASE("one-row, one-column and single-node rasters: every node is an outlet, the seeds are the catchment",
          "[hydrology][upstream][degenerate]") {
    constexpr std::array<std::array<std::size_t, 2>, 4> shapes{{{1, 1}, {1, 7}, {7, 1}, {2, 2}}};
    for (const auto& [rows, cols] : shapes) {
        Dem d{rows, cols, std::vector<float>(rows * cols)};
        for (std::size_t i = 0; i < d.z.size(); ++i)
            d.z[i] = static_cast<float>(i % 3); // ties and a rise, deliberately
        Mask seed = none(d);
        seed[d.z.size() / 2] = 1;
        const UpstreamOutcome out = flood(d, seed);
        INFO(rows << " x " << cols);
        CHECK(out.mask == seed);
        CHECK(out.touches_edge);
    }
}

TEST_CASE("invariants on random DEMs with pits, flats and NoData", "[hydrology][upstream][property]") {
    std::mt19937 rng{22006};
    std::uniform_int_distribution<std::size_t> size(1, 24);
    std::uniform_int_distribution<int> level(0, 4); // few levels: flats, ties, pits
    std::bernoulli_distribution hole(0.04);
    for (int trial = 0; trial < 400; ++trial) {
        Dem d{size(rng), size(rng), {}};
        d.z.resize(d.rows * d.cols);
        const bool with_nodata = trial % 2 == 1;
        for (auto& v : d.z)
            v = with_nodata && hole(rng) ? (trial % 4 == 1 ? kNaN : kSentinel)
                                         : static_cast<float>(level(rng));
        const std::optional<float> nodata =
            with_nodata ? std::optional<float>{kSentinel} : std::nullopt;
        const std::size_t n = d.rows * d.cols;
        const Mask s1 = random_seeds(n, 0.03, rng);
        const Mask s2 = blob(d, rng);
        Mask both(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            both[i] = (s1[i] || s2[i]) ? 1 : 0;

        INFO("trial " << trial << ", " << d.rows << " x " << d.cols);
        const UpstreamOutcome a = flood(d, s1, nodata);
        const UpstreamOutcome b = flood(d, s2, nodata);
        const UpstreamOutcome ab = flood(d, both, nodata);
        check_invariants(d, s1, a, nodata);
        check_invariants(d, s2, b, nodata);
        check_invariants(d, both, ab, nodata);

        // Deterministic: the same input twice gives the same mask and counts.
        const UpstreamOutcome again = flood(d, s1, nodata);
        CHECK(again.mask == a.mask);
        CHECK(again.nodes_in == a.nodes_in);

        // The flood's order does not depend on the labels, so catchments add:
        // the catchment of a union of seed sets is the union of catchments.
        Mask united(n, 0);
        for (std::size_t i = 0; i < n; ++i)
            united[i] = (a.mask[i] || b.mask[i]) ? 1 : 0;
        CHECK(ab.mask == united);
    }
}

TEST_CASE("an owning Raster and a RasterView over the same buffer give the same flood",
          "[hydrology][upstream][concept]") {
    std::mt19937 rng{22007};
    const Dem d = no_pit_dem(13, 17, rng);
    std::vector<double> z64(d.z.begin(), d.z.end());
    const RasterView<double> view{grid(d.rows, d.cols), z64.data(), std::nullopt};
    const Mask seed = blob(d, rng);
    const UpstreamOutcome from_view = upstream(view, std::span<const std::uint8_t>{seed});
    const UpstreamOutcome from_raster = flood(d, seed);
    CHECK(from_view.mask == from_raster.mask);
    CHECK(from_view.nodes_in == from_raster.nodes_in);
}
