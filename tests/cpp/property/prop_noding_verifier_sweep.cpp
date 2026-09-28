// Property suite for increment 16b-0: the verifier's pair search by sort and
// sweep (docs/increments/16b-terrain-polygons.md, R8, I6 and "Pinned by the red
// suite (16b-0)").
//
// INVARIANT-CRITICAL, mutation round required. Three layers, each able to fail
// where the others cannot:
//
//  A. sweep_box_pairs, the pinned kernel-free pair search, against
//     verifier_oracle::brute_force_box_pairs: EXACTLY the pairs of closed boxes
//     that intersect, each reported once, no cell-cell pair. This is the only
//     layer that can see a sweep which drops boxes touching at one x: on-grid
//     input never gives a violating pair whose boxes meet only at an x (two
//     segments sharing just that x meet at an endpoint, which classify calls
//     Touching, and a padded cell's x-bounds sit at half-integer grid offsets).
//  B. NodedPslgBuilder::build<K>() against brute_force_guarantee_14 (I6): the
//     same status on random and adversarial candidates, built by hand as in
//     test_noding_noded_pslg_builder.cpp, plus candidates made by adding one
//     edge or one node to a real noder output.
//  C. node<K> unchanged: its output on integer fixtures at a dyadic spacing,
//     digested over integers only, equals the digest recorded from today's
//     brute-force verifier (6b4fcb9) before any production change; and every
//     Ok output of node<K> on random polylines passes the brute force.
//
// The existing 14(a) and 14(b) refusal fixtures stay in
// test_noding_noded_pslg_builder.cpp, unchanged.
//
// DefaultKernel only, for the reason prop_noding_no_crossings.cpp gives.
// Random inputs draw raw std::mt19937_64 words, never a std::*_distribution,
// whose output the standard does not fix across libraries.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>

#include "noding_generators.h"
#include "verifier_oracle.hpp"

#include <terrain/core/bbox.hpp>
#include <terrain/core/edge_properties.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/node.hpp>
#include <terrain/noding/noded_pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <span>
#include <utility>
#include <vector>

using terrain::Box2;
using terrain::Chain;
using terrain::ChainRole;
using terrain::EdgeProperties;
using terrain::GridPoint;
using terrain::Point2;
using terrain::SnapGrid;
using terrain::noding::NodedPslgBuilder;
using terrain::noding::NodeOptions;
using terrain::noding::NodeStatus;
using terrain::noding::sweep_box_pairs;
using terrain::pred::DefaultKernel;
using verifier_oracle::BoxPairs;
using verifier_oracle::brute_force_box_pairs;

namespace {

// ===========================================================================
// A. sweep_box_pairs
// ===========================================================================

struct SweepRun {
    BoxPairs pairs;
    std::size_t calls{0};
    bool repeated{false};      // a pair reported twice
    bool malformed{false};     // i == j, or an index out of range
    bool completed{false};     // the return value
};

// Runs the sweep, recording every call. stop_after == 0 never stops; k > 0
// makes the k-th call (of either kind) return false.
[[nodiscard]] SweepRun run_sweep(std::span<const Box2> edges, std::span<const Box2> cells,
                                 std::size_t stop_after = 0) {
    SweepRun run;
    auto keep_going = [&run, stop_after] { return stop_after == 0 || run.calls < stop_after; };
    run.completed = sweep_box_pairs(
        edges, cells,
        [&](std::size_t i, std::size_t j) {
            ++run.calls;
            if (i == j || i >= edges.size() || j >= edges.size()) {
                run.malformed = true;
            }
            if (!run.pairs.edge_edge.emplace(std::min(i, j), std::max(i, j)).second) {
                run.repeated = true;
            }
            return keep_going();
        },
        [&](std::size_t e, std::size_t c) {
            ++run.calls;
            if (e >= edges.size() || c >= cells.size()) {
                run.malformed = true;
            }
            if (!run.pairs.edge_cell.emplace(e, c).second) {
                run.repeated = true;
            }
            return keep_going();
        });
    return run;
}

void require_exact_pairs(const std::vector<Box2>& edges, const std::vector<Box2>& cells) {
    const BoxPairs want = brute_force_box_pairs(edges, cells);
    const SweepRun got = run_sweep(edges, cells);
    INFO(edges.size() << " edge boxes, " << cells.size() << " cell boxes; oracle "
                      << want.edge_edge.size() << " edge pairs, " << want.edge_cell.size()
                      << " edge-cell pairs; sweep " << got.pairs.edge_edge.size() << " and "
                      << got.pairs.edge_cell.size());
    REQUIRE(got.completed);
    REQUIRE_FALSE(got.malformed);
    REQUIRE_FALSE(got.repeated);
    REQUIRE(got.pairs.edge_edge == want.edge_edge);
    REQUIRE(got.pairs.edge_cell == want.edge_cell);
    REQUIRE(got.calls == want.edge_edge.size() + want.edge_cell.size());
}

[[nodiscard]] Box2 box(double x0, double y0, double x1, double y1) {
    return Box2{Point2{x0, y0}, Point2{x1, y1}};
}

// Random boxes on a small integer lattice, so that equal xmin, xmax == another's
// xmin, zero width, zero height and corner contact are common rather than rare.
// A tenth are long: they span most of the range, which is the case that keeps
// an item active across every other arrival.
[[nodiscard]] std::vector<Box2> random_boxes(std::mt19937_64& rng, std::size_t n,
                                             std::uint64_t range, double scale,
                                             double offset) {
    std::vector<Box2> out;
    out.reserve(n);
    auto extent = [&rng, range]() -> std::uint64_t {
        switch (rng() % 10) {
            case 0: return 0;                   // zero width or height
            case 1: return range - rng() % 3;   // long
            default: return rng() % 4;
        }
    };
    for (std::size_t k = 0; k < n; ++k) {
        const std::uint64_t x0 = rng() % range;
        const std::uint64_t y0 = rng() % range;
        const std::uint64_t x1 = x0 + extent();
        const std::uint64_t y1 = y0 + extent();
        out.push_back(box(static_cast<double>(x0) * scale + offset,
                          static_cast<double>(y0) * scale + offset,
                          static_cast<double>(x1) * scale + offset,
                          static_cast<double>(y1) * scale + offset));
    }
    return out;
}

TEST_CASE("sweep_box_pairs: nothing to pair", "[16b-0][sweep]") {
    const std::vector<Box2> none;
    const std::vector<Box2> one{box(0, 0, 1, 1)};
    const std::vector<Box2> cells{box(0, 0, 1, 1), box(0, 0, 1, 1), box(0.5, 0.5, 2, 2)};

    SECTION("no boxes at all") { require_exact_pairs(none, none); }
    SECTION("one edge box") { require_exact_pairs(one, none); }
    SECTION("cells only: overlapping cells are never paired with each other") {
        const SweepRun run = run_sweep(none, cells);
        REQUIRE(run.completed);
        REQUIRE(run.calls == 0);
    }
}

TEST_CASE("sweep_box_pairs: closed boxes touching in one coordinate are a pair",
          "[16b-0][sweep]") {
    // Each fixture is given in both index orders: which of the two arrives first
    // in the sweep must not matter.
    auto both_orders = [](Box2 a, Box2 b) {
        require_exact_pairs({a, b}, {});
        require_exact_pairs({b, a}, {});
        require_exact_pairs({a}, {b});
        require_exact_pairs({b}, {a});
    };

    SECTION("an edge's xmax equal to the other's xmin") {
        both_orders(box(0, 0, 10, 1), box(10, 0, 20, 1));
    }
    SECTION("two zero-width boxes at the same x, overlapping in y") {
        both_orders(box(5, 0, 5, 4), box(5, 2, 5, 6));
    }
    SECTION("a zero-width box on the other's right face") {
        both_orders(box(0, 0, 10, 10), box(10, 3, 10, 4));
    }
    SECTION("ymax equal to the other's ymin, x-intervals overlapping") {
        both_orders(box(0, 0, 10, 5), box(3, 5, 4, 9));
    }
    SECTION("corners only") { both_orders(box(0, 0, 1, 1), box(1, 1, 2, 2)); }
    SECTION("a point box on a corner") { both_orders(box(0, 0, 1, 1), box(1, 0, 1, 0)); }
    SECTION("two coincident point boxes") { both_orders(box(3, 3, 3, 3), box(3, 3, 3, 3)); }
    SECTION("the same at huge magnitude") {
        const double big = 0x1p52;
        both_orders(box(-big, big, big, big), box(big, big, 2 * big, 2 * big));
    }
}

TEST_CASE("sweep_box_pairs: boxes one ulp apart, or apart in only one coordinate, are not "
          "a pair",
          "[16b-0][sweep]") {
    const double after10 = std::nextafter(10.0, 20.0);
    auto neither_order = [](Box2 a, Box2 b) {
        require_exact_pairs({a, b}, {});
        require_exact_pairs({b, a}, {});
        require_exact_pairs({a}, {b});
        require_exact_pairs({b}, {a});
    };

    SECTION("one ulp apart in x") { neither_order(box(0, 0, 10, 1), box(after10, 0, 20, 1)); }
    SECTION("one ulp apart in y") { neither_order(box(0, 0, 1, 10), box(0, after10, 1, 20)); }
    // An item kept active past its xmax, or a sweep that skips the y-test, pairs
    // these.
    SECTION("x-disjoint, y-equal") { neither_order(box(0, 0, 1, 1), box(2, 0, 3, 1)); }
    SECTION("y-disjoint under a long box's x-range") {
        neither_order(box(0, 0, 100, 1), box(50, 2, 51, 3));
    }
}

TEST_CASE("sweep_box_pairs: many boxes at one x", "[16b-0][sweep]") {
    // Twenty zero-width boxes on x = 5, end to end and overlapping, and twenty
    // cells straddling x = 5. Equal xmin is a tie in the sort key everywhere.
    std::vector<Box2> edges;
    std::vector<Box2> cells;
    for (int k = 0; k < 20; ++k) {
        edges.push_back(box(5, k, 5, k + 1 + (k % 3)));
        cells.push_back(box(5 - (k % 2), k - 0.5, 5 + (k % 2), k + 0.5));
    }
    require_exact_pairs(edges, cells);
    std::ranges::reverse(edges);
    require_exact_pairs(edges, cells);
}

TEST_CASE("sweep_box_pairs: a long box spanning every other one", "[16b-0][sweep]") {
    std::vector<Box2> edges{box(0, 0, 1000, 1000)};
    std::vector<Box2> cells;
    for (int k = 0; k < 50; ++k) {
        edges.push_back(box(20 * k, 20 * k, 20 * k + 5, 20 * k + 5));
        cells.push_back(box(20 * k + 5, 20 * k + 5, 20 * k + 10, 20 * k + 10));
    }
    edges.push_back(box(-5, 1000, 1005, 1000));  // along the top face
    require_exact_pairs(edges, cells);
}

TEST_CASE("sweep_box_pairs: an edge arriving after a cell and a cell arriving after an "
          "edge are both paired",
          "[16b-0][sweep]") {
    // The two directions are separate code paths in R8: an arriving edge against
    // active cells, an arriving cell against active edges.
    require_exact_pairs({box(5, 0, 9, 1)}, {box(0, 0, 6, 1)});   // cell first
    require_exact_pairs({box(0, 0, 6, 1)}, {box(5, 0, 9, 1)});   // edge first
    require_exact_pairs({box(5, 0, 9, 1)}, {box(5, 0, 6, 1)});   // same xmin
    require_exact_pairs({box(6, 0, 9, 1)}, {box(0, 0, 6, 1)});   // cell xmax == edge xmin
    require_exact_pairs({box(0, 0, 6, 1)}, {box(6, 0, 9, 1)});   // edge xmax == cell xmin
}

TEST_CASE("sweep_box_pairs reports exactly the brute force's pairs on random boxes",
          "[16b-0][sweep][random]") {
    struct Frame {
        double scale;
        double offset;
    };
    // Unit, a non-dyadic scale, huge positive and huge negative coordinates.
    const Frame frame = GENERATE(Frame{1.0, 0.0}, Frame{0.1, 0.0}, Frame{1.0, 0x1p51},
                                 Frame{0.1, -0x1p49});
    const std::uint64_t seed = GENERATE(range<std::uint64_t>(1, 41));
    std::mt19937_64 rng{seed};
    const std::uint64_t range_ = 4 + rng() % 40;
    const std::size_t n_edges = rng() % 60;
    const std::size_t n_cells = rng() % 60;
    const std::vector<Box2> edges = random_boxes(rng, n_edges, range_, frame.scale, frame.offset);
    const std::vector<Box2> cells = random_boxes(rng, n_cells, range_, frame.scale, frame.offset);
    INFO("seed " << seed << ", range " << range_);
    require_exact_pairs(edges, cells);
}

TEST_CASE("sweep_box_pairs stops at the first callback that returns false",
          "[16b-0][sweep]") {
    std::mt19937_64 rng{7};
    const std::vector<Box2> edges = random_boxes(rng, 30, 10, 1.0, 0.0);
    const std::vector<Box2> cells = random_boxes(rng, 30, 10, 1.0, 0.0);
    const std::size_t total = run_sweep(edges, cells).calls;
    REQUIRE(total > 20);  // the probe must have something to stop

    for (std::size_t k = 1; k <= total; ++k) {
        INFO("stop at call " << k << " of " << total);
        const SweepRun run = run_sweep(edges, cells, k);
        REQUIRE_FALSE(run.completed);
        REQUIRE(run.calls == k);
    }
}

TEST_CASE("sweep_box_pairs is a pure function: the same calls in the same order",
          "[16b-0][sweep]") {
    std::mt19937_64 rng{11};
    const std::vector<Box2> edges = random_boxes(rng, 40, 12, 1.0, 0.0);
    const std::vector<Box2> cells = random_boxes(rng, 40, 12, 1.0, 0.0);
    auto sequence = [&] {
        std::vector<std::pair<std::size_t, std::size_t>> seq;
        (void)sweep_box_pairs(
            std::span<const Box2>{edges}, std::span<const Box2>{cells},
            [&](std::size_t i, std::size_t j) { seq.emplace_back(i, j); return true; },
            [&](std::size_t e, std::size_t c) { seq.emplace_back(e, 1000 + c); return true; });
        return seq;
    };
    REQUIRE(sequence() == sequence());
}

// ===========================================================================
// B. NodedPslgBuilder::build<K>() against the brute force (I6)
// ===========================================================================

struct Candidate {
    SnapGrid grid{0.1};
    std::vector<Point2> vertices;
    std::vector<Chain> chains;
    std::vector<std::uint32_t> chain_indices;
    std::vector<EdgeProperties> edge_properties;
    std::vector<std::uint32_t> node_of_input_vertex;

    // Appends a chain over existing vertex ids, keeping 16's structure.
    void add_chain(const std::vector<std::uint32_t>& ids, ChainRole role = ChainRole::Breakline,
                   EdgeProperties props = EdgeProperties::bit(0)) {
        chains.push_back(Chain{static_cast<std::uint32_t>(chain_indices.size()),
                               static_cast<std::uint32_t>(ids.size()), role, props});
        chain_indices.insert(chain_indices.end(), ids.begin(), ids.end());
        const std::size_t edges = terrain::is_closed(role) ? ids.size() : ids.size() - 1;
        edge_properties.insert(edge_properties.end(), edges, props);
    }

    // Returns the id of g's vertex, adding it if it is new.
    std::uint32_t vertex(GridPoint g) {
        const Point2 p = grid.world(g);
        for (std::uint32_t k = 0; k < vertices.size(); ++k) {
            if (grid.snap(vertices[k]) == g) {
                return k;
            }
        }
        vertices.push_back(p);
        return static_cast<std::uint32_t>(vertices.size() - 1);
    }

    [[nodiscard]] NodeStatus built() const {
        return NodedPslgBuilder{grid, vertices, chains, chain_indices, edge_properties,
                                node_of_input_vertex}
            .build<DefaultKernel>()
            .status;
    }

    [[nodiscard]] NodeStatus oracle() const {
        return verifier_oracle::brute_force_guarantee_14<DefaultKernel>(grid, vertices, chains,
                                                                        chain_indices);
    }
};

// The first edge of chain `chain` in the flat edge numbering.
[[nodiscard]] std::size_t base_edge_offset(const Candidate& c, std::size_t chain) {
    std::size_t e = 0;
    for (std::size_t k = 0; k < chain; ++k) {
        const Chain& ch = c.chains[k];
        e += terrain::is_closed(ch.role) ? ch.count : ch.count - 1;
    }
    return e;
}

// The status must be the oracle's, and the fixture's own expectation, so a
// fixture that stops meaning what its name says fails here rather than quietly
// comparing two agreeing answers to the wrong question.
void require_verdict(const Candidate& c, NodeStatus expected) {
    REQUIRE(c.oracle() == expected);
    REQUIRE(c.built() == expected);
}

enum class Layout { Free, AxisLines, Diagonals };

// A random candidate in grid space, then mapped to the world through
// world(g + origin). Structure (16), 11, 12 and 13 hold by construction, so the
// only status build() can return is 14's.
[[nodiscard]] Candidate random_candidate(std::uint64_t seed, double spacing,
                                         std::int64_t origin, Layout layout) {
    std::mt19937_64 rng{seed};
    Candidate c;
    c.grid = SnapGrid{spacing};
    const std::int64_t w = 4 + static_cast<std::int64_t>(rng() % 12);

    auto pick = [&]() -> GridPoint {
        const auto a = static_cast<std::int64_t>(rng() % static_cast<std::uint64_t>(w));
        const auto b = static_cast<std::int64_t>(rng() % static_cast<std::uint64_t>(w));
        switch (layout) {
            case Layout::Free: return {a, b};
            case Layout::AxisLines: return rng() % 2 ? GridPoint{w / 2, b} : GridPoint{a, w / 3};
            case Layout::Diagonals: return rng() % 2 ? GridPoint{a, a} : GridPoint{a, w - 1 - a};
        }
        return {a, b};
    };

    const std::size_t n_vertices = 3 + rng() % 10;
    for (std::size_t k = 0; k < n_vertices; ++k) {
        const GridPoint g = pick();
        (void)c.vertex(GridPoint{g.ix + origin, g.iy + origin});
    }
    const auto n = static_cast<std::uint32_t>(c.vertices.size());
    if (n < 3) {
        return c;
    }

    const std::size_t n_chains = 1 + rng() % 4;
    for (std::size_t k = 0; k < n_chains; ++k) {
        const bool closed = rng() % 5 == 0;
        const std::size_t len = (closed ? 3 : 2) + rng() % 3;
        std::vector<std::uint32_t> ids{static_cast<std::uint32_t>(rng() % n)};
        while (ids.size() < len) {
            const auto next = static_cast<std::uint32_t>(rng() % n);
            if (next != ids.back()) {
                ids.push_back(next);
            }
        }
        if (closed && ids.back() == ids.front()) {
            ids.back() = (ids.front() + 1) % n == ids[ids.size() - 2]
                             ? (ids.front() + 2) % n
                             : (ids.front() + 1) % n;
        }
        if (closed && (ids.back() == ids.front() || ids.back() == ids[ids.size() - 2])) {
            continue;  // n == 3 corner cases: skip rather than emit a zero-length edge
        }
        c.add_chain(ids, closed ? ChainRole::Hole : ChainRole::Breakline,
                    EdgeProperties::bit(static_cast<unsigned>(k % 3)));
    }
    if (c.chains.empty()) {
        c.add_chain({0, 1});
    }
    return c;
}

TEST_CASE("I6: build() returns the brute force's status on random candidates",
          "[16b-0][builder][random]") {
    struct Setting {
        double spacing;
        std::int64_t origin;
    };
    // A non-dyadic spacing, a dyadic one, and grids pushed next to kMaxGridIndex
    // where a world coordinate's ulp is a large part of a cell.
    const Setting s = GENERATE(Setting{0.1, 0}, Setting{0.125, -40}, Setting{1.0, 0},
                               Setting{0.1, (std::int64_t{1} << 50) + 12345},
                               Setting{1.0, terrain::kMaxGridIndex - 64},
                               Setting{0.3, -(terrain::kMaxGridIndex - 64)});
    const Layout layout = GENERATE(Layout::Free, Layout::AxisLines, Layout::Diagonals);

    std::size_t ok = 0;
    std::size_t refused = 0;
    for (std::uint64_t seed = 1; seed <= 300; ++seed) {
        const Candidate c = random_candidate(seed, s.spacing, s.origin, layout);
        const NodeStatus want = c.oracle();
        INFO("seed " << seed << ", spacing " << s.spacing << ", origin " << s.origin);
        REQUIRE(c.built() == want);
        (want == NodeStatus::Ok ? ok : refused) += 1;
    }
    // The probe must be able to fail both ways: a family that is all Ok cannot
    // catch a sweep that misses pairs, and one that is all refused cannot catch
    // a sweep that invents them.
    INFO(ok << " Ok, " << refused << " refused");
    REQUIRE(ok >= 15);
    REQUIRE(refused >= 15);
}

// Real noder output is the "nearly noded" case: every pair is legal until one
// thing is added, and then exactly the pairs that touch the addition matter.
[[nodiscard]] Candidate from_noded(const terrain::NodedPslg& p) {
    Candidate c;
    c.grid = p.grid();
    c.vertices.assign(p.vertices().begin(), p.vertices().end());
    c.chains.assign(p.chains().begin(), p.chains().end());
    c.chain_indices.assign(p.chain_indices().begin(), p.chain_indices().end());
    c.edge_properties.assign(p.edge_properties().begin(), p.edge_properties().end());
    c.node_of_input_vertex.assign(p.node_of_input_vertex().begin(),
                                  p.node_of_input_vertex().end());
    return c;
}

TEST_CASE("I6: build() returns the brute force's status on noder output with one edge or "
          "one node added",
          "[16b-0][builder][random]") {
    std::size_t ok = 0;
    std::size_t refused = 0;
    for (std::uint64_t seed = 1; seed <= 40; ++seed) {
        terrain::testing::PolylineOptions opts;
        opts.breakline_count = 6;
        opts.reach = seed % 2 ? 100.0 : 15.0;
        const terrain::Pslg pslg =
            terrain::testing::to_pslg(terrain::testing::random_polylines(seed, opts));
        const auto noded = terrain::noding::node<DefaultKernel>(pslg, NodeOptions{0.1});
        REQUIRE(noded.ok());
        const Candidate base = from_noded(*noded.pslg);
        REQUIRE(base.oracle() == NodeStatus::Ok);
        REQUIRE(base.built() == NodeStatus::Ok);

        std::mt19937_64 rng{seed * 7919};
        for (int trial = 0; trial < 12; ++trial) {
            Candidate c = base;
            const auto n = static_cast<std::uint32_t>(c.vertices.size());
            if (trial % 3 == 0) {
                // An edge between two existing nodes.
                const auto a = static_cast<std::uint32_t>(rng() % n);
                const auto b = static_cast<std::uint32_t>(rng() % n);
                if (a == b) {
                    continue;
                }
                c.add_chain({a, b});
            } else {
                // A lone node at, or one cell beside, a random edge's midpoint:
                // grazes, near misses and T-junctions for 14(b).
                const std::size_t e = rng() % c.edge_properties.size();
                std::size_t chain = 0;
                while (chain + 1 < c.chains.size() &&
                       e >= base_edge_offset(c, chain + 1)) {
                    ++chain;
                }
                const std::size_t k = e - base_edge_offset(c, chain);
                const Chain& ch = c.chains[chain];
                const std::uint32_t u = c.chain_indices[ch.begin + k];
                const std::uint32_t v = c.chain_indices[ch.begin + (k + 1) % ch.count];
                const GridPoint gu = c.grid.snap(c.vertices[u]);
                const GridPoint gv = c.grid.snap(c.vertices[v]);
                const GridPoint mid{(gu.ix + gv.ix) / 2 + static_cast<std::int64_t>(rng() % 3) - 1,
                                    (gu.iy + gv.iy) / 2 + static_cast<std::int64_t>(rng() % 3) - 1};
                const auto before = c.vertices.size();
                (void)c.vertex(mid);
                if (c.vertices.size() == before) {
                    continue;
                }
            }
            const NodeStatus want = c.oracle();
            INFO("seed " << seed << ", trial " << trial);
            REQUIRE(c.built() == want);
            (want == NodeStatus::Ok ? ok : refused) += 1;
        }
    }
    INFO(ok << " Ok, " << refused << " refused");
    REQUIRE(ok >= 20);
    REQUIRE(refused >= 20);
}

// --- Named adversarial candidates ------------------------------------------
//
// Built on grid indices. Each asserts the oracle's verdict AND a stated one.

[[nodiscard]] Candidate on_grid(double spacing = 0.1) {
    Candidate c;
    c.grid = SnapGrid{spacing};
    return c;
}

TEST_CASE("I6 fixtures: duplicate edges in both orientations are legal", "[16b-0][builder]") {
    Candidate c = on_grid();
    const auto a = c.vertex({0, 0});
    const auto b = c.vertex({10, 3});
    const auto d = c.vertex({20, 0});
    c.add_chain({a, b, d});
    c.add_chain({d, b, a}, ChainRole::Breakline, EdgeProperties::bit(1));
    c.add_chain({a, b}, ChainRole::Breakline, EdgeProperties::bit(2));
    require_verdict(c, NodeStatus::Ok);
}

TEST_CASE("I6 fixtures: collinear partial overlaps are refused", "[16b-0][builder]") {
    SECTION("vertical, the two boxes sharing one x (zero width)") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({5, 0}), c.vertex({5, 4})});
        c.add_chain({c.vertex({5, 2}), c.vertex({5, 6})});
        require_verdict(c, NodeStatus::NotConverged);
    }
    SECTION("horizontal") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({0, 5}), c.vertex({4, 5})});
        c.add_chain({c.vertex({2, 5}), c.vertex({6, 5})});
        require_verdict(c, NodeStatus::NotConverged);
    }
    SECTION("diagonal, one containing the other") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({0, 0}), c.vertex({9, 9})});
        c.add_chain({c.vertex({3, 3}), c.vertex({6, 6})});
        require_verdict(c, NodeStatus::NotConverged);
    }
}

TEST_CASE("I6 fixtures: vertical edges stacked end to end at one x are legal",
          "[16b-0][builder]") {
    Candidate c = on_grid();
    std::vector<std::uint32_t> ids;
    for (std::int64_t k = 0; k <= 12; k += 3) {
        ids.push_back(c.vertex({7, k}));
    }
    for (std::size_t k = 0; k + 1 < ids.size(); ++k) {
        c.add_chain({ids[k], ids[k + 1]});
    }
    // And three more zero-width columns beside it, touching nothing.
    for (std::int64_t x = 9; x <= 11; ++x) {
        c.add_chain({c.vertex({x, 0}), c.vertex({x, 12})});
    }
    require_verdict(c, NodeStatus::Ok);
}

TEST_CASE("I6 fixtures: a T-junction is refused, a node two cells off is not",
          "[16b-0][builder]") {
    SECTION("an endpoint in the other edge's interior") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({0, 0}), c.vertex({20, 0})});
        c.add_chain({c.vertex({10, 0}), c.vertex({10, 8})});
        require_verdict(c, NodeStatus::NotConverged);
    }
    SECTION("a lone node on a vertical edge") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({4, 0}), c.vertex({4, 20})});
        (void)c.vertex({4, 11});
        require_verdict(c, NodeStatus::NotConverged);
    }
    SECTION("a lone node two cells off") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({0, 0}), c.vertex({20, 0})});
        (void)c.vertex({10, 2});
        require_verdict(c, NodeStatus::Ok);
    }
}

TEST_CASE("I6 fixtures: an edge through a cell corner meets that cell", "[16b-0][builder]") {
    // (0,0)-(1,1) in grid units passes through (1/2, 1/2), the shared corner of
    // the cells of (0,0), (1,0), (0,1) and (1,1): closed cells, so (1,0) and
    // (0,1) are met.
    Candidate c = on_grid(0.125);
    c.add_chain({c.vertex({0, 0}), c.vertex({1, 1})});
    (void)c.vertex({1, 0});
    require_verdict(c, NodeStatus::NotConverged);
}

TEST_CASE("I6 fixtures: a node whose cell box meets an edge's box only at a corner",
          "[16b-0][builder]") {
    // The edge (0,0)-(4,2) has box [0,4]x[0,2]; the cell of (5,3) is
    // [4.5,5.5]x[2.5,3.5], and padded by one spacing its box overlaps the
    // edge's box in the corner [3.5,4]x[1.5,2] only. The segment does not reach
    // the cell: a candidate pair for the sweep, Ok for the predicate.
    Candidate c = on_grid(0.125);
    c.add_chain({c.vertex({0, 0}), c.vertex({4, 2})});
    (void)c.vertex({5, 3});
    require_verdict(c, NodeStatus::Ok);
}

TEST_CASE("I6 fixtures: one long edge spanning every other box", "[16b-0][builder]") {
    SECTION("with a node on it far along the sweep: refused") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({0, 0}), c.vertex({1000, 1000})});
        for (std::int64_t k = 0; k < 40; ++k) {
            c.add_chain({c.vertex({25 * k + 3, 25 * k}), c.vertex({25 * k + 5, 25 * k})});
        }
        require_verdict(c, NodeStatus::Ok);
        (void)c.vertex({997, 997});
        require_verdict(c, NodeStatus::NotConverged);
    }
    SECTION("crossed by the last edge in x order: refused") {
        Candidate c = on_grid();
        c.add_chain({c.vertex({0, 0}), c.vertex({1000, 0})});
        for (std::int64_t k = 1; k < 40; ++k) {
            c.add_chain({c.vertex({25 * k, 3}), c.vertex({25 * k + 2, 9})});
        }
        require_verdict(c, NodeStatus::Ok);
        c.add_chain({c.vertex({999, -3}), c.vertex({999, 3})});
        require_verdict(c, NodeStatus::NotConverged);
    }
}

TEST_CASE("I6 fixtures: a crossing next to kMaxGridIndex", "[16b-0][builder]") {
    const std::int64_t o = terrain::kMaxGridIndex - 32;
    Candidate c = on_grid(1.0);
    c.add_chain({c.vertex({o, o}), c.vertex({o + 20, o + 20})});
    require_verdict(c, NodeStatus::Ok);
    c.add_chain({c.vertex({o, o + 20}), c.vertex({o + 20, o})});
    require_verdict(c, NodeStatus::NotConverged);

    Candidate t = on_grid(1.0);
    t.add_chain({t.vertex({-o, -o}), t.vertex({-o + 30, -o})});
    (void)t.vertex({-o + 17, -o});
    require_verdict(t, NodeStatus::NotConverged);
}

// ===========================================================================
// C. node<K> unchanged
// ===========================================================================

class Fnv {
public:
    void add(std::uint64_t v) noexcept {
        for (int k = 0; k < 8; ++k) {
            h_ ^= (v >> (8 * k)) & 0xffu;
            h_ *= 1099511628211ull;
        }
    }
    [[nodiscard]] std::uint64_t value() const noexcept { return h_; }

private:
    std::uint64_t h_{14695981039346656037ull};
};

// Integers only: the status, every node's GridPoint, the chains, the flat
// indices, the masks and node_of_input_vertex. Doubles are left out so FMA
// contraction and libm cannot move the value while the noder's decisions stand.
[[nodiscard]] std::uint64_t noded_digest(const terrain::noding::NodeOutcome& out) {
    Fnv f;
    f.add(static_cast<std::uint64_t>(out.status));
    if (!out.ok()) {
        return f.value();
    }
    const terrain::NodedPslg& p = *out.pslg;
    f.add(p.vertices().size());
    for (const Point2& v : p.vertices()) {
        const GridPoint g = p.grid().snap(v);
        f.add(static_cast<std::uint64_t>(g.ix));
        f.add(static_cast<std::uint64_t>(g.iy));
    }
    for (const Chain& ch : p.chains()) {
        f.add(ch.begin);
        f.add(ch.count);
        f.add(static_cast<std::uint64_t>(ch.role));
    }
    for (const std::uint32_t i : p.chain_indices()) {
        f.add(i);
    }
    for (const EdgeProperties e : p.edge_properties()) {
        f.add(e.bits());
    }
    for (const std::uint32_t i : p.node_of_input_vertex()) {
        f.add(i);
    }
    return f.value();
}

struct Lines {
    std::vector<Point2> vertices;
    std::vector<std::vector<std::uint32_t>> runs;
    std::vector<ChainRole> roles;

    std::uint32_t at(double x, double y) {
        vertices.push_back(Point2{x, y});
        return static_cast<std::uint32_t>(vertices.size() - 1);
    }
    void line(std::vector<std::pair<double, double>> pts,
              ChainRole role = ChainRole::Breakline) {
        std::vector<std::uint32_t> run;
        for (const auto& [x, y] : pts) {
            run.push_back(at(x, y));
        }
        runs.push_back(std::move(run));
        roles.push_back(role);
    }
    [[nodiscard]] terrain::Pslg pslg() const {
        terrain::PslgBuilder b{vertices};
        for (std::size_t c = 0; c < runs.size(); ++c) {
            b.add_chain(runs[c], roles[c], EdgeProperties::bit(static_cast<unsigned>(c % 4)));
        }
        auto r = std::move(b).build<DefaultKernel>();
        REQUIRE(r.ok());
        return std::move(*r.pslg);
    }
};

[[nodiscard]] Lines square(double m) {
    Lines l;
    l.line({{0, 0}, {m, 0}, {m, m}, {0, m}}, ChainRole::Outer);
    return l;
}

[[nodiscard]] Lines fixture(int which) {
    Lines l = square(40);
    switch (which) {
        case 0:  // a hash: two by two crossing lines
            l.line({{10, 2}, {10, 38}});
            l.line({{30, 2}, {30, 38}});
            l.line({{2, 10}, {38, 10}});
            l.line({{2, 30}, {38, 30}});
            break;
        case 1:  // a star of six lines through one point
            for (int k = 0; k < 6; ++k) {
                const double dx = std::array{1, 2, 3, 1, -1, 5}[k];
                const double dy = std::array{0, 1, 3, 3, 2, -2}[k];
                l.line({{20 - 3 * dx, 20 - 3 * dy}, {20 + 3 * dx, 20 + 3 * dy}});
            }
            break;
        case 2:  // collinear partial overlaps on one line, and on the outer ring
            l.line({{0, 5}, {10, 5}});
            l.line({{5, 5}, {15, 5}});
            l.line({{7, 5}, {20, 5}});
            l.line({{3, 0}, {17, 0}});
            break;
        case 3:  // a road along a river, reversed, and a T-junction onto it
            l.line({{2, 20}, {20, 21}, {38, 20}});
            l.line({{38, 20}, {20, 21}, {2, 20}});
            l.line({{20, 21}, {20, 35}});
            l.line({{11, 20.5}, {11, 3}});
            break;
        case 4:  // many verticals at one x, overlapping and crossed
            for (int k = 0; k < 6; ++k) {
                l.line({{17, 2.0 + 5 * k}, {17, 9.0 + 5 * k}});
            }
            l.line({{3, 3}, {37, 37}});
            break;
        default: {  // 40 random polylines from raw mt19937_64 words
            l = square(1000);
            std::mt19937_64 rng{16};
            for (int k = 0; k < 40; ++k) {
                std::vector<std::pair<double, double>> pts;
                for (int v = 0; v < 3; ++v) {
                    // Two draws as arguments are unsequenced (clang: x first, GCC: y first); fixed as clang's.
                    const auto x = static_cast<double>(1 + rng() % 998);
                    const auto y = static_cast<double>(1 + rng() % 998);
                    pts.emplace_back(x, y);
                }
                l.line(pts);
            }
        }
    }
    return l;
}

TEST_CASE("node<K>'s output is today's, digest for digest", "[16b-0][noder]") {
    // Recorded from 6b4fcb9's node<K> (the brute-force verifier) with this
    // digest, before any production change on the branch. Dyadic spacing and
    // integer or half-integer input: nothing here depends on rounding.
    constexpr std::uint64_t kRecorded[] = {
        0x80d6a6a851692c2eull, 0xb71d5b267cb92eacull, 0x9632a57ca8ae58bdull,
        0x1facb9f200887f70ull, 0x055c808b84bf1a80ull, 0xe9be0099ef219da5ull,
    };
    const int which = GENERATE(range(0, 6));
    const auto out = terrain::noding::node<DefaultKernel>(fixture(which).pslg(),
                                                           NodeOptions{0.125});
    INFO("fixture " << which << ", status " << static_cast<int>(out.status) << ", digest 0x"
                    << std::hex << noded_digest(out));
    REQUIRE(noded_digest(out) == kRecorded[which]);
}

TEST_CASE("node<K>'s Ok output passes the brute force", "[16b-0][noder][random]") {
    for (std::uint64_t seed = 1; seed <= 30; ++seed) {
        terrain::testing::PolylineOptions opts;
        opts.breakline_count = 10;
        opts.reach = seed % 3 == 0 ? 100.0 : 20.0;
        const auto out = terrain::noding::node<DefaultKernel>(
            terrain::testing::to_pslg(terrain::testing::random_polylines(seed, opts)),
            NodeOptions{seed % 2 ? 0.1 : 0.125});
        INFO("seed " << seed);
        REQUIRE(out.ok());
        REQUIRE(from_noded(*out.pslg).oracle() == NodeStatus::Ok);
    }
}

}  // namespace

