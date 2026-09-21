// Unit tests for terrain/noding/node.hpp -- the driver, on named fixtures.
//
// node.hpp's SECOND test file, and the one place Gate B's one-file-per-header
// rule is relaxed in writing (docs/increments/05b-noder-driver.md, "What is
// worth testing"): the property file and the fixture file split along a
// different axis than the one Gate B protects.
//
// The mutation round is prop_noding_no_crossings.cpp's, but four of the
// increment's mutants can only be killed here, because each needs a fixture
// whose parts are known by name rather than generated:
//
//   mutant  3  the split pass driven by classify instead of segment_meets_cell
//   mutant  4  segment_meets_cell queried only within a segment's own chain
//   mutant  7  the property union replaced by "take the first contributor"
//   mutant 11  arc order taken by node id instead of by dot product
//
// THE TWO FIXTURES THAT MATTER MOST ARE THE FIRST TWO, and not the degenerate
// ones. The user's case is a road crossing a river, or entering a forest;
// `rasputin draw not-noded` currently draws zero triangles for exactly the
// crossing fixture below, and the road-along-a-river fixture is the one input
// the whole property-set apparatus exists for.
//
// Every coordinate is either derived through grid.world() or is one of the
// gallery fixture's own numbers. The arithmetic of each degenerate fixture was
// checked against the shipped 5a predicates before it was written down, and each
// case restates the check it passed as its own REQUIRE, so the fixture cannot
// quietly stop being the case it claims to be.

#include <catch2/catch_test_macros.hpp>

#include "../property/noding_generators.h"

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/noding/node.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using terrain::ChainRole;
using terrain::EdgeProperties;
using terrain::GridPoint;
using terrain::NodedPslg;
using terrain::Point2;
using terrain::Pslg;
using terrain::PslgBuilder;
using terrain::PslgBuildResult;
using terrain::Segment2;
using terrain::SnapGrid;
using terrain::noding::NodeOptions;
using terrain::noding::NodeOutcome;
using terrain::noding::NodeStatus;
using terrain::noding::SegmentRelation;
using terrain::noding::classify;
using terrain::noding::node;
using terrain::noding::segment_meets_cell;
using terrain::on_segment;
using terrain::pred::DefaultKernel;

namespace {

const EdgeProperties kRoad = EdgeProperties::bit(0);
const EdgeProperties kRiver = EdgeProperties::bit(1);

[[nodiscard]] Pslg built(PslgBuilder builder) {
    PslgBuildResult result = std::move(builder).build<DefaultKernel>();
    REQUIRE(result.ok());
    return std::move(*result.pslg);
}

// A square ring in world units, counterclockwise, as the outer chain.
[[nodiscard]] std::vector<Point2> square(double lo_x, double lo_y, double hi_x, double hi_y) {
    return {Point2{lo_x, lo_y}, Point2{hi_x, lo_y}, Point2{hi_x, hi_y}, Point2{lo_x, hi_y}};
}

// The node id of a grid point in the output, or the vertex count if absent.
[[nodiscard]] std::uint32_t id_of(const NodedPslg& p, const GridPoint& g) {
    const Point2 want = p.grid().world(g);
    for (std::uint32_t i = 0; i < p.vertices().size(); ++i) {
        if (p.vertices()[i].x == want.x && p.vertices()[i].y == want.y) {
            return i;
        }
    }
    return static_cast<std::uint32_t>(p.vertices().size());
}

[[nodiscard]] std::size_t count_of(const NodedPslg& p, const Point2& want) {
    std::size_t n = 0;
    for (const Point2& v : p.vertices()) {
        if (v.x == want.x && v.y == want.y) {
            ++n;
        }
    }
    return n;
}

[[nodiscard]] std::vector<std::uint32_t> run_of(const NodedPslg& p, std::size_t c) {
    const auto idx = p.indices_of(c);
    return std::vector<std::uint32_t>{idx.begin(), idx.end()};
}

}  // namespace

// ---------------------------------------------------------------------------
// NodeOptions. spacing is first and has no default; max_rounds does.
// ---------------------------------------------------------------------------

TEST_CASE("NodeOptions is spelled NodeOptions{spacing} and caps rounds at four",
          "[noding][node]") {
    const NodeOptions o{0.1};
    REQUIRE(o.spacing == 0.1);

    // Four rather than one because a snap-induced cascade is real, and four
    // rather than sixteen because the audit measured it to be rarer than the
    // theory allows. At any value the outcome is either the same mesh or
    // NotConverged, never a different mesh -- which is why this one has a
    // default and spacing does not.
    REQUIRE(o.max_rounds == 4);

    STATIC_REQUIRE(std::is_aggregate_v<NodeOptions>);
    STATIC_REQUIRE(std::is_same_v<decltype(NodeOptions::spacing), double>);
}

TEST_CASE("a default-constructed NodeOutcome from the driver's own type is NotRun",
          "[noding][node]") {
    const NodeOutcome out;
    REQUIRE(out.status == NodeStatus::NotRun);
    REQUIRE_FALSE(out.ok());
}

// ---------------------------------------------------------------------------
// THE CROSSING. The increment's reason for existing, and the exact input
// `rasputin draw not-noded` currently renders with zero triangles.
// ---------------------------------------------------------------------------

namespace {

constexpr double kOriginX = 430000.0;
constexpr double kOriginY = 6900000.0;

[[nodiscard]] Point2 utm(double x, double y) { return Point2{kOriginX + x, kOriginY + y}; }

// The `not-noded` gallery fixture: a 700 m square with two breaklines crossing
// at local (350, 350), a point neither chain names.
[[nodiscard]] Pslg crossing_pslg() {
    PslgBuilder b;
    b.add_chain(square(kOriginX, kOriginY, kOriginX + 700.0, kOriginY + 700.0), ChainRole::Outer);
    const std::vector<Point2> road{utm(100.0, 100.0), utm(600.0, 600.0)};
    const std::vector<Point2> river{utm(100.0, 600.0), utm(600.0, 100.0)};
    b.add_chain(road, ChainRole::Breakline, kRoad);
    b.add_chain(river, ChainRole::Breakline, kRiver);
    return built(std::move(b));
}

}  // namespace

TEST_CASE("two crossing breaklines are split at the crossing, which becomes one node",
          "[noding][node][crossing]") {
    const Pslg in = crossing_pslg();
    const SnapGrid grid{0.1};

    // The fixture is what it claims to be: the two chains genuinely cross, and
    // the crossing point is not a vertex of either.
    REQUIRE(classify<DefaultKernel>(in.edge(1, 0), in.edge(2, 0)) == SegmentRelation::Crossing);

    const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{0.1});
    REQUIRE(out.status == NodeStatus::Ok);
    REQUIRE(out.ok());
    const NodedPslg& p = *out.pslg;

    // The acceptance criterion 05b-noder-driver.md states for 5b and 5c: the
    // vertex array contains the snapped image of (350, 350) + ORIGIN, and
    // contains it ONCE.
    const Point2 meeting = p.grid().world(grid.snap(utm(350.0, 350.0)));
    REQUIRE(count_of(p, meeting) == 1);

    // MUTANT 4: segment_meets_cell queried only for the nodes of the segment's
    // own chain. The crossing node belongs to neither chain's input vertices --
    // it is manufactured in step 2 -- and both chains must be split at it.
    const std::uint32_t junction = id_of(p, grid.snap(utm(350.0, 350.0)));
    REQUIRE(junction < p.vertices().size());

    const std::vector<std::uint32_t> road = run_of(p, 1);
    const std::vector<std::uint32_t> river = run_of(p, 2);
    REQUIRE(road.size() == 3);
    REQUIRE(river.size() == 3);
    REQUIRE(road[1] == junction);
    REQUIRE(river[1] == junction);

    // Guarantee 16: three chains in, three chains out, roles and property sets
    // index-for-index. The outer ring is untouched and is still four nodes.
    REQUIRE(p.chains().size() == 3);
    REQUIRE(p.chains()[0].role == ChainRole::Outer);
    REQUIRE(p.chains()[1].properties == kRoad);
    REQUIRE(p.chains()[2].properties == kRiver);
    REQUIRE(p.edge_count(0) == 4);

    // Each half of each breakline inherits its own chain's set and nothing else:
    // the two chains cross, they do not run together, so no edge is both.
    for (std::size_t k = 0; k < p.edge_count(1); ++k) {
        REQUIRE(p.edge_properties()[p.edge_base(1) + k] == kRoad);
    }
    for (std::size_t k = 0; k < p.edge_count(2); ++k) {
        REQUIRE(p.edge_properties()[p.edge_base(2) + k] == kRiver);
    }
}

TEST_CASE("the crossing fixture is what increment 4 refused", "[noding][node][crossing]") {
    // Recorded as a fixture-level fact rather than as a CDT call, because
    // triangulate does not take a NodedPslg until 5c. What increment 4 said
    // about this input is that its constrained edges are intersecting; what 5b
    // says is that after noding they are not.
    const Pslg in = crossing_pslg();
    const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{0.1});
    REQUIRE(out.ok());

    const NodedPslg& p = *out.pslg;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        for (std::size_t k = 0; k < p.edge_count(c); ++k) {
            for (std::size_t d = c; d < p.chains().size(); ++d) {
                for (std::size_t l = (d == c ? k + 1 : 0); l < p.edge_count(d); ++l) {
                    REQUIRE(classify<DefaultKernel>(p.edge(c, k), p.edge(d, l)) !=
                            SegmentRelation::Crossing);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// THE T-JUNCTION, in the relation the producer uses.
// ---------------------------------------------------------------------------

TEST_CASE("a T-junction is found by proximity, on input that classifies Disjoint",
          "[noding][node][t_junction]") {
    // MUTANT 3: the split pass driven by classify<K> rather than by
    // segment_meets_cell<K> -- the provisional design 5a's audit refuted. This
    // fixture is its kill and the three REQUIREs below are why: the pair is
    // Disjoint, so a classify-driven pass never looks at it; the intruding node's
    // closed cell nevertheless meets the host; and the node is not on the host by
    // any exact reading, so an on_segment-driven pass misses it too.
    const SnapGrid grid{0.1};
    const Segment2 host{grid.world(GridPoint{10, 10}), grid.world(GridPoint{14, 12})};
    const Segment2 spur{grid.world(GridPoint{11, 10}), grid.world(GridPoint{11, 5})};

    REQUIRE(classify<DefaultKernel>(host, spur) == SegmentRelation::Disjoint);
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, host, GridPoint{11, 10}));
    REQUIRE_FALSE(on_segment<DefaultKernel>(host, grid.world(GridPoint{11, 10})));

    PslgBuilder b;
    b.add_chain(square(-10.0, -10.0, 10.0, 10.0), ChainRole::Outer);
    b.add_chain(std::vector<Point2>{host.a, host.b}, ChainRole::Breakline, kRoad);
    b.add_chain(std::vector<Point2>{spur.a, spur.b}, ChainRole::Breakline, kRiver);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{0.1});
    REQUIRE(out.status == NodeStatus::Ok);
    const NodedPslg& p = *out.pslg;

    const std::vector<std::uint32_t> split_host = run_of(p, 1);
    REQUIRE(split_host.size() == 3);
    REQUIRE(split_host[1] == id_of(p, GridPoint{11, 10}));

    // The spur is not split: nothing meets its cells but its own endpoints.
    REQUIRE(run_of(p, 2).size() == 2);

    // The deformation is real and is bounded: the split node is half a cell off
    // the host's own line. Snap rounding moves geometry, and asserting it does
    // not is how increment 5's first oracle came to be red on correct output.
    const Point2 moved = p.grid().world(GridPoint{11, 10});
    REQUIRE_FALSE(on_segment<DefaultKernel>(host, moved));
}

// ---------------------------------------------------------------------------
// THE LATTICE-CORNER TIE. And it does not do what the design says it does.
// ---------------------------------------------------------------------------

TEST_CASE("the unit diagonal grazes both flanking cells and the tie is exact",
          "[noding][node][lattice]") {
    // The arithmetic half, written only in terms of shipped 5a predicates, so it
    // stands whatever node<K> turns out to do. segment_meets_cell uses CLOSED
    // cells and answers true on a grazed corner (intersect.hpp:203-213; 5a
    // mutant 13 pins it), so guarantee 14(b) REQUIRES the diagonal to be split
    // at both flanking nodes, and their arc keys are exactly equal.
    const SnapGrid grid{0.1};
    const Segment2 diagonal{grid.world(GridPoint{0, 0}), grid.world(GridPoint{1, 1})};

    REQUIRE(segment_meets_cell<DefaultKernel>(grid, diagonal, GridPoint{1, 0}));
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, diagonal, GridPoint{0, 1}));

    const auto arc = [&](const GridPoint& g) {
        const Point2 d{diagonal.b.x - diagonal.a.x, diagonal.b.y - diagonal.a.y};
        const Point2 w = grid.world(g);
        return (w.x - diagonal.a.x) * d.x + (w.y - diagonal.a.y) * d.y;
    };
    REQUIRE(arc(GridPoint{1, 0}) == arc(GridPoint{0, 1}));

    // The tie is broken by GridPoint's defaulted operator<=>, so (0,1) precedes
    // (1,0) and the output chain is (0,0) -> (0,1) -> (1,0) -> (1,1).
    REQUIRE(GridPoint{0, 1} < GridPoint{1, 0});

    // AND THAT CHAIN VIOLATES GUARANTEE 14(b). The middle edge (0,1)-(1,0) is
    // the anti-diagonal through the same lattice corner, so it grazes the closed
    // cells of (0,0) and (1,1) -- neither of which is an endpoint of it.
    //
    // An earlier revision of 05b-noder-driver.md, "The unit diagonal through a
    // lattice corner", asserted the opposite -- that guarantee 14 still holds
    // here and @tester should check it by hand rather than assume. Checked by
    // hand, and it does not; the design now records that correction and keeps
    // the refuted bullet as the shape of the error. What it checked was 14(a),
    // over edge PAIRS, and 14(a) does hold: the two non-consecutive edges are
    // parallel and disjoint. It never asks 14(b) -- a (node, edge) claim -- of
    // the chain's own endpoints against its own middle edge. This case pins that
    // arithmetic, and the case below is its consequence.
    const Segment2 middle{grid.world(GridPoint{0, 1}), grid.world(GridPoint{1, 0})};
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, middle, GridPoint{0, 0}));
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, middle, GridPoint{1, 1}));
}

TEST_CASE("a constraint through a lattice corner flanked by nodes does not converge",
          "[noding][node][lattice]") {
    // The consequence of the case above. Round 1 routes the diagonal through both
    // grazed cells; the resulting middle edge grazes the cells of the diagonal's
    // own endpoints, so round 2 must split it at them; that reintroduces both
    // endpoints into the middle of the chain and the sequence does not settle.
    //
    // NotConverged is the right answer and it has an actionable lever behind it
    // -- a finer or coarser spacing moves the corner off the lattice. What is
    // NOT right is Ok, and a driver that reports Ok here has either skipped the
    // verification pass or is running a verifier that agrees with it.
    //
    // DO NOT DELETE THIS CASE AS "A TEST OF A BUG WE ARE GOING TO FIX". Its
    // expected status is increment 5d's to change, but the fixture is mutant 5's
    // killer -- the loop condition spelled "until a round produces no new node"
    // instead of "until the verification passes". Under that spelling the tie
    // converges to Ok at the end of round 1, because (1,0) and (0,1) are already
    // nodes from the flanking chains, so round 1 splits an edge at two EXISTING
    // nodes and manufactures none; under the correct spelling it runs to the cap.
    // 05b-noder-driver.md, "What the lattice-corner fixture is for", records this
    // as the only fixture in the increment on which the two spellings differ in
    // STATUS rather than in edge count, so deleting it costs mutant 5 the only
    // killer that fails loudly.
    const SnapGrid grid{0.1};
    PslgBuilder b;
    b.add_chain(square(-10.0, -10.0, 10.0, 10.0), ChainRole::Outer);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{0, 0}), grid.world(GridPoint{1, 1})},
                ChainRole::Breakline, kRoad);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{1, 0}), grid.world(GridPoint{5, -4})},
                ChainRole::Breakline);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{0, 1}), grid.world(GridPoint{-4, 5})},
                ChainRole::Breakline);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{0.1});
    REQUIRE_FALSE(out.ok());
    REQUIRE_FALSE(out.pslg.has_value());
    REQUIRE(out.status == NodeStatus::NotConverged);
    REQUIRE_FALSE(out.message.empty());
}

// ---------------------------------------------------------------------------
// THE ROAD ALONG A RIVER. Mutant 7, and the one input the property-set
// apparatus exists for.
// ---------------------------------------------------------------------------

TEST_CASE("a road running along a river gives every shared edge both properties",
          "[noding][node][road_river]") {
    // MUTANT 7: the property union replaced by "take the first contributor" or
    // by "take the chain's own set". Guarantee 15's oracle is one-directional and
    // cannot see a spuriously ADDED property, so it cannot see a REPLACED one
    // either when the replacement happens to be right; this fixture is the thing
    // standing between a lost property and a silent wrong answer.
    //
    // The two chains are four millimetres apart, which at a decimetre grid is
    // less than half a cell, so both snap onto the same two nodes and carry the
    // SAME node-id pair. Guarantee 14(a) as amended permits exactly that, and
    // the pre-amendment clause refused it.
    const SnapGrid grid{0.1};
    REQUIRE(grid.snap(Point2{0.004, 0.004}) == GridPoint{0, 0});
    REQUIRE(grid.snap(Point2{10.004, 0.004}) == GridPoint{100, 0});

    PslgBuilder b;
    b.add_chain(square(-10.0, -10.0, 20.0, 20.0), ChainRole::Outer);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{0, 0}), grid.world(GridPoint{100, 0})},
                ChainRole::Breakline, kRiver);
    b.add_chain(std::vector<Point2>{Point2{0.004, 0.004}, Point2{10.004, 0.004}},
                ChainRole::Breakline, kRoad);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{0.1});
    REQUIRE(out.status == NodeStatus::Ok);
    const NodedPslg& p = *out.pslg;

    const std::vector<std::uint32_t> river = run_of(p, 1);
    const std::vector<std::uint32_t> road = run_of(p, 2);
    REQUIRE(river.size() == 2);
    REQUIRE(road.size() == 2);
    REQUIRE(std::min(river[0], river[1]) == std::min(road[0], road[1]));
    REQUIRE(std::max(river[0], river[1]) == std::max(road[0], road[1]));

    // THE ASSERTION MUTANT 7 DIES ON, and it is made on the ROAD's flat
    // position, not on the river's: "take the chain's own set" gets the river's
    // entry right by accident and the road's wrong.
    const EdgeProperties on_road = p.edge_properties()[p.edge_base(2) + 0];
    const EdgeProperties on_river = p.edge_properties()[p.edge_base(1) + 0];
    REQUIRE(on_road.contains(kRoad));
    REQUIRE(on_road.contains(kRiver));
    REQUIRE(on_river.contains(kRoad));
    REQUIRE(on_river.contains(kRiver));

    // Union, not replacement: the two entries agree, which "take the first
    // contributor" also satisfies -- so the contains() pair above is what does
    // the work and this is the corroboration.
    REQUIRE(on_road == on_river);

    // Guarantee 16: the CHAINS are not merged. Two chains in, two out, each
    // keeping its own set. What they come to share is node ids, and therefore
    // per-edge entries, not identity.
    REQUIRE(p.chains().size() == 3);
    REQUIRE(p.chains()[1].properties == kRiver);
    REQUIRE(p.chains()[2].properties == kRoad);
}

// ---------------------------------------------------------------------------
// ARC ORDER.
// ---------------------------------------------------------------------------

TEST_CASE("split points are ordered along the segment, not by node id",
          "[noding][node][arc_order]") {
    // MUTANT 11. The segment runs from high indices to low, so its arc order is
    // strictly DECREASING in node id -- node ids are NodeSet's lexicographic
    // order of grid points. A driver that sorted split points by id would emit
    // the chain reversed through its own middle, and the doubling-back edges
    // would overlap, which guarantee 14(a) forbids.
    //
    // (The design says "any segment running down-and-right"; measured, the
    // direction that actually inverts GridPoint's lexicographic order is
    // decreasing x. Recorded here rather than left to the next reader.)
    const SnapGrid grid{0.1};
    PslgBuilder b;
    b.add_chain(square(-10.0, -10.0, 20.0, 20.0), ChainRole::Outer);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{100, 100}), grid.world(GridPoint{0, 0})},
                ChainRole::Breakline, kRoad);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{70, 70}), grid.world(GridPoint{70, 40})},
                ChainRole::Breakline);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{30, 30}), grid.world(GridPoint{30, 0})},
                ChainRole::Breakline);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{0.1});
    REQUIRE(out.status == NodeStatus::Ok);
    const NodedPslg& p = *out.pslg;

    const std::vector<std::uint32_t> run = run_of(p, 1);
    REQUIRE(run.size() == 4);
    REQUIRE(run[0] == id_of(p, GridPoint{100, 100}));
    REQUIRE(run[1] == id_of(p, GridPoint{70, 70}));
    REQUIRE(run[2] == id_of(p, GridPoint{30, 30}));
    REQUIRE(run[3] == id_of(p, GridPoint{0, 0}));

    // And the ids really do run the other way, so the fixture is the case it
    // claims to be rather than one a node-id sort would also pass.
    REQUIRE(run[0] > run[1]);
    REQUIRE(run[1] > run[2]);
    REQUIRE(run[2] > run[3]);
}

// ---------------------------------------------------------------------------
// The failure statuses, one named case each. The table in
// 05b-noder-driver.md, "Degeneracy and failure policy", is the specification.
// ---------------------------------------------------------------------------

TEST_CASE("an invalid spacing is refused before anything is snapped", "[noding][node][status]") {
    const Pslg in = crossing_pslg();

    // Spelled over four values rather than one: is_valid_spacing is `> 0.0 &&
    // <= max()`, and the `!= 0.0` spelling it warns against accepts every
    // negative spacing, whose grid runs backwards.
    for (const double bad : {0.0, -0.0, -1.0, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
        const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{bad});
        REQUIRE(out.status == NodeStatus::InvalidSnapSpacing);
        REQUIRE_FALSE(out.ok());
        REQUIRE_FALSE(out.message.empty());
    }
}

TEST_CASE("a coordinate past kMaxGridIndex is a status, not an assertion",
          "[noding][node][status]") {
    // MUTANT 14: the can_snap admission pass dropped. Without it snap() aborts
    // in a Debug build and wraps an index in Release, and the noder produces a
    // mesh with a vertex somewhere else entirely.
    //
    // The offending vertex is UNREFERENCED, which is legal in a Pslg and is
    // exactly the case increment 3's stage 2 exists for: step 0 passes over
    // pslg.vertices(), not over the chains.
    PslgBuilder b;
    b.add_chain(square(-10.0, -10.0, 10.0, 10.0), ChainRole::Outer);
    b.append_vertices(std::vector<Point2>{Point2{1e300, 0.0}});

    const SnapGrid grid{0.1};
    REQUIRE_FALSE(grid.can_snap(Point2{1e300, 0.0}));

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{0.1});
    REQUIRE(out.status == NodeStatus::CoordinateOutOfRange);
    REQUIRE_FALSE(out.ok());
}

TEST_CASE("a hole smaller than one cell collapses and is refused", "[noding][node][status]") {
    // MUTANT 9: the RingCollapsed check dropped. All three vertices snap to the
    // same grid point, so the ring has one distinct node and is not a ring.
    // The lever is a FINER spacing, which is the opposite of what
    // CoordinateOutOfRange wants -- 05-noder.md risk 2 showing through the enum.
    const SnapGrid grid{1.0};
    for (const Point2& p : {Point2{5.0, 5.0}, Point2{5.0, 5.01}, Point2{5.01, 5.0}}) {
        REQUIRE(grid.snap(p) == GridPoint{5, 5});
    }

    PslgBuilder b;
    b.add_chain(square(-50.0, -50.0, 50.0, 50.0), ChainRole::Outer);
    // Clockwise, as a hole must be before snapping.
    b.add_chain(std::vector<Point2>{Point2{5.0, 5.0}, Point2{5.0, 5.01}, Point2{5.01, 5.0}},
                ChainRole::Hole);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{1.0});
    REQUIRE(out.status == NodeStatus::RingCollapsed);
    REQUIRE_FALSE(out.ok());
    REQUIRE_FALSE(out.message.empty());
}

TEST_CASE("a breakline smaller than one cell collapses to one node and is Ok",
          "[noding][node][status]") {
    // The open counterpart of the case above, and the outcome is the OPPOSITE:
    // Ok, by design (05b-noder-driver.md, "An open chain that collapses to one
    // node is Ok, and it is pinned"). The count < 3 refusal is is_closed(role)'s
    // alone. A ring of under three nodes has no interior and no winding, so
    // orientation<K> has nothing to answer; an open chain of one node is asked
    // no question it cannot answer, and an over-coarse spacing is the caller's
    // instruction rather than the noder's bug. The builder-side pin is
    // test_noding_noded_pslg_builder.cpp's pair of count-guard cases; this one
    // pins the end-to-end outcome, over arrays the noder really produced.
    //
    // THE CHAIN COUNT IS THE LOAD-BEARING ASSERTION. Guarantee 16 promises
    // index-for-index chain order with the input, which is what lets a
    // Python-side per-chain attribute array stay valid without a map, and it is
    // the mitigation for half of 05-noder.md risk 5. Dropping a collapsed chain
    // -- the obvious alternative to accepting it -- breaks exactly that, and a
    // silently dropped chain is the failure 16 names.
    const SnapGrid grid{1.0};
    const std::vector<Point2> breakline{Point2{5.0, 5.0}, Point2{5.0, 5.2}, Point2{5.1, 5.3}};
    for (const Point2& p : breakline) {
        REQUIRE(grid.snap(p) == GridPoint{5, 5});
    }

    PslgBuilder b;
    b.add_chain(square(-50.0, -50.0, 50.0, 50.0), ChainRole::Outer);
    b.add_chain(breakline, ChainRole::Breakline, kRiver);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{1.0});
    REQUIRE(out.status == NodeStatus::Ok);
    REQUIRE(out.ok());

    const NodedPslg& p = *out.pslg;
    REQUIRE(p.chains().size() == 2);
    REQUIRE(p.chains()[1].role == ChainRole::Breakline);
    REQUIRE(p.chains()[1].properties == kRiver);

    // One node, no edge, and the empty range every downstream loop handles.
    REQUIRE(run_of(p, 1).size() == 1);
    REQUIRE(run_of(p, 1)[0] == id_of(p, GridPoint{5, 5}));
    REQUIRE(p.edge_count(1) == 0);
    REQUIRE(p.edge_base(1) == p.edge_base(2));
}

TEST_CASE("a hole whose winding flips under snapping is refused, never reversed",
          "[noding][node][status]") {
    // THIS CASE ASSERTS NonSimpleRing, NOT RingDegenerateAfterSnap, AND THAT IS
    // NOW THE DOCUMENTED BEHAVIOUR RATHER THAN AN OBSERVATION AWAITING A RULING.
    // 05b-noder-driver.md, "`RingDegenerateAfterSnap` is a self-check, and risk
    // 3's mitigation moves", demotes that status to a self-check with a name,
    // exactly like MalformedOutput and CdtStatus::NotNoded: it stays in the enum
    // and the check stays in the code, but it is not a diagnosis anybody is owed
    // a fixture for, because no input can reach it.
    //
    // The sliver below does flip: clockwise before snapping and counterclockwise
    // after, measured against the shipped orientation<K> at spacing 1.0, where
    // (0, 0), (10, 0.4), (20, 0.6) becomes (0, 0), (10, 0), (20, 1). The noder
    // never gets to ask. A winding flip requires the ring's signed area to fall
    // below the rounding scale, which puts a vertex within half a cell of the
    // edge opposite it -- exactly the grazing condition segment_meets_cell
    // tests. Guarantee 14(b) then REQUIRES that edge to be split there, the
    // split repeats a node id, and NonSimpleRing fires one stage before ring
    // re-validation ever runs.
    //
    // The argument is structural, not an accident of these numbers: thin enough
    // to flip IS near enough to graze. Four candidate slivers were measured and
    // all four grazed; the residual -- a ring simple before snapping, simple and
    // node-repeat-free after it, and reversed -- was settled by search. See
    // "The bounded ring search, and the intermediate result that looks like a
    // refutation": a bounded exhaustive enumeration of 4-gons and 5-gons, run by
    // two independently written programs with different enumerations, simplicity
    // filters and lattices, reports FLIPPED == grazed and SURVIVORS == 0 in every
    // block. That is evidence and not a theorem -- two ring sizes, small lattices
    // -- and the design says so.
    //
    // MUTANT 8 (the winding re-check dropped, or a reversed ring silently
    // repaired) IS THEREFORE STRUCK: it has no killer input anywhere in this
    // increment, and a mutant with no killer is not a weaker mutant, it is a
    // claim the round would have shipped as covered. Its live half is mutant 10,
    // which kills the same silent wrong mesh through the node repeat that does
    // occur. 05-noder.md risk 3 -- a reversed hole meshed as an island with no
    // diagnostic anywhere -- is mitigated by guarantee 14(b) plus the node-repeat
    // check, one stage earlier than that document places it; the winding re-check
    // is defence in depth behind them. Mutant 8 is reinstated the moment the
    // search turns up a survivor at any ring size.
    PslgBuilder b;
    b.add_chain(square(-50.0, -50.0, 50.0, 50.0), ChainRole::Outer);
    b.add_chain(std::vector<Point2>{Point2{0.0, 0.0}, Point2{10.0, 0.4}, Point2{20.0, 0.6}},
                ChainRole::Hole);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{1.0});
    REQUIRE(out.status == NodeStatus::NonSimpleRing);
    REQUIRE_FALSE(out.ok());

    // What the fixture pins beyond the status itself: the hole is NOT silently
    // reversed and handed back as a mesh. 05-noder.md risk 3's failure mode is a
    // hole meshed as an island with no diagnostic, and an engaged pslg here
    // would be exactly that.
    REQUIRE_FALSE(out.pslg.has_value());
    REQUIRE_FALSE(out.message.empty());
}

TEST_CASE("a figure-eight ring is refused for revisiting a node", "[noding][node][status]") {
    // MUTANT 10: NonSimpleRing checked by a segment-pair scan that only finds
    // interior crossings. After noding a ring self-intersects IFF it repeats a
    // node id -- 14(a) has already excluded interior crossings that are not at
    // nodes -- so the check is a sort over std::uint32_t and the quadratic
    // geometric scan is the wrong implementation as well as the slower one.
    //
    // The bow-tie's crossing at (5, 5) becomes a node in step 2, both of its
    // edges are split there in step 4, and the closed chain then visits that
    // node twice. Which is the only form the case can take.
    PslgBuilder b;
    b.add_chain(std::vector<Point2>{Point2{0.0, 0.0}, Point2{10.0, 10.0}, Point2{10.0, 0.0},
                                    Point2{0.0, 10.0}},
                ChainRole::Outer);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{0.1});
    REQUIRE(out.status == NodeStatus::NonSimpleRing);
    REQUIRE_FALSE(out.ok());
}

// ---------------------------------------------------------------------------
// The loop, and the two mutants about how it ends.
// ---------------------------------------------------------------------------

TEST_CASE("input that needs a second round is not reported Ok after the first",
          "[noding][node][rounds]") {
    // MUTANTS 5 AND 6 together.
    //
    //   5  the loop condition spelled "until a round produces no new node"
    //      instead of "until the verification passes". A round can split an edge
    //      at an EXISTING node, so the first spelling reports converged in a
    //      state that violates 14(b).
    //   6  max_rounds treated as a success -- returning the last candidate
    //      instead of NotConverged.
    //
    // This is a SEARCH rather than a hand fixture, and the reason is measured:
    // 05-noder.md risk 9 reports that the near-parallel generator produced no
    // surviving crossings at 0.1 m, so a decimetre cascade is hard to write down
    // by hand. A COARSE grid deforms far more, and the search below finds the
    // cascade instead of guessing at it. The assertion is that at least one
    // input in the seed range needs more than one round and converges in eight:
    // an existence claim, which is exactly what kills both mutants.
    int needs_more_than_one = 0;
    for (std::uint64_t seed = 0; seed < 32; ++seed) {
        const Pslg in = terrain::testing::to_pslg(
            terrain::testing::random_polylines(seed, terrain::testing::PolylineOptions{}));

        const NodeOutcome capped = node<DefaultKernel>(in, NodeOptions{5.0, 1});
        const NodeOutcome full = node<DefaultKernel>(in, NodeOptions{5.0, 8});

        // Whatever the cap, a refusal never hands back a pslg. Mutant 6's own
        // assertion, made on every seed rather than only on the interesting one.
        if (capped.status != NodeStatus::Ok) {
            REQUIRE_FALSE(capped.pslg.has_value());
        }

        if (capped.status == NodeStatus::NotConverged && full.status == NodeStatus::Ok) {
            ++needs_more_than_one;
        }
    }
    REQUIRE(needs_more_than_one >= 1);
}

TEST_CASE("a larger cap never changes the answer, only whether there is one",
          "[noding][node][rounds]") {
    // max_rounds is a cap on a loop rather than a parameter of the answer: at any
    // value the outcome is either the same mesh or NotConverged, never a
    // different mesh. That is the whole reason it may have a default while
    // spacing may not.
    const Pslg in = crossing_pslg();
    const NodeOutcome few = node<DefaultKernel>(in, NodeOptions{0.1, 1});
    const NodeOutcome many = node<DefaultKernel>(in, NodeOptions{0.1, 16});

    REQUIRE(few.status == NodeStatus::Ok);
    REQUIRE(many.status == NodeStatus::Ok);
    REQUIRE(few.pslg->vertices().size() == many.pslg->vertices().size());
    for (std::size_t i = 0; i < few.pslg->vertices().size(); ++i) {
        REQUIRE(few.pslg->vertices()[i].x == many.pslg->vertices()[i].x);
        REQUIRE(few.pslg->vertices()[i].y == many.pslg->vertices()[i].y);
    }
}

TEST_CASE("input that is already noded comes back unchanged", "[noding][node][rounds]") {
    // The path where the noder must do nothing at all, which is the common case
    // on real breakline data and the one a driver that always splits something
    // would fail. The square's own corners are on the grid, nothing crosses, and
    // no node's cell meets a foreign edge.
    const SnapGrid grid{0.1};
    PslgBuilder b;
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{0, 0}), grid.world(GridPoint{100, 0}),
                                    grid.world(GridPoint{100, 100}), grid.world(GridPoint{0, 100})},
                ChainRole::Outer);
    b.add_chain(std::vector<Point2>{grid.world(GridPoint{20, 20}), grid.world(GridPoint{60, 60})},
                ChainRole::Breakline, kRoad);

    const NodeOutcome out = node<DefaultKernel>(built(std::move(b)), NodeOptions{0.1});
    REQUIRE(out.status == NodeStatus::Ok);
    const NodedPslg& p = *out.pslg;

    REQUIRE(p.vertices().size() == 6);
    REQUIRE(p.edge_count(0) == 4);
    REQUIRE(p.edge_count(1) == 1);
    REQUIRE(p.edge_properties().size() == 5);
    REQUIRE(p.edge_properties()[p.edge_base(1)] == kRoad);
    REQUIRE(p.edge_properties()[p.edge_base(0)] == EdgeProperties{});
}

TEST_CASE("node_of_input_vertex maps every input vertex onto its node",
          "[noding][node]") {
    // Guarantee 9 does not survive noding -- per-vertex index identity is lost,
    // because node ids are NodeSet's sorted order. This map is what a caller
    // holding an input vertex index uses instead, and 5c's Python surface
    // depends on it.
    const Pslg in = crossing_pslg();
    const NodeOutcome out = node<DefaultKernel>(in, NodeOptions{0.1});
    REQUIRE(out.ok());
    const NodedPslg& p = *out.pslg;

    REQUIRE(p.node_of_input_vertex().size() == in.vertices().size());
    for (std::size_t i = 0; i < in.vertices().size(); ++i) {
        const std::uint32_t n = p.node_of_input_vertex()[i];
        REQUIRE(n < p.vertices().size());

        // The node an input vertex maps to is its own snapped image, bitwise.
        const Point2 want = p.grid().snapped(in.vertices()[i]);
        REQUIRE(p.vertices()[n].x == want.x);
        REQUIRE(p.vertices()[n].y == want.y);
    }
}
