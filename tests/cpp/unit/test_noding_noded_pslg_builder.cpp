// Unit tests for terrain/noding/noded_pslg_builder.hpp -- the verifier.
//
// INVARIANT-CRITICAL, mutation round (docs/increments/05b-noder-driver.md,
// "What is worth testing"). This is the only suite in increment 5b that can
// fail when the *verifier* is broken, and that is the whole reason the verifier
// is a separate type rather than a pass inside node<K>'s loop:
//
//   A verification pass reachable only through the noder only ever sees inputs
//   the noder produced, so a verifier that agrees with its producer is
//   .claude/REQUIRED-READING.md's "do not write X is verified by Y until Y has
//   been run against a broken X".
//
// EVERY CANDIDATE IN THIS FILE IS BUILT BY HAND. Nothing here calls node<K>,
// and noding/node.hpp is deliberately not included: the moment this suite needs
// the driver to construct an input, it has stopped being able to fail.
//
// Coordinates are derived from GridPoints through grid.world(), never written as
// decimal literals, so guarantee 11 holds by construction in the fixtures that
// are not about violating it. At a non-dyadic spacing a literal that looks like
// world(g) usually is not one (core/snap_grid.hpp, and commit e412a43).
//
// ---------------------------------------------------------------------------
// TWO READINGS OF THE DESIGN ARE PINNED HERE AND BOTH ARE FLAGGED AS READINGS.
//
// (1) THE STATUS FOR A GUARANTEE-14 REFUSAL. 05b-noder-driver.md says the
//     *driver* maps a 14 refusal to NotConverged at the cap and to another
//     round below it, and maps 11/12/13/16 straight through as MalformedOutput.
//     It never names the status build() itself returns for 14. The driver
//     cannot distinguish the two cases unless build() does, so the only
//     coherent assignment is the one asserted below: NotConverged for a 14
//     refusal, MalformedOutput for 11/12/13/16. If the design means something
//     else it has to say so, and this file is where it will fail first.
//
// (2) THE BUILDER IS NOT THE PSLG VALIDATOR. The design lists what build()
//     checks -- 11, 12, 13, 14(a), 14(b), 16 -- and NoOuterChain is not among
//     them, so the fixtures here are breakline-only. Pinned as its own case.
// ---------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/noding/noded_pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <cstddef>
#include <cstdint>
#include <set>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

using terrain::Chain;
using terrain::ChainRole;
using terrain::EdgeProperties;
using terrain::GridPoint;
using terrain::Point2;
using terrain::Segment2;
using terrain::SnapGrid;
using terrain::noding::NodedPslgBuilder;
using terrain::noding::NodeOutcome;
using terrain::noding::NodeStatus;
using terrain::noding::describe;
using terrain::noding::segment_meets_cell;
using terrain::on_segment;
using terrain::pred::DefaultKernel;

namespace {

constexpr double kSpacing = 0.1;

// Reading (1) above, given a name so that a design ruling changing it changes
// one line of this file rather than nine.
constexpr NodeStatus kNotNoded = NodeStatus::NotConverged;
constexpr NodeStatus kOurBug = NodeStatus::MalformedOutput;

const EdgeProperties kRoad = EdgeProperties::bit(0);
const EdgeProperties kRiver = EdgeProperties::bit(1);

template <class B>
concept BuildableAsLvalue = requires(B& b) { b.template build<DefaultKernel>(); };

template <class B>
concept BuildableAsRvalue = requires(B&& b) { std::move(b).template build<DefaultKernel>(); };

// A hand-built candidate, in the order NodedPslgBuilder's constructor takes it.
struct Candidate {
    SnapGrid grid{kSpacing};
    std::vector<Point2> vertices;
    std::vector<Chain> chains;
    std::vector<std::uint32_t> chain_indices;
    std::vector<EdgeProperties> edge_properties;
    std::vector<std::uint32_t> node_of_input_vertex;

    [[nodiscard]] NodeOutcome verify() const {
        return NodedPslgBuilder{grid,
                                vertices,
                                chains,
                                chain_indices,
                                edge_properties,
                                node_of_input_vertex}
            .build<DefaultKernel>();
    }
};

[[nodiscard]] std::vector<Point2> world_of(const SnapGrid& grid,
                                           const std::vector<GridPoint>& nodes) {
    std::vector<Point2> out;
    out.reserve(nodes.size());
    for (const GridPoint& g : nodes) {
        out.push_back(grid.world(g));
    }
    return out;
}

// Builds a candidate from node positions and one index run per open breakline
// chain, filling the per-edge property array densely from the chains' own sets.
// The array's *values* are not this suite's business -- build() does not check
// guarantee 15 (05b-noder-driver.md, "Why guarantee 15 is not the builder's"),
// and prop_noding_no_crossings.cpp owns it -- but its *length* is, because a
// length that disagrees with edge_base(chains().size()) is malformed output.
[[nodiscard]] Candidate breaklines(const std::vector<GridPoint>& nodes,
                                   const std::vector<std::vector<std::uint32_t>>& runs,
                                   const std::vector<EdgeProperties>& props) {
    Candidate c;
    c.vertices = world_of(c.grid, nodes);
    for (std::size_t i = 0; i < runs.size(); ++i) {
        const auto begin = static_cast<std::uint32_t>(c.chain_indices.size());
        c.chain_indices.insert(c.chain_indices.end(), runs[i].begin(), runs[i].end());
        c.chains.push_back(Chain{begin, static_cast<std::uint32_t>(runs[i].size()),
                                 ChainRole::Breakline, props[i]});
        for (std::size_t k = 0; k + 1 < runs[i].size(); ++k) {
            c.edge_properties.push_back(props[i]);
        }
    }
    for (std::uint32_t n = 0; n < nodes.size(); ++n) {
        c.node_of_input_vertex.push_back(n);
    }
    return c;
}

// Two parallel breaklines half a metre apart. Every guarantee holds.
[[nodiscard]] Candidate valid_candidate() {
    return breaklines({GridPoint{0, 0}, GridPoint{0, 50}, GridPoint{100, 0}, GridPoint{100, 50}},
                      {{0, 2}, {1, 3}}, {kRoad, kRiver});
}

}  // namespace

// ---------------------------------------------------------------------------
// The outcome type and the status vocabulary.
// ---------------------------------------------------------------------------

TEST_CASE("a default-constructed NodeOutcome is NotRun and carries no pslg",
          "[noding][noded_pslg_builder]") {
    const NodeOutcome out;
    REQUIRE(out.status == NodeStatus::NotRun);
    REQUIRE_FALSE(out.ok());
    REQUIRE_FALSE(out.pslg.has_value());
    REQUIRE(out.message.empty());
}

TEST_CASE("describe is total over NodeStatus and gives every status its own words",
          "[noding][noded_pslg_builder]") {
    constexpr NodeStatus all[] = {
        NodeStatus::Ok,          NodeStatus::NotRun,
        NodeStatus::InvalidSnapSpacing, NodeStatus::CoordinateOutOfRange,
        NodeStatus::RingCollapsed,      NodeStatus::RingDegenerateAfterSnap,
        NodeStatus::NonSimpleRing,      NodeStatus::NotConverged,
        NodeStatus::MalformedOutput};

    std::set<std::string_view> seen;
    for (const NodeStatus s : all) {
        const std::string_view text = describe(s);
        REQUIRE_FALSE(text.empty());
        seen.insert(text);
    }
    REQUIRE(seen.size() == 9);

    // constexpr, as declared: usable in a static assertion and therefore in a
    // switch a compiler can fold.
    STATIC_REQUIRE(!describe(NodeStatus::Ok).empty());
}

TEST_CASE("ok() is exactly pslg engagement, on both arms", "[noding][noded_pslg_builder]") {
    const NodeOutcome good = valid_candidate().verify();
    REQUIRE(good.status == NodeStatus::Ok);
    REQUIRE(good.ok());
    REQUIRE(good.pslg.has_value());

    // A successful build says nothing, and a refusal always says something.
    // Spelled as two REQUIREs rather than as `ok() != message.empty()`: that
    // one-liner is the form increment 4's design got inverted (04-cdt.md:346,
    // corrected in tests/cpp/unit/test_cdt_backend_seam.cpp).
    REQUIRE(good.message.empty());

    Candidate broken = valid_candidate();
    broken.vertices.push_back(broken.vertices[0]);  // guarantee 12
    const NodeOutcome bad = broken.verify();
    REQUIRE_FALSE(bad.ok());
    REQUIRE_FALSE(bad.pslg.has_value());
    REQUIRE_FALSE(bad.message.empty());
}

TEST_CASE("build is rvalue-ref-qualified, so a builder cannot be verified twice",
          "[noding][noded_pslg_builder]") {
    // The same shape as PslgBuilder::build: the candidate arrays are moved into
    // the NodedPslg, so "build twice and get two objects sharing nothing" is a
    // compile error rather than a subtle question about what the second sees.
    STATIC_REQUIRE(!BuildableAsLvalue<NodedPslgBuilder>);
    STATIC_REQUIRE(BuildableAsRvalue<NodedPslgBuilder>);
}

TEST_CASE("the builder is not the Pslg validator", "[noding][noded_pslg_builder]") {
    // Breakline-only, no outer chain, no winding contract anywhere. The design
    // lists what build() checks and NoOuterChain is not on it; a builder that
    // refused this would be re-running increment 3's stage 1 on data that has
    // already been through it. Reading (2) at the head of this file.
    const NodeOutcome out = valid_candidate().verify();
    REQUIRE(out.status == NodeStatus::Ok);
}

// ---------------------------------------------------------------------------
// Guarantee 11 -- every coordinate is world(g), within range.
// ---------------------------------------------------------------------------

TEST_CASE("an off-grid coordinate is refused as malformed output",
          "[noding][noded_pslg_builder]") {
    Candidate c = valid_candidate();
    // Three tenths of a cell off the lattice: snapped(v) != v bitwise, which is
    // guarantee 11's own checkable spelling.
    c.vertices[2].x += 0.03;
    REQUIRE(c.grid.snapped(c.vertices[2]).x != c.vertices[2].x);

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
    REQUIRE_FALSE(out.ok());
}

TEST_CASE("a coordinate past kMaxGridIndex is refused, not asserted on",
          "[noding][noded_pslg_builder]") {
    // MUTANT 14's verifier-side twin. can_snap must be consulted BEFORE
    // snapped(), because snap() debug-asserts can_snap and cell_min/cell_max
    // debug-assert the index range: a verifier that checks 11 by comparing
    // snapped(v) == v without the range test first aborts here in a Debug build
    // and silently wraps an index in Release. Both are worse than a status.
    Candidate c = valid_candidate();
    c.vertices[1] = Point2{1e300, 0.0};
    REQUIRE_FALSE(c.grid.can_snap(c.vertices[1]));

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
    REQUIRE_FALSE(out.ok());
}

// ---------------------------------------------------------------------------
// Guarantees 12 and 13.
// ---------------------------------------------------------------------------

TEST_CASE("two equal vertices are refused even when neither is referenced twice",
          "[noding][noded_pslg_builder]") {
    Candidate c = valid_candidate();
    c.vertices.push_back(c.vertices[0]);  // an unreferenced duplicate
    c.node_of_input_vertex.push_back(0);

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
}

TEST_CASE("a chain that repeats a consecutive index is a zero-length edge",
          "[noding][noded_pslg_builder]") {
    Candidate c = breaklines(
        {GridPoint{0, 0}, GridPoint{0, 50}, GridPoint{100, 0}, GridPoint{100, 50}},
        {{0, 0, 2}, {1, 3}}, {kRoad, kRiver});

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
}

// ---------------------------------------------------------------------------
// Guarantee 14(a) -- no two edges cross in their interiors, and the amendment.
// ---------------------------------------------------------------------------

TEST_CASE("two edges crossing at a point that is no node are refused",
          "[noding][noded_pslg_builder]") {
    // The increment's reason for existing, handed to the verifier directly:
    // two diagonals meeting at (5, 5) in world units, which is a node of
    // neither chain.
    Candidate c = breaklines(
        {GridPoint{0, 0}, GridPoint{0, 100}, GridPoint{100, 0}, GridPoint{100, 100}},
        {{0, 3}, {1, 2}}, {kRoad, kRiver});

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kNotNoded);
    REQUIRE_FALSE(out.ok());
}

TEST_CASE("a partial collinear overlap is refused", "[noding][noded_pslg_builder]") {
    // The amended 14(a)'s other arm: collinear, sharing a sub-segment of
    // positive length, with UNEQUAL node-id pairs.
    //
    // This fixture cannot isolate 14(a) and nothing can: the inner endpoint of
    // a partial overlap lies on the other edge, so its cell meets an edge it is
    // not an endpoint of and 14(b) is violated too. A partial overlap is
    // therefore always a double refusal. Recorded rather than worked around --
    // the test's job is that it is refused, not which clause names it.
    Candidate c = breaklines(
        {GridPoint{0, 0}, GridPoint{50, 0}, GridPoint{100, 0}, GridPoint{150, 0}},
        {{0, 2}, {1, 3}}, {kRoad, kRiver});

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kNotNoded);
}

TEST_CASE("a road running along a river is accepted: duplicate edges are legal",
          "[noding][noded_pslg_builder]") {
    // THE ACCEPTANCE, and it is as load-bearing as any refusal in this file.
    // 05-noder.md's unamended 14(a) required Disjoint or Touching for every
    // candidate pair, and classify on two identical segments returns
    // Overlapping -- so the pre-amendment clause refuses the one input
    // guarantee 15 exists for. 05b-noder-driver.md, "Guarantee 14, amended".
    Candidate c = breaklines({GridPoint{0, 0}, GridPoint{100, 0}}, {{0, 1}, {0, 1}},
                             {kRoad, kRiver});
    // Faithful to guarantee 15 even though build() does not check it.
    c.edge_properties = {kRoad | kRiver, kRoad | kRiver};

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == NodeStatus::Ok);
    REQUIRE(out.ok());

    // Chain identity survives the duplication: two chains in, two chains out,
    // each keeping its own set. Guarantee 16, and the reason the dedup is over
    // the edge set rather than over chains.
    REQUIRE(out.pslg->chains().size() == 2);
    REQUIRE(out.pslg->chains()[0].properties == kRoad);
    REQUIRE(out.pslg->chains()[1].properties == kRiver);
}

// ---------------------------------------------------------------------------
// Guarantee 14(b) -- the hot-pixel clause, and the one mutant only this suite
// can kill.
// ---------------------------------------------------------------------------

TEST_CASE("a live T-junction is refused, and proximity is the relation that finds it",
          "[noding][noded_pslg_builder]") {
    // MUTANT 13: the verifier's 14(b) returning true unconditionally, or asking
    // on_segment instead of segment_meets_cell. This fixture is the only thing
    // in the increment that kills it, because no property over the noder's own
    // output can -- a noder and a verifier sharing the wrong relation agree.
    //
    // The host runs from world (0, 0) to world (0.4, 0.2); the intruding node
    // sits at world (0.1, 0), where the host's own y is 0.05. The node is NOT on
    // the host by any exact reading. Its closed cell spans y in [-0.05, 0.05]
    // and x in [0.05, 0.15], and the host passes through that cell between
    // y = 0.025 and y = 0.05. That is what snap rounding means by a T-junction,
    // and it is 93 % of the real ones (05-noder.md, guarantee 14).
    const SnapGrid grid{kSpacing};
    const Segment2 host{grid.world(GridPoint{0, 0}), grid.world(GridPoint{4, 2})};
    const GridPoint intruder{1, 0};

    REQUIRE(segment_meets_cell<DefaultKernel>(grid, host, intruder));
    REQUIRE_FALSE(on_segment<DefaultKernel>(host, grid.world(intruder)));

    // Node ids in lexicographic grid order: (0,0)=0, (1,-10)=1, (1,0)=2, (4,2)=3.
    Candidate c = breaklines(
        {GridPoint{0, 0}, GridPoint{1, -10}, GridPoint{1, 0}, GridPoint{4, 2}},
        {{0, 3}, {2, 1}}, {kRoad, kRiver});

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kNotNoded);
    REQUIRE_FALSE(out.ok());
}

TEST_CASE("14(b) is asked of every node against every edge, not within a chain",
          "[noding][noded_pslg_builder]") {
    // The same geometry with the intruder moved out of the host's chain is
    // already the case above -- the intruder is chain 1's. This case is its
    // complement: the intruder is the host's OWN chain's far endpoint, so a
    // verifier that skips same-chain pairs to save time accepts it.
    const SnapGrid grid{kSpacing};
    const Segment2 host{grid.world(GridPoint{0, 0}), grid.world(GridPoint{4, 2})};
    REQUIRE(segment_meets_cell<DefaultKernel>(grid, host, GridPoint{1, 0}));

    // One chain: (0,0) -> (4,2) -> (1,0). Its second edge is legitimate, but
    // node (1,0) is not an endpoint of the FIRST edge and its cell meets it.
    // Node ids lexicographic: (0,0)=0, (1,0)=1, (4,2)=2.
    Candidate c = breaklines({GridPoint{0, 0}, GridPoint{1, 0}, GridPoint{4, 2}},
                             {{0, 2, 1}}, {kRoad});

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kNotNoded);
}

// ---------------------------------------------------------------------------
// Guarantee 16 -- the chain table's internal consistency. The *correspondence*
// to the input Pslg is node<K>'s to establish and test_noding_node.cpp's to
// check; what the builder can see is that the slices tile the flat buffer.
// ---------------------------------------------------------------------------

TEST_CASE("chain slices that overlap are refused", "[noding][noded_pslg_builder]") {
    Candidate c = valid_candidate();
    c.chains[0].count = 3;  // now covers indices 0..2, and chain 1 begins at 2

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
}

TEST_CASE("a gap between chain slices is refused", "[noding][noded_pslg_builder]") {
    Candidate c = valid_candidate();
    c.chain_indices.push_back(0);  // an index no chain claims
    c.chains[1].begin = 3;
    c.chains[1].count = 2;

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
}

TEST_CASE("a chain slice running past the flat buffer is refused",
          "[noding][noded_pslg_builder]") {
    Candidate c = valid_candidate();
    c.chains[1].count = 3;

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
}

TEST_CASE("a node id past the vertex buffer is refused", "[noding][noded_pslg_builder]") {
    Candidate c = valid_candidate();
    c.chain_indices[3] = 99;

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
}

TEST_CASE("an edge_properties array whose length is not edge_base(n) is refused",
          "[noding][noded_pslg_builder]") {
    // Guarantee 15's *values* are out of the builder's reach. Its *length* is
    // not: a dense array is defined by its length, and a short one makes every
    // caller's edge_base(c) + k a read past the end.
    Candidate c = valid_candidate();
    c.edge_properties.pop_back();

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
}

// ---------------------------------------------------------------------------
// THE ASYMMETRY OF THE count < 3 GUARD, PINNED FROM BOTH SIDES.
//
// 05b-noder-driver.md, "An open chain that collapses to one node is Ok, and it
// is pinned": the refusal belongs to is_closed(role) ALONE, and an open chain of
// one node is accepted BY DESIGN rather than tolerated. A ring of under three
// nodes has no interior and no winding, so orientation<K> has nothing to answer;
// an open chain of one node is asked no question it cannot answer, and its
// edge_count is the empty range every loop already handles.
//
// Both halves are needed and neither is redundant -- but NOT for the reason the
// first draft of this comment gave, and the measurement is worth carrying.
// Widening the guard to a bare `count < 3` was said to be free; it is not. It
// fails 26 registered cases across four suites, because every open chain in
// `valid_candidate()` has `count == 2` and the widening refuses them all. The
// mutant that really was unpinned is `|| count < 2` -- how a maintainer would
// actually spell "an open chain of one node ought to be refused" -- and it
// fails exactly two cases, the acceptance case below and its end-to-end
// counterpart in test_noding_node.cpp, and nothing else. That "and nothing
// else" is the evidence the acceptance was uncovered before them. Deleting the
// guard outright looks like the same correction from the other end, and until
// 797144c the closed arm it guards had no test at all. The two cases below are one
// pair and should be read as one.
// ---------------------------------------------------------------------------

TEST_CASE("an open chain of one node is accepted, contributing no edge",
          "[noding][noded_pslg_builder]") {
    // A third breakline whose every vertex snapped to the same node: count == 1,
    // zero edges, five metres clear of the other two so that no other guarantee
    // can be the reason for the verdict. Arrays the noder did not produce, which
    // is the whole reason the builder is a separate type.
    Candidate c = breaklines({GridPoint{0, 0}, GridPoint{0, 50}, GridPoint{100, 0},
                              GridPoint{100, 50}, GridPoint{50, 500}},
                             {{0, 2}, {1, 3}, {4}}, {kRoad, kRiver, kRoad});

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == NodeStatus::Ok);
    REQUIRE(out.message.empty());

    const auto& p = *out.pslg;
    REQUIRE(p.chains().size() == 3);
    REQUIRE(p.edge_count(2) == 0);
    REQUIRE(p.edge_base(2) == p.edge_base(3));
    REQUIRE(p.edge_properties().size() == p.edge_base(p.chains().size()));
}

TEST_CASE("a closed chain of under three nodes is refused, which is the other half",
          "[noding][noded_pslg_builder]") {
    Candidate c = breaklines({GridPoint{0, 0}, GridPoint{0, 50}}, {{0, 1}}, {kRoad});
    c.chains[0].role = ChainRole::Outer;  // the same two nodes, now a ring

    const NodeOutcome out = c.verify();
    REQUIRE(out.status == kOurBug);
    // Named, because the same candidate also has a closed chain's edge count
    // disagreeing with its property array and either would refuse it.
    REQUIRE(out.message.find("closed chain") != std::string::npos);
}

// ---------------------------------------------------------------------------
// The accepted candidate's own arithmetic, so that an acceptance is not merely
// the absence of a refusal.
// ---------------------------------------------------------------------------

TEST_CASE("an accepted candidate hands back exactly the arrays it was given",
          "[noding][noded_pslg_builder]") {
    const Candidate c = valid_candidate();
    const NodeOutcome out = c.verify();
    REQUIRE(out.status == NodeStatus::Ok);

    const auto& p = *out.pslg;
    REQUIRE(p.vertices().size() == c.vertices.size());
    REQUIRE(p.chains().size() == c.chains.size());
    REQUIRE(p.chain_indices().size() == c.chain_indices.size());
    REQUIRE(p.edge_properties().size() == c.edge_properties.size());
    REQUIRE(p.grid().spacing() == kSpacing);

    for (std::size_t i = 0; i < c.vertices.size(); ++i) {
        REQUIRE(p.vertices()[i].x == c.vertices[i].x);
        REQUIRE(p.vertices()[i].y == c.vertices[i].y);
    }
    REQUIRE(p.edge_base(p.chains().size()) == p.edge_properties().size());
}
