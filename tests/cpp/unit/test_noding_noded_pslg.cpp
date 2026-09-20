// Unit tests for terrain/core/noded_pslg.hpp -- increment 5b's output type.
//
// NOT invariant-critical, and NO mutation round. The reason is
// docs/increments/05b-noder-driver.md, "What is worth testing": this is an
// accessor suite over a value type whose failure mode is a typo and is loud.
// The one accessor with arithmetic in it, edge_base, is additionally checked
// against edge_count's running sum in prop_noding_no_crossings.cpp, on noder
// output rather than on a hand-built candidate, which is the stronger check.
//
// THIS FILE INCLUDES noding/noded_pslg_builder.hpp AND THAT IS NOT A GATE B
// VIOLATION. NodedPslg's constructor is private with exactly one friend, so the
// builder is the only route to an instance; the suite tests the accessors and
// the builder is the fixture. Nothing here asserts a builder refusal --
// test_noding_noded_pslg_builder.cpp owns every one of those.
//
// Every candidate below is built from GridPoints through grid.world(), never
// from decimal literals. That is not a style preference: guarantee 11 says
// every coordinate is world(g), and a literal that looks like world(g) at a
// non-dyadic spacing usually is not one (snap_grid.hpp's header comment, and
// commit e412a43). Deriving the coordinates makes the fixtures satisfy 11 by
// construction instead of by luck.

#include <catch2/catch_test_macros.hpp>

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/noded_pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

using terrain::Chain;
using terrain::ChainRole;
using terrain::EdgeProperties;
using terrain::GridPoint;
using terrain::NodedPslg;
using terrain::Point2;
using terrain::Pslg;
using terrain::Segment2;
using terrain::SnapGrid;
using terrain::noding::NodedPslgBuilder;
using terrain::noding::NodeOutcome;
using terrain::noding::NodeStatus;
using terrain::pred::DefaultKernel;

namespace {

// Decimal on purpose: the spacing at which world() is not an affine map, so a
// suite that accidentally depends on exact arithmetic says so here rather than
// in the property file.
constexpr double kSpacing = 0.1;

// The reference candidate: a counterclockwise outer square four cells on a
// side, plus one breakline diagonal well inside it. Node ids are the
// lexicographic order of the grid points, which is what NodeSet produces.
//
//   id 0 = (  0,   0)      id 3 = ( 80,  80)
//   id 1 = (  0, 100)      id 4 = (100,   0)
//   id 2 = ( 20,  20)      id 5 = (100, 100)
//
// No node's cell meets an edge it is not an endpoint of, and no two edges cross:
// the diagonal runs from (2.0, 2.0) to (8.0, 8.0) in world units inside a square
// spanning [0, 10]. Guarantee 14 therefore holds and build() must accept it.
struct Reference {
    SnapGrid grid{kSpacing};
    std::vector<GridPoint> nodes{GridPoint{0, 0},   GridPoint{0, 100}, GridPoint{20, 20},
                                 GridPoint{80, 80}, GridPoint{100, 0}, GridPoint{100, 100}};
    std::vector<Chain> chains{Chain{0, 4, ChainRole::Outer, EdgeProperties::bit(0)},
                              Chain{4, 2, ChainRole::Breakline, EdgeProperties::bit(1)}};
    std::vector<std::uint32_t> chain_indices{0, 4, 5, 1, 2, 3};

    // Five flat edge positions: four ring edges then one breakline edge.
    std::vector<EdgeProperties> edge_properties{
        EdgeProperties::bit(0), EdgeProperties::bit(0), EdgeProperties::bit(0),
        EdgeProperties::bit(0), EdgeProperties::bit(1)};

    // Six input vertices in this fiction, one per node, in input order
    // (outer ring first, then the breakline) rather than node order -- the
    // whole point of the map is that the two orders differ.
    std::vector<std::uint32_t> node_of_input_vertex{0, 4, 5, 1, 2, 3};

    [[nodiscard]] std::vector<Point2> vertices() const {
        std::vector<Point2> out;
        out.reserve(nodes.size());
        for (const GridPoint& g : nodes) {
            out.push_back(grid.world(g));
        }
        return out;
    }

    [[nodiscard]] NodeOutcome build() const {
        return NodedPslgBuilder{grid, vertices(), chains, chain_indices, edge_properties,
                                node_of_input_vertex}
            .build<DefaultKernel>();
    }
};

[[nodiscard]] NodedPslg reference_pslg() {
    NodeOutcome out = Reference{}.build();
    REQUIRE(out.status == NodeStatus::Ok);
    REQUIRE(out.ok());
    return std::move(*out.pslg);
}

}  // namespace

// ---------------------------------------------------------------------------
// The type's shape. These are compile-time and are the architectural claim:
// docs/increments/05b-noder-driver.md, "NodedPslg is not a subclass of Pslg and
// there is no conversion", restated because the CDT signature change at 5c is
// what makes un-noded input unrepresentable, and a silent conversion would undo
// the entire increment.
// ---------------------------------------------------------------------------

TEST_CASE("NodedPslg is not reachable from a Pslg by any conversion", "[noding][noded_pslg]") {
    STATIC_REQUIRE(!std::is_base_of_v<Pslg, NodedPslg>);
    STATIC_REQUIRE(!std::is_convertible_v<Pslg, NodedPslg>);
    STATIC_REQUIRE(!std::is_constructible_v<NodedPslg, const Pslg&>);
    STATIC_REQUIRE(!std::is_constructible_v<Pslg, const NodedPslg&>);
}

TEST_CASE("NodedPslg has no public constructor and is copyable and movable",
          "[noding][noded_pslg]") {
    STATIC_REQUIRE(!std::is_default_constructible_v<NodedPslg>);
    STATIC_REQUIRE(std::is_copy_constructible_v<NodedPslg>);
    STATIC_REQUIRE(std::is_move_constructible_v<NodedPslg>);
}

TEST_CASE("every NodedPslg accessor is const and the read-only ones are noexcept",
          "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();

    STATIC_REQUIRE(noexcept(p.vertices()));
    STATIC_REQUIRE(noexcept(p.chains()));
    STATIC_REQUIRE(noexcept(p.chain_indices()));
    STATIC_REQUIRE(noexcept(p.indices_of(0)));
    STATIC_REQUIRE(noexcept(p.edge_count(0)));
    STATIC_REQUIRE(noexcept(p.edge(0, 0)));
    STATIC_REQUIRE(noexcept(p.grid()));
    STATIC_REQUIRE(noexcept(p.edge_properties()));
    STATIC_REQUIRE(noexcept(p.node_of_input_vertex()));
    STATIC_REQUIRE(noexcept(p.edge_base(0)));

    // ring() is deliberately NOT noexcept, for the reason core/pslg.hpp gives
    // for its own: IndexedRing's constructor rechecks invariants the compiler
    // cannot see are already established, and a noexcept resting on that is a
    // std::terminate waiting for the argument to stop being true.
    STATIC_REQUIRE(!noexcept(p.ring(0)));
}

// ---------------------------------------------------------------------------
// The buffers, unchanged from what the builder consumed.
// ---------------------------------------------------------------------------

TEST_CASE("the vertex buffer is the candidate's, coordinate for coordinate",
          "[noding][noded_pslg]") {
    const Reference ref;
    const NodedPslg p = reference_pslg();
    const std::vector<Point2> expected = ref.vertices();

    REQUIRE(p.vertices().size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        // Bitwise, not approximate. Guarantee 11 is a statement about equality
        // of doubles and nothing weaker would mean anything.
        REQUIRE(p.vertices()[i].x == expected[i].x);
        REQUIRE(p.vertices()[i].y == expected[i].y);
    }
}

TEST_CASE("the chain table survives with its roles and property sets", "[noding][noded_pslg]") {
    const Reference ref;
    const NodedPslg p = reference_pslg();

    REQUIRE(p.chains().size() == 2);
    REQUIRE(p.chains()[0].role == ChainRole::Outer);
    REQUIRE(p.chains()[1].role == ChainRole::Breakline);
    REQUIRE(p.chains()[0].properties == ref.chains[0].properties);
    REQUIRE(p.chains()[1].properties == ref.chains[1].properties);
}

TEST_CASE("chain_indices is the flat buffer and indices_of is its slice",
          "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();

    REQUIRE(p.chain_indices().size() == 6);

    const auto outer = p.indices_of(0);
    REQUIRE(outer.size() == 4);
    REQUIRE(outer[0] == 0);
    REQUIRE(outer[1] == 4);
    REQUIRE(outer[2] == 5);
    REQUIRE(outer[3] == 1);

    const auto breakline = p.indices_of(1);
    REQUIRE(breakline.size() == 2);
    REQUIRE(breakline[0] == 2);
    REQUIRE(breakline[1] == 3);

    // The slices partition the flat buffer with no gap and no overlap. That is
    // guarantee 16's structural half and the builder enforces it; asserted here
    // because a caller reading chain_indices() directly depends on it.
    REQUIRE(p.chains()[0].begin == 0);
    REQUIRE(p.chains()[0].begin + p.chains()[0].count == p.chains()[1].begin);
    REQUIRE(p.chains()[1].begin + p.chains()[1].count == p.chain_indices().size());
}

// ---------------------------------------------------------------------------
// The two accessors with a role branch, and the one with arithmetic.
// ---------------------------------------------------------------------------

TEST_CASE("edge_count counts the closing edge for a ring and not for a breakline",
          "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();

    REQUIRE(p.edge_count(0) == 4);  // four distinct vertices, four edges
    REQUIRE(p.edge_count(1) == 1);  // two vertices, one open edge
}

TEST_CASE("edge is directed and the ring's last edge is its closing edge",
          "[noding][noded_pslg]") {
    const Reference ref;
    const NodedPslg p = reference_pslg();
    const std::vector<Point2> v = ref.vertices();

    const Segment2 first = p.edge(0, 0);
    REQUIRE(first.a.x == v[0].x);
    REQUIRE(first.a.y == v[0].y);
    REQUIRE(first.b.x == v[4].x);
    REQUIRE(first.b.y == v[4].y);

    // (k + 1) % n unconditionally: edge(0, 3) runs from the last stored index
    // back to the first. This is the accessor's whole reason for existing --
    // see core/pslg.hpp's comment on the same method.
    const Segment2 closing = p.edge(0, 3);
    REQUIRE(closing.a.x == v[1].x);
    REQUIRE(closing.a.y == v[1].y);
    REQUIRE(closing.b.x == v[0].x);
    REQUIRE(closing.b.y == v[0].y);

    const Segment2 open = p.edge(1, 0);
    REQUIRE(open.a.x == v[2].x);
    REQUIRE(open.b.x == v[3].x);
}

TEST_CASE("edge_base is the prefix sum of edge_count with one extra entry",
          "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();

    REQUIRE(p.edge_base(0) == 0);
    REQUIRE(p.edge_base(1) == 4);
    REQUIRE(p.edge_base(2) == 5);

    // The convention chain_indices already uses: one entry past the end, equal
    // to the total. Without it every caller holding (c, k) recomputes a prefix
    // sum, which is the derived arithmetic that is right in four places and
    // wrong in the fifth (05b-noder-driver.md, "Two additions").
    REQUIRE(p.edge_base(p.chains().size()) == p.edge_properties().size());

    std::size_t running = 0;
    for (std::size_t c = 0; c < p.chains().size(); ++c) {
        REQUIRE(p.edge_base(c) == running);
        running += p.edge_count(c);
    }
    REQUIRE(running == p.edge_properties().size());
}

TEST_CASE("edge_properties is dense and index-aligned with the flat enumeration",
          "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();

    // Dense, not a sparse override set: one entry per output edge, and the
    // empty set is a legal value meaning *unclassified*, never *absent*
    // (05b-noder-driver.md, "The array is dense").
    REQUIRE(p.edge_properties().size() == 5);

    for (std::size_t k = 0; k < p.edge_count(0); ++k) {
        REQUIRE(p.edge_properties()[p.edge_base(0) + k] == EdgeProperties::bit(0));
    }
    REQUIRE(p.edge_properties()[p.edge_base(1) + 0] == EdgeProperties::bit(1));
}

TEST_CASE("edge_properties is an EdgeProperties span and not a bare word",
          "[noding][noded_pslg]") {
    // 05b-noder-driver.md, "No bare word per edge": a std::uint32_t array here
    // would reinstate `if (edge_is_river)` one level below the field that closed
    // it, at a place where nothing re-validates. The type is the enforcement and
    // this static assertion is the tripwire on it.
    using Element = std::remove_cvref_t<decltype(reference_pslg().edge_properties()[0])>;
    STATIC_REQUIRE(std::is_same_v<Element, EdgeProperties>);
    STATIC_REQUIRE(!std::is_convertible_v<EdgeProperties, std::uint32_t>);
    STATIC_REQUIRE(!std::is_convertible_v<EdgeProperties, bool>);
}

TEST_CASE("node_of_input_vertex maps input positions onto node ids", "[noding][noded_pslg]") {
    const Reference ref;
    const NodedPslg p = reference_pslg();

    REQUIRE(p.node_of_input_vertex().size() == ref.node_of_input_vertex.size());
    for (std::size_t i = 0; i < ref.node_of_input_vertex.size(); ++i) {
        const std::uint32_t node = p.node_of_input_vertex()[i];
        REQUIRE(node == ref.node_of_input_vertex[i]);
        REQUIRE(node < p.vertices().size());
    }
}

TEST_CASE("grid reports the spacing the candidate was snapped on", "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();
    REQUIRE(p.grid().spacing() == kSpacing);

    // Guarantee 11, checkable bitwise and checked that way.
    for (const Point2& v : p.vertices()) {
        const Point2 again = p.grid().snapped(v);
        REQUIRE(again.x == v.x);
        REQUIRE(again.y == v.y);
    }
}

TEST_CASE("ring is available for a closed chain and is the chain's vertices",
          "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();

    REQUIRE_NOTHROW(p.ring(0));
    const auto r = p.ring(0);
    REQUIRE(r.size() == 4);
    REQUIRE(r.vertex(0).x == p.vertices()[0].x);
    REQUIRE(r.vertex(2).x == p.vertices()[5].x);
}

// ---------------------------------------------------------------------------
// Guarantee 10's analogue: a const NodedPslg needs no synchronisation. Same
// shape as prop_pslg_invariants.cpp's concurrent read, and the same reason --
// 5c releases the GIL around the noder and the renderer reads the result from
// whichever thread it lands on.
// ---------------------------------------------------------------------------

TEST_CASE("a const NodedPslg is readable from many threads with no synchronisation",
          "[noding][noded_pslg]") {
    const NodedPslg p = reference_pslg();

    // The reduction is fixed-order within each thread, so every thread must
    // produce the same double bit for bit. Computed once here first, so the
    // comparison is against a value and not against the threads' consensus --
    // eight threads agreeing on a wrong answer would otherwise pass.
    const auto reduce = [&p]() {
        double x_sum = 0.0;
        for (std::size_t c = 0; c < p.chains().size(); ++c) {
            for (std::size_t k = 0; k < p.edge_count(c); ++k) {
                x_sum += p.edge(c, k).a.x + p.edge(c, k).b.y;
            }
        }
        return x_sum;
    };
    const double expected = reduce();

    constexpr int thread_count = 8;
    std::atomic<int> agreements{0};
    std::vector<std::thread> threads;
    threads.reserve(thread_count);

    for (int t = 0; t < thread_count; ++t) {
        threads.emplace_back([&reduce, &agreements, expected]() {
            if (reduce() == expected) {
                agreements.fetch_add(1);
            }
        });
    }
    for (std::thread& th : threads) {
        th.join();
    }
    REQUIRE(agreements.load() == thread_count);
}
