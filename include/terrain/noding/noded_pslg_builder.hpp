#pragma once

// The verifier, and the only producer of a NodedPslg.
//
// WHY THIS IS A TYPE RATHER THAN A PASS INSIDE node<K>'s LOOP. A verification
// pass reachable only through the noder only ever sees inputs the noder
// produced, so a verifier that agrees with its producer is
// .claude/REQUIRED-READING.md's "do not write X is verified by Y until Y has
// been run against a broken X". As a separate type with a public entry point it
// can be handed a hand-built candidate with a live T-junction in it and be
// REQUIRED to refuse, which is what
// tests/cpp/unit/test_noding_noded_pslg_builder.cpp does and what no property
// over the noder's own output could.
//
// HOW TO READ NodeStatus, and it belongs here rather than in a suite comment:
// NodeStatus IS THE STATUS OF AN OUTCOME, NOT OF A CHECK. build() reports the
// status the driver would report if this were the last round. Under that
// reading the driver's mapping at the cap is the identity -- it forwards rather
// than translates -- and the only place the name reads oddly is the unit suite,
// where a hand-built T-junction candidate comes back NotConverged having
// converged nothing.
//
//   * A refusal of guarantee 14, EITHER CLAUSE  -> NotConverged. Snap rounding
//     can genuinely create a crossing that was not in the input, so this is a
//     diagnosis about the data, and the driver's answer below the cap is another
//     round. The lever is the spacing.
//   * A refusal of 11, 12, 13 or 16 -> MalformedOutput. Those the driver
//     establishes by construction, so a refusal means the driver is broken. A
//     self-check with a name, exactly parallel to CdtStatus::MalformedInput.
//   * A candidate violating BOTH reports MalformedOutput: the self-check
//     outranks the diagnosis, because a builder handed malformed arrays has no
//     basis for saying anything about convergence.
//
// WHY GUARANTEE 15 IS NOT CHECKED HERE. Its oracle must be built from the INPUT
// Pslg, and the builder does not have it. Re-deriving 15 from the arrays the
// driver just handed over is the dedup restating itself. 15 is checked in
// prop_noding_no_crossings.cpp, against the input, and nowhere else -- so a
// NodedPslg promises 11-14 and 16 by construction and 15 only as far as that
// suite reaches. Stated rather than papered over.
//
// What the builder CAN see of 15 is its array's LENGTH, which is what makes the
// array dense, and a short one is a read past the end at every caller's
// edge_base(c) + k. That is checked, as part of 16's structural half.
//
// THE BUILDER IS NOT THE PSLG VALIDATOR. It does not ask for an outer chain, it
// does not check ring winding and it does not re-run increment 3's stages: its
// input has already been through them. Ring re-validation after snapping is
// node<K>'s, because the statuses it raises (RingCollapsed, NonSimpleRing) are
// diagnoses about the spacing rather than claims about these arrays.

#include <terrain/core/bbox.hpp>
#include <terrain/core/edge_properties.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <format>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace terrain::noding {

enum class NodeStatus : int {
    Ok = 0,
    NotRun,
    InvalidSnapSpacing,
    CoordinateOutOfRange,
    RingCollapsed,
    RingDegenerateAfterSnap,
    NonSimpleRing,
    NotConverged,
    MalformedOutput,
};

// Total over the enumeration, one distinct sentence each. The lever the reader
// should reach for is part of the text, because three of these rows point at the
// spacing and two of them point in opposite directions.
[[nodiscard]] constexpr std::string_view describe(NodeStatus s) noexcept {
    switch (s) {
        case NodeStatus::Ok:
            return "the constraint set was noded";
        case NodeStatus::NotRun:
            return "the noder has not been run";
        case NodeStatus::InvalidSnapSpacing:
            return "the snap spacing must be finite and strictly positive";
        case NodeStatus::CoordinateOutOfRange:
            return "a coordinate is too large for the snap grid; coarsen the "
                   "spacing or re-project";
        case NodeStatus::RingCollapsed:
            return "a ring collapsed to fewer than three nodes; use a finer spacing";
        case NodeStatus::RingDegenerateAfterSnap:
            return "a ring's winding did not survive snapping; self-check, no "
                   "input is known to reach it";
        case NodeStatus::NonSimpleRing:
            return "a ring visits the same node twice; repair the polygon upstream";
        case NodeStatus::NotConverged:
            return "the constraint set did not settle within the round cap; "
                   "use a finer spacing";
        case NodeStatus::MalformedOutput:
            return "the noder produced output violating its own guarantees; this "
                   "is our bug";
    }
    return "unknown noder status";
}

struct NodeOutcome {
    std::optional<NodedPslg> pslg;  // engaged iff status == Ok
    NodeStatus status{NodeStatus::NotRun};
    std::string message;

    [[nodiscard]] bool ok() const noexcept { return pslg.has_value(); }
};

namespace detail {

[[nodiscard]] inline NodeOutcome refuse(NodeStatus status, std::string message) {
    return NodeOutcome{std::nullopt, status, std::move(message)};
}

// The "nothing to report" outcome of an internal check: no pslg, no message.
// Not a success -- only build() constructs one of those.
[[nodiscard]] inline NodeOutcome pass() {
    return NodeOutcome{std::nullopt, NodeStatus::Ok, std::string{}};
}

// The two verdicts of the builder, as named functions, because which one a
// clause maps to is the whole content of "how to read NodeStatus" above and is
// not a thing to re-derive at twelve call sites.
[[nodiscard]] inline NodeOutcome malformed(std::string message) {
    return refuse(NodeStatus::MalformedOutput, std::move(message));
}

[[nodiscard]] inline NodeOutcome not_noded(std::string message) {
    return refuse(NodeStatus::NotConverged, std::move(message));
}

}  // namespace detail

// The pair search behind guarantee 14: every pair of CLOSED boxes that
// intersect, found by sort and sweep on x (increment 16b-0, R8). Kernel-free:
// it compares box bounds and nothing else, and what a pair MEANS is the
// callbacks' business.
//
// on_edge_pair(i, j) is called exactly once for every unordered pair of edge
// boxes that meet (i != j, either order); on_edge_cell(e, c) exactly once for
// every edge box and cell box that meet. Two cells are never paired. A callback
// returning false stops the sweep at once, and the function returns false;
// otherwise it returns true. The same input gives the same calls in the same
// order.
//
// The sweep: items sorted by (xmin, kind, index), a total order, so the sort is
// deterministic. Each arriving item is tested against the ACTIVE items whose
// xmax >= its xmin (closed), and those whose xmax is below it are dropped on the
// way past, for good: every later item starts at or right of this one. A pair is
// therefore reported once, when its later-sorted member arrives, and only if
// the earlier one still spans that x -- which, with the y-test, is exactly the
// closed-box intersection.
//
// COST. Typical inputs are O(n log n + k) for k reported pairs plus the
// x-overlapping pairs whose y-intervals miss. The worst case is still
// quadratic: many long east-west boxes all overlap in x. It is the pair search
// only; it shares no code with the driver's BroadPhase, and must not (see
// check_guarantee_14).
template <class OnEdgePair, class OnEdgeCell>
[[nodiscard]] bool sweep_box_pairs(std::span<const Box2> edges, std::span<const Box2> cells,
                                   OnEdgePair&& on_edge_pair, OnEdgeCell&& on_edge_cell) {
    struct Item {
        double xmin;
        bool is_cell;  // edges before cells on an x tie; either order is complete
        std::size_t index;
    };
    std::vector<Item> items;
    items.reserve(edges.size() + cells.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
        items.push_back(Item{edges[i].lo().x, false, i});
    }
    for (std::size_t c = 0; c < cells.size(); ++c) {
        items.push_back(Item{cells[c].lo().x, true, c});
    }
    std::ranges::sort(items, [](const Item& a, const Item& b) {
        return std::tie(a.xmin, a.is_cell, a.index) < std::tie(b.xmin, b.is_cell, b.index);
    });

    // Reports the members of `active` meeting `b`, dropping (swap-and-pop) those
    // that end left of it. False as soon as `report` returns false.
    const auto scan = [](std::vector<std::size_t>& active, std::span<const Box2> boxes,
                         const Box2& b, auto&& report) {
        for (std::size_t k = 0; k < active.size();) {
            const Box2& a = boxes[active[k]];
            if (a.hi().x < b.lo().x) {
                active[k] = active.back();
                active.pop_back();
                continue;
            }
            if (a.lo().y <= b.hi().y && b.lo().y <= a.hi().y && !report(active[k])) {
                return false;
            }
            ++k;
        }
        return true;
    };

    std::vector<std::size_t> active_edges;
    std::vector<std::size_t> active_cells;
    for (const Item& t : items) {
        if (t.is_cell) {
            const auto meet = [&](std::size_t e) { return on_edge_cell(e, t.index); };
            if (!scan(active_edges, edges, cells[t.index], meet)) {
                return false;
            }
            active_cells.push_back(t.index);
            continue;
        }
        const auto pair = [&](std::size_t e) { return on_edge_pair(e, t.index); };
        const auto meet = [&](std::size_t c) { return on_edge_cell(t.index, c); };
        if (!scan(active_edges, edges, edges[t.index], pair) ||
            !scan(active_cells, cells, edges[t.index], meet)) {
            return false;
        }
        active_edges.push_back(t.index);
    }
    return true;
}

class NodedPslgBuilder {
public:
    // Takes the candidate arrays by value and consumes them. The kernel is a
    // parameter of build(), not of the class, exactly as PslgBuilder::build is.
    NodedPslgBuilder(SnapGrid grid, std::vector<Point2> vertices, std::vector<Chain> chains,
                     std::vector<std::uint32_t> chain_indices,
                     std::vector<EdgeProperties> edge_properties,
                     std::vector<std::uint32_t> node_of_input_vertex)
        : grid_{grid},
          vertices_{std::move(vertices)},
          chains_{std::move(chains)},
          chain_indices_{std::move(chain_indices)},
          edge_properties_{std::move(edge_properties)},
          node_of_input_vertex_{std::move(node_of_input_vertex)} {}

    // Rvalue-ref-qualified, the same shape as PslgBuilder::build: the arrays are
    // moved into the NodedPslg, so "build twice and get two objects sharing
    // nothing" is a compile error rather than a subtle question about what the
    // second one sees.
    //
    // CHECK ORDER IS NOT ARBITRARY. The structural checks come first because
    // every later one indexes through them; 11 before 12 because 12's node-id
    // reasoning is snap()'s and snap() has can_snap as a precondition; 14 last
    // because it is the expensive one and because its verdict is outranked by
    // any of the others.
    template <pred::GeometryKernel K>
    [[nodiscard]] NodeOutcome build() && {
        if (NodeOutcome bad = check_structure(); bad.status != NodeStatus::Ok) {
            return bad;
        }
        if (NodeOutcome bad = check_grid(); bad.status != NodeStatus::Ok) {
            return bad;
        }
        if (NodeOutcome bad = check_distinct_and_nondegenerate(); bad.status != NodeStatus::Ok) {
            return bad;
        }
        if (NodeOutcome bad = check_guarantee_14<K>(); bad.status != NodeStatus::Ok) {
            return bad;
        }

        NodedPslg out;
        out.grid_ = grid_;
        out.vertices_ = std::move(vertices_);
        out.chains_ = std::move(chains_);
        out.chain_indices_ = std::move(chain_indices_);
        out.edge_properties_ = std::move(edge_properties_);
        out.node_of_input_vertex_ = std::move(node_of_input_vertex_);
        out.edge_base_ = std::move(edge_base_);
        return NodeOutcome{std::move(out), NodeStatus::Ok, std::string{}};
    }

private:
    // Guarantee 16's structural half plus guarantee 15's length. Everything the
    // later checks index through is established here, which is why it is first:
    // a chain slice running past the flat buffer must be a status and not a read
    // past the end.
    [[nodiscard]] NodeOutcome check_structure() {
        edge_base_.assign(chains_.size() + 1, 0);
        std::size_t running = 0;
        std::size_t edges = 0;

        for (std::size_t c = 0; c < chains_.size(); ++c) {
            const Chain& ch = chains_[c];
            const auto begin = static_cast<std::size_t>(ch.begin);
            const auto count = static_cast<std::size_t>(ch.count);
            if (begin != running || count == 0 || begin + count > chain_indices_.size()) {
                return detail::malformed(std::format("chain {} does not tile the buffer", c));
            }
            if (is_closed(ch.role) && count < 3) {
                return detail::malformed(std::format("closed chain {} has under 3 nodes", c));
            }
            running = begin + count;
            edge_base_[c] = edges;
            edges += is_closed(ch.role) ? count : count - 1;
        }
        edge_base_[chains_.size()] = edges;

        if (running != chain_indices_.size()) {
            return detail::malformed("the chain slices leave a gap in the index buffer");
        }
        for (const std::uint32_t i : chain_indices_) {
            if (i >= vertices_.size()) {
                return detail::malformed(std::format("chain index {} is out of range", i));
            }
        }
        for (const std::uint32_t i : node_of_input_vertex_) {
            if (i >= vertices_.size()) {
                return detail::malformed(std::format("node id {} is out of range", i));
            }
        }
        if (edge_properties_.size() != edges) {
            return detail::malformed(std::format("edge_properties has {} entries, not {}",
                                                 edge_properties_.size(), edges));
        }
        return detail::pass();
    }

    // Guarantee 11. can_snap is consulted BEFORE snapped(), and the order is
    // load-bearing rather than tidy: snap() debug-asserts can_snap and
    // cell_min/cell_max debug-assert the index range, so the reversed spelling
    // aborts in a Debug build and silently wraps an index in Release. Both are
    // worse than a status.
    [[nodiscard]] NodeOutcome check_grid() {
        for (std::size_t i = 0; i < vertices_.size(); ++i) {
            const Point2& v = vertices_[i];
            if (!grid_.can_snap(v)) {
                return detail::malformed(std::format("vertex {} is out of grid range", i));
            }
            const Point2 again = grid_.snapped(v);
            if (again.x != v.x || again.y != v.y) {
                return detail::malformed(std::format("vertex {} is not on the grid", i));
            }
        }
        return detail::pass();
    }

    // Guarantees 12 and 13. 12 is decided on grid indices rather than on doubles:
    // after 11 the two are the same question, and std::int64_t has no NaN, so the
    // sort is total and the equality exact.
    [[nodiscard]] NodeOutcome check_distinct_and_nondegenerate() {
        std::vector<GridPoint> keys;
        keys.reserve(vertices_.size());
        for (const Point2& v : vertices_) {
            keys.push_back(grid_.snap(v));
        }
        std::ranges::sort(keys);
        if (std::ranges::adjacent_find(keys) != keys.end()) {
            return detail::malformed("two vertices share a grid point");
        }

        for (std::size_t c = 0; c < chains_.size(); ++c) {
            const std::span<const std::uint32_t> idx = indices_of(c);
            const std::size_t n = idx.size();
            for (std::size_t k = 0; k < edge_count(c); ++k) {
                if (idx[k] == idx[(k + 1) % n]) {
                    return detail::malformed(std::format("chain {} edge {} is zero-length", c, k));
                }
            }
        }
        return detail::pass();
    }

    struct Edge {
        Segment2 segment;
        std::uint32_t lo, hi;  // the UNDIRECTED node-id pair, which is 14(a)'s key
    };

    // Guarantee 14, over the (edge, edge) and (edge, node) pairs whose closed
    // boxes meet, found by sweep_box_pairs. Complete: two segments that cross or
    // overlap share a point, so their boxes meet, and a segment meeting a closed
    // cell meets any box containing the cell. Box bounds are min/max of the
    // endpoints' doubles, so exact. NOT the driver's BroadPhase: the index is the
    // driver's, and a verification pass sharing it would be blind wherever the
    // index is. The sweep is a different algorithm and shares no code with it.
    // The tests on each pair are the brute force's, unchanged; only which
    // violating pair the message names first may differ from it.
    template <pred::GeometryKernel K>
    [[nodiscard]] NodeOutcome check_guarantee_14() {
        std::vector<Edge> edges;
        std::vector<Box2> edge_boxes;
        edges.reserve(edge_base_.back());
        edge_boxes.reserve(edge_base_.back());
        for (std::size_t c = 0; c < chains_.size(); ++c) {
            const std::span<const std::uint32_t> idx = indices_of(c);
            const std::size_t n = idx.size();
            for (std::size_t k = 0; k < edge_count(c); ++k) {
                const std::uint32_t u = idx[k];
                const std::uint32_t v = idx[(k + 1) % n];
                edges.push_back(Edge{Segment2{vertices_[u], vertices_[v]},
                                     std::min(u, v), std::max(u, v)});
                edge_boxes.push_back(bbox(edges.back().segment));
            }
        }

        // Node g's closed cell, padded by one spacing on every side: a false
        // positive only costs a segment_meets_cell call, and the padding leaves
        // no question of how the cell's bounds round. Clamped to the finite
        // range, where Box2 lives; the edges are finite (guarantee 11), so the
        // clamp changes no answer.
        constexpr double kBig = std::numeric_limits<double>::max();
        const double pad = grid_.spacing();
        std::vector<Box2> cell_boxes;
        cell_boxes.reserve(vertices_.size());
        for (const Point2& v : vertices_) {
            const GridPoint g = grid_.snap(v);
            const Point2 lo = grid_.cell_min(g);
            const Point2 hi = grid_.cell_max(g);
            cell_boxes.emplace_back(
                Point2{std::max(lo.x - pad, -kBig), std::max(lo.y - pad, -kBig)},
                Point2{std::min(hi.x + pad, kBig), std::min(hi.y + pad, kBig)});
        }

        std::string why;
        // 14(a), amended: Disjoint or Touching, OR Overlapping with EQUAL node-id
        // pairs. Duplicate edges are legal -- a road noded along a river yields
        // two chains carrying the same edge, and classify on two identical
        // segments returns Overlapping. A PARTIAL overlap is still a violation,
        // and the clause is stated over node ids rather than over collinearity so
        // that the decision stays on integers.
        const auto edge_pair = [&](std::size_t i, std::size_t j) {
            const Edge& a = edges[i];
            const Edge& b = edges[j];
            const SegmentRelation r = classify<K>(a.segment, b.segment);
            const bool same_pair = a.lo == b.lo && a.hi == b.hi;
            if (r == SegmentRelation::Crossing ||
                (r == SegmentRelation::Overlapping && !same_pair)) {
                why = std::format("edges ({}, {}) and ({}, {}) meet off a node", a.lo, a.hi,
                                  b.lo, b.hi);
                return false;
            }
            return true;
        };
        // 14(b), the hot-pixel clause. segment_meets_cell, NEVER on_segment:
        // world(g) is not affine at a non-dyadic spacing, so a snapped node is
        // NEAR a segment and not ON it, and an exact-incidence spelling here
        // blesses the T-junctions the split pass missed instead of catching them.
        const auto edge_cell = [&](std::size_t i, std::size_t g) {
            const Edge& e = edges[i];
            if (g == e.lo || g == e.hi ||
                !segment_meets_cell<K>(grid_, e.segment, grid_.snap(vertices_[g]))) {
                return true;
            }
            why = std::format("node {}'s cell meets edge ({}, {}), which it does not end", g,
                              e.lo, e.hi);
            return false;
        };
        if (!sweep_box_pairs(edge_boxes, cell_boxes, edge_pair, edge_cell)) {
            return detail::not_noded(std::move(why));
        }
        return detail::pass();
    }

    [[nodiscard]] std::span<const std::uint32_t> indices_of(std::size_t c) const noexcept {
        const Chain& ch = chains_[c];
        return std::span<const std::uint32_t>{chain_indices_}.subspan(ch.begin, ch.count);
    }

    [[nodiscard]] std::size_t edge_count(std::size_t c) const noexcept {
        const Chain& ch = chains_[c];
        const auto n = static_cast<std::size_t>(ch.count);
        return is_closed(ch.role) ? n : n - 1;
    }

    SnapGrid grid_;
    std::vector<Point2> vertices_;
    std::vector<Chain> chains_;
    std::vector<std::uint32_t> chain_indices_;
    std::vector<EdgeProperties> edge_properties_;
    std::vector<std::uint32_t> node_of_input_vertex_;
    std::vector<std::size_t> edge_base_;
};

}  // namespace terrain::noding
