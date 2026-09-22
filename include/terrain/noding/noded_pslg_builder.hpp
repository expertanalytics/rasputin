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
#include <optional>
#include <span>
#include <string>
#include <string_view>
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

    // Guarantee 14, brute force over every (edge, edge) and every (node, edge)
    // pair. NO BROAD PHASE: the index is the driver's, and a verification pass
    // sharing it would be blind wherever the index is. The candidates are small
    // by construction -- a noded constraint set, not a mesh.
    template <pred::GeometryKernel K>
    [[nodiscard]] NodeOutcome check_guarantee_14() {
        std::vector<Edge> edges;
        edges.reserve(edge_base_.back());
        for (std::size_t c = 0; c < chains_.size(); ++c) {
            const std::span<const std::uint32_t> idx = indices_of(c);
            const std::size_t n = idx.size();
            for (std::size_t k = 0; k < edge_count(c); ++k) {
                const std::uint32_t u = idx[k];
                const std::uint32_t v = idx[(k + 1) % n];
                edges.push_back(Edge{Segment2{vertices_[u], vertices_[v]},
                                     std::min(u, v), std::max(u, v)});
            }
        }

        // 14(a), amended: Disjoint or Touching, OR Overlapping with EQUAL node-id
        // pairs. Duplicate edges are legal -- a road noded along a river yields
        // two chains carrying the same edge, and classify on two identical
        // segments returns Overlapping. A PARTIAL overlap is still a violation,
        // and the clause is stated over node ids rather than over collinearity so
        // that the decision stays on integers.
        for (std::size_t i = 0; i < edges.size(); ++i) {
            for (std::size_t j = i + 1; j < edges.size(); ++j) {
                const SegmentRelation r = classify<K>(edges[i].segment, edges[j].segment);
                const bool same_pair = edges[i].lo == edges[j].lo && edges[i].hi == edges[j].hi;
                if (r == SegmentRelation::Crossing ||
                    (r == SegmentRelation::Overlapping && !same_pair)) {
                    return detail::not_noded(
                        std::format("edges ({}, {}) and ({}, {}) meet off a node", edges[i].lo,
                                    edges[i].hi, edges[j].lo, edges[j].hi));
                }
            }
        }

        // 14(b), the hot-pixel clause. segment_meets_cell, NEVER on_segment:
        // world(g) is not affine at a non-dyadic spacing, so a snapped node is
        // NEAR a segment and not ON it, and an exact-incidence spelling here
        // blesses the T-junctions the split pass missed instead of catching them.
        for (const Edge& e : edges) {
            for (std::uint32_t g = 0; g < vertices_.size(); ++g) {
                if (g == e.lo || g == e.hi) {
                    continue;
                }
                if (segment_meets_cell<K>(grid_, e.segment, grid_.snap(vertices_[g]))) {
                    return detail::not_noded(std::format(
                        "node {}'s cell meets edge ({}, {}), which it does not end", g, e.lo,
                        e.hi));
                }
            }
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
