#pragma once

// THE FROZEN ORACLE for increment 16b-0 (docs/increments/16b-terrain-polygons.md,
// R8 and I6): "the brute force stays, as the test oracle, in
// tests/cpp/support/".
//
// brute_force_guarantee_14 is NodedPslgBuilder::check_guarantee_14 from
// include/terrain/noding/noded_pslg_builder.hpp COPIED at 6b4fcb9 (the 16b
// rulings commit), before any production change on the branch. The edits are
// the ones a free function needs: the builder's members become parameters,
// indices_of and edge_count are inlined, and the result is the STATUS alone
// (Ok or NotConverged), because R8 lets the sweep name a different violating
// pair first and no test asserts wording. Every (edge, edge) and every
// (edge, node) pair is tested with the same classify<K> and segment_meets_cell<K>
// the builder uses: the producer's relation, never exact incidence.
//
// brute_force_box_pairs is the other half of the oracle: the pairs of CLOSED
// boxes that intersect, by testing every pair. It is the specification
// sweep_box_pairs is held to.
//
// FROZEN: no test may change this file to agree with new code.
//
// Test-only. Nothing under include/ may include it.

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/noding/noded_pslg_builder.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <set>
#include <span>
#include <utility>
#include <vector>

namespace verifier_oracle {

using terrain::Box2;
using terrain::Chain;
using terrain::Point2;
using terrain::Segment2;
using terrain::SnapGrid;
using terrain::noding::NodeStatus;
using terrain::noding::SegmentRelation;

// Guarantee 14 by brute force, over arrays that already satisfy 11, 12, 13 and
// 16 (the caller's generator establishes them; the builder checks them first).
template <terrain::pred::GeometryKernel K>
[[nodiscard]] NodeStatus brute_force_guarantee_14(const SnapGrid& grid,
                                                  std::span<const Point2> vertices,
                                                  std::span<const Chain> chains,
                                                  std::span<const std::uint32_t> chain_indices) {
    struct Edge {
        Segment2 segment;
        std::uint32_t lo, hi;
    };
    std::vector<Edge> edges;
    for (const Chain& ch : chains) {
        const auto idx = chain_indices.subspan(ch.begin, ch.count);
        const std::size_t n = idx.size();
        const std::size_t count = terrain::is_closed(ch.role) ? n : n - 1;
        for (std::size_t k = 0; k < count; ++k) {
            const std::uint32_t u = idx[k];
            const std::uint32_t v = idx[(k + 1) % n];
            edges.push_back(
                Edge{Segment2{vertices[u], vertices[v]}, std::min(u, v), std::max(u, v)});
        }
    }

    for (std::size_t i = 0; i < edges.size(); ++i) {
        for (std::size_t j = i + 1; j < edges.size(); ++j) {
            const SegmentRelation r =
                terrain::noding::classify<K>(edges[i].segment, edges[j].segment);
            const bool same_pair = edges[i].lo == edges[j].lo && edges[i].hi == edges[j].hi;
            if (r == SegmentRelation::Crossing ||
                (r == SegmentRelation::Overlapping && !same_pair)) {
                return NodeStatus::NotConverged;
            }
        }
    }

    for (const Edge& e : edges) {
        for (std::uint32_t g = 0; g < vertices.size(); ++g) {
            if (g == e.lo || g == e.hi) {
                continue;
            }
            if (terrain::noding::segment_meets_cell<K>(grid, e.segment,
                                                        grid.snap(vertices[g]))) {
                return NodeStatus::NotConverged;
            }
        }
    }
    return NodeStatus::Ok;
}

[[nodiscard]] inline bool closed_boxes_meet(const Box2& a, const Box2& b) noexcept {
    return a.lo().x <= b.hi().x && b.lo().x <= a.hi().x && a.lo().y <= b.hi().y &&
           b.lo().y <= a.hi().y;
}

// The unordered edge pairs are stored as (min, max).
struct BoxPairs {
    std::set<std::pair<std::size_t, std::size_t>> edge_edge;
    std::set<std::pair<std::size_t, std::size_t>> edge_cell;  // (edge, cell)

    friend bool operator==(const BoxPairs&, const BoxPairs&) = default;
};

[[nodiscard]] inline BoxPairs brute_force_box_pairs(std::span<const Box2> edges,
                                                    std::span<const Box2> cells) {
    BoxPairs out;
    for (std::size_t i = 0; i < edges.size(); ++i) {
        for (std::size_t j = i + 1; j < edges.size(); ++j) {
            if (closed_boxes_meet(edges[i], edges[j])) {
                out.edge_edge.emplace(i, j);
            }
        }
        for (std::size_t c = 0; c < cells.size(); ++c) {
            if (closed_boxes_meet(edges[i], cells[c])) {
                out.edge_cell.emplace(i, c);
            }
        }
    }
    return out;
}

}  // namespace verifier_oracle
