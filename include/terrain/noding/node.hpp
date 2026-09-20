#pragma once

// The noder's driver: a valid Pslg in, a verified NodedPslg or a status out.
//
// node<K> IS A PURE FUNCTION of (pslg, options). No statics, no caches, no
// global state, no allocation outside its own locals. Two calls on the same
// input from two threads produce bit-identical output, because every ordering in
// it is either GridPoint's total order or the input's chain order, and never the
// order a thread happened to find something in. That is what makes it safe to
// release the GIL around at 5c.
//
// THE SIX STEPS, and what each one is allowed to decide:
//
//   0. Admission. is_valid_spacing -> InvalidSnapSpacing; one pass over
//      pslg.vertices() with SnapGrid::can_snap -> CoordinateOutOfRange. This is
//      the one place kMaxGridIndex is compared against, and after it snap() is
//      branch-free.
//   1. Broad phase over the current segment set, each segment's bbox the query.
//   2. classify<K> per candidate pair. Crossing -> crossing_point<K> contributes
//      a node. Touching and Overlapping contribute nothing new. THE ONLY STEP
//      THAT MANUFACTURES A COORDINATE.
//   3. NodeSet over every snapped input vertex -- referenced or not, which is
//      what makes node_of_input_vertex total -- every current chain position and
//      every constructed crossing. Node ids are its sorted order and are
//      therefore a function of the SET.
//   4. Broad phase again, segments against nodes, queried by the node's
//      cell_min/cell_max BOX and never by world(g). Each segment is split at
//      every node whose closed cell it meets, in arc order along the segment,
//      ties broken by GridPoint's <=>.
//   5. Reassemble chains and merge. The chain structure is preserved; the dedup
//      is over the edge SET and what it produces is the per-edge property array.
//   6. Verify by handing the arrays to NodedPslgBuilder::build<K>().
//
// "NO NEW NODE" IS NOT THE FIXPOINT AND IS NOT THE LOOP CONDITION. A round can
// split an edge at an EXISTING node -- a new edge, no new node -- so a loop
// spelled that way reports converged in a state that violates guarantee 14(b).
// The loop ends when the VERIFICATION passes, and at the cap it ends with
// NotConverged.
//
// WHAT NotConverged DOES NOT ALWAYS MEAN. Its printed lever is a finer spacing,
// and that cures the ordinary case -- snap rounding creating a crossing that was
// not in the input. It does not reliably cure a constraint edge passing exactly
// through a lattice corner with nodes on both flanking cells: refining moves the
// lattice, so whether those cells still hold nodes is re-rolled rather than
// fixed. segment_meets_cell uses CLOSED cells, so such an edge must be split at
// both flanking nodes, and the resulting middle edge grazes the cells of the
// original endpoints -- the orbit has period two and Ok is unreachable at any
// cap. That is a known, named defect, out of scope here and specified in
// docs/increments/05d-corner-graze.md, whose fix changes 5a's shipped predicate.
// 5b ships NotConverged for it deliberately; a driver that answers Ok there has
// either skipped the verification or is running a verifier that agrees with it.
//
// Header-only and templated on K, exactly as increment 3's validator is; nothing
// in src/noding/.

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/noding/broad_phase.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/noding/node_set.hpp>
#include <terrain/noding/noded_pslg_builder.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <format>
#include <map>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace terrain::noding {

struct NodeOptions {
    // NO DEFAULT, EVER, and first so that NodeOptions{0.1} is the only cheap
    // spelling. There is no default spacing anywhere in C++, for the same reason
    // K has no default in ring.hpp and pslg_builder.hpp: the right value is a
    // policy question about the data, and three of the failure statuses point at
    // it, two of them in opposite directions. A default in cli.py at 5c is the
    // composition root making a policy choice and is a different thing.
    double spacing;

    // This one DOES have a default, because it is a cap on a loop rather than a
    // parameter of the answer: at any value the outcome is either the same mesh
    // or NotConverged, never a different mesh. Four rather than one because a
    // snap-induced cascade is real, and four rather than sixteen because it was
    // measured to be rarer than the theory allows.
    std::uint32_t max_rounds{4};
};

namespace detail {

// One segment of the current working set, with the chain and edge it belongs to.
struct OwnedSegment {
    Segment2 segment;
    std::size_t chain;
    std::size_t edge;
};

// Step 4's sort key. The dot product of (world(g) - s.a) with (s.b - s.a),
// which orders the split points ALONG THE SEGMENT rather than by node id -- the
// two disagree whenever the segment's ix decreases, and a driver that sorted by
// id would emit the chain reversed through its own middle. Exact ties -- an edge
// through a lattice corner, where two flanking nodes project identically -- are
// broken by GridPoint's defaulted, therefore lexicographic, operator<=>, which
// is a total order and so makes the result a function of the node set alone.
struct ArcKey {
    double t{};
    GridPoint g{};
    std::uint32_t id{};

    // A free comparator rather than operator<=>: ranges::less demands
    // totally_ordered, and a defaulted operator== over all three members would
    // not agree with an ordering that deliberately ignores the id.
    [[nodiscard]] static bool before(const ArcKey& a, const ArcKey& b) noexcept {
        return a.t != b.t ? a.t < b.t : a.g < b.g;
    }
};

[[nodiscard]] inline double arc_of(const Segment2& s, const Point2& w) noexcept {
    return (w.x - s.a.x) * (s.b.x - s.a.x) + (w.y - s.a.y) * (s.b.y - s.a.y);
}

[[nodiscard]] inline std::size_t edge_count_of(ChainRole role, std::size_t n) noexcept {
    if (is_closed(role)) {
        return n;
    }
    return n > 0 ? n - 1 : 0;
}

// The dedup key: the UNDIRECTED node-id pair. Integers, which is what keeps the
// merge independent of which way a snap rounded.
[[nodiscard]] inline std::pair<std::uint32_t, std::uint32_t> edge_key(
    const std::vector<std::uint32_t>& run, std::size_t k) noexcept {
    const std::uint32_t u = run[k];
    const std::uint32_t v = run[(k + 1) % run.size()];
    return {std::min(u, v), std::max(u, v)};
}

// Ring re-validation, between steps 5 and 6, and the three questions are asked
// in the order the design fixes: distinct node count, then winding, then
// simplicity.
//
// The simplicity test is a SORT OVER std::uint32_t and not a quadratic
// segment-pair scan, and the equivalence that licenses that is worth stating: on
// output satisfying 14(a), a closed chain self-intersects IFF it visits the same
// node twice, because 14(a) has already excluded interior crossings that are not
// at nodes. The obvious implementation is both the slower one and the wrong one
// -- it misses a ring that revisits a node without any edge pair crossing.
//
// RingDegenerateAfterSnap IS A SELF-CHECK WITH A NAME, like MalformedOutput and
// CdtStatus::NotNoded, and no fixture is owed for it: a winding flip needs the
// snapped signed area to cross zero, which puts a vertex within half a cell of
// the edge opposite it -- exactly segment_meets_cell -- so 14(b) requires that
// edge to be split there, the split repeats a node id, and NonSimpleRing fires
// first. Measured by a bounded exhaustive search over 4-gons and 5-gons: every
// winding flip grazes, no block has a survivor. The check stays because deleting
// a cheap check whose job is to be unreachable is how a later change makes it
// reachable in silence.
template <pred::GeometryKernel K>
[[nodiscard]] NodeOutcome revalidate_rings(const Pslg& in, double spacing,
                                           std::span<const Point2> vertices,
                                           const std::vector<std::vector<std::uint32_t>>& runs) {
    const auto distinct_count = [](const std::vector<std::uint32_t>& run) {
        std::vector<std::uint32_t> sorted = run;
        std::ranges::sort(sorted);
        const auto dup = std::ranges::unique(sorted);
        return static_cast<std::size_t>(std::distance(sorted.begin(), dup.begin()));
    };

    // The message carries the EXTENT, which is what a diagnostics vector would
    // have carried and the only thing a single status loses. There is no vector,
    // because every failure here is snap-induced and they share one lever.
    std::size_t collapsed = 0;
    std::size_t first_collapsed = 0;
    for (std::size_t c = 0; c < runs.size(); ++c) {
        if (is_closed(in.chains()[c].role) && distinct_count(runs[c]) < 3) {
            if (collapsed++ == 0) {
                first_collapsed = c;
            }
        }
    }
    if (collapsed != 0) {
        return refuse(NodeStatus::RingCollapsed,
                      std::format("{} ring(s) collapsed to under 3 nodes at spacing {}; first is "
                                  "chain {}",
                                  collapsed, spacing, first_collapsed));
    }

    for (std::size_t c = 0; c < runs.size(); ++c) {
        const ChainRole role = in.chains()[c].role;
        if (!is_closed(role)) {
            continue;
        }
        const pred::Orientation want = role == ChainRole::Outer
                                           ? pred::Orientation::CounterClockwise
                                           : pred::Orientation::Clockwise;
        // Re-evaluated on the SNAPPED ring, and NEVER silently reversed.
        const IndexedRing ring{vertices, std::span<const std::uint32_t>{runs[c]}};
        if (orientation<K>(ring) != want) {
            return refuse(NodeStatus::RingDegenerateAfterSnap,
                          std::format("chain {}'s winding did not survive spacing {}", c, spacing));
        }
    }

    for (std::size_t c = 0; c < runs.size(); ++c) {
        if (!is_closed(in.chains()[c].role)) {
            continue;
        }
        if (distinct_count(runs[c]) != runs[c].size()) {
            return refuse(NodeStatus::NonSimpleRing,
                          std::format("chain {} visits a node twice at spacing {}", c, spacing));
        }
    }
    return pass();
}

}  // namespace detail

// Precondition: none beyond pslg being a valid Pslg, which its type asserts.
template <pred::GeometryKernel K>
[[nodiscard]] NodeOutcome node(const Pslg& pslg, const NodeOptions& options) {
    // --- Step 0: admission -------------------------------------------------
    if (!is_valid_spacing(options.spacing)) {
        return detail::refuse(
            NodeStatus::InvalidSnapSpacing,
            std::format("snap spacing {} is not finite and strictly positive", options.spacing));
    }
    const SnapGrid grid{options.spacing};

    std::vector<GridPoint> input_keys;
    input_keys.reserve(pslg.vertices().size());
    for (std::size_t i = 0; i < pslg.vertices().size(); ++i) {
        if (!grid.can_snap(pslg.vertices()[i])) {
            return detail::refuse(NodeStatus::CoordinateOutOfRange,
                                  std::format("input vertex {} cannot be snapped at spacing {}", i,
                                              options.spacing));
        }
        input_keys.push_back(grid.snap(pslg.vertices()[i]));
    }

    // The working geometry. Round 1 runs on the INPUT coordinates, because a
    // crossing of the input is what the user asked about; from round 2 on it runs
    // on world(node id), which is where snap-induced crossings live.
    std::vector<std::vector<Point2>> pts(pslg.chains().size());
    std::vector<std::vector<EdgeProperties>> props(pslg.chains().size());
    for (std::size_t c = 0; c < pslg.chains().size(); ++c) {
        for (const std::uint32_t i : pslg.indices_of(c)) {
            pts[c].push_back(pslg.vertices()[i]);
        }
        props[c].assign(pslg.edge_count(c), pslg.chains()[c].properties);
    }

    NodeOutcome last;
    for (std::uint32_t round = 0; round < options.max_rounds; ++round) {
        // --- Step 1: the segment set and its index -------------------------
        std::vector<detail::OwnedSegment> segs;
        std::vector<Segment2> raw;
        for (std::size_t c = 0; c < pts.size(); ++c) {
            const std::size_t n = pts[c].size();
            const std::size_t ec = detail::edge_count_of(pslg.chains()[c].role, n);
            for (std::size_t k = 0; k < ec; ++k) {
                segs.push_back(detail::OwnedSegment{Segment2{pts[c][k], pts[c][(k + 1) % n]}, c, k});
                raw.push_back(segs.back().segment);
            }
        }
        const BroadPhase index{raw};

        // --- Step 2: classify, and construct exactly at a crossing ---------
        std::vector<GridPoint> keys = input_keys;
        for (const std::vector<Point2>& chain : pts) {
            for (const Point2& p : chain) {
                keys.push_back(grid.snap(p));
            }
        }
        for (std::size_t i = 0; i < raw.size(); ++i) {
            index.for_each_candidate(bbox(raw[i]), [&](std::uint32_t j) {
                if (j <= i) {
                    return;
                }
                if (classify<K>(raw[i], raw[j]) == SegmentRelation::Crossing) {
                    keys.push_back(crossing_point<K>(grid, raw[i], raw[j]));
                }
            });
        }

        // --- Step 3: the node table ----------------------------------------
        const NodeSet nodes{std::move(keys)};

        // --- Step 4: segments against nodes, queried by the CELL BOX -------
        std::vector<std::vector<detail::ArcKey>> met(raw.size());
        for (std::uint32_t n = 0; n < nodes.size(); ++n) {
            const GridPoint& g = nodes[n];
            const Box2 cell{grid.cell_min(g), grid.cell_max(g)};
            index.for_each_candidate(cell, [&](std::uint32_t i) {
                if (segment_meets_cell<K>(grid, raw[i], g)) {
                    met[i].push_back(detail::ArcKey{detail::arc_of(raw[i], grid.world(g)), g, n});
                }
            });
        }

        // The split sequence of each segment: its own two endpoint nodes first
        // and last, everything else in arc order between them. Forcing the
        // endpoints keeps the reassembled chain continuous; in every
        // non-degenerate case they are already the extremes of the order.
        std::vector<std::vector<std::uint32_t>> split(raw.size());
        for (std::size_t i = 0; i < raw.size(); ++i) {
            std::ranges::sort(met[i], detail::ArcKey::before);
            const std::uint32_t a = nodes.id_of(grid.snap(raw[i].a));
            const std::uint32_t b = nodes.id_of(grid.snap(raw[i].b));
            split[i].push_back(a);
            if (a == b) {
                continue;  // the segment collapsed under snapping
            }
            for (const detail::ArcKey& key : met[i]) {
                if (key.id != a && key.id != b) {
                    split[i].push_back(key.id);
                }
            }
            split[i].push_back(b);
        }

        // --- Step 5: reassemble, then merge over the edge set ---------------
        std::vector<std::vector<std::uint32_t>> runs(pts.size());
        std::vector<std::vector<EdgeProperties>> new_props(pts.size());
        for (std::size_t s = 0; s < segs.size(); ++s) {
            const std::size_t c = segs[s].chain;
            for (const std::uint32_t id : split[s]) {
                if (runs[c].empty()) {
                    runs[c].push_back(id);
                } else if (runs[c].back() != id) {
                    // Zero-length edges -- two consecutive positions with the
                    // same node id -- collapse here and nowhere else.
                    runs[c].push_back(id);
                    new_props[c].push_back(props[c][segs[s].edge]);
                }
            }
        }
        for (std::size_t c = 0; c < runs.size(); ++c) {
            if (is_closed(pslg.chains()[c].role) && runs[c].size() > 1 &&
                runs[c].front() == runs[c].back()) {
                runs[c].pop_back();
            }
        }

        std::vector<Point2> vertices;
        vertices.reserve(nodes.size());
        for (const GridPoint& g : nodes.points()) {
            vertices.push_back(grid.world(g));
        }

        // --- The step 05-noder.md's list does not have: ring re-validation --
        if (NodeOutcome bad =
                detail::revalidate_rings<K>(pslg, options.spacing, vertices, runs);
            bad.status != NodeStatus::Ok) {
            return bad;
        }

        // The merge, at the node-id edge-key dedup and nowhere else. On noded
        // output the integer key IS the overlap relation -- amended 14(a) says
        // two edges that overlap at all have equal node-id pairs -- so a
        // collinear-overlap classification here would re-derive, in floating
        // point, what the key already decides. The reduce is over an UNORDERED
        // collection and union's commutativity, associativity and idempotence are
        // the whole licence for that.
        std::map<std::pair<std::uint32_t, std::uint32_t>, EdgeProperties> merged;
        for (std::size_t c = 0; c < runs.size(); ++c) {
            for (std::size_t k = 0; k < new_props[c].size(); ++k) {
                const auto key = detail::edge_key(runs[c], k);
                merged[key] = merged[key] | new_props[c][k];
            }
        }

        // --- Step 6: verify -------------------------------------------------
        std::vector<Chain> chains;
        std::vector<std::uint32_t> flat;
        std::vector<EdgeProperties> edge_properties;
        chains.reserve(runs.size());
        for (std::size_t c = 0; c < runs.size(); ++c) {
            chains.push_back(Chain{static_cast<std::uint32_t>(flat.size()),
                                   static_cast<std::uint32_t>(runs[c].size()),
                                   pslg.chains()[c].role, pslg.chains()[c].properties});
            flat.insert(flat.end(), runs[c].begin(), runs[c].end());
            for (std::size_t k = 0; k < new_props[c].size(); ++k) {
                edge_properties.push_back(merged[detail::edge_key(runs[c], k)]);
            }
        }

        std::vector<std::uint32_t> node_of_input_vertex;
        node_of_input_vertex.reserve(input_keys.size());
        for (const GridPoint& g : input_keys) {
            node_of_input_vertex.push_back(nodes.id_of(g));
        }

        last = NodedPslgBuilder{grid,
                                std::move(vertices),
                                std::move(chains),
                                std::move(flat),
                                std::move(edge_properties),
                                std::move(node_of_input_vertex)}
                   .template build<K>();
        if (last.status == NodeStatus::Ok || last.status == NodeStatus::MalformedOutput) {
            return last;
        }

        // A guarantee-14 refusal below the cap: iterate from step 1 on the SPLIT
        // segments, which are now on the grid, so the next round sees exactly the
        // geometry the verification just rejected.
        for (std::size_t c = 0; c < runs.size(); ++c) {
            pts[c].clear();
            for (const std::uint32_t id : runs[c]) {
                pts[c].push_back(grid.world(nodes[id]));
            }
            props[c] = new_props[c];
        }
    }

    return detail::refuse(
        NodeStatus::NotConverged,
        std::format("the constraint set did not settle in {} rounds at spacing {}: {}",
                    options.max_rounds, options.spacing,
                    last.message.empty() ? std::string{describe(NodeStatus::NotConverged)}
                                         : last.message));
}

}  // namespace terrain::noding
