#pragma once

// The noded planar straight-line graph: the constraint set after snap rounding,
// with every crossing and every hot-pixel incidence resolved into a node.
//
// A NodedPslg exists only if it passed NodedPslgBuilder's verification, so every
// module above this one -- the CDT wrapper at 5c, refinement later -- may rely on
// the guarantees below WITHOUT re-checking. The private constructor with a single
// friend is what makes that true rather than hoped for, exactly as core/pslg.hpp
// does it.
//
// On top of Pslg's guarantees 1-8 and 10 (9 is replaced by
// node_of_input_vertex), and numbered as docs/increments/05b-noder-driver.md
// numbers them:
//  11. Every coordinate is grid().world(g) for some GridPoint within
//      kMaxGridIndex. Checkable bitwise: grid().snapped(v) == v.
//  12. No two vertices are equal.
//  13. No edge is zero-length; no chain repeats a consecutive index.
//  14. (a) No two edges cross in their interiors -- and two edges that overlap
//      collinearly have EQUAL node-id pairs, which is what makes a road noded
//      along a river legal. (b) No node's cell meets an edge it is not an
//      endpoint of: the hot-pixel clause, strictly stronger than exact
//      incidence.
//  15. edge_properties() is DENSE and index-aligned with the flat edge
//      enumeration (c, k) at offset edge_base(c) + k, each entry the UNION over
//      the property sets of every input chain that contributed geometry to that
//      edge. 15 is established by the driver and checked by its property suite,
//      not by the builder -- see noding/noded_pslg_builder.hpp.
//  16. Chain order, roles and property sets are preserved index-for-index with
//      the input Pslg. Only begin and count change.
//
// THIS HEADER IS KERNEL-FREE. Like Pslg it is a type; the validation that
// establishes its promise is a separate, kernel-templated header. It includes
// nothing from predicates/ and nothing from noding/ beyond the forward
// declaration of its one friend.
//
// NodedPslg IS NOT A SUBCLASS OF Pslg AND THERE IS NO CONVERSION EITHER WAY.
// That is the whole architectural point: once cdt::triangulate takes a
// NodedPslg, un-noded input is unrepresentable at the entry point rather than
// diagnosed inside it.
//
// It does, however, have the same SHAPE as Pslg -- vertices, chains,
// chain_indices, indices_of -- because the Python renderer's scene join is
// structural (viz/protocols.py's PslgLike) and must be fed the noded graph, not
// the input one. Two independent reasons for one decision.
//
// Like Pslg, NodedPslg STORES NO SPAN AND NO ITERATOR INTO ITSELF: every view is
// computed on demand, which is what makes the defaulted copy and move correct.

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace terrain {

namespace noding {
class NodedPslgBuilder;
}  // namespace noding

class NodedPslg {
public:
    NodedPslg(const NodedPslg&) = default;
    NodedPslg& operator=(const NodedPslg&) = default;
    NodedPslg(NodedPslg&&) = default;
    NodedPslg& operator=(NodedPslg&&) = default;
    ~NodedPslg() = default;

    [[nodiscard]] std::span<const Point2> vertices() const noexcept { return vertices_; }
    [[nodiscard]] std::span<const Chain> chains() const noexcept { return chains_; }
    [[nodiscard]] std::span<const std::uint32_t> chain_indices() const noexcept {
        return chain_indices_;
    }

    // Sub-span of chain_indices(). Precondition: c < chains().size().
    [[nodiscard]] std::span<const std::uint32_t> indices_of(std::size_t c) const noexcept {
        assert(c < chains_.size());
        const Chain& ch = chains_[c];
        return std::span<const std::uint32_t>{chain_indices_}.subspan(ch.begin, ch.count);
    }

    // Precondition: is_closed(chains()[c].role). NOT noexcept, for the reason
    // core/pslg.hpp gives for its own ring(): IndexedRing rechecks invariants the
    // compiler cannot see are already established, and a noexcept resting on that
    // is a std::terminate waiting for the argument to stop being true.
    [[nodiscard]] IndexedRing ring(std::size_t c) const {
        assert(is_closed(chains_[c].role));
        return IndexedRing{vertices(), indices_of(c)};
    }

    // The only place the roles differ. Precondition: c < chains().size().
    [[nodiscard]] std::size_t edge_count(std::size_t c) const noexcept {
        assert(c < chains_.size());
        const Chain& ch = chains_[c];
        const auto n = static_cast<std::size_t>(ch.count);
        return is_closed(ch.role) ? n : n - 1;
    }

    // Directed, and (k + 1) % n unconditionally, so edge(c, edge_count(c) - 1) is
    // the closing edge of a ring and the last segment of a breakline. Same ruling
    // and same reasoning as core/pslg.hpp's edge().
    [[nodiscard]] Segment2 edge(std::size_t c, std::size_t k) const noexcept {
        assert(k < edge_count(c));
        const std::span<const std::uint32_t> idx = indices_of(c);
        const std::size_t n = idx.size();
        return Segment2{vertices_[idx[k]], vertices_[idx[(k + 1) % n]]};
    }

    [[nodiscard]] SnapGrid grid() const noexcept { return grid_; }

    // Guarantee 15's array. DENSE -- one entry per output edge -- and not a
    // sparse override set: after noding there is no single source value for an
    // exception to override, because an output edge can descend from two input
    // chains at once. The empty set is a legal value meaning *unclassified*.
    [[nodiscard]] std::span<const EdgeProperties> edge_properties() const noexcept {
        return edge_properties_;
    }

    // Guarantee 9's replacement: input vertex position -> node id. Total over the
    // input's vertex array, unreferenced vertices included.
    [[nodiscard]] std::span<const std::uint32_t> node_of_input_vertex() const noexcept {
        return node_of_input_vertex_;
    }

    // The offset of chain c's first edge in the flat enumeration. One entry past
    // the end, the same convention chain_indices already uses:
    // edge_base(chains().size()) == edge_properties().size(). It exists so that
    // no caller holding (c, k) recomputes a prefix sum -- the kind of derived
    // arithmetic that is right in four places and wrong in the fifth.
    // Precondition: c <= chains().size().
    [[nodiscard]] std::size_t edge_base(std::size_t c) const noexcept {
        assert(c < edge_base_.size());
        return edge_base_[c];
    }

private:
    // Private, with exactly one friend: the only way to hold a NodedPslg is to
    // have passed verification. Friending the class rather than node<K> keeps
    // this header free of pred::GeometryKernel, which is the point of it being
    // kernel-free at all.
    NodedPslg() = default;
    friend class noding::NodedPslgBuilder;

    std::vector<Point2> vertices_;
    std::vector<std::uint32_t> chain_indices_;
    std::vector<Chain> chains_;
    std::vector<EdgeProperties> edge_properties_;
    std::vector<std::uint32_t> node_of_input_vertex_;
    std::vector<std::size_t> edge_base_;
    SnapGrid grid_{1.0};
};

}  // namespace terrain
