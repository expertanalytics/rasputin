#pragma once

// The constraint edge set and the per-triangle constrained-edge mask.
//
// The mask is OURS to compute, from the Pslg and the output triangles, with no
// backend involvement: detria exposes no public per-edge query, and getTopology
// is documented "mostly used for testing" and returns an internal half-edge
// type, so using it would drag the library's topology vocabulary across the
// seam for one bit per edge. Computing it here is the good outcome -- it makes
// the mask testable without a triangulation at all.
//
// A SORTED std::vector<std::uint64_t> WITH BINARY SEARCH, NOT A HASH SET. The
// key is exact and small, so hashing buys nothing; one allocation of known size
// beats a rehashing container; the structure is read-only after construction
// and therefore shareable by every refinement thread without synchronisation,
// per parallel_refinement.md, which also asks for contiguity by name. Sorting
// dedups for free, which matters because a breakline may legitimately repeat an
// edge. This header includes no <unordered_map> and no <unordered_set>.
//
// This is not "dedup" in the sense increment 3 prohibited. That prohibition is
// about COORDINATE keys, which need a snap grid that does not exist until
// increment 5. These keys are vertex INDICES: exact integers, no tolerance
// question.
//
// THIS WORKS ONLY BECAUSE THE BACKEND NEVER SPLITS A CONSTRAINT EDGE. detria
// inserts no Steiner points -- a vertex lying exactly on a constraint is a hard
// PointOnConstrainedEdge failure, not a split -- so every constraint edge
// appears in the output as exactly one triangulation edge and a pair lookup is
// exact. That assumption belongs to DetriaBackend, not to the seam: a backend
// that splits edges (poly2tri does insert points) would need the mask computed
// by walking the split chain. It is stated here because this is the code that
// silently becomes wrong when the backend changes.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/pslg.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

namespace terrain::cdt {

// Every constraint edge of a Pslg, as an unordered index pair, sorted once.
class ConstraintEdgeSet {
public:
    explicit ConstraintEdgeSet(const Pslg& pslg) {
        std::size_t total = 0;
        for (std::size_t c = 0; c < pslg.chains().size(); ++c) {
            total += pslg.edge_count(c);
        }
        // The pre-pass above is what makes the "one allocation of known size"
        // rationale at the top of this header true, rather than redundant work:
        // without it the sorted vector would grow and reallocate while filling.
        keys_.reserve(total);

        for (std::size_t c = 0; c < pslg.chains().size(); ++c) {
            const std::span<const std::uint32_t> idx = pslg.indices_of(c);
            const std::size_t n = idx.size();
            // edge_count(c) is the authority on how many edges a chain has: a
            // ring's closing edge is included and an open breakline's is not,
            // which is the whole difference between the two roles. The modulus
            // then spells the closing pair without a role branch.
            for (std::size_t k = 0; k < pslg.edge_count(c); ++k) {
                keys_.push_back(key(idx[k], idx[(k + 1) % n]));
            }
        }

        std::ranges::sort(keys_);
        const auto dup = std::ranges::unique(keys_);
        keys_.erase(dup.begin(), dup.end());
    }

    [[nodiscard]] bool contains(std::uint32_t a, std::uint32_t b) const noexcept {
        return std::ranges::binary_search(keys_, key(a, b));
    }

    [[nodiscard]] std::size_t size() const noexcept { return keys_.size(); }

    // (min(a, b) << 32) | max(a, b). Exposed so a test can pin the symmetry
    // rather than infer it: key(a, b) == key(b, a), for every pair.
    [[nodiscard]] static constexpr std::uint64_t key(std::uint32_t a, std::uint32_t b) noexcept {
        const std::uint32_t lo = a < b ? a : b;
        const std::uint32_t hi = a < b ? b : a;
        return (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint64_t>(hi);
    }

private:
    std::vector<std::uint64_t> keys_;
};

// Bit e is set iff the edge (v[e], v[(e + 1) % 3]) is a constraint edge.
//
// NOT the "edge e is opposite vertex e" convention CGAL uses, which a reader
// with that background will assume. That one maps edge e to (v[e+1], v[e+2]) --
// our bit e+1 -- so the two are A ROTATION OF EACH OTHER and agree on every
// triangle with zero or three constrained edges. Only a triangle with exactly
// one constrained edge tells them apart, which is why the suite pins the
// convention on exactly that shape.
[[nodiscard]] inline std::uint8_t constrained_mask(const ConstraintEdgeSet& edges,
                                                   const TriangleIndices& t) noexcept {
    std::uint8_t mask = 0;
    for (std::size_t e = 0; e < 3; ++e) {
        if (edges.contains(t[e], t[(e + 1) % 3])) {
            mask |= static_cast<std::uint8_t>(1u << e);
        }
    }
    return mask;
}

}  // namespace terrain::cdt
