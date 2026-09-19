#pragma once

// The planar straight-line graph: the constraint set, checked once and never
// again.
//
// A Pslg exists only if it passed PslgBuilder's validation, so every module
// above this one -- the noder, the CDT wrapper, refinement -- may rely on the
// guarantees listed below WITHOUT re-checking. That list is the entire point of
// the type; the private constructor is what makes it true rather than hoped
// for.
//
// What a valid Pslg guarantees:
//   1. Every coordinate in vertices() is finite, including unreferenced ones.
//   2. Every value in chain_indices() is < vertices().size().
//   3. chains() is non-empty and at least one chain has role Outer.
//   4. Every closed chain has count >= 3; every breakline has count >= 2.
//   5. No closed chain stores its closure: its first and last vertices are
//      distinct POINTS.
//   6. Every Outer ring is counterclockwise and every Hole ring is clockwise
//      under the kernel that validated it. Neither is collinear.
//   7. indices_of(0..n) partition chain_indices() in chain order, contiguously,
//      with no gap and no overlap.
//   8. ring(c) on a closed chain never throws.
//   9. The vertex buffer is element-wise equal to what the builder accumulated:
//      no dedup, no reordering, no reversal, no insertion.
//  10. A const Pslg is safe for concurrent read from any number of threads.
//
// And what it explicitly does NOT promise, and downstream must not assume: that
// any chain is simple, that two chains are disjoint, that a hole lies inside an
// outer ring, that vertices are distinct, or that every vertex is referenced.
// THEREFORE A VALID PSLG IS NOT A VALID CDT INPUT. It is a valid *noder* input:
// it asserts everything that can be decided without constructing a point, and
// nothing that cannot. Increment 5's NodedPslg carries the rest.
//
// Pslg owns all three buffers and is immutable after construction. There is no
// mutator, no non-const accessor, no lazy cache and no mutable member, which is
// what makes guarantee 10 hold with no synchronisation.
//
// Pslg stores NO SPAN AND NO ITERATOR INTO ITSELF. Every view is computed on
// demand from the vectors, and that is what makes the defaulted copy and move
// correct: a cached IndexedRing member would make a copied Pslg point at the
// original's buffers, and that bug survives every test that never copies. Do
// not add one.

#include <terrain/core/edge_properties.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace terrain {

enum class ChainRole : std::uint8_t { Outer, Hole, Breakline };

[[nodiscard]] constexpr bool is_closed(ChainRole r) noexcept {
    return r == ChainRole::Outer || r == ChainRole::Hole;
}

// `count` counts DISTINCT vertices. For Outer and Hole closure is implied and
// the repeated last index is not stored, exactly as Ring specifies. For
// Breakline the chain is open and `count` is the vertex count of an open
// polyline.
//
// role's default member initializer is Breakline and that is deliberate.
// Breakline is the role that promises the least -- no winding contract, no
// implied closure, no interior -- so a Chain that reaches a reader without
// having been through the validator cannot claim it is a domain boundary.
// Aggregate initialization of a missing enum member would otherwise give
// ChainRole(0). If the enumerator order is ever changed, this NSDMI is what
// keeps the property; it is not decoration.
//
// properties is the chain's feature set, permitted on every role -- a wide
// river or a lake is legitimately an area feature, and so is a walled
// enclosure -- and never validated, because it is data, not structure. It is
// NAMED `properties` rather than `edge_properties`: a Chain has no edges of its
// own that outlive it, and the relation is inheritance downward, every output
// edge a chain contributes geometry to receiving that chain's set.
struct Chain {
    std::uint32_t begin{};
    std::uint32_t count{};
    ChainRole role{ChainRole::Breakline};
    EdgeProperties properties{};

    friend constexpr bool operator==(const Chain&, const Chain&) = default;
};

class PslgBuilder;

class Pslg {
public:
    Pslg(const Pslg&) = default;
    Pslg& operator=(const Pslg&) = default;
    Pslg(Pslg&&) = default;
    Pslg& operator=(Pslg&&) = default;

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

    // Zero-copy view of a closed chain: two spans, no allocation, no gather.
    // Precondition: is_closed(chains()[c].role).
    //
    // NOT noexcept, even though it cannot throw. IndexedRing's constructor
    // rechecks size, closure and the two boundary indices, and every one of
    // those was established by the validator -- but the compiler cannot know
    // that, and a noexcept resting on an argument the compiler cannot see is a
    // std::terminate waiting for the argument to stop being true. Those checks
    // being provably unreachable here is the cleanest evidence that the
    // validator's stage order is right.
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

    // Directed, and (k + 1) % n UNCONDITIONALLY, with no role branch:
    // edge(c, edge_count(c) - 1) is therefore the closing edge of a ring and
    // the last segment of a breakline. That is the whole reason this accessor
    // exists -- neither the noder's broad phase nor the CDT wrapper's
    // setConstrainedEdge loop should write (k + 1) % count itself and get the
    // open case wrong.
    //
    // Precondition: k < edge_count(c). THE MODULUS IS NOT THERE TO SERVICE AN
    // OUT-OF-RANGE k. It is there so that a caller who violates the
    // precondition gets a wrong answer instead of undefined behaviour: a bare
    // k + 1 on the open arm reads the flat index buffer out of bounds at
    // edge(c, count - 1), and no in-contract call can tell the two spellings
    // apart. A distinction no test can enforce is not a design freedom.
    [[nodiscard]] Segment2 edge(std::size_t c, std::size_t k) const noexcept {
        assert(k < edge_count(c));
        const std::span<const std::uint32_t> idx = indices_of(c);
        const std::size_t n = idx.size();
        return Segment2{vertices_[idx[k]], vertices_[idx[(k + 1) % n]]};
    }

private:
    // Private, with PslgBuilder a friend: the only way to hold a Pslg is to
    // have passed validation. A hand-built std::vector<Chain> whose begin/count
    // do not match the flat buffer has no public path to this type.
    Pslg(std::vector<Point2> vertices, std::vector<std::uint32_t> chain_indices,
         std::vector<Chain> chains)
        : vertices_{std::move(vertices)},
          chain_indices_{std::move(chain_indices)},
          chains_{std::move(chains)} {}

    friend class PslgBuilder;

    std::vector<Point2> vertices_;
    std::vector<std::uint32_t> chain_indices_;
    std::vector<Chain> chains_;
};

// What deliberately does not exist here: there is NO closed_indices_of(c).
//
// A closed span -- a ring's indices with its first index repeated at the end --
// has no producer in this project and no consumer either. The CDT wrapper was
// expected to need one, because detria's addOutline and addHole take a span
// that bounds a region; measured against the vendored header, an OPEN span is
// closed by detria itself, so indices_of(c) goes straight in. Increment 4
// deleted the closed_index_buffer_size that existed to size the scratch buffer
// for the caller that is now never going to be written.
//
// So the rule survives in its stronger form: no closing index is materialised
// anywhere in the project, and there is therefore no accessor whose result
// would have to be allocated. If a re-pinned backend ever stops auto-closing --
// the backend suite's characterisation test is what would say so -- the fix is
// one scratch buffer inside the wrapper, sub-spans handed out from it, and
// still nothing here.

}  // namespace terrain
