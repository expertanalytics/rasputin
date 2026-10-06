#pragma once

// IndexedRing, the non-owning ring view, and the algorithms over it.
//
// A ring is a cyclic sequence of *distinct* vertices with closure implied,
// never stored. The constructor rejects a stored closure: if both encodings
// were accepted every ring would have two spellings and every off-by-one ring
// bug would live in the gap between them.
//
// The ring is a *view*. Increment 3's PSLG owns one vertex buffer and one
// flat index buffer, and hands the algorithms sub-spans of them; an owning ring
// would force a gather-and-copy per ring, which is both an allocation per ring
// and the exact setup for detria's use-after-free -- addOutline stores a span
// and does not copy, so a ring buffer built inside a loop dangles at
// triangulate(). The rvalue-vector constructors are therefore = delete'd, which
// turns the most likely spelling of that mistake into a compile error. The
// residual hazard -- a view over a vector that is *reallocated* while the view
// lives -- is unguardable in C++ and is the caller's to avoid.
//
// IndexedRing has no conversion operator to std::span<const Point2>: a ring is
// its chain of indices, and a span conversion would hand out the whole vertex
// buffer as if it were the ring.
//
// K is the first template parameter on every kernel-using algorithm, so call
// sites read orientation<DefaultKernel>(ring). There is no default for K:
// defaulting it would make this header include default_kernel.hpp and drag the
// compiled terrain_predicates target into every consumer of a pure-header core
// type.

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace terrain {

namespace detail {

// The degeneracy policy of the increment, in one place.
// Finiteness is deliberately *not* checked -- it is a precondition of the
// vertex buffer, established once at the PSLG and Python boundaries, and
// scanning n coordinates on every construction is the wrong cost for a view
// that may be built per query.
inline void check_ring_size(std::size_t n) {
    if (n < 3) {
        throw std::invalid_argument("terrain: a ring needs at least three vertices");
    }
}

// Called only after check_ring_size, which is what makes the first and last
// vertices safe to name.
inline void check_ring_closure(const Point2& first, const Point2& last) {
    if (first == last) {
        throw std::invalid_argument(
            "terrain: a ring's closure is implied and must not be stored; "
            "drop the repeated last vertex");
    }
}

}  // namespace detail

// The zero-copy model: vertex(i) is vertices[chain[i]], with the chain slice
// already sub-spanned by the caller. This is what increment 3 hands down.
class IndexedRing {
public:
    IndexedRing(std::span<const Point2> vertices, std::span<const std::uint32_t> chain)
        : vertices_{vertices}, chain_{chain} {
        detail::check_ring_size(chain.size());
        // Exactly the two indices the closure check below dereferences. The
        // full per-index range check is O(n) and belongs to the PSLG's one-time
        // validation, not to a view that may be built per query -- but the one
        // constructor documented as not validating indices must still not
        // perform two unchecked out-of-range reads of its own.
        if (chain.front() >= vertices.size() || chain.back() >= vertices.size()) {
            throw std::invalid_argument(
                "terrain: a ring's first and last chain indices must be in range");
        }
        // Closure is a question about the points, not about the indices: two
        // distinct indices onto coincident vertices close the ring just as
        // surely as one index used twice.
        detail::check_ring_closure(vertices[chain.front()], vertices[chain.back()]);
    }

    IndexedRing(std::vector<Point2>&&, std::span<const std::uint32_t>) = delete;
    IndexedRing(std::span<const Point2>, std::vector<std::uint32_t>&&) = delete;

    [[nodiscard]] std::size_t size() const noexcept { return chain_.size(); }

    // Precondition: i < size(), and chain_[i] < vertices_.size() -- the latter
    // validated once by the PSLG, not here.
    [[nodiscard]] const Point2& vertex(std::size_t i) const noexcept {
        return vertices_[chain_[i]];
    }

private:
    std::span<const Point2> vertices_;
    std::span<const std::uint32_t> chain_;
};

namespace detail {

// Local to this header on purpose: point.hpp gets no ordering and no
// operator<=>. The comparison is only meaningful under the finiteness
// precondition a ring carries and Point2 does not, so putting it on the type
// would make it available exactly where it is unsafe.
[[nodiscard]] constexpr bool lexicographically_before(const Point2& a, const Point2& b) noexcept {
    return a.y < b.y || (a.y == b.y && a.x < b.x);
}

// Any convex-hull vertex is convex, so the tie-break here is non-normative:
// min-y/max-y and min-x/max-x are all equally correct and the choice is
// unobservable. Lowest index wins ties, because the comparison is strict.
[[nodiscard]] inline std::size_t extreme_vertex(const IndexedRing& r) {
    std::size_t best = 0;
    for (std::size_t i = 1; i < r.size(); ++i) {
        if (lexicographically_before(r.vertex(i), r.vertex(best))) {
            best = i;
        }
    }
    return best;
}

}  // namespace detail

// EXACT, and it sums nothing: one kernel call on the extreme vertex and its
// neighbours, so there is no cancellation to lose the answer to. Exact only for
// a simple ring, which is accepted -- a non-simple ring has no well-defined
// winding to report.
//
// When the first triple is collinear -- a repeated extreme vertex, or a
// collinear run through it -- prev and next advance INDEPENDENTLY, never in
// lockstep. A lockstep walk returns Collinear for the proper triangle
// {A, A, B, C}, where every symmetric pair about the extreme vertex is
// collinear and the walk exhausts. The nested walk below returns Collinear only
// when every vertex is collinear with the extreme one, which is exactly the
// all-collinear ring the policy specifies it for.
template <pred::GeometryKernel K>
[[nodiscard]] pred::Orientation orientation(const IndexedRing& r) {
    const std::size_t n = r.size();
    const std::size_t e = detail::extreme_vertex(r);
    const Point2& v = r.vertex(e);

    for (std::size_t back = 1; back < n; ++back) {
        const Point2& prev = r.vertex((e + n - back) % n);
        for (std::size_t fwd = 1; fwd < n; ++fwd) {
            const pred::Orientation o = K::orient2d(prev, v, r.vertex((e + fwd) % n));
            if (o != pred::Orientation::Collinear) {
                return o;
            }
        }
    }
    return pred::Orientation::Collinear;
}

}  // namespace terrain
