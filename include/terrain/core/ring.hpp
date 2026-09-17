#pragma once

// The Ring concept, its two non-owning models, and the algorithms over them.
//
// A ring is a cyclic sequence of *distinct* vertices with closure implied,
// never stored. Both models reject a stored closure at construction: if both
// encodings were accepted every ring would have two spellings and every
// off-by-one ring bug would live in the gap between them.
//
// Both models are *views*. Increment 3's PSLG owns one vertex buffer and one
// flat index buffer, and hands the algorithms sub-spans of them; an owning ring
// would force a gather-and-copy per ring, which is both an allocation per ring
// and the exact setup for detria's use-after-free -- addOutline stores a span
// and does not copy, so a ring buffer built inside a loop dangles at
// triangulate(). The rvalue-vector constructors are therefore = delete'd, which
// turns the most likely spelling of that mistake into a compile error. The
// residual hazard -- a view over a vector that is *reallocated* while the view
// lives -- is unguardable in C++ and is the caller's to avoid.
//
// PointRing's span constructor is `explicit`, and neither model has a
// conversion operator to std::span<const Point2>. That is load-bearing, not
// stylistic: bbox.hpp's free bounding_box(std::span<const Point2>) and the
// bounding_box(const R&) template below coexist unambiguously only because no
// ring model converts to a span. Do not relax either.
//
// K is the first template parameter on every kernel-using algorithm, so call
// sites read point_in_ring<DefaultKernel>(ring, p) with R deduced. There is no
// default for K: defaulting it would make this header include
// default_kernel.hpp and drag the compiled terrain_predicates target into every
// consumer of a pure-header core type.

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

// Exactly two obligations, so a future third model owes nothing else. size()
// counts distinct vertices and never the implied closure; vertex(i) has the
// precondition i < size() and returns a reference, so the algorithms iterate
// without copying.
template <typename R>
concept Ring = requires(const R& r, std::size_t i) {
    { r.size() } -> std::same_as<std::size_t>;
    { r.vertex(i) } -> std::same_as<const Point2&>;
};

namespace detail {

// Shared by both models: the degeneracy policy of the increment, in one place.
// Finiteness is deliberately *not* checked -- it is a precondition of the
// vertex buffer, established once at the PSLG and Python boundaries, and
// scanning n coordinates on every construction is the wrong cost for a view
// that may be built per query. all_finite() is how a caller asks.
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

// The contiguous model: tests, Python ingress, a draped DEM border.
class PointRing {
public:
    explicit PointRing(std::span<const Point2> vertices) : vertices_{vertices} {
        detail::check_ring_size(vertices.size());
        detail::check_ring_closure(vertices.front(), vertices.back());
    }

    // A view over a temporary dangles on the next line. Compile error, not a
    // sanitizer finding.
    PointRing(std::vector<Point2>&&) = delete;

    [[nodiscard]] std::size_t size() const noexcept { return vertices_.size(); }

    // Precondition: i < size().
    [[nodiscard]] const Point2& vertex(std::size_t i) const noexcept { return vertices_[i]; }

private:
    std::span<const Point2> vertices_;
};

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

// Directed, and wrapping at the last vertex: edge(r, size() - 1) is the closing
// edge. A free function rather than a concept requirement, so the concept stays
// two lines.
template <Ring R>
[[nodiscard]] Segment2 edge(const R& r, std::size_t i) {
    return Segment2{r.vertex(i), r.vertex((i + 1) % r.size())};
}

// For increment 3's validator, and for tests. O(n), and the only finiteness
// question this header answers.
template <Ring R>
[[nodiscard]] bool all_finite(const R& r) noexcept {
    for (std::size_t i = 0; i < r.size(); ++i) {
        const Point2& p = r.vertex(i);
        if (!std::isfinite(p.x) || !std::isfinite(p.y)) {
            return false;
        }
    }
    return true;
}

// Never the empty box: a ring has at least three vertices. May be degenerate in
// an axis, and is silently corrupt on a non-finite ring -- see all_finite.
template <Ring R>
[[nodiscard]] Box2 bounding_box(const R& r) {
    Box2 box;
    for (std::size_t i = 0; i < r.size(); ++i) {
        box.expand(r.vertex(i));
    }
    return box;
}

// APPROXIMATE, kernel-free, and its SIGN MUST NEVER DRIVE A TOPOLOGY DECISION.
// Over a sliver at UTM33 magnitudes the sum cancels to noise, and taking its
// sign is how a hole gets classified as an outline. orientation<K> is the
// instrument for that question and it does not sum. This is the magnitude, and
// only the magnitude.
//
// Shoelace translated to vertex(0), which is what keeps the error at a relative
// scale instead of the absolute one an untranslated sum produces at UTM33
// northings. The two terms a fan over vertex(0) would contribute zero to are
// simply not visited.
template <Ring R>
[[nodiscard]] double signed_area(const R& r) noexcept {
    const Point2& origin = r.vertex(0);
    double twice = 0.0;
    for (std::size_t i = 1; i + 1 < r.size(); ++i) {
        twice += cross(r.vertex(i) - origin, r.vertex(i + 1) - origin);
    }
    return 0.5 * twice;
}

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
template <Ring R>
[[nodiscard]] std::size_t extreme_vertex(const R& r) {
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
template <pred::GeometryKernel K, Ring R>
[[nodiscard]] pred::Orientation orientation(const R& r) {
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

// The underlying values are the sign convention, as in pred::Orientation.
enum class PointInRing : int {
    Outside = -1,
    Boundary = 0,
    Inside = 1,
};

// Three-valued and exact, even-odd, and independent of winding.
//
// Three-valued rather than bool because CDT hole removal, DEM-border draping
// and the noder all need Boundary distinguished; folding it into either bucket
// is how boundary vertices get deleted along with their hole. Even-odd rather
// than nonzero-winding because even-odd needs no consistent orientation, and a
// self-intersecting ring -- accepted and not detected here -- has none.
//
// This function performs NO ARITHMETIC OF ITS OWN. Every numeric decision is a
// K::orient2d call or a comparison between two coordinates the caller supplied:
// no division, no constructed x-intersection, no accumulation. The observable
// consequence is that translating the whole problem by an exactly representable
// UTM33 offset cannot change an answer.
//
// Boundary is tested per edge and returns early, before the parity update, so a
// point on an edge or on a vertex is Boundary whatever the parity would have
// said. The parity rule is half-open in y -- edge (u, v) counts when
// (u.y <= p.y) != (v.y <= p.y) -- which counts a local extremum zero times and
// a horizontal edge on the ray not at all. Zero-length edges likewise
// contribute nothing, both endpoints falling the same side of the half-open
// test, so they need no special case.
template <pred::GeometryKernel K, Ring R>
[[nodiscard]] PointInRing point_in_ring(const R& r, const Point2& p) {
    const std::size_t n = r.size();
    bool inside = false;

    for (std::size_t i = 0; i < n; ++i) {
        const Point2& u = r.vertex(i);
        const Point2& v = r.vertex((i + 1) % n);

        if (on_segment<K>(Segment2{u, v}, p)) {
            return PointInRing::Boundary;
        }
        const bool below_u = u.y <= p.y;
        if (below_u == (v.y <= p.y)) {
            continue;
        }
        // The edge straddles the ray. It crosses to the +x side of p when p
        // lies left of an upward edge, or right of a downward one -- decided by
        // the kernel, never by a computed intersection.
        const pred::Orientation side = below_u ? pred::Orientation::CounterClockwise
                                               : pred::Orientation::Clockwise;
        if (K::orient2d(u, v, p) == side) {
            inside = !inside;
        }
    }
    return inside ? PointInRing::Inside : PointInRing::Outside;
}

}  // namespace terrain
