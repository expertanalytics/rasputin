#pragma once

// The directed 2D segment and its one predicate.
//
// Segment2's equality is *ordered*: {a, b} != {b, a}. Ring edges are directed
// -- the parity rule in point_in_ring and every future constraint edge depend
// on which way an edge runs -- so a segment comparing equal to its own reverse
// would make "the same edge" an ambiguous phrase in the one module where it
// must not be. An unordered comparison, if ever needed, gets its own name.
//
// on_segment is the only predicate here, and it is exact: collinearity comes
// from the kernel, betweenness from closed comparisons between coordinates the
// caller supplied. There is no division, no parameter t, no constructed
// intersection point and no tolerance. Everything constructive -- intersection
// points, distances, projections -- is deliberately absent until the noder
// exists, because construction rounds and the rounding target is the snap grid.
//
// This header includes the predicate *concept* headers but never
// default_kernel.hpp: naming a default kernel here would drag the compiled
// terrain_predicates target into every consumer of a pure-header core type.
// Callers name their kernel.

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <algorithm>

namespace terrain {

struct Segment2 {
    Point2 a{};
    Point2 b{};

    friend constexpr bool operator==(const Segment2&, const Segment2&) = default;
};

[[nodiscard]] constexpr Segment2 reversed(const Segment2& s) noexcept {
    return Segment2{s.b, s.a};
}

// Exact, not approximate: a segment one ulp long is not degenerate.
[[nodiscard]] constexpr bool is_degenerate(const Segment2& s) noexcept {
    return s.a == s.b;
}

// Validating, like the free bounding_box it delegates to: throws
// std::invalid_argument on a non-finite endpoint.
[[nodiscard]] inline Box2 bbox(const Segment2& s) {
    const Point2 points[] = {s.a, s.b};
    return bounding_box(std::span<const Point2>{points});
}

// Point on closed segment, exact. K is the first template parameter so call
// sites read on_segment<DefaultKernel>(s, p); there is no default for it, for
// the include-dependency reason in the header comment.
//
// Collinearity is the kernel's answer; betweenness is a closed comparison
// against the segment's own coordinate range. No arithmetic of this function's
// own appears anywhere, so under an exact kernel the result is exact for every
// finite input.
//
// A degenerate segment needs no branch: a zero-length segment makes every point
// collinear, so the betweenness comparisons collapse to p == s.a on their own.
// An explicit special case here is how the degenerate answer gets subtly wrong.
template <pred::GeometryKernel K>
[[nodiscard]] bool on_segment(const Segment2& s, const Point2& p) {
    if (K::orient2d(s.a, s.b, p) != pred::Orientation::Collinear) {
        return false;
    }
    return std::min(s.a.x, s.b.x) <= p.x && p.x <= std::max(s.a.x, s.b.x) &&
           std::min(s.a.y, s.b.y) <= p.y && p.y <= std::max(s.a.y, s.b.y);
}

}  // namespace terrain
