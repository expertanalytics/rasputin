#pragma once

// Test-only: the exact point-in-ring classifier the CDT domain property
// (prop_cdt_invariants.cpp, "no triangle lies in a hole or outside the
// domain") uses as its oracle. It moved here verbatim from
// terrain/core/ring.hpp when production stopped calling it (C++ audit PR B,
// docs/increments/cpp-audit.md section 7); only the namespace and the ring
// parameter, now `const IndexedRing&`, changed.
//
// Call it qualified, terrain::test::point_in_ring<K>(r, p): an unqualified
// call would also look in namespace terrain through the ring's type.

#include <terrain/core/point.hpp>
#include <terrain/core/ring.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <cstddef>

namespace terrain::test {

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
template <pred::GeometryKernel K>
[[nodiscard]] PointInRing point_in_ring(const IndexedRing& r, const Point2& p) {
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

}  // namespace terrain::test
