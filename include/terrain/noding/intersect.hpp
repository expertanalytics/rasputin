#pragma once

// Pairwise segment intersection for the noder: classification, then
// construction, never together -- plus the hot-pixel predicate.
//
// The split is the design. `classify` is a PREDICATE: it divides nothing,
// constructs nothing and has no tolerance, so under an exact kernel its answer
// is exact for every finite input, which is what makes it usable to decide
// topology. `crossing_point` is a CONSTRUCTION: it divides, it rounds, and its
// result is a GridPoint rather than a Point2 -- the type-level statement that
// the rounding target is the grid and that no un-snapped constructed coordinate
// may escape into the pipeline.
//
// Disjoint, Touching and Overlapping need no construction at all. Input
// vertices are snapped BEFORE pairwise testing, so every endpoint is already a
// grid point and every split point for a touch or an overlap is one of the four
// endpoints. Construction happens once per crossing pair and nowhere else.
//
// This header includes no `core/pslg.hpp` and knows nothing of chains: 5a is
// everything whose correctness is a question about numbers.

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/core/snap_grid.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace terrain::noding {

// The exact classification of two CLOSED segments. Exhaustive: every pair of
// finite segments is exactly one of these.
enum class SegmentRelation : std::uint8_t {
    Disjoint,     // no common point
    Touching,     // exactly one common point, and it is an endpoint of at least one
    Crossing,     // exactly one common point, interior to BOTH
    Overlapping,  // collinear, sharing a sub-segment of positive length
};

// Exact and construction-free: signs from K, betweenness from closed
// comparisons between coordinates the caller supplied. Total on every pair,
// including degenerate (zero-length) segments, which get no special case --
// segment.hpp's ruling, with its reason intact: "an explicit special case here
// is how the degenerate answer gets subtly wrong".
template <pred::GeometryKernel K>
[[nodiscard]] SegmentRelation classify(const Segment2& s, const Segment2& t) {
    using pred::Orientation;

    const Orientation o1 = K::orient2d(s.a, s.b, t.a);
    const Orientation o2 = K::orient2d(s.a, s.b, t.b);
    const Orientation o3 = K::orient2d(t.a, t.b, s.a);
    const Orientation o4 = K::orient2d(t.a, t.b, s.b);

    // The collinear arm. It must MEASURE the shared extent rather than return
    // Overlapping on sight: two collinear segments sharing exactly one endpoint
    // -- consecutive edges of a ring, the most common pair in the whole input
    // -- are Touching, and an arm that said Overlapping would merge every ring
    // edge with its neighbour.
    //
    // Both axes, not a dominant axis. Choosing "the axis along which the
    // segments vary" needs a branch that is wrong for a vertical pair, for a
    // horizontal pair and for a zero-length segment; intersecting both
    // intervals costs four comparisons and has no degenerate case.
    if (o1 == Orientation::Collinear && o2 == Orientation::Collinear &&
        o3 == Orientation::Collinear && o4 == Orientation::Collinear) {
        const double xlo = std::max(std::min(s.a.x, s.b.x), std::min(t.a.x, t.b.x));
        const double xhi = std::min(std::max(s.a.x, s.b.x), std::max(t.a.x, t.b.x));
        const double ylo = std::max(std::min(s.a.y, s.b.y), std::min(t.a.y, t.b.y));
        const double yhi = std::min(std::max(s.a.y, s.b.y), std::max(t.a.y, t.b.y));

        if (xlo > xhi || ylo > yhi) {
            return SegmentRelation::Disjoint;
        }
        if (xlo == xhi && ylo == yhi) {
            return SegmentRelation::Touching;
        }
        return SegmentRelation::Overlapping;
    }

    // The four-orientation arm: the segments meet at exactly one point. It is
    // Crossing iff that point is interior to both, which is exactly the case
    // where none of the four orientations is Collinear; if any is, the meeting
    // point is an endpoint of one of them.
    if (o1 != o2 && o3 != o4) {
        const bool any_collinear =
            o1 == Orientation::Collinear || o2 == Orientation::Collinear ||
            o3 == Orientation::Collinear || o4 == Orientation::Collinear;
        return any_collinear ? SegmentRelation::Touching : SegmentRelation::Crossing;
    }

    // The endpoint arm is a TOTALITY FALLBACK, not a discriminator, and this
    // comment is why it does not get deleted. Under an exact K it never returns
    // Touching -- reaching it needs o1 == o2 or o3 == o4, and both cases are
    // already decided above -- so its only reachable answer here is Disjoint.
    // It stays because classify is a template on any conforming K, including
    // FastKernel, whose mis-signed or wrongly-Collinear orientation is exactly
    // what drops a genuine meeting into this arm; without it the answer there
    // is the silent wrong one rather than the conservative one. It also keeps
    // exhaustiveness a property of the code rather than of a case analysis in a
    // document.
    if (on_segment<K>(s, t.a) || on_segment<K>(s, t.b) || on_segment<K>(t, s.a) ||
        on_segment<K>(t, s.b)) {
        return SegmentRelation::Touching;
    }
    return SegmentRelation::Disjoint;
}

// THE ONLY CONSTRUCTIVE FUNCTION IN THIS PROJECT.
// Precondition: classify<K>(s, t) == SegmentRelation::Crossing. Debug-asserted;
// K is taken for that assert and for nothing else, since a kernel decides signs
// and this function decides a value.
//
// THE CLAMP IS MANDATORY AND IS NOT A TOLERANCE. `den` is near zero for a
// near-parallel crossing and the computed u can then be off by orders of
// magnitude, putting p arbitrarily far from either segment. The true
// intersection provably lies in the intersection of the two segments'
// coordinate ranges, so clamping into that box bounds the damage to a region
// that provably contains the answer, using only comparisons between coordinates
// the caller supplied.
//
// The clamp is also what ESTABLISHES snap's precondition: after it, p lies
// inside a box whose corners are caller coordinates that already passed
// can_snap in the driver's one-pass check, so can_snap(p) holds by monotonicity
// of |p| / spacing. Remove the clamp and this function acquires an unstated
// precondition on its inputs' magnitudes -- and, in release, can index past
// kMaxGridIndex and break world() injectivity.
//
// What the clamp does NOT buy: it bounds the error, it does not make it small.
// For a near-parallel pair the clamped point can be anywhere in a long thin
// overlap region, and snapping then commits to whatever cell that is. That is
// the standing contract of snap rounding, not a defect of this function.
template <pred::GeometryKernel K>
[[nodiscard]] GridPoint crossing_point(const SnapGrid& grid, const Segment2& s,
                                       const Segment2& t) {
    assert(classify<K>(s, t) == SegmentRelation::Crossing);

    const Point2 d1 = s.b - s.a;
    const Point2 d2 = t.b - t.a;
    const double den = cross(d1, d2);  // nonzero: Crossing implies non-parallel
    const double u = cross(t.a - s.a, d2) / den;

    const double xlo = std::max(std::min(s.a.x, s.b.x), std::min(t.a.x, t.b.x));
    const double xhi = std::min(std::max(s.a.x, s.b.x), std::max(t.a.x, t.b.x));
    const double ylo = std::max(std::min(s.a.y, s.b.y), std::min(t.a.y, t.b.y));
    const double yhi = std::min(std::max(s.a.y, s.b.y), std::max(t.a.y, t.b.y));

    return grid.snap(Point2{std::clamp(s.a.x + u * d1.x, xlo, xhi),
                            std::clamp(s.a.y + u * d1.y, ylo, yhi)});
}

// The hot-pixel question: does the closed cell of g meet segment s?
//
// THIS IS NOT on_segment, AND THE DIFFERENCE IS THE POINT. on_segment asks
// EXACT INCIDENCE; snap rounding is defined on PROXIMITY -- each segment is
// routed through every hot pixel it passes through, not through every pixel
// whose centre it exactly contains. Those differ on the great majority of
// snapped input, because world(g) = g*spacing is not an affine map at a
// non-dyadic spacing, so a design that detects T-junctions through on_segment
// reports Disjoint for most real ones, never splits the host, and then agrees
// with its own error when the verification pass asks the same question.
//
// Separating-axis between a segment and an axis-aligned box: the candidate axes
// are the box's two face normals (step 1) and the segment's normal (step 2).
// Both halves are necessary and neither is sufficient -- dropping step 1 admits
// any cell the segment's LINE crosses however far along it, and dropping step 2
// admits any cell in the segment's bounding box, which for a long diagonal is
// most of the domain.
//
// Exact, total, division-free, tolerance-free, and expressible in the kernel
// exactly as it stands. No precondition beyond g being in range for
// cell_min/cell_max.
template <pred::GeometryKernel K>
[[nodiscard]] bool segment_meets_cell(const SnapGrid& grid, const Segment2& s,
                                      const GridPoint& g) {
    const Point2 lo = grid.cell_min(g);
    const Point2 hi = grid.cell_max(g);

    // 1. The box's face normals. Closed comparisons, no arithmetic.
    if (std::max(s.a.x, s.b.x) < lo.x || std::min(s.a.x, s.b.x) > hi.x) {
        return false;
    }
    if (std::max(s.a.y, s.b.y) < lo.y || std::min(s.a.y, s.b.y) > hi.y) {
        return false;
    }

    // 2. Does the LINE through s separate the cell? It does only when all four
    // corners lie strictly on the same side.
    //
    // "All the same NON-Collinear value", never "all strictly positive or all
    // strictly negative". A corner exactly on the line makes one orientation
    // Collinear, the cell is grazed, and the answer is true; written as a
    // strict sign test the predicate loses exactly the incidences it exists to
    // find. A zero-length s makes every orientation Collinear, which is why it
    // needs no branch of its own: the box test above has already decided it.
    const Point2 corners[] = {lo, Point2{hi.x, lo.y}, hi, Point2{lo.x, hi.y}};
    const pred::Orientation first = K::orient2d(s.a, s.b, corners[0]);
    if (first == pred::Orientation::Collinear) {
        return true;
    }
    for (std::size_t i = 1; i < 4; ++i) {
        if (K::orient2d(s.a, s.b, corners[i]) != first) {
            return true;
        }
    }
    return false;
}

}  // namespace terrain::noding
