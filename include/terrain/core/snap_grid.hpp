#pragma once

// The snap grid: the lattice every coordinate in the noded output lies on, and
// the integral key coincidence is decided by.
//
// One scalar, anchored at zero. There is no offset member and there must never
// be one: the anchor is the entire reason `snap` is a pure function of `p`, and
// an input-dependent anchor would make adding one vertex re-snap every other
// one, turning four parallel-fors into a global reduction
// (`parallel_refinement.md:105-118`). The signature is the enforcement -- there
// is no argument through which an anchor could arrive.
//
// What snapping buys is NOT accuracy. `world(g)` is `index * spacing` in plain
// double and is not exact at a non-dyadic spacing; commit e412a43 settled that
// and nothing here re-derives it. What it buys is the one property the noder
// needs:
//
//     Two points snap to the same GridPoint if and only if their snapped world
//     coordinates are bit-identical doubles, because world() is a pure function
//     of the index.
//
// That is a statement about equality, not about accuracy, and it is what makes
// dedup legal. It needs `world` to be injective on the permitted index range,
// which is where kMaxGridIndex comes from -- see below.
//
// THIS HEADER IS KERNEL-FREE and must stay so. The cell *geometry* lives here;
// the cell *predicate* (`segment_meets_cell`) is noding's, because the question
// it answers is about a segment and this type knows nothing about segments.

#include <terrain/core/point.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

namespace terrain {

// The dedup key. Integral on purpose: std::int64_t has no NaN, so the ordering
// below is total and equality is exact. point.hpp's std::hash<Point2> is
// prohibited as a dedup key here (`parallel_refinement.md:126-128`) -- it
// cannot find a NaN and it hashes near-coincident points apart.
struct GridPoint {
    std::int64_t ix{};
    std::int64_t iy{};

    friend constexpr bool operator==(const GridPoint&, const GridPoint&) = default;

    // Defaulted, therefore lexicographic in (ix, iy) and total. The ordering is
    // not geometric and means nothing spatially; it exists so NodeSet can sort,
    // and so node ids are a function of the SET of points rather than of the
    // order a thread happened to find them in.
    friend constexpr auto operator<=>(const GridPoint&, const GridPoint&) = default;
};

// |index| must not exceed this, on either axis.
//
// For a normalized double, ulp(x) <= x * 2^-52, and each computed fl(k*s) sits
// within k*s*2^-53 of its exact value, so two neighbouring lattice coordinates
// stay distinct doubles while k + 1 < 2^52. 2^51 therefore holds a full factor
// of two of margin over a bound that is itself a bound, and it is nowhere near
// tight: a 1 mm grid over Web Mercator's full extent reaches 2^34.2.
inline constexpr std::int64_t kMaxGridIndex = std::int64_t{1} << 51;

// Finite and strictly positive, and nothing else -- a subnormal spacing is
// admitted, because whether a grid that fine is a sane thing to ask for is the
// driver's policy question rather than this type's.
//
// Free and constexpr for the reason increment 3 pulled out
// detail::sizes_fit_u32: a check that can only be exercised through the object
// that asserts on it is a check that rots. The driver calls this and returns a
// status; SnapGrid asserts it.
//
// Spelled `> 0.0`, never `!= 0.0`. Both reject +0.0 and both reject -0.0, since
// `-0.0 == 0.0` in IEEE. What `!= 0.0` wrongly ACCEPTS is every negative
// spacing, -inf included, whose grid runs backwards. Not NaN and not +inf: the
// `<= max()` conjunct rejects those whichever way the first one is spelled.
[[nodiscard]] constexpr bool is_valid_spacing(double spacing) noexcept {
    return spacing > 0.0 && spacing <= std::numeric_limits<double>::max();
}

class SnapGrid {
public:
    // Precondition: is_valid_spacing(spacing). Debug-asserted.
    explicit constexpr SnapGrid(double spacing) noexcept : spacing_{spacing} {
        assert(is_valid_spacing(spacing));
    }

    [[nodiscard]] constexpr double spacing() const noexcept { return spacing_; }

    // True iff snap(p) is representable and world() is injective there.
    // Precondition: p is finite -- Pslg guarantee 1, not re-checked.
    //
    // This is the driver's one-pass admission check, run once over the vertex
    // buffer, and the only place a coordinate is compared against
    // kMaxGridIndex. After that pass snap() is branch-free in the parallel-for.
    [[nodiscard]] bool can_snap(const Point2& p) const noexcept {
        constexpr double limit = static_cast<double>(kMaxGridIndex);
        return std::abs(p.x) / spacing_ <= limit && std::abs(p.y) / spacing_ <= limit;
    }

    // Precondition: can_snap(p). Debug-asserted, unchecked in release.
    //
    // std::llround, never std::rint or std::nearbyint: the latter two honour the
    // DYNAMIC rounding mode, so the same input can snap to different cells in
    // two threads or two builds. llround rounds half away from zero
    // unconditionally.
    //
    // Divide by the spacing; do not multiply by a precomputed reciprocal.
    // p.x * (1.0 / spacing_) is two roundings where p.x / spacing_ is one, and
    // the two disagree on inputs that are easy to find by search and impossible
    // to predict by reading.
    [[nodiscard]] GridPoint snap(const Point2& p) const noexcept {
        assert(can_snap(p));
        return GridPoint{static_cast<std::int64_t>(std::llround(p.x / spacing_)),
                         static_cast<std::int64_t>(std::llround(p.y / spacing_))};
    }

    // Total, and a pure function of g and spacing(). No precondition.
    [[nodiscard]] Point2 world(const GridPoint& g) const noexcept {
        return Point2{static_cast<double>(g.ix) * spacing_,
                      static_cast<double>(g.iy) * spacing_};
    }

    // world(snap(p)). Same precondition as snap.
    [[nodiscard]] Point2 snapped(const Point2& p) const noexcept { return world(snap(p)); }

    // The closed cell (hot pixel) of g: [(ix-1/2)s, (ix+1/2)s] x likewise in y.
    // Precondition: |ix|, |iy| <= kMaxGridIndex. Debug-asserted.
    //
    // ONE rounding, computed in index space: (2*ix +/- 1) * (spacing_ * 0.5),
    // never world(g) +/- spacing_ * 0.5. Halving is exact, 2*ix +/- 1 is exact
    // (kMaxGridIndex is 2^51 and std::int64_t holds 2^52 with ten bits to
    // spare), so a corner is a single fl() of an exact product and therefore a
    // pure function of g and spacing() in precisely the sense world() is. The
    // world-space spelling is two roundings on top of an already-rounded
    // coordinate, and the cells then fail to abut.
    [[nodiscard]] Point2 cell_min(const GridPoint& g) const noexcept {
        assert(in_range(g));
        const double half = spacing_ * 0.5;
        return Point2{static_cast<double>(2 * g.ix - 1) * half,
                      static_cast<double>(2 * g.iy - 1) * half};
    }

    [[nodiscard]] Point2 cell_max(const GridPoint& g) const noexcept {
        assert(in_range(g));
        const double half = spacing_ * 0.5;
        return Point2{static_cast<double>(2 * g.ix + 1) * half,
                      static_cast<double>(2 * g.iy + 1) * half};
    }

private:
    [[nodiscard]] static constexpr bool in_range(const GridPoint& g) noexcept {
        return g.ix >= -kMaxGridIndex && g.ix <= kMaxGridIndex && g.iy >= -kMaxGridIndex &&
               g.iy <= kMaxGridIndex;
    }

    double spacing_{};
};

}  // namespace terrain
