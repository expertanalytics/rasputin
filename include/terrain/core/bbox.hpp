#pragma once

// The axis-aligned bounding box, and the free bounding_box over a point span.
//
// The load-bearing decision here is that a default-constructed Box2 is the
// EMPTY box -- lo at +inf, hi at -inf -- and not a degenerate box at the
// origin. That is what makes
//
//     Box2 b; for (const Point2& p : pts) b.expand(p);
//
// a correct fold: an empty input yields an empty box rather than one that
// silently contains (0, 0), and no caller has to special-case the first
// point. It is also why Box2 is a class rather than an aggregate -- an
// aggregate {Point2 lo; Point2 hi;} value-initialises to a box at the origin.
//
// The sentinel corners are chosen so that min/max against them is the identity,
// which is what makes expand(Box2) a no-op for an empty argument in either
// direction, and therefore what makes a fold over per-chunk partial boxes
// correct.
//
// There are no tolerances in this header and none may be added. Every
// comparison Box2 makes is between coordinates the caller supplied.
//
// This header deliberately depends on nothing from predicates/: Box2 is a
// filter over coordinates, never a topology decision, so it needs no kernel.

#include <terrain/core/point.hpp>

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <span>
#include <stdexcept>

namespace terrain {

class Box2 {
public:
    // The empty box. See the header comment: this is the identity of the union,
    // not a box at the origin.
    constexpr Box2() noexcept = default;

    // Throws std::invalid_argument on a non-finite corner or on lo > hi in
    // either component.
    //
    // The finiteness test is explicit rather than implied by `lo > hi`: every
    // comparison against NaN is false, so a range test alone accepts a NaN
    // corner, and the box then answers every containment query with a confident
    // `false`. RasterGeometry::cell_of documents the same hazard.
    //
    // This also means the empty sentinel is unreachable through this
    // constructor -- (+inf, +inf), (-inf, -inf) is rejected as non-finite. The
    // empty box is spelled Box2{} or Box2::empty(), and nothing else.
    Box2(const Point2& lo, const Point2& hi) : lo_{lo}, hi_{hi} {
        if (!std::isfinite(lo.x) || !std::isfinite(lo.y) ||
            !std::isfinite(hi.x) || !std::isfinite(hi.y)) {
            throw std::invalid_argument("terrain::Box2: corners must be finite");
        }
        if (lo.x > hi.x || lo.y > hi.y) {
            throw std::invalid_argument("terrain::Box2: lo must not exceed hi in either component");
        }
    }

    [[nodiscard]] static constexpr Box2 empty() noexcept { return Box2{}; }

    [[nodiscard]] constexpr bool is_empty() const noexcept {
        return lo_.x > hi_.x || lo_.y > hi_.y;
    }

    [[nodiscard]] constexpr const Point2& lo() const noexcept { return lo_; }
    [[nodiscard]] constexpr const Point2& hi() const noexcept { return hi_; }

    // Zero on the empty box, where hi - lo would be -inf.
    [[nodiscard]] constexpr double width() const noexcept {
        return is_empty() ? 0.0 : hi_.x - lo_.x;
    }

    [[nodiscard]] constexpr double height() const noexcept {
        return is_empty() ? 0.0 : hi_.y - lo_.y;
    }

    // Precondition: !is_empty().
    //
    // Computed as lo/2 + hi/2, NOT (lo + hi)/2. A box with finite, legal
    // corners can have a corner sum that overflows to infinity, and a box
    // symmetric about the origin hides that because lo + hi cancels to exactly
    // zero. Halving first is exact -- division by two is a lossless exponent
    // decrement on every finite double, subnormals included.
    [[nodiscard]] constexpr Point2 center() const noexcept {
        return lo_ / 2.0 + hi_ / 2.0;
    }

    // Precondition: p is finite.
    //
    // expand is the inner loop of every bounding-box computation in the engine,
    // so it validates nothing and is noexcept. Finiteness is established once,
    // at the PSLG and Python boundaries. A NaN reaching here corrupts the box
    // silently; that is the accepted cost of not branching per vertex. The free
    // bounding_box below is the validating entry point.
    constexpr void expand(const Point2& p) noexcept {
        lo_.x = std::min(lo_.x, p.x);
        lo_.y = std::min(lo_.y, p.y);
        hi_.x = std::max(hi_.x, p.x);
        hi_.y = std::max(hi_.y, p.y);
    }

    // Precondition: b is finite OR empty.
    //
    // Note the "or empty": the empty box's corners are +-inf, and expanding by
    // an empty sub-box must be a no-op. That falls out of min/max against the
    // sentinels with no branch, and it is exactly what makes a fold over
    // per-chunk partial boxes correct when some chunks were empty.
    constexpr void expand(const Box2& b) noexcept {
        lo_.x = std::min(lo_.x, b.lo_.x);
        lo_.y = std::min(lo_.y, b.lo_.y);
        hi_.x = std::max(hi_.x, b.hi_.x);
        hi_.y = std::max(hi_.y, b.hi_.y);
    }

    // Closed and exact: one ulp outside is outside. A non-finite query is not
    // contained -- for NaN because every comparison against it is false, for
    // an infinity because it lies outside any finite box. The empty box
    // contains no point, since +inf <= p is false for every finite p.
    [[nodiscard]] constexpr bool contains(const Point2& p) const noexcept {
        return lo_.x <= p.x && p.x <= hi_.x && lo_.y <= p.y && p.y <= hi_.y;
    }

    // Closed and reflexive. The empty box is the identity of the union and so a
    // subset of every box, including of itself; that is the early return, and
    // it cannot be left to the coordinate comparisons, which would answer false
    // for the +-inf sentinels.
    [[nodiscard]] constexpr bool contains(const Box2& b) const noexcept {
        if (b.is_empty()) {
            return true;
        }
        return lo_.x <= b.lo_.x && b.hi_.x <= hi_.x && lo_.y <= b.lo_.y && b.hi_.y <= hi_.y;
    }

    // Closed: boxes touching at an edge or a corner intersect. The empty box
    // intersects nothing, including itself -- it has no points to share.
    [[nodiscard]] constexpr bool intersects(const Box2& b) const noexcept {
        if (is_empty() || b.is_empty()) {
            return false;
        }
        return lo_.x <= b.hi_.x && b.lo_.x <= hi_.x && lo_.y <= b.hi_.y && b.lo_.y <= hi_.y;
    }

    friend constexpr bool operator==(const Box2&, const Box2&) = default;

private:
    static constexpr double inf_ = std::numeric_limits<double>::infinity();

    Point2 lo_{inf_, inf_};
    Point2 hi_{-inf_, -inf_};
};

// The validating entry point. Unlike expand, this is what a caller reaches for
// with data of unknown provenance, so it throws std::invalid_argument on a
// non-finite coordinate. An empty span yields the empty box.
[[nodiscard]] inline Box2 bounding_box(std::span<const Point2> points) {
    Box2 box;
    for (const Point2& p : points) {
        if (!std::isfinite(p.x) || !std::isfinite(p.y)) {
            throw std::invalid_argument("terrain::bounding_box: every point must be finite");
        }
        box.expand(p);
    }
    return box;
}

}  // namespace terrain

namespace std {

template <>
struct formatter<terrain::Box2> {
    // Same reasoning as point.hpp's specialisation: a spec that parses and is
    // then discarded is a silent lie about the output, so reject what we do not
    // honour.
    constexpr auto parse(format_parse_context& ctx) {
        auto it = ctx.begin();
        if (it != ctx.end() && *it != '}')
            throw format_error("terrain::Box2 does not accept a format spec");
        return it;
    }

    // The empty box renders as "Box2(empty)" rather than through its corners.
    // Printing "Box2(Point2(inf, inf), Point2(-inf, -inf))" hands a reader what
    // reads as a bug report rather than as the identity element it is.
    template <typename FormatContext>
    auto format(const terrain::Box2& b, FormatContext& ctx) const {
        if (b.is_empty()) {
            return std::format_to(ctx.out(), "Box2(empty)");
        }
        return std::format_to(ctx.out(), "Box2({}, {})", b.lo(), b.hi());
    }
};

}  // namespace std
