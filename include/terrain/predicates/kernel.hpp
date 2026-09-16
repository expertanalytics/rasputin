#pragma once

// The geometry kernels: the objects the triangulation asks geometric questions
// of.
//
// `FastKernel` evaluates the determinants in plain double arithmetic.
// `FilteredKernel<E>` evaluates the same determinants, checks the result
// against a static error bound, and defers to the exact backend `E` only when
// the sign is not certain. Both are empty types with static member functions:
// no instance state, no caching or memoizing of orientations at any scope, no
// lazy initialisation and no mutable globals, so every predicate is callable
// concurrently from any number of threads.

#include <terrain/core/point.hpp>
#include <terrain/predicates/exact.hpp>
#include <terrain/predicates/orientation.hpp>

#include <cmath>

namespace terrain::pred {

// Qualified static calls, for the same reason as `ExactPredicates`: a model
// must be usable without an instance.
//
// This is a semantic contract, not just a signature. Every model answers the
// same two questions with the same meaning -- in particular `incircle` is total
// and insensitive to the order of its first three points. A model is permitted
// to be *wrong* near degeneracy; it is not permitted to answer a different
// question.
template <typename K>
concept GeometryKernel = requires(Point2 a, Point2 b, Point2 c, Point2 d) {
    { K::orient2d(a, b, c) } -> std::same_as<Orientation>;
    { K::incircle(a, b, c, d) } -> std::same_as<Incircle>;
};

namespace detail {

// The naive 3x3 lifted determinant, translated so d is the origin. Positive iff
// d lies inside the circle through a, b, c when a, b, c is counterclockwise.
[[nodiscard]] inline double incircle_det(const Point2& a, const Point2& b,
                                         const Point2& c, const Point2& d) noexcept {
    const double adx = a.x - d.x;
    const double ady = a.y - d.y;
    const double bdx = b.x - d.x;
    const double bdy = b.y - d.y;
    const double cdx = c.x - d.x;
    const double cdy = c.y - d.y;

    const double alift = adx * adx + ady * ady;
    const double blift = bdx * bdx + bdy * bdy;
    const double clift = cdx * cdx + cdy * cdy;

    return alift * (bdx * cdy - cdx * bdy)
         + blift * (cdx * ady - adx * cdy)
         + clift * (adx * bdy - bdx * ady);
}

// Shewchuk's permanent for the incircle determinant: the same expression with
// every subtraction of products replaced by a sum of magnitudes. Kept private
// to the filter -- it is an implementation detail of how `incircle_bound_a` is
// applied, not part of any interface.
[[nodiscard]] inline double incircle_permanent(const Point2& a, const Point2& b,
                                               const Point2& c, const Point2& d) noexcept {
    const double adx = a.x - d.x;
    const double ady = a.y - d.y;
    const double bdx = b.x - d.x;
    const double bdy = b.y - d.y;
    const double cdx = c.x - d.x;
    const double cdy = c.y - d.y;

    const double alift = adx * adx + ady * ady;
    const double blift = bdx * bdx + bdy * bdy;
    const double clift = cdx * cdx + cdy * cdy;

    return (std::fabs(bdx * cdy) + std::fabs(cdx * bdy)) * alift
         + (std::fabs(cdx * ady) + std::fabs(adx * cdy)) * blift
         + (std::fabs(adx * bdy) + std::fabs(bdx * ady)) * clift;
}

}  // namespace detail

// Unfiltered, no fallback. It exists so the suite can demonstrate the wrong
// answers that motivate `FilteredKernel`, and as the baseline the filter's
// overhead is measured against. Do not use it for topology decisions.
struct FastKernel {
    // Expressed through `terrain::cross` rather than as a second, independently
    // written determinant, so the two can never drift apart.
    [[nodiscard]] static Orientation orient2d(const Point2& a, const Point2& b,
                                              const Point2& c) noexcept {
        return orientation_of_sign(cross(b - a, c - a));
    }

    // Normalizes exactly as `FilteredKernel::incircle` does, using its own
    // unfiltered `orient2d` -- being unfiltered end to end is its defining
    // property, so it never consults an exact orientation. Without this the two
    // kernels would answer different questions and `GeometryKernel` would be a
    // name for a signature rather than a contract.
    [[nodiscard]] static Incircle incircle(const Point2& a, const Point2& b,
                                           const Point2& c, const Point2& d) noexcept {
        switch (orient2d(a, b, c)) {
            case Orientation::Collinear:
                return Incircle::Cocircular;
            case Orientation::CounterClockwise:
                return incircle_of_sign(detail::incircle_det(a, b, c, d));
            case Orientation::Clockwise:
                break;
        }
        return incircle_of_sign(detail::incircle_det(a, c, b, d));
    }
};

// The kernel the rest of the engine should use. The filter decides *how* an
// answer is reached, never *what* it is: on every input, `FilteredKernel<E>`
// returns what `E` would.
template <ExactPredicates E>
struct FilteredKernel {
    [[nodiscard]] static Orientation orient2d(const Point2& a, const Point2& b,
                                              const Point2& c) {
        const double detleft = (b.x - a.x) * (c.y - a.y);
        const double detright = (b.y - a.y) * (c.x - a.x);
        const double det = detleft - detright;
        const double permanent = std::fabs(detleft) + std::fabs(detright);

        if (std::fabs(det) > orient2d_bound_a * permanent) {
            return orientation_of_sign(det);
        }
        return E::orient2d(a, b, c);
    }

    // Total and precondition-free, unlike the backend's `incircle_ccw`.
    //
    // The lifted determinant's sign is uninterpretable without knowing the
    // orientation of a, b, c, so the orientation is settled *first*; there is no
    // filter-first arrangement of this predicate. That costs roughly 1.2-1.3x on
    // the filter-passing path, which is the price of a predicate that cannot be
    // called wrongly.
    //
    // The orientation used is the *filtered* one, never `E::orient2d`: the
    // filter answers it outright on well-separated input, so the backend is not
    // touched by either predicate there.
    //
    // For a clockwise triple, b and c are swapped and there is no `reversed`
    // anywhere. Swapping makes the triple counterclockwise, and the circle
    // through three points does not depend on their order, so the reordered
    // call already returns the right answer; negating it would invert a correct
    // one. An earlier draft did exactly that, inherited from a version whose
    // backend returned a sign-carrying double.
    [[nodiscard]] static Incircle incircle(const Point2& a, const Point2& b,
                                           const Point2& c, const Point2& d) {
        switch (orient2d(a, b, c)) {
            // Three collinear points define a degenerate circle and every point
            // lies on it. This is a definition rather than an approximation, so
            // it is answered here: no lifted determinant is evaluated, and the
            // backend is not called at all -- there is no counterclockwise
            // triple it could legally be handed.
            case Orientation::Collinear:
                return Incircle::Cocircular;
            case Orientation::CounterClockwise:
                return incircle_ccw(a, b, c, d);
            case Orientation::Clockwise:
                break;
        }
        return incircle_ccw(a, c, b, d);
    }

private:
    // Precondition: a, b, c is counterclockwise. Private, because that
    // precondition is established by `incircle` and by nothing else.
    [[nodiscard]] static Incircle incircle_ccw(const Point2& a, const Point2& b,
                                               const Point2& c, const Point2& d) {
        const double det = detail::incircle_det(a, b, c, d);
        const double permanent = detail::incircle_permanent(a, b, c, d);

        if (std::fabs(det) > incircle_bound_a * permanent) {
            return incircle_of_sign(det);
        }
        return E::incircle_ccw(a, b, c, d);
    }
};

}  // namespace terrain::pred
