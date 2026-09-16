#pragma once

// Normalizing a triple to counterclockwise, using the reference oracle.
//
// `ExactPredicates::incircle_ccw` carries a precondition, and with a real
// backend behind it that precondition is not advisory: detria's `incircle`
// asserts it under `#ifndef NDEBUG` and its assert handler raises SIGTRAP. Any
// test that calls a backend's `incircle_ccw` *directly* -- bypassing
// `FilteredKernel::incircle`, which is the thing that normally normalizes --
// must therefore establish the precondition itself, or it is not testing the
// backend, it is arranging for the Debug CI job to die on a signal.
//
// The normalization here deliberately goes through `RefExact`, the independent
// oracle, and never through the kernel or the backend under test. A helper that
// asked the subject under test which way its own arguments turn would make
// every downstream agreement test vacuous on exactly the inputs where the two
// implementations disagree.

#include <exact_reference.hpp>

#include <terrain/core/point.hpp>
#include <terrain/predicates/orientation.hpp>

#include <optional>

namespace terrain::test {

struct CcwTriple {
    Point2 a;
    Point2 b;
    Point2 c;
};

// Returns the triple in counterclockwise order, or nullopt if it is collinear.
//
// Collinear input has no counterclockwise ordering at all, so there is nothing
// to hand a backend and the caller must skip it -- which is the same conclusion
// `FilteredKernel::incircle` reaches when it answers Cocircular without
// consulting the backend.
[[nodiscard]] inline std::optional<CcwTriple> as_ccw(const Point2& a, const Point2& b,
                                                     const Point2& c) {
    switch (RefExact::orient2d(a, b, c)) {
        case pred::Orientation::CounterClockwise:
            return CcwTriple{a, b, c};
        case pred::Orientation::Clockwise:
            // Swap, never negate the eventual answer: the circle through three
            // points does not depend on their order. See the note on
            // `reversed` in orientation.hpp.
            return CcwTriple{a, c, b};
        case pred::Orientation::Collinear:
            break;
    }
    return std::nullopt;
}

}  // namespace terrain::test
