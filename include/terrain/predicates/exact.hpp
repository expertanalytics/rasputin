#pragma once

// The exact-arithmetic backend interface, and the static error bounds the
// filtered kernel checks against it.
//
// A model of `ExactPredicates` answers the two determinant questions without
// approximation. It returns classifications rather than sign-carrying doubles:
// an exact backend has no residual magnitude worth exposing, and handing back a
// double invites callers to negate it, scale it, or compare it against a
// tolerance -- all of which are meaningless on an already-exact answer, and one
// of which (negation) is how the `incircle` normalization bug in kernel.hpp
// came about.
//
// No first-party exact arithmetic in this module, or in anything downstream of
// it, may use `__int128`. GCC 16 rejects it under `-Wpedantic -Werror`, which
// is this tree's warning posture and the compiler on the ubuntu CI job. Use
// limb-based or expansion-based arithmetic instead.

#include <terrain/core/point.hpp>
#include <terrain/predicates/orientation.hpp>

#include <concepts>
#include <limits>

namespace terrain::pred {

// Spelled with qualified static calls -- `E::orient2d(...)`, not
// `e.orient2d(...)` -- so every model is required to be usable without an
// instance. Predicates are pure functions; nothing should have to construct one
// to ask a question.
template <typename E>
concept ExactPredicates = requires(Point2 a, Point2 b, Point2 c, Point2 d) {
    { E::orient2d(a, b, c) } -> std::same_as<Orientation>;
    // The `_ccw` suffix is the precondition: a, b, c must be counterclockwise.
    // Callers normalize; see `FilteredKernel::incircle`.
    { E::incircle_ccw(a, b, c, d) } -> std::same_as<Incircle>;
};

// Half an ulp: the relative error bound of one correctly rounded double
// operation.
inline constexpr double epsilon = std::numeric_limits<double>::epsilon() / 2.0;

// Shewchuk's static (level A) error bounds for the 2D orientation and incircle
// determinants, in units of the corresponding permanent. If the absolute value
// of the naive determinant exceeds the bound times the permanent, its sign is
// certain and the exact backend is not needed.
//
// These two constants are public *only as a transcription check*: a typo in a
// magic constant of this kind does not change behaviour visibly -- it silently
// widens or disables the exact fallback -- so the suite pins them with
// STATIC_REQUIRE against independently written expressions. They are not a
// tuning knob and not an extension point; nothing outside the filter should
// branch on them.
inline constexpr double orient2d_bound_a = (3.0 + 16.0 * epsilon) * epsilon;
inline constexpr double incircle_bound_a = (10.0 + 96.0 * epsilon) * epsilon;

// The permanent expressions themselves stay private to the filter in
// kernel.hpp. They are an implementation detail of how the bound is applied,
// not part of any interface.

}  // namespace terrain::pred
