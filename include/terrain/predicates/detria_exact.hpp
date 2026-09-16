#pragma once

// The exact-arithmetic backend behind `DefaultKernel`, implemented on top of
// the vendored detria predicates (`lib/detria`, MIT by election).
//
// This header deliberately does **not** include detria.hpp, and no header under
// include/ ever may. detria.hpp is 4500 lines, pulls in <iostream>, <sstream>
// and <csignal>, and its Debug assertions raise SIGTRAP; nothing in the engine
// should compile against all that just to ask which way three points turn. The
// declarations below are enough, and the definitions live in the single
// translation unit src/predicates/detria_exact.cpp.
//
// The rule is enforced structurally, not by convention: lib/detria is on the
// `terrain_predicates` target's include path PRIVATE, so an include from here
// would not even resolve for a consumer; tools/check_detria_boundary.py fails
// the governance job if one is added; and the unit suite carries an
// `#ifdef DETRIA_HPP_INCLUDED -> #error` guard.
//
// Because the definitions are out of line, this backend is not inlinable across
// the boundary. That is affordable: `FilteredKernel` answers well-separated
// input from double arithmetic alone and reaches the backend only near
// degeneracy, where a call is negligible against the expansion arithmetic it is
// about to do.

#include <terrain/core/point.hpp>
#include <terrain/predicates/orientation.hpp>

namespace terrain::pred {

// Stateless, static, and thread-safe: detria's robust predicates keep every
// floating-point expansion on the stack, compute their error bounds as
// constexpr, and have no lazily initialised tables -- unlike Shewchuk's
// original predicates.c, whose `exactinit()` globals were the reason the ports
// that wrap it were rejected. See lib/detria/README.md.
struct DetriaExact {
    [[nodiscard]] static Orientation orient2d(const Point2& a, const Point2& b,
                                              const Point2& c) noexcept;

    // Precondition: a, b, c is counterclockwise. This one is enforced by the
    // backend itself -- detria's `incircle` asserts it under `#ifndef NDEBUG`
    // and its assert handler calls std::raise(SIGTRAP) -- so a violation in a
    // Debug build kills the process on a signal rather than returning a wrong
    // answer. `FilteredKernel::incircle` establishes the precondition; nothing
    // else should call this directly without doing the same.
    [[nodiscard]] static Incircle incircle_ccw(const Point2& a, const Point2& b,
                                               const Point2& c, const Point2& d) noexcept;
};

}  // namespace terrain::pred
