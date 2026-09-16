#pragma once

// Test-only reference implementations of the exact predicate backend.
//
// These exist so the predicates suite can be written and run before
// `terrain::pred::DetriaExact` (and the vendored detria.hpp behind it) exist,
// and so it keeps a second, independently derived oracle afterwards. A filter
// checked only against the backend it wraps proves nothing about the backend.
//
// `RefExact` is deliberately unsophisticated: no adaptivity, no expansions, no
// error bounds. It restricts itself to integer-valued coordinates and evaluates
// the determinants in exact fixed-width integer arithmetic, so its answers are
// correct by construction rather than by a numerical argument. Integer-valued
// coordinates are also the family that matters most operationally, since that
// is what a snap grid produces.
//
// `__int128` is not used: GCC rejects it under `-Wpedantic -Werror`, which this
// tree compiles with. The 128-bit arithmetic below is four 32-bit limbs.

#include <terrain/core/point.hpp>
#include <terrain/predicates/orientation.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace terrain::test {

// ---------------------------------------------------------------------------
// Minimal unsigned 128-bit magnitude: four base-2^32 limbs, little endian.
// ---------------------------------------------------------------------------

struct U128 {
    std::array<std::uint64_t, 4> limb{};
};

// Adds `value` (any 64-bit quantity) starting at limb `idx`, propagating carry.
// Overflow throws rather than wrapping: an oracle that silently returns a wrong
// sign is worse than no oracle.
inline void u128_add_at(U128& r, std::size_t idx, std::uint64_t value) {
    while (value != 0) {
        if (idx >= r.limb.size()) {
            throw std::overflow_error("RefExact: 128-bit reference arithmetic overflowed");
        }
        const std::uint64_t sum = r.limb[idx] + (value & 0xFFFFFFFFULL);
        r.limb[idx] = sum & 0xFFFFFFFFULL;
        value = (value >> 32) + (sum >> 32);
        ++idx;
    }
}

[[nodiscard]] inline U128 u128_from_u64(std::uint64_t v) {
    U128 r{};
    r.limb[0] = v & 0xFFFFFFFFULL;
    r.limb[1] = v >> 32;
    return r;
}

[[nodiscard]] inline U128 u128_mul_u64(const U128& a, std::uint64_t b) {
    const std::array<std::uint64_t, 2> bs{b & 0xFFFFFFFFULL, b >> 32};
    U128 r{};
    for (std::size_t i = 0; i < a.limb.size(); ++i) {
        for (std::size_t j = 0; j < bs.size(); ++j) {
            const std::uint64_t product = a.limb[i] * bs[j];
            if (product == 0) {
                continue;
            }
            u128_add_at(r, i + j, product);
        }
    }
    return r;
}

[[nodiscard]] inline int u128_cmp(const U128& a, const U128& b) {
    for (std::size_t i = a.limb.size(); i-- > 0;) {
        if (a.limb[i] != b.limb[i]) {
            return a.limb[i] < b.limb[i] ? -1 : 1;
        }
    }
    return 0;
}

[[nodiscard]] inline U128 u128_add(const U128& a, const U128& b) {
    U128 r = a;
    for (std::size_t i = 0; i < b.limb.size(); ++i) {
        u128_add_at(r, i, b.limb[i]);
    }
    return r;
}

// Precondition: a >= b.
[[nodiscard]] inline U128 u128_sub(const U128& a, const U128& b) {
    U128 r{};
    std::uint64_t borrow = 0;
    for (std::size_t i = 0; i < a.limb.size(); ++i) {
        const std::uint64_t lhs = a.limb[i];
        const std::uint64_t rhs = b.limb[i] + borrow;
        if (lhs >= rhs) {
            r.limb[i] = lhs - rhs;
            borrow = 0;
        } else {
            r.limb[i] = (lhs + (1ULL << 32)) - rhs;
            borrow = 1;
        }
    }
    return r;
}

// ---------------------------------------------------------------------------
// Signed 128-bit value as sign plus magnitude.
// ---------------------------------------------------------------------------

struct I128 {
    int sign{0};  // -1, 0 or +1
    U128 mag{};
};

[[nodiscard]] inline I128 i128_normalized(int sign, const U128& mag) {
    const U128 zero{};
    if (u128_cmp(mag, zero) == 0) {
        return I128{};
    }
    return I128{sign, mag};
}

[[nodiscard]] inline I128 i128_from_i64(std::int64_t v) {
    if (v == 0) {
        return I128{};
    }
    const std::uint64_t mag = v < 0 ? (~static_cast<std::uint64_t>(v) + 1ULL)
                                    : static_cast<std::uint64_t>(v);
    return I128{v < 0 ? -1 : 1, u128_from_u64(mag)};
}

[[nodiscard]] inline I128 i128_mul_i64(const I128& a, std::int64_t b) {
    if (a.sign == 0 || b == 0) {
        return I128{};
    }
    const std::uint64_t mag = b < 0 ? (~static_cast<std::uint64_t>(b) + 1ULL)
                                    : static_cast<std::uint64_t>(b);
    return i128_normalized(b < 0 ? -a.sign : a.sign, u128_mul_u64(a.mag, mag));
}

[[nodiscard]] inline I128 i128_mul(std::int64_t a, std::int64_t b) {
    return i128_mul_i64(i128_from_i64(a), b);
}

[[nodiscard]] inline I128 i128_add(const I128& a, const I128& b) {
    if (a.sign == 0) {
        return b;
    }
    if (b.sign == 0) {
        return a;
    }
    if (a.sign == b.sign) {
        return i128_normalized(a.sign, u128_add(a.mag, b.mag));
    }
    const int c = u128_cmp(a.mag, b.mag);
    if (c == 0) {
        return I128{};
    }
    return c > 0 ? i128_normalized(a.sign, u128_sub(a.mag, b.mag))
                 : i128_normalized(b.sign, u128_sub(b.mag, a.mag));
}

[[nodiscard]] inline I128 i128_neg(const I128& a) {
    return I128{-a.sign, a.mag};
}

[[nodiscard]] inline I128 i128_sub(const I128& a, const I128& b) {
    return i128_add(a, i128_neg(b));
}

// ---------------------------------------------------------------------------
// RefExact: exact predicates over integer-valued coordinates.
// ---------------------------------------------------------------------------

// `orient2d` needs products of coordinate differences (2 factors), `incircle`
// needs triple products against a lifted coordinate (~4 factors), so the two
// admit different coordinate ranges before the 128-bit accumulator overflows.
// Both bounds comfortably cover UTM33 easting/northing magnitudes.
inline constexpr double ref_orient2d_coord_bound = 4503599627370496.0;  // 2^52
inline constexpr double ref_incircle_coord_bound = 33554432.0;          // 2^25

[[nodiscard]] inline std::int64_t exact_integer(double v, double bound, const char* what) {
    if (!(std::isfinite(v)) || v != std::floor(v) || std::abs(v) > bound) {
        throw std::invalid_argument(
            std::string{"RefExact: "} + what + " requires integer-valued coordinates within its bound");
    }
    return static_cast<std::int64_t>(v);
}

struct RefExact {
    [[nodiscard]] static pred::Orientation orient2d(const Point2& a, const Point2& b, const Point2& c) {
        const std::int64_t ax = exact_integer(a.x, ref_orient2d_coord_bound, "orient2d");
        const std::int64_t ay = exact_integer(a.y, ref_orient2d_coord_bound, "orient2d");
        const std::int64_t bx = exact_integer(b.x, ref_orient2d_coord_bound, "orient2d");
        const std::int64_t by = exact_integer(b.y, ref_orient2d_coord_bound, "orient2d");
        const std::int64_t cx = exact_integer(c.x, ref_orient2d_coord_bound, "orient2d");
        const std::int64_t cy = exact_integer(c.y, ref_orient2d_coord_bound, "orient2d");

        const I128 det = i128_sub(i128_mul(bx - ax, cy - ay), i128_mul(by - ay, cx - ax));
        return static_cast<pred::Orientation>(det.sign);
    }

    // Mathematically total: the 3x3 lifted determinant is defined for any four
    // points. The `_ccw` suffix documents the caller's obligation, and
    // `CcwCheckingExact` below is what actually polices it.
    [[nodiscard]] static pred::Incircle incircle_ccw(const Point2& a, const Point2& b,
                                               const Point2& c, const Point2& d) {
        const auto coord = [](double v) {
            return exact_integer(v, ref_incircle_coord_bound, "incircle_ccw");
        };
        const std::int64_t dx = coord(d.x);
        const std::int64_t dy = coord(d.y);

        const std::int64_t adx = coord(a.x) - dx;
        const std::int64_t ady = coord(a.y) - dy;
        const std::int64_t bdx = coord(b.x) - dx;
        const std::int64_t bdy = coord(b.y) - dy;
        const std::int64_t cdx = coord(c.x) - dx;
        const std::int64_t cdy = coord(c.y) - dy;

        // Each lift is at most 2^53 and each 2x2 cross at most 2^53, so both fit
        // exactly in int64; only the triple products need the 128-bit path.
        const std::int64_t alift = adx * adx + ady * ady;
        const std::int64_t blift = bdx * bdx + bdy * bdy;
        const std::int64_t clift = cdx * cdx + cdy * cdy;

        // det = | adx ady alift ; bdx bdy blift ; cdx cdy clift |, expanded
        // along the first row. Positive iff d lies inside the circle through
        // a, b, c when a, b, c is counterclockwise.
        const I128 t0 = i128_mul_i64(i128_sub(i128_mul(bdy, clift), i128_mul(blift, cdy)), adx);
        const I128 t1 = i128_mul_i64(i128_sub(i128_mul(bdx, clift), i128_mul(blift, cdx)), ady);
        const I128 t2 = i128_mul_i64(i128_sub(i128_mul(bdx, cdy), i128_mul(bdy, cdx)), alift);

        const I128 det = i128_add(i128_sub(t0, t1), t2);
        return static_cast<pred::Incircle>(det.sign);
    }
};


// ---------------------------------------------------------------------------
// Instrumenting wrappers.
// ---------------------------------------------------------------------------

// Counts how often the filter actually falls back. A filter that never defers
// to the exact backend is wrong, and one that always defers is pointless; only
// a call count distinguishes the two from the outside.
//
// The counters are static because `ExactPredicates` requires static member
// functions, so the wrapper has nowhere else to put them. That makes this type
// unusable from concurrent tests -- use `RefExact` directly there.
template <typename E>
struct CountingExact {
    static inline int orient2d_calls = 0;
    static inline int incircle_calls = 0;

    static void reset() {
        orient2d_calls = 0;
        incircle_calls = 0;
    }

    [[nodiscard]] static pred::Orientation orient2d(const Point2& a, const Point2& b, const Point2& c) {
        ++orient2d_calls;
        return E::orient2d(a, b, c);
    }

    [[nodiscard]] static pred::Incircle incircle_ccw(const Point2& a, const Point2& b,
                                               const Point2& c, const Point2& d) {
        ++incircle_calls;
        return E::incircle_ccw(a, b, c, d);
    }
};

// Stands in for detria's debug-build precondition check.
//
// detria's `incircle` asserts under `#ifndef NDEBUG` that its first three
// arguments are counterclockwise, and its assert handler raises SIGTRAP. In the
// asan+ubsan Debug CI job an un-normalized call would therefore kill the
// process on a signal instead of failing a test, which is unreadable as a
// failure report. This wrapper records the same violation as an ordinary flag
// so the suite can assert on it, and it must stay in the suite after
// `DetriaExact` lands: it is the only thing standing between a regression in
// `FilteredKernel::incircle`'s normalization and an unreproducible CI death.
template <typename E>
struct CcwCheckingExact {
    static inline int violations = 0;
    static inline int incircle_calls = 0;

    static void reset() {
        violations = 0;
        incircle_calls = 0;
    }

    [[nodiscard]] static pred::Orientation orient2d(const Point2& a, const Point2& b, const Point2& c) {
        return E::orient2d(a, b, c);
    }

    [[nodiscard]] static pred::Incircle incircle_ccw(const Point2& a, const Point2& b,
                                               const Point2& c, const Point2& d) {
        ++incircle_calls;
        if (E::orient2d(a, b, c) != pred::Orientation::CounterClockwise) {
            ++violations;
        }
        return E::incircle_ccw(a, b, c, d);
    }
};

}  // namespace terrain::test
