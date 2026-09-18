#pragma once

// The per-edge feature property set: one opaque 32-bit word, and no feature
// name anywhere in terrain::.
//
// A constraint edge may carry several properties at once -- at a coarse
// resolution the same line segment is both a road and a river -- so the carrier
// is a SET and the merge is a union. Union is commutative, associative and
// idempotent, which is what lets the noder reduce over an UNORDERED set of
// contributing chains: the answer may not depend on the order the broad phase
// happens to visit buckets in, and a priority scheme has none of the three.
//
// WHAT THIS TYPE DELIBERATELY DOES NOT HAVE, and must never grow:
//
//   * A named member, an enumerator, a `River`. The mapping from bit position
//     to feature name is a boundary concern and lives in one Pydantic model,
//     src_python/tin_engine/features.py, alongside CRS metadata and everything
//     else the Python layer keeps out of C++. A refinement or hydrology policy
//     that wants "refine harder near rivers" is HANDED a mask by its caller and
//     is never compiled against a vocabulary. That is the cost of the firewall.
//   * Any conversion to or from bool, int or std::uint32_t, in either
//     direction. bits() is the one way out and bit()/operator| the one way in.
//     Without that, `add_chain(idx, role, true)` would keep compiling with a
//     changed meaning, and `if (chain.properties)` would compile at all.
//
// 32 is a ceiling named rather than discovered: one word per edge, trivially
// parallel-reducible, no allocation. The vocabulary it must hold is LINEAR
// features -- river, road, railway, coastline, contour, wall, ditch -- which is
// under ten, against the 23 and 44 class land-cover enumerations in legacy/,
// and land cover stays face-based. Widening to 64 is a one-line change to a
// type nobody can pattern-match on, because it has no named members.
//
// See docs/increments/07-edge-properties.md for the rulings above.

#include <cassert>
#include <cstdint>

namespace terrain {

class EdgeProperties {
public:
    static constexpr unsigned kMaxProperties = 32;

    // The empty set: an unclassified constraint. Legal everywhere, and the
    // default -- not an error, not a diagnostic, not a warning.
    constexpr EdgeProperties() = default;

    // The singleton {i}.
    //
    // Precondition: i < kMaxProperties. Debug-asserted, unchecked in release --
    // SnapGrid::can_snap/snap's shape (snap_grid.hpp:94-96, :114). Admission
    // happens once where untrusted data arrives (the mask range check in the
    // bindings); everything after it is branch-free.
    [[nodiscard]] static constexpr EdgeProperties bit(unsigned i) noexcept {
        assert(i < kMaxProperties);
        return EdgeProperties{FromBits{}, std::uint32_t{1} << i};
    }

    [[nodiscard]] constexpr bool empty() const noexcept { return bits_ == 0u; }

    // Superset, NOT "shares a bit": a.contains(b) iff every member of b is a
    // member of a. The intersection spelling would make a singleton contain the
    // pair it belongs to.
    [[nodiscard]] constexpr bool contains(EdgeProperties other) const noexcept {
        return (bits_ & other.bits_) == other.bits_;
    }

    [[nodiscard]] constexpr std::uint32_t bits() const noexcept { return bits_; }

    [[nodiscard]] friend constexpr EdgeProperties operator|(EdgeProperties a,
                                                            EdgeProperties b) noexcept {
        return EdgeProperties{FromBits{}, a.bits_ | b.bits_};
    }

    [[nodiscard]] friend constexpr bool operator==(const EdgeProperties&,
                                                   const EdgeProperties&) = default;

private:
    // Private, and two-argument, so that no conversion from a word exists even
    // by accident: std::is_constructible_v<EdgeProperties, std::uint32_t> is
    // false on both counts.
    struct FromBits {};

    constexpr EdgeProperties(FromBits, std::uint32_t bits) noexcept : bits_{bits} {}

    std::uint32_t bits_{};
};

}  // namespace terrain
