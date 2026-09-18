#pragma once

// The noder's node table: every distinct grid point in the noded input, sorted
// once, with node ids as indices into that sorted sequence.
//
// A sorted std::vector with a binary search, not a hash map, for the reasons
// increment 4 gave for ConstraintEdgeSet and which apply unchanged: the key is
// exact and small so hashing buys nothing; one allocation of known size beats a
// rehashing container; the structure is read-only after construction and is
// therefore shareable by every thread without synchronisation; and sorting
// dedups for free. This header includes no <unordered_map> and no
// <unordered_set>, deliberately.
//
// NODE IDS ARE IN LEXICOGRAPHIC GRID ORDER, NOT FIRST-APPEARANCE ORDER, and
// that is a ruling with a consequence. Lexicographic order is a function of the
// SET of points and of nothing else -- not of input order, not of which thread
// found a crossing first, not of the order the broad phase visited buckets in
// -- so the node numbering, and therefore every downstream index in the mesh,
// is reproducible across runs and across thread counts. First-appearance order
// would make the numbering depend on scheduling.
//
// Dedup here is exact integer equality and nothing else. Whether two nearby
// points are coincident was decided upstream, by SnapGrid; this type only
// collapses keys that are already equal.

#include <terrain/core/snap_grid.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <span>
#include <utility>
#include <vector>

namespace terrain::noding {

class NodeSet {
public:
    // Sorts and uniques. Takes by value and consumes.
    explicit NodeSet(std::vector<GridPoint> points) : points_{std::move(points)} {
        std::ranges::sort(points_);
        const auto duplicates = std::ranges::unique(points_);
        points_.erase(duplicates.begin(), duplicates.end());
    }

    [[nodiscard]] std::size_t size() const noexcept { return points_.size(); }

    // Sorted and unique.
    [[nodiscard]] std::span<const GridPoint> points() const noexcept { return points_; }

    // Precondition: g was among the constructor's arguments. Debug-asserted.
    [[nodiscard]] std::uint32_t id_of(const GridPoint& g) const noexcept {
        const auto it = std::ranges::lower_bound(points_, g);
        assert(it != points_.end() && *it == g);
        return static_cast<std::uint32_t>(std::distance(points_.begin(), it));
    }

    [[nodiscard]] const GridPoint& operator[](std::uint32_t id) const noexcept {
        assert(id < points_.size());
        return points_[id];
    }

private:
    std::vector<GridPoint> points_;
};

}  // namespace terrain::noding
