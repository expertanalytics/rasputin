#pragma once

// The noder's broad phase: a uniform bucket index over segment bounding boxes,
// with exactly one query.
//
// THE CONTRACT IS ONE-DIRECTIONAL AND IS NOT WEAKENED ANYWHERE:
//
//     for_each_candidate(q, f) invokes f once for every indexed segment whose
//     CLOSED bounding box intersects q, and may invoke it for others.
//
// False positives are free; a false negative is a bug in a function whose entire
// job is not to have them -- the noder would then never compare a pair that
// crosses, and a verification pass sharing the index would agree with it.
//
// THERE IS NO Point2 OVERLOAD AND THERE MUST NEVER BE ONE. The caller that would
// reach for it is the node query in node.hpp's step 4, and the right argument
// there is the node's CELL BOX, never world(g): a cell straddling a bucket
// boundary loses every segment on the far side when queried by its centre, and a
// T-junction with them.
//
// BUCKETING IS BY BOUNDING BOX, NOT BY EXACT TRAVERSAL, and that is what makes
// the postcondition provable in one line rather than by case analysis: a segment
// is inserted into every bucket its closed bbox overlaps, so a query box that
// meets that bbox shares at least one bucket with it, both being rasterised by
// the same overlap rule over the same lattice. An exact DDA traversal would be
// tighter and would cost that sentence.
//
// BUCKET SIZING IS BY SEGMENT COUNT AND DOMAIN EXTENT, NEVER BY
// SnapGrid::spacing(). k = max(1, ceil(sqrt(n))) buckets per axis over the
// segment set's bounding box, so occupancy is O(1) for uniformly distributed
// segments and the index is O(n) in memory. A 100 km domain at a 5 cm spacing
// would want 4e12 buckets the other way. THIS HEADER THEREFORE INCLUDES NO
// core/snap_grid.hpp, and that is grep-checkable.
//
// IT IS ALSO KERNEL-FREE: it compares bounding boxes and never asks an
// orientation, which is why prop_noding_broad_phase registers with the plain
// add_terrain_test helper and a broad phase that reached for a predicate would
// fail to LINK rather than pass quietly.
//
// The index is immutable after construction and const-queried. That is what lets
// the split pass and the verification pass share one instance, and what lets a
// future parallel-for over segments query it without synchronisation.

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

namespace terrain::noding {

class BroadPhase {
public:
    // Copies each segment's bounding box; the span need not outlive the index.
    // Throws std::invalid_argument on a non-finite endpoint, via bbox().
    explicit BroadPhase(std::span<const Segment2> segments) {
        boxes_.reserve(segments.size());
        for (const Segment2& s : segments) {
            boxes_.push_back(bbox(s));
            extent_.expand(boxes_.back());
        }

        const auto n = static_cast<double>(segments.size());
        const auto k = static_cast<std::size_t>(std::max(1.0, std::ceil(std::sqrt(n))));
        // A zero-width or zero-height extent gives one bucket on that axis, which
        // is what keeps the division below out of the query's hot path without a
        // branch of its own.
        nx_ = extent_.width() > 0.0 ? k : 1;
        ny_ = extent_.height() > 0.0 ? k : 1;

        buckets_.resize(nx_ * ny_);
        for (std::size_t i = 0; i < boxes_.size(); ++i) {
            const Cells c = cells_of(boxes_[i]);
            for (std::size_t y = c.y0; y <= c.y1; ++y) {
                for (std::size_t x = c.x0; x <= c.x1; ++x) {
                    buckets_[y * nx_ + x].push_back(static_cast<std::uint32_t>(i));
                }
            }
        }
    }

    [[nodiscard]] std::size_t size() const noexcept { return boxes_.size(); }

    // F is invoked as f(std::uint32_t segment_index), at most once per segment.
    //
    // The at-most-once is not part of the contract a caller may lean on for
    // correctness -- "may visit others" already permits repeats -- but it is
    // cheap here and it keeps the caller's pair enumeration from growing with
    // the number of buckets a long segment spans. It falls out of emitting a
    // segment only from the FIRST bucket its range shares with the query's,
    // which needs no per-query state and so leaves the object const and
    // thread-safe.
    template <class F>
    void for_each_candidate(Box2 q, F&& f) const {
        if (q.is_empty() || buckets_.empty()) {
            return;
        }
        const Cells qc = cells_of(q);
        const std::size_t qx0 = qc.x0;
        const std::size_t qy0 = qc.y0;
        for (std::size_t y = qy0; y <= qc.y1; ++y) {
            for (std::size_t x = qx0; x <= qc.x1; ++x) {
                for (const std::uint32_t i : buckets_[y * nx_ + x]) {
                    const Cells s = cells_of(boxes_[i]);
                    if (x == std::max(s.x0, qx0) && y == std::max(s.y0, qy0)) {
                        f(i);
                    }
                }
            }
        }
    }

private:
    struct Cells {
        std::size_t x0, y0, x1, y1;
    };

    [[nodiscard]] std::size_t axis(double v, double lo, double span, std::size_t n) const noexcept {
        if (n == 1) {
            return 0;
        }
        const double t = std::floor((v - lo) / span * static_cast<double>(n));
        if (!(t > 0.0)) {
            return 0;
        }
        const auto limit = static_cast<double>(n - 1);
        return t >= limit ? n - 1 : static_cast<std::size_t>(t);
    }

    // Closed on both ends, so a box exactly on a bucket boundary occupies both
    // buckets; that is the overlap rule the postcondition's one-line argument
    // rests on.
    [[nodiscard]] Cells cells_of(const Box2& b) const noexcept {
        const double lx = extent_.lo().x;
        const double ly = extent_.lo().y;
        const double wx = extent_.width();
        const double wy = extent_.height();
        return Cells{axis(b.lo().x, lx, wx, nx_), axis(b.lo().y, ly, wy, ny_),
                     axis(b.hi().x, lx, wx, nx_), axis(b.hi().y, ly, wy, ny_)};
    }

    std::vector<Box2> boxes_;
    std::vector<std::vector<std::uint32_t>> buckets_;
    Box2 extent_;
    std::size_t nx_{1};
    std::size_t ny_{1};
};

}  // namespace terrain::noding
