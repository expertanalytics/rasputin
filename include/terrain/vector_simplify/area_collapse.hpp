#pragma once

// Area-preserving segment collapse of a simple ring to a horizontal tolerance
// (docs/increments/22-auto-catchment.md, "The reduction"; after Kronenfeld,
// Stanislawski, Buttenfield and Brockmeyer 2020, APSC). Each collapse replaces
// the two vertices B, C of a chain A-B-C-D by one point E on the line where
// A-E-D encloses the area A-B-C-D did, so the ring's area is kept up to
// rounding. A collapse is made only if every fine vertex it stands for stays
// within the tolerance of the new edges, E within it of the fine ring, the
// ring stays simple (exact kernel) and every keep-point stays inside.
//
// Locality: every collapse reads four vertices, their fine ranges and a grid
// query; the heap is global and the loop serial, ordered by (deviation, id),
// so the output depends on nothing but the input.

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <queue>
#include <span>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace terrain::vector_simplify {

enum class ReduceStatus : std::uint8_t { Ok, InvalidTolerance, NotCounterClockwise, TooFewVertices };

struct ReduceCounts {
    std::size_t collinear{};  // vertices dropped by the collinear pass
    std::size_t collapses{};  // B, C -> E replacements made
    std::size_t rejected_crossing{}, rejected_seed{}, rejected_tolerance{};
};

struct ReduceOutcome {
    std::vector<Point2> ring; // open: the first vertex not repeated
    ReduceStatus status{ReduceStatus::Ok};
    ReduceCounts counts{};
};

namespace detail {

[[nodiscard]] inline double segment_distance(const Point2& p, const Point2& a, const Point2& b) {
    const Point2 d = b - a;
    const double len2 = dot(d, d);
    const double t = len2 > 0.0 ? std::clamp(dot(p - a, d) / len2, 0.0, 1.0) : 0.0;
    const Point2 q = a + t * d;
    return std::hypot(p.x - q.x, p.y - q.y);
}

// The current ring: a doubly linked list over nodes. Edge u -> next[u] stands
// for `count[u]` fine vertices starting at fine index `start[u]`.
struct Ring {
    std::vector<Point2> p;
    std::vector<std::size_t> prev, next, start, count, id;
    std::vector<std::uint32_t> version;
    std::vector<bool> alive;

    std::size_t add(Point2 q, std::size_t ident) {
        p.push_back(q);
        prev.push_back(0);
        next.push_back(0);
        start.push_back(0);
        count.push_back(0);
        id.push_back(ident);
        version.push_back(0);
        alive.push_back(true);
        return p.size() - 1;
    }
};

// Current edges by grid bucket, keyed by their start node. Entries are never
// removed: a stale one is filtered by the caller against the live ring.
class EdgeGrid {
public:
    EdgeGrid(double side, Point2 origin) : side_{side}, origin_{origin} {}

    void insert(std::size_t u, const Point2& a, const Point2& b) {
        each_cell(a, b, a, [&](std::uint64_t key) { cells_[key].push_back(u); });
    }

    template <typename F>
    void query(const Point2& a, const Point2& b, const Point2& c, F&& f) const {
        each_cell(a, b, c, [&](std::uint64_t key) {
            if (const auto it = cells_.find(key); it != cells_.end())
                for (const std::size_t u : it->second)
                    f(u);
        });
    }

private:
    template <typename F>
    void each_cell(const Point2& a, const Point2& b, const Point2& c, F&& f) const {
        const auto cell = [&](double v, double o) {
            return static_cast<std::int64_t>(std::floor((v - o) / side_));
        };
        const auto x0 = cell(std::min({a.x, b.x, c.x}), origin_.x);
        const auto x1 = cell(std::max({a.x, b.x, c.x}), origin_.x);
        const auto y0 = cell(std::min({a.y, b.y, c.y}), origin_.y);
        const auto y1 = cell(std::max({a.y, b.y, c.y}), origin_.y);
        for (auto ix = x0; ix <= x1; ++ix)
            for (auto iy = y0; iy <= y1; ++iy)
                f((static_cast<std::uint64_t>(ix) << 32) ^ static_cast<std::uint64_t>(iy & 0xffffffff));
    }

    double side_;
    Point2 origin_;
    std::unordered_map<std::uint64_t, std::vector<std::size_t>> cells_;
};

struct Candidate {
    Point2 e;
    std::size_t split; // fine vertices [0, split) of the range go to A-E, the rest to E-D
    double deviation;
};

// The deviation of placing E for the chain whose range is `fine` from `s`,
// `n` vertices long: the best split of the range between A-E and E-D, and E's
// distance to the fine segments in the range.
[[nodiscard]] inline Candidate deviation_of(std::span<const Point2> fine, std::size_t s,
                                            std::size_t n, const Point2& a, const Point2& e,
                                            const Point2& d) {
    const std::size_t size = fine.size();
    const auto at = [&](std::size_t k) -> const Point2& { return fine[(s + k) % size]; };
    std::vector<double> suffix(n + 1, 0.0);
    for (std::size_t k = n; k-- > 0;)
        suffix[k] = std::max(suffix[k + 1], segment_distance(at(k), e, d));
    Candidate best{e, 0, suffix[0]};
    double prefix = 0.0;
    for (std::size_t m = 1; m <= n; ++m) {
        prefix = std::max(prefix, segment_distance(at(m - 1), a, e));
        if (const double v = std::max(prefix, suffix[m]); v < best.deviation)
            best = Candidate{e, m, v};
    }
    double to_fine = std::numeric_limits<double>::infinity();
    for (std::size_t k = 0; k < std::max<std::size_t>(n, 1); ++k)
        to_fine = std::min(to_fine, segment_distance(e, at(k), at(k + 1)));
    best.deviation = std::max(best.deviation, to_fine);
    return best;
}

} // namespace detail

template <pred::GeometryKernel K>
[[nodiscard]] ReduceOutcome reduce_ring(std::span<const Point2> fine, double tolerance,
                                        std::span<const Point2> keep) {
    using pred::Orientation;
    ReduceOutcome out;
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        out.status = ReduceStatus::InvalidTolerance;
        return out;
    }
    const std::size_t n = fine.size();
    if (n < 4) {
        out.status = ReduceStatus::TooFewVertices;
        return out;
    }
    double twice = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        twice += cross(fine[i] - fine[0], fine[(i + 1) % n] - fine[0]);
    if (!(twice > 0.0)) {
        out.status = ReduceStatus::NotCounterClockwise;
        return out;
    }

    // 1. The collinear pass: drop every vertex exactly between its neighbours.
    std::vector<std::size_t> kept;
    for (std::size_t i = 0; i < n; ++i) {
        const Point2& a = fine[(i + n - 1) % n];
        const Point2& b = fine[(i + 1) % n];
        const bool between = std::min(a.x, b.x) <= fine[i].x && fine[i].x <= std::max(a.x, b.x)
                             && std::min(a.y, b.y) <= fine[i].y && fine[i].y <= std::max(a.y, b.y);
        if (K::orient2d(a, fine[i], b) != Orientation::Collinear || !between)
            kept.push_back(i);
    }
    out.counts.collinear = n - kept.size();
    detail::Ring ring;
    for (std::size_t k = 0; k < kept.size(); ++k) {
        const std::size_t u = ring.add(fine[kept[k]], kept[k]);
        ring.start[u] = kept[k];
        ring.count[u] = (kept[(k + 1) % kept.size()] + n - kept[k]) % n;
    }
    const std::size_t m0 = kept.size();
    for (std::size_t u = 0; u < m0; ++u) {
        ring.next[u] = (u + 1) % m0;
        ring.prev[u] = (u + m0 - 1) % m0;
    }
    std::size_t size = m0, next_id = n;
    const auto emit = [&] {
        std::size_t u = 0;
        while (!ring.alive[u])
            ++u;
        for (std::size_t k = 0; k < size; ++k, u = ring.next[u])
            out.ring.push_back(ring.p[u]);
        return out;
    };
    if (tolerance == 0.0 || size <= 4)
        return emit();

    // The grid of current edges: bucket side the tolerance, at least the mean
    // fine edge, so a query meets few buckets.
    double perimeter = 0.0;
    Point2 low = fine[0];
    for (std::size_t i = 0; i < n; ++i) {
        const Point2 d = fine[(i + 1) % n] - fine[i];
        perimeter += std::hypot(d.x, d.y);
        low = Point2{std::min(low.x, fine[i].x), std::min(low.y, fine[i].y)};
    }
    detail::EdgeGrid grid{std::max(tolerance, perimeter / static_cast<double>(n)), low};
    for (std::size_t u = 0; u < m0; ++u)
        grid.insert(u, ring.p[u], ring.p[ring.next[u]]);

    // 2.-4. Candidates: E on line A-B or on line C-D (ties to A-B), or the foot
    // of B-C's midpoint when both are parallel to A-D.
    const auto candidate = [&](std::size_t b) -> std::optional<detail::Candidate> {
        const std::size_t a = ring.prev[b], c = ring.next[b], d = ring.next[c];
        const Point2 A = ring.p[a], vb = ring.p[b] - A, vc = ring.p[c] - A, vd = ring.p[d] - A;
        const double area = cross(vb, vc) + cross(vc, vd); // A-E-D must enclose this
        const std::size_t s = ring.start[a], len = ring.count[a] + ring.count[b] + ring.count[c];
        std::optional<detail::Candidate> best;
        const auto consider = [&](Point2 e) {
            if (!std::isfinite(e.x) || !std::isfinite(e.y))
                return;
            const detail::Candidate cand = detail::deviation_of(fine, s, len, A, A + e, ring.p[d]);
            if (!best || cand.deviation < best->deviation)
                best = cand;
        };
        if (const double bd = cross(vb, vd); bd != 0.0)
            consider(vb * (area / bd));
        if (const double cd = cross(vc, vd); cd != 0.0)
            consider(vc + (vd - vc) * (1.0 - area / cd));
        if (!best && dot(vd, vd) > 0.0) {
            const Point2 mid = (vb + vc) * 0.5;
            consider(mid + Point2{-vd.y, vd.x} * ((cross(mid, vd) - area) / dot(vd, vd)));
        }
        return best;
    };

    using Entry = std::tuple<double, std::size_t, std::size_t, std::uint32_t>; // dev, id, node, version
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> heap;
    std::vector<std::optional<detail::Candidate>> pending(m0);
    const auto evaluate = [&](std::size_t b) {
        ++ring.version[b];
        pending[b] = candidate(b);
        if (pending[b] && pending[b]->deviation <= tolerance)
            heap.emplace(pending[b]->deviation, ring.id[b], b, ring.version[b]);
        else if (pending[b])
            ++out.counts.rejected_tolerance;
    };
    for (std::size_t u = 0; u < m0; ++u)
        evaluate(u);

    std::vector<std::size_t> stamp(m0, 0);
    std::size_t round = 0;
    while (!heap.empty() && size > 4) {
        const auto [dev, ident, b, version] = heap.top();
        heap.pop();
        if (!ring.alive[b] || ring.version[b] != version)
            continue;
        const std::size_t a = ring.prev[b], c = ring.next[b], d = ring.next[c];
        const detail::Candidate cand = *pending[b];
        const Point2 A = ring.p[a], E = cand.e, D = ring.p[d];

        // 5. Simple: A-E and E-D against every current edge near them.
        const Segment2 ae{A, E}, ed{E, D};
        bool crosses = noding::classify<K>(ae, ed) != noding::SegmentRelation::Touching;
        ++round;
        grid.query(A, E, D, [&](std::size_t u) {
            if (crosses || !ring.alive[u] || stamp[u] == round || u == a || u == b || u == c)
                return;
            stamp[u] = round;
            const std::size_t v = ring.next[u];
            const Segment2 edge{ring.p[u], ring.p[v]};
            const auto want = [&](bool shared) {
                return shared ? noding::SegmentRelation::Touching : noding::SegmentRelation::Disjoint;
            };
            crosses = noding::classify<K>(ae, edge) != want(v == a)
                      || noding::classify<K>(ed, edge) != want(u == d);
        });
        if (crosses) {
            ++out.counts.rejected_crossing;
            continue;
        }
        // Keep-points: winding zero round A-B-C-D-E-A, and on neither new edge.
        const std::array<Point2, 5> loop{A, ring.p[b], ring.p[c], D, E};
        bool lost = false;
        for (const Point2& k : keep) {
            int winding = 0;
            for (std::size_t i = 0; i < 5; ++i) {
                const Point2& p = loop[i];
                const Point2& q = loop[(i + 1) % 5];
                const Orientation o = K::orient2d(p, q, k);
                winding += p.y <= k.y && q.y > k.y && o == Orientation::CounterClockwise;
                winding -= p.y > k.y && q.y <= k.y && o == Orientation::Clockwise;
            }
            const Segment2 at{k, k};
            lost = lost || winding != 0
                   || noding::classify<K>(at, ae) != noding::SegmentRelation::Disjoint
                   || noding::classify<K>(at, ed) != noding::SegmentRelation::Disjoint;
        }
        if (lost) {
            ++out.counts.rejected_seed;
            continue;
        }

        // 6. Apply, then re-evaluate the candidates whose four vertices changed.
        const std::size_t len = ring.count[a] + ring.count[b] + ring.count[c];
        const std::size_t e = ring.add(E, next_id++);
        pending.emplace_back();
        stamp.push_back(0);
        ring.alive[b] = ring.alive[c] = false;
        ring.next[a] = e;
        ring.prev[e] = a;
        ring.next[e] = d;
        ring.prev[d] = e;
        ring.count[a] = cand.split;
        ring.start[e] = (ring.start[a] + cand.split) % n;
        ring.count[e] = len - cand.split;
        grid.insert(a, A, E);
        grid.insert(e, E, D);
        --size;
        ++out.counts.collapses;
        if (size > 4)
            for (const std::size_t u : {ring.prev[a], a, e, d})
                evaluate(u);
    }
    return emit();
}

} // namespace terrain::vector_simplify
