#pragma once

// Land-cover borders simplified within a band, every face's area kept
// (docs/increments/32-landcover-simplify.md, sections 6 to 8). The coverage is
// cut at its junctions into borders (de Berg, van Kreveld and Schirra 1998);
// each border is simplified by area-preserving segment collapse (Kronenfeld et
// al. 2020, as reduce_ring), so both neighbours of a border keep their area. A
// collapse is made only if (i) the new edges cross or touch no edge of any
// border except at their shared ends, (ii) no vertex of any border lies in the
// swept region or on a new edge, and the anchored check holds: E has an anchor
// F_E on the source border, F_A <= F_E <= F_D in order along it, |E - F_E| <=
// band, and every source vertex between the anchors is within the band of its
// new edge. That bounds the Hausdorff distance both ways (section 7). With a
// clearance (section 15.2), E is placed at least that far from A and D, and a
// collapse is refused if a vertex comes closer than it to a new edge it does
// not end, or E to an edge other than A-B, B-C, C-D.
//
// Locality: every collapse reads four nodes, a stretch of its source border and
// a grid query; the heap is global and the loop serial, ordered by (deviation,
// node), so the output depends on nothing but the input.

#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/noding/intersect.hpp>
#include <terrain/predicates/kernel.hpp>
#include <terrain/predicates/orientation.hpp>
#include <terrain/vector_simplify/area_collapse.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>
#include <queue>
#include <span>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace terrain::vector_simplify {

enum class BorderStatus : std::uint8_t { Ok, InvalidBand, BadRings, InvalidClearance };

struct BorderCounts {
    std::size_t junctions{}, borders{}, fixed_borders{};
    std::size_t collinear{}, collapses{};
    std::size_t rejected_crossing{}, rejected_side{};
    std::size_t rejected_clearance{}, skipped_placements{}; // section 15.2
};

struct BorderOutcome {
    std::vector<Point2> points;             // every ring, open, back to back
    std::vector<std::uint64_t> ring_starts; // ring k is [starts[k], starts[k+1])
    BorderStatus status{BorderStatus::Ok};
    BorderCounts counts{};
};

namespace detail {

inline constexpr std::size_t no_node = std::numeric_limits<std::size_t>::max();

// A position on a border's source polyline: segment i, parameter t in [0, 1].
struct Anchor {
    std::size_t i;
    double t;
};

// One border: its source polyline (closed: first vertex not repeated), its end
// nodes (open) or a live node (closed), and its live node count.
struct Border {
    std::vector<Point2> source;
    bool closed, fixed;
    std::size_t first, last, size;
};

struct Placed {
    Point2 e;
    Anchor f;
    double deviation;
};

// The anchored check for E between A and D (anchors fa, fd): over every
// source segment from fa's to fd's, E's nearest point there in order between
// fa and fd, the deviation is the largest of |E - F_E| and the distances of
// the source vertices between fa and F_E to A-E and between F_E and fd to
// E-D; the smallest such deviation.
[[nodiscard]] inline std::optional<Placed> anchored(const Border& g, Anchor fa, Anchor fd,
                                                    const Point2& a, const Point2& e,
                                                    const Point2& d) {
    const std::size_t m = g.source.size(), segments = g.closed ? m : m - 1;
    std::size_t span = g.closed ? (fd.i + m - fa.i) % m : fd.i - fa.i;
    if (g.closed && span == 0 && fd.t < fa.t)
        span = m;
    const auto at = [&](std::size_t r) -> const Point2& { return g.source[(fa.i + r) % m]; };
    std::vector<double> suffix(span + 2, 0.0);
    for (std::size_t r = span; r >= 1; --r)
        suffix[r] = std::max(suffix[r + 1], segment_distance(at(r), e, d));
    std::optional<Placed> best;
    double prefix = 0.0;
    for (std::size_t r = 0; r <= span; ++r) {
        if (r > 0)
            prefix = std::max(prefix, segment_distance(at(r), a, e));
        const std::size_t j = g.closed ? (fa.i + r) % m : fa.i + r;
        if (j >= segments)
            continue;
        const Point2 p = at(r), q = at(r + 1) - p;
        const double len2 = dot(q, q);
        const double t = len2 > 0.0 ? std::clamp(dot(e - p, q) / len2, 0.0, 1.0) : 0.0;
        if ((r == 0 && t < fa.t) || (r == span && t > fd.t))
            continue; // the anchor order, F_A <= F_E <= F_D
        const Point2 f = p + t * q;
        const double v = std::max({prefix, suffix[r + 1], std::hypot(e.x - f.x, e.y - f.y)});
        if (!best || v < best->deviation)
            best = Placed{e, Anchor{j, t}, v};
    }
    return best;
}

template <pred::GeometryKernel K>
[[nodiscard]] int winding(std::span<const Point2> loop, const Point2& k) {
    using pred::Orientation;
    int w = 0;
    for (std::size_t i = 0; i < loop.size(); ++i) {
        const Point2& p = loop[i];
        const Point2& q = loop[(i + 1) % loop.size()];
        const Orientation o = K::orient2d(p, q, k);
        w += p.y <= k.y && q.y > k.y && o == Orientation::CounterClockwise;
        w -= p.y > k.y && q.y <= k.y && o == Orientation::Clockwise;
    }
    return w;
}

} // namespace detail

template <pred::GeometryKernel K>
[[nodiscard]] BorderOutcome simplify_borders(std::span<const Point2> points,
                                             std::span<const std::uint64_t> ring_starts,
                                             double band, double clearance = 0.0) {
    using noding::SegmentRelation;
    using detail::no_node;
    BorderOutcome out;
    if (!std::isfinite(band) || band < 0.0) {
        out.status = BorderStatus::InvalidBand;
        return out;
    }
    if (!std::isfinite(clearance) || clearance < 0.0) {
        out.status = BorderStatus::InvalidClearance;
        return out;
    }
    bool bad = ring_starts.empty() || ring_starts.front() != 0 || ring_starts.back() != points.size();
    for (std::size_t k = 0; !bad && k + 1 < ring_starts.size(); ++k)
        bad = ring_starts[k + 1] < ring_starts[k] || ring_starts[k + 1] - ring_starts[k] < 3;
    for (const Point2& p : points)
        bad = bad || !std::isfinite(p.x) || !std::isfinite(p.y);
    if (bad) {
        out.status = BorderStatus::BadRings;
        return out;
    }
    out.points.assign(points.begin(), points.end());
    out.ring_starts.assign(ring_starts.begin(), ring_starts.end());
    if (band == 0.0)
        return out;

    // a. Vertices by exact coordinates (-0.0 is 0.0), each at its first
    // occurrence's coordinates; edges with the rings that use them.
    const std::size_t n = points.size(), rings = ring_starts.size() - 1;
    const auto xy = [&](std::size_t i) { return std::pair{points[i].x + 0.0, points[i].y + 0.0}; };
    std::vector<std::size_t> order(n), vid(n), ring_of(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::stable_sort(order, {}, xy);
    std::vector<Point2> vertex;
    for (std::size_t k = 0; k < n; ++k) {
        if (k == 0 || xy(order[k]) != xy(order[k - 1]))
            vertex.push_back(points[order[k]]);
        vid[order[k]] = vertex.size() - 1;
    }
    const std::uint64_t nv = vertex.size();
    for (std::size_t k = 0; k < rings; ++k)
        std::fill(ring_of.begin() + static_cast<std::ptrdiff_t>(ring_starts[k]),
                  ring_of.begin() + static_cast<std::ptrdiff_t>(ring_starts[k + 1]), k);
    const auto succ = [&](std::size_t i) { return i + 1 == ring_starts[ring_of[i] + 1] ? ring_starts[ring_of[i]] : i + 1; };
    const auto edge_key = [&](std::size_t u, std::size_t v) { return std::min(u, v) * nv + std::max(u, v); };
    std::vector<std::pair<std::uint64_t, std::size_t>> uses(n);
    for (std::size_t i = 0; i < n; ++i)
        uses[i] = {edge_key(vid[i], vid[succ(i)]), ring_of[i]};
    std::ranges::sort(uses);
    std::vector<std::uint64_t> keys;
    std::vector<std::vector<std::size_t>> users;
    for (const auto& [e, r] : uses) {
        if (keys.empty() || keys.back() != e) {
            keys.push_back(e);
            users.emplace_back();
        }
        users.back().push_back(r);
    }
    const auto edge_of = [&](std::size_t u, std::size_t v) {
        return static_cast<std::size_t>(std::ranges::lower_bound(keys, edge_key(u, v)) - keys.begin());
    };
    // b. Junctions: other than two distinct incident edges, or two whose rings
    // differ. A shared edge is used by exactly two distinct rings.
    std::vector<std::size_t> degree(nv, 0), incident(2 * nv, 0);
    for (std::size_t e = 0; e < keys.size(); ++e)
        for (const std::uint64_t v : {keys[e] / nv, keys[e] % nv})
            if (degree[v]++ < 2)
                incident[2 * v + degree[v] - 1] = e;
    std::vector<bool> junction(nv);
    for (std::size_t v = 0; v < nv; ++v) {
        junction[v] = degree[v] != 2 || users[incident[2 * v]] != users[incident[2 * v + 1]]
                      || keys[incident[2 * v]] / nv == keys[incident[2 * v]] % nv;
        out.counts.junctions += junction[v];
    }
    const auto shared = [&](std::size_t u, std::size_t v) {
        const auto& who = users[edge_of(u, v)];
        return who.size() == 2 && who[0] != who[1];
    };

    // Borders and their nodes; each ring as a list of (border, forward).
    std::vector<detail::Border> borders;
    std::vector<Point2> p;
    std::vector<std::size_t> prev, next, owner;
    std::vector<detail::Anchor> anchor;
    std::vector<std::uint32_t> version;
    std::vector<bool> alive;
    const auto add_node = [&](Point2 q, std::size_t b, detail::Anchor f) {
        p.push_back(q);
        prev.push_back(no_node);
        next.push_back(no_node);
        owner.push_back(b);
        anchor.push_back(f);
        version.push_back(0);
        alive.push_back(true);
        return p.size() - 1;
    };
    std::unordered_map<std::uint64_t, std::pair<std::size_t, bool>> by_edge; // directed first edge
    const auto directed = [&](std::size_t u, std::size_t v) { return std::uint64_t{u} * nv + v; };
    const auto border_of = [&](const std::vector<std::size_t>& ids, bool closed) {
        const std::size_t m = ids.size();
        if (const auto it = by_edge.find(directed(ids[0], ids[1])); it != by_edge.end())
            return it->second;
        const std::size_t b = borders.size();
        detail::Border g{{}, closed, !shared(ids[0], ids[1]), p.size(), 0, m};
        for (std::size_t k = 0; k < m; ++k) {
            g.source.push_back(vertex[ids[k]]);
            const std::size_t u = add_node(vertex[ids[k]], b, {k, 0.0});
            if (k > 0) {
                next[u - 1] = u;
                prev[u] = u - 1;
            }
        }
        g.last = p.size() - 1;
        if (closed) {
            next[g.last] = g.first;
            prev[g.first] = g.last;
        }
        out.counts.fixed_borders += g.fixed;
        borders.push_back(std::move(g));
        by_edge[directed(ids[0], ids[1])] = {b, true};
        by_edge[directed(ids[closed ? 0 : m - 1], ids[closed ? m - 1 : m - 2])] = {b, false};
        return std::pair{b, true};
    };
    std::vector<std::vector<std::pair<std::size_t, bool>>> ring_borders(rings);
    for (std::size_t k = 0; k < rings; ++k) {
        const std::size_t s = ring_starts[k], len = ring_starts[k + 1] - s;
        const auto id = [&](std::size_t r) { return vid[s + r % len]; };
        std::vector<std::size_t> at;
        for (std::size_t r = 0; r < len; ++r)
            if (junction[id(r)])
                at.push_back(r);
        if (at.empty()) { // one closed border, started at its smallest vertex
            std::size_t low = 0;
            for (std::size_t r = 1; r < len; ++r)
                low = id(r) < id(low) ? r : low;
            std::vector<std::size_t> ids;
            for (std::size_t r = 0; r < len; ++r)
                ids.push_back(id(low + r));
            ring_borders[k].push_back(border_of(ids, true));
            continue;
        }
        for (std::size_t j = 0; j < at.size(); ++j) {
            const std::size_t to = j + 1 < at.size() ? at[j + 1] : at[0] + len;
            std::vector<std::size_t> ids;
            for (std::size_t r = at[j]; r <= to; ++r)
                ids.push_back(id(r));
            ring_borders[k].push_back(border_of(ids, false));
        }
    }
    out.counts.borders = borders.size();

    // The collinear pass: drop every inner vertex exactly between its
    // neighbours, on borders that may move.
    for (std::size_t u = 0; u < p.size(); ++u) {
        detail::Border& g = borders[owner[u]];
        if (g.fixed || prev[u] == no_node || next[u] == no_node || (g.closed && g.size <= 3))
            continue;
        const Point2 &a = p[prev[u]], &b = p[next[u]];
        const bool between = std::min(a.x, b.x) <= p[u].x && p[u].x <= std::max(a.x, b.x)
                             && std::min(a.y, b.y) <= p[u].y && p[u].y <= std::max(a.y, b.y);
        if (K::orient2d(a, p[u], b) != pred::Orientation::Collinear || !between)
            continue;
        alive[u] = false;
        next[prev[u]] = next[u];
        prev[next[u]] = prev[u];
        g.first = g.first == u ? next[u] : g.first;
        --g.size;
        ++out.counts.collinear;
    }

    // The grid of current edges: bucket side the band, at least the mean edge.
    double total = 0.0;
    std::size_t edges = 0;
    Point2 low = p.empty() ? Point2{} : p[0];
    for (std::size_t u = 0; u < p.size(); ++u) {
        low = Point2{std::min(low.x, p[u].x), std::min(low.y, p[u].y)};
        if (alive[u] && next[u] != no_node) {
            total += std::hypot(p[next[u]].x - p[u].x, p[next[u]].y - p[u].y);
            ++edges;
        }
    }
    detail::EdgeGrid grid{std::max(band, total / static_cast<double>(std::max<std::size_t>(edges, 1))), low};
    for (std::size_t u = 0; u < p.size(); ++u)
        if (alive[u] && next[u] != no_node)
            grid.insert(u, p[u], p[next[u]]);

    // c. Candidates: E on line A-B or on line C-D, or the foot of B-C's
    // midpoint when both are parallel to A-D; the smaller anchored deviation.
    const auto candidate = [&](std::size_t b) -> std::optional<detail::Placed> {
        const detail::Border& g = borders[owner[b]];
        const std::size_t a = prev[b], c = next[b];
        if (g.fixed || (g.closed && g.size <= 4) || a == no_node || c == no_node || next[c] == no_node)
            return std::nullopt;
        const std::size_t d = next[c];
        const Point2 A = p[a], vb = p[b] - A, vc = p[c] - A, vd = p[d] - A;
        const double area = cross(vb, vc) + cross(vc, vd); // A-E-D must enclose this
        std::optional<detail::Placed> best;
        const auto consider = [&](Point2 e) {
            if (!std::isfinite(e.x) || !std::isfinite(e.y))
                return;
            if (std::hypot(e.x, e.y) < clearance || std::hypot(e.x - vd.x, e.y - vd.y) < clearance) {
                ++out.counts.skipped_placements; // too close to A or D
                return;
            }
            const auto placed = detail::anchored(g, anchor[a], anchor[d], A, A + e, p[d]);
            if (placed && (!best || placed->deviation < best->deviation))
                best = placed;
        };
        if (const double bd = cross(vb, vd); bd != 0.0)
            consider(vb * (area / bd));
        if (const double cd = cross(vc, vd); cd != 0.0)
            consider(vc + (vd - vc) * (1.0 - area / cd));
        if (!best && dot(vd, vd) > 0.0 && cross(vb, vd) == 0.0 && cross(vc, vd) == 0.0) {
            const Point2 mid = (vb + vc) * 0.5;
            consider(mid + Point2{-vd.y, vd.x} * ((cross(mid, vd) - area) / dot(vd, vd)));
        }
        return best;
    };

    using Entry = std::tuple<double, std::size_t, std::uint32_t>; // deviation, node, version
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> heap;
    std::vector<std::optional<detail::Placed>> pending(p.size());
    const auto evaluate = [&](std::size_t b) {
        if (b == no_node)
            return;
        ++version[b];
        pending[b] = candidate(b);
        if (pending[b] && pending[b]->deviation <= band)
            heap.emplace(pending[b]->deviation, b, version[b]);
    };
    for (std::size_t u = 0; u < p.size(); ++u)
        if (alive[u])
            evaluate(u);

    std::vector<std::size_t> stamp(p.size(), 0);
    std::size_t round = 0;
    while (!heap.empty()) {
        const auto [dev, b, ver] = heap.top();
        heap.pop();
        detail::Border& g = borders[owner[b]];
        if (!alive[b] || version[b] != ver || (g.closed && g.size <= 4))
            continue;
        const std::size_t a = prev[b], c = next[b], d = next[c];
        const detail::Placed cand = *pending[b];
        const Point2 A = p[a], E = cand.e, D = p[d];
        const std::array<Point2, 5> loop{A, p[b], p[c], D, E};
        Point2 lo = A, hi = A;
        for (const Point2& q : loop) {
            lo = Point2{std::min(lo.x, q.x - clearance), std::min(lo.y, q.y - clearance)};
            hi = Point2{std::max(hi.x, q.x + clearance), std::max(hi.y, q.y + clearance)};
        }
        // (i) A-E and E-D cross or touch no current edge but at a shared end;
        // (ii) no vertex lies in the swept region or on a new edge; the
        // clearance (near), in floating point, does not end the scan. The scan
        // skips the edge from A, so A against E-D and D against A-E here.
        const Segment2 ae{A, E}, ed{E, D};
        bool crosses = noding::classify<K>(ae, ed) != SegmentRelation::Touching, side = false;
        bool near = detail::segment_distance(A, E, D) < clearance || detail::segment_distance(D, A, E) < clearance;
        ++round;
        grid.query(lo, hi, hi, [&](std::size_t u) {
            if (crosses || side || !alive[u] || next[u] == no_node || stamp[u] == round)
                return;
            stamp[u] = round;
            if (u == a || u == b || u == c)
                return;
            const std::size_t v = next[u];
            const Segment2 edge{p[u], p[v]};
            const auto want = [&](const Point2& end) {
                return p[u] == end || p[v] == end ? SegmentRelation::Touching : SegmentRelation::Disjoint;
            };
            crosses = noding::classify<K>(ae, edge) != want(A) || noding::classify<K>(ed, edge) != want(D);
            near = near || detail::segment_distance(E, p[u], p[v]) < clearance;
            for (const std::size_t w : {u, v}) {
                if (crosses || w == b || w == c)
                    continue;
                // A's copies only from A-E, D's only from E-D (15.2).
                near = near || (p[w] != A && detail::segment_distance(p[w], A, E) < clearance)
                       || (p[w] != D && detail::segment_distance(p[w], E, D) < clearance);
                if (p[w] == A || p[w] == D)
                    continue;
                const Segment2 at{p[w], p[w]};
                side = side || detail::winding<K>(loop, p[w]) != 0
                       || noding::classify<K>(at, ae) != SegmentRelation::Disjoint
                       || noding::classify<K>(at, ed) != SegmentRelation::Disjoint;
            }
        });
        if (crosses || side || near) {
            ++(crosses ? out.counts.rejected_crossing : side ? out.counts.rejected_side : out.counts.rejected_clearance);
            continue;
        }
        // Apply, then re-evaluate the candidates whose four nodes changed.
        const std::size_t e = add_node(E, owner[b], cand.f);
        pending.emplace_back();
        stamp.push_back(0);
        alive[b] = alive[c] = false;
        next[a] = e;
        prev[e] = a;
        next[e] = d;
        prev[d] = e;
        g.first = g.first == b || g.first == c ? e : g.first;
        grid.insert(a, A, E);
        grid.insert(e, E, D);
        --g.size;
        ++out.counts.collapses;
        for (const std::size_t u : {prev[a], a, e, d})
            evaluate(u);
    }

    // d. Rings rebuilt from their borders: same rings, order and orientation.
    out.points.clear();
    out.ring_starts.assign(1, 0);
    for (const auto& list : ring_borders) {
        for (const auto& [id, forward] : list) {
            const detail::Border& g = borders[id];
            const std::size_t start = forward || g.closed ? g.first : g.last;
            const std::size_t stop = g.closed ? g.first : forward ? g.last : g.first;
            std::size_t u = start;
            do {
                out.points.push_back(p[u]);
                u = forward ? next[u] : prev[u];
            } while (u != stop);
        }
        out.ring_starts.push_back(out.points.size());
    }
    return out;
}

} // namespace terrain::vector_simplify
