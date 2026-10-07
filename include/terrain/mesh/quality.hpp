#pragma once

// A minimum-angle pass over a LatticeMesh, run once on the start mesh before
// DEM refinement (docs/increments/20-start-quality.md, R2 to R6, R10).
//
// A triangle is bad when its circumradius-to-shortest-edge ratio, in the world
// frame (col * dx, -(row * dy)), exceeds 1 / (2 sin theta): its minimum angle
// is below theta. The ratio is a double proposal only; nothing topological
// depends on it. For a bad triangle the pass inserts the DEM node nearest its
// circumcentre (R3), located by a visibility walk with the exact orient_sign,
// and inserted with split_inside or split_edge -- the latter exactly when the
// node lies on an edge, constrained or not (R6) -- then legalised around.
//
// Skips, in the order the code tests them (R4's list, with the vertex test
// after the walk, which is what finds it): circumradius below the floor
// sqrt(dx^2 + dy^2), so the node is within R / 2 of the centre (R5); centre
// outside the node rectangle (never clamped); the node rejected by the caller's
// validity callable (a NoData node, which trim would remove; the fix of
// 20-start-quality.md); the walk crossing a constrained edge or leaving the
// mesh; the node already a vertex; the node on a frozen edge of the triangle it
// lies in (docs/increments/23-basin-scale.md, N4). The walk is bounded by
// triangle_count() steps.
//
// Constraint feet (docs/increments/20c-soft-quality.md, R2). With
// constraint_feet, a located node closer than min(dx, dy) / 2 to a constraint
// (constraint_foot.hpp) goes in as its foot when that is a Hit the validity
// callable accepts, and is skipped (skipped_near_line) for any other answer
// but None. The callable is asked about the foot as a MeshVertex.
//
// The soft criterion (20c R7). With min_gain_deg >= 0, a candidate (node or
// foot) goes in only if the worst angle of the triangles it would make, capped
// at theta, is at least the worst of those it would replace plus the gain,
// less kGainSlackDeg; otherwise skipped_no_gain. The triangles come from
// quality_cavity, read-only, with Lawson's own flip test, so they are what
// legalise_around then writes. The line split (20c R8). With the criterion
// and constraint_feet on, a walk blocked by a constrained, non-frozen edge
// with a triangle beyond splits that edge at the node's foot, if the foot is
// at least min(dx, dy) / 2 from both ends, fits, and is valid; not judged by
// R7; otherwise skipped_blocked.
//
// Serial and deterministic: the queue key is (ratio descending, slot
// ascending), and an entry whose slot no longer holds its three vertices is
// stale and dropped. It ends because every insertion is a new node of a finite
// lattice and entries are added only for written slots (R10).
//
// Depends on core, predicates, lattice_mesh.hpp and lawson.hpp; knows no
// raster and reads no height. The grid's rows and cols come in as numbers, and
// whether a node holds data comes in as a callable.

#include <terrain/core/point.hpp>
#include <terrain/mesh/constraint_foot.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <concepts>
#include <cstdint>
#include <limits>
#include <numbers>
#include <optional>
#include <queue>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

namespace terrain::mesh {

// R7's slack s, degrees: a candidate is accepted when new >= old + gain - s,
// so a foot that keeps its edge's end angle exactly is not decided by its
// last bit (20c R7, "Why the slack").
inline constexpr double kGainSlackDeg = 1e-6;

struct QualityOptions {
    double min_angle_deg = 0.0;  // <= 0 (or NaN): the pass does nothing
    std::size_t rows = 0;        // the node rectangle, [0, cols-1] x [0, rows-1]
    std::size_t cols = 0;
    bool constraint_feet = false;  // R2.6: a node near a constraint goes in as its foot
    double min_gain_deg = -1.0;    // 20c R7's gain P; negative (or NaN): off, and R8 with it
};

struct QualityOutcome {
    std::size_t inserted = 0;         // DEM nodes added
    std::size_t skipped_floor = 0;    // circumradius below the floor
    std::size_t skipped_outside = 0;  // circumcentre outside the node rectangle
    std::size_t skipped_vertex = 0;   // the snapped node is already a vertex
    std::size_t skipped_blocked = 0;  // the walk met a constrained edge or left the mesh
    std::size_t walk_bound_hits = 0;  // the walk took triangle_count() steps
    std::size_t skipped_frozen = 0;   // the snapped node lies on a frozen edge
    std::size_t skipped_void = 0;     // the snapped node is not valid (NoData)
    std::size_t feet = 0;               // feet inserted instead of the node, not in `inserted`
    std::size_t skipped_near_line = 0;  // near a constraint with no usable foot
    std::size_t skipped_no_gain = 0;    // R7: the candidate would lower the worst angle
    std::size_t line_splits = 0;        // R8: lines split at a node's foot, in no count above
};

// The default validity callable: every node holds data.
struct AllNodesValid {
    constexpr bool operator()(const MeshVertex&) const noexcept { return true; }
};

namespace detail {

struct QualityShape {
    double ratio;   // circumradius / shortest edge
    double radius;  // circumradius, world units
    Point2 centre;  // circumcentre, world frame
};

[[nodiscard]] inline QualityShape quality_shape(const LatticeMesh& m, std::uint32_t t,
                                                const LatticeFrame& f) {
    const Point2 a = f.at(m.corner(t, 0)), b = f.at(m.corner(t, 1)), c = f.at(m.corner(t, 2));
    const double bx = b.x - a.x, by = b.y - a.y, cx = c.x - a.x, cy = c.y - a.y;
    const double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    const double d = 2.0 * (bx * cy - by * cx);
    const double ux = (cy * b2 - by * c2) / d, uy = (bx * c2 - cx * b2) / d;
    const double radius = std::hypot(ux, uy);
    const double shortest = std::min({std::sqrt(b2), std::sqrt(c2), std::hypot(c.x - b.x, c.y - b.y)});
    return QualityShape{radius / shortest, radius, Point2{a.x + ux, a.y + uy}};
}

struct QualityEntry {
    double ratio;
    std::uint32_t slot;
    TriangleIndices tri;

    // Lower priority: smaller ratio, then larger slot.
    friend bool operator<(const QualityEntry& a, const QualityEntry& b) noexcept {
        return a.ratio < b.ratio || (a.ratio == b.ratio && a.slot > b.slot);
    }
};

// R7's prediction: the slots inserting p would replace and the triangles it
// would make, p joined to each boundary edge, counter-clockwise. `on` 0..2:
// p goes in on t's edge `on` and the triangle across seeds too; 3: inside t.
// Grown across unconstrained edges by quad_flips, Lawson's own test.
struct QualityCavity {
    std::vector<std::uint32_t> removed;
    std::vector<std::array<MeshVertex, 3>> created;
};

// Into c, whose buffers are reused from call to call.
template <pred::GeometryKernel K>
void quality_cavity(const LatticeMesh& m, std::uint32_t t, unsigned on, MeshVertex p, const LatticeFrame& f,
                    QualityCavity& c) {
    c.removed.assign(1, t);
    c.created.clear();
    if (on < 3 && m.neighbours(t)[on] != kNoNeighbour)
        c.removed.push_back(m.neighbours(t)[on]);
    // A few slots: a loop, not std::find, which libc++ hands to wmemchr.
    const auto in = [&](std::uint32_t u) { return std::ranges::any_of(c.removed, [u](auto r) { return r == u; }); };
    for (std::size_t i = 0; i < c.removed.size(); ++i)
        for (unsigned k = 0; k < 3; ++k) {
            const std::uint32_t r = c.removed[i], u = m.neighbours(r)[k];
            if (u == kNoNeighbour || in(u) || m.is_constrained(r, k))
                continue;
            unsigned j = 0;
            while (m.triangles()[u][j] != m.triangles()[r][(k + 1) % 3])
                ++j;
            if (quad_flips<K>(m.corner(r, k), m.corner(r, (k + 1) % 3), p, m.corner(u, (j + 2) % 3), f))
                c.removed.push_back(u);
        }
    for (const auto r : c.removed)
        for (unsigned k = 0; k < 3; ++k)
            if (const auto u = m.neighbours(r)[k]; u == kNoNeighbour ? !(r == t && k == on) : !in(u))
                c.created.push_back({m.corner(r, k), m.corner(r, (k + 1) % 3), p});
}

template <pred::GeometryKernel K>
[[nodiscard]] QualityCavity quality_cavity(const LatticeMesh& m, std::uint32_t t, unsigned on, MeshVertex p,
                                           const LatticeFrame& f) {
    QualityCavity c;
    quality_cavity<K>(m, t, on, p, f, c);
    return c;
}

// The smallest angle of (a, b, c), degrees.
[[nodiscard]] inline double smallest_angle_deg(Point2 a, Point2 b, Point2 c) {
    const std::array<Point2, 3> v{a, b, c};
    double best = std::numbers::pi;
    for (std::size_t k = 0; k < 3; ++k) {
        const double ux = v[(k + 1) % 3].x - v[k].x, uy = v[(k + 1) % 3].y - v[k].y;
        const double wx = v[(k + 2) % 3].x - v[k].x, wy = v[(k + 2) % 3].y - v[k].y;
        best = std::min(best, std::atan2(std::abs(ux * wy - uy * wx), ux * wx + uy * wy));
    }
    return best * 180.0 / std::numbers::pi;
}

// R7's constants and scratch: theta, the gain, theta's cosine and sine, the cavity.
struct GainJudge {
    double theta, gain, cos_t = std::cos(theta * std::numbers::pi / 180.0),
                        sin_t = std::sin(theta * std::numbers::pi / 180.0);
    QualityCavity c{};

    // min(theta, smallest_angle_deg(a, b, c)), with the atan2 calls skipped
    // when every corner is above theta by a margin: corner angle phi > theta
    // iff |cross| cos - dot sin = r sin(phi - theta) > 0, and the margin,
    // 1e-9 of |cross| + |dot| >= |u||w|, is far above the rounding of these
    // terms, of atan2 and of the degree conversion, so the skip leaves the
    // value smallest_angle_deg's min with theta gives.
    [[nodiscard]] double capped(Point2 a, Point2 b, Point2 c) const {
        const std::array<Point2, 3> v{a, b, c};
        for (std::size_t k = 0; k < 3; ++k) {
            const double ux = v[(k + 1) % 3].x - v[k].x, uy = v[(k + 1) % 3].y - v[k].y;
            const double wx = v[(k + 2) % 3].x - v[k].x, wy = v[(k + 2) % 3].y - v[k].y;
            const double cross = std::abs(ux * wy - uy * wx), dot = ux * wx + uy * wy;
            if (!(cross * cos_t - dot * sin_t > 1e-9 * (cross + std::abs(dot))))
                return std::min(theta, smallest_angle_deg(a, b, c));
        }
        return theta;
    }
};

// R7: inserting p at (t, on) leaves the worst angle, capped at theta, no lower
// than before plus gain, less the slack. The worst before only lowers the
// bar, so (by monotone rounding) the answer is yes once the worst after
// clears it at theta, or at any worst-before-so-far.
template <pred::GeometryKernel K>
[[nodiscard]] bool pays(const LatticeMesh& m, std::uint32_t t, unsigned on, MeshVertex p, const LatticeFrame& f,
                        GainJudge& j) {
    quality_cavity<K>(m, t, on, p, f, j.c);
    double old_w = j.theta, new_w = j.theta;
    for (const auto& x : j.c.created)
        new_w = std::min(new_w, j.capped(f.at(x[0]), f.at(x[1]), f.at(x[2])));
    for (const auto r : j.c.removed) {
        if (new_w >= old_w + j.gain - kGainSlackDeg)
            return true;
        old_w = std::min(old_w, j.capped(f.at(m.corner(r, 0)), f.at(m.corner(r, 1)), f.at(m.corner(r, 2))));
    }
    return new_w >= old_w + j.gain - kGainSlackDeg;
}

}  // namespace detail

template <pred::GeometryKernel K, std::predicate<const MeshVertex&> Valid = AllNodesValid>
QualityOutcome improve(LatticeMesh& m, const LatticeFrame& f, const QualityOptions& o,
                       const Valid& valid = {}) {
    QualityOutcome out;
    if (!(o.min_angle_deg > 0.0) || o.rows == 0 || o.cols == 0)
        return out;
    const double bound = 1.0 / (2.0 * std::sin(o.min_angle_deg * std::numbers::pi / 180.0));
    const double floor = std::hypot(f.dx, f.dy);

    std::priority_queue<detail::QualityEntry> queue;
    const auto offer = [&](std::uint32_t t) {
        if (const double r = detail::quality_shape(m, t, f).ratio; r > bound)
            queue.push({r, t, m.triangles()[t]});
    };
    for (std::uint32_t t = 0; t < m.triangle_count(); ++t)
        offer(t);

    const bool judged = o.min_gain_deg >= 0.0, split_lines = judged && o.constraint_feet;
    const double half_cell = std::min(f.dx, f.dy) / 2.0;  // delta_q, R2.2
    detail::GainJudge judge{o.min_angle_deg, o.min_gain_deg};
    std::vector<std::uint32_t> written;
    // p into t (`on` 3) or onto t's edge `on`, legalised; written slots offered.
    const auto insert = [&](std::uint32_t t, unsigned on, MeshVertex p) {
        const auto before = static_cast<std::uint32_t>(m.triangle_count());
        written.assign({t, before, before + 1});
        std::uint32_t q = 0;
        if (on == 3) {
            q = m.split_inside(t, p);
        } else {
            const auto u = m.neighbours(t)[on];
            q = m.split_edge(t, on, p);
            if (u == kNoNeighbour)
                written.pop_back();
            else
                written.push_back(u);
        }
        const std::vector<std::uint32_t> seeds = written;
        legalise_around<K>(m, q, std::span<const std::uint32_t>{seeds}, f,
                           [&](std::uint32_t w) { written.push_back(w); });
        std::sort(written.begin(), written.end());
        written.erase(std::unique(written.begin(), written.end()), written.end());
        for (const auto w : written)
            offer(w);
    };
    while (!queue.empty()) {
        const detail::QualityEntry e = queue.top();
        queue.pop();
        if (m.triangles()[e.slot] != e.tri)
            continue;  // stale: the slot was rewritten, and its new triangle was offered
        const auto s = detail::quality_shape(m, e.slot, f);
        if (s.radius < floor) {
            ++out.skipped_floor;
            continue;
        }
        const double col = s.centre.x / f.dx, row = -s.centre.y / f.dy;
        if (!(col >= 0.0 && row >= 0.0 && col <= static_cast<double>(o.cols - 1)
              && row <= static_cast<double>(o.rows - 1))) {
            ++out.skipped_outside;
            continue;
        }
        const LatticeVertex node{static_cast<std::uint32_t>(std::round(row)),
                                 static_cast<std::uint32_t>(std::round(col))};
        if (!valid(node)) {
            ++out.skipped_void;
            continue;
        }
        const MeshVertex p{node};

        // Visibility walk: cross the first edge with p strictly beyond it.
        std::uint32_t t = e.slot;
        std::size_t steps = 0;
        unsigned hit = 0;  // the edge a blocked walk stopped at
        bool blocked = false, located = false;
        while (!blocked && !located) {
            if (++steps > m.triangle_count()) {
                ++out.walk_bound_hits;
                break;
            }
            located = true;
            for (unsigned k = 0; k < 3; ++k)
                if (orient_sign(m.corner(t, k), m.corner(t, (k + 1) % 3), p) < 0) {
                    located = false;
                    const auto u = m.neighbours(t)[k];
                    blocked = u == kNoNeighbour || m.is_constrained(t, k);
                    std::tie(t, hit) = std::pair{blocked ? t : u, k};
                    break;
                }
        }
        if (blocked) {  // R8: split the line at the node's foot, or skip
            std::optional<FootSearch> s;
            if (split_lines && m.neighbours(t)[hit] != kNoNeighbour)
                s = detail::foot_on(m, t, hit, p, std::numeric_limits<double>::infinity(), half_cell, f);
            if (s && s->status == FootStatus::Hit && valid(s->at)) {
                insert(t, hit, s->at);
                ++out.line_splits;
            } else {
                ++out.skipped_blocked;
            }
        }
        if (!located)
            continue;

        std::array<int, 3> side{};
        for (unsigned k = 0; k < 3; ++k)
            side[k] = orient_sign(m.corner(t, k), m.corner(t, (k + 1) % 3), p);
        const auto zeros = std::count(side.begin(), side.end(), 0);
        if (zeros > 1) {
            ++out.skipped_vertex;
            continue;
        }
        const auto on = static_cast<unsigned>(std::find(side.begin(), side.end(), 0) - side.begin());
        if (zeros == 1 && m.is_frozen(t, on)) {
            ++out.skipped_frozen;
            continue;
        }
        std::uint32_t owner = t;
        unsigned edge = zeros == 0 ? 3 : on;
        MeshVertex at = p;
        bool foot = false;
        if (o.constraint_feet) {
            const auto s = constraint_foot(m, t, p, half_cell, f);
            if (s.status != FootStatus::None) {
                if (s.status != FootStatus::Hit || !valid(s.at)) {
                    ++out.skipped_near_line;
                    continue;
                }
                foot = true;
                std::tie(owner, edge, at) = std::tuple{s.owner, s.edge, s.at};
            }
        }
        if (judged && !detail::pays<K>(m, owner, edge, at, f, judge)) {
            ++out.skipped_no_gain;
            continue;
        }
        insert(owner, edge, at);
        ++(foot ? out.feet : out.inserted);
    }
    return out;
}

}  // namespace terrain::mesh
