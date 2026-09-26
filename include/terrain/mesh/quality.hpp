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
// Skips, in R4's order: circumradius below the floor sqrt(dx^2 + dy^2), so the
// node is within R / 2 of the centre (R5); centre outside the node rectangle
// (never clamped); the node already a vertex; the walk crossing a constrained
// edge or leaving the mesh. The walk is bounded by triangle_count() steps.
//
// Serial and deterministic: the queue key is (ratio descending, slot
// ascending), and an entry whose slot no longer holds its three vertices is
// stale and dropped. It ends because every insertion is a new node of a finite
// lattice and entries are added only for written slots (R10).
//
// Depends on core, predicates, lattice_mesh.hpp and lawson.hpp; knows no
// raster and reads no height. The grid's rows and cols come in as numbers.

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/predicates/kernel.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <queue>
#include <span>
#include <vector>

namespace terrain::mesh {

struct QualityOptions {
    double min_angle_deg = 0.0;  // <= 0 (or NaN): the pass does nothing
    std::size_t rows = 0;        // the node rectangle, [0, cols-1] x [0, rows-1]
    std::size_t cols = 0;
};

struct QualityOutcome {
    std::size_t inserted = 0;         // DEM nodes added
    std::size_t skipped_floor = 0;    // circumradius below the floor
    std::size_t skipped_outside = 0;  // circumcentre outside the node rectangle
    std::size_t skipped_vertex = 0;   // the snapped node is already a vertex
    std::size_t skipped_blocked = 0;  // the walk met a constrained edge or left the mesh
    std::size_t walk_bound_hits = 0;  // the walk took triangle_count() steps
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

}  // namespace detail

template <pred::GeometryKernel K>
QualityOutcome improve(LatticeMesh& m, const LatticeFrame& f, const QualityOptions& o) {
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

    std::vector<std::uint32_t> written;
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
        const MeshVertex p{node};

        // Visibility walk: cross the first edge with p strictly beyond it.
        std::uint32_t t = e.slot;
        std::size_t steps = 0;
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
                    t = blocked ? t : u;
                    break;
                }
        }
        if (blocked)
            ++out.skipped_blocked;
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
        const auto before = static_cast<std::uint32_t>(m.triangle_count());
        written.assign({t, before, before + 1});
        std::uint32_t q = 0;
        if (zeros == 0) {
            q = m.split_inside(t, node);
        } else {
            const auto k = static_cast<unsigned>(std::find(side.begin(), side.end(), 0) - side.begin());
            const auto u = m.neighbours(t)[k];
            q = m.split_edge(t, k, node);
            if (u == kNoNeighbour)
                written.pop_back();
            else
                written.push_back(u);
        }
        ++out.inserted;
        const std::vector<std::uint32_t> seeds = written;
        legalise_around<K>(m, q, std::span<const std::uint32_t>{seeds}, f,
                           [&](std::uint32_t w) { written.push_back(w); });
        std::sort(written.begin(), written.end());
        written.erase(std::unique(written.begin(), written.end()), written.end());
        for (const auto w : written)
            offer(w);
    }
    return out;
}

}  // namespace terrain::mesh
