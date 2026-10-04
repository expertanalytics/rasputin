#pragma once

// The seam pass (docs/increments/23-basin-scale.md, "The seam protocol" step 2,
// N8-N12): the vertices a frozen seam edge gets before either piece refines.
//
// The edge is ordered so a < b by (x, y), whichever order it was given in, and
// handed to the edge strip's generator (constraint_points.hpp) as vertices
// {a, b}, so its check points and their parameters s run from a and do not
// depend on the order given. Each end's z is vertex_z at its lattice position,
// none on NoData.
//
// One-dimensional greedy, Douglas-Peucker in the vertical. The polyline starts
// as (a, b). A piece with an end that has no z has no lerp: while it holds a
// check point, the one nearest its invalid end is inserted (N10; the points
// are sorted by s, so that is its first or last, and each invalid end costs at
// most one). Then, while some check point on a valid piece is more than
// `tolerance` from the lerp between the piece's ends (in s, as the strip's
// `along`), the worst goes in, ties to the smallest s. Each piece keeps its
// own worst point in a queue, so the global worst is the queue's top; a piece
// is split only at that point. It ends: each point is inserted at most once.
//
// An inserted point is output at (x_min + col dx, y_max - row dy) of its
// lattice position, with the check point's own z and s. max_error is the
// largest error left on a valid piece. Pure: two pieces that call it with one
// raster and one edge get the same output bit for bit (K4).

#include <terrain/core/point.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/scan.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <queue>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace terrain::refinement {

struct SeamPoint {
    Point2 at;  // world
    double z;   // the check point's vertex_z
    double s;   // its parameter from a
};

struct SeamOutcome {
    Point2 a, b;                     // the edge's ends, a < b by (x, y)
    std::optional<double> z_a, z_b;  // vertex_z at each end; none on NoData
    std::vector<SeamPoint> points;   // inserted, s strictly increasing
    std::size_t check_points = 0;    // the generator's points for the edge
    std::size_t no_data = 0;         // check points it dropped for a NoData stencil
    double max_error = 0.0;          // over check points on valid pieces, at the end
};

template <raster::RasterSource R>
[[nodiscard]] SeamOutcome refine_seam(const R& dem, Point2 a, Point2 b, double tolerance) {
    const raster::RasterGeometry& g = dem.geometry();
    if (!std::isfinite(tolerance) || tolerance < 0.0)
        throw std::invalid_argument("refine_seam: tolerance must be finite and >= 0");
    if (a == b)
        throw std::invalid_argument("refine_seam: the edge's two ends are the same point");
    if (!g.cell_of(a) || !g.cell_of(b))
        throw std::invalid_argument("refine_seam: an end is outside the DEM's node rectangle");
    if (std::tie(b.x, b.y) < std::tie(a.x, a.y))
        std::swap(a, b);
    const std::array<Point2, 2> ends{a, b};
    const std::array<std::array<std::uint32_t, 2>, 1> edge{{{0, 1}}};
    const ConstraintCheckPoints cp = constraint_check_points(dem, std::span<const Point2>{ends},
                                                             std::span<const std::array<std::uint32_t, 2>>{edge});
    const std::span<const ConstraintPoint> pts = cp.on_edge(0);
    const auto n = static_cast<std::ptrdiff_t>(pts.size());
    SeamOutcome out{a, b, vertex_z(dem, detail::lattice_position(g, a)), vertex_z(dem, detail::lattice_position(g, b)),
                    {}, cp.size(), cp.no_data(), 0.0};

    // Knot i in [-1, n]: an end (-1 is a, n is b) or an inserted check point.
    const auto s_of = [&](std::ptrdiff_t i) { return i < 0 ? 0.0 : i == n ? 1.0 : pts[static_cast<std::size_t>(i)].s; };
    const auto z_of = [&](std::ptrdiff_t i) {
        return i < 0 ? out.z_a : i == n ? out.z_b : std::optional<double>{pts[static_cast<std::size_t>(i)].z};
    };
    std::vector<char> in(pts.size(), 0);
    if (n > 0 && !out.z_a)
        in.front() = 1;  // N10: the point nearest each invalid end
    if (n > 0 && !out.z_b)
        in.back() = 1;

    struct Worst {
        double error;
        std::ptrdiff_t at, lo, hi;  // the point, and its piece's knots
        bool operator<(const Worst& o) const { return error < o.error || (error == o.error && at > o.at); }
    };
    std::priority_queue<Worst> queue;
    const auto offer = [&](std::ptrdiff_t lo, std::ptrdiff_t hi) {  // the piece (lo, hi)'s worst point
        const auto zl = z_of(lo), zh = z_of(hi);
        if (!zl || !zh || hi - lo < 2)
            return;
        const double sl = s_of(lo), sh = s_of(hi);
        Worst w{-1.0, -1, lo, hi};
        for (std::ptrdiff_t i = lo + 1; i < hi; ++i) {
            const ConstraintPoint& p = pts[static_cast<std::size_t>(i)];
            const double e = std::abs(p.z - (*zl + (p.s - sl) / (sh - sl) * (*zh - *zl)));
            if (e > w.error)
                w = Worst{e, i, lo, hi};
        }
        queue.push(w);
    };
    std::ptrdiff_t lo = -1;
    for (std::ptrdiff_t i = 0; i <= n; ++i)
        if (i == n || in[static_cast<std::size_t>(i)] != 0) {
            offer(lo, i);
            lo = i;
        }
    while (!queue.empty() && queue.top().error > tolerance) {
        const Worst w = queue.top();
        queue.pop();
        in[static_cast<std::size_t>(w.at)] = 1;
        offer(w.lo, w.at);
        offer(w.at, w.hi);
    }
    out.max_error = queue.empty() ? 0.0 : queue.top().error;
    for (std::size_t i = 0; i < pts.size(); ++i)
        if (in[i] != 0)
            out.points.push_back(SeamPoint{
                Point2{g.x_min() + pts[i].at.col * g.delta_x(), g.y_max() - pts[i].at.row * g.delta_y()},
                pts[i].z, pts[i].s});
    return out;
}

}  // namespace terrain::refinement
