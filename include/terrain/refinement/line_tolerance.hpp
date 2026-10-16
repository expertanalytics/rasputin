#pragma once

// A vertical tolerance that varies with distance to named lines
// (docs/increments/33-feature-tolerance.md, sections 3 and 4.2, pins in 9.1).
//
// A tolerance policy says what each triangle is allowed. UniformTolerance is
// today's single number. LineTolerance is the ramp t(d) of section 3 at the
// triangle's distance d(T) to a set of segments: the smallest distance from
// any point of the closed triangle to any segment, less the margin the lines
// were simplified by. t never decreases with d and d(T) <= d(n) for every
// node n of T, so a triangle within t(d(T)) holds every node to its own t.
//
// The distance is a threshold, not a topology decision: plain doubles, in the
// lattice-metre frame (col dx, -row dy), so no UTM offset enters it. It is a
// pure function of one triangle and one segment, so the result does not
// depend on the thread count or on the index's visiting order.
//
// The search (4.2): noding::BroadPhase is queried with the triangle's box
// grown by g = E/16, E/8, E/4, E/2, E, E + margin. Every segment within g of
// the triangle has its box in the grown box, so once the best found is at
// most g it is the minimum. With nothing within E + margin the distance is
// infinite (9.1, A): its ramp value is F, also for a step ramp.

#include <terrain/core/bbox.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/segment.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/noding/broad_phase.hpp>
#include <terrain/raster/geometry.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace terrain::refinement {

struct ToleranceRamp {  // metres; 0 <= near <= far, 0 <= start <= end, all finite
    double near, far, start, end;
    double margin = 0.0;  // subtracted from every distance, >= 0

    // Section 3's t(max(0, d - margin)); a step (start == end) divides by nothing.
    [[nodiscard]] double at(double distance) const noexcept {
        const double x = std::max(0.0, distance - margin);
        if (x <= start)
            return near;
        if (x >= end)
            return far;
        return near + (far - near) * (x - start) / (end - start);
    }
};

// lowest() and highest() bound at() over every triangle; refine uses them to
// skip a query (4.4).
template <class P>
concept TolerancePolicy = requires(const P& p, const mesh::LatticeMesh& m, std::uint32_t t) {
    { p.lowest() } -> std::convertible_to<double>;
    { p.highest() } -> std::convertible_to<double>;
    { p.at(m, t) } -> std::convertible_to<double>;
};

struct UniformTolerance {  // today's behaviour
    double value;
    [[nodiscard]] double lowest() const noexcept { return value; }
    [[nodiscard]] double highest() const noexcept { return value; }
    [[nodiscard]] double at(const mesh::LatticeMesh&, std::uint32_t) const noexcept { return value; }
};

class LineTolerance {
public:
    // nullopt with one plain sentence in `why` for a ramp value that breaks its
    // bound or is not finite, or a non-finite coordinate. Zero segments is a
    // field too: every triangle gets far. Segments are x0 y0 x1 y1, world.
    [[nodiscard]] static std::optional<LineTolerance> make(const raster::RasterGeometry& g,
                                                           std::span<const std::array<double, 4>> segments,
                                                           ToleranceRamp ramp, std::string& why) {
        const std::pair<const char*, double> values[] = {
            {"near", ramp.near}, {"far", ramp.far}, {"start", ramp.start}, {"end", ramp.end}, {"margin", ramp.margin}};
        for (const auto& [name, v] : values)
            if (!std::isfinite(v)) {
                why = std::string{name} + " must be finite";
                return std::nullopt;
            }
        if (ramp.near < 0.0 || ramp.near > ramp.far)
            why = "near must be >= 0 and at most far";
        else if (ramp.start < 0.0 || ramp.start > ramp.end)
            why = "start must be >= 0 and at most end";
        else if (ramp.margin < 0.0)
            why = "margin must be >= 0";
        if (!why.empty())
            return std::nullopt;
        std::vector<Segment2> frame;
        frame.reserve(segments.size());
        for (const auto& s : segments) {
            if (!std::ranges::all_of(s, [](double v) { return std::isfinite(v); })) {
                why = "every segment coordinate must be finite";
                return std::nullopt;
            }
            frame.push_back(Segment2{Point2{s[0] - g.x_min(), s[1] - g.y_max()},
                                     Point2{s[2] - g.x_min(), s[3] - g.y_max()}});
        }
        return LineTolerance{g, ramp, std::move(frame)};
    }

    [[nodiscard]] double lowest() const noexcept { return segments_.empty() ? ramp_.far : ramp_.near; }
    [[nodiscard]] double highest() const noexcept { return ramp_.far; }
    [[nodiscard]] const raster::RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] const ToleranceRamp& ramp() const noexcept { return ramp_; }

    [[nodiscard]] double at(const mesh::LatticeMesh& m, std::uint32_t t) const { return ramp_.at(search(m, t)); }

    // d(T) to the segments as given (margin not subtracted), capped at end + margin.
    [[nodiscard]] double distance(const mesh::LatticeMesh& m, std::uint32_t t) const {
        return std::min(search(m, t), ramp_.end + ramp_.margin);
    }

private:
    LineTolerance(const raster::RasterGeometry& g, ToleranceRamp ramp, std::vector<Segment2> segments)
        : geometry_{g}, ramp_{ramp}, segments_{std::move(segments)}, index_{segments_} {}

    static double orient(Point2 a, Point2 b, Point2 c) noexcept {
        return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    }

    static double point_segment(Point2 p, Point2 a, Point2 b) noexcept {
        const double ux = b.x - a.x, uy = b.y - a.y, len2 = ux * ux + uy * uy;
        const double s = len2 > 0.0 ? std::clamp(((p.x - a.x) * ux + (p.y - a.y) * uy) / len2, 0.0, 1.0) : 0.0;
        return std::hypot(p.x - (a.x + s * ux), p.y - (a.y + s * uy));
    }

    // Closed segments ab and pq meet; four zero orientations are collinear, decided by the boxes.
    static bool meet(Point2 a, Point2 b, Point2 p, Point2 q) noexcept {
        const double o1 = orient(a, b, p), o2 = orient(a, b, q), o3 = orient(p, q, a), o4 = orient(p, q, b);
        if (o1 == 0.0 && o2 == 0.0 && o3 == 0.0 && o4 == 0.0)
            return std::max(std::min(a.x, b.x), std::min(p.x, q.x)) <= std::min(std::max(a.x, b.x), std::max(p.x, q.x))
                   && std::max(std::min(a.y, b.y), std::min(p.y, q.y)) <= std::min(std::max(a.y, b.y), std::max(p.y, q.y));
        return !((o1 > 0.0 && o2 > 0.0) || (o1 < 0.0 && o2 < 0.0) || (o3 > 0.0 && o4 > 0.0) || (o3 < 0.0 && o4 < 0.0));
    }

    // The closed CCW triangle v to the segment s: 0 when an end is inside or s meets an edge.
    static double pair(const std::array<Point2, 3>& v, const Segment2& s) noexcept {
        const auto inside = [&](Point2 p) {
            return orient(v[0], v[1], p) >= 0.0 && orient(v[1], v[2], p) >= 0.0 && orient(v[2], v[0], p) >= 0.0;
        };
        if (inside(s.a) || inside(s.b))
            return 0.0;
        double d = std::numeric_limits<double>::infinity();
        for (unsigned k = 0; k < 3; ++k) {
            const Point2 a = v[k], b = v[(k + 1) % 3];
            if (meet(a, b, s.a, s.b))
                return 0.0;
            d = std::min({d, point_segment(s.a, a, b), point_segment(s.b, a, b), point_segment(a, s.a, s.b)});
        }
        return d;
    }

    // The uncapped d(T), infinity when no segment lies within end + margin.
    [[nodiscard]] double search(const mesh::LatticeMesh& m, std::uint32_t t) const {
        constexpr double inf = std::numeric_limits<double>::infinity();
        if (segments_.empty())
            return inf;
        std::array<Point2, 3> v;
        for (unsigned k = 0; k < 3; ++k) {
            const mesh::MeshVertex c = m.corner(t, k);
            v[k] = Point2{c.col * geometry_.delta_x(), -(c.row * geometry_.delta_y())};
        }
        const double x0 = std::min({v[0].x, v[1].x, v[2].x}), x1 = std::max({v[0].x, v[1].x, v[2].x});
        const double y0 = std::min({v[0].y, v[1].y, v[2].y}), y1 = std::max({v[0].y, v[1].y, v[2].y});
        const double reach = ramp_.end + ramp_.margin;
        double best = inf;
        for (const double g : {ramp_.end / 16, ramp_.end / 8, ramp_.end / 4, ramp_.end / 2, ramp_.end, reach}) {
            index_.for_each_candidate(Box2{Point2{x0 - g, y0 - g}, Point2{x1 + g, y1 + g}},
                                      [&](std::uint32_t i) { best = std::min(best, pair(v, segments_[i])); });
            if (best <= g)
                return best;
        }
        return inf;
    }

    raster::RasterGeometry geometry_;
    ToleranceRamp ramp_;
    std::vector<Segment2> segments_;  // in the lattice-metre frame (x - x_min, y - y_max)
    noding::BroadPhase index_;
};

}  // namespace terrain::refinement
