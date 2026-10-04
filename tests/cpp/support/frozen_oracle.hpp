#pragma once

// Test-only fixtures and oracles for increment 23b
// (docs/increments/23-basin-scale.md, "The seam protocol", K1 to K4 and
// "Tests @tester can write red", FE1 to FE5 and SP1 to SP4): frozen edges
// and the seam pass. Shared by property/prop_refinement_frozen.cpp and
// property/prop_refinement_seam.cpp.
//
// Independence. Nothing here includes scan.hpp, refine.hpp, refine_points.hpp,
// quality.hpp or seam.hpp. The oracles read the OUTPUT (world vertices, z,
// triangles, constraint edges and masks) and the START the run was given, map
// world points to the lattice by the affine map written in strip_oracle.hpp,
// and decide "on a frozen edge" with the exact orientation on (col, -row): the
// producer's predicate, never its records (computational-geometry skill, §3).
//
// Nothing here asserts: the oracles RETURN their findings, so each can be
// planted and shown to fail (FZ0, unit/test_refinement_frozen_oracle.cpp).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include "strip_oracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace frozen_oracle {

using strip_oracle::Edges;
using strip_oracle::Lat;
using strip_oracle::Mesh;
using strip_oracle::Start;
using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;

inline constexpr std::uint32_t kOutline = 1u;
inline constexpr std::uint32_t kFeature = 2u;
inline constexpr std::uint32_t kSeam = 4u;
inline constexpr std::uint32_t kNoEdgeHasThis = 1u << 31;  // a frozen mask that meets no fixture edge

// ---------------------------------------------------------------------------
// Exact incidence in the lattice frame (col, -row)
// ---------------------------------------------------------------------------

inline Point2 frame(Lat p) { return Point2{p.col, -p.row}; }

inline int orient(Lat a, Lat b, Lat p) {
    return static_cast<int>(terrain::pred::DefaultKernel::orient2d(frame(a), frame(b), frame(p)));
}

// p strictly inside the open segment (a, b): exactly collinear, inside the
// closed box, and neither end.
inline bool on_open(Lat a, Lat b, Lat p) {
    if (p == a || p == b || orient(a, b, p) != 0) return false;
    return std::min(a.col, b.col) <= p.col && p.col <= std::max(a.col, b.col) && std::min(a.row, b.row) <= p.row
        && p.row <= std::max(a.row, b.row);
}

// p within `slack` cells of the open segment, its projection strictly inside
// it by more than `slack` from either end.
inline bool near_open(Lat a, Lat b, Lat p, double slack) {
    const auto [t, d] = strip_oracle::param_dist(a, b, p);
    const double len = std::hypot(b.col - a.col, b.row - a.row);
    return d <= slack && t * len > slack && (1.0 - t) * len > slack;
}

// ---------------------------------------------------------------------------
// K2: frozen means frozen
// ---------------------------------------------------------------------------

struct FrozenFindings {
    std::size_t frozen_edges = 0;  // start edges whose mask meets `frozen` (so a sweep cannot pass vacuously)
    std::size_t missing = 0;       // a frozen start edge that is not an output constraint edge between the same two vertices
    std::size_t wrong_mask = 0;    // ... present, with another mask
    std::size_t moved = 0;         // a frozen edge's end not output where the start had it
    std::size_t on = 0;            // output vertices exactly on a frozen edge's open segment
    std::size_t near = 0;          // output vertices within slack of one, not on it (a split by rounding)
};

inline FrozenFindings frozen_findings(const RasterGeometry& g, std::span<const Point2> start_vertices,
                                      const Edges& start_edges, const std::vector<std::uint32_t>& start_masks,
                                      std::uint32_t frozen, const Mesh& out) {
    FrozenFindings f;
    const double slack = strip_oracle::on_edge_slack(g);
    std::vector<Lat> lv;
    for (const Point2 v : out.vertices) lv.push_back(strip_oracle::lat(g, v));
    for (std::size_t k = 0; k < start_edges.size(); ++k) {
        if ((start_masks[k] & frozen) == 0) continue;
        ++f.frozen_edges;
        const auto [i, j] = std::minmax(start_edges[k][0], start_edges[k][1]);
        for (const auto v : {i, j})
            if (v >= out.vertices.size() || !(out.vertices[v] == start_vertices[v])) ++f.moved;
        bool found = false;
        for (std::size_t e = 0; e < out.edges.size(); ++e)
            // minmax returns a pair of references; libstdc++ (GCC 13) has no
            // == between pair<const T&, const T&> and pair<T, T>, so compare values.
            if (const auto [lo, hi] = std::minmax(out.edges[e][0], out.edges[e][1]); lo == i && hi == j) {
                found = true;
                if (out.masks[e] != start_masks[k]) ++f.wrong_mask;
            }
        if (!found) ++f.missing;
        const Lat a = strip_oracle::lat(g, start_vertices[i]), b = strip_oracle::lat(g, start_vertices[j]);
        for (std::size_t v = 0; v < lv.size(); ++v) {
            if (v == i || v == j) continue;
            if (on_open(a, b, lv[v]))
                ++f.on;
            else if (near_open(a, b, lv[v], slack))
                ++f.near;
        }
    }
    return f;
}

inline bool clean(const FrozenFindings& f) {
    return f.missing == 0 && f.wrong_mask == 0 && f.moved == 0 && f.on == 0 && f.near == 0;
}

// ---------------------------------------------------------------------------
// The tolerance oracle with the seam's nodes left to the seam pass
// ---------------------------------------------------------------------------

// The frozen start segments in the lattice.
inline std::vector<std::array<Lat, 2>> frozen_segments(const RasterGeometry& g, std::span<const Point2> vertices,
                                                        const Edges& edges, const std::vector<std::uint32_t>& masks,
                                                        std::uint32_t frozen) {
    std::vector<std::array<Lat, 2>> out;
    for (std::size_t k = 0; k < edges.size(); ++k)
        if ((masks[k] & frozen) != 0)
            out.push_back({strip_oracle::lat(g, vertices[edges[k][0]]), strip_oracle::lat(g, vertices[edges[k][1]])});
    return out;
}

// strip_oracle::node_findings, except that a DEM node lying exactly on the
// open segment of a frozen edge is not checked (the design's FE2: "a
// brute-force oracle over all nodes not on frozen edges"). Every other valid
// node in every closed output triangle with three valid vertices is measured
// against that triangle's plane, recomputed from the output's z.
inline strip_oracle::NodeFindings node_findings_off_frozen(const Raster<float>& dem, const Mesh& out, double tol,
                                                           const std::vector<std::array<Lat, 2>>& frozen) {
    using terrain::pred::DefaultKernel;
    using terrain::pred::Orientation;
    const RasterGeometry& g = dem.geometry();
    std::vector<Point2> fp;
    double zmax = 0.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        fp.push_back(frame(strip_oracle::lat(g, out.vertices[i])));
        if (out.valid[i]) zmax = std::max(zmax, std::abs(out.z[i]));
    }
    strip_oracle::NodeFindings f;
    for (const auto& tri : out.triangles) {
        const Point2 a = fp[tri[0]], b = fp[tri[1]], c = fp[tri[2]];
        if (DefaultKernel::orient2d(a, b, c) != Orientation::CounterClockwise) {
            ++f.not_ccw;
            continue;
        }
        if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
        const double two_a = strip_oracle::cross(a, b, c);
        const auto lo_c = static_cast<std::int64_t>(std::ceil(std::min({a.x, b.x, c.x})));
        const auto hi_c = static_cast<std::int64_t>(std::floor(std::max({a.x, b.x, c.x})));
        const auto lo_r = static_cast<std::int64_t>(std::ceil(-std::max({a.y, b.y, c.y})));
        const auto hi_r = static_cast<std::int64_t>(std::floor(-std::min({a.y, b.y, c.y})));
        for (std::int64_t r = std::max<std::int64_t>(lo_r, 0); r <= hi_r; ++r)
            for (std::int64_t col = std::max<std::int64_t>(lo_c, 0); col <= hi_c; ++col) {
                const Point2 p{static_cast<double>(col), -static_cast<double>(r)};
                if (DefaultKernel::orient2d(a, b, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(b, c, p) == Orientation::Clockwise
                    || DefaultKernel::orient2d(c, a, p) == Orientation::Clockwise)
                    continue;
                const auto ur = static_cast<std::size_t>(r), uc = static_cast<std::size_t>(col);
                if (p == a || p == b || p == c || strip_oracle::nodata(dem, ur, uc)) continue;
                const Lat node{static_cast<double>(col), static_cast<double>(r)};
                if (std::any_of(frozen.begin(), frozen.end(),
                                [&](const auto& s) { return on_open(s[0], s[1], node); }))
                    continue;
                const double plane = (strip_oracle::cross(p, b, c) * out.z[tri[0]]
                                      + strip_oracle::cross(a, p, c) * out.z[tri[1]]
                                      + strip_oracle::cross(a, b, p) * out.z[tri[2]])
                                   / two_a;
                const double err = std::abs(strip_oracle::at(dem, ur, uc) - plane);
                f.worst = std::max(f.worst, err);
                if (err > tol + 1e-9 * std::max(1.0, zmax)) ++f.over;
            }
    }
    return f;
}

// ---------------------------------------------------------------------------
// A domain cut by one seam
// ---------------------------------------------------------------------------

enum class Side { Both, Left, Right };

// The rectangle [left, left + w] x [0, h] in (col, row), cut by a seam from
// T = (top, 0) to B = (bottom, h), with `inner` (lattice positions strictly inside the
// seam, ordered from T to B) as its vertices: the seam pass's points, already
// in the start (the design's step 3, "The start slice, with fans").
//
// Vertices, Both: 0 TL, 1 T, 2 TR, 3 BR, 4 B, 5 BL,
// then the inner points in order. Left is the polygon TL, BL, B, inner
// reversed, T; Right is T, inner, B, BR, TR; each is two triangles with the
// one holding the seam fanned around its far corner (BL on the left, BR on the
// right), so the fans are the design's "(a, p1, c), (p1, p2, c), ..., (pk, b, c)".
// Left and Right number their vertices in first-use order, as a slice does.
//
// Outline edges carry kOutline, seam sub-edges `seam_mask` (kSeam, or
// kSeam | kFeature for "a seam along a feature edge").
inline Start seam_domain(const RasterGeometry& g, double w, double h, double top, double bottom,
                         const std::vector<Lat>& inner, Side side, std::uint32_t seam_mask = kSeam,
                         double left = 0.0) {
    const std::vector<Lat> named{{left, 0.0}, {top, 0.0}, {left + w, 0.0}, {left + w, h}, {bottom, h}, {left, h}};
    enum : std::uint32_t { TL, T, TR, BR, B, BL };
    std::vector<std::uint32_t> chain{T};  // the seam from T to B, as named indices (inner from 6)
    for (std::uint32_t i = 0; i < inner.size(); ++i) chain.push_back(6 + i);
    chain.push_back(B);
    auto position = [&](std::uint32_t v) { return v < 6 ? named[v] : inner[v - 6]; };

    std::vector<terrain::TriangleIndices> named_tris;
    std::vector<std::pair<std::array<std::uint32_t, 2>, std::uint32_t>> named_edges;
    if (side != Side::Right) {
        named_tris.push_back({TL, BL, T});
        for (std::size_t i = chain.size() - 1; i > 0; --i)  // B -> ... -> T, around BL
            named_tris.push_back({chain[i], chain[i - 1], BL});
        named_edges.push_back({{TL, BL}, kOutline});
        named_edges.push_back({{BL, B}, kOutline});
        named_edges.push_back({{T, TL}, kOutline});
    }
    if (side != Side::Left) {
        named_tris.push_back({T, BR, TR});
        for (std::size_t i = 0; i + 1 < chain.size(); ++i)  // T -> ... -> B, around BR
            named_tris.push_back({chain[i], chain[i + 1], BR});
        named_edges.push_back({{B, BR}, kOutline});
        named_edges.push_back({{BR, TR}, kOutline});
        named_edges.push_back({{TR, T}, kOutline});
    }
    for (std::size_t i = 0; i + 1 < chain.size(); ++i) named_edges.push_back({{chain[i], chain[i + 1]}, seam_mask});

    // Renumber in first-use order over the triangles.
    std::vector<std::uint32_t> index(6 + inner.size(), std::numeric_limits<std::uint32_t>::max());
    std::vector<Point2> xy;
    std::vector<terrain::TriangleIndices> tris;
    for (const auto& t : named_tris) {
        terrain::TriangleIndices out{};
        for (unsigned k = 0; k < 3; ++k) {
            if (index[t[k]] == std::numeric_limits<std::uint32_t>::max()) {
                index[t[k]] = static_cast<std::uint32_t>(xy.size());
                const Lat p = position(t[k]);
                xy.push_back(strip_oracle::world(g, p.col, p.row));
            }
            out[k] = index[t[k]];
        }
        tris.push_back(out);
    }
    Start s;
    for (const auto& [e, m] : named_edges) {
        s.edges.push_back({std::min(index[e[0]], index[e[1]]), std::max(index[e[0]], index[e[1]])});
        s.masks.push_back(m);
    }
    const auto n = tris.size();
    s.mesh = IndexedMesh2{std::move(xy), std::move(tris), std::vector<std::uint8_t>(n, 0)};
    return s;
}

// The output vertices on the seam (T, B and everything between them, exactly
// on it or within the lattice slack: a seam-pass point off a general seam's
// line by rounding is still its vertex), sorted from T to B, as (world point,
// z) pairs: what K4 compares between two pieces, bit for bit.
inline std::vector<std::pair<Point2, double>> seam_sequence(const RasterGeometry& g, const Mesh& out, Lat t,
                                                            Lat b) {
    const double slack = strip_oracle::on_edge_slack(g);
    std::vector<std::pair<double, std::pair<Point2, double>>> on;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const Lat p = strip_oracle::lat(g, out.vertices[i]);
        if (p == t || p == b || on_open(t, b, p) || near_open(t, b, p, slack))
            on.push_back({strip_oracle::param_dist(t, b, p).first, {out.vertices[i], out.z[i]}});
    }
    std::sort(on.begin(), on.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
    std::vector<std::pair<Point2, double>> seq;
    for (const auto& e : on) seq.push_back(e.second);
    return seq;
}

// ---------------------------------------------------------------------------
// The seam pass's polyline
// ---------------------------------------------------------------------------

// A polyline along one segment in the lattice: its vertices at increasing
// parameter along a -> b, each with its z. Linear between neighbours.
struct Polyline {
    std::vector<double> t;
    std::vector<double> z;
    std::vector<bool> valid;

    // The value at parameter s, and whether both ends of its piece are valid.
    [[nodiscard]] std::optional<double> at(double s) const {
        auto hi = std::upper_bound(t.begin(), t.end(), s);
        if (hi == t.begin()) hi = std::next(hi);
        if (hi == t.end()) hi = std::prev(hi);
        const auto j = static_cast<std::size_t>(hi - t.begin()), i = j - 1;
        if (s == t[i]) return valid[i] ? std::optional<double>{z[i]} : std::nullopt;
        if (s == t[j]) return valid[j] ? std::optional<double>{z[j]} : std::nullopt;
        if (!valid[i] || !valid[j]) return std::nullopt;
        return z[i] + (s - t[i]) / (t[j] - t[i]) * (z[j] - z[i]);
    }
};

struct SeamFindings {
    std::size_t checked = 0;
    std::size_t on_void = 0;  // a check point on a piece with an invalid end
    std::size_t over = 0;     // |z - polyline| over tolerance + 1e-9 relative
    double worst = 0.0;
};

// Every oracle check point (strip_oracle::ruled_points over the one edge)
// against the polyline: the ruled points are generated independently of the
// producer's generator, from the 15f rulings the seam pass reuses.
inline SeamFindings seam_findings(Lat a, Lat b, const Polyline& line,
                                  const std::vector<strip_oracle::OraclePoint>& pts, double tol) {
    SeamFindings f;
    for (const auto& p : pts) {
        ++f.checked;
        const double s = strip_oracle::param_dist(a, b, p.at).first;
        const auto v = line.at(s);
        if (!v) {
            ++f.on_void;
            continue;
        }
        const double err = std::abs(p.z - *v);
        f.worst = std::max(f.worst, err);
        if (err > tol + 1e-9 * std::max(1.0, std::abs(p.z))) ++f.over;
    }
    return f;
}

}  // namespace frozen_oracle
