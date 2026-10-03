#pragma once

// Test-only fixtures and oracles for increment 15f-2
// (docs/increments/15f-edge-strip.md, D4, "The guarantee, and its wording" E1
// to E8, and "Tests for @tester" ES1 to ES11): the edge strip in the
// refinement loop. Shared by property/prop_refinement_edge_strip.cpp (the red
// suite, which calls refine_strip and refine_points(..., strip)) and
// unit/test_refinement_strip_oracle.cpp (which shows, on refine()'s output,
// that every oracle here can fail).
//
// Independence. Nothing here includes constraint_points.hpp,
// refine_points.hpp or scan.hpp. The strip oracle generates its own crossings
// and midpoints from the START mesh's constraint edges (in the world, mapped to
// the lattice by the affine map written here), takes z from the raster's own
// values by bilinear interpolation written here, finds the OUTPUT constraint
// edge holding each point by its own projection, and interpolates linearly
// along it. It borrows the producer's predicates (exact orientation, the
// lattice frame, vertex_z's NoData rule: a point is not a check point when any
// corner of the cell holding it is NoData) and never its records (no store
// order, no s, no sub-edge map).
//
// Frame. Lattice (col, row), world x = x_min + col dx, y = y_max - row dy, so
// rows grow downward. exact_geometry() takes x_min = y_max = 0, dx = 2, dy = 1:
// every lattice position maps to the world and back bit for bit, so the
// oracles recover the producer's lattice positions exactly from the output's
// world points (inserted vertices are output at x_min + col dx, y_max - row dy,
// D4 and refine.hpp). Non-square cells, so a row/col swap cannot pass. The
// producer's Delaunay frame is (col dx, -(row dy)).
//
// Nothing here asserts: the oracles RETURN their findings, so each can be
// planted and shown to fail (unit/test_refinement_strip_oracle.cpp).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace strip_oracle {

using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::raster::Raster;
using terrain::raster::RasterGeometry;
using Edges = std::vector<std::array<std::uint32_t, 2>>;

inline constexpr std::size_t kN = 33;      // nodes a side in the property fixtures
inline constexpr double kOnEdge = 1e-9;    // cells: "distance under 1e-9 cells" (ES2)

inline RasterGeometry exact_geometry(std::size_t cols = kN, std::size_t rows = kN) {
    return RasterGeometry{0.0, 0.0, 2.0, 1.0, cols, rows};
}

struct Lat {
    double col;
    double row;
    friend bool operator==(const Lat&, const Lat&) = default;
    friend auto operator<=>(const Lat&, const Lat&) = default;
};

inline Lat lat(const RasterGeometry& g, Point2 p) {
    return Lat{(p.x - g.x_min()) / g.delta_x(), (g.y_max() - p.y) / g.delta_y()};
}
inline Point2 world(const RasterGeometry& g, double col, double row) {
    return Point2{g.x_min() + col * g.delta_x(), g.y_max() - row * g.delta_y()};
}

// ---------------------------------------------------------------------------
// Heights, written from Q14's words and vertex_z's documented rule
// ---------------------------------------------------------------------------

inline bool nodata(const Raster<float>& dem, std::size_t r, std::size_t c) { return dem.is_nodata({r, c}); }
inline double at(const Raster<float>& dem, std::size_t r, std::size_t c) {
    return static_cast<double>(dem.value_at({r, c}));
}

// A node's own value at a node; otherwise bilinear over the cell holding the
// point (floor, the last cell on the far sides), nullopt when any of that
// cell's four corners is NoData.
inline std::optional<double> height(const Raster<float>& dem, Lat p) {
    const RasterGeometry& g = dem.geometry();
    if (p.col == std::floor(p.col) && p.row == std::floor(p.row)) {
        const auto r = static_cast<std::size_t>(p.row), c = static_cast<std::size_t>(p.col);
        if (nodata(dem, r, c)) return std::nullopt;
        return at(dem, r, c);
    }
    const std::size_t r0 = std::min(static_cast<std::size_t>(p.row), g.rows() - 2);
    const std::size_t c0 = std::min(static_cast<std::size_t>(p.col), g.cols() - 2);
    for (const auto& [r, c] : {std::pair{r0, c0}, {r0, c0 + 1}, {r0 + 1, c0}, {r0 + 1, c0 + 1}})
        if (nodata(dem, r, c)) return std::nullopt;
    const double ty = p.row - static_cast<double>(r0), tx = p.col - static_cast<double>(c0);
    return at(dem, r0, c0) * (1 - tx) * (1 - ty) + at(dem, r0, c0 + 1) * tx * (1 - ty)
         + at(dem, r0 + 1, c0) * (1 - tx) * ty + at(dem, r0 + 1, c0 + 1) * tx * ty;
}

// ---------------------------------------------------------------------------
// The ruled check points, generated independently
// ---------------------------------------------------------------------------

struct OraclePoint {
    Lat at;
    double z;
    std::size_t edge;  // index into the start's edge list
};

// Q14 (a) as ruled, extended, with the default that an edge's ends count as
// neighbours: every crossing of the open edge with a column or row line, and
// the midpoint of each two neighbours of (end, crossings..., end). Nothing
// recursive: Q1 (the every-point form) is not ruled. Crossings closer than
// 1e-12 in the parameter are one (a node crossing seen from both families);
// the producer snaps those to the node exactly, so the oracle's copy is within
// rounding of it. A point whose cell touches NoData is not a check point.
inline std::vector<OraclePoint> ruled_points(const Raster<float>& dem, std::span<const Point2> vertices,
                                             const Edges& edges) {
    const RasterGeometry& g = dem.geometry();
    std::vector<OraclePoint> out;
    for (std::size_t k = 0; k < edges.size(); ++k) {
        const std::uint32_t i0 = std::min(edges[k][0], edges[k][1]), i1 = std::max(edges[k][0], edges[k][1]);
        const Lat a = lat(g, vertices[i0]), b = lat(g, vertices[i1]);
        std::vector<std::pair<double, Lat>> list;
        for (double c = std::floor(std::min(a.col, b.col)) + 1; c < std::max(a.col, b.col); ++c) {
            const double t = (c - a.col) / (b.col - a.col);
            list.push_back({t, Lat{c, a.row + t * (b.row - a.row)}});
        }
        for (double r = std::floor(std::min(a.row, b.row)) + 1; r < std::max(a.row, b.row); ++r) {
            const double t = (r - a.row) / (b.row - a.row);
            list.push_back({t, Lat{a.col + t * (b.col - a.col), r}});
        }
        std::sort(list.begin(), list.end());
        std::vector<std::pair<double, Lat>> chain{{0.0, a}};
        for (const auto& c : list)
            if (c.first - chain.back().first > 1e-12 && c.first < 1.0) chain.push_back(c);
        chain.push_back({1.0, b});
        auto keep = [&](Lat p) {
            if (const auto z = height(dem, p)) out.push_back({p, *z, k});
        };
        for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
            if (i > 0) keep(chain[i].second);
            const Lat l = chain[i].second, r = chain[i + 1].second;
            keep(Lat{(l.col + r.col) / 2, (l.row + r.row) / 2});
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// The outcome, as the oracles read it
// ---------------------------------------------------------------------------

// The fields every refinement outcome has; templated callers pass
// RefineOutcome or PointRefineOutcome.
struct Mesh {
    std::vector<Point2> vertices;
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    std::vector<terrain::TriangleIndices> triangles;
    Edges edges;
    std::vector<std::uint32_t> masks;
};

template <class Outcome>
Mesh mesh_of(const Outcome& o) {
    return Mesh{o.vertices, o.z, o.valid, o.triangles, o.edges, o.masks};
}

// ---------------------------------------------------------------------------
// E1: every ruled point within tolerance of the output constraint line
// ---------------------------------------------------------------------------

struct StripFindings {
    std::size_t checked = 0;
    std::size_t unfiled = 0;   // on no output constraint edge
    std::size_t on_void = 0;   // on an output constraint edge with an invalid end
    std::size_t over = 0;      // |z - mesh value| over tolerance + slack
    double worst = 0.0;        // the largest error over filed, non-void points
    std::vector<Lat> over_at;  // where each `over` point is
};

inline std::pair<double, double> param_dist(Lat a, Lat b, Lat p) {
    const double ux = b.col - a.col, uy = b.row - a.row, len2 = ux * ux + uy * uy;
    const double t = ((p.col - a.col) * ux + (p.row - a.row) * uy) / len2;
    const double d = std::abs((p.col - a.col) * uy - (p.row - a.row) * ux) / std::sqrt(len2);
    return {t, d};
}

inline StripFindings strip_findings(const RasterGeometry& g, const std::vector<OraclePoint>& pts,
                                    const Mesh& out, double tol) {
    std::vector<Lat> lv;
    for (const Point2 v : out.vertices) lv.push_back(lat(g, v));
    StripFindings f;
    for (const OraclePoint& p : pts) {
        ++f.checked;
        bool filed = false, is_void = false;
        double err = std::numeric_limits<double>::infinity();
        for (const auto& e : out.edges) {
            const Lat a = lv[e[0]], b = lv[e[1]];
            const auto [t, d] = param_dist(a, b, p.at);
            if (d > kOnEdge || t < -1e-12 || t > 1.0 + 1e-12) continue;
            filed = true;
            if (!(out.valid[e[0]] && out.valid[e[1]])) {
                is_void = true;
                continue;
            }
            const double s = std::clamp(t, 0.0, 1.0);
            err = std::min(err, std::abs(p.z - (out.z[e[0]] + s * (out.z[e[1]] - out.z[e[0]]))));
        }
        if (!filed) {
            ++f.unfiled;
            continue;
        }
        if (std::isinf(err)) {
            f.on_void += is_void ? 1 : 0;
            continue;
        }
        f.worst = std::max(f.worst, err);
        if (err > tol + 1e-9 * std::max(1.0, std::abs(p.z))) {
            ++f.over;
            f.over_at.push_back(p.at);
        }
    }
    return f;
}

// ---------------------------------------------------------------------------
// The DEM-node guarantee (increment 14, E2) and the triangle orientation
// ---------------------------------------------------------------------------

struct NodeFindings {
    std::size_t not_ccw = 0;  // output triangles not strictly counter-clockwise
    std::size_t over = 0;     // (node, triangle) pairs over tolerance + slack
    double worst = 0.0;
};

inline double cross(Point2 a, Point2 b, Point2 c) { return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x); }

// Every valid DEM node in every CLOSED output triangle with three valid
// vertices, against that triangle's plane recomputed from the output's z.
// Membership is the exact orientation on (col, -row); the plane is the
// barycentric weights in the same frame.
inline NodeFindings node_findings(const Raster<float>& dem, const Mesh& out, double tol) {
    using terrain::pred::DefaultKernel;
    using terrain::pred::Orientation;
    const RasterGeometry& g = dem.geometry();
    std::vector<Point2> fp;
    double zmax = 0.0;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const Lat l = lat(g, out.vertices[i]);
        fp.push_back(Point2{l.col, -l.row});
        if (out.valid[i]) zmax = std::max(zmax, std::abs(out.z[i]));
    }
    NodeFindings f;
    for (const auto& tri : out.triangles) {
        const Point2 a = fp[tri[0]], b = fp[tri[1]], c = fp[tri[2]];
        if (DefaultKernel::orient2d(a, b, c) != Orientation::CounterClockwise) {
            ++f.not_ccw;
            continue;
        }
        if (!(out.valid[tri[0]] && out.valid[tri[1]] && out.valid[tri[2]])) continue;
        const double two_a = cross(a, b, c);
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
                if (p == a || p == b || p == c || nodata(dem, ur, uc)) continue;
                const double plane = (cross(p, b, c) * out.z[tri[0]] + cross(a, p, c) * out.z[tri[1]]
                                      + cross(a, b, p) * out.z[tri[2]]) / two_a;
                const double err = std::abs(at(dem, ur, uc) - plane);
                f.worst = std::max(f.worst, err);
                if (err > tol + 1e-9 * std::max(1.0, zmax)) ++f.over;
            }
    }
    return f;
}

// ---------------------------------------------------------------------------
// Constrained Delaunay, in the producer's frame (col dx, -(row dy))
// ---------------------------------------------------------------------------

// The number of (interior non-constraint edge, side) pairs whose apex lies
// strictly inside the other triangle's circumcircle, by the exact incircle.
inline std::size_t delaunay_violations(const RasterGeometry& g, const Mesh& out) {
    using terrain::pred::DefaultKernel;
    using terrain::pred::Incircle;
    auto lf = [&](std::uint32_t i) {
        const Lat l = lat(g, out.vertices[i]);
        return Point2{l.col * g.delta_x(), -(l.row * g.delta_y())};
    };
    std::set<std::pair<std::uint32_t, std::uint32_t>> constrained;
    for (const auto& e : out.edges) constrained.insert(std::minmax(e[0], e[1]));
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> sides;
    for (std::size_t t = 0; t < out.triangles.size(); ++t)
        for (unsigned k = 0; k < 3; ++k)
            sides[std::minmax(out.triangles[t][k], out.triangles[t][(k + 1) % 3])].push_back(t);
    std::size_t bad = 0;
    for (const auto& [e, ts] : sides) {
        if (ts.size() > 2) ++bad;  // not a manifold triangulation
        if (ts.size() != 2 || constrained.contains(e)) continue;
        for (unsigned s = 0; s < 2; ++s) {
            const auto& tri = out.triangles[ts[s]];
            std::uint32_t apex = 0;
            for (const auto x : out.triangles[ts[1 - s]])
                if (x != e.first && x != e.second) apex = x;
            if (DefaultKernel::incircle(lf(tri[0]), lf(tri[1]), lf(tri[2]), lf(apex)) == Incircle::Inside) ++bad;
        }
    }
    return bad;
}

// ---------------------------------------------------------------------------
// E5, E6, E8: what the run may add, and where the constraints go
// ---------------------------------------------------------------------------

struct ShapeFindings {
    std::size_t moved_start = 0;     // a start vertex not output as given
    std::size_t stray = 0;           // a new vertex neither a DEM node nor strictly inside a start constraint edge
    std::size_t wrong_z = 0;         // a new vertex whose z is not the DEM's there (1e-9 relative)
    std::size_t off_node_new = 0;    // new vertices that are not DEM nodes (strip insertions)
    std::size_t coincident = 0;      // pairs of output vertices at one position
    std::size_t unplaced_edge = 0;   // an output constraint edge on no start constraint edge
    std::size_t wrong_mask = 0;      // ... or carrying another mask
    std::size_t broken_chain = 0;    // a start constraint edge whose pieces do not chain 0 -> 1
};

inline bool is_node(Lat p) { return p.col == std::floor(p.col) && p.row == std::floor(p.row); }

// `start_vertices`, `start_edges`, `start_masks` describe the mesh the run
// started from. A new vertex at a node is a DEM node (inserted by the rescan,
// or a strip point at a node crossing); any other new vertex must lie strictly
// inside a start constraint edge (a strip point, E5: never inside a triangle).
// Either way its z is the DEM's there (E6: a strip vertex carries its point's
// own z, which is vertex_z).
inline ShapeFindings shape_findings(const Raster<float>& dem, std::span<const Point2> start_vertices,
                                    const Edges& start_edges, const std::vector<std::uint32_t>& start_masks,
                                    const Mesh& out) {
    const RasterGeometry& g = dem.geometry();
    ShapeFindings f;
    const std::size_t n0 = start_vertices.size();
    std::vector<std::array<Lat, 2>> segs;
    for (const auto& e : start_edges) segs.push_back({lat(g, start_vertices[e[0]]), lat(g, start_vertices[e[1]])});
    std::set<Lat> seen;
    for (std::size_t i = 0; i < out.vertices.size(); ++i) {
        const Lat p = lat(g, out.vertices[i]);
        if (!seen.insert(p).second) ++f.coincident;
        if (i < n0) {
            if (!(out.vertices[i] == start_vertices[i])) ++f.moved_start;
            continue;
        }
        if (!is_node(p)) {
            ++f.off_node_new;
            bool on = false;
            for (const auto& s : segs) {
                const auto [t, d] = param_dist(s[0], s[1], p);
                on = on || (d <= kOnEdge && t > 0.0 && t < 1.0);
            }
            if (!on) ++f.stray;
        }
        const auto z = height(dem, p);
        if (!z || !out.valid[i] || std::abs(out.z[i] - *z) > 1e-9 * std::max(1.0, std::abs(*z))) ++f.wrong_z;
    }
    std::vector<std::vector<std::pair<double, double>>> pieces(segs.size());
    for (std::size_t k = 0; k < out.edges.size(); ++k) {
        const Lat p = lat(g, out.vertices[out.edges[k][0]]), q = lat(g, out.vertices[out.edges[k][1]]);
        std::optional<std::size_t> s;
        for (std::size_t i = 0; i < segs.size() && !s; ++i) {
            const auto [tp, dp] = param_dist(segs[i][0], segs[i][1], p);
            const auto [tq, dq] = param_dist(segs[i][0], segs[i][1], q);
            if (dp <= kOnEdge && dq <= kOnEdge && std::min(tp, tq) >= -kOnEdge && std::max(tp, tq) <= 1 + kOnEdge)
                s = i;
        }
        if (!s) {
            ++f.unplaced_edge;
            continue;
        }
        if (out.masks[k] != start_masks[*s]) ++f.wrong_mask;
        pieces[*s].push_back(std::minmax(param_dist(segs[*s][0], segs[*s][1], p).first,
                                         param_dist(segs[*s][0], segs[*s][1], q).first));
    }
    for (auto& ps : pieces) {
        std::sort(ps.begin(), ps.end());
        bool ok = !ps.empty() && std::abs(ps.front().first) <= kOnEdge && std::abs(ps.back().second - 1.0) <= kOnEdge;
        for (std::size_t j = 1; ok && j < ps.size(); ++j) ok = std::abs(ps[j].first - ps[j - 1].second) <= kOnEdge;
        if (!ok) ++f.broken_chain;
    }
    return f;
}

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

// A start mesh in world coordinates with its constraint edges.
struct Start {
    IndexedMesh2 mesh;
    Edges edges;
    std::vector<std::uint32_t> masks;
};

// A jittered 5 x 5 vertex grid over a 33 x 33 node DEM, two triangles per
// quad (diagonal top-left to bottom-right), nominal lattice positions
// 2 + 7 i, so the outline is inside the node rectangle and off it.
// Constrained chains, each with its own mask:
//   1: the outline, jittered off-grid except its two diagonal corners;
//   2: the middle row of vertices, row jitter 0: a chain ALONG grid line 16;
//   4: the second column of vertices, column jitter 0: a chain along column 9;
//   8: the main diagonal, no jitter: nodes (2,2) ... (30,30), so each of its
//      edges passes THROUGH six DEM nodes.
// Jitter is a multiple of 1/64 cell in [-1.5, 1.5], so every position is
// dyadic and exact in exact_geometry().
inline Start jittered(const RasterGeometry& g, std::uint32_t seed) {
    std::mt19937 gen{seed};
    constexpr std::uint32_t n = 5;
    auto jit = [&] { return (static_cast<double>(gen() % 193u) - 96.0) / 64.0; };
    std::vector<Point2> xy;
    for (std::uint32_t i = 0; i < n; ++i)
        for (std::uint32_t j = 0; j < n; ++j) {
            double dr = jit(), dc = jit();
            if (i == 2) dr = 0.0;
            if (j == 1) dc = 0.0;
            if (i == j) dr = dc = 0.0;
            xy.push_back(world(g, 2.0 + 7.0 * j + dc, 2.0 + 7.0 * i + dr));
        }
    auto at = [](std::uint32_t i, std::uint32_t j) { return i * n + j; };
    std::vector<terrain::TriangleIndices> tris;
    for (std::uint32_t i = 0; i + 1 < n; ++i)
        for (std::uint32_t j = 0; j + 1 < n; ++j) {
            const auto tl = at(i, j), tr = at(i, j + 1), bl = at(i + 1, j), br = at(i + 1, j + 1);
            tris.push_back({tl, bl, br});
            tris.push_back({tl, br, tr});
        }
    Start s;
    s.mesh = IndexedMesh2{std::move(xy), std::move(tris), std::vector<std::uint8_t>(2 * (n - 1) * (n - 1), 0)};
    auto add = [&](std::uint32_t a, std::uint32_t b, std::uint32_t mask) {
        s.edges.push_back({std::min(a, b), std::max(a, b)});
        s.masks.push_back(mask);
    };
    for (std::uint32_t j = 0; j + 1 < n; ++j) add(at(0, j), at(0, j + 1), 1);
    for (std::uint32_t i = 0; i + 1 < n; ++i) add(at(i, n - 1), at(i + 1, n - 1), 1);
    for (std::uint32_t j = 0; j + 1 < n; ++j) add(at(n - 1, j), at(n - 1, j + 1), 1);
    for (std::uint32_t i = 0; i + 1 < n; ++i) add(at(i, 0), at(i + 1, 0), 1);
    for (std::uint32_t j = 0; j + 1 < n; ++j) add(at(2, j), at(2, j + 1), 2);
    for (std::uint32_t i = 0; i + 1 < n; ++i) add(at(i, 1), at(i + 1, 1), 4);
    for (std::uint32_t i = 0; i + 1 < n; ++i) add(at(i, i), at(i + 1, i + 1), 8);
    return s;
}

enum class Terrain { Rough, Smooth, SmoothWithHole, Plane };

// Seeded DEMs over g. mt19937's raw output is specified by the standard; the
// distributions are not, so values are scaled by hand. Every value is a float
// of few bits. SmoothWithHole adds a 3 x 3 NoData patch whose position depends
// on the seed, somewhere the outline or a chain may cross.
inline Raster<float> terrain_dem(const RasterGeometry& g, Terrain kind, std::uint32_t seed) {
    std::mt19937 gen{seed * 7919u + 17u};
    const std::size_t rows = g.rows(), cols = g.cols();
    std::vector<float> v(rows * cols);
    const double p1 = static_cast<double>(gen() % 628u) / 100.0, p2 = static_cast<double>(gen() % 628u) / 100.0;
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < cols; ++c) {
            const double x = static_cast<double>(c), y = static_cast<double>(r);
            double z = 0.0;
            switch (kind) {
            case Terrain::Rough: z = static_cast<double>(gen() % 4096u) / 64.0; break;
            case Terrain::Smooth:
            case Terrain::SmoothWithHole:
                z = std::round((50.0 + 20.0 * std::sin(0.7 * x + p1) + 15.0 * std::cos(0.55 * y + p2)) * 64.0) / 64.0;
                break;
            case Terrain::Plane: z = 100.0 + 0.5 * x + 0.25 * y; break;
            }
            v[r * cols + c] = static_cast<float>(z);
        }
    if (kind == Terrain::SmoothWithHole) {
        const std::size_t r0 = 1 + gen() % (rows - 4), c0 = 1 + gen() % (cols - 4);
        for (std::size_t r = r0; r < r0 + 3; ++r)
            for (std::size_t c = c0; c < c0 + 3; ++c) v[r * cols + c] = -9999.0f;
        return Raster<float>{g, std::move(v), -9999.0f};
    }
    return Raster<float>{g, std::move(v)};
}

// ES1: Surprise 3 in miniature. A 9 x 9 node DEM, flat 0 except nodes
// (row 4, col 3) and (row 4, col 4), which are 10. A quadrilateral domain
// whose bottom side runs along row 3.5 from col 0.5 to col 7.5, between grid
// rows 3 and 4: every DEM node inside the domain is 0 and so is every vertex
// z, so refine inserts nothing, yet the bottom side crosses columns 3 and 4
// at z (0 + 10) / 2 = 5. Vertex order v0 (7.5, 3.5), v1 (7.25, 0.5),
// v2 (0.75, 0.5), v3 (0.5, 3.5): counter-clockwise in world, fanned from v0,
// so the bottom side v3 -> v0 runs from the HIGHER index to the lower in the
// only triangle that holds it (ES5). The top corners are pulled in so the
// four are not cocircular.
inline Raster<float> bump_dem() {
    const auto g = exact_geometry(9, 9);
    std::vector<float> z(81, 0.0f);
    z[4 * 9 + 3] = 10.0f;
    z[4 * 9 + 4] = 10.0f;
    return Raster<float>{g, std::move(z)};
}

inline std::vector<std::array<double, 2>> bump_ring() { return {{7.5, 3.5}, {7.25, 0.5}, {0.75, 0.5}, {0.5, 3.5}}; }

// The ring above, fanned from vertex 0; every side constrained, the bottom
// side with mask 2 and the others with mask 1. With `below`, a second
// quadrilateral under the bottom side (to row 6.5) makes that side interior.
inline Start bump_start(bool below = false) {
    const auto g = exact_geometry(9, 9);
    std::vector<Point2> xy;
    for (const auto& [c, r] : bump_ring()) xy.push_back(world(g, c, r));
    std::vector<terrain::TriangleIndices> tris{{0, 1, 2}, {0, 2, 3}};
    Start s;
    s.edges = {{0, 1}, {1, 2}, {2, 3}, {0, 3}};
    s.masks = {1, 1, 1, 2};
    if (below) {
        xy.push_back(world(g, 0.75, 6.5));  // 4
        xy.push_back(world(g, 7.25, 6.5));  // 5
        tris.push_back({3, 4, 5});
        tris.push_back({3, 5, 0});
        s.edges = {{0, 1}, {1, 2}, {2, 3}, {0, 3}, {3, 4}, {4, 5}, {0, 5}};
        s.masks = {1, 1, 1, 2, 1, 1, 1};
    }
    s.mesh = IndexedMesh2{std::move(xy), std::move(tris), std::vector<std::uint8_t>(below ? 4 : 2, 0)};
    return s;
}

// The z and valid a run starts from: the oracle's own heights at the start
// vertices (vertex_z's rule), 0 and invalid where it refuses.
inline std::pair<std::vector<double>, std::vector<std::uint8_t>> start_z(const Raster<float>& dem,
                                                                         std::span<const Point2> vertices) {
    std::vector<double> z;
    std::vector<std::uint8_t> valid;
    for (const Point2 v : vertices) {
        const auto h = height(dem, lat(dem.geometry(), v));
        z.push_back(h.value_or(0.0));
        valid.push_back(h ? 1 : 0);
    }
    return {z, valid};
}

}  // namespace strip_oracle
