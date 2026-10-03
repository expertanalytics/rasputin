#pragma once

// The final check (docs/increments/15c-geographic-dem.md, D5): greedy insertion
// against stored check points, starting from phase 1's mesh. refine.hpp's loop
// with a different scanner; refine.hpp itself is not edited (J1).
//
// The store's geometry() is the frame and the buckets, nothing more: phase 2
// reads no raster value. A start vertex keeps the z it is given (NaN where it
// is not valid); an inserted vertex is a check point at its stored position,
// output at (x_min + col dx, y_max - row dy), with the point's own z.
//
// Scan (parallel, read-only, one result per triangle). For each cell row the
// triangle meets, the column range it covers in that band, widened by a cell
// each side, bounds the stored points tested; membership in the CLOSED
// triangle is three exact orient_sign calls. A point equal to a corner is
// skipped. Error is |z - plane| with scan.hpp's off-node double expression
// (the largest corner difference for a sliver whose 2A rounds to <= 0); the
// worst point wins by strictly larger error, so ties go to the first in store
// order. A void triangle (a NaN corner) takes the point nearest a void corner,
// and counts its points in `uncovered`.
//
// Split (serial, triangle-index order), as refine: split_inside, or split_edge
// for a point on an edge, skipped this round when the neighbour across it was
// already touched; legalise_around follows. The start is legalised once first,
// as refine does. Terminates: each insertion is a stored point not yet a
// vertex, and the store is finite.
//
// Coincident points. A stored point equal to a start vertex is a corner of
// every triangle holding it, so it is never scanned nor inserted. After the
// loop each start vertex looks it up in its own cell; matches are counted in
// `coincident`, their largest |z - vertex z| in `coincident_max_error`.
//
// The edge strip (docs/increments/15f-edge-strip.md, D4 and L1-L9). One loop,
// detail::point_loop, scans up to three sets per triangle: the store's points
// (refine_points), the strip's points, and the DEM's nodes (refine_strip, only
// in a slot this run has written, L3). Strip points are filed by constrained
// sub-edge, not found by membership (F1): a triangle scans the strip sub-edges
// it owns (lower vertex index to higher, or no triangle across), and the mesh
// value there is linear in s between the sub-edge's ends. Across sets the
// first void result wins, else the strictly largest error, ties to the earlier
// set (source, strip, DEM). A strip point goes in by split_edge guarded by
// strip_fits (foot_fits, and the split's new edges locally Delaunay); refused,
// it is marked and never named again, and its triangle stays active (L1).
// Every split of a strip sub-edge replaces its record by two (step 4); a point
// of another set lands on the strip point at its exact position, if any, which
// is then consumed (L5). After the loop every strip point is measured against
// its final sub-edge (step 6).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/constraint_points.hpp>
#include <terrain/refinement/refine.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace terrain::refinement {

struct PointRefineOptions {
    double tolerance = 0.0;  // metres, finite and >= 0
    unsigned threads = 0;    // 0: hardware concurrency
};

// RefineOutcome, with max_error over check points, plus the coincident points
// and the edge strip's figures (15f-edge-strip.md, D4 "Outcome", L4, L6).
struct PointRefineOutcome : RefineOutcome {
    std::size_t coincident = 0;         // stored points equal to a start vertex
    double coincident_max_error = 0.0;  // their largest |z - vertex z|, valid vertices only
    std::size_t strip_points = 0;       // strip.size(), 0 without a strip
    std::size_t strip_inserted = 0;     // a subset of `inserted`, carved ones included
    double strip_max_error = 0.0;       // at the end, refused points excluded
    std::size_t strip_refused = 0;
    double strip_refused_max_error = 0.0;
    std::size_t nodes_inserted = 0;  // refine_strip: DEM nodes, a subset of `inserted`
};

namespace detail {

enum class PointSet : std::uint8_t { Source, Strip, Dem };

// One triangle's scan: the winner across sets (for the split), and the
// non-strip set's own largest error and every set's uncovered count (L4).
struct PointScan {
    double max_error = 0.0;                  // the non-strip set's own; 0 when void
    std::optional<mesh::MeshVertex> point;  // the winner, or a void carve point
    double z = 0.0;                          // its z
    double error = 0.0;                      // its error
    NodeLocation where = NodeLocation::Inside;
    bool is_void = false;       // the winner is a carve point (source: the triangle is void)
    std::size_t uncovered = 0;  // points of every set in void triangles or on void sub-edges
    PointSet set = PointSet::Source;
    std::size_t strip_index = 0;  // a strip winner's flag index
    double s = 0.0;               // a strip winner's parameter

    // Offered in set order: the first void candidate wins, else the strictly larger error.
    void offer(const PointScan& c) {
        if ((point && is_void) || !(c.is_void || !point || c.error > error))
            return;
        const double own = max_error;
        const std::size_t n = uncovered;
        *this = c;
        max_error = own;
        uncovered = n;
    }
};

template <class Store>
[[nodiscard]] PointScan scan_points(const Store& points, const mesh::LatticeMesh& m,
                                    std::span<const double> zt, std::uint32_t t) {
    using mesh::MeshVertex;
    const raster::RasterGeometry& g = points.geometry();
    const auto& tri = m.triangles()[t];
    const std::array<MeshVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
    const std::array<double, 3> zv{zt[tri[0]], zt[tri[1]], zt[tri[2]]};
    PointScan r;
    r.is_void = std::isnan(zv[0]) || std::isnan(zv[1]) || std::isnan(zv[2]);
    const Point2 f0 = v[0].frame(), f1 = v[1].frame(), f2 = v[2].frame();
    const double two_a = (f1.x - f0.x) * (f2.y - f0.y) - (f1.y - f0.y) * (f2.x - f0.x);
    auto value = [&](unsigned k, MeshVertex p) {
        const Point2 a = v[k].frame(), b = v[(k + 1) % 3].frame(), q = p.frame();
        return (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x);
    };
    double nearest = std::numeric_limits<double>::infinity();
    auto visit = [&](MeshVertex p, float zf) {
        if (p == v[0] || p == v[1] || p == v[2])
            return;
        for (unsigned k = 0; k < 3; ++k)
            if (mesh::orient_sign(v[k], v[(k + 1) % 3], p) < 0)
                return;
        const auto z = static_cast<double>(zf);
        if (r.is_void) {
            ++r.uncovered;
            double d = std::numeric_limits<double>::infinity();
            for (unsigned k = 0; k < 3; ++k)
                if (std::isnan(zv[k]))
                    d = std::min(d, (p.row - v[k].row) * (p.row - v[k].row) + (p.col - v[k].col) * (p.col - v[k].col));
            if (d < nearest) {
                nearest = d;
                r.point = p;
                r.z = z;
            }
            return;
        }
        const double err =
            two_a > 0.0 ? std::abs(z - (value(1, p) * zv[0] + value(2, p) * zv[1] + value(0, p) * zv[2]) / two_a)
                        : std::max({std::abs(z - zv[0]), std::abs(z - zv[1]), std::abs(z - zv[2])});
        if (err > r.max_error) {
            r.max_error = err;
            r.point = p;
            r.z = z;
        }
    };

    // Per cell row, the triangle's column range in the band [b, b + 1]: its
    // vertices in the band and its edges' crossings of the band's two lines.
    const std::size_t last_row = last_cell(g.rows()), last_col = last_cell(g.cols());
    const auto [rlo, rhi] = std::minmax({v[0].row, v[1].row, v[2].row});
    const auto b_hi = std::min(static_cast<std::size_t>(rhi), last_row);
    for (std::size_t b = std::min(static_cast<std::size_t>(rlo), last_row); b <= b_hi; ++b) {
        const auto top = static_cast<double>(b), bottom = top + 1.0;
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (unsigned k = 0; k < 3; ++k) {
            const MeshVertex a = v[k], c = v[(k + 1) % 3];
            if (a.row >= top && a.row <= bottom) {
                lo = std::min(lo, a.col);
                hi = std::max(hi, a.col);
            }
            for (const double y : {top, bottom})
                if ((a.row - y) * (c.row - y) < 0.0) {
                    const double x = a.col + (c.col - a.col) * (y - a.row) / (c.row - a.row);
                    lo = std::min(lo, x);
                    hi = std::max(hi, x);
                }
        }
        if (!(lo <= hi))
            continue;
        const auto c0 = static_cast<std::size_t>(std::max(std::floor(lo) - 1.0, 0.0));
        const auto c1 = std::min(static_cast<std::size_t>(std::floor(hi) + 1.0), last_col);
        points.for_each_in(b, c0, c1, visit);
    }
    if (r.point)
        r.where = mesh::orient_sign(v[0], v[1], *r.point) == 0   ? NodeLocation::Edge0
                  : mesh::orient_sign(v[1], v[2], *r.point) == 0 ? NodeLocation::Edge1
                  : mesh::orient_sign(v[2], v[0], *r.point) == 0 ? NodeLocation::Edge2
                                                                 : NodeLocation::Inside;
    return r;
}

// A constrained sub-edge of strip edge k, ends a and b in the order of s
// (s_a < s_b, a nearer P0), keyed by its unordered vertex pair (D4 step 1, L5).
struct SubEdge {
    std::size_t k;
    std::uint32_t a, b;
    double s_a, s_b;
};
using SubEdges = std::unordered_map<std::uint64_t, SubEdge>;

[[nodiscard]] inline std::uint64_t edge_key(std::uint32_t u, std::uint32_t v) {
    return std::uint64_t{std::min(u, v)} << 32 | std::max(u, v);
}

// The mesh value along a sub-edge at s: z_a + sigma (z_b - z_a), and an end's
// own z at its own s (step 6, L5). The scan and step 6 share it, so they agree.
[[nodiscard]] inline double along(const SubEdge& e, std::span<const double> zt, double s) {
    if (s == e.s_a || s == e.s_b)
        return zt[s == e.s_a ? e.a : e.b];
    return zt[e.a] + (s - e.s_a) / (e.s_b - e.s_a) * (zt[e.b] - zt[e.a]);
}

// The strip's candidate in t, offered to r: over the sub-edges t owns, the
// unrefused points strictly inside, the worst by strictly larger error (lower
// edge index, then lower s, on ties); on a void sub-edge the point nearest an
// invalid end, its points counted in `uncovered` (D4 step 2).
inline void scan_strip(const ConstraintCheckPoints& strip, std::span<const std::size_t> offset,
                       std::span<const char> refused, const SubEdges& subs, const mesh::LatticeMesh& m,
                       std::span<const double> zt, std::uint32_t t, PointScan& r) {
    const auto& tri = m.triangles()[t];
    PointScan best;
    double nearest = std::numeric_limits<double>::infinity();
    for (unsigned e = 0; e < 3; ++e) {
        const std::uint32_t u = tri[e], v = tri[(e + 1) % 3];
        if (!m.is_constrained(t, e) || (u > v && m.neighbours(t)[e] != mesh::kNoNeighbour))
            continue;
        const auto it = subs.find(edge_key(u, v));
        if (it == subs.end())
            continue;
        const SubEdge& se = it->second;
        const auto pts = strip.on_edge(se.k);
        const bool void_a = std::isnan(zt[se.a]), void_b = std::isnan(zt[se.b]);
        auto j = static_cast<std::size_t>(
            std::partition_point(pts.begin(), pts.end(), [&](const ConstraintPoint& p) { return p.s <= se.s_a; })
            - pts.begin());
        for (; j < pts.size() && pts[j].s < se.s_b; ++j) {
            const ConstraintPoint& p = pts[j];
            if (refused[offset[se.k] + j] != 0)
                continue;
            PointScan c;
            c.set = PointSet::Strip;
            c.point = p.at;
            c.z = p.z;
            c.where = static_cast<NodeLocation>(e + 1);
            c.strip_index = offset[se.k] + j;
            c.s = p.s;
            if (void_a || void_b) {
                ++r.uncovered;
                const double inf = std::numeric_limits<double>::infinity();
                const double d = std::min(void_a ? p.s - se.s_a : inf, void_b ? se.s_b - p.s : inf);
                c.is_void = true;
                if (d < nearest) {
                    nearest = d;
                    best = c;
                }
            } else if (c.error = std::abs(p.z - along(se, zt, p.s)); !best.is_void && c.error > best.error) {
                best = c;
            }
        }
    }
    if (best.point)
        r.offer(best);
}

// D4 step 4: a split at vertex q of the constrained edge (u, v) replaces its
// sub-edge record by two. s_q is a strip point's own s; for any other point,
// that of the strip point at exactly q's position (consumed, L5), else q's
// projection onto the sub-edge mapped into [s_a, s_b].
inline void cut(SubEdges& subs, const ConstraintCheckPoints& strip, const mesh::LatticeMesh& m,
                std::uint32_t u, std::uint32_t v, std::uint32_t q, std::optional<double> s) {
    const auto it = subs.find(edge_key(u, v));
    if (it == subs.end())
        return;
    const SubEdge se = it->second;
    const mesh::MeshVertex p = m.vertices()[q], a = m.vertices()[se.a], b = m.vertices()[se.b];
    for (const ConstraintPoint& c : strip.on_edge(se.k))
        if (!s && c.s > se.s_a && c.s < se.s_b && c.at == p)
            s = c.s;
    if (!s) {
        const double dc = b.col - a.col, dr = b.row - a.row;
        const double sigma =
            std::clamp(((p.col - a.col) * dc + (p.row - a.row) * dr) / (dc * dc + dr * dr), 0.0, 1.0);
        s = se.s_a + sigma * (se.s_b - se.s_a);
    }
    subs.erase(it);
    subs[edge_key(se.a, q)] = SubEdge{se.k, se.a, q, se.s_a, *s};
    subs[edge_key(q, se.b)] = SubEdge{se.k, q, se.b, *s, se.s_b};
}

// foot_fits, and the split's new edges q-c (and q-d) locally Delaunay. A point
// a hair off its sub-edge near an end can lie outside t's circle, and
// legalise_around tests only the edges opposite q.
inline bool strip_fits(const mesh::LatticeMesh& m, std::uint32_t t, unsigned e, mesh::MeshVertex f,
                       const mesh::LatticeFrame& fr) {
    using K = pred::DefaultKernel;
    if (!foot_fits(m, t, e, f))
        return false;
    // Whether edge (p, q) of the counter-clockwise (p, q, r) must flip, s across it (lawson.hpp's must_flip).
    const auto flips = [&](mesh::MeshVertex p, mesh::MeshVertex q, mesh::MeshVertex r, mesh::MeshVertex s) {
        const Point2 a = fr.at(p), b = fr.at(q), c = fr.at(r), d = fr.at(s);
        if (K::orient2d(a, b, c) == pred::Orientation::CounterClockwise)
            return K::incircle(a, b, c, d) == pred::Incircle::Inside;
        return K::orient2d(b, a, d) == pred::Orientation::CounterClockwise
            && K::incircle(b, a, d, c) == pred::Incircle::Inside;
    };
    const mesh::MeshVertex a = m.corner(t, e), b = m.corner(t, (e + 1) % 3), c = m.corner(t, (e + 2) % 3);
    if (flips(f, c, a, b))
        return false;
    const std::uint32_t u = m.neighbours(t)[e];
    if (u == mesh::kNoNeighbour)
        return true;
    unsigned k = 0;
    while (m.neighbours(u)[k] != t)
        ++k;
    return !flips(f, m.corner(u, (k + 2) % 3), b, a);
}

// L12: the constrained edge of t whose line p lies within 1e-10 lattice units
// of, (col, row) Euclidean, with p's projection strictly inside it; the
// nearest, ties to the lower edge index.
inline std::optional<unsigned> near_constraint(const mesh::LatticeMesh& m, std::uint32_t t, mesh::MeshVertex p) {
    std::optional<unsigned> best;
    double best_d = 1e-10;
    for (unsigned e = 0; e < 3; ++e) {
        const mesh::MeshVertex a = m.corner(t, e), b = m.corner(t, (e + 1) % 3);
        const double dc = b.col - a.col, dr = b.row - a.row, len2 = dc * dc + dr * dr;
        const double sigma = ((p.col - a.col) * dc + (p.row - a.row) * dr) / len2;
        const double d = std::abs(dc * (p.row - a.row) - dr * (p.col - a.col)) / std::sqrt(len2);
        if (m.is_constrained(t, e) && sigma > 0.0 && sigma < 1.0 && (d < best_d || (!best && d == best_d))) {
            best = e;
            best_d = d;
        }
    }
    return best;
}

struct NoSet {};  // point_loop without a store, or without a DEM

[[nodiscard]] inline PointRefineOutcome point_refusal(RefineOutcome r) {
    PointRefineOutcome out;
    static_cast<RefineOutcome&>(out) = std::move(r);
    return out;
}

// The loop of refine_points and refine_strip (D4). `points` (a store) or `dem`
// may be null, `strip` too; refusals (2) to (4) of L2, the tolerance being the
// caller's.
template <class Store, class R>
[[nodiscard]] PointRefineOutcome point_loop(const std::string& name, const raster::RasterGeometry& g,
                                            const Store* points, const R* dem, const ConstraintCheckPoints* strip,
                                            const IndexedMesh2& start, std::span<const double> z,
                                            std::span<const std::uint8_t> valid,
                                            std::span<const std::array<std::uint32_t, 2>> edges,
                                            std::span<const std::uint32_t> masks, const PointRefineOptions& options) {
    constexpr bool has_store = !std::is_same_v<Store, NoSet>, has_dem = !std::is_same_v<R, NoSet>;
    SubEdges subs;
    std::vector<std::size_t> offset{0};
    if (strip) {
        if (!(strip->geometry() == g))
            throw std::logic_error(name + ": the edge strip's raster geometry differs from the run's");
        std::vector<std::uint64_t> given;
        for (const auto& e : edges)
            given.push_back(edge_key(e[0], e[1]));
        std::sort(given.begin(), given.end());
        for (std::size_t k = 0; k < strip->edge_count(); ++k) {
            const auto [p0, p1] = strip->edge(k);
            if (!std::binary_search(given.begin(), given.end(), edge_key(p0, p1)))
                throw std::logic_error(name + ": strip edge (" + std::to_string(p0) + ", " + std::to_string(p1)
                                       + ") is not a constraint edge of the start mesh");
            subs[edge_key(p0, p1)] = SubEdge{k, p0, p1, 0.0, 1.0};
            offset.push_back(offset.back() + strip->on_edge(k).size());
        }
    }
    auto built = to_lattice(g, start, edges, masks);
    if (auto* refused = std::get_if<RefineOutcome>(&built))
        return point_refusal(std::move(*refused));
    auto& m = std::get<mesh::LatticeMesh>(built);

    PointRefineOutcome out;
    out.strip_points = strip ? strip->size() : 0;
    const std::size_t n0 = start.vertices().size();
    std::vector<double> zt(n0, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < n0 && i < z.size() && i < valid.size(); ++i)
        if (valid[i] != 0)
            zt[i] = z[i];

    const auto frame = mesh::lattice_frame(g.delta_x(), g.delta_y(), g.rows(), g.cols());
    std::vector<char> written(m.triangle_count(), 0), refused(offset.back(), 0);  // L3, L1
    out.flips = mesh::legalise_all<pred::DefaultKernel>(m, frame, [&](std::uint32_t s) { written[s] = 1; });
    std::vector<PointScan> results;
    mesh::FlipStack flip_stack;
    std::vector<std::uint32_t> active(m.triangle_count());
    for (std::uint32_t t = 0; t < active.size(); ++t)
        active[t] = t;

    // Source, strip, DEM, in that order (D4, "Combining").
    const auto scan_one = [&](std::uint32_t t) {
        PointScan r;
        if constexpr (has_store) {
            r = scan_points(*points, m, zt, t);
            r.error = r.max_error;
        }
        if (strip)
            scan_strip(*strip, offset, refused, subs, m, zt, t, r);
        if constexpr (has_dem)
            if (written[t] != 0) {
                const ScanResult d = scan(*dem, m, t);
                r.max_error = d.max_error;
                r.uncovered += d.uncovered;
                if (d.node) {
                    PointScan c;
                    c.point = mesh::MeshVertex{*d.node};
                    c.z = vertex_z(*dem, *c.point).value_or(std::numeric_limits<double>::quiet_NaN());
                    c.error = d.max_error;
                    c.where = d.where;
                    c.is_void = d.is_void;
                    c.set = PointSet::Dem;
                    r.offer(c);
                }
            }
        return r;
    };

    using clock = std::chrono::steady_clock;
    const auto since = [](clock::time_point t0) {
        return std::chrono::duration<double>(clock::now() - t0).count();
    };
    while (true) {
        ++out.rounds;
        results.resize(m.triangle_count());
        auto t0 = clock::now();
        parallel_util::for_each_block(active.size(), options.threads, parallel_util::BlockSchedule{},
                                      [&](std::size_t begin, std::size_t end) {
                                          for (std::size_t i = begin; i < end; ++i)
                                              results[active[i]] = scan_one(active[i]);
                                      });
        out.scan_seconds += since(t0);
        t0 = clock::now();
        std::vector<char> touched(m.triangle_count(), 0);
        std::vector<std::uint32_t> skipped;
        bool any = false;
        for (const std::uint32_t t : active) {
            const PointScan& r = results[t];
            if (!r.point || !(r.is_void || r.error > options.tolerance))
                continue;
            any = true;
            if (touched[t] != 0)
                continue;
            const auto before = static_cast<std::uint32_t>(m.triangle_count());
            std::array<std::uint32_t, 4> seeds{t, before, before + 1, before + 1};
            std::size_t n_seeds = 3;
            std::uint32_t q = 0;
            // L12: with a strip, a point Inside but a hair off a constrained edge goes in on it.
            std::optional<unsigned> edge;
            if (r.where != NodeLocation::Inside)
                edge = static_cast<unsigned>(r.where) - 1;
            else if (strip)
                edge = near_constraint(m, t, *r.point);
            const bool near = edge && r.where == NodeLocation::Inside;
            if (edge) {
                const std::uint32_t u = m.neighbours(t)[*edge];
                if (u != mesh::kNoNeighbour && touched[u] != 0) {
                    skipped.push_back(t);
                    continue;
                }
                if ((r.set == PointSet::Strip || near) && !strip_fits(m, t, *edge, *r.point, frame)) {
                    if (near) {  // L12: back to split_inside, nothing refused
                        edge.reset();
                    } else {  // L1
                        refused[r.strip_index] = 1;
                        ++out.strip_refused;
                        skipped.push_back(t);
                        continue;
                    }
                }
            }
            if (!edge) {
                q = m.split_inside(t, *r.point);
            } else {
                const auto e = *edge;
                const std::uint32_t u = m.neighbours(t)[e];
                const auto& tri = m.triangles()[t];
                const std::uint32_t ea = tri[e], eb = tri[(e + 1) % 3];
                const bool constrained = m.is_constrained(t, e);
                q = m.split_edge(t, e, *r.point);
                if (strip && constrained)
                    cut(subs, *strip, m, ea, eb, q,
                        r.set == PointSet::Strip ? std::optional<double>{r.s} : std::nullopt);
                n_seeds = u != mesh::kNoNeighbour ? 4 : 2;
                if (n_seeds == 4)
                    seeds[3] = u;
            }
            zt.push_back(r.z);
            touched.resize(m.triangle_count(), 1);
            touched[t] = 1;
            if (n_seeds == 4)
                touched[seeds[3]] = 1;
            out.flips += mesh::legalise_around<pred::DefaultKernel>(
                m, q, std::span<const std::uint32_t>{seeds.data(), n_seeds}, frame, flip_stack,
                [&](std::uint32_t s) { touched[s] = 1; });
            ++out.inserted;
            out.carved += r.is_void ? 1 : 0;
            out.strip_inserted += r.set == PointSet::Strip ? 1 : 0;
            out.nodes_inserted += r.set == PointSet::Dem ? 1 : 0;
        }
        out.split_seconds += since(t0);
        if (!any)
            break;
        written.resize(m.triangle_count(), 0);
        for (std::size_t s = 0; s < touched.size(); ++s)
            written[s] = static_cast<char>(written[s] | touched[s]);
        rebuild_active(touched, skipped, active);
    }

    for (const PointScan& r : results) {
        out.uncovered += r.uncovered;
        out.max_error = std::max(out.max_error, r.max_error);
    }
    // Step 6: every strip point against its final sub-edge, refused ones apart.
    for (const auto& [key, se] : subs) {
        if (std::isnan(zt[se.a]) || std::isnan(zt[se.b]))
            continue;
        const auto pts = strip->on_edge(se.k);
        for (std::size_t j = 0; j < pts.size(); ++j)
            if (pts[j].s >= se.s_a && pts[j].s <= se.s_b) {
                double& worst = refused[offset[se.k] + j] != 0 ? out.strip_refused_max_error : out.strip_max_error;
                worst = std::max(worst, std::abs(pts[j].z - along(se, zt, pts[j].s)));
            }
    }
    for (std::size_t i = 0; i < m.vertices().size(); ++i) {
        const mesh::MeshVertex v = m.vertices()[i];
        if constexpr (has_store)
            if (i < n0) {
                const auto row = std::min(static_cast<std::size_t>(v.row), last_cell(g.rows()));
                const auto col = std::min(static_cast<std::size_t>(v.col), last_cell(g.cols()));
                points->for_each_in(row, col, col, [&](mesh::MeshVertex p, float pz) {
                    if (p != v)
                        return;
                    ++out.coincident;
                    if (!std::isnan(zt[i]))
                        out.coincident_max_error =
                            std::max(out.coincident_max_error, std::abs(static_cast<double>(pz) - zt[i]));
                });
            }
        out.vertices.push_back(i < n0 ? start.vertices()[i]
                                      : Point2{g.x_min() + v.col * g.delta_x(), g.y_max() - v.row * g.delta_y()});
        out.z.push_back(std::isnan(zt[i]) ? 0.0 : zt[i]);
        out.valid.push_back(std::isnan(zt[i]) ? 0 : 1);
    }
    out.triangles.assign(m.triangles().begin(), m.triangles().end());
    std::tie(out.edges, out.masks) = m.constraint_edges();
    return out;
}

}  // namespace detail

// Store: CheckPoints, or a test double with geometry(), frozen() and for_each_in.
// An unfrozen store is a programming error (std::logic_error), as add after
// freeze is: RefineStatus has no value for it, and refine.hpp is not edited.
// With a strip (the reprojected path), its points join the loop (D4); a strip
// on another geometry, or with an edge that is not a constraint edge of the
// start, is std::logic_error (L2).
template <class Store>
[[nodiscard]] PointRefineOutcome refine_points(const Store& points, const IndexedMesh2& start,
                                               std::span<const double> z, std::span<const std::uint8_t> valid,
                                               std::span<const std::array<std::uint32_t, 2>> edges,
                                               std::span<const std::uint32_t> masks,
                                               const PointRefineOptions& options,
                                               const ConstraintCheckPoints* strip = nullptr) {
    if (!std::isfinite(options.tolerance) || options.tolerance < 0.0)
        return detail::point_refusal(detail::refusal(RefineStatus::InvalidTolerance,
                                                     "refine_points: tolerance must be finite and >= 0"));
    if (!points.frozen())
        throw std::logic_error("refine_points: the check-point store is not frozen");
    return detail::point_loop("refine_points", points.geometry(), &points,
                              static_cast<const detail::NoSet*>(nullptr), strip, start, z, valid, edges, masks,
                              options);
}

// The projected path (D4, F2): the strip's points, and the DEM's nodes in every
// triangle the run has written (L3), with refine's own scan. Refusals as
// refine_points', against dem.geometry() (L2).
template <raster::RasterSource R>
[[nodiscard]] PointRefineOutcome refine_strip(const R& dem, const ConstraintCheckPoints& strip,
                                              const IndexedMesh2& start, std::span<const double> z,
                                              std::span<const std::uint8_t> valid,
                                              std::span<const std::array<std::uint32_t, 2>> edges,
                                              std::span<const std::uint32_t> masks,
                                              const PointRefineOptions& options) {
    if (!std::isfinite(options.tolerance) || options.tolerance < 0.0)
        return detail::point_refusal(detail::refusal(RefineStatus::InvalidTolerance,
                                                     "refine_strip: tolerance must be finite and >= 0"));
    return detail::point_loop("refine_strip", dem.geometry(), static_cast<const detail::NoSet*>(nullptr), &dem,
                              &strip, start, z, valid, edges, masks, options);
}

}  // namespace terrain::refinement
