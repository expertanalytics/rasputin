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
//
// Frozen edges (docs/increments/23-basin-scale.md, N1, N2, N6, N7, N16). With
// frozen_mask, a stored point on a frozen edge of its triangle (after L14's
// corner test; with a strip, also one within r(g) of it, projection strictly
// inside) is never named. It is counted once in `on_frozen`, by the triangle
// that owns the edge when it lies exactly on it (15f D4's rule), with its
// error against the edge's linear z at its projection. L12 never takes a
// frozen edge, the DEM rescan is refine's scan, and a frozen strip edge is a
// programming error (std::logic_error).
//
// Constraint feet (docs/increments/20c-soft-quality.md, R4). With
// constraint_feet, a source point or DEM node that L12 left Inside and that
// lies within delta_p of a constraint (mesh::constraint_foot) goes in as its
// foot, once; if its error stays above the tolerance it goes in later as
// itself (feet_fallback). Feet are counted in `inserted`.
//
// A tolerance policy (docs/increments/33-feature-tolerance.md, 4.4): each
// triangle's error is compared with the policy's allowed error, asked only
// between lowest() and highest(), as refine does; options.tolerance is not
// read by the overloads that take one.

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
#include <terrain/refinement/strip_scan.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace terrain::refinement {

struct PointRefineOptions {
    double tolerance = 0.0;  // metres, finite and >= 0
    unsigned threads = 0;    // 0: hardware concurrency
    std::uint32_t frozen_mask = 0;  // 23b: edges whose mask meets it are never split
    bool constraint_feet = false;   // 20c R4.6: a point near a constraint goes in as its foot first
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
    std::size_t on_frozen = 0;          // stored points on a frozen edge, each once (N6)
    double on_frozen_max_error = 0.0;   // their largest error, edges with two valid ends only
    std::size_t feet_fallback = 0;      // 20c R4.5: footed points later inserted as themselves
};

namespace detail {

template <class Store>
[[nodiscard]] PointScan scan_points(const Store& points, const mesh::LatticeMesh& m,
                                    std::span<const double> zt, std::uint32_t t, double radius = 0.0) {
    using mesh::MeshVertex;
    const raster::RasterGeometry& g = points.geometry();
    const auto& tri = m.triangles()[t];
    const std::array<MeshVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
    const std::array<double, 3> zv{zt[tri[0]], zt[tri[1]], zt[tri[2]]};
    PointScan r;
    r.is_void = std::isnan(zv[0]) || std::isnan(zv[1]) || std::isnan(zv[2]);
    const bool frozen = m.is_frozen(t, 0) || m.is_frozen(t, 1) || m.is_frozen(t, 2);
    const Point2 f0 = v[0].frame(), f1 = v[1].frame(), f2 = v[2].frame();
    const double two_a = (f1.x - f0.x) * (f2.y - f0.y) - (f1.y - f0.y) * (f2.x - f0.x);
    auto value = [&](unsigned k, MeshVertex p) {
        const Point2 a = v[k].frame(), b = v[(k + 1) % 3].frame(), q = p.frame();
        return (b.x - a.x) * (q.y - a.y) - (b.y - a.y) * (q.x - a.x);
    };
    double nearest = std::numeric_limits<double>::infinity();
    auto visit = [&](MeshVertex p, float zf) {
        for (unsigned k = 0; k < 3; ++k)  // a corner, or within `radius` of one (L14)
            if (p == v[k] || (radius > 0.0 && std::hypot(p.col - v[k].col, p.row - v[k].row) <= radius))
                return;
        for (unsigned k = 0; k < 3; ++k)
            if (mesh::orient_sign(v[k], v[(k + 1) % 3], p) < 0)
                return;
        const auto z = static_cast<double>(zf);
        if (const auto e = frozen ? frozen_edge_at(m, t, p, radius) : std::nullopt) {
            const unsigned f = (*e + 1) % 3;
            if (mesh::orient_sign(v[*e], v[f], p) != 0 || tri[*e] < tri[f] || m.neighbours(t)[*e] == mesh::kNoNeighbour) {
                ++r.on_frozen;
                const double dc = v[f].col - v[*e].col, dr = v[f].row - v[*e].row;
                const double sigma = ((p.col - v[*e].col) * dc + (p.row - v[*e].row) * dr) / (dc * dc + dr * dr);
                if (!std::isnan(zv[*e]) && !std::isnan(zv[f]))
                    r.frozen_error = std::max(r.frozen_error, std::abs(z - (zv[*e] + sigma * (zv[f] - zv[*e]))));
            }
            return;
        }
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

struct NoSet {};  // point_loop without a store, or without a DEM

[[nodiscard]] inline PointRefineOutcome point_refusal(RefineOutcome r) {
    PointRefineOutcome out;
    static_cast<RefineOutcome&>(out) = std::move(r);
    return out;
}

// The loop of refine_points and refine_strip (D4). `points` (a store) or `dem`
// may be null, `strip` too; refusals (2) to (4) of L2, the tolerance being the
// caller's.
template <class Store, class R, TolerancePolicy P>
[[nodiscard]] PointRefineOutcome point_loop(const std::string& name, const raster::RasterGeometry& g,
                                            const Store* points, const R* dem, const ConstraintCheckPoints* strip,
                                            const IndexedMesh2& start, std::span<const double> z,
                                            std::span<const std::uint8_t> valid,
                                            std::span<const std::array<std::uint32_t, 2>> edges,
                                            std::span<const std::uint32_t> masks, const PointRefineOptions& options,
                                            const P& policy) {
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
        for (std::size_t i = 0; i < edges.size() && i < masks.size(); ++i)  // N16, after check (3)
            if ((masks[i] & options.frozen_mask) != 0 && subs.contains(edge_key(edges[i][0], edges[i][1])))
                throw std::logic_error(name + ": strip edge (" + std::to_string(std::min(edges[i][0], edges[i][1]))
                                       + ", " + std::to_string(std::max(edges[i][0], edges[i][1])) + ") is frozen");
    }
    auto built = to_lattice(g, start, edges, masks);
    if (auto* refused = std::get_if<RefineOutcome>(&built))
        return point_refusal(std::move(*refused));
    auto& m = std::get<mesh::LatticeMesh>(built);
    m.set_frozen_mask(options.frozen_mask);

    PointRefineOutcome out;
    out.strip_points = strip ? strip->size() : 0;
    const double radius = strip ? coincidence_radius(g) : 0.0;  // L14, L16; 0 keeps 15c's run
    const std::size_t n0 = start.vertices().size();
    std::vector<double> zt(n0, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < n0 && i < z.size() && i < valid.size(); ++i)
        if (valid[i] != 0)
            zt[i] = z[i];

    const auto frame = mesh::lattice_frame(g.delta_x(), g.delta_y(), g.rows(), g.cols());
    std::vector<char> written(m.triangle_count(), 0), refused(offset.back(), 0);  // L3, L1
    out.flips = mesh::legalise_all<pred::DefaultKernel>(m, frame, [&](std::uint32_t s) { written[s] = 1; });
    std::vector<PointScan> results;
    [[maybe_unused]] Allowed<P> allowed;
    mesh::FlipStack flip_stack;
    std::vector<std::uint32_t> active(m.triangle_count());
    for (std::uint32_t t = 0; t < active.size(); ++t)
        active[t] = t;

    // 20c R4: feet on constraints, each point footed once (by exact position),
    // delta_p half the smaller cell side. R4.3, the foot's z: with a raster,
    // vertex_z; else linear in s between the strip points bracketing the foot
    // on its sub-edge (the sub-edge's end where there is none), or linear
    // between the edge's ends where it has no strip record. NaN refuses it.
    std::set<std::pair<double, double>> footed;
    const double delta_p = std::min(g.delta_x(), g.delta_y()) / 2.0;
    const auto foot_z = [&](const mesh::FootSearch& s) {
        if constexpr (has_dem) {
            return vertex_z(*dem, s.at).value_or(std::numeric_limits<double>::quiet_NaN());
        } else {
            const auto sigma = [&](std::uint32_t a, std::uint32_t b) {
                const mesh::MeshVertex va = m.vertices()[a], vb = m.vertices()[b];
                const double dc = vb.col - va.col, dr = vb.row - va.row;
                return ((s.at.col - va.col) * dc + (s.at.row - va.row) * dr) / (dc * dc + dr * dr);
            };
            const auto& tri = m.triangles()[s.owner];
            const std::uint32_t ea = tri[s.edge], eb = tri[(s.edge + 1) % 3];
            const auto it = strip ? subs.find(edge_key(ea, eb)) : subs.end();
            if (it == subs.end())
                return zt[ea] + sigma(ea, eb) * (zt[eb] - zt[ea]);
            const SubEdge& se = it->second;
            const double sf = se.s_a + sigma(se.a, se.b) * (se.s_b - se.s_a);
            double s0 = se.s_a, z0 = zt[se.a], s1 = se.s_b, z1 = zt[se.b];
            for (const ConstraintPoint& c : strip->on_edge(se.k)) {
                if (c.s > se.s_a && c.s <= sf)
                    std::tie(s0, z0) = std::pair{c.s, c.z};
                else if (c.s > sf && c.s < s1)
                    std::tie(s1, z1) = std::pair{c.s, c.z};
            }
            return s1 > s0 ? z0 + (sf - s0) / (s1 - s0) * (z1 - z0) : z0;
        }
    };

    // Source, strip, DEM, in that order (D4, "Combining").
    const auto scan_one = [&](std::uint32_t t) {
        PointScan r;
        if constexpr (has_store) {
            r = scan_points(*points, m, zt, t, radius);
            r.error = r.max_error;
        }
        if (strip)
            scan_strip(*strip, offset, refused, subs, m, zt, t, r);
        if constexpr (has_dem)
            if (written[t] != 0) {
                const ScanResult d = scan(*dem, m, t, radius);
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
        if constexpr (varies<P>)
            allowed.resize(m.triangle_count());
        auto t0 = clock::now();
        parallel_util::for_each_block(active.size(), options.threads, parallel_util::BlockSchedule{},
                                      [&](std::size_t begin, std::size_t end) {
                                          for (std::size_t i = begin; i < end; ++i) {
                                              results[active[i]] = scan_one(active[i]);
                                              if constexpr (varies<P>)
                                                  allowed[active[i]] =
                                                      allowed_at(policy, m, active[i], results[active[i]].error);
                                          }
                                      });
        out.scan_seconds += since(t0);
        t0 = clock::now();
        std::vector<char> touched(m.triangle_count(), 0);
        std::vector<std::uint32_t> skipped;
        bool any = false;
        for (const std::uint32_t t : active) {
            const PointScan& r = results[t];
            double limit = policy.lowest();
            if constexpr (varies<P>)
                limit = allowed[t];
            if (!r.point || !(r.is_void || r.error > limit))
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
                edge = near_constraint(m, t, *r.point, radius);
            const bool near = edge && r.where == NodeLocation::Inside;
            // 20c R4: a source point or DEM node L12 left Inside, not yet footed, may go in as its foot.
            std::uint32_t owner = t;
            mesh::MeshVertex p = *r.point;
            double pz = r.z;
            std::size_t refused_feet = 0;
            std::optional<mesh::FootSearch> foot;
            const bool was_footed = r.set != PointSet::Strip && footed.contains({p.col, p.row});
            if (options.constraint_feet && !edge && !r.is_void && r.set != PointSet::Strip && !was_footed) {
                const auto s = mesh::constraint_foot(m, t, p, delta_p, frame);
                const double fz = s.status == mesh::FootStatus::Hit ? foot_z(s) : 0.0;
                const bool fits = s.status == mesh::FootStatus::Hit && !std::isnan(fz)
                                  && (!strip || strip_fits(m, s.owner, s.edge, s.at, frame));
                if ((foot = usable(s, fits, refused_feet)))
                    std::tie(owner, edge, p, pz) = std::tuple{foot->owner, foot->edge, foot->at, fz};
            }
            if (edge) {
                const std::uint32_t u = m.neighbours(owner)[*edge];
                if (touched[owner] != 0 || (u != mesh::kNoNeighbour && touched[u] != 0)) {
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
                const std::uint32_t u = m.neighbours(owner)[e];
                const auto& tri = m.triangles()[owner];
                const std::uint32_t ea = tri[e], eb = tri[(e + 1) % 3];
                const bool constrained = m.is_constrained(owner, e);
                q = m.split_edge(owner, e, p);
                if (strip && constrained)
                    cut(subs, *strip, m, ea, eb, q,
                        r.set == PointSet::Strip ? std::optional<double>{r.s} : std::nullopt);
                seeds[0] = owner;
                n_seeds = u != mesh::kNoNeighbour ? 4 : 2;
                if (n_seeds == 4)
                    seeds[3] = u;
            }
            if (owner != t)  // R4.2, as R3: t is unchanged and its point still a candidate
                skipped.push_back(t);
            zt.push_back(pz);
            touched.resize(m.triangle_count(), 1);
            touched[owner] = 1;
            if (n_seeds == 4)
                touched[seeds[3]] = 1;
            out.flips += mesh::legalise_around<pred::DefaultKernel>(
                m, q, std::span<const std::uint32_t>{seeds.data(), n_seeds}, frame, flip_stack,
                [&](std::uint32_t s) { touched[s] = 1; });
            ++out.inserted;
            out.carved += r.is_void ? 1 : 0;
            out.strip_inserted += r.set == PointSet::Strip ? 1 : 0;
            out.nodes_inserted += r.set == PointSet::Dem && !foot ? 1 : 0;
            out.feet += foot ? 1 : 0;
            out.feet_fallback += was_footed ? 1 : 0;
            out.feet_refused += refused_feet;
            if (foot)
                footed.insert({r.point->col, r.point->row});
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
        out.on_frozen += r.on_frozen;
        out.on_frozen_max_error = std::max(out.on_frozen_max_error, r.frozen_error);
    }
    out.max_error_near = out.max_error;
    if constexpr (varies<P>) {
        out.max_error_near = 0.0;
        for (std::uint32_t t = 0; t < results.size(); ++t)
            if (results[t].max_error > out.max_error_near && policy.at(m, t) <= policy.lowest())
                out.max_error_near = results[t].max_error;
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
    const auto coincide = [&](std::size_t i, double pz) {
        ++out.coincident;
        if (!std::isnan(zt[i]))
            out.coincident_max_error = std::max(out.coincident_max_error, std::abs(pz - zt[i]));
    };
    std::vector<std::uint64_t> node_vertices;  // the vertices that are DEM nodes, as row * cols + col
    if constexpr (has_dem) {
        for (const mesh::MeshVertex v : m.vertices())
            if (v.is_node())
                node_vertices.push_back(static_cast<std::uint64_t>(v.row) * g.cols() + static_cast<std::uint64_t>(v.col));
        std::sort(node_vertices.begin(), node_vertices.end());
    }
    for (std::size_t i = 0; i < m.vertices().size(); ++i) {
        const mesh::MeshVertex v = m.vertices()[i];
        // Coincident points: 15c's pass, a stored point equal to a start vertex;
        // with a strip (L14), also any stored point or DEM node within the
        // radius of a vertex it is not, which the scan skipped. A stored point
        // within the radius of two vertices is counted twice, and a stored
        // point that is itself a vertex (inserted) is counted against a strip
        // vertex inserted within the radius of it. Both are rounding-scale
        // events; the figure is a report, not an invariant.
        if constexpr (has_store)
            if (i < n0 || radius > 0.0) {
                const auto cell = [&](double x, std::size_t n) {
                    return std::min(static_cast<std::size_t>(std::max(x, 0.0)), last_cell(n));
                };
                for (std::size_t b = cell(v.row - radius, g.rows()); b <= cell(v.row + radius, g.rows()); ++b)
                    points->for_each_in(b, cell(v.col - radius, g.cols()), cell(v.col + radius, g.cols()),
                                        [&](mesh::MeshVertex p, float pz) {
                                            const bool same = p == v;
                                            if ((same && i < n0) || (!same && radius > 0.0
                                                                     && std::hypot(p.col - v.col, p.row - v.row) <= radius))
                                                coincide(i, static_cast<double>(pz));
                                        });
            }
        if constexpr (has_dem) {
            const mesh::MeshVertex n{std::round(v.col), std::round(v.row)};
            const auto key = static_cast<std::uint64_t>(n.row) * g.cols() + static_cast<std::uint64_t>(n.col);
            if (!v.is_node() && std::hypot(v.col - n.col, v.row - n.row) <= radius
                && !std::binary_search(node_vertices.begin(), node_vertices.end(), key))
                if (const auto nz = vertex_z(*dem, n))
                    coincide(i, *nz);
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
template <class Store, TolerancePolicy P>
[[nodiscard]] PointRefineOutcome refine_points(const Store& points, const IndexedMesh2& start,
                                               std::span<const double> z, std::span<const std::uint8_t> valid,
                                               std::span<const std::array<std::uint32_t, 2>> edges,
                                               std::span<const std::uint32_t> masks,
                                               const PointRefineOptions& options, const P& policy,
                                               const ConstraintCheckPoints* strip = nullptr) {
    if (detail::bad_policy(policy))
        return detail::point_refusal(detail::refusal(RefineStatus::InvalidTolerance,
                                                     "refine_points: tolerance must be finite and >= 0"));
    if (!points.frozen())
        throw std::logic_error("refine_points: the check-point store is not frozen");
    return detail::point_loop("refine_points", points.geometry(), &points,
                              static_cast<const detail::NoSet*>(nullptr), strip, start, z, valid, edges, masks,
                              options, policy);
}

template <class Store>
[[nodiscard]] PointRefineOutcome refine_points(const Store& points, const IndexedMesh2& start,
                                               std::span<const double> z, std::span<const std::uint8_t> valid,
                                               std::span<const std::array<std::uint32_t, 2>> edges,
                                               std::span<const std::uint32_t> masks,
                                               const PointRefineOptions& options,
                                               const ConstraintCheckPoints* strip = nullptr) {
    return refine_points(points, start, z, valid, edges, masks, options, UniformTolerance{options.tolerance}, strip);
}

// The projected path (D4, F2): the strip's points, and the DEM's nodes in every
// triangle the run has written (L3), with refine's own scan. Refusals as
// refine_points', against dem.geometry() (L2).
template <raster::RasterSource R, TolerancePolicy P>
[[nodiscard]] PointRefineOutcome refine_strip(const R& dem, const ConstraintCheckPoints& strip,
                                              const IndexedMesh2& start, std::span<const double> z,
                                              std::span<const std::uint8_t> valid,
                                              std::span<const std::array<std::uint32_t, 2>> edges,
                                              std::span<const std::uint32_t> masks,
                                              const PointRefineOptions& options, const P& policy) {
    if (detail::bad_policy(policy))
        return detail::point_refusal(detail::refusal(RefineStatus::InvalidTolerance,
                                                     "refine_strip: tolerance must be finite and >= 0"));
    return detail::point_loop("refine_strip", dem.geometry(), static_cast<const detail::NoSet*>(nullptr), &dem,
                              &strip, start, z, valid, edges, masks, options, policy);
}

template <raster::RasterSource R>
[[nodiscard]] PointRefineOutcome refine_strip(const R& dem, const ConstraintCheckPoints& strip,
                                              const IndexedMesh2& start, std::span<const double> z,
                                              std::span<const std::uint8_t> valid,
                                              std::span<const std::array<std::uint32_t, 2>> edges,
                                              std::span<const std::uint32_t> masks,
                                              const PointRefineOptions& options) {
    return refine_strip(dem, strip, start, z, valid, edges, masks, options, UniformTolerance{options.tolerance});
}

}  // namespace terrain::refinement
