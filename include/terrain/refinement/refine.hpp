#pragma once

// Adaptive refinement of a start mesh against a DEM, to a sup-norm tolerance
// (docs/increments/14-adaptive-refinement.md, R5 to R8).
//
// Each round scans the active triangles in parallel (read-only, one result
// slot per triangle), then splits the ones that did not converge serially, in
// triangle-index order. The split phase depends only on the scan results and
// the index order, never on the thread count, so the output is bit-identical
// for any `threads`.
//
// A triangle converges when its error is at most `tolerance`; a void triangle
// (a NoData vertex) when its node set holds no valid node. Each split inserts
// a valid DEM node that is not yet a vertex, and the lattice is finite, so
// the loop terminates for every tolerance, 0 included.
//
// Delaunay insertion (docs/increments/14b-delaunay-insertion.md, R1 to R5).
// The start mesh is legalised once, and each split is followed by Lawson
// legalisation around the new vertex, in the serial phase. Every slot a split
// or a flip writes is `touched` and rescanned next round; an unwritten slot
// still holds the triangle its last scan measured, so its result stays exact.
// The output is constrained Delaunay in the frame (col * dx, -(row * dy)).
//
// Off-node start vertices (docs/increments/16-domain-polygon.md, R2). A start
// vertex may lie anywhere in the node rectangle; one that is a node bit for
// bit is a node, as before, and any other gets fractional (col, row), computed
// once here. It keeps its world point as given for the output, and its z is
// bilinear there (R0), or it is invalid where bilinear refuses.
//
// A quality start (docs/increments/20-start-quality.md, R1). With
// min_angle_deg > 0, mesh::improve runs once after legalise_all and before the
// first scan, and inserts DEM nodes until the start meets the angle or says why
// not. Its nodes are counted in quality_inserted, not in inserted.
//
// Constraint feet (docs/increments/20b-min-insertion-distance.md, R2 to R6).
// With constraint_feet, a worst node N closer than eps(N) to a constrained edge
// of its triangle is replaced by its foot F on that edge, once per node; N may
// still go in later if its error stays above tolerance (R5). An inserted
// off-node vertex is output at (x_min + col dx, y_max - row dy) with vertex_z.

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
#include <terrain/mesh/quality.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace terrain::refinement {

enum class RefineStatus : std::uint8_t { Ok, OutsideGrid, NotCounterClockwise, InvalidTolerance };

struct RefineOptions {
    double tolerance = 0.0;  // metres, finite and >= 0
    unsigned threads = 0;    // 0: hardware concurrency
    double min_angle_deg = 0.0;  // the start-quality pass; 0 (or NaN) is off
    bool constraint_feet = false;  // 20b R9: feet on constraint segments
};

struct RefineOutcome {
    RefineStatus status = RefineStatus::Ok;
    std::string message;  // empty on Ok

    std::vector<Point2> vertices;  // world: a start vertex as given, else RasterGeometry::node
    std::vector<double> z;         // value_at for a node, else bilinear; 0.0 where there is none
    std::vector<std::uint8_t> valid;
    std::vector<TriangleIndices> triangles;
    std::vector<std::array<std::uint32_t, 2>> edges;  // constraint edges, each once
    std::vector<std::uint32_t> masks;                 // one per edge

    std::size_t rounds = 0;     // scans run; 1 when nothing needed a split
    std::size_t inserted = 0;   // vertices added
    std::size_t flips = 0;      // Lawson flips, the start mesh's included
    double max_error = 0.0;     // over triangles with three valid vertices
    std::size_t uncovered = 0;  // valid nodes still inside void triangles
    std::size_t carved = 0;     // inserts that split a void triangle, a subset of `inserted`
    std::size_t quality_inserted = 0;  // start-quality nodes, not in `inserted`
    std::size_t quality_skipped = 0;   // start-quality skips, every reason summed
    std::size_t feet = 0;              // feet inserted, a subset of `inserted`
    std::size_t feet_refused = 0;      // 20b R2 step 4: a foot refused, N inserted instead

    // Wall seconds on the calling thread, steady_clock (17-mesh-stats.md R6).
    // for_each_chunk joins its workers before returning, so no worker reads a
    // clock. Not part of the determinism guarantee.
    double legalise_seconds = 0.0;  // the start mesh's legalise_all
    double scan_seconds = 0.0;      // every round's parallel scan, summed
    double split_seconds = 0.0;     // every round's serial split + flip phase, summed
    double quality_seconds = 0.0;   // the start-quality pass

    [[nodiscard]] bool ok() const noexcept { return status == RefineStatus::Ok; }
};

namespace detail {

[[nodiscard]] inline RefineOutcome refusal(RefineStatus status, std::string message) {
    RefineOutcome out;
    out.status = status;
    out.message = std::move(message);
    return out;
}

// The lattice mesh for `start`, or the refusal. A vertex is a node iff
// RasterGeometry::node of its rounded (row, col) gives it back bit for bit;
// otherwise it is off-node, and refused only outside the node rectangle.
[[nodiscard]] inline std::variant<mesh::LatticeMesh, RefineOutcome> to_lattice(
    const raster::RasterGeometry& g, const IndexedMesh2& start,
    std::span<const std::array<std::uint32_t, 2>> edges, std::span<const std::uint32_t> masks) {
    std::vector<mesh::MeshVertex> lattice;
    lattice.reserve(start.vertices().size());
    for (std::size_t i = 0; i < start.vertices().size(); ++i) {
        const Point2 p = start.vertices()[i];
        if (!g.cell_of(p))
            return refusal(RefineStatus::OutsideGrid, "refine: start vertex " + std::to_string(i)
                                                          + " is outside the DEM's node rectangle");
        const double col = std::clamp((p.x - g.x_min()) / g.delta_x(), 0.0,
                                      static_cast<double>(g.cols() - 1));
        const double row = std::clamp((g.y_max() - p.y) / g.delta_y(), 0.0,
                                      static_cast<double>(g.rows() - 1));
        const raster::CellIndex node{static_cast<std::size_t>(std::round(row)),
                                     static_cast<std::size_t>(std::round(col))};
        lattice.push_back(g.node(node) == p ? mesh::MeshVertex{static_cast<double>(node.col),
                                                               static_cast<double>(node.row)}
                                            : mesh::MeshVertex{col, row});
    }

    std::map<std::pair<std::uint32_t, std::uint32_t>, std::uint32_t> constraint;
    for (std::size_t i = 0; i < edges.size() && i < masks.size(); ++i)
        constraint[std::minmax(edges[i][0], edges[i][1])] = masks[i];
    std::vector<std::uint8_t> bits(start.triangle_count(), 0);
    std::vector<std::array<std::uint32_t, 3>> edge_masks(start.triangle_count(), {0, 0, 0});
    for (std::size_t t = 0; t < start.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto& tri = start.triangles()[t];
            if (const auto it = constraint.find(std::minmax(tri[k], tri[(k + 1) % 3]));
                it != constraint.end()) {
                bits[t] |= static_cast<std::uint8_t>(1u << k);
                edge_masks[t][k] = it->second;
            }
        }

    auto m = mesh::LatticeMesh::build(std::move(lattice), {start.triangles().begin(),
                                                           start.triangles().end()},
                                      std::move(bits), std::move(edge_masks));
    if (!m)
        return refusal(RefineStatus::NotCounterClockwise,
                       "refine: a start triangle is not counter-clockwise, or the start mesh "
                       "is not a valid triangulation");
    return std::move(*m);
}

[[nodiscard]] inline bool needs_split(const ScanResult& r, double tolerance) noexcept {
    return r.node.has_value() && (r.is_void || r.max_error > tolerance);
}

// R3: eps(n) = clamp(tol / G, floor, cap), G the largest bilinear slope bound
// over the valid cells sharing n; the cap where G is 0 (flat, or no cell).
template <raster::RasterSource R>
[[nodiscard]] double foot_epsilon(const R& dem, mesh::LatticeVertex n, double tol) {
    const raster::RasterGeometry& g = dem.geometry();
    const double dx = g.delta_x(), dy = g.delta_y();
    double grad = 0.0;
    for (std::uint32_t r = n.row > 0 ? n.row - 1 : 0; r <= n.row && r + 1 < g.rows(); ++r)
        for (std::uint32_t c = n.col > 0 ? n.col - 1 : 0; c <= n.col && c + 1 < g.cols(); ++c) {
            const auto z00 = vertex_z(dem, mesh::LatticeVertex{r, c}),
                       z01 = vertex_z(dem, mesh::LatticeVertex{r, c + 1}),
                       z10 = vertex_z(dem, mesh::LatticeVertex{r + 1, c}),
                       z11 = vertex_z(dem, mesh::LatticeVertex{r + 1, c + 1});
            if (!z00 || !z01 || !z10 || !z11)
                continue;
            grad = std::max(grad, std::hypot(std::max(std::abs(*z01 - *z00), std::abs(*z11 - *z10)) / dx,
                                             std::max(std::abs(*z10 - *z00), std::abs(*z11 - *z01)) / dy));
        }
    const double cap = std::min(dx, dy) / 2.0, floor = std::min(dx, dy) / 100.0;
    return grad > 0.0 ? std::clamp(tol / grad, floor, cap) : cap;
}

struct Foot {
    unsigned edge;
    mesh::MeshVertex at;
};

// R2 steps 1 and 2: on the first constrained edge of t that n is closer to
// than eps (in world distance), the foot, or nothing when that foot is within
// eps of an end or no such edge exists. n exactly on an edge is skipped.
template <raster::RasterSource R>
[[nodiscard]] std::optional<Foot> foot_of(const R& dem, const mesh::LatticeMesh& m,
                                          std::uint32_t t, mesh::LatticeVertex n, double tol) {
    const double dx = dem.geometry().delta_x(), dy = dem.geometry().delta_y();
    std::optional<double> eps;
    for (unsigned e = 0; e < 3; ++e) {
        const mesh::MeshVertex a = m.corner(t, e), b = m.corner(t, (e + 1) % 3), p{n};
        if (!m.is_constrained(t, e) || mesh::orient_sign(a, b, p) == 0)
            continue;
        const double ux = (b.col - a.col) * dx, uy = (b.row - a.row) * dy;
        const double px = (p.col - a.col) * dx, py = (p.row - a.row) * dy;
        const double s = std::clamp((px * ux + py * uy) / (ux * ux + uy * uy), 0.0, 1.0);
        if (!eps)
            eps = foot_epsilon(dem, n, tol);
        if (std::hypot(px - s * ux, py - s * uy) >= *eps)
            continue;
        const double len = std::hypot(ux, uy);
        if (s * len < *eps || (1.0 - s) * len < *eps)
            return std::nullopt;
        return Foot{e, {a.col + s * (b.col - a.col), a.row + s * (b.row - a.row)}};
    }
    return std::nullopt;
}

// R2 step 4: every child of splitting t's edge e (and its neighbour's) at f
// is strictly counter-clockwise.
[[nodiscard]] inline bool foot_fits(const mesh::LatticeMesh& m, std::uint32_t t, unsigned e,
                                    mesh::MeshVertex f) {
    const mesh::MeshVertex a = m.corner(t, e), b = m.corner(t, (e + 1) % 3),
                           c = m.corner(t, (e + 2) % 3);
    if (mesh::orient_sign(a, f, c) <= 0 || mesh::orient_sign(f, b, c) <= 0)
        return false;
    const std::uint32_t u = m.neighbours(t)[e];
    if (u == mesh::kNoNeighbour)
        return true;
    unsigned k = 0;
    while (m.neighbours(u)[k] != t)
        ++k;
    const mesh::MeshVertex d = m.corner(u, (k + 2) % 3);
    return mesh::orient_sign(b, f, d) > 0 && mesh::orient_sign(f, a, d) > 0;
}

}  // namespace detail

template <raster::RasterSource R>
[[nodiscard]] RefineOutcome refine(const R& dem, const IndexedMesh2& start,
                                   std::span<const std::array<std::uint32_t, 2>> edges,
                                   std::span<const std::uint32_t> masks,
                                   const RefineOptions& options) {
    if (!std::isfinite(options.tolerance) || options.tolerance < 0.0)
        return detail::refusal(RefineStatus::InvalidTolerance,
                               "refine: tolerance must be finite and >= 0");
    const raster::RasterGeometry& g = dem.geometry();
    auto built = detail::to_lattice(g, start, edges, masks);
    if (auto* refused = std::get_if<RefineOutcome>(&built))
        return std::move(*refused);
    auto& m = std::get<mesh::LatticeMesh>(built);

    const mesh::LatticeFrame frame{g.delta_x(), g.delta_y()};
    RefineOutcome out;
    using clock = std::chrono::steady_clock;
    const auto since = [](clock::time_point t0) {
        return std::chrono::duration<double>(clock::now() - t0).count();
    };
    auto t0 = clock::now();
    out.flips = mesh::legalise_all<pred::DefaultKernel>(m, frame, [](std::uint32_t) {});
    out.legalise_seconds = since(t0);
    if (options.min_angle_deg > 0.0) {
        t0 = clock::now();
        const auto q = mesh::improve<pred::DefaultKernel>(
            m, frame, mesh::QualityOptions{options.min_angle_deg, g.rows(), g.cols()});
        out.quality_inserted = q.inserted;
        out.quality_skipped = q.skipped_floor + q.skipped_outside + q.skipped_vertex
                            + q.skipped_blocked + q.walk_bound_hits;
        out.quality_seconds = since(t0);
    }
    std::vector<ScanResult> results;
    std::set<std::pair<std::uint32_t, std::uint32_t>> footed;  // (row, col), R2 step 5
    std::vector<std::uint32_t> active(m.triangle_count());
    for (std::uint32_t t = 0; t < active.size(); ++t)
        active[t] = t;

    while (true) {
        ++out.rounds;
        results.resize(m.triangle_count());
        t0 = clock::now();
        parallel_util::for_each_chunk(active.size(), options.threads,
                                      [&](std::size_t begin, std::size_t end) {
                                          for (std::size_t i = begin; i < end; ++i)
                                              results[active[i]] = scan(dem, m, active[i]);
                                      });
        out.scan_seconds += since(t0);
        t0 = clock::now();

        // `touched` is every slot a split or a flip wrote this round. A marked triangle
        // skipped because its neighbour was touched is itself unchanged and
        // still unconverged, so it stays active rather than being forgotten.
        std::vector<char> touched(m.triangle_count(), 0);
        std::vector<std::uint32_t> skipped;
        bool any = false;
        for (const std::uint32_t t : active) {
            const ScanResult& r = results[t];
            if (!detail::needs_split(r, options.tolerance))
                continue;
            any = true;
            if (touched[t] != 0)
                continue;
            const auto before = static_cast<std::uint32_t>(m.triangle_count());
            std::array<std::uint32_t, 4> seeds{t, before, before + 1, before + 1};
            std::size_t n_seeds = 3;
            std::uint32_t q = 0;
            std::optional<unsigned> edge;
            if (r.where != NodeLocation::Inside)
                edge = static_cast<unsigned>(r.where) - 1;
            mesh::MeshVertex p = *r.node;
            bool foot = false, refused = false;
            if (options.constraint_feet && !r.is_void && !footed.contains({r.node->row, r.node->col}))
                if (const auto f = detail::foot_of(dem, m, t, *r.node, options.tolerance)) {
                    refused = !vertex_z(dem, f->at) || !detail::foot_fits(m, t, f->edge, f->at);
                    foot = !refused;
                    if (foot) {
                        edge = f->edge;
                        p = f->at;
                    }
                }
            if (!edge) {
                q = m.split_inside(t, *r.node);
            } else {
                const auto e = *edge;
                const std::uint32_t u = m.neighbours(t)[e];
                if (u != mesh::kNoNeighbour && touched[u] != 0) {
                    skipped.push_back(t);
                    continue;
                }
                q = m.split_edge(t, e, p);
                n_seeds = u != mesh::kNoNeighbour ? 4 : 2;
                if (n_seeds == 4)
                    seeds[3] = u;
            }
            touched.resize(m.triangle_count(), 1);
            touched[t] = 1;
            if (n_seeds == 4)
                touched[seeds[3]] = 1;
            out.flips += mesh::legalise_around<pred::DefaultKernel>(
                m, q, std::span<const std::uint32_t>{seeds.data(), n_seeds}, frame,
                [&](std::uint32_t s) { touched[s] = 1; });
            ++out.inserted;
            out.carved += r.is_void ? 1 : 0;
            out.feet += foot ? 1 : 0;
            out.feet_refused += refused ? 1 : 0;
            if (foot)
                footed.insert({r.node->row, r.node->col});
        }
        out.split_seconds += since(t0);
        if (!any)
            break;
        active.clear();
        for (std::uint32_t t = 0; t < touched.size(); ++t)
            if (touched[t] != 0)
                active.push_back(t);
        active.insert(active.end(), skipped.begin(), skipped.end());
        std::sort(active.begin(), active.end());
        active.erase(std::unique(active.begin(), active.end()), active.end());
    }

    // By the stopping rule a void triangle holds no valid node, so `uncovered`
    // sums zeros unless that rule changes; it is reported so a change shows.
    for (const ScanResult& r : results) {
        if (r.is_void)
            out.uncovered += r.uncovered;
        else
            out.max_error = std::max(out.max_error, r.max_error);
    }
    for (std::size_t i = 0; i < m.vertices().size(); ++i) {
        const mesh::MeshVertex v = m.vertices()[i];
        const bool node = v.is_node();
        const raster::CellIndex c{node ? static_cast<std::size_t>(v.row) : 0,
                                  node ? static_cast<std::size_t>(v.col) : 0};
        const bool given = i < start.vertices().size();
        const Point2 p = given ? start.vertices()[i]
                               : Point2{g.x_min() + v.col * g.delta_x(), g.y_max() - v.row * g.delta_y()};
        const std::optional<double> z =
            !given || (node && g.node(c) == p) ? vertex_z(dem, v) : raster::bilinear(dem, p);
        out.vertices.push_back(p);
        out.z.push_back(z.value_or(0.0));
        out.valid.push_back(z ? 1 : 0);
    }
    out.triangles.assign(m.triangles().begin(), m.triangles().end());
    std::tie(out.edges, out.masks) = m.constraint_edges();
    return out;
}

}  // namespace terrain::refinement
