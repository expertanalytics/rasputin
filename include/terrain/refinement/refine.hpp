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
// of its triangle or of a neighbour (mesh::constraint_foot, 20c R3) is replaced
// by its foot F on that edge, once per node; N may still go in later if its
// error stays above tolerance (R5). When F is on a neighbour's edge, N's
// triangle stays active. The quality start foots too (20c R2), counted in
// quality_feet. An inserted off-node vertex is output at
// (x_min + col dx, y_max - row dy) with vertex_z.
//
// Frozen edges (docs/increments/23-basin-scale.md, "Refine with the seam
// frozen", N1-N5). An edge whose mask meets frozen_mask gets no vertex: the
// scan skips the nodes on it, no foot is taken on it, and the quality pass
// skips a node on it. frozen_mask 0 is the run as before, bit for bit (K1).

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/constraint_foot.hpp>
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
#include <iterator>
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
    std::uint32_t frozen_mask = 0;  // 23b: edges whose mask meets it are never split
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
    std::size_t quality_feet = 0;      // start-quality feet (20c R2), in neither count above
    std::size_t feet = 0;              // feet inserted, a subset of `inserted`
    std::size_t feet_refused = 0;      // 20b R2 step 4: a foot refused, N inserted instead

    // Wall seconds on the calling thread, steady_clock (17-mesh-stats.md R6).
    // for_each_block joins its workers before returning, so no worker reads a
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

// The lattice position of a world point in the node rectangle (the caller
// checks): the node when RasterGeometry::node of the rounded (row, col) gives
// p back bit for bit, otherwise the clamped fractional (col, row). to_lattice
// and the edge strip's generator (constraint_points.hpp) share it, so the
// generator's edge ends are the loop's vertices bit for bit.
[[nodiscard]] inline mesh::MeshVertex lattice_position(const raster::RasterGeometry& g, Point2 p) {
    const double col =
        std::clamp((p.x - g.x_min()) / g.delta_x(), 0.0, static_cast<double>(g.cols() - 1));
    const double row =
        std::clamp((g.y_max() - p.y) / g.delta_y(), 0.0, static_cast<double>(g.rows() - 1));
    const raster::CellIndex node{static_cast<std::size_t>(std::round(row)),
                                 static_cast<std::size_t>(std::round(col))};
    return g.node(node) == p
               ? mesh::MeshVertex{static_cast<double>(node.col), static_cast<double>(node.row)}
               : mesh::MeshVertex{col, row};
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
        lattice.push_back(lattice_position(g, p));
    }

    // The constraint lookup is a sorted vector, not a std::map, for the
    // reason LatticeMesh::build's table is flat. The sort is stable, so the
    // last entry of an edge listed twice is the later mask, as the map's
    // assignment gave; an edge found sets its bit even when its mask is 0.
    using Key = std::pair<std::uint32_t, std::uint32_t>;
    using Entry = std::pair<Key, std::uint32_t>;  // ((min, max), mask)
    std::vector<Entry> constraint;
    for (std::size_t i = 0; i < edges.size() && i < masks.size(); ++i)
        constraint.emplace_back(std::minmax(edges[i][0], edges[i][1]), masks[i]);
    std::ranges::stable_sort(constraint, {}, &Entry::first);
    std::vector<std::uint8_t> bits(start.triangle_count(), 0);
    std::vector<std::array<std::uint32_t, 3>> edge_masks(start.triangle_count(), {0, 0, 0});
    for (std::size_t t = 0; t < start.triangle_count(); ++t)
        for (unsigned k = 0; k < 3; ++k) {
            const auto& tri = start.triangles()[t];
            const Key key = std::minmax(tri[k], tri[(k + 1) % 3]);
            if (const auto it = std::ranges::upper_bound(constraint, key, {}, &Entry::first);
                it != constraint.begin() && std::prev(it)->first == key) {
                bits[t] |= static_cast<std::uint8_t>(1u << k);
                edge_masks[t][k] = std::prev(it)->second;
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

// The next round's active set: every touched slot and every skipped one,
// ascending and once each. `skipped` is ascending (filled while walking the
// ascending `active`) and below touched.size(), so one linear merge with the
// slot walk gives what collecting, sorting and deduplicating gave (21a, QW3).
inline void rebuild_active(std::span<const char> touched, std::span<const std::uint32_t> skipped,
                           std::vector<std::uint32_t>& active) {
    active.clear();
    auto s = skipped.begin();
    for (std::uint32_t t = 0; t < touched.size(); ++t) {
        bool in = touched[t] != 0;
        for (; s != skipped.end() && *s == t; ++s)
            in = true;
        if (in)
            active.push_back(t);
    }
}

[[nodiscard]] inline bool needs_split(const ScanResult& r, double tolerance) noexcept {
    return r.node.has_value() && (r.is_void || r.max_error > tolerance);
}

[[nodiscard]] inline double foot_cap(const raster::RasterGeometry& g) noexcept {
    return std::min(g.delta_x(), g.delta_y()) / 2.0;
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
    const double cap = foot_cap(g), floor = std::min(dx, dy) / 100.0;
    return grad > 0.0 ? std::clamp(tol / grad, floor, cap) : cap;
}

// A Hit of `found` that may go in: a z there, else counted in `refused`.
// Anything but a Hit is nothing; NotCounterClockwise counts as refused (R5).
[[nodiscard]] inline std::optional<mesh::FootSearch> usable(const mesh::FootSearch& found, bool has_z,
                                                            std::size_t& refused) {
    if (found.status == mesh::FootStatus::Hit && has_z)
        return found;
    refused += found.status == mesh::FootStatus::NotCounterClockwise || found.status == mesh::FootStatus::Hit ? 1 : 0;
    return std::nullopt;
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
    m.set_frozen_mask(options.frozen_mask);

    const auto frame = mesh::lattice_frame(g.delta_x(), g.delta_y(), g.rows(), g.cols());
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
        // The pass never inserts a NoData node or a foot without a z: trim would remove it.
        const auto q = mesh::improve<pred::DefaultKernel>(
            m, frame,
            mesh::QualityOptions{options.min_angle_deg, g.rows(), g.cols(), options.constraint_feet},
            [&](const mesh::MeshVertex& v) { return vertex_z(dem, v).has_value(); });
        out.quality_inserted = q.inserted;
        out.quality_feet = q.feet;
        out.quality_skipped = q.skipped_floor + q.skipped_outside + q.skipped_vertex
                            + q.skipped_blocked + q.walk_bound_hits + q.skipped_frozen
                            + q.skipped_void + q.skipped_near_line;
        out.quality_seconds = since(t0);
    }
    std::vector<ScanResult> results;
    std::set<std::pair<std::uint32_t, std::uint32_t>> footed;  // (row, col), R2 step 5
    const double foot_cap = detail::foot_cap(g);
    mesh::FlipStack flip_stack;  // one buffer for every legalise_around
    std::vector<std::uint32_t> active(m.triangle_count());
    for (std::uint32_t t = 0; t < active.size(); ++t)
        active[t] = t;

    while (true) {
        ++out.rounds;
        results.resize(m.triangle_count());
        t0 = clock::now();
        parallel_util::for_each_block(active.size(), options.threads, parallel_util::BlockSchedule{},
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
            std::uint32_t owner = t;
            std::size_t refused = 0;
            std::optional<mesh::FootSearch> foot;
            // eps <= cap, so a search to cap that finds nothing is the eps search's answer;
            // eps (DEM reads) is computed only when an edge is within cap.
            auto s = options.constraint_feet && !r.is_void && mesh::foot_reachable(m, t)
                         ? mesh::constraint_foot(m, t, p, foot_cap, frame)
                         : mesh::FootSearch{};
            if (s.status != mesh::FootStatus::None && !footed.contains({r.node->row, r.node->col})) {
                if (const double eps = detail::foot_epsilon(dem, *r.node, options.tolerance); eps < foot_cap)
                    s = mesh::constraint_foot(m, t, p, eps, frame);
                foot = detail::usable(s, s.status == mesh::FootStatus::Hit && vertex_z(dem, s.at), refused);
                if (foot)
                    std::tie(owner, edge, p) = std::tuple{foot->owner, foot->edge, foot->at};
            }
            if (!edge) {
                q = m.split_inside(t, *r.node);
            } else {
                const auto e = *edge;
                const std::uint32_t u = m.neighbours(owner)[e];
                if (touched[owner] != 0 || (u != mesh::kNoNeighbour && touched[u] != 0)) {
                    skipped.push_back(t);
                    continue;
                }
                q = m.split_edge(owner, e, p);
                seeds[0] = owner;
                n_seeds = u != mesh::kNoNeighbour ? 4 : 2;
                if (n_seeds == 4)
                    seeds[3] = u;
            }
            if (owner != t)  // R3: t is unchanged and unconverged; rescan it
                skipped.push_back(t);
            touched.resize(m.triangle_count(), 1);
            touched[owner] = 1;
            if (n_seeds == 4)
                touched[seeds[3]] = 1;
            out.flips += mesh::legalise_around<pred::DefaultKernel>(
                m, q, std::span<const std::uint32_t>{seeds.data(), n_seeds}, frame, flip_stack,
                [&](std::uint32_t s) { touched[s] = 1; });
            ++out.inserted;
            out.carved += r.is_void ? 1 : 0;
            out.feet += foot ? 1 : 0;
            out.feet_refused += refused;
            if (foot)
                footed.insert({r.node->row, r.node->col});
        }
        out.split_seconds += since(t0);
        if (!any)
            break;
        detail::rebuild_active(touched, skipped, active);
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
