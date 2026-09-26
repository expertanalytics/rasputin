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

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/lawson.hpp>
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

    // Wall seconds on the calling thread, steady_clock (17-mesh-stats.md R6).
    // for_each_chunk joins its workers before returning, so no worker reads a
    // clock. Not part of the determinism guarantee.
    double legalise_seconds = 0.0;  // the start mesh's legalise_all
    double scan_seconds = 0.0;      // every round's parallel scan, summed
    double split_seconds = 0.0;     // every round's serial split + flip phase, summed

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
    std::vector<ScanResult> results;
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
            if (r.where == NodeLocation::Inside) {
                q = m.split_inside(t, *r.node);
            } else {
                const auto e = static_cast<unsigned>(r.where) - 1;
                const std::uint32_t u = m.neighbours(t)[e];
                if (u != mesh::kNoNeighbour && touched[u] != 0) {
                    skipped.push_back(t);
                    continue;
                }
                q = m.split_edge(t, e, *r.node);
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
        const Point2 p = i < start.vertices().size() ? start.vertices()[i] : g.node(c);
        const std::optional<double> z =
            node && g.node(c) == p ? vertex_z(dem, v) : raster::bilinear(dem, p);
        out.vertices.push_back(p);
        out.z.push_back(z.value_or(0.0));
        out.valid.push_back(z ? 1 : 0);
    }
    out.triangles.assign(m.triangles().begin(), m.triangles().end());
    std::tie(out.edges, out.masks) = m.constraint_edges();
    return out;
}

}  // namespace terrain::refinement
