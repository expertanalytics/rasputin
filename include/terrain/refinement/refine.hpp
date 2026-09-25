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

#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/refinement/scan.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace terrain::refinement {

enum class RefineStatus : std::uint8_t { Ok, OffLattice, NotCounterClockwise, InvalidTolerance };

struct RefineOptions {
    double tolerance = 0.0;  // metres, finite and >= 0
    unsigned threads = 0;    // 0: hardware concurrency
};

struct RefineOutcome {
    RefineStatus status = RefineStatus::Ok;
    std::string message;  // empty on Ok

    std::vector<Point2> vertices;  // world, via RasterGeometry::node
    std::vector<double> z;         // the node's value; 0.0 where it is NoData
    std::vector<std::uint8_t> valid;
    std::vector<TriangleIndices> triangles;
    std::vector<std::array<std::uint32_t, 2>> edges;  // constraint edges, each once
    std::vector<std::uint32_t> masks;                 // one per edge

    std::size_t rounds = 0;     // scans run; 1 when nothing needed a split
    std::size_t inserted = 0;   // vertices added
    double max_error = 0.0;     // over triangles with three valid vertices
    std::size_t uncovered = 0;  // valid nodes still inside void triangles

    [[nodiscard]] bool ok() const noexcept { return status == RefineStatus::Ok; }
};

namespace detail {

[[nodiscard]] inline RefineOutcome refusal(RefineStatus status, std::string message) {
    RefineOutcome out;
    out.status = status;
    out.message = std::move(message);
    return out;
}

// The lattice mesh for `start`, or the refusal. A vertex is on the lattice iff
// RasterGeometry::node of its rounded (row, col) gives it back bit for bit.
[[nodiscard]] inline std::variant<mesh::LatticeMesh, RefineOutcome> to_lattice(
    const raster::RasterGeometry& g, const IndexedMesh2& start,
    std::span<const std::array<std::uint32_t, 2>> edges, std::span<const std::uint32_t> masks) {
    std::vector<mesh::LatticeVertex> lattice;
    lattice.reserve(start.vertices().size());
    for (std::size_t i = 0; i < start.vertices().size(); ++i) {
        const Point2 p = start.vertices()[i];
        const double col = std::round((p.x - g.x_min()) / g.delta_x());
        const double row = std::round((g.y_max() - p.y) / g.delta_y());
        const bool in_grid = col >= 0.0 && row >= 0.0 && col < static_cast<double>(g.cols())
                          && row < static_cast<double>(g.rows());
        const mesh::LatticeVertex v{in_grid ? static_cast<std::uint32_t>(row) : 0u,
                                    in_grid ? static_cast<std::uint32_t>(col) : 0u};
        const Point2 back = g.node({v.row, v.col});
        if (!in_grid || back.x != p.x || back.y != p.y)
            return refusal(RefineStatus::OffLattice,
                           "refine: start vertex " + std::to_string(i) + " is not a DEM node");
        lattice.push_back(v);
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

    RefineOutcome out;
    std::vector<ScanResult> results;
    std::vector<std::uint32_t> active(m.triangle_count());
    for (std::uint32_t t = 0; t < active.size(); ++t)
        active[t] = t;

    while (true) {
        ++out.rounds;
        results.resize(m.triangle_count());
        parallel_util::for_each_chunk(active.size(), options.threads,
                                      [&](std::size_t begin, std::size_t end) {
                                          for (std::size_t i = begin; i < end; ++i)
                                              results[active[i]] = scan(dem, m, active[i]);
                                      });

        // `touched` is every slot a split wrote this round. A marked triangle
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
            if (r.where == NodeLocation::Inside) {
                m.split_inside(t, *r.node);
            } else {
                const auto e = static_cast<unsigned>(r.where) - 1;
                const std::uint32_t u = m.neighbours(t)[e];
                if (u != mesh::kNoNeighbour && touched[u] != 0) {
                    skipped.push_back(t);
                    continue;
                }
                m.split_edge(t, e, *r.node);
                if (u != mesh::kNoNeighbour)
                    touched[u] = 1;
            }
            touched[t] = 1;
            touched.resize(m.triangle_count(), 1);
            ++out.inserted;
        }
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
    for (const mesh::LatticeVertex v : m.vertices()) {
        const raster::CellIndex c{v.row, v.col};
        const bool valid = !dem.is_nodata(c);
        out.vertices.push_back(g.node(c));
        out.z.push_back(valid ? static_cast<double>(dem.value_at(c)) : 0.0);
        out.valid.push_back(valid ? 1 : 0);
    }
    out.triangles.assign(m.triangles().begin(), m.triangles().end());
    std::tie(out.edges, out.masks) = m.constraint_edges();
    return out;
}

}  // namespace terrain::refinement
