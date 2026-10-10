#pragma once

// A vertical tolerance that follows the slope of the DEM at every node
// (docs/increments/34-slope-tolerance.md, sections 3, 4.3 and 4.5).
//
// SlopeRamp is section 3's t(s): N from END degrees up, F up to START, linear
// between; tested in that order, so a step (START = END) gives N from END up.
// SlopeTolerance holds every node's class (raster::steepness) and two
// 256-entry tables: allowed(c) = t(c / 2) for classes 0 to 180 and F for 181
// to 255 (255 is a NoData node, whose error is never measured), and
// weight(c) = 1 / allowed(c). Immutable after make, so any number of threads
// may read it.
//
// Sloped<P> is a per-triangle policy (33's TolerancePolicy) with the per-node
// slope beside it: lowest(), highest() and at() are the triangle part's, and
// the scans read `slope`. Every DEM node is held to allowed(its class), and a
// check point or strip point to allowed(cell_class(its position)): the largest
// class of the valid corners of the cell CheckPoints files it in, 0 when none
// is valid.

#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/mesh/row_spans.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/steepness.hpp>
#include <terrain/refinement/check_points.hpp>
#include <terrain/refinement/line_tolerance.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace terrain::refinement {

struct SlopeRamp {  // metres and degrees; 0 < near <= far, 0 <= start <= end < 90, all finite
    double near, far, start_deg, end_deg;

    [[nodiscard]] double at(double slope_deg) const noexcept {
        if (slope_deg >= end_deg)
            return near;
        if (slope_deg <= start_deg)
            return far;
        return far + (near - far) * (slope_deg - start_deg) / (end_deg - start_deg);
    }
};

class SlopeTolerance {
public:
    // nullopt with one plain sentence in `why` for a ramp value that breaks
    // its bound or is not finite. Computes steepness(dem, threads).
    template <raster::RasterSource R>
    [[nodiscard]] static std::optional<SlopeTolerance> make(const R& dem, SlopeRamp ramp, unsigned threads,
                                                            std::string& why) {
        const std::pair<const char*, double> values[] = {
            {"near", ramp.near}, {"far", ramp.far}, {"start", ramp.start_deg}, {"end", ramp.end_deg}};
        for (const auto& [name, v] : values)
            if (!std::isfinite(v)) {
                why = std::string{name} + " must be finite";
                return std::nullopt;
            }
        if (ramp.near <= 0.0)
            why = "near must be above 0";
        else if (ramp.near > ramp.far)
            why = "near must be at most far";
        else if (ramp.start_deg < 0.0)
            why = "start must be >= 0";
        else if (ramp.start_deg > ramp.end_deg)
            why = "start must be at most end";
        else if (ramp.end_deg >= 90.0)
            why = "end must be below 90 degrees";
        if (!why.empty())
            return std::nullopt;
        return SlopeTolerance{dem.geometry(), ramp, raster::steepness(dem, threads)};
    }

    [[nodiscard]] const raster::RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] std::span<const std::uint8_t> row(std::size_t r) const noexcept {
        return std::span<const std::uint8_t>{classes_}.subspan(r * geometry_.cols(), geometry_.cols());
    }
    // The largest valid corner class of the cell CheckPoints files p in,
    // (floor(row), floor(col)) clamped to the last cell; 0 if none is valid.
    [[nodiscard]] std::uint8_t cell_class(mesh::MeshVertex p) const noexcept {
        const std::size_t r0 = std::min(static_cast<std::size_t>(p.row), last_cell(geometry_.rows()));
        const std::size_t c0 = std::min(static_cast<std::size_t>(p.col), last_cell(geometry_.cols()));
        std::uint8_t best = 0;
        for (std::size_t r = r0; r <= r0 + 1 && r < geometry_.rows(); ++r)
            for (std::size_t c = c0; c <= c0 + 1 && c < geometry_.cols(); ++c)
                if (const std::uint8_t k = classes_[r * geometry_.cols() + c]; k != raster::kNoDataClass)
                    best = std::max(best, k);
        return best;
    }
    [[nodiscard]] double allowed(std::uint8_t c) const noexcept { return allowed_[c]; }
    [[nodiscard]] double weight(std::uint8_t c) const noexcept { return weight_[c]; }
    [[nodiscard]] std::array<std::size_t, 256> histogram() const {
        std::array<std::size_t, 256> h{};
        for (const std::uint8_t c : classes_)
            ++h[c];
        return h;
    }
    [[nodiscard]] double near() const noexcept { return ramp_.near; }
    [[nodiscard]] double far() const noexcept { return ramp_.far; }

private:
    SlopeTolerance(const raster::RasterGeometry& g, SlopeRamp ramp, std::vector<std::uint8_t> classes)
        : geometry_{g}, ramp_{ramp}, classes_{std::move(classes)} {
        for (std::size_t c = 0; c < allowed_.size(); ++c) {
            allowed_[c] = c < raster::kSteepnessClasses ? ramp.at(static_cast<double>(c) / 2.0) : ramp.far;
            weight_[c] = 1.0 / allowed_[c];
        }
    }

    raster::RasterGeometry geometry_;
    SlopeRamp ramp_;
    std::vector<std::uint8_t> classes_;
    std::array<double, 256> allowed_{}, weight_{};  // entry 255 (NoData): F
};

// A per-triangle policy with the per-node slope beside it (4.3). P may be a
// const reference, so a caller need not copy a LineTolerance.
template <TolerancePolicy P>
struct Sloped {
    P triangle;
    const SlopeTolerance* slope;
    [[nodiscard]] double lowest() const { return triangle.lowest(); }
    [[nodiscard]] double highest() const { return triangle.highest(); }
    [[nodiscard]] double at(const mesh::LatticeMesh& m, std::uint32_t t) const { return triangle.at(m, t); }
};

// A policy's triangle part (itself unless Sloped), and whether it is Sloped.
template <class P>
struct SlopeParts {
    using Triangle = P;
    static constexpr bool sloped = false;
};
template <class P>
struct SlopeParts<Sloped<P>> {
    using Triangle = std::remove_cvref_t<P>;
    static constexpr bool sloped = true;
};
template <class P>
inline constexpr bool is_sloped = SlopeParts<P>::sloped;

struct SlopeNodes {
    std::size_t valid = 0;      // valid DEM nodes inside the mesh
    std::size_t tightened = 0;  // of them, held to less than F by their slope
};

// Section 4.5: the valid nodes in the closed node sets of m's triangles, each
// once: a node on an edge by the lower-index triangle of the two across it, a
// vertex from the vertex list. Only a span's ends, and a row through a vertex
// (where an edge may be horizontal), can hold a vertex or an edge node.
// Integer sums per fixed block of slots, so any thread count gives the same.
[[nodiscard]] inline SlopeNodes count_slope_nodes(const mesh::LatticeMesh& m, const SlopeTolerance& s,
                                                  unsigned threads) {
    constexpr std::size_t b = 4096;
    std::vector<SlopeNodes> part(m.triangle_count() / b + 1);
    const auto add = [&](SlopeNodes& to, std::uint8_t c) {
        to.valid += c != raster::kNoDataClass ? 1 : 0;
        to.tightened += c != raster::kNoDataClass && s.allowed(c) < s.far() ? 1 : 0;
    };
    parallel_util::for_each_block(m.triangle_count(), threads, parallel_util::BlockSchedule{b},
                                  [&](std::size_t begin, std::size_t end) {
        SlopeNodes& to = part[begin / b];
        for (auto t = static_cast<std::uint32_t>(begin); t < end; ++t) {
            const std::array<mesh::MeshVertex, 3> v{m.corner(t, 0), m.corner(t, 1), m.corner(t, 2)};
            mesh::for_each_row_span(v, [&](mesh::RowSpan span) {
                const auto row = s.row(span.row);
                const bool vertex_row = v[0].row == span.row || v[1].row == span.row || v[2].row == span.row;
                for (std::uint32_t c = span.c0; c <= span.c1; ++c) {
                    if (vertex_row || c == span.c0 || c == span.c1) {
                        const mesh::MeshVertex p{mesh::LatticeVertex{span.row, c}};
                        if (p == v[0] || p == v[1] || p == v[2])
                            continue;
                        unsigned e = 0;
                        while (e < 3 && mesh::orient_sign(v[e], v[(e + 1) % 3], p) != 0)
                            ++e;
                        if (e < 3 && m.neighbours(t)[e] != mesh::kNoNeighbour && m.neighbours(t)[e] < t)
                            continue;
                    }
                    add(to, row[c]);
                }
            });
        }
    });
    SlopeNodes total;
    for (const mesh::MeshVertex v : m.vertices())
        if (const auto n = v.as_node())
            add(total, s.row(n->row)[n->col]);
    for (const SlopeNodes& p : part) {
        total.valid += p.valid;
        total.tightened += p.tightened;
    }
    return total;
}

}  // namespace terrain::refinement
