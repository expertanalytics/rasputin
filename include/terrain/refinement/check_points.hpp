#pragma once

// The check-point store of the final check (docs/increments/15c-geographic-dem.md,
// D4): points given as world (x, y) in the target CRS with a float z, filed by
// lattice cell of the target grid, and read back by cell row and column range.
//
// Numbers only: the grid is a RasterGeometry, and no CRS, path or raster value
// reaches this file.
//
// Filing. A point goes to cell (floor(row), floor(col)) of its fractional
// lattice position, the last cell for a point on the far edges. A point outside
// the node rectangle, or with a non-finite coordinate or z, is dropped and
// counted in `outside`.
//
// Packing, 16 bytes per point: the cell's column, the two offsets inside the
// cell as float, and z. The row is the bucket the point sits in. A float offset
// resolves 2^-24 of a cell. The check point IS the stored position: callers see
// only (cell + offset) as a double, so membership, error and an inserted vertex
// all use the same number. An offset that rounds up to a whole cell moves the
// point to the next cell's corner (offset 0), so every point lies in the closed
// cell it is filed in.
//
// freeze() sorts each row by (col, row offset, col offset, z), so iteration
// order does not depend on the order of add calls, and keeps the first of
// points at one stored position, counting the rest in `duplicates`. After it
// the store is immutable and any number of threads may read it.

#include <terrain/core/point.hpp>
#include <terrain/mesh/lattice_mesh.hpp>
#include <terrain/raster/geometry.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace terrain::refinement {

// The index of the last cell along an axis of n nodes (n >= 1).
[[nodiscard]] constexpr std::size_t last_cell(std::size_t n) noexcept { return std::max<std::size_t>(n, 2) - 2; }

class CheckPoints {
public:
    explicit CheckPoints(const raster::RasterGeometry& g) : geometry_{g}, rows_(last_cell(g.rows()) + 1) {}

    [[nodiscard]] const raster::RasterGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] bool frozen() const noexcept { return frozen_; }
    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] std::size_t duplicates() const noexcept { return duplicates_; }
    [[nodiscard]] std::size_t outside() const noexcept { return outside_; }

    void add(std::span<const Point2> xy, std::span<const float> z) {
        if (frozen_)
            throw std::logic_error("CheckPoints::add: the store is frozen");
        if (xy.size() != z.size())
            throw std::invalid_argument("CheckPoints::add: xy and z differ in length");
        const raster::RasterGeometry& g = geometry_;
        const auto max_col = static_cast<double>(g.cols() - 1), max_row = static_cast<double>(g.rows() - 1);
        for (std::size_t i = 0; i < xy.size(); ++i) {
            const double col = (xy[i].x - g.x_min()) / g.delta_x(), row = (g.y_max() - xy[i].y) / g.delta_y();
            // Written so that NaN fails every test: a NaN coordinate or z is outside.
            if (!(col >= 0.0 && col <= max_col && row >= 0.0 && row <= max_row && std::isfinite(z[i]))) {
                ++outside_;
                continue;
            }
            const auto [r, dr] = file(row, g.rows());
            const auto [c, dc] = file(col, g.cols());
            rows_[r].push_back(Packed{c, dc, dr, z[i]});
        }
    }

    void freeze() {
        if (frozen_)
            return;
        for (auto& row : rows_) {
            std::ranges::sort(row, {}, [](const Packed& p) { return std::tuple{p.col, p.dr, p.dc, p.z}; });
            const auto same = [](const Packed& a, const Packed& b) {
                return a.col == b.col && a.dr == b.dr && a.dc == b.dc;
            };
            const auto end = std::unique(row.begin(), row.end(), same);  // keeps the first
            duplicates_ += static_cast<std::size_t>(row.end() - end);
            row.erase(end, row.end());
            row.shrink_to_fit();
            size_ += row.size();
        }
        frozen_ = true;
    }

    // f(MeshVertex p, float z) for every point in cell row `row` whose cell
    // column is in [c0, c1] (inclusive), in stored order. p is the stored
    // fractional (col, row) position.
    template <class F>
    void for_each_in(std::size_t row, std::size_t c0, std::size_t c1, F&& f) const {
        if (row >= rows_.size())
            return;
        const auto& r = rows_[row];
        auto it = std::ranges::lower_bound(r, c0, {}, [](const Packed& p) { return std::size_t{p.col}; });
        for (; it != r.end() && it->col <= c1; ++it)
            f(mesh::MeshVertex{static_cast<double>(it->col) + static_cast<double>(it->dc),
                               static_cast<double>(row) + static_cast<double>(it->dr)},
              it->z);
    }

private:
    struct Packed {
        std::uint32_t col;
        float dc, dr;  // offsets inside the cell, [0, 1]
        float z;
    };
    static_assert(sizeof(Packed) == 16);

    // The cell along one axis of n nodes and the offset in it.
    static std::pair<std::uint32_t, float> file(double f, std::size_t n) {
        const std::size_t last = last_cell(n);
        auto c = std::min(static_cast<std::size_t>(f), last);
        auto o = static_cast<float>(f - static_cast<double>(c));
        if (o >= 1.0f && c < last) {
            ++c;
            o = 0.0f;
        }
        return {static_cast<std::uint32_t>(c), o};
    }

    raster::RasterGeometry geometry_;
    std::vector<std::vector<Packed>> rows_;  // one per cell row
    std::size_t size_ = 0, duplicates_ = 0, outside_ = 0;
    bool frozen_ = false;
};

}  // namespace terrain::refinement
