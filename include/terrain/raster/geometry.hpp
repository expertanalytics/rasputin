#pragma once

#include <terrain/core/point.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>

namespace terrain::raster {

// Row and column are named, not positional. The legacy code passed a bare
// two-element array between coordinates_to_indices and extract_buffer_values,
// and the two call sites disagreed about which slot was which -- a transposed
// read that only failed to crash on square rasters. Named members make that
// particular mistake unwriteable.
struct CellIndex {
    std::size_t row{};
    std::size_t col{};

    friend constexpr bool operator==(const CellIndex&, const CellIndex&) = default;
};

// North-up, grid-registered (pixel-is-point) raster geometry.
//
// Row 0 lies at y_max and y decreases as the row index grows; delta_y is
// stored positive. The extent spans (n - 1) spacings, not n: the samples are
// points on a grid, not areal cells. Both conventions are load-bearing --
// inverting the first flips every DEM, and changing the second shifts every
// extent by half a cell.
class RasterGeometry {
public:
    constexpr RasterGeometry(double x_min, double y_max, double delta_x, double delta_y,
                             std::size_t cols, std::size_t rows)
        : x_min_{x_min}, y_max_{y_max}, delta_x_{delta_x}, delta_y_{delta_y},
          cols_{cols}, rows_{rows} {
        if (!(delta_x > 0.0) || !(delta_y > 0.0))
            throw std::invalid_argument("RasterGeometry: delta_x and delta_y must be positive");
        if (cols == 0 || rows == 0)
            throw std::invalid_argument("RasterGeometry: cols and rows must both be non-zero");
    }

    [[nodiscard]] constexpr double x_min() const noexcept { return x_min_; }
    [[nodiscard]] constexpr double y_max() const noexcept { return y_max_; }
    [[nodiscard]] constexpr double delta_x() const noexcept { return delta_x_; }
    [[nodiscard]] constexpr double delta_y() const noexcept { return delta_y_; }
    [[nodiscard]] constexpr std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] constexpr std::size_t rows() const noexcept { return rows_; }

    [[nodiscard]] constexpr double x_max() const noexcept {
        return x_min_ + static_cast<double>(cols_ - 1) * delta_x_;
    }
    [[nodiscard]] constexpr double y_min() const noexcept {
        return y_max_ - static_cast<double>(rows_ - 1) * delta_y_;
    }

    [[nodiscard]] constexpr std::size_t size() const noexcept { return cols_ * rows_; }

    // Row-major. std::size_t throughout: the legacy used int and overflowed
    // past 2^31 cells.
    [[nodiscard]] constexpr std::size_t linear_index(const CellIndex& c) const noexcept {
        return c.row * cols_ + c.col;
    }

    [[nodiscard]] constexpr Point2 node(const CellIndex& c) const noexcept {
        return Point2{x_min_ + static_cast<double>(c.col) * delta_x_,
                      y_max_ - static_cast<double>(c.row) * delta_y_};
    }

    // Index of the node at or before p. Returns nullopt outside the domain --
    // the legacy clamped silently, turning every out-of-domain query into a
    // plausible-looking edge sample.
    //
    // The finiteness test must come first and must be explicit. Every
    // comparison against NaN is false, so a range check alone does not reject
    // it: control reaches the saturating path, which answers NaN with a
    // confident CellIndex{0,0}, and bilinear then returns an engaged optional
    // carrying NaN. Infinity happens to fall out of the range check correctly,
    // but relying on that is what hid the NaN hole.
    [[nodiscard]] std::optional<CellIndex> cell_of(const Point2& p) const noexcept {
        if (!std::isfinite(p.x) || !std::isfinite(p.y))
            return std::nullopt;
        if (p.x < x_min_ || p.x > x_max() || p.y > y_max_ || p.y < y_min())
            return std::nullopt;
        return clamped_cell_of(p);
    }

    // The legacy clamping behaviour, kept for callers that genuinely want a
    // saturating lookup -- but they now have to ask for it by name.
    //
    // Precondition: p is finite. This function has no failure channel, so a
    // non-finite coordinate yields an arbitrary in-range index. Callers that
    // cannot guarantee finiteness must use cell_of.
    [[nodiscard]] CellIndex clamped_cell_of(const Point2& p) const noexcept {
        return CellIndex{axis_index((y_max_ - p.y) / delta_y_, rows_),
                         axis_index((p.x - x_min_) / delta_x_, cols_)};
    }

    // Cell whose 2x2 node neighbourhood is fully in range, for interpolation.
    // This is the fix for the legacy out-of-bounds read: get_indices clamped
    // to n-1 and the interpolant then read index+1 unconditionally.
    [[nodiscard]] std::optional<CellIndex> bilinear_cell_of(const Point2& p) const noexcept {
        if (rows_ < 2 || cols_ < 2)
            return std::nullopt;
        auto c = cell_of(p);
        if (!c)
            return std::nullopt;
        c->row = std::min(c->row, rows_ - 2);
        c->col = std::min(c->col, cols_ - 2);
        return c;
    }

    // Strictly interior, with the legacy's relative epsilon. Boundary draping
    // relies on the strictness: an edge lying on the raster border must be
    // recognised and skipped, or it is re-emitted as a duplicate constraint
    // along the whole DEM boundary.
    [[nodiscard]] bool contains_strict(const Point2& p) const noexcept {
        if (!std::isfinite(p.x) || !std::isfinite(p.y))
            return false;
        const double eps = boundary_epsilon();
        return p.x > x_min_ + eps && p.x < x_max() - eps
            && p.y > y_min() + eps && p.y < y_max_ - eps;
    }

    // Tolerance for "lies on the boundary", scaled to coordinate magnitude
    // rather than to cell size alone.
    //
    // The legacy used sqrt(dx^2 + dy^2) * 1e-10. At UTM33 northings around
    // 7.9e6 with 10 m cells that is ~1.4e-9, while one ulp at that magnitude
    // is ~9.3e-10 -- so the tolerance was roughly 1.5 ulp, indistinguishable
    // from rounding noise at exactly the coordinates this project targets.
    // Boundary points arrive from clipping against this rectangle, so they
    // land within a few ulp of it; the tolerance has to sit well above that
    // and still far below one cell.
    [[nodiscard]] double boundary_epsilon() const noexcept {
        const double extent = std::max({std::abs(x_min_), std::abs(x_max()),
                                        std::abs(y_min()), std::abs(y_max_)});
        const double cell = std::sqrt(delta_x_ * delta_x_ + delta_y_ * delta_y_);
        return std::max(extent, cell) * 1e-12;
    }

private:
    [[nodiscard]] static constexpr std::size_t axis_index(double f, std::size_t n) noexcept {
        if (!(f > 0.0))
            return 0;
        const auto limit = static_cast<double>(n - 1);
        return f >= limit ? n - 1 : static_cast<std::size_t>(f);
    }

    double x_min_, y_max_, delta_x_, delta_y_;
    std::size_t cols_, rows_;
};

} // namespace terrain::raster
