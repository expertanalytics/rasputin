#pragma once

// A triangle's DEM node set as row spans (docs/increments/18-row-span-scan.md,
// R1 to R1c): for each lattice row the triangle's closed interior meets, the
// inclusive column interval [c0, c1] of its nodes, vertices excluded. The union
// of the spans is exactly the set a box walk finds with three exact orientation
// tests per node, so a caller walks contiguous memory instead of testing.
//
// Per row, edge k from a to b bounds the column: orient(a, b, (c, r)) is affine
// in c with slope b.row - a.row, so its exact sign is monotone along the row. A
// rising edge gives a lower bound, a falling one an upper bound, and a
// horizontal one either empties the row, leaves it alone, or lies on it
// (flat_edge). Node-only triangles take the bound by exact integer division
// (R1a). Otherwise a double estimate is corrected with the exact predicate,
// normally one step each way (R1b). A node vertex on the row is always an end
// of its interval and is removed there (R1c).
//
// Pure lattice geometry: it sees three vertices and no raster or frame.

#include <terrain/mesh/lattice_mesh.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>

namespace terrain::mesh {

struct RowSpan {
    std::uint32_t row;
    std::uint32_t c0, c1;   // inclusive, c0 <= c1; empty rows are never reported
    std::int8_t flat_edge;  // k when edge k is horizontal and lies on this row, else -1
};

// Floor and ceiling of n / d for any signs; C++'s / truncates toward zero.
[[nodiscard]] constexpr std::int64_t floor_div(std::int64_t n, std::int64_t d) noexcept {
    const std::int64_t q = n / d;
    return (n % d != 0 && ((n % d < 0) != (d < 0))) ? q - 1 : q;
}
[[nodiscard]] constexpr std::int64_t ceil_div(std::int64_t n, std::int64_t d) noexcept {
    return -floor_div(-n, d);
}

// Calls f(RowSpan) once per non-empty row, rows ascending. v is counter-
// clockwise. Every product below is of two differences bounded by the grid's
// dimensions, and rows * cols < 2^48 for any DEM held in memory, so no int64
// overflows (R1a).
template <std::invocable<RowSpan> F>
void for_each_row_span(const std::array<MeshVertex, 3>& v, F&& f) {
    const bool nodes = v[0].is_node() && v[1].is_node() && v[2].is_node();
    const auto [rlo, rhi] = std::minmax({v[0].row, v[1].row, v[2].row});
    const auto [clo, chi] = std::minmax({v[0].col, v[1].col, v[2].col});
    const auto box_lo = static_cast<std::int64_t>(std::ceil(clo));
    const auto box_hi = static_cast<std::int64_t>(std::floor(chi));
    const auto row_lo = static_cast<std::uint32_t>(std::ceil(rlo));
    const auto row_hi = static_cast<std::uint32_t>(std::floor(rhi));

    for (std::uint32_t r = row_lo; r <= row_hi; ++r) {
        std::int64_t c0 = box_lo, c1 = box_hi;
        std::int8_t flat = -1;
        for (unsigned k = 0; k < 3; ++k) {
            const MeshVertex a = v[k], b = v[(k + 1) % 3];
            if (a.row == b.row) {
                // No column bound: orient is exactly (b.col - a.col)(a.row - r),
                // one sign for the whole row, on nodes or not.
                const int sign = ((b.col > a.col) - (b.col < a.col)) * ((a.row > r) - (a.row < r));
                if (sign < 0)
                    c1 = c0 - 1;
                else if (sign == 0)
                    flat = static_cast<std::int8_t>(k);
                continue;
            }
            const bool lower = b.row > a.row;
            std::int64_t c;
            if (nodes) {
                const auto ac = static_cast<std::int64_t>(a.col), ar = static_cast<std::int64_t>(a.row);
                const std::int64_t dr = static_cast<std::int64_t>(b.row) - ar;
                const std::int64_t n = (static_cast<std::int64_t>(b.col) - ac) * (r - ar);
                c = ac + (lower ? ceil_div(n, dr) : floor_div(n, dr));
            } else {
                // Estimate, clamped to the box widened by one (NaN to its low
                // end), then walked to the exact bound: s is non-decreasing in c
                // for a lower bound and non-increasing for an upper one.
                auto s = [&](std::int64_t col) {
                    return orient_sign(a, b, MeshVertex{static_cast<double>(col), static_cast<double>(r)});
                };
                const double x = a.col + (b.col - a.col) * (r - a.row) / (b.row - a.row);
                const double e = lower ? std::ceil(x) : std::floor(x);
                const double lo = static_cast<double>(box_lo - 1), hi = static_cast<double>(box_hi + 1);
                c = static_cast<std::int64_t>(e == e ? std::clamp(e, lo, hi) : lo);
                if (lower) {
                    while (c <= box_hi && s(c) < 0) ++c;
                    while (c - 1 >= box_lo && s(c - 1) >= 0) --c;
                } else {
                    while (c >= box_lo && s(c) < 0) --c;
                    while (c + 1 <= box_hi && s(c + 1) >= 0) ++c;
                }
            }
            if (lower)
                c0 = std::max(c0, c);
            else
                c1 = std::min(c1, c);
        }
        // A node vertex on this row is an end of its interval, never inside.
        for (const MeshVertex& p : v)
            if (p.is_node() && p.row == r) {
                const auto pc = static_cast<std::int64_t>(p.col);
                if (c0 == pc) ++c0;
                if (c1 == pc) --c1;
            }
        if (c0 <= c1)
            f(RowSpan{r, static_cast<std::uint32_t>(c0), static_cast<std::uint32_t>(c1), flat});
    }
}

}  // namespace terrain::mesh
