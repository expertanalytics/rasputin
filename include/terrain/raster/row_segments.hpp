#pragma once

// One DEM row's columns [c0, c1] as contiguous segments, ascending
// (docs/increments/18-row-span-scan.md, R6). A single grid is one segment of
// row(i), chosen at compile time; a source that cannot hand out a whole row
// (a later multi-tile mosaic) provides row_segments(i, c0, c1, f) and is
// forwarded to. Segments partition [c0, c1] and share the source's one
// nodata() (Ola's C2), so a walker indexing by first_col + j visits the same
// cells in the same order however the row is cut.

#include <terrain/raster/raster.hpp>

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <span>

namespace terrain::raster {

template <typename T>
struct RowSegment {
    std::span<const T> values;  // cells first_col .. first_col + size - 1
    std::uint32_t first_col;    // global lattice column of values[0]
};

template <RasterSource R, std::invocable<RowSegment<typename R::value_type>> F>
void for_each_row_segment(const R& dem, std::size_t row, std::uint32_t c0, std::uint32_t c1,
                          F&& f) {
    if constexpr (requires { dem.row_segments(row, c0, c1, f); })
        dem.row_segments(row, c0, c1, f);
    else
        f(RowSegment<typename R::value_type>{dem.row(row).subspan(c0, c1 - c0 + 1), c0});
}

}  // namespace terrain::raster
