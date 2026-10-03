// Increment 15f-2, L16 (docs/increments/15f-edge-strip.md): the coincidence
// radius of L12 and L14 scales with the lattice,
//
//   detail::coincidence_radius(const raster::RasterGeometry& g)
//       = max(1e-10, 64 * ulp(M)),  M = max(cols, rows) - 1,
//       ulp(M) = std::nextafter(M, inf) - M.
//
// The expected values are written by hand as powers of two, never computed
// by the formula under test, and the helper is read with no mesh or loop, so
// a wrong radius shows here before any refinement case runs.

#include <catch2/catch_test_macros.hpp>

#include <terrain/raster/geometry.hpp>
#include <terrain/refinement/strip_scan.hpp>

#include <cmath>
#include <cstddef>

using terrain::raster::RasterGeometry;
using terrain::refinement::detail::coincidence_radius;

namespace {

RasterGeometry lattice(std::size_t cols, std::size_t rows) { return RasterGeometry{0.0, 0.0, 2.0, 1.0, cols, rows}; }

}  // namespace

TEST_CASE("L16: coincidence_radius is the 1e-10 floor up to 8,191 nodes across, then 64 ulps of the extent",
          "[edge_strip][L16]") {
    SECTION("M = 8, the ES lattices: the floor, exactly") {
        CHECK(coincidence_radius(lattice(9, 9)) == 1e-10);
        CHECK(coincidence_radius(lattice(9, 2)) == 1e-10);
    }
    SECTION("M = 50,296, the whole basin as one window: 64 * 2^-37, whichever axis is longer") {
        CHECK(coincidence_radius(lattice(50297, 41332)) == std::ldexp(64.0, -37));
        CHECK(coincidence_radius(lattice(41332, 50297)) == std::ldexp(64.0, -37));
    }
    SECTION("M = 2^20, a corridor: 64 * 2^-32") {
        CHECK(coincidence_radius(lattice((std::size_t{1} << 20) + 1, 2)) == std::ldexp(64.0, -32));
    }
    SECTION("the floor's edge: M = 8,190 keeps 1e-10, M = 16,384 is 64 * 2^-38") {
        CHECK(coincidence_radius(lattice(8191, 9)) == 1e-10);
        CHECK(coincidence_radius(lattice(16385, 9)) == std::ldexp(64.0, -38));
    }
}
