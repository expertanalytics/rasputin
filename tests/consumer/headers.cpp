// Every promised header of the installed core, compiled with only the
// imported target's settings: a RasterView<float> over a 2 x 2 array, sampled
// at a node and at the cell's centre. Exit 0 when both values are right.

#include <terrain/core/point.hpp>
#include <terrain/parallel_util/chunks.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/raster.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/raster/view.hpp>

#include <array>
#include <cstdio>
#include <optional>

static_assert(__cplusplus >= 202002L, "rasputin::terrain_headers must ask for C++20");

#ifdef CONSUMER_EXPECT_HARDENING
#if defined(_LIBCPP_VERSION)
#if !defined(_LIBCPP_HARDENING_MODE) || _LIBCPP_HARDENING_MODE != _LIBCPP_HARDENING_MODE_FAST
#error "rasputin::terrain_headers must carry the core's libc++ hardening mode"
#endif
#elif defined(__GLIBCXX__)
#ifndef _GLIBCXX_ASSERTIONS
#error "rasputin::terrain_headers must carry the core's _GLIBCXX_ASSERTIONS"
#endif
#endif
#endif

int main() {
    // Row 0 at y_max = 20, rows 10 apart; columns 10 apart from x_min = 0.
    const std::array<float, 4> nodes{1.0F, 2.0F, 3.0F, 4.0F};
    const terrain::raster::RasterGeometry geometry{0.0, 20.0, 10.0, 10.0, 2, 2};
    const terrain::raster::RasterView<float> view{geometry, nodes.data(), std::nullopt};

    const std::optional<double> node = terrain::raster::bilinear(view, terrain::Point2{10.0, 20.0});
    const std::optional<double> centre = terrain::raster::bilinear(view, terrain::Point2{5.0, 15.0});
    // Exact in binary: a node reads itself, the centre is the mean 2.5.
    if (node != 2.0 || centre != 2.5) {
        std::fprintf(stderr, "consumer: sampled %g and %g, expected 2 and 2.5\n",
                     node.value_or(-1.0), centre.value_or(-1.0));
        return 1;
    }
    return 0;
}
