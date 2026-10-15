// The installed binding header in a module of its own, as the private package
// uses it: its own tag, so its raster view is a distinct C++ type from the
// core's, bound module-local.

#include <rasputin_bindings/numpy_raster.hpp>

#include <pybind11/pybind11.h>

#include <terrain/raster/view.hpp>

#include <type_traits>
#include <variant>

namespace py = pybind11;

struct ConsumerRasterTag;
using ConsumerRasterView = TaggedRasterView<ConsumerRasterTag>;

// The record the core binds today, with the tag parameter unused in its body.
static_assert(std::is_same_v<decltype(ConsumerRasterView::array), py::object>);
static_assert(std::is_same_v<decltype(ConsumerRasterView::view),
                             std::variant<terrain::raster::RasterView<float>,
                                          terrain::raster::RasterView<double>>>);
static_assert(!std::is_same_v<ConsumerRasterView, TaggedRasterView<struct OtherTag>>);

PYBIND11_MODULE(consumer_ext, m) {
    py::class_<ConsumerRasterView>(m, "RasterView", py::module_local());
}
