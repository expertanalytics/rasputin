// The HBV snow kernel (include/terrain/snow/hbv_snow.hpp) for Python: arrays
// in, read-only arrays out, the GIL released during a step. No file and no
// CRS cross here.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <terrain/snow/hbv_snow.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;
using namespace terrain::snow;

namespace {

using F64 = py::array_t<double, py::array::c_style | py::array::forcecast>;

F64 vector_arg(const py::object& values, const char* name) {
    auto a = F64::ensure(values);
    if (!a || a.ndim() != 1)
        throw py::value_error(std::format("{} must be a one-dimensional float64 array", name));
    return a;
}

std::span<const double> span_of(const F64& a) { return {a.data(), static_cast<std::size_t>(a.size())}; }

// A read-only array over the state's own buffer, which `owner` keeps alive.
py::array_t<double> view(const py::object& owner, std::span<const double> v,
                         std::vector<py::ssize_t> shape) {
    static constexpr double kEmpty[1]{};  // numpy copies when handed a null pointer
    std::vector<py::ssize_t> strides(shape.size(), sizeof(double));
    if (shape.size() == 2)
        strides[0] = shape[1] * static_cast<py::ssize_t>(sizeof(double));
    py::array_t<double> a{std::move(shape), std::move(strides), v.empty() ? kEmpty : v.data(), owner};
    py::detail::array_proxy(a.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
    return a;
}

const char* status_name(StepStatus s) {
    switch (s) {
        case StepStatus::ok: return "ok";
        case StepStatus::invalid_parameters: return "invalid_parameters";
        case StepStatus::size_mismatch: return "size_mismatch";
        case StepStatus::non_finite_forcing: return "non_finite_forcing";
        case StepStatus::negative_precipitation: return "negative_precipitation";
        case StepStatus::bad_step_length: return "bad_step_length";
    }
    return "unknown";
}

// The first triangle whose temperature, precipitation or anomaly is not finite.
std::size_t first_non_finite(const StepForcing& f) {
    for (std::size_t i = 0; i < f.temperature.size(); ++i)
        if (!std::isfinite(f.temperature[i]) || !std::isfinite(f.precipitation[i]) ||
            (!f.radiation_anomaly.empty() && !std::isfinite(f.radiation_anomaly[i])))
            return i;
    return f.temperature.size();
}

}  // namespace

void bind_snow(py::module_& m) {
    py::enum_<SunCorrection>(m, "SunCorrection", "How the radiation anomaly corrects melt.")
        .value("none", SunCorrection::none)
        .value("melt_factor", SunCorrection::melt_factor)
        .value("threshold", SunCorrection::threshold);

    py::class_<SnowDistribution>(m, "SnowDistribution", R"doc(
The relative snow depth of N equal-area bins within a triangle. factors must
be 1 to 64 finite, non-negative, ascending values of mean 1 within 1e-9, else
ValueError; they are kept divided by their mean.
)doc")
        .def(py::init([](const py::object& factors) {
                 const auto a = vector_arg(factors, "factors");
                 return SnowDistribution{span_of(a)};
             }),
             py::arg("factors"))
        .def_property_readonly(
            "factors",
            [](const SnowDistribution& d) {
                py::array_t<double> a(static_cast<py::ssize_t>(d.factors().size()));
                std::ranges::copy(d.factors(), a.mutable_data());
                return a;
            },
            "The factors, float64, mean 1.");

    py::class_<SnowState>(m, "SnowState", R"doc(
The snow pack of n triangles of N bins each: solid sp and liquid sw (n, N) in
mm, and per triangle swe, sca and the last step's outflow (n,). The arrays are
read-only views that change with each step.
)doc")
        .def_static("uniform", &SnowState::uniform, py::arg("triangles"), py::arg("distribution"),
                    py::arg("initial_swe_mm"),
                    "Every bin f_k * initial_swe_mm of ice, no water; a negative or "
                    "non-finite initial SWE is a ValueError.")
        .def_property_readonly("sp", [](const py::object& self) {
            const auto& s = self.cast<const SnowState&>();
            return view(self, s.sp(), {static_cast<py::ssize_t>(s.triangles()), static_cast<py::ssize_t>(s.bins())});
        })
        .def_property_readonly("sw", [](const py::object& self) {
            const auto& s = self.cast<const SnowState&>();
            return view(self, s.sw(), {static_cast<py::ssize_t>(s.triangles()), static_cast<py::ssize_t>(s.bins())});
        })
        .def_property_readonly("swe", [](const py::object& self) {
            const auto& s = self.cast<const SnowState&>();
            return view(self, s.swe(), {static_cast<py::ssize_t>(s.triangles())});
        })
        .def_property_readonly("sca", [](const py::object& self) {
            const auto& s = self.cast<const SnowState&>();
            return view(self, s.sca(), {static_cast<py::ssize_t>(s.triangles())});
        })
        .def_property_readonly("outflow", [](const py::object& self) {
            const auto& s = self.cast<const SnowState&>();
            return view(self, s.outflow(), {static_cast<py::ssize_t>(s.triangles())});
        })
        .def(
            "step",
            [](SnowState& s, double tx, double ts, double cx, double lw, double cfr, SunCorrection mode,
               double r, double d, const py::object& temperature, const py::object& precipitation,
               const py::object& radiation_anomaly, double dt_seconds, unsigned threads) {
                const auto t = vector_arg(temperature, "temperature");
                const auto p = vector_arg(precipitation, "precipitation");
                const auto di = vector_arg(radiation_anomaly, "radiation_anomaly");
                const StepForcing f{span_of(t), span_of(p), span_of(di), dt_seconds};
                StepStatus status;
                {
                    const py::gil_scoped_release unlocked;  // the arrays are held by locals
                    status = step({tx, ts, cx, lw, cfr}, {mode, r, d}, s, f, threads);
                }
                if (status == StepStatus::non_finite_forcing)
                    throw py::value_error(std::format(
                        "snow step refused, non_finite_forcing: a forcing value at triangle {} is not finite",
                        first_non_finite(f)));
                if (status != StepStatus::ok)
                    throw py::value_error(std::format("snow step refused: {}", status_name(status)));
            },
            py::kw_only(), py::arg("tx"), py::arg("ts"), py::arg("cx"), py::arg("lw"), py::arg("cfr"),
            py::arg("mode"), py::arg("r"), py::arg("d"), py::arg("temperature"), py::arg("precipitation"),
            py::arg("radiation_anomaly"), py::arg("dt_seconds"), py::arg("threads"), R"doc(
One HBV step of dt_seconds for every triangle. temperature (degC) and
precipitation (mm over the step) are of length n; radiation_anomaly (W/m2) of
length n under a correction and empty under none. A refused input is a
ValueError naming the status (and, for a non-finite value, the first triangle)
and leaves the state as it was. threads only splits the work (0: all cores).
Releases the GIL.
)doc");
}
