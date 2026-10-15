// The sun's position, the receiver check, the horizon search and the instant
// evaluator (include/terrain/radiation/) for Python: times cross as int64 UTC
// milliseconds since the Unix epoch, angles in degrees, heights as float32
// RasterViews. No time zone or CRS crosses here; radiation/time.py turns aware
// times into UTC milliseconds, radiation/frame.py builds the receivers.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "numpy_raster.hpp"

#include <terrain/radiation/evaluate.hpp>
#include <terrain/radiation/horizon.hpp>
#include <terrain/radiation/receivers.hpp>
#include <terrain/radiation/solar_position.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <optional>
#include <span>
#include <string>

namespace py = pybind11;

struct CoreRasterTag;  // the view _core binds (core.cpp)
using BoundRasterView = TaggedRasterView<CoreRasterTag>;
using namespace terrain::radiation;

namespace {

utc_ms at(std::int64_t ms) { return utc_ms{std::chrono::milliseconds{ms}}; }

using Doubles = py::array_t<double, py::array::c_style | py::array::forcecast>;

std::span<const double> view(const Doubles& a) { return {a.data(), static_cast<std::size_t>(a.size())}; }

// The float32 view inside a bound raster, or a TypeError naming float32.
const terrain::raster::RasterView<float>& float32(const BoundRasterView& r, const char* what) {
    if (const auto* v = std::get_if<terrain::raster::RasterView<float>>(&r.view))
        return *v;
    throw py::type_error(std::string{"horizon_search: the "} + what + " raster must be float32; convert it with "
                         "astype(numpy.float32)");
}

struct SunPositions {
    py::array_t<double> elevation, azimuth, distance;
};

struct InstantArrays {
    py::array lit, cos_i, b, sin_e;
};

// A (rows, cols) NumPy copy of one of the evaluator's time-major outputs.
template <class T>
py::array_t<T> table(const std::vector<T>& v, py::ssize_t rows, py::ssize_t cols) {
    py::array_t<T> a({rows, cols});
    std::copy(v.begin(), v.end(), a.mutable_data());
    return a;
}

}  // namespace

void bind_radiation(py::module_& m) {
    m.def("delta_t_seconds", [](std::int64_t t) { return delta_t_seconds(at(t)); }, py::arg("time_utc_ms"),
          "The default delta T (TT - UT) in seconds: Espenak and Meeus's polynomials, "
          "the long-term parabola outside 1986-2050.");

    py::class_<SolarPosition>(m, "SolarPosition", "e0_deg (no refraction), e_deg (refracted) and azimuth_deg.")
        .def_readonly("e0_deg", &SolarPosition::e0_deg)
        .def_readonly("e_deg", &SolarPosition::e_deg)
        .def_readonly("azimuth_deg", &SolarPosition::azimuth_deg);
    m.def(
        "solar_position",
        [](std::int64_t t, double lat, double lon, double h, double dt, double p, double temp) {
            return solar_position(at(t), lat, lon, h, dt, p, temp);
        },
        py::arg("time_utc_ms"), py::arg("lat_deg"), py::arg("lon_deg"), py::arg("height_m"), py::arg("delta_t_s"),
        py::arg("pressure_hpa"), py::arg("temperature_c"), R"doc(
NREL's SPA at one time and place: the elevation without refraction e0_deg,
with refraction e_deg (held at its -1.8 degree value below that), and the
azimuth from true north, clockwise.
)doc");

    py::class_<SunPositions>(m, "SunPositions", "elevation, azimuth (degrees) and distance (AU), float64 arrays.")
        .def_readonly("elevation", &SunPositions::elevation)
        .def_readonly("azimuth", &SunPositions::azimuth)
        .def_readonly("distance", &SunPositions::distance);
    m.def(
        "sun_positions",
        [](const py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>& times, double lat,
           double lon, double h) {
            const py::ssize_t n = times.ndim() == 1 ? times.shape(0) : -1;
            if (n < 0)
                throw py::value_error("times_utc_ms must be a one-dimensional int64 array");
            SunPositions out{py::array_t<double>(n), py::array_t<double>(n), py::array_t<double>(n)};
            const std::int64_t* t = times.data();
            double *e = out.elevation.mutable_data(), *a = out.azimuth.mutable_data(), *r = out.distance.mutable_data();
            {
                const py::gil_scoped_release unlocked;  // the arrays are held by locals
                for (py::ssize_t i = 0; i < n; ++i) {
                    const GeocentricSun g = geocentric_sun(at(t[i]), delta_t_seconds(at(t[i])));
                    const TopocentricSun s = topocentric_sun(g, lat, lon, h);
                    e[i] = s.e0_deg;
                    a[i] = s.azimuth_deg;
                    r[i] = g.r_au;
                }
            }
            return out;
        },
        py::arg("times_utc_ms"), py::arg("lat"), py::arg("lon"), py::arg("height_m"), R"doc(
The sun at one place for many UTC times (int64 milliseconds) at the default
delta T: the elevation without refraction and the azimuth from true north,
clockwise, in degrees, and the sun's distance in AU. Releases the GIL.
)doc");

    m.def(
        "check_receivers",
        [](const Doubles& x, const Doubles& y, const Doubles& z_obs, const Doubles& lat, const Doubles& lon,
           const Doubles& jacobian, const Doubles& normal, const Doubles& exclusion_radius) {
            require_receivers(ReceiverSet{view(x), view(y), view(z_obs), view(lat), view(lon), view(jacobian),
                                          view(normal), view(exclusion_radius)});
        },
        py::arg("x"), py::arg("y"), py::arg("z_obs"), py::arg("lat_deg"), py::arg("lon_deg"), py::arg("jacobian"),
        py::arg("normal"), py::arg("exclusion_radius"), R"doc(
The core's check of a receiver set: None, or a ValueError naming the first
failing receiver (a number that is not finite, a latitude beyond 89.9
degrees, a negative exclusion radius) or arrays of unequal length. jacobian is
(n, 2, 2) and normal (n, 3), or the same numbers flat.
)doc");

    m.def(
        "horizon_search",
        [](const Doubles& x, const Doubles& y, const Doubles& z_obs, const Doubles& lat, const Doubles& lon,
           const Doubles& jacobian, const Doubles& normal, const Doubles& exclusion_radius, const BoundRasterView& near,
           const BoundRasterView* far, std::size_t n_directions, double curvature, double d_near, double d_far,
           unsigned threads) {
            const auto& near_view = float32(near, "near");
            const auto far_view = far ? std::optional{float32(*far, "far")} : std::nullopt;
            const ReceiverSet rc{view(x), view(y), view(z_obs), view(lat), view(lon), view(jacobian), view(normal),
                                 view(exclusion_radius)};
            HorizonTable t;
            {
                const py::gil_scoped_release unlocked;  // the arrays are held by the arguments
                t = horizon_search(rc, near_view, far_view, {n_directions, curvature, d_near, d_far}, threads);
            }
            const auto n = static_cast<py::ssize_t>(rc.size()), N = static_cast<py::ssize_t>(n_directions);
            py::array_t<std::int16_t> horizon({n, N});
            py::array_t<std::uint8_t> reach({n, N});
            std::copy(t.horizon.begin(), t.horizon.end(), horizon.mutable_data());
            std::copy(t.reach.begin(), t.reach.end(), reach.mutable_data());
            return py::make_tuple(horizon, reach);
        },
        py::arg("x"), py::arg("y"), py::arg("z_obs"), py::arg("lat_deg"), py::arg("lon_deg"), py::arg("jacobian"),
        py::arg("normal"), py::arg("exclusion_radius"), py::arg("near"), py::arg("far").none(true), py::kw_only(),
        py::arg("n_directions"), py::arg("curvature"), py::arg("d_near"), py::arg("d_far"), py::arg("threads") = 0,
        R"doc(
The horizon of each receiver in n_directions true azimuths (direction k at
k 360 / n_directions), searched on the near raster to d_near and on the far
raster from there to d_far; with far None the search ends at d_near. Earth
curvature and refraction enter as curvature = (1 - k) / (2 R) per metre.
Returns (horizon, reach), both (n, n_directions): horizon int16 in 0.01
degrees, -9000 where no sample was found; reach uint8, the whole km searched
before the first sample off the raster or on NoData, 255 when none was.
Both rasters must be float32. Refusals are a ValueError naming the field or
the receiver. Releases the GIL.
)doc");

    py::class_<InstantArrays>(m, "InstantFactors",
                              "lit (uint8), cos_i, b and sin_e (float32), each (times, receivers).")
        .def_readonly("lit", &InstantArrays::lit)
        .def_readonly("cos_i", &InstantArrays::cos_i)
        .def_readonly("b", &InstantArrays::b)
        .def_readonly("sin_e", &InstantArrays::sin_e);
    m.def(
        "evaluate_instant",
        [](const Doubles& lat, const Doubles& lon, const Doubles& normal, const py::array& horizon,
           const py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>& times, double limb,
           double pressure, double temperature, std::size_t n_directions, unsigned threads) {
            if (!horizon.dtype().is(py::dtype::of<std::int16_t>()))
                throw py::type_error("evaluate_instant: horizon must be int16, in 0.01 degrees");
            const auto h = py::array_t<std::int16_t, py::array::c_style>::ensure(horizon);
            const LitReceivers r{view(lat), view(lon), view(normal),
                                 {h.data(), static_cast<std::size_t>(h.size())}, n_directions};
            const std::span<const std::int64_t> t{times.data(), static_cast<std::size_t>(times.size())};
            InstantFactors f;
            {
                const py::gil_scoped_release unlocked;  // the arrays are held by the arguments and h
                f = evaluate_instant(r, t, SunSettings{limb, pressure, temperature}, threads);
            }
            const auto nt = static_cast<py::ssize_t>(t.size()), n = static_cast<py::ssize_t>(lat.size());
            return InstantArrays{table(f.lit, nt, n), table(f.cos_i, nt, n), table(f.b, nt, n),
                                 table(f.sin_e, nt, n)};
        },
        py::arg("lat_deg"), py::arg("lon_deg"), py::arg("normal"), py::arg("horizon"), py::arg("times_utc_ms"),
        py::kw_only(), py::arg("limb_deg") = 0.0, py::arg("pressure_hpa") = 1010.0, py::arg("temperature_c") = 10.0,
        py::arg("n_directions") = 64, py::arg("threads") = 0, R"doc(
Per UTC time (int64 milliseconds) and receiver: lit, the sun above the
horizon at its azimuth (raised by limb_deg) and in front of the receiver's
plane; cos_i, the normal against the apparent sun, unclipped; b = lit
max(0, cos_i); sin_e, the sine of the apparent elevation. normal is (n, 3) in
true east, north, up; horizon (n, n_directions) int16 in 0.01 degrees. Each
output is (times, n), time-major. Refraction at pressure_hpa and
temperature_c; a pressure of 0 is none. A ValueError names the field or the
receiver; a horizon that is not int16 is a TypeError. Releases the GIL.
)doc");
}
