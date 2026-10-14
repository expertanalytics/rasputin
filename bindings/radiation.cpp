// The sun's position and the receiver check (include/terrain/radiation/) for
// Python: times cross as int64 UTC milliseconds since the Unix epoch, angles in
// degrees. No time zone or CRS crosses here; radiation/time.py turns aware
// times into UTC milliseconds, radiation/frame.py builds the receivers.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <terrain/radiation/receivers.hpp>
#include <terrain/radiation/solar_position.hpp>

#include <chrono>
#include <cstdint>
#include <span>
#include <string>

namespace py = pybind11;
using namespace terrain::radiation;

namespace {

utc_ms at(std::int64_t ms) { return utc_ms{std::chrono::milliseconds{ms}}; }

using Doubles = py::array_t<double, py::array::c_style | py::array::forcecast>;

std::span<const double> view(const Doubles& a) { return {a.data(), static_cast<std::size_t>(a.size())}; }

struct SunPositions {
    py::array_t<double> elevation, azimuth, distance;
};

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
            const ReceiverCheck c = check_receivers(ReceiverSet{view(x), view(y), view(z_obs), view(lat), view(lon),
                                                                view(jacobian), view(normal), view(exclusion_radius)});
            const std::string at = "receiver " + std::to_string(c.index);
            switch (c.status) {
            case ReceiverStatus::ok:
                return;
            case ReceiverStatus::size_mismatch:
                throw py::value_error("the receiver arrays differ in length: x, y, z_obs, lat_deg, lon_deg and "
                                      "exclusion_radius need n values, jacobian 4 n and normal 3 n");
            case ReceiverStatus::non_finite:
                throw py::value_error(at + " has a number that is not finite");
            case ReceiverStatus::polar:
                throw py::value_error(at + " is at a latitude beyond 89.9 degrees, where east and north are "
                                           "undefined");
            case ReceiverStatus::negative_exclusion_radius:
                throw py::value_error(at + " has a negative exclusion radius");
            }
        },
        py::arg("x"), py::arg("y"), py::arg("z_obs"), py::arg("lat_deg"), py::arg("lon_deg"), py::arg("jacobian"),
        py::arg("normal"), py::arg("exclusion_radius"), R"doc(
The core's check of a receiver set: None, or a ValueError naming the first
failing receiver (a number that is not finite, a latitude beyond 89.9
degrees, a negative exclusion radius) or arrays of unequal length. jacobian is
(n, 2, 2) and normal (n, 3), or the same numbers flat.
)doc");
}
