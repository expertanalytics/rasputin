#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "solar_position.h"
#include "shade.h"

namespace py = pybind11;

using namespace std::chrono;

void bind_solars(py::module& m) {
     m.def("timestamp_solar_position", [] (const double timestamp, const double geographic_latitude, const double geographic_longitude, const double masl) {
            const auto tp = sys_days{1970y /January / 1} + seconds(long(std::round(timestamp)));
            return rasputin::solar_position::time_point_solar_position(tp, geographic_latitude, geographic_longitude, masl,
                                                                     rasputin::solar_position::collectors::azimuth_and_elevation(),
                                                                     rasputin::solar_position::delta_t_calculator::coarse_timestamp_calc());
         }, "Compute azimuth and elevation of sun for given UTC timestamp.")
     .def("calendar_solar_position", [] (unsigned int year,unsigned int month, double day, const double geographic_latitude, const double geographic_longitude, const double masl) {
            return rasputin::solar_position::calendar_solar_position(year, month, day, geographic_latitude, geographic_longitude, masl,
                                                                     rasputin::solar_position::collectors::azimuth_and_elevation(),
                                                                     rasputin::solar_position::delta_t_calculator::coarse_date_calc());
         }, "Compute azimuth and elevation of sun for given UT calendar coordinate.")
     .def("solar_elevation_correction", &rasputin::solar_position::corrected_solar_elevation, "Correct elevation based on pressure and temperature.");
}

void bind_shades(py::module& m) {
    m.def("compute_shadow", (std::vector<int> (*)(const rasputin::Mesh &, const rasputin::point3 &))&rasputin::compute_shadow, "Compute shadows for given topocentric sun position.")
     .def("compute_shadow", (std::vector<int> (*)(const rasputin::Mesh &, const double, const double))&rasputin::compute_shadow, "Compute shadows for given azimuth and elevation.")
     .def("compute_shadows", &rasputin::compute_shadows, "Compute shadows for a series of times and ray directions.")
     .def("shade", [] (const rasputin::Mesh& mesh, const double timestamp) {
         const auto secs = seconds(int(std::round(timestamp)));
         const auto millisecs = milliseconds(int(round(1000*fmod(timestamp, 1))));
         const auto tp = sys_days{1970y / January / 1} + secs + millisecs;
         const auto shade_vec = rasputin::shade(mesh, tp);
         py::array_t<bool> result(shade_vec.size());
         auto info = result.request();
         auto *data = static_cast<bool*>(info.ptr);
         std::copy(std::begin(shade_vec), std::end(shade_vec), data);
         return result;
     });
}
