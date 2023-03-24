#include <pybind11/pybind11.h>

#include "solar_position.h"

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
