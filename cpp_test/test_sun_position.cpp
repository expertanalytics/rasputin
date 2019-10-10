#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include <catch2/catch.hpp>
#include <solar_position.h>
#include <cmath>
//#include "date/date.h"
#include <chrono>
#include <ctime>

namespace rasputin::test_utils {

auto fixed_cal_delta_t_calc() {
    return [] (const unsigned int, const unsigned int, const unsigned int) {
        return 67.0;
    };
}

auto fixed_time_point_delta_t_calc() {
    return [] (std::chrono::system_clock::time_point) {
        return 67.0;
    };
}
} // namespace rasputin::test_utils

TEST_CASE("Reference example test", "[reference]") {
    const unsigned int year = 2003;
    const unsigned int month = 10;
    const double day = 17.0 + 12.0/24.0 + 30.0/(60.0*24) + 30.0/(60*60*24) + 7.0/24.0;
    const double lat = 39.742476;
    const double lon = -105.1786;
    const double masl = 1830.14;
    const double P = 820;
    const double T = 11;
    using namespace rasputin::test_utils;
    using namespace rasputin::solar_position;
    const auto [Phi, e0] = calendar_solar_position(year, month, day, lat, lon, masl,
                                                   collectors::azimuth_and_elevation(),
                                                   fixed_cal_delta_t_calc());
    REQUIRE(abs(Phi - 194.34024) < 1.0e-4);
    const auto [e, Theta] = corrected_solar_elevation(e0, P, T);
    REQUIRE(abs(Theta - 50.11162) < 1.0e-4);
}

TEST_CASE("JD test 1", "[jd1]") {
    const unsigned int year = 2000;
    const unsigned int month = 1;
    const double day = 1.5;
    using namespace rasputin::solar_position;
    REQUIRE(jd_from_cal(year, month, day) == 2451545.0);
}

TEST_CASE("UTC cal test 1", "[utccalc]") {
    using namespace rasputin::test_utils;
    using namespace rasputin::solar_position;

    const double lat = 39.742476;
    const double lon = -105.1786;
    const double masl = 1830.14;
    const double P = 820;
    const double T = 11;
    const unsigned int year = 2000;
    const unsigned int mon = 1;
    const unsigned int dom = 1;
    const unsigned int h = 19;

    using namespace std::chrono;
    //using namespace date;
    const auto [Phi0, elevation0] = calendar_solar_position(year, mon, dom + h/24.0, lat, lon, masl, 
                                                            collectors::azimuth_and_elevation(), 
                                                            fixed_cal_delta_t_calc());
    const auto tp = sys_days{month(mon)/dom/year} + hours(h);
    const auto [Phi1, elevation1] = time_point_solar_position(tp, lat, lon, masl, 
                                                              collectors::azimuth_and_elevation(),
                                                              fixed_time_point_delta_t_calc());
    REQUIRE(elevation0 == elevation1);
    REQUIRE(Phi0 == Phi1);
}

TEST_CASE("DeltaT test", "[DT]") {
    using namespace rasputin::test_utils;

    const double lat = 39.742476;
    const double lon = -105.1786;
    const double masl = 1830.14;
    const double P = 820;
    const double T = 11;
    const unsigned int year = 2000;
    const unsigned int mon = 1;
    const unsigned int dom = 1;
    const unsigned int h = 19;

    using namespace std::chrono;
    //using namespace date;
    using namespace rasputin::solar_position;
    auto tp = sys_days{month(mon)/dom/year};
    const auto result = time_point_solar_position(tp, lat, lon, masl, 
                                                  collectors::azimuth_and_elevation(),
                                                  delta_t_calculator::coarse_timestamp_calc());
}

