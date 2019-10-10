enable_testing()
find_package(Catch2 REQUIRED)
add_executable(foo test_sun_position.cpp)
target_link_libraries(foo Catch2::Catch2)
target_include_directories(foo PRIVATE "../src/rasputin")
include(CTest)
include(Catch)
catch_discover_tests(foo)

