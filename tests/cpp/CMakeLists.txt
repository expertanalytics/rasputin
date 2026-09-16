include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.6.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

add_executable(test_point test_point.cpp)
target_link_libraries(test_point PRIVATE
    terrain_headers
    Catch2::Catch2WithMain
)
target_compile_options(test_point PRIVATE
    -Wall -Wextra -Wpedantic -Werror
)

catch_discover_tests(test_point)

add_executable(test_raster test_raster.cpp)
target_link_libraries(test_raster PRIVATE
    terrain_headers
    Catch2::Catch2WithMain
)
target_compile_options(test_raster PRIVATE
    -Wall -Wextra -Wpedantic -Werror
)

catch_discover_tests(test_raster)

# Characterisation tests for the legacy solar-position header. These pin the
# UTC/Julian-day behaviour the post-CGAL port must reproduce, and guard removal
# of the lib/date fallback. The legacy header carries pre-existing
# unused-parameter warnings, so that one flag is waived -- -Werror stays on.
add_executable(test_solar_position test_solar_position.cpp)
target_include_directories(test_solar_position PRIVATE
    "${CMAKE_SOURCE_DIR}/legacy"
)
target_link_libraries(test_solar_position PRIVATE
    Catch2::Catch2WithMain
)
target_compile_options(test_solar_position PRIVATE
    -Wall -Wextra -Wpedantic -Werror -Wno-unused-parameter
)

catch_discover_tests(test_solar_position)
