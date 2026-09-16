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

# Every suite gets the same warning posture: the tests compile under exactly the
# flags the core does, so a warning introduced by a header is caught by whoever
# includes it first.
function(add_terrain_test name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE
        terrain_headers
        Catch2::Catch2WithMain
    )
    target_compile_options(${name} PRIVATE
        -Wall -Wextra -Wpedantic -Werror
    )
    catch_discover_tests(${name})
endfunction()

add_terrain_test(test_point unit/test_point.cpp)
add_terrain_test(test_raster unit/test_raster.cpp)

# Characterisation tests for the legacy solar-position header. These pin the
# UTC/Julian-day behaviour the post-CGAL port must reproduce, and guard removal
# of the lib/date fallback. The legacy header carries pre-existing
# unused-parameter warnings, so that one flag is waived -- -Werror stays on.
# It also does not link terrain_headers, so it does not go through the helper.
add_executable(test_solar_position unit/test_solar_position.cpp)
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
