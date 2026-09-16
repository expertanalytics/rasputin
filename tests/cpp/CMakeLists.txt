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

# Test-only helpers: reference exact predicates and the deterministic
# adversarial point families. Not part of the shipped headers, so they live
# under tests/ and are visible only to test targets.
add_library(terrain_test_support INTERFACE)
target_include_directories(terrain_test_support INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/support"
)

# Every suite gets the same warning posture: the tests compile under exactly the
# flags the core does, so a warning introduced by a header is caught by whoever
# includes it first.
function(add_terrain_test name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE
        terrain_headers
        terrain_test_support
        Catch2::Catch2WithMain
    )
    target_compile_options(${name} PRIVATE
        -Wall -Wextra -Wpedantic -Werror
    )
    catch_discover_tests(${name})
endfunction()

add_terrain_test(test_point unit/test_point.cpp)
add_terrain_test(test_raster unit/test_raster.cpp)

# Predicates, increment 1a: the orientation vocabulary and the two kernels.
add_terrain_test(test_predicates_orientation unit/test_predicates_orientation.cpp)
add_terrain_test(test_predicates_fast_kernel unit/test_predicates_fast_kernel.cpp)
add_terrain_test(test_predicates_filtered_kernel unit/test_predicates_filtered_kernel.cpp)
add_terrain_test(prop_predicates_symmetry property/prop_predicates_symmetry.cpp)

# The concurrency tests in the filtered-kernel and detria-exact suites are the
# only things in the tree that start threads.
find_package(Threads REQUIRED)
target_link_libraries(test_predicates_filtered_kernel PRIVATE Threads::Threads)

# Predicates, increment 1b: the detria-backed exact backend.
#
# These do not build yet -- terrain/predicates/detria_exact.hpp,
# terrain/predicates/default_kernel.hpp and the terrain_predicates library that
# compiles src/predicates/detria_exact.cpp do not exist. That is the intended
# red step of the TDD protocol, not a broken checkout.
#
# `terrain_predicates` is the compiled target these suites expect: a small
# static library whose only translation unit is src/predicates/detria_exact.cpp,
# with lib/detria on its include path PRIVATE. Keeping detria's include
# directory off this target's INTERFACE is what makes the "detria.hpp must not
# escape one translation unit" rule enforceable rather than aspirational -- if
# detria_exact.hpp ever includes it, these suites fail to resolve the include.
# The preprocessor guard at the top of test_predicates_detria_exact.cpp covers
# the case where a future build change does put it on the path anyway.
function(add_terrain_backend_test name)
    add_terrain_test(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE terrain_predicates)
endfunction()

add_terrain_backend_test(test_predicates_detria_exact unit/test_predicates_detria_exact.cpp)
add_terrain_backend_test(test_predicates_default_kernel unit/test_predicates_default_kernel.cpp)
add_terrain_backend_test(prop_predicates_detria_agreement property/prop_predicates_detria_agreement.cpp)

# The concurrency test in the detria-exact suite starts threads, and the forked
# probe in the default-kernel suite needs no extra linkage on POSIX but is
# compiled out elsewhere.
target_link_libraries(test_predicates_detria_exact PRIVATE Threads::Threads)

# Core, increment 2: the 2D geometry value types.
#
# These do not build yet -- terrain/core/bbox.hpp, terrain/core/segment.hpp and
# terrain/core/ring.hpp do not exist. That is the intended red step of the TDD
# protocol, not a broken checkout.
#
# test_bbox names no kernel, so it registers with add_terrain_test. Everything
# else names DefaultKernel, which pulls in detria_exact.hpp and therefore needs
# the compiled terrain_predicates target linked -- add_terrain_backend_test.
add_terrain_test(test_bbox unit/test_bbox.cpp)
add_terrain_backend_test(test_segment unit/test_segment.cpp)
add_terrain_backend_test(test_ring unit/test_ring.cpp)
add_terrain_backend_test(prop_ring_invariants property/prop_ring_invariants.cpp)

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
