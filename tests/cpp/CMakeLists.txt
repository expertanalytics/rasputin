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
# test_bbox names no kernel, so it registers with add_terrain_test. Everything
# else names DefaultKernel, which pulls in detria_exact.hpp and therefore needs
# the compiled terrain_predicates target linked -- add_terrain_backend_test.
add_terrain_test(test_bbox unit/test_bbox.cpp)
add_terrain_backend_test(test_segment unit/test_segment.cpp)
add_terrain_backend_test(test_ring unit/test_ring.cpp)
add_terrain_backend_test(prop_ring_invariants property/prop_ring_invariants.cpp)

# Core, increment 3: the PSLG.
#
# All three name a kernel -- the builder suite instantiates build<K>() and both
# others build under DefaultKernel -- so none of them uses the plain helper.
add_terrain_backend_test(test_pslg unit/test_pslg.cpp)
add_terrain_backend_test(test_pslg_builder unit/test_pslg_builder.cpp)
add_terrain_backend_test(prop_pslg_invariants property/prop_pslg_invariants.cpp)

# The concurrent-read property asserts guarantee 10 -- that a const Pslg needs
# no synchronisation -- by having eight threads reduce over the same one and
# requiring they agree.
target_link_libraries(prop_pslg_invariants PRIVATE Threads::Threads)

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

# CDT, increment 4: the constrained Delaunay triangulation.
#
# `terrain_cdt` is the compiled target these suites expect: a static library
# whose only translation unit is src/cdt/detria_backend.cpp, with lib/detria on
# its include path PRIVATE -- the same shape, for the same reason, as
# terrain_predicates. terrain_predicates is linked too because every suite here
# builds its Pslg under DefaultKernel.
function(add_terrain_cdt_test name)
    add_terrain_test(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE terrain_cdt terrain_predicates)
endfunction()

add_terrain_cdt_test(test_cdt_detria_backend unit/test_cdt_detria_backend.cpp)
add_terrain_cdt_test(test_cdt_constrained_edges unit/test_cdt_constrained_edges.cpp)
add_terrain_cdt_test(prop_cdt_invariants property/prop_cdt_invariants.cpp)

# test_cdt_indexed_mesh names no kernel and no backend: IndexedMesh2 is a value
# type in core/, so the plain helper is right for it.
add_terrain_test(test_cdt_indexed_mesh unit/test_cdt_indexed_mesh.cpp)

# THE SEAM SUITE USES THE PLAIN HELPER DELIBERATELY. THAT IS THE TEST, NOT AN
# OVERSIGHT. It links terrain_headers and Catch2 and nothing else, so the one
# symbol terrain_cdt contributes -- DetriaBackend::triangulate -- is absent from
# this link line, and a generic entry point that ever calls the shipped backend
# fails to LINK here. Registering it with add_terrain_cdt_test would leave every
# assertion passing and the proof gone. It builds its Pslg under the
# header-only FastKernel so it needs terrain_predicates no more than it needs
# terrain_cdt; see the comment at the top of the file.
add_terrain_test(test_cdt_backend_seam unit/test_cdt_backend_seam.cpp)
