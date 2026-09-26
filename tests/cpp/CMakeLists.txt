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

# Noding, increment 5a: the noder's numeric half -- the snap grid, the pairwise
# intersection classification and construction, the hot-pixel predicate and the
# node set.
#
# THE LINK LINES ARE PART OF THE SEAM BETWEEN 5a AND 5b, and are checkable
# rather than argued: no target here links terrain_cdt, and no 5a header
# includes core/pslg.hpp. 5a is everything whose correctness is a question about
# numbers; the moment one of these suites needs a Chain, the split has been
# violated and the build says so.
#
# Registration follows docs/increments/05-noder.md, "Files and LOC":
# test_noding_node_set names no kernel and uses the plain helper; the other
# three use add_terrain_backend_test, which links the compiled
# terrain_predicates target that DefaultKernel needs.
add_terrain_backend_test(test_noding_snap_rounding unit/test_noding_snap_rounding.cpp)
add_terrain_backend_test(test_noding_pairwise_intersection unit/test_noding_pairwise_intersection.cpp)
add_terrain_backend_test(prop_noding_snap_invariants property/prop_noding_snap_invariants.cpp)

add_terrain_test(test_noding_node_set unit/test_noding_node_set.cpp)

# Edge properties, increment 7: the opaque per-edge property set.
#
# The plain helper, deliberately: EdgeProperties names no kernel, includes no
# predicate header and is constexpr throughout, so it links terrain_headers and
# Catch2 and nothing else. test_cdt_indexed_mesh is registered the same way for
# the same reason -- a value type in core/ needs neither terrain_predicates nor
# terrain_cdt, and a link line that pulled either in would hide a dependency
# this header must not grow.
add_terrain_test(test_edge_properties unit/test_edge_properties.cpp)

# Noding, increment 5b: the noder's topology half -- the output type, the broad
# phase, the verifier and the driver.
#
# FIVE SUITES AGAINST FOUR HEADERS, which is Gate B's mapping plus its one
# declared relaxation (docs/increments/05b-noder-driver.md, "What is worth
# testing"): node.hpp gets a property file and a fixture file, split along a
# different axis than the one Gate B protects. Everything else is one test file
# per production header.
#
# No target here links terrain_cdt. Nothing in 5b knows the CDT exists; 5c is
# where that changes, and the link line is where it will show.
add_terrain_backend_test(test_noding_noded_pslg unit/test_noding_noded_pslg.cpp)

# The concurrent-read case in the NodedPslg suite is 5b's analogue of guarantee
# 10: a const NodedPslg needs no synchronisation, which is what lets 5c release
# the GIL around node() and read the result from whichever thread it lands on.
target_link_libraries(test_noding_noded_pslg PRIVATE Threads::Threads)

# The verifier. Invariant-critical: the only suite in 5b that can fail when the
# VERIFIER is broken, which is the whole reason NodedPslgBuilder is a separate
# type rather than a pass inside node<K>'s loop.
add_terrain_backend_test(test_noding_noded_pslg_builder unit/test_noding_noded_pslg_builder.cpp)

# The broad phase. Invariant-critical, and it USES THE PLAIN HELPER
# DELIBERATELY, which is the same kind of claim the CDT seam suite above makes.
# broad_phase.hpp is kernel-free -- it compares bounding boxes and never asks an
# orientation -- so terrain_predicates is absent from this link line and a broad
# phase that ever reaches for a predicate fails to LINK here rather than passing
# quietly. Registering it with add_terrain_backend_test would leave every
# assertion passing and the proof gone.
add_terrain_test(prop_noding_broad_phase property/prop_noding_broad_phase.cpp)

# The driver's property suite. Invariant-critical. The file name is the one
# testing.md's layout block reserves for 5b; `grep -n prop_noding_no_crossings
# testing.md` resolves it.
add_terrain_backend_test(prop_noding_no_crossings property/prop_noding_no_crossings.cpp)

# node.hpp's SECOND file: the driver's hand-built fixtures. This is Gate B's one
# declared relaxation -- a property file and a fixture file for the same header
# split along a different axis than the one Gate B protects, and both would move
# to a 5b-ii together if a seam is ever cut.
add_terrain_backend_test(test_noding_node unit/test_noding_node.cpp)

# Raster view, increment 12: RasterView<T>, row access on RasterSource, and the
# batch sampler (docs/increments/12-dem-to-mesh.md, R2 and R4).
add_terrain_test(test_raster_view unit/test_raster_view.cpp)

# Adaptive refinement, increment 14 (docs/increments/14-adaptive-refinement.md,
# "Tests for @tester"). test_mesh_lattice_split and test_refinement_scan are the
# invariant-critical pair. All three link Threads: the loop and chunk suites
# start threads through for_each_chunk. The pair are backend tests because
# lattice_mesh.hpp includes default_kernel.hpp (increment 16's orient_sign),
# which needs terrain_predicates.
add_terrain_test(test_refinement_chunks unit/test_refinement_chunks.cpp)
add_terrain_backend_test(test_mesh_lattice_split unit/test_mesh_lattice_split.cpp)
add_terrain_backend_test(test_refinement_scan unit/test_refinement_scan.cpp)
target_link_libraries(test_refinement_chunks PRIVATE Threads::Threads)
target_link_libraries(test_refinement_scan PRIVATE Threads::Threads)
target_link_libraries(test_mesh_lattice_split PRIVATE Threads::Threads)

# Delaunay insertion, increment 14b (docs/increments/14b-delaunay-insertion.md,
# "Tests for @tester"). test_mesh_lawson and prop_refinement_refine's T3 oracle
# are the invariant-critical pair; both need terrain_predicates (DefaultKernel).
add_terrain_backend_test(test_mesh_lawson unit/test_mesh_lawson.cpp)
target_link_libraries(test_mesh_lawson PRIVATE Threads::Threads)

# Off-node start vertices, increment 16 (docs/increments/16-domain-polygon.md,
# "Tests for @tester"). test_refinement_scan_offnode and prop_refinement_refine
# are the invariant-critical pair.
add_terrain_backend_test(test_refinement_scan_offnode unit/test_refinement_scan_offnode.cpp)
target_link_libraries(test_refinement_scan_offnode PRIVATE Threads::Threads)
add_terrain_backend_test(prop_refinement_refine property/prop_refinement_refine.cpp)
target_link_libraries(prop_refinement_refine PRIVATE Threads::Threads)

# Mesh statistics, increment 17 (docs/increments/17-mesh-stats.md, "Tests for
# @tester"): refine's phase seconds and carving counter. Not invariant-critical.
add_terrain_backend_test(test_refinement_refine unit/test_refinement_refine.cpp)
target_link_libraries(test_refinement_refine PRIVATE Threads::Threads)

# The row-span scan, increment 18 (docs/increments/18-row-span-scan.md, "Tests
# for @tester"). test_mesh_row_spans (T1) and prop_refinement_scan_equivalence
# (T2, T2-seg, C2's concept change) are the invariant-critical pair; both
# compare against the frozen box walk in support/scan_oracle.hpp.
add_terrain_backend_test(test_mesh_row_spans unit/test_mesh_row_spans.cpp)
add_terrain_backend_test(prop_refinement_scan_equivalence
                         property/prop_refinement_scan_equivalence.cpp)

# A quality start, increment 20 (docs/increments/20-start-quality.md, "Tests for
# @tester"). test_mesh_quality (Q1 to Q6, Q8, Q9) is the invariant-critical
# suite; prop_refinement_quality re-runs 14b's T3 and T6 with the pass on and
# holds Q7, Q8 and Q10 through refine. Both run in the TSan job.
add_terrain_backend_test(test_mesh_quality unit/test_mesh_quality.cpp)
add_terrain_backend_test(prop_refinement_quality property/prop_refinement_quality.cpp)
target_link_libraries(test_mesh_quality PRIVATE Threads::Threads)
target_link_libraries(prop_refinement_quality PRIVATE Threads::Threads)
