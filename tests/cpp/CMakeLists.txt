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
#
# Optional `COST <value>`: a ctest scheduling hint put on every Catch2 case of
# the target, so `ctest --parallel` starts the long ones first. Without it the
# suite is registered exactly as before. The helpers below forward ${ARGN}, so
# they accept it too.
function(add_terrain_test name)
    cmake_parse_arguments(PARSE_ARGV 1 arg "" "COST" "")
    add_executable(${name} ${arg_UNPARSED_ARGUMENTS})
    target_link_libraries(${name} PRIVATE
        terrain_headers
        terrain_test_support
        Catch2::Catch2WithMain
    )
    target_compile_options(${name} PRIVATE
        -Wall -Wextra -Wpedantic -Werror
    )
    if(DEFINED arg_COST)
        catch_discover_tests(${name} PROPERTIES COST ${arg_COST})
    else()
        catch_discover_tests(${name})
    endif()
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
# invariant-critical pair. Both link Threads but start none, so neither runs in
# the TSan job (h17 §4c, docs/increments/h17-ci-test-time.md). Both are backend
# tests because lattice_mesh.hpp includes default_kernel.hpp (increment 16's
# orient_sign), which needs terrain_predicates.
add_terrain_backend_test(test_mesh_lattice_split unit/test_mesh_lattice_split.cpp)
add_terrain_backend_test(test_refinement_scan unit/test_refinement_scan.cpp)
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
# holds Q7, Q8 and Q10 through refine. prop_refinement_quality starts threads
# through refine and runs in the TSan job; test_mesh_quality starts none, so it
# is not in it (h17 §4c).
add_terrain_backend_test(test_mesh_quality unit/test_mesh_quality.cpp)
add_terrain_backend_test(prop_refinement_quality property/prop_refinement_quality.cpp)
target_link_libraries(test_mesh_quality PRIVATE Threads::Threads)
target_link_libraries(prop_refinement_quality PRIVATE Threads::Threads)

# The pass skips NoData nodes (docs/increments/20-start-quality.md, "Fix: the
# start-quality pass skips NoData nodes", Q-V1): improve's validity callable
# and skipped_void, in a target of their own so test_mesh_quality's
# default-call guard builds without them. Not mutation-critical.
add_terrain_backend_test(test_mesh_quality_void unit/test_mesh_quality_void.cpp)
target_link_libraries(test_mesh_quality_void PRIVATE Threads::Threads)

# Constraint feet, increment 20b (docs/increments/20b-min-insertion-distance.md,
# "Tests for @tester"). test_refinement_constraint_feet (F1 to F5, plus 14b's
# T6 with the rule on) is the invariant-critical suite; it runs in the TSan job.
add_terrain_backend_test(test_refinement_constraint_feet
                         property/prop_refinement_constraint_feet.cpp)
target_link_libraries(test_refinement_constraint_feet PRIVATE Threads::Threads)

# Parallel refine quick wins, increment 21a (docs/increments/21-parallel-refine.md,
# section 3 and "Pinned by the red suite (21a)"). Bit-identical: 14's T6 and
# 18's golden digests are the output oracle and stay as they are.
# test_refinement_chunks_dynamic (QW1, the dynamic block scheduler) is the
# invariant-critical suite; test_refinement_active (QW3) and
# test_mesh_lawson_stack (legalise_around's caller-owned stack) compare with
# today's code. The first two run in the TSan job; test_mesh_lawson_stack
# starts no threads, so it does not (h17 §4c).

add_terrain_test(test_refinement_chunks_dynamic unit/test_refinement_chunks_dynamic.cpp)
target_link_libraries(test_refinement_chunks_dynamic PRIVATE Threads::Threads)
add_terrain_backend_test(test_refinement_active unit/test_refinement_active.cpp)
target_link_libraries(test_refinement_active PRIVATE Threads::Threads)
add_terrain_backend_test(test_mesh_lawson_stack unit/test_mesh_lawson_stack.cpp)
target_link_libraries(test_mesh_lawson_stack PRIVATE Threads::Threads)

# The integer incircle, increment 21b (docs/increments/21-parallel-refine.md,
# section 3 QW2 and "Pinned by the red suite (21b)"). Bit-identical: 14's T6
# and 18's golden digests are the output oracle and stay as they are.
# test_mesh_lattice_incircle is the invariant-critical suite: agreement with
# DetriaExact on the frame doubles wherever lattice_incircle answers, refusal
# wherever QW2's conditions fail, and must_flip's decisions as today's. It
# starts no threads, so it is not in the TSan job.

add_terrain_backend_test(test_mesh_lattice_incircle unit/test_mesh_lattice_incircle.cpp)

# The verifier's pair search by sort and sweep, increment 16b-0
# (docs/increments/16b-terrain-polygons.md, R8, I6 and "Pinned by the red suite
# (16b-0)"). prop_noding_verifier_sweep is the invariant-critical suite:
# sweep_box_pairs reports exactly the brute force's box pairs, build() returns
# the brute-force verifier's status (tests/cpp/support/verifier_oracle.hpp,
# frozen at 6b4fcb9), and node<K>'s output digests are today's. It starts no
# threads, so it is not in the TSan job.

add_terrain_backend_test(prop_noding_verifier_sweep property/prop_noding_verifier_sweep.cpp)

# The fine catchment's flood, increment 22 PR 1
# (docs/increments/22-auto-catchment.md, "Flow and membership: one flood" and
# "The red suites"). Not named invariant-critical; the descent oracle and the
# invariants on random DEMs are the defence. It starts no threads.
add_terrain_test(test_hydrology_upstream unit/test_hydrology_upstream.cpp)

# Accumulation from the same flood, increment 29 PR 1
# (docs/increments/29-nve-reference-catchments.md, "Accumulation, from the same
# flood (C++)" and "The red suites"). INVARIANT-CRITICAL: the exact oracle
# against upstream, node for node; mutation target is accumulate's on_reach
# callback and its reverse pass. It starts no threads.
add_terrain_test(test_hydrology_accumulate unit/test_hydrology_accumulate.cpp)

# The outline reduction, increment 22 PR 2 (docs/increments/22-auto-catchment.md,
# "The reduction" and "The red suites"): area-preserving segment collapse with
# the exact kernel. A backend test because it names DefaultKernel. Not named
# invariant-critical; it starts no threads.
add_terrain_backend_test(test_area_collapse unit/test_area_collapse.cpp)

# The final check, increment 15c-1 (docs/increments/15c-geographic-dem.md, D4,
# D5 and "Tests for @tester"): the check-point store (CP1), refine_points on
# placed points (RP1, RP2, RP5 to RP8, the scan's tie-break, the unfrozen-store
# refusal, and LatticeMesh::split_inside taking an off-node MeshVertex), and its
# property suite (RP3's tolerance and Delaunay oracles, RP5's edge splits with
# work on both sides of an interior edge, RP4's thread and add-order sweep).
# The record names RP3 and the membership and tie-break tests (RP2, RP4, RP5)
# invariant-critical. The two refine_points suites start threads; both run in
# the TSan job.
add_terrain_test(test_refinement_check_points unit/test_refinement_check_points.cpp)
# Increment 15e fix 4 (docs/increments/15e-memory-fixes.md): the store's slab
# arena, its reserved_points() bound and rows across chunk boundaries. Its own
# target, so a compile failure in it leaves every other suite built and run.
add_terrain_test(test_refinement_check_points_arena unit/test_refinement_check_points_arena.cpp)
add_terrain_backend_test(test_refinement_refine_points unit/test_refinement_refine_points.cpp)
target_link_libraries(test_refinement_refine_points PRIVATE Threads::Threads)
add_terrain_backend_test(prop_refinement_refine_points property/prop_refinement_refine_points.cpp)
target_link_libraries(prop_refinement_refine_points PRIVATE Threads::Threads)

# The edge strip's generator and store, increment 15f-1
# (docs/increments/15f-edge-strip.md, D2, D3 and "Tests for @tester", CC1 to
# CC9), and detail::lattice_position, extracted from refine's to_lattice with no
# change in behaviour (D2 step 1). Two targets, one per header under test, so
# lattice_position's suite stays free of the strip. Not invariant-critical (the
# mutation suite is 15f-2's ES2, ES3 and ES5); neither starts threads, but
# refine.hpp includes the chunk scheduler, so both link Threads.
# test_refinement_lattice_position also pins to_lattice's constraint lookup
# (increment 15f-4, B4).
add_terrain_backend_test(test_refinement_constraint_points
                         unit/test_refinement_constraint_points.cpp)
target_link_libraries(test_refinement_constraint_points PRIVATE Threads::Threads)
add_terrain_backend_test(test_refinement_lattice_position
                         unit/test_refinement_lattice_position.cpp)
target_link_libraries(test_refinement_lattice_position PRIVATE Threads::Threads)

# The edge strip in the refinement loop, increment 15f-2
# (docs/increments/15f-edge-strip.md, D4, E1 to E8 and "Tests for @tester",
# ES1 to ES16). Two targets, so the oracles stay buildable on their own:
# test_refinement_strip_oracle shows, against refine() alone, that every
# oracle in support/strip_oracle.hpp can fail (ES1's defect half among them);
# prop_refinement_edge_strip is the suite for refine_strip and
# refine_points(..., strip). ES2, ES3 and ES5 are the invariant-critical
# cases. Both start threads (refine's scheduler; ES4 sweeps 1, 2 and 8).
add_terrain_backend_test(test_refinement_strip_oracle unit/test_refinement_strip_oracle.cpp)
target_link_libraries(test_refinement_strip_oracle PRIVATE Threads::Threads)
# COST 100: its ES9 case runs 97-171 s under the sanitizers; without a cost
# ctest --parallel would start it near the end (h14 section 2c).
add_terrain_backend_test(prop_refinement_edge_strip property/prop_refinement_edge_strip.cpp
                         COST 100)
target_link_libraries(prop_refinement_edge_strip PRIVATE Threads::Threads)

# L16 (docs/increments/15f-edge-strip.md): detail::coincidence_radius, by
# hand-written powers of two; a helper with no mesh, so it has its own target.
add_terrain_backend_test(test_refinement_coincidence_radius unit/test_refinement_coincidence_radius.cpp)
target_link_libraries(test_refinement_coincidence_radius PRIVATE Threads::Threads)

# Frozen edges and the seam pass, increment 23b
# (docs/increments/23-basin-scale.md, "The seam protocol", K1 to K4 and "Tests
# @tester can write red", FE1 to FE5, SP1 to SP4). Each its own target, so a
# compile failure in one leaves the others built and run:
#   test_mesh_frozen              LatticeMesh's frozen mask and the quality
#                                 pass's skip (FE4);
#   test_refinement_scan_frozen   the scan's skip of nodes on a frozen edge (FE2);
#   test_refinement_frozen_oracle the oracles of support/frozen_oracle.hpp,
#                                 shown able to fail against today's refine;
#   prop_refinement_frozen        refine, refine_points and refine_strip with
#                                 frozen edges (FE1 to FE5, K1, K2, §3D);
#   prop_refinement_seam          refine_seam, and a piece refined with its seam
#                                 frozen after it (SP1 to SP4, K3, K4).
# Invariant-critical (mutation runs at green): FE2 to FE5 (test_mesh_frozen's
# FE4 cases, test_refinement_scan_frozen, prop_refinement_frozen) and SP1, SP2
# (prop_refinement_seam). The refine suites start threads.
add_terrain_backend_test(test_mesh_frozen unit/test_mesh_frozen.cpp)
add_terrain_backend_test(test_refinement_scan_frozen unit/test_refinement_scan_frozen.cpp)
add_terrain_backend_test(test_refinement_frozen_oracle unit/test_refinement_frozen_oracle.cpp)
target_link_libraries(test_refinement_frozen_oracle PRIVATE Threads::Threads)
add_terrain_backend_test(prop_refinement_frozen property/prop_refinement_frozen.cpp)
target_link_libraries(prop_refinement_frozen PRIVATE Threads::Threads)
add_terrain_backend_test(prop_refinement_seam property/prop_refinement_seam.cpp)
target_link_libraries(prop_refinement_seam PRIVATE Threads::Threads)

# The release-hardening guard, increment 24
# (docs/increments/24-release-hardening.md, sections 5 and 12, T1 to T3). The
# suite is told what RASPUTIN_HARDENING asked for, 1 or 0, and refuses to
# compile without it; its static_assert compares that with what the standard
# library reports, and its forked case checks that an out-of-range read stops
# the process. POSIX only (fork), which both CI platforms are.
#
# Registered by hand rather than with add_terrain_test (same links, same
# warnings) for one property: under RASPUTIN_HARDENING=OFF the forked case
# SKIPs, and a Catch2 run whose only case skipped exits 4, which CTest counts as
# a failure unless SKIP_RETURN_CODE says otherwise. Catch2 v3.6.0's
# catch_discover_tests does not set it, so the unchecked CI leg would go red on
# a correct build.
add_executable(test_build_hardening unit/test_build_hardening.cpp)
target_link_libraries(test_build_hardening PRIVATE
    terrain_headers
    terrain_test_support
    Catch2::Catch2WithMain
)
target_compile_options(test_build_hardening PRIVATE -Wall -Wextra -Wpedantic -Werror)
target_compile_definitions(test_build_hardening PRIVATE
    RASPUTIN_EXPECT_HARDENING=$<BOOL:${RASPUTIN_HARDENING}>)
# A latent hole: Catch2 exits with the number of failed assertions, so a run
# with exactly four failures would also exit 4 and read as skipped. Not
# reachable here: each case runs as its own ctest test and none can fail more
# than one assertion (T1 and T3 have one check each, T2 one CHECK after REQUIREs
# that stop the case). A case that can fail four times needs this revisited.
catch_discover_tests(test_build_hardening PROPERTIES SKIP_RETURN_CODE 4)

# Increment 27 (docs/increments/27-node-sampling.md, S5): RasterGeometry::node_at
# against refine's own node test. Its own target, because it includes refine.hpp
# and so needs the backend target and Threads, which test_raster does not.
# refine.hpp names DefaultKernel, hence the backend target.
add_terrain_backend_test(test_raster_node_at unit/test_raster_node_at.cpp)
target_link_libraries(test_raster_node_at PRIVATE Threads::Threads)

# The foot rule on every insertion path, increment 20c PR 20c-1
# (docs/increments/20c-soft-quality.md, R1 to R6 and "Tests @tester writes red
# first", 20c-1). The design's invariant-critical suite test_constraint_foot,
# in four targets so a compile failure in one leaves the others built and run:
#   test_constraint_foot          CF1, the helper constraint_foot on
#                                 hand-built meshes; starts no threads;
#   test_constraint_foot_quality  CF2, the quality start's foot, through
#                                 improve and through refine, with 14b's T3
#                                 and T6 re-run with the rule on;
#   test_constraint_foot_refine   CF3, refinement's search through a
#                                 neighbour, and the features-start off switch;
#   test_constraint_foot_final    CF4, the final check's foot on both paths,
#                                 its bound and fallback, T3 and T6, and the
#                                 off switch.
# The last three start threads through refine's scheduler; their T6 cases are
# the TSan job's (the design: "T6 is the TSan job's").
add_terrain_backend_test(test_constraint_foot unit/test_constraint_foot.cpp)
add_terrain_backend_test(test_constraint_foot_quality unit/test_constraint_foot_quality.cpp)
target_link_libraries(test_constraint_foot_quality PRIVATE Threads::Threads)
add_terrain_backend_test(test_constraint_foot_refine property/prop_constraint_foot_refine.cpp)
target_link_libraries(test_constraint_foot_refine PRIVATE Threads::Threads)
add_terrain_backend_test(test_constraint_foot_final property/prop_constraint_foot_final.cpp)
target_link_libraries(test_constraint_foot_final PRIVATE Threads::Threads)

# The soft criterion and the line split in the quality start, increment 20c
# PR 20c-2 (docs/increments/20c-soft-quality.md, R7 and R8 and "Tests @tester
# writes red first", 20c-2). Two targets, so a compile failure in one leaves
# the other built and run:
#   test_quality_gain         the design's invariant-critical suite: T-P1 (the
#                             predicted cavity is what legalise_around writes),
#                             T-P2 (no accepted candidate lowers the worst
#                             angle) and LS1 (R8's split), on improve; starts
#                             no threads;
#   test_quality_gain_refine  T-P3 (gain -1 is 20c-1's output, digests
#                             recorded from 41bda81a), determinism over
#                             threads, and the QA rules' section D oracles
#                             with the gain on, through refine. It starts
#                             threads; its determinism case is the TSan job's.
add_terrain_backend_test(test_quality_gain unit/test_quality_gain.cpp)
add_terrain_backend_test(test_quality_gain_refine property/prop_quality_gain_refine.cpp)
target_link_libraries(test_quality_gain_refine PRIVATE Threads::Threads)

# A vertical tolerance that varies with distance to named lines, increment 33
# (docs/increments/33-feature-tolerance.md, sections 3, 4 and 9). Two targets:
#   test_refinement_line_tolerance   the field alone: the ramp (test 4), the
#                                    triangle-to-segment distance and its
#                                    search (test 6), make's refusals (test 10);
#                                    starts no threads;
#   prop_refinement_line_tolerance   the field in refine, refine_points and
#                                    refine_strip: tests 1, 2, 3, 5, 7, 8 and 9,
#                                    with the QA rules' section D oracles. It
#                                    starts threads; test 9 is the TSan job's.
# Tests 1 to 8 are the design's invariant-critical suite.
add_terrain_backend_test(test_refinement_line_tolerance unit/test_refinement_line_tolerance.cpp)
add_terrain_backend_test(prop_refinement_line_tolerance property/prop_refinement_line_tolerance.cpp)
target_link_libraries(prop_refinement_line_tolerance PRIVATE Threads::Threads)

# A vertical tolerance that follows the slope, increment 34
# (docs/increments/34-slope-tolerance.md, sections 3, 4 and 9). Three targets:
#   test_raster_steepness            test 2: every node's class against the
#                                    test's own Horn; starts threads (1 and 8);
#   test_refinement_slope_tolerance  the ramp (test 3), cell_class beside NoData
#                                    (test 6 (b), M9), make's refusals (test 9),
#                                    the histogram; starts no threads;
#   prop_refinement_slope_tolerance  the slope in refine, refine_strip and
#                                    refine_points: tests 1, 4, 5, 6, 7 and 8 and
#                                    the laziness count (M8), with the QA rules'
#                                    section D oracles. It starts threads; test 7
#                                    is the TSan job's.
# Tests 2, 4, 5, 6 and 8 are the design's invariant-critical suite. The files
# are guarded on terrain/raster/steepness.hpp and
# terrain/refinement/slope_tolerance.hpp: without them each compiles to one
# failing case that names the missing header (test 1 runs today), and the tree
# still builds.
add_terrain_test(test_raster_steepness unit/test_raster_steepness.cpp)
target_link_libraries(test_raster_steepness PRIVATE Threads::Threads)
add_terrain_backend_test(test_refinement_slope_tolerance unit/test_refinement_slope_tolerance.cpp)
add_terrain_backend_test(prop_refinement_slope_tolerance property/prop_refinement_slope_tolerance.cpp)
target_link_libraries(prop_refinement_slope_tolerance PRIVATE Threads::Threads)
