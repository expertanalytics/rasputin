# Increment 137: kill record on green a1d6c1d9

Written by `@tester` on 2026-11-04 after code review round 1; copied verbatim into the branch by the main session.


Scratch copy of a1d6c1d9 (tools/scratch_copy.py), own CMake Release build
with hardening on: `_core`, `prop_statistics_height_shares`,
`test_land_cover_overlap`. Before each mutant: old binaries deleted, rebuilt,
`__pycache__` purged, PYTHONDONTWRITEBYTECODE=1, `_core` checked to load from
the copy. Suite per mutant: the five 137 Python files, `test_land_cover_overlap.py`,
and the two C++ binaries. Baseline: 214 Python passed, both C++ green.
Three mutants first failed to compile under -Werror (unused variable) and
were rewritten so they compile; those reruns are the ones counted.

Two real survivors, each killed by a new test in its own commit:
- aeeda8b9 `test_dem_height_shares.py::TestClassEdges::test_a_residual_past_n_plus_1_metres_stays_inside_the_classes` kills W1.
- 483ce7a9 `test_dem_heights.py::TestHalfCount::test_a_tolerance_on_an_edge_whose_quotient_rounds_low` kills H1.

## Table, one row per design target (hpp = height_shares.hpp, py = dem_heights.py, ov = overlap.hpp)

| # | target | mutant | where | result (Python kills + C++ kills) |
|---|---|---|---|---|
| 1 | lower class edge | `<` to `<=` | hpp:71 | survives, equivalent: the raise step (hpp:73) moves k back up |
| 2 | upper class edge | `>=` to `>` | hpp:73 | killed, 1: TestClassEdges::test_a_residual_on_an_edge_whose_quotient_rounds_low |
| 2b | (classes not centred, `floor(r)`) | floor(q), edges k w, (k+1) w | hpp:65,71,73 | killed, 19 |
| 3 | division by w | `rho / w` to `rho` | hpp:65 | survived the red suite; killed by aeeda8b9, 1 |
| 3b | division by w | `rho / w` to `rho * w` | hpp:65 | killed, 1: TestClassEdges::test_the_width_divides_the_residual |
| 4 | clamp to +-half | `std::clamp(k, ...)` to `k` | hpp:75 | killed, 21 |
| 4b | (early clamps removed) | drop hpp:66-69 | hpp:66 | survives, equivalent on arm64: guards only an int overflow of q |
| 5 | pixel-centre evaluation | residual at the piece centroid | hpp:148 | killed, 10 |
| 6 | held-node test | `< -1e-9` to `<= 0` (strict, no margin) | hpp:130 | killed, 8 |
| 6b | held-node test | `< -1e-9` to `< 0` (no margin) | hpp:130 | killed, 3 |
| 6c | (span over all pieces) | held test bypassed | hpp:153 | killed, 12 |
| 7 | fold to the span's ends | low and high ends swapped | hpp:167-170 | killed, 8 |
| 7b | fold | no fold | hpp:167-170 | killed, 11 |
| 8 | no-node class | class 0 | hpp:159 | killed, 5 |
| 8b | no-node class | pieces' spread | hpp:136,151,159 | killed, 4 |
| 9 | NoData branch | NoData counted in class 0 | hpp:144 | killed, 12 + 3 C++ |
| 9b | NoData branch | NoData ignored | hpp:143 | killed, 12 + 3 C++ |
| 10 | denominator | shares over the valid weight | hpp:174 | killed, 9 + 3 C++ |
| 10b | denominator | NoData share over the valid weight | hpp:163 | killed, 9 + 3 C++ |
| (extra) | row_area ignored | `w = area` | hpp:141 | killed, 2 + 1 C++ |
| 11 | half_count | raise while `<`, not `<=` | py:46 | survived the red suite; killed by 483ce7a9, 1 |
| 11b | half_count | lower while `>=` | py:44 | survives, equivalent: the raise step repairs it |
| 11c | half_count | floor without + 1/2 | py:43 | survives, equivalent: the raise step repairs it |
| 11d | half_count | `n >= MAX_HALF` refused | py:48 | killed, 1 |
| 11e | half_count | n + 1 | py:53 | killed, 31 |
| 12 | edges_m | j from -half | py:58 | killed, 30 + 12 errors |
| 12b | edges_m | j - 1/2 | py:58 | killed, 18 |
| 12c | edges_m | width dropped | py:58 | killed, 5 |
| 13 | whole-pixel test | `full_begin` ceil to floor | ov:215 | killed, 45 + 6 C++ |
| 13b | whole-pixel test | `full_end` floor to ceil | ov:216 | killed, 41 + 5 C++ |
| 13c | whole-pixel test | `c < full_end` to `<=` | ov:222 | killed, 52 + 12 C++ |
| 13d | whole-pixel test | no whole pixels (all clipped) | ov:222 | survives, equivalent: a clipped unit square has area exactly 1 and its centre is held |
| 14 | piece area | cut piece counted 1 | ov:236 | killed, 52 + 14 C++ |
| 14b | piece area | whole pixel 0.5 | ov:236 | killed, 47 + 6 C++ |

## Every killing test, per mutant (from the reruns; W1 and H1 shown with the new tests)

### A1_cut_counted_whole
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestClassEdges::test_6_the_plane_runs_through_the_written_z
- test_dem_height_shares.py::TestClassEdges::test_a_residual_on_an_edge_whose_quotient_rounds_low
- test_dem_height_shares.py::TestClassEdges::test_a_width_of_1_216_keeps_its_own_edges
- test_dem_height_shares.py::TestClassEdges::test_every_edge_and_both_clamps_where_the_design_says
- test_dem_height_shares.py::TestClassEdges::test_the_width_divides_the_residual
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_on_the_edge_is_held
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_within_the_margin_is_held
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestRefusals::test_a_float32_dem_is_taken
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_height_shares.py::TestRowAreas::test_two_triangles_over_rows_of_area_one_and_three
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[ne-sw]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[nw-se]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- test_dem_height_shares.py::TestStepsAtTheBinding::test_without_a_sentinel_the_sentinel_value_is_a_height
- test_dem_heights_off.py::test_class_areas_keeps_its_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[28]
- test_land_cover_overlap.py::TestPerCellSums::test_each_triangles_pieces_sum_to_its_area
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[26]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[27]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[28]
- test_land_cover_overlap.py::TestRefusals::test_the_window_edge_is_inside
- C++ prop_statistics_height_shares: a float DEM and its sentinel
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
- C++ test_land_cover_overlap: 1. a triangle wholly inside one cell: one entry, that class, its area times the row's
- C++ test_land_cover_overlap: 10. random triangles: every row well formed and summing to the shoelace area
- C++ test_land_cover_overlap: 14. the window is [0, cols] x [0, rows], closed: on the edge is inside, past it is OutsideWindow
- C++ test_land_cover_overlap: 2. a right triangle with legs 2 on a 2 x 2 grid: 1, 0.5, 0.5, the fourth cell absent
- C++ test_land_cover_overlap: 3. vertices on grid corners: no zero-area entry from the cells they touch
- C++ test_land_cover_overlap: 4. an edge lying on a grid line: the row below it gets nothing
- C++ test_land_cover_overlap: 5. a long sliver along a row: ten positive pieces that sum to the triangle's area
- C++ test_land_cover_overlap: 6. a triangle a nanometre across keeps its positive area
- C++ test_land_cover_overlap: 7. rows weigh differently: each piece is multiplied by its row's cell area
- C++ test_land_cover_overlap: 8. the vertex order does not change the areas: a rotation, and the mirror orientation
- C++ test_land_cover_overlap: 9. two triangles on a shared edge: in every cell their pieces sum to the cell's covered area
### A2_whole_half
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestClassEdges::test_6_the_plane_runs_through_the_written_z
- test_dem_height_shares.py::TestClassEdges::test_a_residual_on_an_edge_whose_quotient_rounds_low
- test_dem_height_shares.py::TestClassEdges::test_a_width_of_1_216_keeps_its_own_edges
- test_dem_height_shares.py::TestClassEdges::test_every_edge_and_both_clamps_where_the_design_says
- test_dem_height_shares.py::TestClassEdges::test_the_width_divides_the_residual
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestRefusals::test_a_float32_dem_is_taken
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_height_shares.py::TestRowAreas::test_two_triangles_over_rows_of_area_one_and_three
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[ne-sw]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[nw-se]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- test_dem_height_shares.py::TestStepsAtTheBinding::test_without_a_sentinel_the_sentinel_value_is_a_height
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[28]
- test_land_cover_overlap.py::TestPerCellSums::test_each_triangles_pieces_sum_to_its_area
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[26]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[27]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[28]
- C++ prop_statistics_height_shares: a float DEM and its sentinel
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
- C++ test_land_cover_overlap: 10. random triangles: every row well formed and summing to the shoelace area
- C++ test_land_cover_overlap: 2. a right triangle with legs 2 on a 2 x 2 grid: 1, 0.5, 0.5, the fourth cell absent
- C++ test_land_cover_overlap: 8. the vertex order does not change the areas: a rotation, and the mirror orientation
### C1_no_clamp
- test_cli_mesh_dem_heights.py::TestConvergence::test_below_half_a_class_every_triangle_has_one
- test_cli_mesh_dem_heights.py::TestConvergence::test_every_count_is_at_most_2n_plus_1
- test_cli_mesh_dem_heights.py::TestConvergence::test_the_mean_count_strictly_falls
- test_cli_mesh_dem_heights.py::TestConvergence::test_the_one_class_share_never_falls
- test_cli_mesh_dem_heights.py::TestConvergence::test_the_triangle_counts
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_piece_file_states_the_edges
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestCutRun::test_four_piece_files_and_1163_triangles
- test_cli_mesh_dem_heights.py::TestReprojected::test_969_triangles_and_n_5
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestThreads::test_equal_arrays
- test_cli_mesh_dem_heights.py::TestTriangles::test_the_edges_follow_the_chosen_tolerance
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestClassEdges::test_6_the_plane_runs_through_the_written_z
- test_dem_height_shares.py::TestClassEdges::test_a_pixel_outside_the_triangle_counts_nowhere
- test_dem_height_shares.py::TestClassEdges::test_every_edge_and_both_clamps_where_the_design_says
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_test_4s_scene
- test_dem_height_shares.py::TestRefusals::test_a_float32_dem_is_taken
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### C2_no_early_clamp
### D1_nodata_in_class0
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestAllNoData::test_count_0_first_0_nodata_1
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- test_dem_heights.py::TestAllNoData::test_count_0_first_0_nodata_1_and_counted[-9999.0]
- test_dem_heights.py::TestAllNoData::test_count_0_first_0_nodata_1_and_counted[None]
- C++ prop_statistics_height_shares: a float DEM and its sentinel
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
### D2_nodata_ignored
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestAllNoData::test_count_0_first_0_nodata_1
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- test_dem_heights.py::TestAllNoData::test_count_0_first_0_nodata_1_and_counted[-9999.0]
- test_dem_heights.py::TestAllNoData::test_count_0_first_0_nodata_1_and_counted[None]
- C++ prop_statistics_height_shares: a float DEM and its sentinel
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
### E1_lower_lt_le
### E2_upper_ge_gt
- test_dem_height_shares.py::TestClassEdges::test_a_residual_on_an_edge_whose_quotient_rounds_low
### E2b_not_centred_floor
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestClassEdges::test_6_the_plane_runs_through_the_written_z
- test_dem_height_shares.py::TestClassEdges::test_a_residual_on_an_edge_whose_quotient_rounds_low
- test_dem_height_shares.py::TestClassEdges::test_a_width_of_1_216_keeps_its_own_edges
- test_dem_height_shares.py::TestClassEdges::test_every_edge_and_both_clamps_where_the_design_says
- test_dem_height_shares.py::TestClassEdges::test_the_width_divides_the_residual
- test_dem_height_shares.py::TestNoHeldNode::test_5_the_plane_at_the_pixel_centre_not_the_piece_centroid
- test_dem_height_shares.py::TestRefusals::test_a_float32_dem_is_taken
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_height_shares.py::TestResolution::test_classes_minus_5_to_plus_5_occur
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### F1_wrong_end
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### F2_no_fold
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestTriangles::test_the_edges_follow_the_chosen_tolerance
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_past_the_margin_is_not_held_and_its_piece_folds
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### G1_range_short
- test_cli_mesh_dem_heights.py::TestFormatsAgree::test_ascii_files_carry_the_same_numbers
- test_cli_mesh_dem_heights.py::TestFormatsAgree::test_read_mesh_returns_the_height_table
- test_cli_mesh_dem_heights.py::TestFormatsAgree::test_the_gpkg_keeps_the_zeros_inside_a_span
- test_cli_mesh_dem_heights.py::TestFormatsAgree::test_vtk_line_cells_get_nodata_0_and_no_rows
- test_cli_mesh_dem_heights.py::TestRecord::test_2_no_triangle_is_without_data_and_the_width_is_the_default
- test_cli_mesh_dem_heights.py::TestRecord::test_the_json_record_has_the_entries
- test_cli_mesh_dem_heights.py::TestRecord::test_the_ply_comments_carry_the_entries
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_cli_mesh_dem_heights.py::TestThreads::test_equal_arrays
- test_cli_mesh_dem_heights.py::TestTriangles::test_the_edges_follow_the_chosen_tolerance
- test_dem_heights.py::TestAllNoData::test_count_0_first_0_nodata_1_and_counted[-9999.0]
- test_dem_heights.py::TestAllNoData::test_count_0_first_0_nodata_1_and_counted[None]
- test_dem_heights.py::TestEdges::test_2_5_m_classes_reach_11_25
- test_dem_heights.py::TestEdges::test_edges_of_n_1_at_1_m_are_exact
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[0-1.0]
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[10-1.0]
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[4-1.216]
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[4-2.5]
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_tolerance_sets_n
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestOutside::test_a_held_part_holding_the_mesh_reads_the_sources_indices
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### G2_minus_half
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_piece_file_states_the_edges
- test_cli_mesh_dem_heights.py::TestFormatsAgree::test_the_same_table_and_edges
- test_cli_mesh_dem_heights.py::TestReprojected::test_969_triangles_and_n_5
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_cli_mesh_dem_heights.py::TestTriangles::test_the_edges_follow_the_chosen_tolerance
- test_dem_heights.py::TestEdges::test_2_5_m_classes_reach_11_25
- test_dem_heights.py::TestEdges::test_edges_of_n_1_at_1_m_are_exact
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[0-1.0]
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[10-1.0]
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[4-1.216]
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[4-2.5]
- test_dem_heights.py::TestInTheDemsCrs::test_the_tolerance_sets_n
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
### G3_no_width
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_dem_heights.py::TestEdges::test_2_5_m_classes_reach_11_25
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[4-1.216]
- test_dem_heights.py::TestEdges::test_the_edges_are_one_product_each[4-2.5]
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
### H1_nonstrict_margin
- test_dem_heights.py::TestHalfCount::test_a_tolerance_on_an_edge_whose_quotient_rounds_low
### H2_lower_ge
### H3_max_ge
- test_dem_heights.py::TestHalfCount::test_n[127.4-1-127]
### H4_floor_no_half
### H5_plus_one
- test_cli_mesh_dem_heights.py::TestConvergence::test_below_half_a_class_every_triangle_has_one
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_piece_file_states_the_edges
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestFormatsAgree::test_the_gpkg_keeps_the_zeros_inside_a_span
- test_cli_mesh_dem_heights.py::TestFormatsAgree::test_the_same_table_and_edges
- test_cli_mesh_dem_heights.py::TestReprojected::test_969_triangles_and_n_5
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_cli_mesh_dem_heights.py::TestTriangles::test_the_edges_follow_the_chosen_tolerance
- test_dem_heights.py::TestHalfCount::test_n[0.4-1-0]
- test_dem_heights.py::TestHalfCount::test_n[0.5-1-1]
- test_dem_heights.py::TestHalfCount::test_n[1-1-1]
- test_dem_heights.py::TestHalfCount::test_n[10-1-10]
- test_dem_heights.py::TestHalfCount::test_n[10-2.5-4]
- test_dem_heights.py::TestHalfCount::test_n[127.4-1-127]
- test_dem_heights.py::TestHalfCount::test_n[2-1-2]
- test_dem_heights.py::TestHalfCount::test_n[2.5-1-3]
- test_dem_heights.py::TestHalfCount::test_n[20-1-20]
- test_dem_heights.py::TestHalfCount::test_n[5-1.216-4]
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_tolerance_sets_n
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### N1_strict_no_margin
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_on_the_edge_is_held
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_within_the_margin_is_held
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_height_shares.py::TestResolution::test_the_same_count_for_every_triangle
- test_dem_height_shares.py::TestRowAreas::test_two_triangles_over_rows_of_area_one_and_three
### N2_no_margin
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_within_the_margin_is_held
### N3_span_all_pieces
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_past_the_margin_is_not_held_and_its_piece_folds
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestNoHeldNode::test_the_mean_residual_not_the_pieces_spread
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### O1_full_begin_floor
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_on_the_edge_is_held
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_past_the_margin_is_not_held_and_its_piece_folds
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_within_the_margin_is_held
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[ne-sw]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[nw-se]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- test_dem_height_shares.py::TestStepsAtTheBinding::test_without_a_sentinel_the_sentinel_value_is_a_height
- test_dem_heights_off.py::test_class_areas_keeps_its_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[28]
- test_land_cover_overlap.py::TestPerCellSums::test_each_triangles_pieces_sum_to_its_area
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[26]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[27]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[28]
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
- C++ test_land_cover_overlap: 10. random triangles: every row well formed and summing to the shoelace area
- C++ test_land_cover_overlap: 3. vertices on grid corners: no zero-area entry from the cells they touch
- C++ test_land_cover_overlap: 7. rows weigh differently: each piece is multiplied by its row's cell area
- C++ test_land_cover_overlap: 8. the vertex order does not change the areas: a rotation, and the mirror orientation
- C++ test_land_cover_overlap: 9. two triangles on a shared edge: in every cell their pieces sum to the cell's covered area
### O2_full_end_ceil
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_on_the_edge_is_held
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_within_the_margin_is_held
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[ne-sw]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[nw-se]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- test_dem_heights_off.py::test_class_areas_keeps_its_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[28]
- test_land_cover_overlap.py::TestPerCellSums::test_each_triangles_pieces_sum_to_its_area
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[26]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[27]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[28]
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
- C++ test_land_cover_overlap: 10. random triangles: every row well formed and summing to the shoelace area
- C++ test_land_cover_overlap: 3. vertices on grid corners: no zero-area entry from the cells they touch
- C++ test_land_cover_overlap: 8. the vertex order does not change the areas: a rotation, and the mirror orientation
- C++ test_land_cover_overlap: 9. two triangles on a shared edge: in every cell their pieces sum to the cell's covered area
### O3_full_end_incl
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_1_both_triangles_whichever_the_diagonal[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_18_the_count_follows_the_tolerance[steps.vtk]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_cli_mesh_dem_heights.py::TestTheSourcesWidth::test_20
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestClassEdges::test_6_the_plane_runs_through_the_written_z
- test_dem_height_shares.py::TestClassEdges::test_a_residual_on_an_edge_whose_quotient_rounds_low
- test_dem_height_shares.py::TestClassEdges::test_a_width_of_1_216_keeps_its_own_edges
- test_dem_height_shares.py::TestClassEdges::test_every_edge_and_both_clamps_where_the_design_says
- test_dem_height_shares.py::TestClassEdges::test_the_width_divides_the_residual
- test_dem_height_shares.py::TestFold::test_19_a_low_piece_centred_outside_folds_into_the_low_end
- test_dem_height_shares.py::TestFold::test_19_the_piece_centred_outside_folds_into_class_2
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_on_the_edge_is_held
- test_dem_height_shares.py::TestHeldNodeMargin::test_a_node_within_the_margin_is_held
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestNoHeldNode::test_the_mean_residual_not_the_pieces_spread
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestRefusals::test_a_float32_dem_is_taken
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_height_shares.py::TestRowAreas::test_two_triangles_over_rows_of_area_one_and_three
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[ne-sw]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_1_both_triangles_on_either_diagonal[nw-se]
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- test_dem_height_shares.py::TestStepsAtTheBinding::test_without_a_sentinel_the_sentinel_value_is_a_height
- test_dem_heights_off.py::test_class_areas_keeps_its_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestInTheDemsCrs::test_the_two_triangles_of_test_1
- test_dem_heights.py::TestInTheDemsCrs::test_the_width_is_the_callers
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_class_totals_over_the_mesh_are_the_classes_areas_inside_the_square[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_every_triangles_class_areas_are_geoss_in_cell_units[28]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[26]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[27]
- test_land_cover_overlap.py::TestF2AgainstGeos::test_with_a_cell_area_per_row_each_piece_is_weighted_by_its_rows[28]
- test_land_cover_overlap.py::TestPerCellSums::test_each_triangles_pieces_sum_to_its_area
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[26]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[27]
- test_land_cover_overlap.py::TestPerCellSums::test_the_pieces_of_a_cell_over_all_triangles_sum_to_its_covered_area[28]
- C++ prop_statistics_height_shares: a float DEM and its sentinel
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
- C++ test_land_cover_overlap: 10. random triangles: every row well formed and summing to the shoelace area
- C++ test_land_cover_overlap: 14. the window is [0, cols] x [0, rows], closed: on the edge is inside, past it is OutsideWindow
- C++ test_land_cover_overlap: 2. a right triangle with legs 2 on a 2 x 2 grid: 1, 0.5, 0.5, the fourth cell absent
- C++ test_land_cover_overlap: 3. vertices on grid corners: no zero-area entry from the cells they touch
- C++ test_land_cover_overlap: 5. a long sliver along a row: ten positive pieces that sum to the triangle's area
- C++ test_land_cover_overlap: 6. a triangle a nanometre across keeps its positive area
- C++ test_land_cover_overlap: 7. rows weigh differently: each piece is multiplied by its row's cell area
- C++ test_land_cover_overlap: 8. the vertex order does not change the areas: a rotation, and the mirror orientation
- C++ test_land_cover_overlap: 9. two triangles on a shared edge: in every cell their pieces sum to the cell's covered area
### O4_no_full
### P1_piece_centroid
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_spread_and_class_0_share_are_the_sources
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestNoHeldNode::test_5_the_plane_at_the_pixel_centre_not_the_piece_centroid
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[10.0]
- test_dem_height_shares.py::TestResolution::test_against_geos_and_every_row_sums_to_one[5.0]
- test_dem_heights.py::TestInTheDemsCrs::test_chunks_of_one_triangle_give_the_same_bytes
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
### Q1_share_over_valid
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- C++ prop_statistics_height_shares: a float DEM and its sentinel
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
### Q2_nodata_over_valid
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.ply]
- test_cli_mesh_dem_heights.py::TestSteps::test_2_a_nodata_node_moves_its_pixel_to_the_nodata_share[steps.vtk]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[0.5-0]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.0-10]
- test_dem_height_shares.py::TestAgainstTheOracle::test_the_same_table[1.216-3]
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestOneClassAtHalfZero::test_with_nodata
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nan_pixel_is_nodata_without_a_sentinel
- test_dem_height_shares.py::TestStepsAtTheBinding::test_2_a_nodata_node_weighs_on_the_nodata_share_in_the_denominator
- C++ prop_statistics_height_shares: a float DEM and its sentinel
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
- C++ prop_statistics_height_shares: each pixel's area once over a triangulation, 300 cases
### R1_row_area_ignored
- test_dem_height_shares.py::TestRowAreas::test_two_triangles_over_rows_of_area_one_and_three
- test_dem_heights.py::TestReprojected::test_against_the_oracle_on_the_sources_pixels
- C++ prop_statistics_height_shares: a whole pixel inside a triangle weighs exactly its row area
### W1_no_div
- test_dem_height_shares.py::TestClassEdges::test_a_residual_past_n_plus_1_metres_stays_inside_the_classes
### W2_mul
- test_dem_height_shares.py::TestClassEdges::test_the_width_divides_the_residual
### X1_class0
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_dem_height_shares.py::TestNoHeldNode::test_5_the_plane_at_the_pixel_centre_not_the_piece_centroid
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestNoHeldNode::test_the_mean_is_clamped_to_the_end_class
- test_dem_height_shares.py::TestNoHeldNode::test_the_mean_residual_not_the_pieces_spread
### X2_spread
- test_cli_mesh_dem_heights.py::TestCutRun::test_every_triangle_of_every_piece_against_the_oracle
- test_cli_mesh_dem_heights.py::TestReprojected::test_the_oracle_on_the_sources_pixels
- test_dem_height_shares.py::TestNoHeldNode::test_held_nodes_all_nodata_with_valid_pieces_take_the_mean
- test_dem_height_shares.py::TestNoHeldNode::test_the_mean_residual_not_the_pieces_spread

## The design's extra mutants (`@tester`, after code review round 2; scratch copy of a1d6c1d9 with HEAD's tests)

| # | target | mutant | where | result |
|---|---|---|---|---|
| 15 | resampled grid on the reprojected path | the tile's windows, `mesh_crs` None | cli.py:3710, 3738, 3912 | killed, 2: TestReprojected::test_spread_and_class_0_share_are_the_sources, ::test_the_oracle_on_the_sources_pixels |
| 16 | shared accumulator across threads | `weight` moved outside the per-block lambda | hpp:116 | killed, 6, in 3 of 3 runs: test_dem_height_shares TestThreads::test_threads_1_2_7_and_all_give_equal_bytes; cli TestThreads::test_equal_arrays, TestCutRun::test_every_triangle_of_every_piece_against_the_oracle, TestTriangles::test_the_edges_follow_the_chosen_tolerance, both TestReprojected tests. The C++ prop_statistics_height_shares test does not kill it (3 runs). |
| 17 | flag ignored in `_write_piece` | `dem_heights=None`, `face_heights=None`, no edges comment | cli.py:2779, 2786, 2789 | killed, 3: cli TestThreads::test_equal_arrays, TestCutRun::test_every_piece_file_states_the_edges, ::test_every_triangle_of_every_piece_against_the_oracle |
