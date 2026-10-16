# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/vtol-field/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin_data/banenor_banenettverk/corridor_domain.geojson --tolerance 20 --tolerance-near /Users/skavhaug/projects/rasputin_data/banenor_banenettverk/bergen_line3.geojson 1 --tolerance-ramp 0 3000 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/corridor_ramp.vtk --stats docs/benchmarks/2026-10-09/vtol-field/corridor/corridor.stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 13211 × 28376 (10 m) |
| domain vertices | 1154 (1 ring, 3 holes) |
| start vertices | 1154 |
| start triangles | 1158 |
| output vertices | 794927 |
| output triangles | 1587426 |
| constraint edges | 2432 |
| corridor_ramp.vtk | 57.2 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 28376 columns x 13211 rows, 10 m apart | `dem_grid` |
| The tiles or downloaded blocks used | 6600_1_10m_z33.tif; 6600_2_10m_z33.tif; 6600_3_10m_z33.tif; 6600_4_10m_z33.tif; 6601_1_10m_z33.tif; 6601_2_10m_z33.tif; 6601_3_10m_z33.tif; 6601_4_10m_z33.tif; 6602_3_10m_z33.tif; 6602_4_10m_z33.tif; 66m1_1_10m_z33.tif; 66m1_2_10m_z33.tif; 6700_1_10m_z33.tif; 6700_2_10m_z33.tif; 6700_3_10m_z33.tif; 6700_4_10m_z33.tif; 6701_1_10m_z33.tif; 6701_2_10m_z33.tif; 6701_3_10m_z33.tif; 6701_4_10m_z33.tif; 6702_3_10m_z33.tif; 6702_4_10m_z33.tif; 67m1_1_10m_z33.tif; 67m1_2_10m_z33.tif | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The domain file and shape | corridor_domain.geojson: 1 outline, 3 holes, 1154 vertices | `domain` |
| The domain's coordinate system | EPSG:25833 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Tolerance on the tolerance lines | 1 | `tolerance_near_m` |
| Distance from the lines over which the tolerance rises to its far value | 0 to 3000 | `tolerance_ramp_m` |
| The tolerance lines file, and its segments after simplifying | bergen_line3.geojson: 11877 segments | `tolerance_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.71 % | 0.00524° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 17 | 228 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 20 | `tolerance_m` |
| Largest height error | 19.998650429774443 | `max_error_m` |
| DEM | 6600_1_10m_z33.tif; 6600_2_10m_z33.tif; 6600_3_10m_z33.tif; 6600_4_10m_z33.tif; 6601_1_10m_z33.tif; 6601_2_10m_z33.tif; 6601_3_10m_z33.tif; 6601_4_10m_z33.tif; 6602_3_10m_z33.tif; 6602_4_10m_z33.tif; 66m1_1_10m_z33.tif; 66m1_2_10m_z33.tif; 6700_1_10m_z33.tif; 6700_2_10m_z33.tif; 6700_3_10m_z33.tif; 6700_4_10m_z33.tif; 6701_1_10m_z33.tif; 6701_2_10m_z33.tif; 6701_3_10m_z33.tif; 6701_4_10m_z33.tif; 6702_3_10m_z33.tif; 6702_4_10m_z33.tif; 67m1_1_10m_z33.tif; 67m1_2_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Largest height error where the tolerance is the lines' own | 0.9999643961588447 | `max_error_near_lines_m` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 211710 | `line_points_checked` |
| Their largest difference from the mesh | 19.99876889607185 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 189 | `line_points_inserted` |
| DEM nodes the line check added | 11 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 44 | `refinement_rounds` |
| Points added | 792285 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 1891802 | `edge_flips` |
| Points added to improve the starting mesh | 1287 | `start_quality_points_inserted` |
| Tries skipped while improving it | 215 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 201 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 1088 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 1 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 1154 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.0 % |
| decode | 1.361 | 21.4 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.001 | 0.0 % |
| start mesh: triangulate | 0.001 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.0 % |
| refine | 4.367 | 68.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.003 | 0.1 % |
| refine: scan (parallel) | 3.723 | 58.5 % |
| refine: split + flip (serial) | 0.550 | 8.6 % |
| refine: setup + output | 0.092 | 1.4 % |
| edge strip: generate | 0.008 | 0.1 % |
| edge strip: scan (parallel) | 0.005 | 0.1 % |
| edge strip: split + flip (serial) | 0.008 | 0.1 % |
| trim | 0.056 | 0.9 % |
| write: encode | 0.024 | 0.4 % |
| write: disk | 0.008 | 0.1 % |
| other | 0.525 | 8.2 % |
| **total** | **6.369** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.296 s, not included above.
