# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6901_3_10m_z33.tif --tolerance 10 --tolerance-slope 2 25 35 --out /tmp/perf-104b-meshes/romsdal-slope-plain-r5.vtk --binary --stats docs/benchmarks/2026-10-25/104b/rms/stats/romsdal-slope-plain-r5.md --record docs/benchmarks/2026-10-25/104b/rms/stats/romsdal-slope-plain-r5.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 5051 × 5051 (10 m) |
| start vertices | 16384 |
| start triangles | 32258 |
| output vertices | 1549323 |
| output triangles | 3096202 |
| constraint edges | 2442 |
| romsdal-slope-plain-r5.vtk | 111.5 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 5051 columns x 5051 rows, 10 m apart | `dem_grid` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| What refinement started from | every 40th DEM node | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 10 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 1.80 % | 0.935° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 21 | 1083 | 1 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999908093780277 | `max_error_m` |
| Tolerance on steep ground | 2 | `tolerance_slope_m` |
| Slope where the tolerance starts to tighten, and where it reaches the steep value | 25 to 35 | `tolerance_slope_deg` |
| DEM nodes held to less than the general tolerance because of their slope | 7917630 of 25512601 DEM nodes (31.0 %) | `slope_nodes_tightened` |
| Largest error as a share of the tolerance its slope allows | 1 | `max_error_slope_share_of_tolerance` |
| DEM | 6901_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 38132 | `line_points_checked` |
| Their largest difference from the mesh | 9.955078125 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 174 | `line_points_inserted` |
| DEM nodes the line check added | 28 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 30 | `refinement_rounds` |
| Points added | 1532234 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 3323985 | `edge_flips` |
| Points added to improve the starting mesh | 503 | `start_quality_points_inserted` |
| Tries skipped while improving it | 252 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 0 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 0 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 0 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| decode | 0.067 | 1.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.003 | 0.1 % |
| start mesh: triangulate | 0.028 | 0.8 % |
| start mesh: constraint edges | 0.007 | 0.2 % |
| slope | 0.106 | 2.9 % |
| refine | 2.021 | 56.2 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.0 % |
| refine: scan (parallel) | 0.569 | 15.8 % |
| refine: split + flip (serial) | 1.178 | 32.8 % |
| refine: setup + output | 0.273 | 7.6 % |
| edge strip: generate | 0.001 | 0.0 % |
| edge strip: scan (parallel) | 0.009 | 0.3 % |
| edge strip: split + flip (serial) | 0.024 | 0.7 % |
| trim | 0.119 | 3.3 % |
| write: encode | 0.039 | 1.1 % |
| write: disk | 0.015 | 0.4 % |
| other | 1.155 | 32.1 % |
| **total** | **3.594** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.684 s, not included above.
