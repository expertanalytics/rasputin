# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin/.claude/worktrees/tri-budget-rms/tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --threads 1 --rms --out /tmp/perf-104b-meshes/tile-rms-t1-r2.vtk --binary --stats docs/benchmarks/2026-10-25/104b/rms/stats/tile-rms-t1-r2.md --record docs/benchmarks/2026-10-25/104b/rms/stats/tile-rms-t1-r2.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 5051 × 5051 (10 m) |
| start vertices | 16384 |
| start triangles | 32258 |
| output vertices | 236525 |
| output triangles | 464290 |
| constraint edges | 872 |
| tile-rms-t1-r2.vtk | 16.9 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 5051 columns x 5051 rows, 10 m apart | `dem_grid` |
| Unit of z | metres | `dem_vertical_unit` |
| What refinement started from | every 40th DEM node | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 1 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 2.32 % | 0.63° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 74 | 324 | 197 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 1 | `tolerance_m` |
| Largest height error | 0.9999759197235107 | `max_error_m` |
| DEM | 7908_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 199 | `nodata_vertices_removed` |
| Root-mean-square height error | 0.06749796106317188 | `rms_error_m` |
| DEM nodes inside the mesh the height errors are over | 25233600 | `error_nodes_measured` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 23664 | `line_points_checked` |
| Their largest difference from the mesh | 0.9995651245117188 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 15664 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 0 | `line_points_inserted` |
| DEM nodes the line check added | 0 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 53 | `refinement_rounds` |
| Points added | 219837 | `points_inserted` |
| Of them, into triangles with a NoData corner | 7887 | `points_inserted_on_nodata` |
| Edge swaps | 445675 | `edge_flips` |
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
not in it. Threads used: 1.

| phase | seconds | share |
|---|---|---|
| decode | 0.052 | 5.7 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.003 | 0.3 % |
| start mesh: triangulate | 0.027 | 3.0 % |
| start mesh: constraint edges | 0.007 | 0.7 % |
| refine | 0.472 | 51.4 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.1 % |
| refine: scan (parallel) | 0.351 | 38.3 % |
| refine: split + flip (serial) | 0.093 | 10.1 % |
| refine: setup + output | 0.027 | 2.9 % |
| edge strip: generate | 0.001 | 0.1 % |
| edge strip: scan (parallel) | 0.006 | 0.7 % |
| edge strip: split + flip (serial) | 0.003 | 0.3 % |
| error summary | 0.220 | 23.9 % |
| trim | 0.017 | 1.8 % |
| write: encode | 0.009 | 1.0 % |
| write: disk | 0.003 | 0.3 % |
| other | 0.099 | 10.8 % |
| **total** | **0.918** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.090 s, not included above.
