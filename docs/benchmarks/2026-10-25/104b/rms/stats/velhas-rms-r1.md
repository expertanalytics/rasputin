# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin/.claude/worktrees/tri-budget-rms/tests/fixtures/velhas/anadem_velhas.tif --out-crs EPSG:31983 --bbox 592000 7845500 599000 7854000 --tolerance 0.5 --rms --out /tmp/perf-104b-meshes/velhas-rms-r1.vtk --binary --stats docs/benchmarks/2026-10-25/104b/rms/stats/velhas-rms-r1.md --record docs/benchmarks/2026-10-25/104b/rms/stats/velhas-rms-r1.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 288 × 239 (30 m) |
| start vertices | 7857 |
| start triangles | 15360 |
| output vertices | 49955 |
| output triangles | 99236 |
| constraint edges | 672 |
| velhas-rms-r1.vtk | 3.6 MB |

## Inputs

| item | value | name |
|---|---|---|
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4674 | `dem_crs` |
| Conversion from the DEM's coordinate system | UTM zone 23S (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 30 m square grid in EPSG:31983, 239 columns x 288 rows | `resampled_grid` |
| What refinement started from | every 3rd DEM node | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 10 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.94° | 0.00 % | 1.73 % | 1.09° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 12 | 1 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:31983 | `crs` |
| Tolerance | 0.5 | `tolerance_m` |
| Largest height error | 0.4999983232074783 | `max_error_m` |
| DEM | anadem_velhas.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Root-mean-square height error | 0.2471039538012421 | `rms_error_m` |
| DEM nodes inside the mesh the height errors are over | 54411 | `error_nodes_measured` |
| Largest error against the resampled grid | 0.49993896484375 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 72859 | `dem_nodes_checked` |
| Points that comparison added | 18478 | `dem_check_points_inserted` |
| Passes of that comparison | 7 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 1458 | `line_points_checked` |
| Their largest difference from the mesh | 0.48442084174530464 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 0 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 9 | `refinement_rounds` |
| Points added | 23525 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 25449 | `edge_flips` |
| Points added to improve the starting mesh | 95 | `start_quality_points_inserted` |
| Tries skipped while improving it | 0 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 0 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 0 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 30 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 20 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 0 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| decode | 0.036 | 25.9 % |
| start mesh: build | 0.000 | 0.1 % |
| start mesh: node | 0.001 | 0.9 % |
| start mesh: triangulate | 0.010 | 7.4 % |
| start mesh: constraint edges | 0.003 | 2.3 % |
| refine | 0.013 | 9.1 % |
| refine: legalise start | 0.000 | 0.2 % |
| refine: start quality | 0.000 | 0.2 % |
| refine: scan (parallel) | 0.003 | 2.0 % |
| refine: split + flip (serial) | 0.006 | 4.6 % |
| refine: setup + output | 0.003 | 2.2 % |
| edge strip: generate | 0.000 | 0.1 % |
| check points: project | 0.013 | 9.5 % |
| check points: store | 0.002 | 1.1 % |
| final check: scan (parallel) | 0.008 | 5.7 % |
| final check: split + flip (serial) | 0.010 | 7.5 % |
| error summary | 0.009 | 6.6 % |
| trim | 0.004 | 3.0 % |
| write: encode | 0.005 | 3.7 % |
| write: disk | 0.001 | 0.5 % |
| other | 0.023 | 16.4 % |
| **total** | **0.138** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.020 s, not included above.
