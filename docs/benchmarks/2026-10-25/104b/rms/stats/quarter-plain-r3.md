# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin/.claude/worktrees/tri-budget-rms/tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain /Users/skavhaug/projects/rasputin/.claude/worktrees/tri-budget-rms/docs/benchmarks/2026-09-26/quarter.geojson --tolerance 1 --out /tmp/perf-104b-meshes/quarter-plain-r3.vtk --binary --stats docs/benchmarks/2026-10-25/104b/rms/stats/quarter-plain-r3.md --record docs/benchmarks/2026-10-25/104b/rms/stats/quarter-plain-r3.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 3001 × 3001 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 536 |
| start triangles | 534 |
| output vertices | 214645 |
| output triangles | 428218 |
| constraint edges | 1070 |
| quarter-plain-r3.vtk | 15.5 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 3001 columns x 3001 rows, 10 m apart | `dem_grid` |
| Unit of z | metres | `dem_vertical_unit` |
| The domain file and shape | quarter.geojson: 1 outline, 0 holes, 536 vertices | `domain` |
| The domain's coordinate system | EPSG:25833 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 10 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 0.82 % | 0.396° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 18 | 184 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 1 | `tolerance_m` |
| Largest height error | 0.9999796549479072 | `max_error_m` |
| DEM | 7908_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 23411 | `line_points_checked` |
| Their largest difference from the mesh | 0.9995651245117188 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 1 | `line_points_inserted` |
| DEM nodes the line check added | 0 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 41 | `refinement_rounds` |
| Points added | 213464 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 445657 | `edge_flips` |
| Points added to improve the starting mesh | 644 | `start_quality_points_inserted` |
| Tries skipped while improving it | 0 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 0 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 6 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 235 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| domain read | 0.001 | 0.4 % |
| decode | 0.045 | 12.8 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.000 | 0.1 % |
| start mesh: triangulate | 0.000 | 0.1 % |
| start mesh: constraint edges | 0.001 | 0.2 % |
| refine | 0.165 | 46.9 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.002 | 0.5 % |
| refine: scan (parallel) | 0.045 | 12.8 % |
| refine: split + flip (serial) | 0.099 | 28.2 % |
| refine: setup + output | 0.019 | 5.4 % |
| edge strip: generate | 0.001 | 0.2 % |
| edge strip: scan (parallel) | 0.001 | 0.3 % |
| edge strip: split + flip (serial) | 0.004 | 1.1 % |
| trim | 0.016 | 4.4 % |
| write: encode | 0.009 | 2.6 % |
| write: disk | 0.003 | 0.8 % |
| other | 0.106 | 30.1 % |
| **total** | **0.353** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.081 s, not included above.
