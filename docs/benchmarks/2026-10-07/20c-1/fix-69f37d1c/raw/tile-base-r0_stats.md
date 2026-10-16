# rasputin mesh — statistics

`drive.py --pkg /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c1b/base/build-bench/pkg --json /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/tile-base-r0.json -- mesh --dem /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c1b/meshes/tile-base.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/benchmarks/2026-10-07/20c-1/fix-69f37d1c/raw/tile-base-r0_stats.md`

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
| tile-base.vtk | 16.9 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 5051 columns x 5051 rows, 10 m apart | `dem_grid` |
| Unit of z | metres | `dem_vertical_unit` |
| What refinement started from | every 40th DEM node | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| DEM nodes very close to a line were moved onto it | on | `snap_to_lines` |

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
| Refinement passes | 53 | `refinement_rounds` |
| Points added | 219837 | `points_inserted` |
| Of them, into triangles with a NoData corner | 7887 | `points_inserted_on_nodata` |
| Edge swaps | 445675 | `edge_flips` |
| Points added to improve the starting mesh | 503 | `start_quality_points_inserted` |
| Tries skipped while improving it | 252 | `start_quality_points_skipped` |
| Points moved onto lines | 0 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Starting-mesh vertices not on a DEM node | 0 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| decode | 0.047 | 11.6 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.003 | 0.7 % |
| start mesh: triangulate | 0.029 | 7.2 % |
| start mesh: constraint edges | 0.007 | 1.6 % |
| refine | 0.190 | 46.5 % |
| refine: legalise start | 0.000 | 0.1 % |
| refine: start quality | 0.001 | 0.2 % |
| refine: scan (parallel) | 0.058 | 14.1 % |
| refine: split + flip (serial) | 0.105 | 25.7 % |
| refine: setup + output | 0.026 | 6.4 % |
| edge strip: generate | 0.001 | 0.2 % |
| edge strip: scan (parallel) | 0.001 | 0.3 % |
| edge strip: split + flip (serial) | 0.001 | 0.2 % |
| trim | 0.017 | 4.2 % |
| write: encode | 0.009 | 2.2 % |
| write: disk | 0.003 | 0.7 % |
| other | 0.101 | 24.7 % |
| **total** | **0.408** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.075 s, not included above.
