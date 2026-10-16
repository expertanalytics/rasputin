# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/20-nodata/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin/.claude/worktrees/20-nodata/tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/t20.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/20-nodata/docs/benchmarks/2026-10-04/20-nodata-acceptance/tile_fix.stats.md`

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
| t20.vtk | 16.9 MB |

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
| decode | 0.049 | 15.8 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.003 | 0.8 % |
| start mesh: triangulate | 0.029 | 9.3 % |
| start mesh: constraint edges | 0.008 | 2.6 % |
| refine | 0.192 | 61.1 % |
| refine: legalise start | 0.000 | 0.1 % |
| refine: start quality | 0.001 | 0.2 % |
| refine: scan (parallel) | 0.056 | 17.9 % |
| refine: split + flip (serial) | 0.107 | 34.0 % |
| refine: setup + output | 0.028 | 8.9 % |
| trim | 0.018 | 5.8 % |
| write: encode | 0.010 | 3.1 % |
| write: disk | 0.004 | 1.1 % |
| other | 0.001 | 0.3 % |
| **total** | **0.314** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.082 s, not included above.
