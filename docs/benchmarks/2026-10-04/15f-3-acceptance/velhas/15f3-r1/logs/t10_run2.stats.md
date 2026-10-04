# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983 --tolerance 10 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/velhas/15f3-r1/t10.bin.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/docs/benchmarks/2026-10-04/15f-3-acceptance/velhas/15f3-r1/logs/t10_run2.stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 7346 × 4204 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 358020 |
| output triangles | 708612 |
| constraint edges | 7426 |
| t10.bin.vtk | 25.7 MB |

## Inputs

| item | value | name |
|---|---|---|
| The tiles or downloaded blocks used | anadem_v1_compressed_COG | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4674 | `dem_crs` |
| Conversion from the DEM's coordinate system | UTM zone 23S (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 30 m square grid in EPSG:31983, 4204 columns x 7346 rows | `resampled_grid` |
| The domain file and shape | bho2017_5k_76949_outline_epsg4674.geojson: 1 outline, 0 holes, 7309 vertices | `domain` |
| The domain's coordinate system | EPSG:4674 | `domain_crs` |
| Conversion from the domain's coordinate system | UTM zone 23S (with axis order normalized for visualization) | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| DEM nodes very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.19 % | 0.0525° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 27 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:31983 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999996159862349 | `max_error_m` |
| DEM | anadem-v1 | `dem_source` |
| Credit | Ag\xeancia Nacional de \xc1guas e Saneamento B\xe1sico. (2025). ANADEM: A Digital Terrain Model for South America. Distributed by OpenTopography. https://doi.org/10.5069/G9736P4G. | `dem_credit` |
| Licence | CC BY 4.0 (Creative Commons Attribution 4.0 International) | `licence_note` |
| Please cite | Laipelt, L., et al. (2024). ANADEM: A Digital Terrain Model for South America. Remote Sensing, 16(13), 2321. https://doi.org/10.3390/rs16132321 | `cite` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Largest error against the resampled grid | 9.999976905616563 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 18791066 | `dem_nodes_checked` |
| Points that comparison added | 29406 | `dem_check_points_inserted` |
| Passes of that comparison | 11 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 73511 | `line_points_checked` |
| Their largest difference from the mesh | 9.957514884958755 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 45 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 35 | `refinement_rounds` |
| Points added | 310810 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 789337 | `edge_flips` |
| Points added to improve the starting mesh | 10495 | `start_quality_points_inserted` |
| Tries skipped while improving it | 341 | `start_quality_points_skipped` |
| Points moved onto lines | 72 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Starting-mesh vertices not on a DEM node | 7309 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 0.2 % |
| decode | 1.068 | 27.4 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.2 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.012 | 0.3 % |
| refine | 0.445 | 11.4 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.013 | 0.3 % |
| refine: scan (parallel) | 0.213 | 5.5 % |
| refine: split + flip (serial) | 0.185 | 4.8 % |
| refine: setup + output | 0.034 | 0.9 % |
| edge strip: generate | 0.003 | 0.1 % |
| check points: project | 0.557 | 14.3 % |
| check points: store | 0.421 | 10.8 % |
| final check: scan (parallel) | 0.409 | 10.5 % |
| final check: split + flip (serial) | 0.032 | 0.8 % |
| trim | 0.028 | 0.7 % |
| write: encode | 0.017 | 0.4 % |
| write: disk | 0.005 | 0.1 % |
| other | 0.880 | 22.6 % |
| **total** | **3.898** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.141 s, not included above.
