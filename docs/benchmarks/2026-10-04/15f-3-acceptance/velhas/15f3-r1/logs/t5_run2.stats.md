# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983 --tolerance 5 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/velhas/15f3-r1/t5.bin.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/docs/benchmarks/2026-10-04/15f-3-acceptance/velhas/15f3-r1/logs/t5_run2.stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 7346 × 4204 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 934979 |
| output triangles | 1862154 |
| constraint edges | 7802 |
| t5.bin.vtk | 67.3 MB |

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
| minimum angle | 36.87° | 0.00 % | 0.28 % | 0.0558° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 22 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:31983 | `crs` |
| Tolerance | 5 | `tolerance_m` |
| Largest height error | 4.999999720676556 | `max_error_m` |
| DEM | anadem-v1 | `dem_source` |
| Credit | Ag\xeancia Nacional de \xc1guas e Saneamento B\xe1sico. (2025). ANADEM: A Digital Terrain Model for South America. Distributed by OpenTopography. https://doi.org/10.5069/G9736P4G. | `dem_credit` |
| Licence | CC BY 4.0 (Creative Commons Attribution 4.0 International) | `licence_note` |
| Please cite | Laipelt, L., et al. (2024). ANADEM: A Digital Terrain Model for South America. Remote Sensing, 16(13), 2321. https://doi.org/10.3390/rs16132321 | `cite` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Largest error against the resampled grid | 4.999982561383945 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 18791066 | `dem_nodes_checked` |
| Points that comparison added | 142874 | `dem_check_points_inserted` |
| Passes of that comparison | 12 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 73657 | `line_points_checked` |
| Their largest difference from the mesh | 4.998462104529381 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 275 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 37 | `refinement_rounds` |
| Points added | 774301 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 1885927 | `edge_flips` |
| Points added to improve the starting mesh | 10495 | `start_quality_points_inserted` |
| Tries skipped while improving it | 341 | `start_quality_points_skipped` |
| Points moved onto lines | 218 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Starting-mesh vertices not on a DEM node | 7309 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 0.1 % |
| decode | 0.988 | 15.0 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.1 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.012 | 0.2 % |
| refine | 0.986 | 15.0 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.013 | 0.2 % |
| refine: scan (parallel) | 0.327 | 5.0 % |
| refine: split + flip (serial) | 0.568 | 8.6 % |
| refine: setup + output | 0.079 | 1.2 % |
| edge strip: generate | 0.004 | 0.1 % |
| check points: project | 0.541 | 8.2 % |
| check points: store | 0.420 | 6.4 % |
| final check: scan (parallel) | 0.793 | 12.0 % |
| final check: split + flip (serial) | 0.179 | 2.7 % |
| trim | 0.071 | 1.1 % |
| write: encode | 0.028 | 0.4 % |
| write: disk | 0.048 | 0.7 % |
| other | 2.496 | 37.9 % |
| **total** | **6.586** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.377 s, not included above.
