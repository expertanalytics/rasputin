# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 10 --binary --out /Users/skavhaug/projects/rasputin_scratch/sao_francisco_piece_meshes/26b-fractions/velhas_t10_lcoff.bin.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/mapbiomas-26b/docs/benchmarks/2026-10-12/26b-fractions/raw/velhas_t10_lcoff_r1.stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 328469 |
| output triangles | 649504 |
| constraint edges | 7432 |
| velhas_t10_lcoff.bin.vtk | 23.6 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 4200 columns x 7339 rows, 30 m apart | `dem_grid` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The domain file and shape | bho2017_5k_76949_outline_epsg4674.geojson: 1 outline, 0 holes, 7309 vertices | `domain` |
| The domain's coordinate system | EPSG:4674 | `domain_crs` |
| Conversion from the domain's coordinate system | UTM zone 23S (with axis order normalized for visualization) | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.15 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 13 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:31983 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999969482421875 | `max_error_m` |
| DEM | bho2017_5k_76949_anadem_epsg31983_30m.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 73518 | `line_points_checked` |
| Their largest difference from the mesh | 9.957514884958755 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 44 | `line_points_inserted` |
| DEM nodes the line check added | 7 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 38 | `refinement_rounds` |
| Points added | 313068 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 799011 | `edge_flips` |
| Points added to improve the starting mesh | 8036 | `start_quality_points_inserted` |
| Tries skipped while improving it | 1311 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 1031 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 5 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 74 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 7309 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 0.8 % |
| decode | 0.164 | 17.3 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.9 % |
| start mesh: triangulate | 0.004 | 0.4 % |
| start mesh: constraint edges | 0.009 | 0.9 % |
| refine | 0.532 | 56.1 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.020 | 2.1 % |
| refine: scan (parallel) | 0.281 | 29.7 % |
| refine: split + flip (serial) | 0.199 | 21.0 % |
| refine: setup + output | 0.032 | 3.3 % |
| edge strip: generate | 0.003 | 0.3 % |
| edge strip: scan (parallel) | 0.002 | 0.2 % |
| edge strip: split + flip (serial) | 0.004 | 0.5 % |
| trim | 0.025 | 2.7 % |
| write: encode | 0.008 | 0.9 % |
| write: disk | 0.009 | 0.9 % |
| other | 0.171 | 18.0 % |
| **total** | **0.949** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.144 s, not included above.
