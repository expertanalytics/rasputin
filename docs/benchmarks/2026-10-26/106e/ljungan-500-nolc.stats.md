# rasputin mesh — statistics

`rasputin mesh --dem glo30 --out-crs EPSG:3006 --simplify-domain --triangles 500 --stats docs/benchmarks/2026-10-26/106e/ljungan-500-nolc.stats.md --domain /Users/skavhaug/rasputin-output/opal-se2-2026-10-21/catchments-original/ljungan.geojson --out /Users/skavhaug/projects/rasputin_scratch/perf-106e/ljungan-500-nolc.vtk`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 3923 × 8237 (31 m) |
| domain vertices | 198 (1 ring, 0 holes) |
| start vertices | 198 |
| start triangles | 196 |
| output vertices | 345 |
| output triangles | 483 |
| constraint edges | 205 |
| ljungan-500-nolc.vtk | 25.1 kB |

## Inputs

| item | value | name |
|---|---|---|
| The tiles or downloaded blocks used | Copernicus_DSM_COG_10_N62_00_E012_00_DEM; Copernicus_DSM_COG_10_N62_00_E013_00_DEM; Copernicus_DSM_COG_10_N62_00_E014_00_DEM; Copernicus_DSM_COG_10_N62_00_E015_00_DEM; Copernicus_DSM_COG_10_N62_00_E016_00_DEM; Copernicus_DSM_COG_10_N62_00_E017_00_DEM; Copernicus_DSM_COG_10_N63_00_E012_00_DEM; Copernicus_DSM_COG_10_N63_00_E013_00_DEM; Copernicus_DSM_COG_10_N63_00_E014_00_DEM; Copernicus_DSM_COG_10_N63_00_E015_00_DEM; Copernicus_DSM_COG_10_N63_00_E016_00_DEM; Copernicus_DSM_COG_10_N63_00_E017_00_DEM | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4326 | `dem_crs` |
| Conversion from the DEM's coordinate system | Inverse of SWEREF99 to WGS 84 (1) + SWEREF99 TM (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 31 m square grid in EPSG:3006, 8237 columns x 3923 rows | `resampled_grid` |
| The domain file and shape | ljungan.geojson: 1 outline, 0 holes, 198 vertices | `domain` |
| The domain's coordinate system | EPSG:3006 | `domain_crs` |
| Band the domain outline was simplified within, m | 962.5741548764302 | `domain_band_m` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 10 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 34.00° | 0.62 % | 3.52 % | 0.214° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 4 | 8 | 10 | 0 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:3006 | `crs` |
| Tolerance | 246.8265826656846 | `tolerance_m` |
| Largest height error | 243.9420178559693 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Triangles asked for (about this many) | 500 | `triangles_asked` |
| Triangles at most | 525 | `triangles_cap` |
| Tolerances tried for the triangle budget, and the triangles each gave | 1e+06 m: 354; 107 m: 4139; 26.8 m: 73248; 493 m: 354; 135 m: 2262; 210 m: 647; 247 m: 483 | `budget_trials` |
| Root-mean-square height error | 69.19457180364434 | `rms_error_m` |
| DEM nodes inside the mesh the height errors are over | 14565966 | `error_nodes_measured` |
| Largest error against the resampled grid | 243.95175651129478 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 23093999 | `dem_nodes_checked` |
| Points that comparison added | 0 | `dem_check_points_inserted` |
| Passes of that comparison | 1 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 87513 | `line_points_checked` |
| Their largest difference from the mesh | 227.35302207689983 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 0 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 7 | `refinement_rounds` |
| Points added | 68 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 158 | `edge_flips` |
| Points added to improve the starting mesh | 79 | `start_quality_points_inserted` |
| Tries skipped while improving it | 43 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 27 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 7 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 198 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run per trial, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. The triangle budget ran 7 trials: the refine, edge strip, final check and trim rows sum all of them. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| domain read | 0.095 | 1.4 % |
| domain simplify | 2.534 | 37.0 % |
| decode | 1.632 | 23.8 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.000 | 0.0 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.000 | 0.0 % |
| check points: project | 0.885 | 12.9 % |
| check points: store | 0.465 | 6.8 % |
| refine | 0.457 | 6.7 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.002 | 0.0 % |
| refine: scan (parallel) | 0.431 | 6.3 % |
| refine: split + flip (serial) | 0.020 | 0.3 % |
| refine: setup + output | 0.005 | 0.1 % |
| edge strip: generate | 0.023 | 0.3 % |
| final check: scan (parallel) | 0.622 | 9.1 % |
| final check: split + flip (serial) | 0.003 | 0.0 % |
| trim | 0.004 | 0.1 % |
| error summary | 0.038 | 0.6 % |
| write: encode | 0.004 | 0.1 % |
| write: disk | 0.000 | 0.0 % |
| other | 0.085 | 1.2 % |
| **total** | **6.848** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.000 s, not included above.
