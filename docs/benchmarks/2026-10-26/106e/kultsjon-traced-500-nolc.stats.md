# rasputin mesh — statistics

`rasputin mesh --dem glo30 --out-crs EPSG:3006 --simplify-domain --triangles 500 --stats docs/benchmarks/2026-10-26/106e/kultsjon-traced-500-nolc.stats.md --domain /Users/skavhaug/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon-lake-seeded-check.geojson --out /Users/skavhaug/projects/rasputin_scratch/perf-106e/kultsjon-traced-500-nolc.vtk`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 1246 × 1976 (31 m) |
| domain vertices | 122 (1 ring, 0 holes) |
| start vertices | 122 |
| start triangles | 120 |
| output vertices | 320 |
| output triangles | 502 |
| constraint edges | 136 |
| kultsjon-traced-500-nolc.vtk | 23.6 kB |

## Inputs

| item | value | name |
|---|---|---|
| The tiles or downloaded blocks used | Copernicus_DSM_COG_10_N64_00_E014_00_DEM; Copernicus_DSM_COG_10_N64_00_E015_00_DEM; Copernicus_DSM_COG_10_N65_00_E014_00_DEM; Copernicus_DSM_COG_10_N65_00_E015_00_DEM | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4326 | `dem_crs` |
| Conversion from the DEM's coordinate system | Inverse of SWEREF99 to WGS 84 (1) + SWEREF99 TM (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 31 m square grid in EPSG:3006, 1976 columns x 1246 rows | `resampled_grid` |
| The domain file and shape | kultsjon-lake-seeded-check.geojson: 1 outline, 0 holes, 122 vertices | `domain` |
| The domain's coordinate system | EPSG:4326 | `domain_crs` |
| Band the domain outline was simplified within, m | 281.4083704119474 | `domain_band_m` |
| Conversion from the domain's coordinate system | Inverse of SWEREF99 to WGS 84 (1) + SWEREF99 TM (with axis order normalized for visualization) | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 10 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 34.42° | 0.40 % | 2.99 % | 0.329° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 5 | 8 | 9 | 0 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:3006 | `crs` |
| Tolerance | 152.3235640184785 | `tolerance_m` |
| Largest height error | 152.18614415633806 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Triangles asked for (about this many) | 500 | `triangles_asked` |
| Triangles at most | 525 | `triangles_cap` |
| Tolerances tried for the triangle budget, and the triangles each gave | 1e+06 m: 240; 97.2 m: 1133; 24.3 m: 8479; 210 m: 329; 126 m: 731; 152 m: 502 | `budget_trials` |
| Root-mean-square height error | 44.69166005367105 | `rms_error_m` |
| DEM nodes inside the mesh the height errors are over | 1354746 | `error_nodes_measured` |
| Largest error against the resampled grid | 151.91061405362893 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 2687045 | `dem_nodes_checked` |
| Points that comparison added | 8 | `dem_check_points_inserted` |
| Passes of that comparison | 4 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 18494 | `line_points_checked` |
| Their largest difference from the mesh | 151.83698310202828 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 2 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 8 | `refinement_rounds` |
| Points added | 130 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 255 | `edge_flips` |
| Points added to improve the starting mesh | 60 | `start_quality_points_inserted` |
| Tries skipped while improving it | 28 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 23 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 12 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 122 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run per trial, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. The triangle budget ran 6 trials: the refine, edge strip, final check and trim rows sum all of them. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.2 % |
| domain simplify | 0.170 | 17.5 % |
| decode | 0.429 | 44.2 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.000 | 0.0 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.000 | 0.0 % |
| check points: project | 0.131 | 13.5 % |
| check points: store | 0.055 | 5.7 % |
| refine | 0.064 | 6.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.1 % |
| refine: scan (parallel) | 0.060 | 6.1 % |
| refine: split + flip (serial) | 0.002 | 0.2 % |
| refine: setup + output | 0.001 | 0.1 % |
| edge strip: generate | 0.003 | 0.3 % |
| final check: scan (parallel) | 0.085 | 8.8 % |
| final check: split + flip (serial) | 0.001 | 0.1 % |
| trim | 0.001 | 0.1 % |
| error summary | 0.005 | 0.5 % |
| write: encode | 0.004 | 0.4 % |
| write: disk | 0.000 | 0.0 % |
| other | 0.020 | 2.0 % |
| **total** | **0.972** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.000 s, not included above.
