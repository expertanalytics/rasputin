# rasputin mesh — statistics

`rasputin mesh --dem glo30 --cache /Users/skavhaug/projects/rasputin_data/sweden_glo30_cache --out-crs EPSG:3006 --domain /Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --features /Users/skavhaug/projects/rasputin_data/sweden_corine/lagan_clc2018_3035.geojson --features-crs EPSG:3035 --features-map corine --tolerance 10 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-30d2/out/p1_r2_lagan_geojson.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/buffer-speed-2/docs/benchmarks/2026-10-07/30d-buffer/raw/runs/p1_r2_lagan_geojson_stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 4777 × 3638 (31 m) |
| domain vertices | 56260 (1 ring, 0 holes) |
| start vertices | 191590 |
| start triangles | 325950 |
| output vertices | 462541 |
| output triangles | 867612 |
| constraint edges | 201437 |
| p1_r2_lagan_geojson.vtk | 41.8 MB |

## Inputs

| item | value | name |
|---|---|---|
| The tiles or downloaded blocks used | Copernicus_DSM_COG_10_N56_00_E012_00_DEM; Copernicus_DSM_COG_10_N56_00_E013_00_DEM; Copernicus_DSM_COG_10_N56_00_E014_00_DEM; Copernicus_DSM_COG_10_N57_00_E012_00_DEM; Copernicus_DSM_COG_10_N57_00_E013_00_DEM; Copernicus_DSM_COG_10_N57_00_E014_00_DEM | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4326 | `dem_crs` |
| Conversion from the DEM's coordinate system | Inverse of SWEREF99 to WGS 84 (1) + SWEREF99 TM (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 31 m square grid in EPSG:3006, 3638 columns x 4777 rows | `resampled_grid` |
| The domain file and shape | lagan_mouth_svar2022_3006.geojson: 1 outline, 0 holes, 56260 vertices | `domain` |
| The domain's coordinate system | EPSG:3006 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| Each features file: layer, class map, features, lines, vertices | lagan_clc2018_3035.geojson, class map corine: 2790 features, 4527 lines, 274426 vertices | `features` |
| Each features file's coordinate system | EPSG:3035 | `features_crs` |
| Conversion from each features file's coordinate system | Inverse of Europe Equal Area 2001 + Inverse of SWEREF99 to ETRS89 (1) + SWEREF99 TM (with axis order normalized for visualization) | `features_transform` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.08 % | 4.53 % | 0.00238° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 17 | 235 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:3006 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.99975411631047 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Largest error against the resampled grid | 9.999847412109375 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 12340039 | `dem_nodes_checked` |
| Points that comparison added | 26708 | `dem_check_points_inserted` |
| Passes of that comparison | 12 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 1280115 | `line_points_checked` |
| Their largest difference from the mesh | 9.998990037059428 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 111 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 15 | `refinement_rounds` |
| Points added | 46144 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 107719 | `edge_flips` |
| Points added to improve the starting mesh | 192169 | `start_quality_points_inserted` |
| Tries skipped while improving it | 104936 | `start_quality_points_skipped` |
| Points moved onto lines while improving the starting mesh | 5930 | `start_quality_points_snapped_to_lines` |
| Points refinement moved onto lines | 720 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 793 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 19 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 191530 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.046 | 0.1 % |
| decode | 19.325 | 33.9 % |
| features read | 30.196 | 53.0 % |
| features clip | 3.140 | 5.5 % |
| start mesh: build | 0.006 | 0.0 % |
| start mesh: node | 0.629 | 1.1 % |
| start mesh: triangulate | 0.197 | 0.3 % |
| start mesh: constraint edges | 0.424 | 0.7 % |
| refine | 0.608 | 1.1 % |
| refine: legalise start | 0.009 | 0.0 % |
| refine: start quality | 0.382 | 0.7 % |
| refine: scan (parallel) | 0.053 | 0.1 % |
| refine: split + flip (serial) | 0.040 | 0.1 % |
| refine: setup + output | 0.125 | 0.2 % |
| edge strip: generate | 0.052 | 0.1 % |
| check points: project | 0.477 | 0.8 % |
| check points: store | 0.277 | 0.5 % |
| final check: scan (parallel) | 0.328 | 0.6 % |
| final check: split + flip (serial) | 0.035 | 0.1 % |
| trim | 0.036 | 0.1 % |
| land cover | 0.460 | 0.8 % |
| write: encode | 0.015 | 0.0 % |
| write: disk | 0.023 | 0.0 % |
| other | 0.651 | 1.1 % |
| **total** | **56.927** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.163 s, not included above.
