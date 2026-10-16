# rasputin mesh — statistics

`rasputin mesh --dem glo30 --cache /Users/skavhaug/projects/rasputin_data/sweden_glo30_cache --out-crs EPSG:3006 --domain /Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/ljungan_flasjo_svar2022_3006.geojson --features /Users/skavhaug/projects/rasputin_data/sweden_corine/ljungan_flasjo_clc2018_3035.geojson --features-crs EPSG:3035 --features-map corine --tolerance 10 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-30d2/out/master_r1_ljungan_flasjo_geojson.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/buffer-speed-2/docs/benchmarks/2026-10-07/30d-buffer/raw/runs/master_r1_ljungan_flasjo_geojson_stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 1366 × 2099 (31 m) |
| domain vertices | 23259 (1 ring, 0 holes) |
| start vertices | 50068 |
| start triangles | 76521 |
| output vertices | 97930 |
| output triangles | 172153 |
| constraint edges | 51937 |
| master_r1_ljungan_flasjo_geojson.vtk | 8.8 MB |

## Inputs

| item | value | name |
|---|---|---|
| The tiles or downloaded blocks used | Copernicus_DSM_COG_10_N62_00_E012_00_DEM; Copernicus_DSM_COG_10_N62_00_E013_00_DEM; Copernicus_DSM_COG_10_N63_00_E012_00_DEM; Copernicus_DSM_COG_10_N63_00_E013_00_DEM | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4326 | `dem_crs` |
| Conversion from the DEM's coordinate system | Inverse of SWEREF99 to WGS 84 (1) + SWEREF99 TM (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 31 m square grid in EPSG:3006, 2099 columns x 1366 rows | `resampled_grid` |
| The domain file and shape | ljungan_flasjo_svar2022_3006.geojson: 1 outline, 0 holes, 23259 vertices | `domain` |
| The domain's coordinate system | EPSG:3006 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| Each features file: layer, class map, features, lines, vertices | ljungan_flasjo_clc2018_3035.geojson, class map corine: 449 features, 844 lines, 54268 vertices | `features` |
| Each features file's coordinate system | EPSG:3035 | `features_crs` |
| Conversion from each features file's coordinate system | Inverse of Europe Equal Area 2001 + Inverse of SWEREF99 to ETRS89 (1) + SWEREF99 TM (with axis order normalized for visualization) | `features_transform` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.05 % | 2.48 % | 0.024° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 20 | 41 | 1 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:3006 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999125803873994 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Largest error against the resampled grid | 9.9995005082111 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 3005077 | `dem_nodes_checked` |
| Points that comparison added | 1082 | `dem_check_points_inserted` |
| Passes of that comparison | 9 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 280177 | `line_points_checked` |
| Their largest difference from the mesh | 9.951234763313892 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 26 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 11 | `refinement_rounds` |
| Points added | 4979 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 11500 | `edge_flips` |
| Points added to improve the starting mesh | 40566 | `start_quality_points_inserted` |
| Tries skipped while improving it | 21611 | `start_quality_points_skipped` |
| Points moved onto lines while improving the starting mesh | 1235 | `start_quality_points_snapped_to_lines` |
| Points refinement moved onto lines | 91 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 36 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 1 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 50036 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.019 | 0.2 % |
| decode | 4.749 | 53.2 % |
| features read | 2.852 | 31.9 % |
| features clip | 0.557 | 6.2 % |
| start mesh: build | 0.001 | 0.0 % |
| start mesh: node | 0.103 | 1.2 % |
| start mesh: triangulate | 0.038 | 0.4 % |
| start mesh: constraint edges | 0.088 | 1.0 % |
| refine | 0.094 | 1.0 % |
| refine: legalise start | 0.002 | 0.0 % |
| refine: start quality | 0.054 | 0.6 % |
| refine: scan (parallel) | 0.009 | 0.1 % |
| refine: split + flip (serial) | 0.003 | 0.0 % |
| refine: setup + output | 0.025 | 0.3 % |
| edge strip: generate | 0.010 | 0.1 % |
| check points: project | 0.126 | 1.4 % |
| check points: store | 0.066 | 0.7 % |
| final check: scan (parallel) | 0.040 | 0.4 % |
| final check: split + flip (serial) | 0.002 | 0.0 % |
| trim | 0.007 | 0.1 % |
| land cover | 0.080 | 0.9 % |
| write: encode | 0.003 | 0.0 % |
| write: disk | 0.001 | 0.0 % |
| other | 0.100 | 1.1 % |
| **total** | **8.934** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.030 s, not included above.
