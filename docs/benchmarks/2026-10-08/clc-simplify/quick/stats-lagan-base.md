# rasputin mesh — statistics

`bench.py _child --pkg /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/base/build-bench/pkg --threads 0 -- mesh --dem glo30 --cache /Users/skavhaug/projects/rasputin_data/sweden_glo30_cache --out-crs EPSG:3006 --domain /Users/skavhaug/projects/rasputin_data/sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson --features /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-crs EPSG:3035 --features-map corine --tolerance 10 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/lagan-base.vtk --stats /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/lagan-base.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 4777 × 3638 (31 m) |
| domain vertices | 56260 (1 ring, 0 holes) |
| start vertices | 191725 |
| start triangles | 326646 |
| output vertices | 427798 |
| output triangles | 798554 |
| constraint edges | 211916 |
| lagan-base.vtk | 39.1 MB |

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
| Each features file: layer, class map, features, lines, vertices | U2018_CLC2018_V2020_20u1.gpkg layer U2018_CLC2018_V2020_20u1, class map corine: 23 features, 4297 lines, 275145 vertices | `features` |
| Each features file's coordinate system | EPSG:3035 | `features_crs` |
| Conversion from each features file's coordinate system | Inverse of Europe Equal Area 2001 + Inverse of SWEREF99 to ETRS89 (1) + SWEREF99 TM (with axis order normalized for visualization) | `features_transform` |
| Land-cover borders closer than this made one, narrower gaps filled | 1 | `features_repair_m` |
| Borders between land-cover polygons of one class dropped | on | `features_merge_same_class` |
| Land-cover borders simplified by this much (0 = off) | 0 | `features_tolerance_m` |
| Land-cover borders this close to the outline moved onto it | 5 | `features_outline_snap_m` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 35.63° | 0.00 % | 4.31 % | 0.729° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 85 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:3006 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999930431586677 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Land-cover vertices before and after clean-up | 380826 after the clip, 385640 after clean-up | `land_cover_vertices` |
| Land-cover area inside the outline that the outline rule gave to another polygon, m2 | 16050.3 | `land_cover_area_moved_m2` |
| Largest error against the resampled grid | 9.999993296063082 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 12340039 | `dem_nodes_checked` |
| Points that comparison added | 28944 | `dem_check_points_inserted` |
| Passes of that comparison | 13 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 1284790 | `line_points_checked` |
| Their largest difference from the mesh | 9.998990037059428 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 66 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 17 | `refinement_rounds` |
| Points added | 54195 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 131636 | `edge_flips` |
| Points added to improve the starting mesh | 136181 | `start_quality_points_inserted` |
| Tries skipped while improving it | 110975 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 34540 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 3431 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 13322 | `start_quality_lines_split` |
| Points refinement moved onto lines | 476 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 690 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 12 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 191665 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.042 | 0.2 % |
| decode | 4.567 | 23.3 % |
| features read | 0.960 | 4.9 % |
| features clip | 9.731 | 49.7 % |
| features clip: clean-up | 5.654 | 28.9 % |
| start mesh: build | 0.006 | 0.0 % |
| start mesh: node | 0.662 | 3.4 % |
| start mesh: triangulate | 0.192 | 1.0 % |
| start mesh: constraint edges | 0.416 | 2.1 % |
| refine | 0.658 | 3.4 % |
| refine: legalise start | 0.009 | 0.0 % |
| refine: start quality | 0.431 | 2.2 % |
| refine: scan (parallel) | 0.052 | 0.3 % |
| refine: split + flip (serial) | 0.044 | 0.2 % |
| refine: setup + output | 0.123 | 0.6 % |
| edge strip: generate | 0.052 | 0.3 % |
| check points: project | 0.685 | 3.5 % |
| check points: store | 0.275 | 1.4 % |
| final check: scan (parallel) | 0.303 | 1.5 % |
| final check: split + flip (serial) | 0.034 | 0.2 % |
| trim | 0.032 | 0.2 % |
| land cover | 0.394 | 2.0 % |
| write: encode | 0.015 | 0.1 % |
| write: disk | 0.006 | 0.0 % |
| other | 0.556 | 2.8 % |
| **total** | **19.587** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.150 s, not included above.
