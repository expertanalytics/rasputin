# rasputin mesh — statistics

`rasputin mesh --dem glo30 --cache /Users/skavhaug/projects/rasputin_data/germany_glo30_cache --out-crs EPSG:25832 --domain /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/fused_catchment_glo30.geojson --features /Users/skavhaug/projects/rasputin_data/germany_corine/isar_loisach_clc2018_3035.geojson --features-crs EPSG:3035 --features-map corine --features-tolerance 10 --tolerance 20 --start-min-angle 25 --out /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify/ft10_tol20_a25.vtk --stats /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify/ft10_tol20_a25_stats.md --record /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify/ft10_tol20_a25_record.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 1403 × 1724 (31 m) |
| domain vertices | 515 (1 ring, 0 holes) |
| start vertices | 62790 |
| start triangles | 124554 |
| output vertices | 147192 |
| output triangles | 292730 |
| constraint edges | 66973 |
| ft10_tol20_a25.vtk | 13.8 MB |

## Inputs

| item | value | name |
|---|---|---|
| The tiles or downloaded blocks used | Copernicus_DSM_COG_10_N47_00_E010_00_DEM; Copernicus_DSM_COG_10_N47_00_E011_00_DEM | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4326 | `dem_crs` |
| Conversion from the DEM's coordinate system | Inverse of ETRS89 to WGS 84 (1) + UTM zone 32N (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 31 m square grid in EPSG:25832, 1724 columns x 1403 rows | `resampled_grid` |
| The domain file and shape | fused_catchment_glo30.geojson: 1 outline, 0 holes, 515 vertices | `domain` |
| The domain's coordinate system | EPSG:25832 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| Each features file: layer, class map, features, lines, vertices | isar_loisach_clc2018_3035.geojson, class map corine: 18 features, 1157 lines, 125969 vertices | `features` |
| Each features file's coordinate system | EPSG:3035 | `features_crs` |
| Conversion from each features file's coordinate system | Inverse of Europe Equal Area 2001 + UTM zone 32N (with axis order normalized for visualization) | `features_transform` |
| Land-cover borders closer than this made one, narrower gaps filled | 1 | `features_repair_m` |
| Borders between land-cover polygons of one class dropped | on | `features_merge_same_class` |
| Land-cover borders simplified by this much (0 = off) | 10 | `features_tolerance_m` |
| Land-cover borders this close to the outline moved onto it | 5 | `features_outline_snap_m` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.89° | 0.00 % | 0.64 % | 0.489° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 13 | 5 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25832 | `crs` |
| Tolerance | 20 | `tolerance_m` |
| Largest height error | 19.99982279993128 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Land-cover vertices before and after clean-up | 221171 after the clip, 184903 after clean-up | `land_cover_vertices` |
| Land-cover area inside the outline that the outline rule gave to another polygon, m2 | 7987.3 | `land_cover_area_moved_m2` |
| Largest error against the resampled grid | 19.999589041267882 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 3326503 | `dem_nodes_checked` |
| Points that comparison added | 4286 | `dem_check_points_inserted` |
| Passes of that comparison | 7 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 414297 | `line_points_checked` |
| Their largest difference from the mesh | 19.998724747721326 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 141 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 14 | `refinement_rounds` |
| Points added | 15459 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 35526 | `edge_flips` |
| Points added to improve the starting mesh | 62017 | `start_quality_points_inserted` |
| Tries skipped while improving it | 29268 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 7770 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 659 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 1981 | `start_quality_lines_split` |
| Points refinement moved onto lines | 440 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 171 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 9 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 62790 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.0 % |
| decode | 0.308 | 8.4 % |
| features read | 0.598 | 16.3 % |
| features clip | 1.535 | 41.8 % |
| features clip: clean-up | 1.326 | 36.1 % |
| start mesh: build | 0.002 | 0.1 % |
| start mesh: node | 0.214 | 5.8 % |
| start mesh: triangulate | 0.054 | 1.5 % |
| start mesh: constraint edges | 0.129 | 3.5 % |
| refine | 0.214 | 5.8 % |
| refine: legalise start | 0.003 | 0.1 % |
| refine: start quality | 0.146 | 4.0 % |
| refine: scan (parallel) | 0.016 | 0.4 % |
| refine: split + flip (serial) | 0.010 | 0.3 % |
| refine: setup + output | 0.040 | 1.1 % |
| edge strip: generate | 0.016 | 0.4 % |
| check points: project | 0.140 | 3.8 % |
| check points: store | 0.055 | 1.5 % |
| final check: scan (parallel) | 0.065 | 1.8 % |
| final check: split + flip (serial) | 0.005 | 0.1 % |
| trim | 0.011 | 0.3 % |
| land cover | 0.129 | 3.5 % |
| write: encode | 0.005 | 0.1 % |
| write: disk | 0.003 | 0.1 % |
| other | 0.185 | 5.0 % |
| **total** | **3.671** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.052 s, not included above.
