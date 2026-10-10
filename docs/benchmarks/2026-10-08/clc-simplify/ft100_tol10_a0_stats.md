# rasputin mesh — statistics

`rasputin mesh --dem glo30 --cache /Users/skavhaug/projects/rasputin_data/germany_glo30_cache --out-crs EPSG:25832 --domain /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/fused_catchment_glo30.geojson --features /Users/skavhaug/projects/rasputin_data/germany_corine/isar_loisach_clc2018_3035.geojson --features-crs EPSG:3035 --features-map corine --features-tolerance 100 --tolerance 10 --start-min-angle 0 --out /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify/ft100_tol10_a0.vtk --stats /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify/ft100_tol10_a0_stats.md --record /Users/skavhaug/projects/rasputin_scratch/germany/isar-loisach/clc_simplify/ft100_tol10_a0_record.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 1403 × 1724 (31 m) |
| domain vertices | 515 (1 ring, 0 holes) |
| start vertices | 13629 |
| start triangles | 26315 |
| output vertices | 144050 |
| output triangles | 285950 |
| constraint edges | 27625 |
| ft100_tol10_a0.vtk | 12.7 MB |

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
| Each features file: layer, class map, features, lines, vertices | isar_loisach_clc2018_3035.geojson, class map corine: 18 features, 1060 lines, 27599 vertices | `features` |
| Each features file's coordinate system | EPSG:3035 | `features_crs` |
| Conversion from each features file's coordinate system | Inverse of Europe Equal Area 2001 + UTM zone 32N (with axis order normalized for visualization) | `features_transform` |
| Land-cover borders closer than this made one, narrower gaps filled | 1 | `features_repair_m` |
| Borders between land-cover polygons of one class dropped | on | `features_merge_same_class` |
| Land-cover borders simplified by this much (0 = off) | 100 | `features_tolerance_m` |
| Land-cover borders this close to the outline moved onto it | 5 | `features_outline_snap_m` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 0 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 33.69° | 0.02 % | 2.98 % | 0.0416° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 10 | 14 | 158 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25832 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999952936927002 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Land-cover vertices before and after clean-up | 221171 after the clip, 41066 after clean-up | `land_cover_vertices` |
| Land-cover area inside the outline that the outline rule gave to another polygon, m2 | 8101.9 | `land_cover_area_moved_m2` |
| Largest error against the resampled grid | 9.999885481154365 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 3326503 | `dem_nodes_checked` |
| Points that comparison added | 29610 | `dem_check_points_inserted` |
| Passes of that comparison | 13 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 330992 | `line_points_checked` |
| Their largest difference from the mesh | 9.996563819618132 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 1845 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 20 | `refinement_rounds` |
| Points added | 100811 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 222847 | `edge_flips` |
| Points added to improve the starting mesh | 0 | `start_quality_points_inserted` |
| Tries skipped while improving it | 0 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 0 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 8736 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 2658 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 238 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 13629 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.1 % |
| decode | 0.308 | 10.2 % |
| features read | 0.602 | 19.9 % |
| features clip | 1.312 | 43.4 % |
| features clip: clean-up | 1.100 | 36.4 % |
| start mesh: build | 0.001 | 0.0 % |
| start mesh: node | 0.048 | 1.6 % |
| start mesh: triangulate | 0.010 | 0.3 % |
| start mesh: constraint edges | 0.029 | 0.9 % |
| refine | 0.099 | 3.3 % |
| refine: legalise start | 0.001 | 0.0 % |
| refine: start quality | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.027 | 0.9 % |
| refine: split + flip (serial) | 0.055 | 1.8 % |
| refine: setup + output | 0.015 | 0.5 % |
| edge strip: generate | 0.012 | 0.4 % |
| check points: project | 0.147 | 4.9 % |
| check points: store | 0.057 | 1.9 % |
| final check: scan (parallel) | 0.083 | 2.7 % |
| final check: split + flip (serial) | 0.028 | 0.9 % |
| trim | 0.010 | 0.3 % |
| land cover | 0.121 | 4.0 % |
| write: encode | 0.005 | 0.2 % |
| write: disk | 0.002 | 0.1 % |
| other | 0.149 | 4.9 % |
| **total** | **3.024** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.051 s, not included above.
