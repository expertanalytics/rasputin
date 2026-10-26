# rasputin mesh — statistics

`rasputin mesh --out /Users/skavhaug/projects/rasputin_scratch/perf-106c-out/kultsjon-4000-ft200.gpkg --dem glo30 --out-crs EPSG:3006 --domain /Users/skavhaug/rasputin-output/opal-se2-2026-10-21/catchments/kultsjon.geojson --catchments smhi-svar2022-subcatchments --land-cover clc2018 --triangles 4000 --features-tolerance 200 --stats /Users/skavhaug/projects/rasputin_scratch/perf-106c-out/kultsjon-4000-ft200.stats.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| resampled grid nodes | 1577 × 2153 (31 m) |
| domain vertices | 195 (1 ring, 0 holes) |
| start vertices | 765 |
| start triangles | 1333 |
| output vertices | 2095 |
| output triangles | 3972 |
| constraint edges | 1187 |
| kultsjon-4000-ft200.gpkg | 2.4 MB |

## Inputs

| item | value | name |
|---|---|---|
| The tiles or downloaded blocks used | Copernicus_DSM_COG_10_N64_00_E014_00_DEM; Copernicus_DSM_COG_10_N64_00_E015_00_DEM; Copernicus_DSM_COG_10_N65_00_E014_00_DEM; Copernicus_DSM_COG_10_N65_00_E015_00_DEM | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres | `dem_vertical_unit` |
| The DEM's own coordinate system | EPSG:4326 | `dem_crs` |
| Conversion from the DEM's coordinate system | Inverse of SWEREF99 to WGS 84 (1) + SWEREF99 TM (with axis order normalized for visualization) | `dem_transform` |
| The grid the DEM was interpolated onto | 31 m square grid in EPSG:3006, 2153 columns x 1577 rows | `resampled_grid` |
| The domain file and shape | kultsjon.geojson: 1 outline, 0 holes, 195 vertices | `domain` |
| The domain's coordinate system | EPSG:3006 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 10 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 35.93° | 0.03 % | 1.46 % | 0.535° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 10 | 0 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:3006 | `crs` |
| Tolerance | 100.184504331562 | `tolerance_m` |
| Largest height error | 100.1521242057969 | `max_error_m` |
| DEM | glo30 | `dem_source` |
| Credit | produced using Copernicus WorldDEM-30 (c) DLR e.V. 2010-2014 and (c) Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved | `dem_credit` |
| Licence | Licence for Copernicus DEM instance COP-DEM-GLO-30-F; Art. 6(c): "The organisations in charge of the Copernicus programme by law or by delegation do not incur any liability for any use of the Copernicus WorldDEM-30" | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Land-cover classes: the source and the year | clc2018 | `land_cover_source` |
| Credit for the land-cover classes | CORINE Land Cover 2018 (vector/raster 100 m), Europe, 6-yearly. European Union's Copernicus Land Monitoring Service information, the European Environment Agency's map service https://image.discomap.eea.europa.eu/arcgis/rest/services/Corine/CLC2018_LAEA/MapServer/0 (accessed on 2026-10-17). Modified: the polygons meeting a box, stored as a GeoPackage by rasputin. Not endorsed by the European Union. | `land_cover_credit` |
| Licence of the land-cover classes | Copernicus Land Monitoring Service: full, open and free access (Commission Delegated Regulation (EU) No 1159/2013). Name the source, state any modification, and do not imply endorsement by the European Union. | `land_cover_licence_note` |
| What land_cover_code holds | clc2018 class number; 0 = no data, and every constraint line | `land_cover_codes` |
| Reservoir polygons splitting the cells' lake share | none | `reservoirs` |
| Land-cover polygons that crossed themselves, repaired | 43 | `land_cover_repaired` |
| Sub-catchment units: their borders constraints, their ids the cells' | smhi-svar2022-subcatchments | `catchments` |
| Sub-catchment units whose inner point lies in the domain | 75 | `catchments_units` |
| Band the sub-catchment borders were simplified within, m | 200 | `catchments_band_m` |
| Triangles asked for (about this many) | 4000 | `triangles_asked` |
| Triangles at most | 4200 | `triangles_cap` |
| Tolerances tried for the triangle budget, and the triangles each gave | 1e+06 m: 3498; 88.8 m: 4219; 22.2 m: 17414; 5.55 m: 136165; 171 m: 3564; 123 m: 3727; 100 m: 3972 | `budget_trials` |
| Root-mean-square height error | 25.975378043446067 | `rms_error_m` |
| DEM nodes inside the mesh the height errors are over | 2108226 | `error_nodes_measured` |
| Largest error against the resampled grid | 100.12968737410893 | `resampled_grid_max_error_m` |
| Nodes of the original DEM compared with the mesh | 3671533 | `dem_nodes_checked` |
| Points that comparison added | 11 | `dem_check_points_inserted` |
| Passes of that comparison | 3 | `dem_check_rounds` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 82160 | `line_points_checked` |
| Their largest difference from the mesh | 99.74306466749306 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 1 | `line_points_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 11 | `refinement_rounds` |
| Points added | 236 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 490 | `edge_flips` |
| Points added to improve the starting mesh | 765 | `start_quality_points_inserted` |
| Tries skipped while improving it | 555 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 478 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 15 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 303 | `start_quality_lines_split` |
| Points refinement moved onto lines | 29 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 765 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run per trial, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. The triangle budget ran 7 trials: the refine, edge strip, final check and trim rows sum all of them. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.0 % |
| catchments | 1.314 | 2.8 % |
| decode | 0.428 | 0.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.002 | 0.0 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.003 | 0.0 % |
| check points: project | 0.192 | 0.4 % |
| check points: store | 0.076 | 0.2 % |
| refine | 0.114 | 0.2 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.014 | 0.0 % |
| refine: scan (parallel) | 0.057 | 0.1 % |
| refine: split + flip (serial) | 0.034 | 0.1 % |
| refine: setup + output | 0.009 | 0.0 % |
| edge strip: generate | 0.017 | 0.0 % |
| final check: scan (parallel) | 0.153 | 0.3 % |
| final check: split + flip (serial) | 0.013 | 0.0 % |
| trim | 0.008 | 0.0 % |
| error summary | 0.011 | 0.0 % |
| land cover fractions | 45.113 | 94.8 % |
| write: geopackage | 0.052 | 0.1 % |
| other | 0.114 | 0.2 % |
| **total** | **47.612** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.001 s, not included above.
