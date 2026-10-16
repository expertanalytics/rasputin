# rasputin mesh — statistics

`bench.py _child --pkg /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/base-c193cb1/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/docs/benchmarks/2026-09-29/bygdin-landcover/bygdin_reduced_t20.geojson --tolerance 10 --features /Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/bygdin/base_t10_b1_run2.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/docs/benchmarks/2026-10-04/15f-3-acceptance/bygdin/base_t10_b1_run2.stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 8242 |
| start triangles | 15614 |
| output vertices | 42110 |
| output triangles | 83171 |
| constraint edges | 8882 |
| base_t10_b1_run2.vtk | 3.7 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 3314 columns x 2195 rows, 10 m apart | `dem_grid` |
| The tiles or downloaded blocks used | 6801_2_10m_z33.tif; 6801_3_10m_z33.tif | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The domain file and shape | bygdin_reduced_t20.geojson: 1 outline, 0 holes, 740 vertices | `domain` |
| The domain's coordinate system | EPSG:25833 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| Each features file: layer, class map, features, lines, vertices | corine2018_dtm10_utm33.gpkg layer corine2018, class map corine: 79 features, 208 lines, 15106 vertices | `features` |
| Each features file's coordinate system | EPSG:25833 | `features_crs` |
| Conversion from each features file's coordinate system | none | `features_transform` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| DEM nodes very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.07 % | 1.70 % | 0.00408° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 16 | 31 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999732907715497 | `max_error_m` |
| DEM | 6801_2_10m_z33.tif; 6801_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 17 | `refinement_rounds` |
| Points added | 17879 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 39582 | `edge_flips` |
| Points added to improve the starting mesh | 15989 | `start_quality_points_inserted` |
| Tries skipped while improving it | 2442 | `start_quality_points_skipped` |
| Points moved onto lines | 534 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Starting-mesh vertices not on a DEM node | 8242 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.1 % |
| decode | 0.176 | 9.0 % |
| features read | 0.131 | 6.7 % |
| features clip | 1.475 | 75.5 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.018 | 0.9 % |
| start mesh: triangulate | 0.005 | 0.3 % |
| start mesh: constraint edges | 0.016 | 0.8 % |
| refine | 0.048 | 2.4 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.015 | 0.8 % |
| refine: scan (parallel) | 0.014 | 0.7 % |
| refine: split + flip (serial) | 0.009 | 0.5 % |
| refine: setup + output | 0.009 | 0.4 % |
| trim | 0.003 | 0.2 % |
| land cover | 0.071 | 3.6 % |
| write: encode | 0.001 | 0.1 % |
| write: disk | 0.002 | 0.1 % |
| other | 0.005 | 0.2 % |
| **total** | **1.955** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.015 s, not included above.
