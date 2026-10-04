# rasputin mesh — statistics

`bench.py _child --pkg /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/base-c193cb1/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/docs/benchmarks/2026-09-29/bygdin-landcover/bygdin_reduced_t20.geojson --tolerance 1 --features /Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/af318b82-b96d-429e-956c-0296b113f23b/scratchpad/bygdin/base_t1_b1_run2.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/15f-3/docs/benchmarks/2026-10-04/15f-3-acceptance/bygdin/base_t1_b1_run2.stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 8242 |
| start triangles | 15614 |
| output vertices | 573758 |
| output triangles | 1145434 |
| constraint edges | 14074 |
| base_t1_b1_run2.vtk | 48.5 MB |

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
| minimum angle | 45.00° | 0.10 % | 1.36 % | 0.00408° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 29 | 1377 | 37 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 1 | `tolerance_m` |
| Largest height error | 1 | `max_error_m` |
| DEM | 6801_2_10m_z33.tif; 6801_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 24 | `refinement_rounds` |
| Points added | 549527 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 1152221 | `edge_flips` |
| Points added to improve the starting mesh | 15989 | `start_quality_points_inserted` |
| Tries skipped while improving it | 2442 | `start_quality_points_skipped` |
| Points moved onto lines | 5726 | `points_snapped_to_lines` |
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
| decode | 0.170 | 4.5 % |
| features read | 0.135 | 3.5 % |
| features clip | 1.482 | 38.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.017 | 0.5 % |
| start mesh: triangulate | 0.006 | 0.1 % |
| start mesh: constraint edges | 0.017 | 0.5 % |
| refine | 0.595 | 15.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.016 | 0.4 % |
| refine: scan (parallel) | 0.099 | 2.6 % |
| refine: split + flip (serial) | 0.426 | 11.2 % |
| refine: setup + output | 0.054 | 1.4 % |
| trim | 0.042 | 1.1 % |
| land cover | 1.296 | 34.0 % |
| write: encode | 0.018 | 0.5 % |
| write: disk | 0.022 | 0.6 % |
| other | 0.005 | 0.1 % |
| **total** | **3.808** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.210 s, not included above.
