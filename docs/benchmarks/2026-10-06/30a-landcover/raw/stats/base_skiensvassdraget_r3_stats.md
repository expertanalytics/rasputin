# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin_scratch/norway/skiensvassdraget/skiensvassdraget_outline_nve.geojson --features /Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --tolerance 10 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/30a/base_skiensvassdraget_r3.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/landcover-speed/docs/benchmarks/2026-10-06/30a-landcover/raw/stats/base_skiensvassdraget_r3_stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 13675 × 14653 (10 m) |
| domain vertices | 16007 (1 ring, 0 holes) |
| start vertices | 256482 |
| start triangles | 496210 |
| output vertices | 1523287 |
| output triangles | 3029295 |
| constraint edges | 275572 |
| base_skiensvassdraget_r3.vtk | 133.6 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 14653 columns x 13675 rows, 10 m apart | `dem_grid` |
| The tiles or downloaded blocks used | 6500_1_10m_z33.tif; 6501_1_10m_z33.tif; 6501_4_10m_z33.tif; 6600_1_10m_z33.tif; 6600_2_10m_z33.tif; 6601_1_10m_z33.tif; 6601_2_10m_z33.tif; 6601_3_10m_z33.tif; 6601_4_10m_z33.tif; 6700_2_10m_z33.tif; 6701_2_10m_z33.tif; 6701_3_10m_z33.tif | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The domain file and shape | skiensvassdraget_outline_nve.geojson: 1 outline, 0 holes, 16007 vertices | `domain` |
| The domain's coordinate system | EPSG:25833 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| Each features file: layer, class map, features, lines, vertices | corine2018_dtm10_utm33.gpkg layer corine2018, class map corine: 2628 features, 4138 lines, 484504 vertices | `features` |
| Each features file's coordinate system | EPSG:25833 | `features_crs` |
| Conversion from each features file's coordinate system | none | `features_transform` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| DEM nodes very close to a line were moved onto it | on | `snap_to_lines` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.11 % | 2.09 % | 3.15e-05° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 20 | 733 | 2 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999987707205037 | `max_error_m` |
| DEM | 6500_1_10m_z33.tif; 6501_1_10m_z33.tif; 6501_4_10m_z33.tif; 6600_1_10m_z33.tif; 6600_2_10m_z33.tif; 6601_1_10m_z33.tif; 6601_2_10m_z33.tif; 6601_3_10m_z33.tif; 6601_4_10m_z33.tif; 6700_2_10m_z33.tif; 6701_2_10m_z33.tif; 6701_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 5027869 | `line_points_checked` |
| Their largest difference from the mesh | 9.999998555189222 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 2979 | `line_points_inserted` |
| DEM nodes the line check added | 483 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 26 | `refinement_rounds` |
| Points added | 714270 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 1734812 | `edge_flips` |
| Points added to improve the starting mesh | 549073 | `start_quality_points_inserted` |
| Tries skipped while improving it | 125398 | `start_quality_points_skipped` |
| Points moved onto lines | 14024 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Starting-mesh vertices not on a DEM node | 256482 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.015 | 0.1 % |
| decode | 1.289 | 6.3 % |
| features read | 0.267 | 1.3 % |
| features clip | 5.649 | 27.6 % |
| start mesh: build | 0.008 | 0.0 % |
| start mesh: node | 1.014 | 5.0 % |
| start mesh: triangulate | 0.267 | 1.3 % |
| start mesh: constraint edges | 0.616 | 3.0 % |
| refine | 2.908 | 14.2 % |
| refine: legalise start | 0.013 | 0.1 % |
| refine: start quality | 1.013 | 5.0 % |
| refine: scan (parallel) | 0.831 | 4.1 % |
| refine: split + flip (serial) | 0.775 | 3.8 % |
| refine: setup + output | 0.275 | 1.3 % |
| edge strip: generate | 0.205 | 1.0 % |
| edge strip: scan (parallel) | 0.027 | 0.1 % |
| edge strip: split + flip (serial) | 0.025 | 0.1 % |
| trim | 0.117 | 0.6 % |
| land cover | 6.539 | 32.0 % |
| write: encode | 0.043 | 0.2 % |
| write: disk | 0.074 | 0.4 % |
| other | 1.394 | 6.8 % |
| **total** | **20.458** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.606 s, not included above.
