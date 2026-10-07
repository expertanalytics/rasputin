# rasputin mesh — statistics

`drive.py --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/build-bench/pkg --json /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/benchmarks/2026-10-07/20c-1/raw/numed-nofeet-r7.json -- mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin_scratch/norway/numedalslagen/numedalslagen_outline_nve.geojson --features /Users/skavhaug/projects/rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --tolerance 10 --binary --no-constraint-feet --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a/scratchpad/perf-20c1/meshes/numed-nofeet.vtk --stats /Users/skavhaug/projects/rasputin/.claude/worktrees/soft-quality/docs/benchmarks/2026-10-07/20c-1/raw/numed-nofeet-r7_stats.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 17882 × 15667 (10 m) |
| domain vertices | 14092 (1 ring, 0 holes) |
| start vertices | 149930 |
| start triangles | 284938 |
| output vertices | 659212 |
| output triangles | 1303047 |
| constraint edges | 156567 |
| numed-nofeet.vtk | 58.4 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 15667 columns x 17882 rows, 10 m apart | `dem_grid` |
| The tiles or downloaded blocks used | 6500_1_10m_z33.tif; 6501_1_10m_z33.tif; 6501_4_10m_z33.tif; 6502_4_10m_z33.tif; 6600_1_10m_z33.tif; 6600_2_10m_z33.tif; 6601_1_10m_z33.tif; 6601_2_10m_z33.tif; 6601_3_10m_z33.tif; 6601_4_10m_z33.tif; 6602_3_10m_z33.tif; 6602_4_10m_z33.tif; 6700_2_10m_z33.tif; 6701_2_10m_z33.tif; 6701_3_10m_z33.tif; 6702_3_10m_z33.tif | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | 6601_1_10m_z33.tif and 6602_3_10m_z33.tif disagree at 1804 nodes, by up to 0.043335 m (median 0.0050354 m); 6601_1_10m_z33.tif and 6602_4_10m_z33.tif disagree at 8379 nodes, by up to 1.18677 m (median 0.0448608 m); 6601_2_10m_z33.tif and 6602_3_10m_z33.tif disagree at 56278 nodes, by up to 0.207977 m (median 0.00482178 m); 6601_2_10m_z33.tif and 6602_4_10m_z33.tif disagree at 2026 nodes, by up to 1.18677 m (median 0.0470276 m); 6602_3_10m_z33.tif and 6602_4_10m_z33.tif disagree at 2029 nodes, by up to 1.19373 m (median 0.0465088 m); 6700_2_10m_z33.tif and 6701_3_10m_z33.tif disagree at 31004 nodes, by up to 0.267212 m (median 0.0020752 m) | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The domain file and shape | numedalslagen_outline_nve.geojson: 1 outline, 0 holes, 14092 vertices | `domain` |
| The domain's coordinate system | EPSG:25833 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| Each features file: layer, class map, features, lines, vertices | corine2018_dtm10_utm33.gpkg layer corine2018, class map corine: 1611 features, 2765 lines, 273953 vertices | `features` |
| Each features file's coordinate system | EPSG:25833 | `features_crs` |
| Conversion from each features file's coordinate system | none | `features_transform` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Points very close to a line were moved onto it | off | `snap_to_lines` |

## DEM seams

| tile | tile | nodes | max | median |
|---|---|---|---|---|
| 6601_1_10m_z33.tif | 6602_3_10m_z33.tif | 1804 | 0.043335 | 0.0050354 |
| 6601_1_10m_z33.tif | 6602_4_10m_z33.tif | 8379 | 1.18677 | 0.0448608 |
| 6601_2_10m_z33.tif | 6602_3_10m_z33.tif | 56278 | 0.207977 | 0.00482178 |
| 6601_2_10m_z33.tif | 6602_4_10m_z33.tif | 2026 | 1.18677 | 0.0470276 |
| 6602_3_10m_z33.tif | 6602_4_10m_z33.tif | 2029 | 1.19373 | 0.0465088 |
| 6700_2_10m_z33.tif | 6701_3_10m_z33.tif | 31004 | 0.267212 | 0.0020752 |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.49 % | 3.49 % | 7.22e-06° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 19 | 487 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999971563165786 | `max_error_m` |
| DEM | 6500_1_10m_z33.tif; 6501_1_10m_z33.tif; 6501_4_10m_z33.tif; 6502_4_10m_z33.tif; 6600_1_10m_z33.tif; 6600_2_10m_z33.tif; 6601_1_10m_z33.tif; 6601_2_10m_z33.tif; 6601_3_10m_z33.tif; 6601_4_10m_z33.tif; 6602_3_10m_z33.tif; 6602_4_10m_z33.tif; 6700_2_10m_z33.tif; 6701_2_10m_z33.tif; 6701_3_10m_z33.tif; 6702_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the features data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `features_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, map corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 2941574 | `line_points_checked` |
| Their largest difference from the mesh | 9.999588520580346 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 5173 | `line_points_inserted` |
| DEM nodes the line check added | 385 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Refinement passes | 21 | `refinement_rounds` |
| Points added | 187233 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 441327 | `edge_flips` |
| Points added to improve the starting mesh | 316491 | `start_quality_points_inserted` |
| Tries skipped while improving it | 79591 | `start_quality_points_skipped` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Points refinement moved onto lines | 0 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 149930 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.013 | 0.1 % |
| decode | 1.742 | 15.5 % |
| features read | 0.252 | 2.2 % |
| features clip | 5.735 | 51.0 % |
| start mesh: build | 0.005 | 0.0 % |
| start mesh: node | 0.474 | 4.2 % |
| start mesh: triangulate | 0.141 | 1.3 % |
| start mesh: constraint edges | 0.343 | 3.1 % |
| refine | 1.065 | 9.5 % |
| refine: legalise start | 0.008 | 0.1 % |
| refine: start quality | 0.498 | 4.4 % |
| refine: scan (parallel) | 0.285 | 2.5 % |
| refine: split + flip (serial) | 0.143 | 1.3 % |
| refine: setup + output | 0.131 | 1.2 % |
| edge strip: generate | 0.110 | 1.0 % |
| edge strip: scan (parallel) | 0.016 | 0.1 % |
| edge strip: split + flip (serial) | 0.021 | 0.2 % |
| trim | 0.049 | 0.4 % |
| land cover | 0.676 | 6.0 % |
| write: encode | 0.021 | 0.2 % |
| write: disk | 0.012 | 0.1 % |
| other | 0.569 | 5.1 % |
| **total** | **11.246** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.248 s, not included above.
