# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin/.claude/worktrees/lakes-overrule/docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance 10 --edges all --threads 8 --prefer lake=/Users/skavhaug/projects/rasputin_scratch/temp-check-20261021/nve-lakes-25833.geojson --land-cover clc2018-clean-25833 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/out/byg-clc-master.vtk --record /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/out/byg-clc-master.json --stats /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/out/byg-clc-master.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 4217 |
| start triangles | 7565 |
| output vertices | 30367 |
| output triangles | 59623 |
| constraint edges | 6662 |
| byg-clc-master.vtk | 2.7 MB |

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
| Each source of edges: file, layer, land-cover class file or line kind, features, lines, vertices | clc2018-clean-25833 layer features, classes corine: 7 features, 788 lines, 7958 vertices | `edges` |
| Lines and vertices of all edges, joined | 788 lines, 7307 vertices after joining shared edges and adding crossings | `edges_joined` |
| Each land-cover file's coordinate system | EPSG:25833 | `land_cover_crs` |
| Conversion from each land-cover file's coordinate system | none | `land_cover_transform` |
| Land-cover borders closer than this made one, narrower gaps filled | 1 | `edges_repair_m` |
| Borders between land-cover polygons of one class dropped | on | `edges_merge_same_class` |
| Land-cover borders simplified, each at most this far from its border after the repair and the outline rule, each class's area kept (0 = off) | 50 | `edges_max_shift_m` |
| Land-cover borders this close to the outline moved onto it | 5 | `edges_outline_snap_m` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 15 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 8 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 32.52° | 0.01 % | 3.19 % | 0.766° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 10 | 16 | 38 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.999755859375 | `max_error_m` |
| DEM | 6801_2_10m_z33.tif; 6801_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the classes taken from a preferred source | lake: nve-lakes-25833.geojson | `land_cover_prefer_credit` |
| Credit for the land-cover data | Contains modified CORINE Land Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land Monitoring Service 2018, European Environment Agency (EEA): clipped and re-encoded, produced with funding by the European Union, not endorsed by the EU | `land_cover_notice` |
| What land_cover_code holds | CORINE Land Cover level-3 code, attribute Code_18, legend corine; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Class shares per triangle | not computed (add --shares) | `land_cover_shares` |
| Classes taken from a preferred source: inside the mesh, the land cover's area of the class, the preferred source's, and what went to neighbouring classes | lake from nve-lakes-25833.geojson: 50.1911 km2 replaced, 55.3325 km2 given, 0.5717 km2 refilled | `land_cover_prefer` |
| Land-cover vertices before and after clean-up | 85853 after the clip, 7942 after clean-up | `land_cover_vertices` |
| Land-cover area inside the outline that the outline rule gave to another polygon, m2 | 1905.0 | `land_cover_area_moved_m2` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 201250 | `line_points_checked` |
| Their largest difference from the mesh | 9.990996754614343 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 171 | `line_points_inserted` |
| DEM nodes the line check added | 48 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 1 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 20 | `refinement_rounds` |
| Points added | 23194 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 53558 | `edge_flips` |
| Points added to improve the starting mesh | 1823 | `start_quality_points_inserted` |
| Tries skipped while improving it | 622 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 354 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 36 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 873 | `start_quality_lines_split` |
| Points refinement moved onto lines | 1129 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 5 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 4217 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads used: 8.

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.1 % |
| decode | 0.323 | 9.6 % |
| land cover read | 1.516 | 44.9 % |
| edges | 0.984 | 29.2 % |
| edges: clean-up | 0.983 | 29.1 % |
| start mesh: build | 0.001 | 0.0 % |
| start mesh: node | 0.012 | 0.3 % |
| start mesh: triangulate | 0.003 | 0.1 % |
| start mesh: constraint edges | 0.010 | 0.3 % |
| refine | 0.046 | 1.4 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.005 | 0.1 % |
| refine: scan (parallel) | 0.023 | 0.7 % |
| refine: split + flip (serial) | 0.013 | 0.4 % |
| refine: setup + output | 0.006 | 0.2 % |
| edge strip: generate | 0.008 | 0.2 % |
| edge strip: scan (parallel) | 0.002 | 0.1 % |
| edge strip: split + flip (serial) | 0.002 | 0.1 % |
| trim | 0.003 | 0.1 % |
| land cover | 0.023 | 0.7 % |
| write: encode | 0.002 | 0.1 % |
| write: disk | 0.001 | 0.0 % |
| other | 0.436 | 12.9 % |
| **total** | **3.374** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.014 s, not included above.
