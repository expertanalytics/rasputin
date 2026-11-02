# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin/.claude/worktrees/lakes-overrule/docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson --tolerance 10 --edges all --threads 8 --prefer lake=/Users/skavhaug/projects/rasputin_scratch/temp-check-20261021/nve-lakes-25833.geojson --land-cover ar50-innlandet-clean-25833 --land-cover ar50-vestland-clean-25833 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/out/byg-innvest-master.vtk --record /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/out/byg-innvest-master.json --stats /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/out/byg-innvest-master.md`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 31785 |
| start triangles | 62623 |
| output vertices | 63582 |
| output triangles | 125973 |
| constraint edges | 34073 |
| byg-innvest-master.vtk | 6.1 MB |

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
| Each source of edges: file, layer, land-cover class file or line kind, features, lines, vertices | ar50-innlandet-clean-25833 + ar50-vestland-clean-25833 layer ArealressursFlate, classes ar50: 8 features, 1275 lines, 63817 vertices | `edges` |
| Lines and vertices of all edges, joined | 1275 lines, 62773 vertices after joining shared edges and adding crossings | `edges_joined` |
| Each land-cover file's coordinate system | EPSG:25833; EPSG:25833 | `land_cover_crs` |
| Conversion from each land-cover file's coordinate system | none; none | `land_cover_transform` |
| Area inside the domain where land-cover files overlap, and the file that kept it | none | `land_cover_overlap` |
| Land-cover borders closer than this made one, narrower gaps filled | 0.1 | `edges_repair_m` |
| Borders between land-cover polygons of one class dropped | on | `edges_merge_same_class` |
| Land-cover borders simplified, each at most this far from its border after the repair and the outline rule, each class's area kept (0 = off) | 5 | `edges_max_shift_m` |
| Land-cover borders this close to the outline moved onto it | 0.5 | `edges_outline_snap_m` |
| What refinement started from | the domain outline and the feature lines | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 15 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 8 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 31.55° | 0.01 % | 0.72 % | 0.296° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 13 | 12 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 10 | `tolerance_m` |
| Largest height error | 9.997639973958485 | `max_error_m` |
| DEM | 6801_2_10m_z33.tif; 6801_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Credit for the classes taken from a preferred source | lake: nve-lakes-25833.geojson | `land_cover_prefer_credit` |
| Credit for the land-cover data | Contains data under the Norwegian licence for Open Government Data (NLOD) distributed by NIBIO: AR50 land-resource map 1:50 000, re-encoded | `land_cover_notice` |
| Licence of the land-cover data | Norwegian licence for Open Government Data (NLOD) 1.0, https://data.norge.no/nlod/en/1.0 | `land_cover_licence` |
| What land_cover_code holds | AR50 arealtype code, attribute arealtype, class 30 split by treslag into 3031, 3032, 3033, 3099, class 50 split by vegetasjonsdekke into 5051, 5052, 5054, 5055, 5099, legend ar50; 0 = in no polygon, and every constraint line | `land_cover_codes` |
| Class shares per triangle | not computed (add --shares) | `land_cover_shares` |
| Classes taken from a preferred source: inside the mesh, the land cover's area of the class, the preferred source's, and what went to neighbouring classes | lake from nve-lakes-25833.geojson: 55.1328 km2 replaced, 55.3325 km2 given, 0.7604 km2 refilled | `land_cover_prefer` |
| Land-cover vertices before and after clean-up | 259920 after the clip, 63768 after clean-up | `land_cover_vertices` |
| Land-cover area inside the outline that the outline rule gave to another polygon, m2 | 25.7 | `land_cover_area_moved_m2` |
| DEM nodes on a vertex (to rounding), not compared | 9 | `dem_nodes_at_vertices` |
| Their largest difference | 5.775291356258094e-11 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 372920 | `line_points_checked` |
| Their largest difference from the mesh | 9.968105575050231 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 83 | `line_points_inserted` |
| DEM nodes the line check added | 25 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 1273 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 18 | `refinement_rounds` |
| Points added | 16118 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 39888 | `edge_flips` |
| Points added to improve the starting mesh | 14090 | `start_quality_points_inserted` |
| Tries skipped while improving it | 1972 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 1239 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 137 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 1342 | `start_quality_lines_split` |
| Points refinement moved onto lines | 283 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 2 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 27921 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads used: 8.

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.0 % |
| decode | 0.159 | 2.1 % |
| land cover read | 3.282 | 42.5 % |
| edges | 2.680 | 34.7 % |
| edges: clean-up | 2.678 | 34.7 % |
| start mesh: build | 0.002 | 0.0 % |
| start mesh: node | 0.101 | 1.3 % |
| start mesh: triangulate | 0.026 | 0.3 % |
| start mesh: constraint edges | 0.087 | 1.1 % |
| refine | 0.098 | 1.3 % |
| refine: legalise start | 0.002 | 0.0 % |
| refine: start quality | 0.043 | 0.6 % |
| refine: scan (parallel) | 0.021 | 0.3 % |
| refine: split + flip (serial) | 0.009 | 0.1 % |
| refine: setup + output | 0.022 | 0.3 % |
| edge strip: generate | 0.013 | 0.2 % |
| edge strip: scan (parallel) | 0.003 | 0.0 % |
| edge strip: split + flip (serial) | 0.005 | 0.1 % |
| trim | 0.006 | 0.1 % |
| land cover | 0.055 | 0.7 % |
| write: encode | 0.004 | 0.1 % |
| write: disk | 0.001 | 0.0 % |
| other | 1.193 | 15.5 % |
| **total** | **7.716** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.029 s, not included above.
