# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain /Users/skavhaug/projects/rasputin_data/banenor_banenettverk/section_domain.geojson --tolerance 20 --tolerance-near /Users/skavhaug/projects/rasputin_data/banenor_banenettverk/bergensbanen.geojson 1 --tolerance-ramp 0 3000 --out /tmp/perf-104b-meshes/geilo-al-ramp-plain-r1.vtk --binary --stats docs/benchmarks/2026-10-25/104b/rms/stats/geilo-al-ramp-plain-r1.md --record docs/benchmarks/2026-10-25/104b/rms/stats/geilo-al-ramp-plain-r1.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 1970 × 2301 (10 m) |
| domain vertices | 52 (1 ring, 0 holes) |
| start vertices | 52 |
| start triangles | 50 |
| output vertices | 32188 |
| output triangles | 64162 |
| constraint edges | 212 |
| geilo-al-ramp-plain-r1.vtk | 2.3 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 2301 columns x 1970 rows, 10 m apart | `dem_grid` |
| The tiles or downloaded blocks used | 6701_2_10m_z33.tif; 6701_3_10m_z33.tif | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | none: the tiles agree where they overlap | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The domain file and shape | section_domain.geojson: 1 outline, 0 holes, 52 vertices | `domain` |
| The domain's coordinate system | EPSG:25833 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Tolerance on the tolerance lines | 1 | `tolerance_near_m` |
| Distance from the lines over which the tolerance rises to its far value | 0 to 3000 | `tolerance_ramp_m` |
| The tolerance lines file, and its segments after simplifying | bergensbanen.geojson: 724 segments | `tolerance_lines` |
| Threads used | 10 | `threads` |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.01 % | 0.86 % | 0.578° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 12 | 2 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 20 | `tolerance_m` |
| Largest height error | 19.984358261826856 | `max_error_m` |
| DEM | 6701_2_10m_z33.tif; 6701_3_10m_z33.tif | `dem_source` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Largest height error where the tolerance is the lines' own | 0.9996405707465783 | `max_error_near_lines_m` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 17794 | `line_points_checked` |
| Their largest difference from the mesh | 19.744482010954698 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 6 | `line_points_inserted` |
| DEM nodes the line check added | 0 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 34 | `refinement_rounds` |
| Points added | 32080 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 80205 | `edge_flips` |
| Points added to improve the starting mesh | 50 | `start_quality_points_inserted` |
| Tries skipped while improving it | 21 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 15 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 0 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 52 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 52 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| domain read | 0.001 | 0.4 % |
| decode | 0.140 | 41.7 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.000 | 0.0 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.000 | 0.0 % |
| refine | 0.100 | 29.9 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.083 | 24.8 % |
| refine: split + flip (serial) | 0.014 | 4.1 % |
| refine: setup + output | 0.003 | 1.0 % |
| edge strip: generate | 0.001 | 0.2 % |
| edge strip: scan (parallel) | 0.001 | 0.2 % |
| edge strip: split + flip (serial) | 0.000 | 0.1 % |
| trim | 0.002 | 0.7 % |
| write: encode | 0.004 | 1.2 % |
| write: disk | 0.001 | 0.2 % |
| other | 0.085 | 25.3 % |
| **total** | **0.335** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.014 s, not included above.
