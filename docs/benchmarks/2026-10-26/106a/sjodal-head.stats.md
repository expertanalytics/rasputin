# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/demo-fixes/build-bench/pkg --threads 0 -- mesh --out /Users/skavhaug/projects/rasputin_scratch/perf-106a-out/sjodal-head.gpkg --dem dtm10-norway --domain /Users/skavhaug/projects/rasputin_scratch/perf-106a-out/nedre-sjodalsvatn.geojson --land-cover ar50-innlandet --reservoirs nve-reservoirs --triangles 6000 --stats /Users/skavhaug/projects/rasputin_scratch/perf-106a-out/sjodal-head.stats.json`

bounds checks: on (libc++ fast)

## Sizes

| item | count |
|---|---|
| DEM nodes | 2214 × 4136 (10 m) |
| domain vertices | 1662 (1 ring, 0 holes) |
| start vertices | 1662 |
| start triangles | 1660 |
| output vertices | 3806 |
| output triangles | 5940 |
| constraint edges | 1670 |
| sjodal-head.gpkg | 3.5 MB |

## Inputs

| item | value | name |
|---|---|---|
| DEM grid size and spacing | 4136 columns x 2214 rows, 10 m apart | `dem_grid` |
| The tiles or downloaded blocks used | 6801_2_10m_z33.tif; 6801_3_10m_z33.tif | `dem_tiles` |
| Overlapping tiles that disagree, and by how much | 6801_2_10m_z33.tif and 6801_3_10m_z33.tif disagree at 38902 nodes, by up to 0.446899 m (median 0.00500488 m) | `dem_seams` |
| Unit of z | metres (assumed: the DEM file does not say) | `dem_vertical_unit` |
| The domain file and shape | nedre-sjodalsvatn.geojson: 1 outline, 0 holes, 1662 vertices | `domain` |
| The domain's coordinate system | EPSG:25833 | `domain_crs` |
| Conversion from the domain's coordinate system | none | `domain_transform` |
| What refinement started from | the domain outline | `start_mesh` |
| Starting mesh improved to this smallest angle (0 = off) | 25 | `start_min_angle_deg` |
| Starting mesh: a point is added only if it raises the smallest angle around it by at least this many degrees (negative = always added) | 0 | `start_quality_gain_deg` |
| Points very close to a line were moved onto it | on | `snap_to_lines` |
| Threads used | 10 | `threads` |

## DEM seams

| tile | tile | nodes | max | median |
|---|---|---|---|---|
| 6801_2_10m_z33.tif | 6801_3_10m_z33.tif | 38902 | 0.446899 | 0.00500488 |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 35.42° | 0.00 % | 0.27 % | 3.05° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 5 | 8 | 10 | 0 | 0 |

## Result

| item | value | name |
|---|---|---|
| Coordinate system | EPSG:25833 | `crs` |
| Tolerance | 545.9200599537053 | `tolerance_m` |
| Largest height error | 545.9200599537053 | `max_error_m` |
| DEM | 6801_2_10m_z33.tif; 6801_3_10m_z33.tif | `dem_source` |
| Credit | DTM10 \xa9 Kartverket | `dem_credit` |
| Licence | CC BY 4.0 | `licence_note` |
| Vertices removed on NoData | 0 | `nodata_vertices_removed` |
| Land-cover classes: the source and the year | ar50-innlandet | `land_cover_source` |
| Credit for the land-cover classes | Inneholder data under Norsk lisens for offentlige data (NLOD) tilgjengeliggjort av NIBIO. Kilde: NIBIO, AR50 Arealressurskart 1:50 000, Innlandet, lastet ned 2026-10-19; konvertert til GeoPackage av rasputin. | `land_cover_credit` |
| Licence of the land-cover classes | Norsk lisens for offentlige data (NLOD) 1.0, https://data.norge.no/nlod/no/1.0. | `land_cover_licence_note` |
| What land_cover_code holds | ar50-innlandet class number; 0 = no data, and every constraint line | `land_cover_codes` |
| Reservoir polygons splitting the cells' lake share | nve-reservoirs | `reservoirs` |
| Land-cover polygons that crossed themselves, repaired | 32 | `land_cover_repaired` |
| Reservoir polygons that crossed themselves, repaired | 0 | `reservoirs_repaired` |
| Triangles asked for (about this many) | 6000 | `triangles_asked` |
| Triangles at most | 6300 | `triangles_cap` |
| Tolerances tried for the triangle budget, and the triangles each gave | 1e+06 m: 5940 | `budget_trials` |
| Root-mean-square height error | 113.11096175286805 | `rms_error_m` |
| DEM nodes inside the mesh the height errors are over | 4803976 | `error_nodes_measured` |
| DEM nodes on a vertex (to rounding), not compared | 0 | `dem_nodes_at_vertices` |
| Their largest difference | 0 | `dem_nodes_at_vertices_max_error_m` |
| Points along lines (grid-line crossings and halfway between) compared | 38734 | `line_points_checked` |
| Their largest difference from the mesh | 96.27915635986437 | `line_max_error_m` |
| Points along lines on NoData cells, not compared | 0 | `line_points_on_nodata` |
| Points along lines that could not be added | 0 | `line_points_refused` |
| Their largest difference from the mesh | 0 | `line_points_refused_max_error_m` |
| Points along lines added | 0 | `line_points_inserted` |
| DEM nodes the line check added | 0 | `line_check_dem_nodes_inserted` |
| Points along lines dropped as duplicates | 0 | `line_points_duplicate` |
| Self-check: DEM nodes with data left outside the mesh | 0 | `dem_nodes_outside_mesh` |
| Triangles flat or turned over at the written coordinates; their edges were left out of the last edge-swap check | 0 | `triangles_flat_as_written` |
| Edges that fail the empty-circle test at the written coordinates and could not be swapped | 0 | `edges_not_delaunay_as_written` |
| Refinement passes | 1 | `refinement_rounds` |
| Points added | 0 | `points_inserted` |
| Of them, into triangles with a NoData corner | 0 | `points_inserted_on_nodata` |
| Edge swaps | 0 | `edge_flips` |
| Points added to improve the starting mesh | 2136 | `start_quality_points_inserted` |
| Tries skipped while improving it | 780 | `start_quality_points_skipped` |
| Tries that would not have improved the angles | 692 | `start_quality_points_without_gain` |
| Points moved onto lines while improving the starting mesh | 8 | `start_quality_points_snapped_to_lines` |
| Land-cover and outline lines split to improve angles | 0 | `start_quality_lines_split` |
| Points refinement moved onto lines | 0 | `points_snapped_to_lines` |
| Moves onto lines refused | 0 | `snaps_refused` |
| Points the final check against the DEM moved onto lines | 0 | `final_check_points_snapped_to_lines` |
| Of them, also added where they were | 0 | `final_check_snapped_points_added_anyway` |
| Starting-mesh vertices not on a DEM node | 1662 | `start_vertices_between_dem_nodes` |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run per trial, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. The triangle budget ran 1 trial: the refine, edge strip, final check and trim rows sum all of them. Threads used: 10.

| phase | seconds | share |
|---|---|---|
| domain read | 0.003 | 0.1 % |
| decode | 0.181 | 5.6 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.002 | 0.0 % |
| start mesh: triangulate | 0.001 | 0.0 % |
| start mesh: constraint edges | 0.002 | 0.1 % |
| refine | 0.010 | 0.3 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.006 | 0.2 % |
| refine: scan (parallel) | 0.003 | 0.1 % |
| refine: split + flip (serial) | 0.000 | 0.0 % |
| refine: setup + output | 0.001 | 0.0 % |
| edge strip: generate | 0.001 | 0.0 % |
| edge strip: scan (parallel) | 0.000 | 0.0 % |
| edge strip: split + flip (serial) | 0.000 | 0.0 % |
| trim | 0.000 | 0.0 % |
| error summary | 0.005 | 0.2 % |
| land cover fractions | 2.492 | 77.9 % |
| reservoirs | 0.412 | 12.9 % |
| write: geopackage | 0.074 | 2.3 % |
| other | 0.018 | 0.6 % |
| **total** | **3.200** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.002 s, not included above.
