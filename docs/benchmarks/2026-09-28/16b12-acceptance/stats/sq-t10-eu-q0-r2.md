# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/square48_25833.geojson --tolerance 10 --start-min-angle 0 --features /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/sq-t10-eu-q0.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/sq-t10-eu-q0-r2.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 4801 × 4801 (10 m) |
| domain vertices | 4 (1 ring, 0 holes) |
| start vertices | 51395 |
| start triangles | 102627 |
| output vertices | 106401 |
| output triangles | 212255 |
| constraint edges | 52531 |
| vertices without data dropped | 0 |
| sq-t10-eu-q0.vtk | 12.9 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 28.59° | 3.45 % | 13.63 % | 0.0147° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 11 | 25 | 540 | 9 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999949490582509 m | 23 | 55006 | 0 | 192561 | 0 | 0 | 0 | 359 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.005 | 0.0 % |
| decode | 0.159 | 0.7 % |
| features read | 14.345 | 58.7 % |
| features clip | 8.652 | 35.4 % |
| start mesh: build | 0.002 | 0.0 % |
| start mesh: node | 0.151 | 0.6 % |
| start mesh: triangulate | 0.045 | 0.2 % |
| start mesh: constraint edges | 0.121 | 0.5 % |
| refine | 0.340 | 1.4 % |
| refine: legalise start | 0.004 | 0.0 % |
| refine: start quality | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.229 | 0.9 % |
| refine: split + flip (serial) | 0.040 | 0.2 % |
| refine: setup + output | 0.067 | 0.3 % |
| trim | 0.009 | 0.0 % |
| write: encode | 0.599 | 2.4 % |
| write: disk | 0.002 | 0.0 % |
| other | 0.017 | 0.1 % |
| **total** | **24.447** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.041 s, not included above.
