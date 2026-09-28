# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/square48_25833.geojson --tolerance 1 --start-min-angle 0 --features /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/sq-t1-eu-q0.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/sq-t1-eu-q0-r2.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 4801 × 4801 (10 m) |
| domain vertices | 4 (1 ring, 0 holes) |
| start vertices | 51395 |
| start triangles | 102627 |
| output vertices | 3329882 |
| output triangles | 6655874 |
| constraint edges | 87501 |
| vertices without data dropped | 0 |
| sq-t1-eu-q0.vtk | 393.5 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.04 % | 1.01 % | 0.0306° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 23 | 3479 | 12 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 38 | 3278487 | 0 | 7689701 | 0 | 0 | 0 | 31986 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.005 | 0.0 % |
| decode | 0.135 | 0.3 % |
| features read | 18.415 | 38.4 % |
| features clip | 8.516 | 17.8 % |
| start mesh: build | 0.002 | 0.0 % |
| start mesh: node | 0.148 | 0.3 % |
| start mesh: triangulate | 0.044 | 0.1 % |
| start mesh: constraint edges | 0.110 | 0.2 % |
| refine | 5.174 | 10.8 % |
| refine: legalise start | 0.003 | 0.0 % |
| refine: start quality | 0.000 | 0.0 % |
| refine: scan (parallel) | 1.022 | 2.1 % |
| refine: split + flip (serial) | 3.748 | 7.8 % |
| refine: setup + output | 0.401 | 0.8 % |
| trim | 0.301 | 0.6 % |
| write: encode | 14.973 | 31.3 % |
| write: disk | 0.063 | 0.1 % |
| other | 0.018 | 0.0 % |
| **total** | **47.905** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 1.500 s, not included above.
