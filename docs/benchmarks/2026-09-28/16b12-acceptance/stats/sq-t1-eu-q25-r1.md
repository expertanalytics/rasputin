# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/square48_25833.geojson --tolerance 1 --features /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/sq-t1-eu-q25.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/sq-t1-eu-q25-r1.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 4801 × 4801 (10 m) |
| domain vertices | 4 (1 ring, 0 holes) |
| start vertices | 51395 |
| start triangles | 102627 |
| output vertices | 3378505 |
| output triangles | 6753084 |
| constraint edges | 85862 |
| vertices without data dropped | 0 |
| sq-t1-eu-q25.vtk | 399.2 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.03 % | 0.98 % | 0.00081° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 23 | 3657 | 15 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 30 | 3172537 | 0 | 7127927 | 0 | 154573 | 51555 | 30311 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.0 % |
| decode | 0.142 | 0.3 % |
| features read | 20.190 | 39.1 % |
| features clip | 8.845 | 17.1 % |
| start mesh: build | 0.003 | 0.0 % |
| start mesh: node | 0.156 | 0.3 % |
| start mesh: triangulate | 0.048 | 0.1 % |
| start mesh: constraint edges | 0.118 | 0.2 % |
| refine | 6.156 | 11.9 % |
| refine: legalise start | 0.003 | 0.0 % |
| refine: start quality | 0.221 | 0.4 % |
| refine: scan (parallel) | 1.083 | 2.1 % |
| refine: split + flip (serial) | 4.440 | 8.6 % |
| refine: setup + output | 0.410 | 0.8 % |
| trim | 0.324 | 0.6 % |
| write: encode | 15.423 | 29.9 % |
| write: disk | 0.218 | 0.4 % |
| other | 0.021 | 0.0 % |
| **total** | **51.645** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 1.543 s, not included above.
