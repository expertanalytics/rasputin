# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/square48_25833.geojson --tolerance 1 --features /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/sq-t1-eu-q25-fix.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/sq-t1-eu-q25-fix-r3.md`

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
| sq-t1-eu-q25-fix.vtk | 399.2 MB |

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
| decode | 0.133 | 0.4 % |
| features read | 1.011 | 3.3 % |
| features clip | 8.490 | 27.7 % |
| start mesh: build | 0.002 | 0.0 % |
| start mesh: node | 0.148 | 0.5 % |
| start mesh: triangulate | 0.044 | 0.1 % |
| start mesh: constraint edges | 0.117 | 0.4 % |
| refine | 5.309 | 17.4 % |
| refine: legalise start | 0.003 | 0.0 % |
| refine: start quality | 0.204 | 0.7 % |
| refine: scan (parallel) | 0.886 | 2.9 % |
| refine: split + flip (serial) | 3.828 | 12.5 % |
| refine: setup + output | 0.389 | 1.3 % |
| trim | 0.286 | 0.9 % |
| write: encode | 14.988 | 49.0 % |
| write: disk | 0.055 | 0.2 % |
| other | 0.012 | 0.0 % |
| **total** | **30.597** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 1.474 s, not included above.
