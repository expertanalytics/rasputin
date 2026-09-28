# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/square48_25833.geojson --tolerance 10 --features /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/sq-t10-eu-q25.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/sq-t10-eu-q25-r2.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 4801 × 4801 (10 m) |
| domain vertices | 4 (1 ring, 0 holes) |
| start vertices | 51395 |
| start triangles | 102627 |
| output vertices | 240769 |
| output triangles | 480941 |
| constraint edges | 52572 |
| vertices without data dropped | 0 |
| sq-t10-eu-q25.vtk | 28.1 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.25 % | 4.58 % | 0.0013° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 19 | 219 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999960550447781 m | 17 | 34801 | 0 | 86893 | 0 | 154573 | 51555 | 350 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.001 | 0.0 % |
| decode | 0.131 | 0.5 % |
| features read | 17.609 | 62.1 % |
| features clip | 8.693 | 30.7 % |
| start mesh: build | 0.002 | 0.0 % |
| start mesh: node | 0.161 | 0.6 % |
| start mesh: triangulate | 0.046 | 0.2 % |
| start mesh: constraint edges | 0.114 | 0.4 % |
| refine | 0.419 | 1.5 % |
| refine: legalise start | 0.003 | 0.0 % |
| refine: start quality | 0.209 | 0.7 % |
| refine: scan (parallel) | 0.116 | 0.4 % |
| refine: split + flip (serial) | 0.026 | 0.1 % |
| refine: setup + output | 0.065 | 0.2 % |
| trim | 0.019 | 0.1 % |
| write: encode | 1.132 | 4.0 % |
| write: disk | 0.005 | 0.0 % |
| other | 0.017 | 0.1 % |
| **total** | **28.348** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.090 s, not included above.
