# rasputin mesh — statistics

`rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain docs/benchmarks/2026-09-26/quarter.geojson --tolerance 1 --features tests/fixtures/corine/clc2018_7908_3.gpkg --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/qc-t1-corine.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/qc-t1-corine-r1.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 3001 × 3001 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 5627 |
| start triangles | 10679 |
| output vertices | 221480 |
| output triangles | 441855 |
| constraint edges | 9345 |
| vertices without data dropped | 0 |
| qc-t1-corine.vtk | 24.4 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.08 % | 2.07 % | 0.0395° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 10 | 28 | 746 | 6 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 0.9999987363815306 m | 22 | 204483 | 0 | 404582 | 0 | 11370 | 1326 | 3134 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.1 % |
| decode | 0.091 | 6.5 % |
| features read | 0.020 | 1.4 % |
| features clip | 0.039 | 2.8 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.013 | 0.9 % |
| start mesh: triangulate | 0.004 | 0.3 % |
| start mesh: constraint edges | 0.011 | 0.8 % |
| refine | 0.197 | 14.1 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.011 | 0.8 % |
| refine: scan (parallel) | 0.040 | 2.8 % |
| refine: split + flip (serial) | 0.125 | 8.9 % |
| refine: setup + output | 0.021 | 1.5 % |
| trim | 0.018 | 1.3 % |
| write: encode | 0.990 | 70.9 % |
| write: disk | 0.004 | 0.3 % |
| other | 0.006 | 0.4 % |
| **total** | **1.395** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.077 s, not included above.
