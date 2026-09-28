# rasputin mesh — statistics

`rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain docs/benchmarks/2026-09-26/quarter.geojson --tolerance 10 --features tests/fixtures/corine/clc2018_7908_3.gpkg --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/qc-t10-corine.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/qc-t10-corine-r2.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 3001 × 3001 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 5627 |
| start triangles | 10679 |
| output vertices | 27738 |
| output triangles | 54840 |
| constraint edges | 6035 |
| vertices without data dropped | 0 |
| qc-t10-corine.vtk | 2.9 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.04 % | 1.48 % | 0.0395° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 16 | 6 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.99795150756836 m | 14 | 10741 | 0 | 22907 | 0 | 11370 | 1326 | 287 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.6 % |
| decode | 0.068 | 21.8 % |
| features read | 0.007 | 2.2 % |
| features clip | 0.036 | 11.7 % |
| start mesh: build | 0.000 | 0.1 % |
| start mesh: node | 0.013 | 4.2 % |
| start mesh: triangulate | 0.004 | 1.3 % |
| start mesh: constraint edges | 0.011 | 3.6 % |
| refine | 0.033 | 10.7 % |
| refine: legalise start | 0.000 | 0.1 % |
| refine: start quality | 0.011 | 3.7 % |
| refine: scan (parallel) | 0.010 | 3.2 % |
| refine: split + flip (serial) | 0.006 | 1.8 % |
| refine: setup + output | 0.006 | 1.8 % |
| trim | 0.002 | 0.7 % |
| write: encode | 0.125 | 40.2 % |
| write: disk | 0.002 | 0.7 % |
| other | 0.007 | 2.2 % |
| **total** | **0.311** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.010 s, not included above.
