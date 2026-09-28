# rasputin mesh — statistics

`rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain docs/benchmarks/2026-09-26/quarter.geojson --tolerance 10 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/qc-t10-none.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/qc-t10-none-r3.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 3001 × 3001 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 536 |
| start triangles | 534 |
| output vertices | 16088 |
| output triangles | 31581 |
| constraint edges | 593 |
| vertices without data dropped | 0 |
| qc-t10-none.vtk | 1.5 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 33.69° | 0.01 % | 2.54 % | 0.287° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 10 | 16 | 26 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.998682022094727 m | 27 | 14908 | 0 | 35522 | 0 | 644 | 0 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 1.2 % |
| decode | 0.069 | 44.2 % |
| start mesh: build | 0.000 | 0.1 % |
| start mesh: node | 0.001 | 0.4 % |
| start mesh: triangulate | 0.000 | 0.3 % |
| start mesh: constraint edges | 0.001 | 0.4 % |
| refine | 0.024 | 15.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.6 % |
| refine: scan (parallel) | 0.016 | 10.3 % |
| refine: split + flip (serial) | 0.006 | 3.8 % |
| refine: setup + output | 0.001 | 0.8 % |
| trim | 0.001 | 0.8 % |
| write: encode | 0.053 | 34.0 % |
| write: disk | 0.002 | 1.2 % |
| other | 0.003 | 2.0 % |
| **total** | **0.155** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.006 s, not included above.
