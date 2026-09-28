# rasputin mesh — statistics

`rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain docs/benchmarks/2026-09-26/quarter.geojson --tolerance 1 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/qc-t1-none.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/qc-t1-none-r3.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 3001 × 3001 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 536 |
| start triangles | 534 |
| output vertices | 214644 |
| output triangles | 428217 |
| constraint edges | 1069 |
| vertices without data dropped | 0 |
| qc-t1-none.vtk | 21.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 0.82 % | 0.396° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 18 | 184 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 0.9999796549479072 m | 41 | 213464 | 0 | 445657 | 0 | 644 | 0 | 6 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.2 % |
| decode | 0.067 | 7.0 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.001 | 0.1 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.1 % |
| refine | 0.166 | 17.2 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.1 % |
| refine: scan (parallel) | 0.048 | 5.0 % |
| refine: split + flip (serial) | 0.099 | 10.3 % |
| refine: setup + output | 0.018 | 1.9 % |
| trim | 0.016 | 1.6 % |
| write: encode | 0.706 | 73.1 % |
| write: disk | 0.004 | 0.4 % |
| other | 0.003 | 0.3 % |
| **total** | **0.966** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.072 s, not included above.
