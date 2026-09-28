# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/square48_25833.geojson --tolerance 1 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/sq-t1-none-q25.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/sq-t1-none-q25-r1.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 4801 × 4801 (10 m) |
| domain vertices | 4 (1 ring, 0 holes) |
| start vertices | 4 |
| start triangles | 2 |
| output vertices | 3292586 |
| output triangles | 6581348 |
| constraint edges | 3822 |
| vertices without data dropped | 0 |
| sq-t1-none-q25.vtk | 359.9 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.29 % | 0.217° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 17 | 527 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 61 | 3292582 | 0 | 7657458 | 0 | 0 | 0 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.004 | 0.0 % |
| decode | 0.175 | 1.0 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.000 | 0.0 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.000 | 0.0 % |
| refine | 5.267 | 31.1 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.000 | 0.0 % |
| refine: scan (parallel) | 1.527 | 9.0 % |
| refine: split + flip (serial) | 3.390 | 20.0 % |
| refine: setup + output | 0.350 | 2.1 % |
| trim | 0.308 | 1.8 % |
| write: encode | 11.067 | 65.3 % |
| write: disk | 0.111 | 0.7 % |
| other | 0.003 | 0.0 % |
| **total** | **16.936** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 1.452 s, not included above.
