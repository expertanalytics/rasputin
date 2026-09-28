# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925/6603_4_10m_z33.tif --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/square48_25833.geojson --tolerance 10 --start-min-angle 0 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/sq-t10-none-q0.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/sq-t10-none-q0-r2.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 4801 × 4801 (10 m) |
| domain vertices | 4 (1 ring, 0 holes) |
| start vertices | 4 |
| start triangles | 2 |
| output vertices | 66932 |
| output triangles | 133379 |
| constraint edges | 483 |
| vertices without data dropped | 0 |
| sq-t10-none-q0.vtk | 6.6 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 34.45° | 0.03 % | 1.78 % | 0.258° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 15 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999497794833673 m | 52 | 66928 | 0 | 186432 | 0 | 0 | 0 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.001 | 0.1 % |
| decode | 0.134 | 13.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.000 | 0.0 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.000 | 0.0 % |
| refine | 0.600 | 62.4 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.563 | 58.5 % |
| refine: split + flip (serial) | 0.030 | 3.1 % |
| refine: setup + output | 0.007 | 0.7 % |
| trim | 0.005 | 0.5 % |
| write: encode | 0.217 | 22.5 % |
| write: disk | 0.001 | 0.1 % |
| other | 0.004 | 0.4 % |
| **total** | **0.962** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.026 s, not included above.
