# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin/tests/fixtures/dem_archive/7908_3_10m_z33.tif --tolerance 1 --ascii --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/tile_head.vtk --stats /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/tile_head.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 5051 × 5051 (10 m) |
| start vertices | 16384 |
| start triangles | 32258 |
| output vertices | 236525 |
| output triangles | 464290 |
| constraint edges | 872 |
| vertices without data dropped | 199 |
| tile_head.vtk | 23.2 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 2.32 % | 0.63° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 74 | 324 | 197 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 0.9999759197235107 m | 53 | 219837 | 7887 | 445675 | 0 | 503 | 252 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| decode | 0.055 | 4.8 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.026 | 2.3 % |
| start mesh: triangulate | 0.026 | 2.3 % |
| start mesh: constraint edges | 0.007 | 0.6 % |
| refine | 0.258 | 22.4 % |
| refine: legalise start | 0.003 | 0.2 % |
| refine: start quality | 0.001 | 0.0 % |
| refine: scan (parallel) | 0.073 | 6.3 % |
| refine: split + flip (serial) | 0.142 | 12.3 % |
| refine: setup + output | 0.041 | 3.5 % |
| trim | 0.017 | 1.5 % |
| write: encode | 0.756 | 65.7 % |
| write: disk | 0.005 | 0.4 % |
| other | 0.001 | 0.1 % |
| **total** | **1.151** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.078 s, not included above.
