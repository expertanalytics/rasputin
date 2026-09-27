# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin/tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/paraview/inputs/quarter.geojson --tolerance 1 --ascii --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/quarter_i20b.vtk --stats /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/quarter_i20b.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 5051 × 5051 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 536 |
| start triangles | 534 |
| output vertices | 214644 |
| output triangles | 428217 |
| constraint edges | 1069 |
| vertices without data dropped | 0 |
| quarter_i20b.vtk | 21.8 MB |

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
| decode | 0.055 | 5.4 % |
| domain read | 0.003 | 0.3 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.003 | 0.3 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.1 % |
| refine | 0.241 | 23.7 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.1 % |
| refine: scan (parallel) | 0.066 | 6.5 % |
| refine: split + flip (serial) | 0.142 | 13.9 % |
| refine: setup + output | 0.032 | 3.2 % |
| trim | 0.016 | 1.5 % |
| write: encode | 0.694 | 68.1 % |
| write: disk | 0.005 | 0.5 % |
| other | 0.001 | 0.1 % |
| **total** | **1.018** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.072 s, not included above.
