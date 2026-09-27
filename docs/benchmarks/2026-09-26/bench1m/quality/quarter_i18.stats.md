# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin/tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/paraview/inputs/quarter.geojson --tolerance 1 --ascii --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/quarter_i18.vtk --stats /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/quarter_i18.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 5051 × 5051 (10 m) |
| domain vertices | 536 (1 ring, 0 holes) |
| start vertices | 536 |
| start triangles | 534 |
| output vertices | 214422 |
| output triangles | 427779 |
| constraint edges | 1063 |
| vertices without data dropped | 0 |
| quarter_i18.vtk | 21.7 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.03 % | 0.92 % | 0.0117° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 43 | 216 | 9 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered |
|---|---|---|---|---|---|---|
| 1 m | 0.999997615814209 m | 40 | 213886 | 0 | 452067 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| decode | 0.054 | 4.7 % |
| domain read | 0.003 | 0.3 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.003 | 0.2 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.1 % |
| refine | 0.372 | 32.5 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.203 | 17.8 % |
| refine: split + flip (serial) | 0.138 | 12.0 % |
| refine: setup + output | 0.031 | 2.7 % |
| trim | 0.016 | 1.4 % |
| write: encode | 0.691 | 60.3 % |
| write: disk | 0.004 | 0.4 % |
| other | 0.001 | 0.1 % |
| **total** | **1.145** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.073 s, not included above.
