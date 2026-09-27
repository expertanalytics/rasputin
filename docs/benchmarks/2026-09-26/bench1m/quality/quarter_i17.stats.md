# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin/tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/paraview/inputs/quarter.geojson --tolerance 1 --ascii --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/quarter_i17.vtk --stats /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/runs/quarter_i17.stats.md`

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
| quarter_i17.vtk | 21.7 MB |

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
| decode | 0.057 | 1.4 % |
| domain read | 0.003 | 0.1 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.003 | 0.1 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.0 % |
| refine | 3.300 | 81.0 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: scan (parallel) | 3.125 | 76.7 % |
| refine: split + flip (serial) | 0.144 | 3.5 % |
| refine: setup + output | 0.032 | 0.8 % |
| trim | 0.016 | 0.4 % |
| write: encode | 0.689 | 16.9 % |
| write: disk | 0.005 | 0.1 % |
| other | 0.001 | 0.0 % |
| **total** | **4.075** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.072 s, not included above.
