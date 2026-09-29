# rasputin mesh — statistics

`rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/catchment_t0.geojson --tolerance 1 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/mesh_fine_t1_run2.vtk --stats docs/benchmarks/2026-09-29/bygdin/logs/mesh_fine_t1_run2.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 5513 (1 ring, 0 holes) |
| start vertices | 5513 |
| start triangles | 5511 |
| output vertices | 571651 |
| output triangles | 1137703 |
| constraint edges | 5597 |
| vertices without data dropped | 0 |
| mesh_fine_t1_run2.vtk | 41.1 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.01 % | 0.27 % | 0.0834° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 16 | 80 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 33 | 555349 | 0 | 1193738 | 0 | 10789 | 1841 | 84 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.006 | 0.7 % |
| decode | 0.312 | 33.5 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.006 | 0.7 % |
| start mesh: triangulate | 0.003 | 0.3 % |
| start mesh: constraint edges | 0.007 | 0.7 % |
| refine | 0.483 | 51.9 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.011 | 1.2 % |
| refine: scan (parallel) | 0.100 | 10.7 % |
| refine: split + flip (serial) | 0.326 | 35.1 % |
| refine: setup + output | 0.046 | 4.9 % |
| trim | 0.044 | 4.7 % |
| write: encode | 0.040 | 4.3 % |
| write: disk | 0.026 | 2.8 % |
| other | 0.003 | 0.4 % |
| **total** | **0.931** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.205 s, not included above.
