# rasputin mesh — statistics

`rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/catchment_t0.geojson --tolerance 10 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/mesh_fine_t10_run1.vtk --stats docs/benchmarks/2026-09-29/bygdin/logs/mesh_fine_t10_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 5513 (1 ring, 0 holes) |
| start vertices | 5513 |
| start triangles | 5511 |
| output vertices | 38914 |
| output triangles | 72292 |
| constraint edges | 5534 |
| vertices without data dropped | 0 |
| mesh_fine_t10_run1.vtk | 2.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.01 % | 1.10 % | 0.655° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 13 | 3 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.99996099528471 m | 22 | 22612 | 0 | 54390 | 0 | 10789 | 1841 | 21 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 1.8 % |
| decode | 0.312 | 80.0 % |
| start mesh: build | 0.000 | 0.1 % |
| start mesh: node | 0.006 | 1.6 % |
| start mesh: triangulate | 0.003 | 0.8 % |
| start mesh: constraint edges | 0.007 | 1.8 % |
| refine | 0.043 | 11.0 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.011 | 2.8 % |
| refine: scan (parallel) | 0.018 | 4.5 % |
| refine: split + flip (serial) | 0.009 | 2.3 % |
| refine: setup + output | 0.005 | 1.2 % |
| trim | 0.003 | 0.7 % |
| write: encode | 0.004 | 1.1 % |
| write: disk | 0.002 | 0.5 % |
| other | 0.003 | 0.8 % |
| **total** | **0.390** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.013 s, not included above.
