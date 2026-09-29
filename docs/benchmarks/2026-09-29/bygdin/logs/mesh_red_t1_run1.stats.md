# rasputin mesh — statistics

`rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/catchment_t20_run1.geojson --tolerance 1 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/mesh_red_t1_run1.vtk --stats docs/benchmarks/2026-09-29/bygdin/logs/mesh_red_t1_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 740 |
| start triangles | 738 |
| output vertices | 565511 |
| output triangles | 1129026 |
| constraint edges | 1994 |
| vertices without data dropped | 0 |
| mesh_red_t1_run1.vtk | 40.7 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.03 % | 0.36 % | 0.0461° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 18 | 102 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 36 | 563569 | 0 | 1218749 | 0 | 1202 | 145 | 1254 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.003 | 0.4 % |
| decode | 0.310 | 33.3 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.001 | 0.1 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.001 | 0.1 % |
| refine | 0.497 | 53.3 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.1 % |
| refine: scan (parallel) | 0.105 | 11.3 % |
| refine: split + flip (serial) | 0.347 | 37.2 % |
| refine: setup + output | 0.044 | 4.7 % |
| trim | 0.043 | 4.6 % |
| write: encode | 0.044 | 4.7 % |
| write: disk | 0.030 | 3.2 % |
| other | 0.003 | 0.3 % |
| **total** | **0.933** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.197 s, not included above.
