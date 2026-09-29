# rasputin mesh — statistics

`rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/catchment_t20_run1.geojson --tolerance 10 --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/mesh_red_t10_run1.vtk --stats docs/benchmarks/2026-09-29/bygdin/logs/mesh_red_t10_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 740 |
| start triangles | 738 |
| output vertices | 26576 |
| output triangles | 52198 |
| constraint edges | 952 |
| vertices without data dropped | 0 |
| mesh_red_t10_run1.vtk | 1.9 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 34.43° | 0.01 % | 1.23 % | 0.356° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 12 | 3 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999777071593144 m | 28 | 24634 | 0 | 60095 | 0 | 1202 | 145 | 212 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.004 | 1.0 % |
| decode | 0.306 | 85.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.001 | 0.2 % |
| start mesh: triangulate | 0.000 | 0.1 % |
| start mesh: constraint edges | 0.001 | 0.2 % |
| refine | 0.034 | 9.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.001 | 0.3 % |
| refine: scan (parallel) | 0.020 | 5.5 % |
| refine: split + flip (serial) | 0.011 | 3.1 % |
| refine: setup + output | 0.003 | 0.7 % |
| trim | 0.002 | 0.5 % |
| write: encode | 0.004 | 1.2 % |
| write: disk | 0.001 | 0.3 % |
| other | 0.003 | 0.9 % |
| **total** | **0.356** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.010 s, not included above.
