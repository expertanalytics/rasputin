# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-measure/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 20 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece/velhas76949_t20.bin.vtk --stats docs/benchmarks/2026-10-01/basin-piece/runs/velhas76949_glo30/logs/velhas76949_t20_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 160916 |
| output triangles | 314507 |
| constraint edges | 7323 |
| vertices without data dropped | 0 |
| velhas76949_t20.bin.vtk | 11.6 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.23 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 18 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 20 m | 19.999989197317518 m | 31 | 143112 | 0 | 377249 | 0 | 10495 | 341 | 14 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.009 | 1.6 % |
| decode | 0.253 | 45.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 1.6 % |
| start mesh: triangulate | 0.004 | 0.7 % |
| start mesh: constraint edges | 0.009 | 1.7 % |
| refine | 0.243 | 44.0 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 2.2 % |
| refine: scan (parallel) | 0.142 | 25.8 % |
| refine: split + flip (serial) | 0.070 | 12.8 % |
| refine: setup + output | 0.018 | 3.3 % |
| trim | 0.012 | 2.1 % |
| write: encode | 0.010 | 1.8 % |
| write: disk | 0.002 | 0.3 % |
| other | 0.001 | 0.2 % |
| **total** | **0.551** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.060 s, not included above.
