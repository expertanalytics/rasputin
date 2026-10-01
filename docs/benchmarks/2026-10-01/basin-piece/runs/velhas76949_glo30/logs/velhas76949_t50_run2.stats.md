# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-measure/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 50 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece/velhas76949_t50.bin.vtk --stats docs/benchmarks/2026-10-01/basin-piece/runs/velhas76949_glo30/logs/velhas76949_t50_run2.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 44914 |
| output triangles | 82517 |
| constraint edges | 7309 |
| vertices without data dropped | 0 |
| velhas76949_t50.bin.vtk | 3.2 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 37.06° | 0.01 % | 0.43 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 11 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 50 m | 49.998221261160666 m | 25 | 27110 | 0 | 74625 | 0 | 10495 | 341 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 2.4 % |
| decode | 0.211 | 59.7 % |
| start mesh: build | 0.000 | 0.1 % |
| start mesh: node | 0.009 | 2.4 % |
| start mesh: triangulate | 0.004 | 1.2 % |
| start mesh: constraint edges | 0.010 | 2.7 % |
| refine | 0.102 | 28.7 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.013 | 3.6 % |
| refine: scan (parallel) | 0.070 | 19.9 % |
| refine: split + flip (serial) | 0.012 | 3.3 % |
| refine: setup + output | 0.007 | 1.9 % |
| trim | 0.003 | 0.9 % |
| write: encode | 0.005 | 1.4 % |
| write: disk | 0.001 | 0.2 % |
| other | 0.001 | 0.3 % |
| **total** | **0.354** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.015 s, not included above.
