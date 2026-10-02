# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-anadem/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 10 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece-anadem/velhas76949_t10.bin.vtk --stats docs/benchmarks/2026-10-02/basin-piece-anadem/runs/velhas76949_anadem/logs/velhas76949_t10_run5.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 328614 |
| output triangles | 649845 |
| constraint edges | 7381 |
| vertices without data dropped | 0 |
| velhas76949_t10.bin.vtk | 23.6 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.15 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 23 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999976905616563 m | 35 | 310810 | 0 | 789337 | 0 | 10495 | 341 | 72 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 1.2 % |
| decode | 0.200 | 28.8 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 1.2 % |
| start mesh: triangulate | 0.004 | 0.6 % |
| start mesh: constraint edges | 0.009 | 1.3 % |
| refine | 0.416 | 59.7 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 1.7 % |
| refine: scan (parallel) | 0.200 | 28.8 % |
| refine: split + flip (serial) | 0.171 | 24.6 % |
| refine: setup + output | 0.032 | 4.5 % |
| trim | 0.023 | 3.3 % |
| write: encode | 0.013 | 1.9 % |
| write: disk | 0.014 | 1.9 % |
| other | 0.001 | 0.2 % |
| **total** | **0.697** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.123 s, not included above.
