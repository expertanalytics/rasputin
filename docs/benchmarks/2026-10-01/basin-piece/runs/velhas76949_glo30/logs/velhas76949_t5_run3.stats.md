# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-measure/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 5 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece/velhas76949_t5.bin.vtk --stats docs/benchmarks/2026-10-01/basin-piece/runs/velhas76949_glo30/logs/velhas76949_t5_run3.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 1000480 |
| output triangles | 1993312 |
| constraint edges | 7646 |
| vertices without data dropped | 0 |
| velhas76949_t5.bin.vtk | 72.0 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.16 % | 0.107° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 28 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 5 m | 5 m | 39 | 982676 | 0 | 2373326 | 0 | 10495 | 341 | 337 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 0.5 % |
| decode | 0.199 | 13.2 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.6 % |
| start mesh: triangulate | 0.004 | 0.3 % |
| start mesh: constraint edges | 0.009 | 0.6 % |
| refine | 1.157 | 76.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.8 % |
| refine: scan (parallel) | 0.360 | 23.9 % |
| refine: split + flip (serial) | 0.698 | 46.2 % |
| refine: setup + output | 0.086 | 5.7 % |
| trim | 0.076 | 5.0 % |
| write: encode | 0.033 | 2.2 % |
| write: disk | 0.014 | 0.9 % |
| other | 0.001 | 0.1 % |
| **total** | **1.510** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.399 s, not included above.
