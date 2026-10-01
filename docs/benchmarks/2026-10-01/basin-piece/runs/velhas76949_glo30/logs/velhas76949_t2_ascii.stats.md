# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-measure/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 2 --ascii --out /Users/skavhaug/projects/rasputin_scratch/basin-piece/velhas76949_t2.ascii.vtk --stats docs/benchmarks/2026-10-01/basin-piece/runs/velhas76949_glo30/logs/velhas76949_t2_ascii.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 2872292 |
| output triangles | 5736406 |
| constraint edges | 8176 |
| vertices without data dropped | 0 |
| velhas76949_t2.ascii.vtk | 309.7 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 0.06 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 13 | 82 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 2 m | 2 m | 43 | 2854488 | 0 | 6145955 | 0 | 10495 | 341 | 867 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 0.1 % |
| decode | 0.208 | 1.6 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.1 % |
| start mesh: triangulate | 0.004 | 0.0 % |
| start mesh: constraint edges | 0.009 | 0.1 % |
| refine | 3.400 | 25.7 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.1 % |
| refine: scan (parallel) | 0.665 | 5.0 % |
| refine: split + flip (serial) | 2.471 | 18.7 % |
| refine: setup + output | 0.251 | 1.9 % |
| trim | 0.227 | 1.7 % |
| write: encode | 9.335 | 70.5 % |
| write: disk | 0.042 | 0.3 % |
| other | 0.001 | 0.0 % |
| **total** | **13.244** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 1.198 s, not included above.
