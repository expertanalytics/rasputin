# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-measure/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 1 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece/velhas76949_t1.bin.vtk --stats docs/benchmarks/2026-10-01/basin-piece/runs/velhas76949_glo30/logs/velhas76949_t1_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 5220159 |
| output triangles | 10431955 |
| constraint edges | 8361 |
| vertices without data dropped | 0 |
| velhas76949_t1.bin.vtk | 375.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 0.04 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 114 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 42 | 5202355 | 0 | 10147510 | 0 | 10495 | 341 | 1052 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.009 | 0.1 % |
| decode | 0.209 | 2.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.1 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.009 | 0.1 % |
| refine | 6.317 | 87.8 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.2 % |
| refine: scan (parallel) | 0.985 | 13.7 % |
| refine: split + flip (serial) | 4.886 | 67.9 % |
| refine: setup + output | 0.434 | 6.0 % |
| trim | 0.433 | 6.0 % |
| write: encode | 0.150 | 2.1 % |
| write: disk | 0.053 | 0.7 % |
| other | 0.001 | 0.0 % |
| **total** | **7.194** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 2.176 s, not included above.
