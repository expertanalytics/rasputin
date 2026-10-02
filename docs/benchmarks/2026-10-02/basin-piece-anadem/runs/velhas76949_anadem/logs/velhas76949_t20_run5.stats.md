# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-anadem/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 20 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece-anadem/velhas76949_t20.bin.vtk --stats docs/benchmarks/2026-10-02/basin-piece-anadem/runs/velhas76949_anadem/logs/velhas76949_t20_run5.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 136077 |
| output triangles | 264833 |
| constraint edges | 7319 |
| vertices without data dropped | 0 |
| velhas76949_t20.bin.vtk | 9.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 37.05° | 0.00 % | 0.20 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 14 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 20 m | 19.999762066340054 m | 30 | 118273 | 0 | 311670 | 0 | 10495 | 341 | 10 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 1.8 % |
| decode | 0.202 | 43.0 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 1.8 % |
| start mesh: triangulate | 0.004 | 0.8 % |
| start mesh: constraint edges | 0.009 | 1.9 % |
| refine | 0.212 | 45.1 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 2.6 % |
| refine: scan (parallel) | 0.130 | 27.8 % |
| refine: split + flip (serial) | 0.054 | 11.6 % |
| refine: setup + output | 0.015 | 3.1 % |
| trim | 0.010 | 2.2 % |
| write: encode | 0.007 | 1.5 % |
| write: disk | 0.007 | 1.4 % |
| other | 0.001 | 0.2 % |
| **total** | **0.469** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.049 s, not included above.
