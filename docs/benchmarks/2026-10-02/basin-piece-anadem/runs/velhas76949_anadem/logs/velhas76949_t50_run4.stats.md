# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-anadem/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 50 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece-anadem/velhas76949_t50.bin.vtk --stats docs/benchmarks/2026-10-02/basin-piece-anadem/runs/velhas76949_anadem/logs/velhas76949_t50_run4.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 41076 |
| output triangles | 74841 |
| constraint edges | 7309 |
| vertices without data dropped | 0 |
| velhas76949_t50.bin.vtk | 2.9 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 37.22° | 0.01 % | 0.40 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 7 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 50 m | 49.99881540770707 m | 26 | 23272 | 0 | 63762 | 0 | 10495 | 341 | 0 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 2.4 % |
| decode | 0.205 | 60.3 % |
| start mesh: build | 0.000 | 0.1 % |
| start mesh: node | 0.009 | 2.5 % |
| start mesh: triangulate | 0.004 | 1.2 % |
| start mesh: constraint edges | 0.009 | 2.7 % |
| refine | 0.094 | 27.7 % |
| refine: legalise start | 0.000 | 0.1 % |
| refine: start quality | 0.012 | 3.5 % |
| refine: scan (parallel) | 0.065 | 19.2 % |
| refine: split + flip (serial) | 0.010 | 3.0 % |
| refine: setup + output | 0.006 | 1.8 % |
| trim | 0.003 | 0.9 % |
| write: encode | 0.005 | 1.3 % |
| write: disk | 0.002 | 0.7 % |
| other | 0.001 | 0.3 % |
| **total** | **0.340** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.013 s, not included above.
