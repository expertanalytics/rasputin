# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-anadem/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 1 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece-anadem/velhas76949_t1.bin.vtk --stats docs/benchmarks/2026-10-02/basin-piece-anadem/runs/velhas76949_anadem/logs/velhas76949_t1_run3.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 4353814 |
| output triangles | 8699303 |
| constraint edges | 8323 |
| vertices without data dropped | 0 |
| velhas76949_t1.bin.vtk | 313.4 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.00 % | 0.04 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 91 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 42 | 4336010 | 0 | 8670211 | 0 | 10495 | 341 | 1014 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 0.1 % |
| decode | 0.213 | 3.5 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.1 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.010 | 0.2 % |
| refine | 5.240 | 85.5 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.2 % |
| refine: scan (parallel) | 0.870 | 14.2 % |
| refine: split + flip (serial) | 3.995 | 65.2 % |
| refine: setup + output | 0.363 | 5.9 % |
| trim | 0.371 | 6.1 % |
| write: encode | 0.104 | 1.7 % |
| write: disk | 0.169 | 2.8 % |
| other | 0.001 | 0.0 % |
| **total** | **6.130** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 1.825 s, not included above.
