# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-anadem/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 2 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece-anadem/velhas76949_t2.bin.vtk --stats docs/benchmarks/2026-10-02/basin-piece-anadem/runs/velhas76949_anadem/logs/velhas76949_t2_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 2296371 |
| output triangles | 4584784 |
| constraint edges | 7956 |
| vertices without data dropped | 0 |
| velhas76949_t2.bin.vtk | 165.3 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 37.87° | 0.00 % | 0.06 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 46 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 2 m | 2 m | 41 | 2278567 | 0 | 5024846 | 0 | 10495 | 341 | 647 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.010 | 0.3 % |
| decode | 0.308 | 9.3 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.3 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.009 | 0.3 % |
| refine | 2.685 | 80.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.4 % |
| refine: scan (parallel) | 0.579 | 17.4 % |
| refine: split + flip (serial) | 1.895 | 56.8 % |
| refine: setup + output | 0.199 | 6.0 % |
| trim | 0.186 | 5.6 % |
| write: encode | 0.074 | 2.2 % |
| write: disk | 0.045 | 1.4 % |
| other | 0.002 | 0.1 % |
| **total** | **3.333** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.931 s, not included above.
