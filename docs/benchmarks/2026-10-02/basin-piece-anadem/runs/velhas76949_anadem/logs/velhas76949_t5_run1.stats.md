# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-anadem/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 5 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-piece-anadem/velhas76949_t5.bin.vtk --stats docs/benchmarks/2026-10-02/basin-piece-anadem/runs/velhas76949_anadem/logs/velhas76949_t5_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 792105 |
| output triangles | 1576681 |
| constraint edges | 7527 |
| vertices without data dropped | 0 |
| velhas76949_t5.bin.vtk | 57.0 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.14 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 23 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 5 m | 4.999982561383945 m | 37 | 774301 | 0 | 1885927 | 0 | 10495 | 341 | 218 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.009 | 0.7 % |
| decode | 0.283 | 20.7 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.6 % |
| start mesh: triangulate | 0.004 | 0.3 % |
| start mesh: constraint edges | 0.009 | 0.7 % |
| refine | 0.915 | 66.9 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.013 | 0.9 % |
| refine: scan (parallel) | 0.317 | 23.2 % |
| refine: split + flip (serial) | 0.514 | 37.6 % |
| refine: setup + output | 0.071 | 5.2 % |
| trim | 0.066 | 4.8 % |
| write: encode | 0.036 | 2.7 % |
| write: disk | 0.034 | 2.5 % |
| other | 0.002 | 0.1 % |
| **total** | **1.367** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.311 s, not included above.
