# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-measure/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_glo30_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 10 --ascii --out /Users/skavhaug/projects/rasputin_scratch/basin-piece/velhas76949_t10.ascii.vtk --stats docs/benchmarks/2026-10-01/basin-piece/runs/velhas76949_glo30/logs/velhas76949_t10_ascii.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7339 × 4200 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 407376 |
| output triangles | 807337 |
| constraint edges | 7413 |
| vertices without data dropped | 0 |
| velhas76949_t10.ascii.vtk | 41.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.18 % | 0.113° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 22 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999978027343786 m | 33 | 389572 | 0 | 987771 | 0 | 10495 | 341 | 104 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.009 | 0.4 % |
| decode | 0.204 | 9.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.4 % |
| start mesh: triangulate | 0.004 | 0.2 % |
| start mesh: constraint edges | 0.009 | 0.4 % |
| refine | 0.496 | 24.0 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.6 % |
| refine: scan (parallel) | 0.225 | 10.9 % |
| refine: split + flip (serial) | 0.222 | 10.7 % |
| refine: setup + output | 0.036 | 1.7 % |
| trim | 0.029 | 1.4 % |
| write: encode | 1.303 | 63.0 % |
| write: disk | 0.005 | 0.3 % |
| other | 0.001 | 0.1 % |
| **total** | **2.070** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.157 s, not included above.
