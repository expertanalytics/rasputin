# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15c-2/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983 --tolerance 1 --binary --out /Users/skavhaug/.claude/jobs/85c14e7c/tmp/perf15c2-scr/EPSG31983_t1.bin.vtk --stats docs/benchmarks/2026-10-02/15c-2-acceptance/runs/geo/logs/EPSG31983_t1_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7346 × 4204 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 6943339 |
| output triangles | 13878353 |
| constraint edges | 8323 |
| vertices without data dropped | 0 |
| EPSG31983_t1.bin.vtk | 499.9 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 37.02° | 0.00 % | 1.01 % | 0.000353° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 16 | 211 | 0 |

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
| domain read | 0.007 | 0.0 % |
| decode | 0.982 | 3.3 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.0 % |
| start mesh: triangulate | 0.004 | 0.0 % |
| start mesh: constraint edges | 0.012 | 0.0 % |
| refine | 5.770 | 19.5 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.013 | 0.0 % |
| refine: scan (parallel) | 0.933 | 3.1 % |
| refine: split + flip (serial) | 4.438 | 15.0 % |
| refine: setup + output | 0.387 | 1.3 % |
| check points: project | 0.534 | 1.8 % |
| check points: store | 0.392 | 1.3 % |
| final check: scan (parallel) | 0.000 | 0.0 % |
| final check: split + flip (serial) | 0.000 | 0.0 % |
| trim | 0.663 | 2.2 % |
| write: encode | 0.210 | 0.7 % |
| write: disk | 0.094 | 0.3 % |
| other | 20.967 | 70.7 % |
| **total** | **29.645** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 3.458 s, not included above.
