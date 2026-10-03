# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15c-2/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs EPSG:31983 --tolerance 5 --binary --out /Users/skavhaug/.claude/jobs/85c14e7c/tmp/perf15c2-scr/EPSG31983_t5.bin.vtk --stats docs/benchmarks/2026-10-02/15c-2-acceptance/runs/geo/logs/EPSG31983_t5_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7346 × 4204 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 934748 |
| output triangles | 1861967 |
| constraint edges | 7527 |
| vertices without data dropped | 0 |
| EPSG31983_t5.bin.vtk | 67.3 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.29 % | 0.0303° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 22 | 0 |

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
| domain read | 0.007 | 0.1 % |
| decode | 0.958 | 15.9 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.1 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.012 | 0.2 % |
| refine | 0.966 | 16.1 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.2 % |
| refine: scan (parallel) | 0.335 | 5.6 % |
| refine: split + flip (serial) | 0.545 | 9.1 % |
| refine: setup + output | 0.074 | 1.2 % |
| check points: project | 0.560 | 9.3 % |
| check points: store | 0.376 | 6.3 % |
| final check: scan (parallel) | 0.000 | 0.0 % |
| final check: split + flip (serial) | 0.000 | 0.0 % |
| trim | 0.068 | 1.1 % |
| write: encode | 0.026 | 0.4 % |
| write: disk | 0.010 | 0.2 % |
| other | 3.017 | 50.2 % |
| **total** | **6.012** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.377 s, not included above.
