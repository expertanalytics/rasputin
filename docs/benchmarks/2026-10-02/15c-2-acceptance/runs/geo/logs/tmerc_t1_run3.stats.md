# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15c-2/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs '+proj=tmerc +lat_0=0 +lon_0=-44.1 +k=0.999972 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs' --tolerance 1 --binary --out /Users/skavhaug/.claude/jobs/85c14e7c/tmp/perf15c2-scr/tmerc_t1.bin.vtk --stats docs/benchmarks/2026-10-02/15c-2-acceptance/runs/geo/logs/tmerc_t1_run3.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7336 × 4237 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 6942696 |
| output triangles | 13877141 |
| constraint edges | 8249 |
| vertices without data dropped | 0 |
| tmerc_t1.bin.vtk | 499.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 37.00° | 0.00 % | 0.99 % | 0.000646° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 14 | 230 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 41 | 4336804 | 0 | 8672236 | 0 | 10534 | 334 | 940 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 0.0 % |
| decode | 1.153 | 3.8 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.008 | 0.0 % |
| start mesh: triangulate | 0.004 | 0.0 % |
| start mesh: constraint edges | 0.012 | 0.0 % |
| refine | 5.882 | 19.2 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.0 % |
| refine: scan (parallel) | 0.994 | 3.3 % |
| refine: split + flip (serial) | 4.493 | 14.7 % |
| refine: setup + output | 0.382 | 1.3 % |
| check points: project | 0.577 | 1.9 % |
| check points: store | 0.286 | 0.9 % |
| final check: scan (parallel) | 0.000 | 0.0 % |
| final check: split + flip (serial) | 0.000 | 0.0 % |
| trim | 0.771 | 2.5 % |
| write: encode | 0.219 | 0.7 % |
| write: disk | 0.123 | 0.4 % |
| other | 21.526 | 70.4 % |
| **total** | **30.568** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 3.674 s, not included above.
