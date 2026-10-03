# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15c-2/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs '+proj=tmerc +lat_0=0 +lon_0=-44.1 +k=0.999972 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs' --tolerance 10 --binary --out /Users/skavhaug/.claude/jobs/85c14e7c/tmp/perf15c2-scr/tmerc_t10.bin.vtk --stats docs/benchmarks/2026-10-02/15c-2-acceptance/runs/geo/logs/tmerc_t10_run2.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7336 × 4237 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 358258 |
| output triangles | 709134 |
| constraint edges | 7380 |
| vertices without data dropped | 0 |
| tmerc_t10.bin.vtk | 25.8 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.19 % | 0.0417° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 13 | 17 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999976021902853 m | 33 | 310907 | 0 | 791238 | 0 | 10534 | 334 | 71 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 0.2 % |
| decode | 1.119 | 31.7 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.008 | 0.2 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.012 | 0.3 % |
| refine | 0.420 | 11.9 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.3 % |
| refine: scan (parallel) | 0.207 | 5.9 % |
| refine: split + flip (serial) | 0.170 | 4.8 % |
| refine: setup + output | 0.031 | 0.9 % |
| check points: project | 0.551 | 15.6 % |
| check points: store | 0.283 | 8.0 % |
| final check: scan (parallel) | 0.000 | 0.0 % |
| final check: split + flip (serial) | 0.000 | 0.0 % |
| trim | 0.025 | 0.7 % |
| write: encode | 0.012 | 0.4 % |
| write: disk | 0.004 | 0.1 % |
| other | 1.080 | 30.6 % |
| **total** | **3.524** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.134 s, not included above.
