# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15c-2/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs '+proj=tmerc +lat_0=0 +lon_0=-44.1 +k=0.999972 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs' --tolerance 5 --binary --out /Users/skavhaug/.claude/jobs/85c14e7c/tmp/perf15c2-scr/tmerc_t5.bin.vtk --stats docs/benchmarks/2026-10-02/15c-2-acceptance/runs/geo/logs/tmerc_t5_run2.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7336 × 4237 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 933820 |
| output triangles | 1860101 |
| constraint edges | 7537 |
| vertices without data dropped | 0 |
| tmerc_t5.bin.vtk | 67.2 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.00 % | 0.28 % | 0.0487° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 13 | 38 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 5 m | 5 m | 35 | 772720 | 0 | 1882855 | 0 | 10534 | 334 | 228 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 0.1 % |
| decode | 1.237 | 19.3 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.1 % |
| start mesh: triangulate | 0.004 | 0.1 % |
| start mesh: constraint edges | 0.013 | 0.2 % |
| refine | 1.019 | 15.9 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.013 | 0.2 % |
| refine: scan (parallel) | 0.354 | 5.5 % |
| refine: split + flip (serial) | 0.581 | 9.0 % |
| refine: setup + output | 0.071 | 1.1 % |
| check points: project | 0.571 | 8.9 % |
| check points: store | 0.289 | 4.5 % |
| final check: scan (parallel) | 0.000 | 0.0 % |
| final check: split + flip (serial) | 0.000 | 0.0 % |
| trim | 0.068 | 1.1 % |
| write: encode | 0.027 | 0.4 % |
| write: disk | 0.012 | 0.2 % |
| other | 3.168 | 49.3 % |
| **total** | **6.424** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.380 s, not included above.
