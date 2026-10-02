# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/15c-2/build-bench/pkg --threads 0 -- mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --out-crs '+proj=tmerc +lat_0=0 +lon_0=-44.1 +k=0.999972 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs' --tolerance 20 --binary --out /Users/skavhaug/.claude/jobs/85c14e7c/tmp/perf15c2-scr/tmerc_t20.bin.vtk --stats docs/benchmarks/2026-10-02/15c-2-acceptance/runs/geo/logs/tmerc_t20_run1.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 7336 × 4237 (30 m) |
| domain vertices | 7309 (1 ring, 0 holes) |
| start vertices | 7309 |
| start triangles | 7307 |
| output vertices | 141698 |
| output triangles | 276076 |
| constraint edges | 7318 |
| vertices without data dropped | 0 |
| tmerc_t20.bin.vtk | 10.2 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.97° | 0.00 % | 0.20 % | 0.0828° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 13 | 13 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 20 m | 19.999885979599867 m | 33 | 118143 | 0 | 311797 | 0 | 10534 | 334 | 9 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.007 | 0.2 % |
| decode | 1.198 | 42.4 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.008 | 0.3 % |
| start mesh: triangulate | 0.004 | 0.2 % |
| start mesh: constraint edges | 0.012 | 0.4 % |
| refine | 0.255 | 9.0 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.4 % |
| refine: scan (parallel) | 0.172 | 6.1 % |
| refine: split + flip (serial) | 0.056 | 2.0 % |
| refine: setup + output | 0.014 | 0.5 % |
| check points: project | 0.589 | 20.8 % |
| check points: store | 0.290 | 10.3 % |
| final check: scan (parallel) | 0.000 | 0.0 % |
| final check: split + flip (serial) | 0.000 | 0.0 % |
| trim | 0.010 | 0.4 % |
| write: encode | 0.007 | 0.2 % |
| write: disk | 0.001 | 0.1 % |
| other | 0.442 | 15.7 % |
| **total** | **2.824** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.049 s, not included above.
