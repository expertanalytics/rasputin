# rasputin mesh — statistics

`bench.py _child --pkg /Users/skavhaug/projects/rasputin/.claude/worktrees/basin-anadem/build-bench/pkg --threads 0 -- mesh --dem /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/derived/bho2017_5k_76949_anadem_epsg31983_30m.tif --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_5k_76949_outline_epsg4674.geojson --tolerance 1 --ascii --out /Users/skavhaug/projects/rasputin_scratch/basin-piece-anadem/velhas76949_t1.ascii.vtk --stats docs/benchmarks/2026-10-02/basin-piece-anadem/runs/velhas76949_anadem/logs/velhas76949_t1_ascii.stats.md`

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
| velhas76949_t1.ascii.vtk | 473.0 MB |

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
| domain read | 0.008 | 0.0 % |
| decode | 0.209 | 1.0 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.009 | 0.0 % |
| start mesh: triangulate | 0.004 | 0.0 % |
| start mesh: constraint edges | 0.009 | 0.0 % |
| refine | 5.272 | 26.1 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.012 | 0.1 % |
| refine: scan (parallel) | 0.871 | 4.3 % |
| refine: split + flip (serial) | 4.019 | 19.9 % |
| refine: setup + output | 0.369 | 1.8 % |
| trim | 0.384 | 1.9 % |
| write: encode | 14.152 | 70.1 % |
| write: disk | 0.133 | 0.7 % |
| other | 0.001 | 0.0 % |
| **total** | **20.181** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 1.750 s, not included above.
