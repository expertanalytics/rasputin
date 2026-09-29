# rasputin mesh — statistics

`rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 --domain /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/catchment_t20_run1.geojson --tolerance 10 --features ../rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/run/mesh_feat_t10_run2.vtk --stats docs/benchmarks/2026-09-29/bygdin/logs/mesh_feat_t10_run2.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 8242 |
| start triangles | 15614 |
| output vertices | 42110 |
| output triangles | 83171 |
| constraint edges | 8882 |
| vertices without data dropped | 0 |
| mesh_feat_t10_run2.vtk | 3.3 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.87° | 0.07 % | 1.70 % | 0.00408° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 16 | 31 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999732907715497 m | 17 | 17879 | 0 | 39582 | 0 | 15989 | 2442 | 534 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.003 | 0.2 % |
| decode | 0.308 | 15.6 % |
| features read | 0.131 | 6.6 % |
| features clip | 1.432 | 72.6 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.017 | 0.9 % |
| start mesh: triangulate | 0.005 | 0.3 % |
| start mesh: constraint edges | 0.016 | 0.8 % |
| refine | 0.047 | 2.4 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.015 | 0.8 % |
| refine: scan (parallel) | 0.013 | 0.7 % |
| refine: split + flip (serial) | 0.010 | 0.5 % |
| refine: setup + output | 0.009 | 0.4 % |
| trim | 0.003 | 0.2 % |
| write: encode | 0.005 | 0.3 % |
| write: disk | 0.001 | 0.0 % |
| other | 0.004 | 0.2 % |
| **total** | **1.972** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.015 s, not included above.
