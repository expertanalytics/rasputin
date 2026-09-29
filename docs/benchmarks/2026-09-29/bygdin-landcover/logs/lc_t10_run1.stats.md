# rasputin mesh — statistics

`rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 --domain docs/benchmarks/2026-09-29/bygdin-landcover/bygdin_reduced_t20.geojson --tolerance 10 --features ../rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --binary --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/lc/lc_t10_run1.vtk --stats docs/benchmarks/2026-09-29/bygdin-landcover/logs/lc_t10_run1.stats.md`

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
| lc_t10_run1.vtk | 3.7 MB |

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
| domain read | 0.011 | 0.4 % |
| decode | 0.621 | 24.5 % |
| features read | 0.227 | 9.0 % |
| features clip | 1.498 | 59.1 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.017 | 0.7 % |
| start mesh: triangulate | 0.005 | 0.2 % |
| start mesh: constraint edges | 0.016 | 0.6 % |
| refine | 0.048 | 1.9 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.016 | 0.6 % |
| refine: scan (parallel) | 0.014 | 0.5 % |
| refine: split + flip (serial) | 0.010 | 0.4 % |
| refine: setup + output | 0.009 | 0.3 % |
| trim | 0.004 | 0.2 % |
| land cover | 0.075 | 3.0 % |
| write: encode | 0.002 | 0.1 % |
| write: disk | 0.001 | 0.0 % |
| other | 0.008 | 0.3 % |
| **total** | **2.534** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.016 s, not included above.
