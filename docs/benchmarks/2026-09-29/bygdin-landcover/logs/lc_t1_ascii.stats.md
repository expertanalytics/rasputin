# rasputin mesh — statistics

`rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925 --domain docs/benchmarks/2026-09-29/bygdin-landcover/bygdin_reduced_t20.geojson --tolerance 1 --features ../rasputin_data/corine2018_dtm10_utm33.gpkg --features-layer corine2018 --features-map corine --ascii --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/lc/lc_t1_ascii.vtk --stats docs/benchmarks/2026-09-29/bygdin-landcover/logs/lc_t1_ascii.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2195 × 3314 (10 m) |
| domain vertices | 740 (1 ring, 0 holes) |
| start vertices | 8242 |
| start triangles | 15614 |
| output vertices | 573758 |
| output triangles | 1145434 |
| constraint edges | 14074 |
| vertices without data dropped | 0 |
| lc_t1_ascii.vtk | 69.0 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.10 % | 1.36 % | 0.00408° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 29 | 1377 | 37 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 24 | 549527 | 0 | 1152221 | 0 | 15989 | 2442 | 5726 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.002 | 0.0 % |
| decode | 0.310 | 4.7 % |
| features read | 0.118 | 1.8 % |
| features clip | 1.418 | 21.6 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.017 | 0.3 % |
| start mesh: triangulate | 0.005 | 0.1 % |
| start mesh: constraint edges | 0.016 | 0.2 % |
| refine | 0.536 | 8.2 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.015 | 0.2 % |
| refine: scan (parallel) | 0.094 | 1.4 % |
| refine: split + flip (serial) | 0.377 | 5.7 % |
| refine: setup + output | 0.049 | 0.7 % |
| trim | 0.042 | 0.6 % |
| land cover | 1.256 | 19.1 % |
| write: encode | 2.810 | 42.8 % |
| write: disk | 0.036 | 0.5 % |
| other | 0.005 | 0.1 % |
| **total** | **6.571** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.201 s, not included above.
