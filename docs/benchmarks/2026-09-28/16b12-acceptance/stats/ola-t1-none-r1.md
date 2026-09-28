# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/catchment_4326.geojson --tolerance 1 --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/ola-t1-none.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/ola-t1-none-r1.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2089 × 2007 (10 m) |
| domain vertices | 48 (1 ring, 0 holes) |
| start vertices | 48 |
| start triangles | 46 |
| output vertices | 608791 |
| output triangles | 1216839 |
| constraint edges | 741 |
| vertices without data dropped | 0 |
| ola-t1-none.vtk | 63.6 MB |

## DEM seams

| tile | tile | nodes | max | median |
|---|---|---|---|---|
| 6602_1_10m_z33.tif | 6603_3_10m_z33.tif | 2472 | 0.25087 | 0.0138931 |
| 6602_1_10m_z33.tif | 6603_4_10m_z33.tif | 18977 | 0.376282 | 0.00230408 |
| 6602_2_10m_z33.tif | 6603_3_10m_z33.tif | 59754 | 1.74188 | 0.0246887 |
| 6603_3_10m_z33.tif | 6603_4_10m_z33.tif | 51512 | 2.10982 | 0.0206909 |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 45.00° | 0.02 % | 0.27 % | 0.0116° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 25 | 75 | 1 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 40 | 608706 | 0 | 1338108 | 0 | 37 | 0 | 693 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.008 | 0.2 % |
| decode | 0.795 | 22.6 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.000 | 0.0 % |
| start mesh: triangulate | 0.000 | 0.0 % |
| start mesh: constraint edges | 0.000 | 0.0 % |
| refine | 0.618 | 17.6 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.000 | 0.0 % |
| refine: scan (parallel) | 0.146 | 4.1 % |
| refine: split + flip (serial) | 0.423 | 12.1 % |
| refine: setup + output | 0.049 | 1.4 % |
| trim | 0.054 | 1.6 % |
| write: encode | 2.015 | 57.4 % |
| write: disk | 0.018 | 0.5 % |
| other | 0.003 | 0.1 % |
| **total** | **3.512** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.225 s, not included above.
