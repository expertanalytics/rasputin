# rasputin mesh — statistics

`rasputin mesh --dem /Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925 --domain docs/benchmarks/2026-09-28/16b12-acceptance/catchment/catchment_4326.geojson --tolerance 1 --features /Users/skavhaug/projects/rasputin_data/corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1 --features-map corine --out /private/tmp/claude-501/-Users-skavhaug-projects-rasputin/38487caf-f56e-46a3-b05b-867af1fb1619/scratchpad/p16b12/meshes/ola-t1-eu3035.vtk --stats docs/benchmarks/2026-09-28/16b12-acceptance/stats/ola-t1-eu3035-r1.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 2089 × 2007 (10 m) |
| domain vertices | 48 (1 ring, 0 holes) |
| start vertices | 5642 |
| start triangles | 11169 |
| output vertices | 617990 |
| output triangles | 1235226 |
| constraint edges | 11027 |
| vertices without data dropped | 0 |
| ola-t1-eu3035.vtk | 69.8 MB |

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
| minimum angle | 43.60° | 0.06 % | 0.91 % | 0.0116° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 25 | 731 | 9 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 1 m | 1 m | 29 | 596775 | 0 | 1282174 | 0 | 15573 | 5291 | 5314 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.001 | 0.0 % |
| decode | 0.482 | 1.6 % |
| features read | 17.468 | 58.3 % |
| features clip | 8.585 | 28.6 % |
| start mesh: build | 0.000 | 0.0 % |
| start mesh: node | 0.013 | 0.0 % |
| start mesh: triangulate | 0.004 | 0.0 % |
| start mesh: constraint edges | 0.012 | 0.0 % |
| refine | 0.624 | 2.1 % |
| refine: legalise start | 0.000 | 0.0 % |
| refine: start quality | 0.016 | 0.1 % |
| refine: scan (parallel) | 0.114 | 0.4 % |
| refine: split + flip (serial) | 0.437 | 1.5 % |
| refine: setup + output | 0.056 | 0.2 % |
| trim | 0.047 | 0.2 % |
| write: encode | 2.713 | 9.1 % |
| write: disk | 0.010 | 0.0 % |
| other | 0.007 | 0.0 % |
| **total** | **29.967** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.226 s, not included above.
