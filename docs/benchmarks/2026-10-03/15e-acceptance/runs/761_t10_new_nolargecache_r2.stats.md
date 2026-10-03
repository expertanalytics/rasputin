# rasputin mesh — statistics

`rasputin mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_50k_level3/761_epsg4674.geojson --out-crs 'PROJCRS["undefined",BASEGEOGCRS["SIRGAS 2000",DATUM["Sistema de Referencia Geocentrico para las AmericaS 2000",ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4674]],CONVERSION["unknown",METHOD["Transverse Mercator",ID["EPSG",9807]],PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8801]],PARAMETER["Longitude of natural origin",-42,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["Scale factor at natural origin",0.997548,SCALEUNIT["unity",1],ID["EPSG",8805]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]' --tolerance 10 --binary --out /Users/skavhaug/projects/rasputin_scratch/15e-acceptance/761_t10_new_nolargecache_r2.vtk --stats runs/761_t10_new_nolargecache_r2.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 19377 × 27786 (30 m) |
| domain vertices | 36983 (1 ring, 0 holes) |
| start vertices | 36983 |
| start triangles | 36981 |
| output vertices | 1488537 |
| output triangles | 2940024 |
| constraint edges | 37048 |
| vertices without data dropped | 0 |
| 761_t10_new_nolargecache_r2.vtk | 107.0 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.16° | 0.00 % | 0.41 % | 0.00584° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 17 | 301 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 9.999999395691475 m | 43 | 1314475 | 0 | 3421162 | 0 | 53335 | 1220 | 65 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.035 | 0.1 % |
| decode | 20.777 | 42.1 % |
| start mesh: build | 0.001 | 0.0 % |
| start mesh: node | 0.054 | 0.1 % |
| start mesh: triangulate | 0.023 | 0.0 % |
| start mesh: constraint edges | 0.045 | 0.1 % |
| refine | 4.031 | 8.2 % |
| refine: legalise start | 0.001 | 0.0 % |
| refine: start quality | 0.065 | 0.1 % |
| refine: scan (parallel) | 2.772 | 5.6 % |
| refine: split + flip (serial) | 1.024 | 2.1 % |
| refine: setup + output | 0.170 | 0.3 % |
| check points: project | 9.674 | 19.6 % |
| check points: store | 5.431 | 11.0 % |
| final check: scan (parallel) | 4.340 | 8.8 % |
| final check: split + flip (serial) | 0.129 | 0.3 % |
| trim | 0.137 | 0.3 % |
| write: encode | 0.081 | 0.2 % |
| write: disk | 0.039 | 0.1 % |
| other | 4.511 | 9.1 % |
| **total** | **49.307** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.697 s, not included above.
