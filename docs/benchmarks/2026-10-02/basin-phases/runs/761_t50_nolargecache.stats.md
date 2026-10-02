# rasputin mesh — statistics

`rasputin mesh --dem anadem-v1 --cache /Users/skavhaug/projects/rasputin_data/cache --domain /Users/skavhaug/projects/rasputin_data/sao_francisco_piece/bho2017_50k_level3/761_epsg4674.geojson --out-crs 'PROJCRS["undefined",BASEGEOGCRS["SIRGAS 2000",DATUM["Sistema de Referencia Geocentrico para las AmericaS 2000",ELLIPSOID["GRS 1980",6378137,298.257222101,LENGTHUNIT["metre",1]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433]],ID["EPSG",4674]],CONVERSION["unknown",METHOD["Transverse Mercator",ID["EPSG",9807]],PARAMETER["Latitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8801]],PARAMETER["Longitude of natural origin",-42,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["Scale factor at natural origin",0.997548,SCALEUNIT["unity",1],ID["EPSG",8805]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]' --tolerance 50 --binary --out /Users/skavhaug/projects/rasputin_scratch/basin-phases/761_t50_nlc.vtk --stats runs/761_t50_nolargecache.stats.md`

## Sizes

| item | count |
|---|---|
| DEM nodes | 19377 × 27786 (30 m) |
| domain vertices | 36983 (1 ring, 0 holes) |
| start vertices | 36983 |
| start triangles | 36981 |
| output vertices | 186334 |
| output triangles | 335681 |
| constraint edges | 36985 |
| vertices without data dropped | 0 |
| 761_t50_nlc.vtk | 13.3 MB |

## Quality (plan view, x/y)

| metric | median | < 1° | < 10° | worst |
|---|---|---|---|---|
| minimum angle | 36.42° | 0.01 % | 1.02 % | 0.00584° |

| metric | median | p99 | max | ≥ 12 | ≥ 20 |
|---|---|---|---|---|---|
| vertex degree (triangles) | 6 | 9 | 15 | 102 | 0 |

## Refinement

| tolerance | achieved max error | rounds | inserted | carved | flips | uncovered | quality inserted | quality skipped | feet |
|---|---|---|---|---|---|---|---|---|---|
| 50 m | 49.99978802443161 m | 33 | 94294 | 0 | 254205 | 0 | 53335 | 1220 | 2 |

## Timings

Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it. Threads: 10 (hardware concurrency).

| phase | seconds | share |
|---|---|---|
| domain read | 0.036 | 0.1 % |
| decode | 19.975 | 52.4 % |
| start mesh: build | 0.001 | 0.0 % |
| start mesh: node | 0.054 | 0.1 % |
| start mesh: triangulate | 0.023 | 0.1 % |
| start mesh: constraint edges | 0.048 | 0.1 % |
| refine | 1.119 | 2.9 % |
| refine: legalise start | 0.001 | 0.0 % |
| refine: start quality | 0.065 | 0.2 % |
| refine: scan (parallel) | 0.965 | 2.5 % |
| refine: split + flip (serial) | 0.052 | 0.1 % |
| refine: setup + output | 0.037 | 0.1 % |
| check points: project | 9.343 | 24.5 % |
| check points: store | 5.906 | 15.5 % |
| final check: scan (parallel) | 1.156 | 3.0 % |
| final check: split + flip (serial) | 0.004 | 0.0 % |
| trim | 0.016 | 0.0 % |
| write: encode | 0.015 | 0.0 % |
| write: disk | 0.002 | 0.0 % |
| other | 0.406 | 1.1 % |
| **total** | **38.103** | **100 %** |

Sub-rows sum to their parent and are not added to the total.

Statistics computed in 0.068 s, not included above.
