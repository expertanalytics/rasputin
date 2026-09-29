# Increment 22 acceptance: Bygdin, end to end (@perf, 2026-09-29)

The summary and the verdict are in `../bygdin.md`. This file has the method,
the tables and where every number comes from.

## What was run

- **Commit**: `6757925` on `increment22-autocatchment` (`logs/commit.txt`).
- **Build**: `build-pyext` (Release), rebuilt with
  `cmake --build build-pyext -j --target _core` (exit 0, nothing out of
  date) and copied into `.venv`; the installed `_core` sha256 is in
  `logs/commit.txt`, equal to the build's.
- **Machine**: Apple M1 Max, 10 cores, 32 GB. `mesh` used 10 threads
  (hardware concurrency, as `--stats` reports).
- **Power**: AC throughout. `pmset -g batt` was recorded at the start and end
  of the catchment block, of the timed mesh block, and after the quality runs
  (`logs/pmset_*.txt`): every one reads "AC Power, 80 %, AC attached".
- **Data**: DEM `../rasputin_data/DTM10_UTM33_20260925` (10 m, EPSG:25833).
  Lakes and features: the CORINE 2018 Norway extract
  `../rasputin_data/corine2018_dtm10_utm33.gpkg`, layer `corine2018`.
- **Seed**: `--seed 8.5425 61.3512` (EPSG:4326), as in the design; it falls
  in CORINE's Bygdin polygon (39.936963 km², 399,371 nodes).
- **Reference**: NVE's reservoir catchment ("delfelt") 1187 "BYGDIN",
  `delfeltAreal_km2` = 305.54. Its polygon was fetched on 2026-09-29 (by
  curl, no GDAL) with exactly this command:

  ```sh
  curl -sS --max-time 60 -o nve_delfelt_1187.geojson 'https://gis3.nve.no/map/rest/services/Mapservices/VassdragsreguleringVannkraft/MapServer/8/query?where=delfeltNr%3D1187&outFields=*&returnGeometry=true&outSR=25833&f=geojson'
  ```

  The file as fetched has sha256
  `fe6d298b08573e50e50e711f1be9f70e9459ea823a0626516a539cd62cc4e76c`: one
  Polygon, 1376 coordinates, 305.5393 km² by shapely, valid. As the design
  says, it is **not committed** (see "Data and licence" below). The
  service can change, so a later fetch may not match the sha256. If it
  does not, the overlap figures below apply to the file with this hash,
  not to the new one.

Scripts, both in this directory:

- `run.sh` runs every timed block and writes the raw logs to `logs/`:
  - `catchment` three times at the default tolerance, then once each at
    `--outline-tolerance` 0, 10, 20 and 50;
  - five mesh configurations, each three times with `--binary`;
  - then each once more with `--ascii`, for the quality check.
  Each command runs under `/usr/bin/time -l`, which gives wall time,
  maximum RSS and peak memory footprint at the end of each `.err`.
- `analyse.py` reads the logs and outputs and prints `logs/analysis.txt`.
  That is every table below. The quality check calls `tools/bench.py`'s
  `quality()` and `read_vtk_ascii()` unchanged: the angles and degree as
  `--stats` computes them, the tolerance check, and the constrained Delaunay
  check.

Meshes and catchment GeoJSONs other than the committed one were written to
the session scratchpad and are not committed. To regenerate them, run
`SCRATCH=<dir> bash docs/benchmarks/2026-09-29/bygdin/run.sh` from the
repository root with `.venv` active, fetch NVE's polygon as above, then
`python .../analyse.py <dir> <nve_delfelt_1187.geojson>`.
The committed `bygdin_reduced_t20.geojson` is `catchment_t20_run1.geojson`,
the default output. All four default-tolerance outputs have the same sha256,
`8083fb5b...`.

## 1. `rasputin catchment`

Windows, as every run printed them (`logs/catchment_*.err`; all seven runs
chose the same windows):

| window | x (m) | y (m) | rows x cols | nodes | flood (s, run 1/2/3) | outcome |
|---|---|---|---|---|---|---|
| 1 | 140240-170470 | 6811200-6825570 | 1438 x 3024 | 4.35 M | 0.62 / 0.61 / 0.61 | grown on north, west, south, east |
| 2 | 136250-172700 | 6807210-6829560 | 2236 x 3646 | 8.15 M | 1.22 / 1.22 / 1.22 | grown on north, west |
| 3 | 128580-176700 | 6802360-6837550 | 3520 x 4813 | 16.94 M | 2.84 / 2.84 / 2.89 | contained |

The design expected one growth step, to about 4100 x 3000 nodes (12 M). The
loop took two growth steps and ended at 16.9 M nodes.

The result, identical in every run:

| item | value |
|---|---|
| catchment nodes | 3,049,095 (304.909500 km² of node area) |
| fine outline | 17,812 vertices, 304.909550 km², 1 hole filled (0.000050 km²), 0 rings dropped |
| reduced outline (20 m) | 740 vertices, 304.909550 km², difference -2.98e-07 m² (-9.8e-16 relative) |
| against NVE's 305.54 km² | **-0.206 %** (fine and reduced alike; -0.630 km²) |
| Hausdorff, fine to reduced | 19.92 m at the 20 m tolerance |

Times, the median of three runs at the default tolerance, from each run's
stderr and `/usr/bin/time -l`:

| phase | median (s) | runs |
|---|---|---|
| flood, three windows summed | 4.68 | 4.68, 4.67, 4.72 |
| trace | 0.13 | 0.13, 0.14, 0.11 |
| reduce | 0.04 | 0.04, 0.04, 0.04 |
| the rest of the wall time (start-up, DEM read, seed, write) | 1.68 | by difference |
| **wall** | **6.53** | 6.71, 6.48, 6.53 |
| maximum RSS | 989 MB | 1069, 989, 982 MB |
| peak memory footprint | 920 MB | 923, 920, 861 MB |

The "rest" row is the wall time less the three reported phases. It is not
broken down further: the command does not report a read time, and nothing was
profiled.

**`--outline-tolerance` ladder.** One run each. The vertex count, area and
validity come from the output GeoJSON. The Hausdorff distance was measured
two ways between the fine outline's boundary and the reduced one's:

- shapely `hausdorff_distance(..., densify=0.05)`, as the tests measure it;
- every 1 m sample of one boundary against the other's segments, both ways.

| tolerance (m) | vertices | area (km²) | vs fine (m²) | vs NVE | Hausdorff, densify 0.05 | 1 m: fine to reduced, reduced to fine | valid |
|---|---|---|---|---|---|---|---|
| 0 | 5,513 | 304.909550 | 0 | -0.206 % | 0.00 m | 0.00, 0.00 m | yes |
| 10 | 1,382 | 304.909550 | 0 | -0.206 % | 9.99 m | 9.99, 9.96 m | yes |
| 20 (default) | 740 | 304.909550 | -2.98e-07 | -0.206 % | 19.92 m | 19.92, 19.62 m | yes |
| 50 | 295 | 304.909550 | -1.43e-06 | -0.206 % | 49.85 m | 49.85, 48.06 m | yes |

The traced outline has 17,812 vertices. The Hausdorff distance stayed under
the tolerance at every step of the ladder. The design lists that as measured,
not guaranteed.

## 2. Node overlap with NVE's polygon

This compares the DEM nodes (10 m lattice, x and y multiples of 10) strictly
inside our fine outline with those inside NVE's polygon, using shapely
`contains_xy` over the union's bounds.

| item | nodes |
|---|---|
| ours | 3,049,096 |
| NVE's | 3,055,360 |
| both | 3,028,611 |
| **NVE's nodes in ours** | **99.12 %** |
| **ours in NVE's** | **99.33 %** |
| ours only | 20,485 (2.05 km²) |
| NVE's only | 26,749 (2.67 km²) |

Our count is one more than the command's 3,049,095. The fine outline has
one hole filled, of 50 m² (half a cell, the midpoint diamond around a single
out-node), so one node that is outside the catchment falls inside the filled
ring.

The largest pieces of disagreement, as area and centroid in EPSG:25833, are
in `logs/analysis.txt`:

- ours and not NVE's: 0.255 km² at (163743, 6813036);
- NVE's and not ours: 0.577 km² at (149046, 6826265).

Between the two boundaries, the 1 m-sampled distance is at most 485 m from
ours to NVE's and 782 m from NVE's to ours. Why the divides differ was not
investigated. The prototype's figures in the design (99.1 % and 99.3 %,
-0.21 %) are reproduced by the built code.

## 3. `rasputin mesh --domain`

- **Timed runs**: three per configuration, `--binary`, `--stats`. The
  medians are shown.
- **Refine and total** are from `--stats`. Total is the command body, without
  interpreter start-up.
- **Wall and maximum RSS** are from `/usr/bin/time -l`.
- **Quality** is from the one `--ascii` run per configuration, which has the
  same triangle count as the timed runs.
- "reduced" is the default 740-vertex outline; "fine" is
  `--outline-tolerance 0` (5,513 vertices).

| domain | tol | triangles | vertices | start vertices | refine (s) | total (s) | wall (s) | max RSS (MB) | worst angle | max degree | max error (m) | Delaunay checked / ambiguous / violations |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| reduced | 10 m | 52,198 | 26,576 | 740 | 0.034 | 0.356 | 0.62 | 646 | 0.356° | 12 | 9.9998 | 77,821 / 1,122 / 0 |
| fine | 10 m | 72,292 | 38,914 | 5,513 | 0.043 | 0.388 | 0.66 | 647 | 0.655° | 13 | 10.0000 | 105,671 / 1,546 / 0 |
| reduced | 1 m | 1,129,026 | 565,511 | 740 | 0.497 | 0.927 | 1.38 | 925 | 0.0461° | 18 | 1.0000 | 1,692,542 / 172,872 / 0 |
| fine | 1 m | 1,137,703 | 571,651 | 5,513 | 0.483 | 0.917 | 1.38 | 942 | 0.0834° | 16 | 1.0000 | 1,703,756 / 172,995 / 0 |
| reduced + CORINE | 10 m | 83,171 | 42,110 | 8,242 | 0.046 | 1.972 | 2.24 | 776 | 0.00408° | 16 | 9.9997 | 116,398 / 1,095 / 0 |

Every run is within its tolerance, with 0 valid DEM nodes uncovered and 0
Delaunay violations. Per-run figures are in `logs/analysis.txt`. The
refine-time spread between runs is at most 0.008 s.

**What the reduction saves.**

- At 10 m it saves 27.8 % of the triangles: 52,198 against 72,292.
- At 1 m it saves 0.8 %: 1,129,026 against 1,137,703. At that tolerance the
  interior dominates the triangle count, and the outline does not.
- At both tolerances the reduced outline gives the smaller worst angle:
  0.356° against 0.655° at 10 m, and 0.046° against 0.083° at 1 m. The
  maximum degree is 12 against 13 at 10 m, and 18 against 16 at 1 m. The
  cause was not measured.

**Land cover inside the catchment** (`--features <Norway extract>
--features-layer corine2018 --features-map corine`, 10 m). stderr reports:

- 79 features kept, 342 dropped outside, 35 clipped;
- 15,846 input vertices, 8,242 after noding.

In the mesh:

- 8,882 constraint edges;
- 7,835 edge cells with `land_cover` = 1;
- 2,005 edge cells with `water` = 1.

The worst angle, 0.00408°, comes with the CORINE constraints; without them
the same domain gives 0.356°. Of the 1.97 s total, refine takes 0.046 s. The
rest is the feature read, clip and noding, which were not broken down here
(see the `.stats.md` phase table).

The class areas below come from clipping CORINE to the reduced polygon with
shapely, which is separate from the mesh:

| CORINE code | class | km² | share |
|---|---|---|---|
| 333 | sparsely vegetated areas | 129.093 | 42.34 % |
| 332 | bare rocks | 65.012 | 21.32 % |
| 322 | moors and heathland | 51.849 | 17.00 % |
| 512 | water bodies | 50.191 | 16.46 % |
| 335 | glaciers and perpetual snow | 7.452 | 2.44 % |
| 412 | peat bogs | 1.021 | 0.33 % |
| 142 | sport and leisure facilities | 0.292 | 0.10 % |

The classes cover 304.910 km² of the polygon's 304.910 km².

## No performance baseline

This is the first run of `catchment`, and the first mesh of this domain.
There is no earlier figure on AC or on battery to compare speed against, so
the speed figures above are the baseline for the next increment. The
increment does not change refine. The 1 m benchmark and the thread sweep
were not run; this brief did not ask for them.

## Data and licence

- **NVE's catchment polygon** is © Norges vassdrags- og energidirektorat
  (NVE). It is NVE's reservoir sub-catchment ("delfelt") 1187, "BYGDIN",
  from the map service VassdragsreguleringVannkraft, layer 8
  (https://gis3.nve.no/map/rest/services/Mapservices/VassdragsreguleringVannkraft/MapServer/8).
  It is not a REGINE unit.
  - It is used under the Norwegian Licence for Open Government Data
    (NLOD), which is compatible with CC BY 4.0. NVE asks users to credit
    the rights holder and link to its services:
    https://konto.nve.no/Information.
  - The polygon was committed in `ad3df66` and removed in `1985720`, so it
    stays in the history.
  - To reproduce the overlap figures, fetch it with the curl command under
    "What was run".
  - The main session checked these terms on 2026-09-29.
- **The DEM**, DTM10 (`DTM10_UTM33_20260925`), is © Kartverket.
- **CORINE**, the lake polygon and the `--features` land cover, is credited
  as increment 16b's NOTICE words it (`CORINE_NOTICE` in
  `src_python/tin_engine/feature_input.py`): "Contains modified CORINE Land
  Cover 2018 data (version 2020_20u1), (c) European Union, Copernicus Land
  Monitoring Service 2018, European Environment Agency (EEA): clipped and
  re-encoded, produced with funding by the European Union, not endorsed by
  the EU".

## Files

- `run.sh`, `analyse.py`: the scripts above.
- `logs/`:
  - `*.err`: stderr with `/usr/bin/time -l`;
  - `*.out`: stdout;
  - `*.stats.md`: `mesh --stats`;
  - `exits.txt`: every command's exit status, all 0;
  - `pmset_*.txt`;
  - `commit.txt`;
  - `analysis.txt`.
- `bygdin_reduced_t20.geojson`: the reduced catchment (740 vertices,
  EPSG:25833), which `mesh --domain` reads.
