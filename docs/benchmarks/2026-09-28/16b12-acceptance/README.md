# Increment 16b-1/2 acceptance (@perf, 2026-09-28/29)

**Verdict: NOT ACCEPTED as designed, on one measure.** With Ola's European
GeoPackage (EPSG:3035), reading the candidate features takes **17.5-19.9 s**.
The design admits about 1 s for this read (R5). Everything else passes: the
standard benchmark is unchanged, the quarter circle reproduces M5, both
GeoPackage routes give bit-identical meshes on Ola's catchment, and the CORINE
baseline matches 16b-0's within 0.1 % of the triangles. Short summary:
`../16b12-acceptance.md`.

## What was measured, and on what

- **Branch** `increment16b12-features` at `5f3a522`, the head after `54afc27`.
  **Base**: `origin/master` `14f5fe3` (#107, 16b-0 merged), in a detached
  `git worktree` in the session scratchpad.
  `git diff --stat 14f5fe3 5f3a522 -- include src_cpp bindings src_python`
  lists Python files only (`cli.py`, `chains.py`, `feature_input.py`,
  `io/geopackage.py`, `io/gml.py`, …). No C++ changed.
- **Builds**: bench.py's Release builds, `-O3 -DNDEBUG`, AppleClang 21. The
  `_core` sha256 prefix is `190587e7ba0e` for base and branch in all 8 bench
  runs. That is the same `.so` as 16b-0's. The CLI runs used
  `build-pyext`, rebuilt (`cmake --build build-pyext -j --target _core`,
  exit 0) and copied into `.venv`. It has the same sha256. `tin_engine` was
  loaded from the branch's `src_python`, checked with `__file__`. No C++ was
  patched, so there was no scratch build to sanitize.
- **Machine**: Apple M1 Max (8 P + 2 E), 32 GiB, macOS 27.0, Python 3.14.7.
  `caffeinate -ims` was held for every block.
- **Power: battery, discharging, for every run.** The brief expected AC, but
  `pmset` said battery from the first check. The coordinator confirmed that
  Ola could not plug in. The level fell from 83 % (23:36) to 67 % (00:13).
  Each bench `run.json` carries `pmset` before and after. Every CLI case line
  in `logs/*.log` names the power source. `logs/power.log` has the
  percentages around the CLI blocks. **All comparisons here are battery
  against battery.** There is no AC evidence for 16b-1/2 or for 16b-0.
- The runs crossed midnight. The directory keeps the date the acceptance
  started.

## 1. The quarter circle with the committed extract (M5)

`rasputin mesh --dem tests/fixtures/dem_archive/7908_3_10m_z33.tif --domain
docs/benchmarks/2026-09-26/quarter.geojson --tolerance {1,10} [--features
tests/fixtures/corine/clc2018_7908_3.gpkg --features-map corine] --out X.vtk
--stats Y.md` (`scripts/runs.sh qc`). Each case ran 3 times. Triangles and
quality were identical across the 3 runs, and the times are medians.

| tol | features | start vertices | triangles | M5 triangles | < 1° | M5 | worst | M5 | max degree | achieved max error | `node` | M5 `node` | total |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 1 m | none | 536 | 428 217 | 428 217 | 0.00 % | 0.00 % | 0.396° | 0.40° | 18 | 0.99998 m | — | — | 0.97 s |
| 1 m | CORINE | 5 627 | 441 855 | 441 863 | 0.08 % | 0.08 % | 0.0395° | 0.040° | 28 | 0.9999987 m | 0.013 s | 0.33 s | 1.34 s |
| 10 m | none | 536 | 31 581 | 31 581 | 0.01 % | 0.01 % | 0.287° | 0.29° | 16 | 9.9987 m | — | — | 0.15 s |
| 10 m | CORINE | 5 627 | 54 840 | 54 840 | 0.04 % | 0.04 % | 0.0395° | 0.040° | 16 | 9.9980 m | 0.014 s | 0.33 s | 0.31 s |

- **M5 is reproduced.** The noded start vertices equal M5's 5 627 input
  vertices. Three of the four triangle counts are equal, and 1 m with CORINE
  is 8 fewer (−0.002 %). `node` is 25x faster than M5's figure, which was
  measured before 16b-0.
- **The tolerance holds** in every run, with 0 valid DEM nodes uncovered and
  0 feet refused.
- **stderr**: `60 features kept, 9 dropped outside, 19 clipped, 0 empty
  skipped`, then `10802 input vertices, 5627 noded vertices`. There is no
  `table scanned` line, because the extract has an R-tree.
- **VTK reads the file.** `scripts/check_vtk.py` reads each file with
  `vtkPolyDataReader`, VTK 9.7.0 (`logs/qc.log`). At 1 m with CORINE it reads
  221 480 points, 441 855 triangles and 9 345 lines. The fields present are
  `feature_bits`, `feature_names`, `feature_vocabulary`, `crs`,
  `elevation_source` (with the clause `start domain boundary and features,
  vertex z bilinear`), `domain`, `domain_crs`, `domain_transform`,
  `features` (`clc2018_7908_3.gpkg:U2018_CLC2018_V2020_20u1, map corine, 60
  features, 87 chains, 10266 vertices`), `features_crs` (EPSG:3035),
  `features_transform` and `features_notice` (the Copernicus attribution).
  The cell arrays are `feature_mask`, `land_cover` and `water`.
- **The constraint edges carry the bits**: 8 242 of 9 345 carry `land_cover`
  and 3 905 carry `water`. Every `water` edge also carries `land_cover`. The
  1 103 edges with mask 0 are domain-boundary edges that no feature shares.
  At 10 m, 5 698 of 6 035 edges carry `land_cover` and 768 carry `water`.
  Without features, the file has no `land_cover` or `water` arrays and no
  `features*` fields.

## 2. Ola's case: a 4326 catchment over the DTM10 archive

**The domain** is `catchment/catchment_4326.geojson`, written by
`scripts/make_catchment.py`. It is not a real catchment. There is no
committed 4326 catchment, so a smooth, irregular 48-vertex ring was made in
lon/lat. Its area is **287.3 km²** in EPSG:25833, between 11.2 and 11.6° E
around 59.94° N. It is centred on the corner shared by DTM10 tiles `6602_1`,
`6602_2`, `6603_3` and `6603_4` (300 000, 6 650 000 in 25833). So `--dem
DTM10_UTM33_20260925` (the directory) mosaics 4 tiles and crosses 4 seams. It
also lies inside the 254 tiles that `corine2018_dtm10_utm33.gpkg` covers.

The command was `rasputin mesh --dem ../rasputin_data/DTM10_UTM33_20260925
--domain catchment_4326.geojson --tolerance {1,10}`, with one of:

- **Europe (reproject-first route)**: `--features
  .../U2018_CLC2018_V2020_20u1.gpkg --features-layer U2018_CLC2018_V2020_20u1
  --features-map corine`. That file is EPSG:3035, 8.9 GB, with 2 375 406
  R-tree rows.
- **Norway (same-CRS route)**: `--features corine2018_dtm10_utm33.gpkg
  --features-layer corine2018 --features-map corine`. That file is
  EPSG:25833, with 124 054 rows.

Each case ran 3 times (`scripts/runs.sh ola`; logs `ola-r1.log` and
`ola-r23.log`). Times are medians.

| tol | features | triangles | start vertices | < 1° | worst | max degree | features read | features clip | `node` | refine (scan / split) | encode | total | peak RSS |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 1 m | none | 1 216 839 | 48 | 0.02 % | 0.0116° | 25 | — | — | — | 0.604 (0.144 / 0.410) | 2.02 | 3.36 s | 0.93 GB |
| 1 m | Europe, 3035 | 1 235 226 | 5 642 | 0.06 % | 0.0116° | 25 | **17.90** | **8.50** | 0.013 | 0.624 (0.114 / 0.437) | 2.80 | **30.34 s** | 0.80 GB |
| 1 m | Norway, 25833 | 1 235 226 | 5 642 | 0.06 % | 0.0116° | 25 | 0.894 | 1.269 | 0.013 | 0.611 (0.113 / 0.427) | 2.79 | 6.28 s | 0.89 GB |
| 10 m | none | 30 366 | 48 | 0.14 % | 0.0188° | 12 | — | — | — | 0.044 | 0.05 | 0.58 s | 0.47 GB |
| 10 m | Europe, 3035 | 63 020 | 5 642 | 0.23 % | 0.0078° | 14 | **19.08** | **8.44** | 0.013 | 0.044 | 0.14 | **28.23 s** | 0.47 GB |
| 10 m | Norway, 25833 | 63 020 | 5 642 | 0.23 % | 0.0078° | 14 | 0.875 | 1.291 | 0.013 | 0.042 | 0.14 | 3.09 s | 0.51 GB |

- **Both routes give the same mesh, bit for bit**: points, triangles, lines
  and `feature_mask` hash equal at both tolerances (`logs/route_identity.log`).
  The files differ only in `features`, `features_crs` and
  `features_transform`: `none` for Norway, and `Inverse of Europe Equal Area
  2001 + UTM zone 33N` for Europe.
- **stderr**, the same in all 3 runs:
  - Europe: `87 features kept, 258 dropped outside, 28 clipped, 0 empty
    skipped`.
  - Norway: `87 features kept, 175 dropped outside, 28 clipped, 0 empty
    skipped`.
  - Both: `11316 input vertices, 5642 noded vertices`.
  - No scan line: both files have R-trees. The difference in "dropped
    outside" depends on the index, as R10 says it may.
- The file records `features: …, map corine, 87 features, 155 chains, 11268
  vertices`, `domain_crs EPSG:4326` and `domain_transform Inverse of ETRS89
  to WGS 84 (1) + UTM zone 33N`, and the `dem_tiles` and `dem_seams` fields.
  The tolerance holds, with 0 nodes uncovered and 0 feet refused. At 1 m,
  11 027 constraint edges: 10 275 carry `land_cover` and 979 carry `water`.

### The candidate lookup: measured at 17.5-19.9 s, not ~1 s (the defect)

The design (R5, "Long edges") admits that the widened R-tree query scans the
R-tree's columns: "on Ola's Europe file, for a 50 km box, 1 006 rows in
0.75-1.3 s". Measured on this catchment (`scripts/candidates.py`,
`scripts/query_plan.py`; logs `candidates.log`, `query_plan.log`):

| Europe file, this catchment's query box | rows | time |
|---|---:|---:|
| plain box, ids only | 132 | 0.2 ms |
| widened query, ids only (the R-tree scan the design admits) | 345 | 0.71 s |
| **widened query as built** (`io/geopackage.py` `query_features`: JOIN, ORDER BY pk), blobs fetched | 345 (235 MB) | **17.5-18.7 s** |
| the same rows fetched by primary key (`WHERE pk IN (SELECT id FROM rtree …)`) | 345 (235 MB) | 0.85 s |
| WKB decode of the 345 blobs | | 0.25 s |

- **Cause, measured with `EXPLAIN QUERY PLAN`**: SQLite 3.53.4 plans the
  JOIN as `SCAN t; SCAN r VIRTUAL TABLE INDEX 1`. It walks all 2.4 M rows of
  the feature table and probes the R-tree once per row. The `IN` subquery
  gets `SEARCH t USING INTEGER PRIMARY KEY` and returns the same rows in
  0.85 s, 22x faster. On the Norway file (124 054 rows), the same plan costs
  0.77 s against 0.062 s.
- **The widening also brings in large far-away polygons.** 213 of the 345
  Europe candidates are rows that only the widening admits. They hold
  13.8 M of the 14.7 M candidate vertices, and polygons up to 409 242
  vertices. All of them are then dropped outside: 87 are kept, and the
  plain box already returns 132. On the Norway file, 132 rows only the
  widening admits hold 3.1 M of 3.9 M vertices.
- **`features clip` is 8.4-8.6 s (Europe) and 1.27-1.29 s (Norway).** This
  has not been profiled. It is thought to be the reproject-first transform
  and pre-clip of all 14.7 M (3.9 M) candidate vertices. The design's
  0.09 s per 360 000 points would give about 3.7 s for the transform alone.
- On the Europe file, feature input is 26 s of a 28-30 s run whose mesh
  work is under 4 s. The design admitted about 1 s. This is a defect in the
  branch, reported, not fixed. `src_python/` is outside this remit.

## 3. The 1 m benchmark and thread sweep (`tools/bench.py`, blob `77765b1`)

Four pairs, base then branch, back to back, 23:36-23:55, on battery
(83 % → 75 %). The standard inputs are the quarter circle and the tile of
`7908_3` at 1 m, with no features. Threads 1-20 plus the default (10), and
5 repeats. The commands are in `scripts/pairs.sh` (`PAIRS=4` for the
fourth). Each branch run names its pair's base with `--baseline`.

bench.py's verdicts:

- Pairs 1, 2 and 4: **ACCEPTED**.
- Pair 3: **REGRESSION** on 12 quarter cells (up to +13.4 %, at 1 thread).
- The base runs have no `--baseline`, so they are judged against the newest
  stored comparable ancestor. Base r3 and r4 flagged cells too, although
  they are the same build as base r1.

The fourth pair was run as a tie-break after pair 3.

**Pair 3 is judged noise, not a regression:**

- **The code is identical.** `_core` has the same sha256 in all 8 runs,
  `refine_s` times the C++ `refine` call alone, and both domains' mesh
  sha256 is equal in all 8 runs, so `refine` got the same input.
- The same cells were −1.7 %, −1.1 % and +0.8 % (quarter, 1 thread) in
  pairs 1, 2 and 4.
- The same build moves more than that between its own runs: up to +13.6 %
  for 16b-1/2 (r2 → r3, quarter, 1 thread) and +11.1 % for the base.
- **Pooled over all 42 cells** (median of the 4 base runs against median of
  the 4 branch runs): median **+0.27 %**, range −3.3 % to +3.4 %, and no
  cell above +5 %.

Refine seconds, median of 5 per run (ms):

| domain | threads | base r1 / r2 / r3 / r4 | 16b-1/2 r1 / r2 / r3 / r4 | pair changes (%) | pooled |
|---|---:|---|---|---|---:|
| quarter | default (10) | 167.2 / 164.9 / 163.4 / 166.6 | 167.1 / 166.1 / 173.0 / 170.4 | −0.1, +0.7, +5.8, +2.3 | +1.8 % |
| quarter | 1 | 422.8 / 416.2 / 412.4 / 441.4 | 415.4 / 411.7 / 467.6 / 445.1 | −1.7, −1.1, +13.4, +0.8 | +2.6 % |
| quarter | 2 | 277.9 / 269.4 / 266.4 / 286.6 | 267.0 / 267.3 / 295.5 / 287.1 | −3.9, −0.8, +11.0, +0.2 | +1.3 % |
| quarter | 4 | 201.6 / 199.8 / 198.5 / 208.4 | 199.5 / 199.3 / 208.1 / 206.6 | −1.0, −0.3, +4.9, −0.8 | +1.2 % |
| quarter | 8 | 162.7 / 161.7 / 162.6 / 163.2 | 162.6 / 162.0 / 166.6 / 165.9 | −0.1, +0.2, +2.5, +1.6 | +1.0 % |
| quarter | 10 | 165.0 / 163.7 / 166.0 / 172.7 | 163.9 / 164.1 / 176.7 / 165.5 | −0.6, +0.3, +6.5, −4.2 | −0.4 % |
| quarter | 16 | 168.5 / 165.0 / 163.1 / 170.3 | 162.7 / 164.9 / 171.6 / 169.4 | −3.5, −0.1, +5.2, −0.5 | +0.2 % |
| quarter | 20 | 165.8 / 165.0 / 165.5 / 176.0 | 164.2 / 165.3 / 168.8 / 167.3 | −1.0, +0.2, +2.0, −5.0 | +0.4 % |
| tile | default (10) | 191.9 / 185.2 / 193.3 / 191.8 | 187.3 / 191.7 / 183.1 / 189.0 | −2.4, +3.5, −5.2, −1.5 | −1.9 % |
| tile | 1 | 489.7 / 465.7 / 493.5 / 512.1 | 471.4 / 473.1 / 479.7 / 494.6 | −3.7, +1.6, −2.8, −3.4 | −3.1 % |
| tile | 2 | 318.4 / 301.6 / 329.2 / 333.2 | 300.6 / 306.9 / 319.2 / 321.5 | −5.6, +1.8, −3.0, −3.5 | −3.3 % |
| tile | 4 | 232.3 / 227.1 / 233.5 / 234.4 | 228.8 / 232.1 / 229.5 / 236.8 | −1.5, +2.2, −1.7, +1.0 | −0.9 % |
| tile | 8 | 185.3 / 184.6 / 190.9 / 188.8 | 185.4 / 191.3 / 186.1 / 192.8 | +0.0, +3.6, −2.5, +2.1 | +0.9 % |
| tile | 10 | 187.5 / 184.7 / 187.5 / 194.3 | 184.6 / 186.4 / 186.1 / 189.4 | −1.5, +0.9, −0.8, −2.5 | −0.7 % |
| tile | 16 | 187.9 / 186.3 / 189.6 / 190.7 | 188.0 / 185.0 / 186.9 / 190.4 | +0.1, −0.7, −1.4, −0.2 | −0.7 % |
| tile | 20 | 189.6 / 187.4 / 191.5 / 193.7 | 186.8 / 185.7 / 187.4 / 195.2 | −1.5, −0.9, −2.1, +0.8 | −1.8 % |

**Ceiling** (speed-up over 1 thread: the best, and at 20 threads):

- Base: tile 2.52-2.71x, quarter 2.57-2.70x.
- 16b-1/2: tile 2.56-2.63x, quarter 2.54-2.82x.
- On battery, 16b-0's runs gave 2.5-2.9x (`../16b0-acceptance/README.md`).
- From 8 threads on, the curve is flat in every run. All 8 values are in
  `logs/pairs_table.txt`.

**Quality: identical in all 8 runs.**

| domain | worst angle | max degree | within tolerance | Delaunay violations | mesh sha256 |
|---|---:|---:|---|---:|---|
| quarter | 0.3955° | 18 | yes | 0 of 641 791 | `ccebf96a86c6c5e2…` |
| tile | 0.6296° | 74 | yes | 0 of 692 056 | `11741a81adfa17b3…` |

Both hashes equal 16b-0's.

## 4. The CORINE baseline update: M4's 48 km square

The square is 301 000-349 000 × 6 651 000-6 699 000 in EPSG:25833
(`catchment/square48_25833.geojson`), meshed over tile `6603_4` of
`DTM10_UTM33_20260925`. It now goes **through the real CLI**: `rasputin mesh
--dem …/6603_4_10m_z33.tif --domain square48_25833.geojson --tolerance {1,10}
[--start-min-angle 0] [--features U2018_CLC2018_V2020_20u1.gpkg
--features-layer U2018_CLC2018_V2020_20u1 --features-map corine]`
(`scripts/runs.sh sq`). This is the same Europe file as 16b-0. There were 3
runs per row, triangles identical across them, and the times are medians.
Peak RSS is `/usr/bin/time -l` of the whole `rasputin` process, including
writing the `.vtk` and the `--stats` pass.

| tol | features | start quality | triangles | 16b-0 | start-quality nodes (16b-0) | feet (16b-0) | < 1° | worst | features read + clip | `node` (16b-0) | start quality | refine (scan / split) | total | peak RSS (16b-0) |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|
| 1 m | none | 25° | 6 581 348 | 6 581 348 | 0 | 0 | 0.00 % | 0.217° | — | — | — | 5.27 s (1.53 / 3.39) | 16.9 s | 3.07 GB (1.93) |
| 1 m | none | off | 6 581 348 | 6 581 348 | 0 | 0 | 0.00 % | 0.217° | — | — | — | 5.07 s (1.49 / 3.25) | 16.6 s | 3.03 GB (1.96) |
| 1 m | CORINE | 25° | 6 753 084 | 6 752 517 | 154 573 (154 640) | 30 311 (30 311) | 0.03 % | 0.00081° | 17.7 + 8.6 s | 0.150 s (0.068) | 0.209 s | 5.87 s (1.00 / 4.27) | 48.0 s | 3.06 GB (2.35) |
| 1 m | CORINE | off | 6 655 874 | 6 655 862 | 0 | 31 986 (31 990) | 0.04 % | 0.031° | 18.4 + 8.6 s | 0.148 s (0.067) | — | 5.45 s (1.03 / 3.88) | 47.9 s | 2.88 GB (2.27) |
| 10 m | none | 25° | 133 379 | 133 379 | 0 | 0 | 0.03 % | 0.258° | — | — | — | 0.60 s | 0.98 s | 0.57 GB (0.41) |
| 10 m | none | off | 133 379 | 133 379 | 0 | 0 | 0.03 % | 0.258° | — | — | — | 0.60 s | 0.96 s | 0.57 GB (0.41) |
| 10 m | CORINE | 25° | **480 941** | 480 961 | 154 573 (154 640) | 350 (350) | 0.25 % | 0.0013° | 17.8 + 8.6 s | 0.155 s (0.067) | 0.209 s | 0.40 s | 28.5 s | 0.60 GB (0.69) |
| 10 m | CORINE | off | 212 255 | 212 255 | 0 | 359 (359) | 3.45 % | 0.0147° | 14.3 + 8.7 s | 0.150 s (0.066) | — | 0.26 s | 24.5 s | 0.52 GB (0.62) |

- **The baseline holds through the real CLI.** Without features, the
  triangles are identical to 16b-0's. With CORINE they differ by −20 to
  +567 (at most 0.008 %), and the start-quality nodes by 67 of 154 640. The
  quality-start tripling is unchanged: **480 941 triangles at 10 m, 3.61x
  the featureless 133 379**, and 2.27x the same features with quality off.
  That is still 20c's input.
- **What differs from 16b-0, and why the numbers are not one-to-one:**
  - 16b-0 fed pre-merged, deduplicated lines (layout A, `clc_mesh.py`). The
    CLI feeds each clipped ring as it comes: 103 420 input vertices, noded to
    the same 51 395.
  - So `node` is 0.15 s here, against 0.068 s for layout A. 16b-0 measured
    0.1435 s for layout B.
  - The small triangle differences are thought to come from the different
    chains (layout, clip), not measured.
- **Peak RSS is not comparable to 16b-0.** This run writes the 6.6 M-triangle
  `.vtk` (encode 11-15 s) and computes `--stats`. 16b-0's `drv.py` did
  neither. The featureless 1 m run rose by the same 1.1 GB, so the rise is
  thought to be in the method, not in 16b-1/2.
- **`features read` is 14-18 s here as well**: the same query-plan defect
  as in part 2, on a 48 km box. Of a 48 s run at 1 m, 26 s is feature input.

## Found in passing

- **`tools/bench.py`: a `--label` with a `/` and `--mesh-dir` fails.** The
  quality mesh path is `mesh_dir / f"{label}_{d.name}.vtk"`
  (`tools/bench.py:605`), and its parent directory is never created. The
  child exits with `BadParameter: … is not an existing directory` and exit
  3. The first attempt at 23:33 failed this way and produced no evidence;
  its log was kept in the scratchpad, not here. The workaround was to create
  `<mesh-dir>/16b12-acceptance` first. This is probably the "bench.py
  mesh-dir error" that 16b-0's README mentions. It is reported, not fixed:
  a fix to `tools/` needs a red test from `@tester` first.

## Not measured

- A profile of `features clip` on the Europe file. The share of the
  transform in it is thought to be, not measured.
- Whether the widening's far rows (213 on this catchment) could be excluded
  cheaply. That is a design question.
- Any AC run.
- The 16b-0-style layout-A mesh on this branch. `drv.py` monkeypatches
  `cli._domain_chains`, which no longer exists.

## Files

- `base-14f5fe3[-r2|-r3|-r4]/`, `16b12[-r2|-r3|-r4]/`: bench.py run
  directories (`run.json`, `raw.tsv`, generated `README.md`).
- `catchment/`: `catchment_4326.geojson` (Ola's case) and
  `square48_25833.geojson` (M4's square).
- `stats/`: every CLI run's `--stats` report, `<case>-r<n>.md`.
- `logs/`: `pairs.log` and `pair4.log` (bench); `qc.log`, `ola-r1.log`,
  `ola-r23.log` and `sq.log` (the CLI runs, with stderr, peak RSS and the
  VTK read-back after each case's first run); `candidates.log`,
  `query_plan.log`, `route_identity.log`, `pairs_table.txt`, `aggregate.md`,
  `power.log`.
- `scripts/`: `pairs.sh`, `runs.sh`, `check_vtk.py`, `make_catchment.py`,
  `candidates.py`, `query_plan.py`, `aggregate.py` and `pairs_table.py`. Run
  them from the repository root, except `pairs_table.py`, which runs from
  this directory. They name the session scratchpad for the meshes and the
  base worktree; change `SP` to rerun.
- The meshes are not kept. They stayed in the session scratchpad
  (`p16b12/meshes/`). The largest is the 48 km square at 1 m: 6.7 M
  triangles, 399 MB. To regenerate them, rerun
  `scripts/runs.sh <qc|ola|sq> 1`.
