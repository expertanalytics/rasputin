# Increment 16c acceptance: Bygdin in natural colours (@perf, 2026-09-29)

The summary and the verdict are in `../bygdin-landcover.md`. This file has the
method, the tables and where every number comes from. The design's acceptance
section (`docs/increments/16c-landcover-labels.md`, "Acceptance") names the
directory `16c-bygdin/`; this run was asked to use `bygdin-landcover/`.

## What was run

- **Commit**: `f7d5f14` on `increment16c-landcover-labels` (`logs/commit.txt`).
  The 16c commits (`196147e` to `f7d5f14`) change nothing under `include/`,
  `src/`, `tests/cpp/`, `tools/` or `CMakeLists.txt`
  (`git diff --stat 196147e~1..f7d5f14 -- include src tests/cpp CMakeLists.txt tools`
  is empty), so the bench and thread sweep are not required (design,
  Acceptance).
- **Build**: `cmake --build build-pyext -j --target _core` (exit 0), copied into
  `.venv`; the installed and built `_core` sha256 are equal
  (`logs/commit.txt`). Python is the editable `src_python/` tree.
- **Machine**: Apple M1 Max, 10 cores, 32 GB. `mesh` used 10 threads.
- **Power**: AC throughout. `pmset -g batt` at the start, after the timed
  runs, and at the end (`logs/pmset_*.txt`): each reads "AC Power, 80 %, AC
  attached".
- **Inputs**:
  - domain `bygdin_reduced_t20.geojson`, committed here, copied with
    `git show increment22-autocatchment:docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson`
    (sha256 `8083fb5b...`, the same as 22's);
  - DEM `../rasputin_data/DTM10_UTM33_20260925` (10 m, EPSG:25833);
  - features `../rasputin_data/corine2018_dtm10_utm33.gpkg`, layer
    `corine2018`, `--features-map corine`.

Scripts, all in this directory:

- `run.sh`: three timed `--binary` runs at `--tolerance 10` and at
  `--tolerance 1`, then one `--ascii` run of each for the quality check, then
  `rasputin palette corine --out corine_natural.json`. Every command runs under
  `/usr/bin/time -l`; logs and `--stats` reports go to `logs/`.
- `analyse.py`: every table below (`logs/analysis.txt` is its output). It reads
  the meshes with VTK's `vtkPolyDataReader`. The CORINE reference is read
  straight from the GeoPackage with `sqlite3` and shapely, not through
  `tin_engine.feature_input`, so it does not share code with what it checks.
  The quality check calls `tools/bench.py`'s `quality()` and
  `read_vtk_ascii()` unchanged.
- `render.py`: the two pictures.

Meshes stay out of the repository. To regenerate: from the repository root
with `.venv` active, `SCRATCH=<dir> bash docs/benchmarks/2026-09-29/bygdin-landcover/run.sh`,
then `python .../analyse.py <dir>` and
`python .../render.py <dir>/lc_t10_run1.vtk docs/benchmarks/2026-09-29/bygdin-landcover`.

## 1. The file opens in VTK with `land_cover_code` on every cell

| run | cells | lines | triangles | `land_cover_code` values | nonzero on lines | `FieldData` `land_cover_codes` |
|---|---|---|---|---|---|---|
| 10 m | 92,053 | 8,882 | 83,171 | 92,053 | 0 | present |
| 1 m | 1,159,508 | 14,074 | 1,145,434 | 1,159,508 | 0 | present |

The field text reads `CORINE Land Cover level-3 code, attribute Code_18, map
corine; 0 = in no polygon, and every constraint line`. The active `SCALARS`
stays `feature_mask` (13, ruling 5). The ASCII run's codes equal the binary
run's, cell for cell, at both tolerances.

## 2. Class shares against CORINE clipped to the catchment

79 CORINE polygons meet the domain; clipped to it they cover 304.909550 km²,
the domain's area. The mesh's plan area is 304.909551 km² at both tolerances.
Triangle area (x, y) per code:

| code | class | mesh km² | mesh share | CORINE clipped share | difference (pp) | 22's table |
|---|---|---|---|---|---|---|
| 333 | sparsely vegetated areas | 129.093253 | 42.33821 % | 42.33821 % | -0.000001 | 42.34 % |
| 332 | bare rocks | 65.011569 | 21.32159 % | 21.32159 % | +0.000000 | 21.32 % |
| 322 | moors and heathland | 51.849126 | 17.00476 % | 17.00476 % | +0.000000 | 17.00 % |
| 512 | water bodies | 50.191119 | 16.46099 % | 16.46099 % | +0.000000 | 16.46 % |
| 335 | glaciers and perpetual snow | 7.452202 | 2.44407 % | 2.44407 % | -0.000000 | 2.44 % |
| 412 | peat bogs | 1.020748 | 0.33477 % | 0.33477 % | +0.000000 | 0.33 % |
| 142 | sport and leisure facilities | 0.291534 | 0.09561 % | 0.09561 % | -0.000000 | 0.10 % |

The 10 m and 1 m rows agree to the printed digits; the table is both. The
same seven codes, no triangle with code 0, and the largest difference is
0.000001 percentage points (limit 0.01).

**The check can fail.** Relabelling the 142 region as 333 on the 10 m mesh
moves 333's share to 42.43383 %, 0.096 pp off, and the oracle below then
disagrees on 115 triangles (a one-off probe, not committed).

## 3. I1 (spread) and I2 (oracle)

| run | interior unconstrained edges | I1: edges with differing codes | I2: triangles checked (r > 3 mm) | I2: disagreements |
|---|---|---|---|---|
| 10 m | 116,398 | 0 | 83,171 of 83,171 | 0 |
| 1 m | 1,705,117 | 0 | 1,145,434 of 1,145,434 | 0 |

The oracle is the design's (R1): each triangle's centroid against the CORINE
polygons (read independently, unclipped), smallest area wins; margin
`2 × 1 mm`, so the exact set is `r > 3 mm`, which is every triangle here.
`--stats` stderr, every run: `land cover: 124 regions, 0 outside every
polygon, 0 in more than one, 0 thinner than the snap`.

## 4. Time (recorded, not gated)

Median of three binary runs; phase times from each run's `--stats` report,
wall and RSS from `/usr/bin/time -l`.

| tolerance | triangles | components | union-find rounds | land cover (s) | refine (s) | features clip (s) | total (s) | wall (s) | max RSS |
|---|---|---|---|---|---|---|---|---|---|
| 10 m | 83,171 | 124 | 5 | 0.071 (0.075, 0.071, 0.070) | 0.047 | 1.460 | 2.075 | 2.34 | 719 MB |
| 1 m | 1,145,434 | 124 | 7 | 1.296 (1.303, 1.296, 1.263) | 0.559 | 1.447 | 3.887 | 4.36 | 859 MB |

"Union-find rounds" counts the hook rounds of `landcover.regions`' outer loop:
`analyse.py` runs a copy of that loop with a counter and checks that its
result equals `regions()` (it does, at both tolerances).

At 1 m the `land cover` phase (1.30 s) costs more than `refine` (0.56 s).
The design's R2 names this as the condition for moving `regions` into the
core. Where inside the phase the time goes was not profiled.

## 5. Quality (the ASCII runs, `tools/bench.py`'s `quality()`)

| tolerance | worst angle | max degree | max error | within tolerance | Delaunay checked / ambiguous / violations |
|---|---|---|---|---|---|
| 10 m | 0.00408° | 16 | 9.99973 m | yes | 116,398 / 1,095 / 0 |
| 1 m | 0.00408° | 29 | 1.0 m | yes | 1,705,117 / 162,591 / 0 |

The 10 m figures equal 22's "reduced + CORINE" row (83,171 triangles, 0.00408°,
degree 16). 22 has no 1 m run with CORINE features; its 1 m run without them
had max degree 18. Labelling does not touch the mesh, so these are the
refine's figures on this input, not 16c's.

## 6. The pictures and the ParaView preset

- `bygdin_landcover_oblique.png`: the 10 m mesh from the south-east, z
  exaggerated 2×; `bygdin_landcover_top.png`: map view. `render.py` draws the
  triangles only (the constraint lines, code 0, are left out), coloured through
  an indexed `vtkLookupTable` built from `tin_engine.palettes.CORINE_NATURAL`,
  with a legend of the seven codes present.
- `corine_natural.json`: written by `rasputin palette corine --out`.
  **In ParaView**: colour by `land_cover_code` (cells), then Colour Map
  Editor → Choose Preset → Import `corine_natural.json`, select "rasputin
  CORINE natural", Apply (tick *Interpret Values As Categories* if it is not
  on).
- Not done here (design item 6, manual, for Ola): the ParaView import itself,
  its version, and whether the preset switches to categorical by itself.
