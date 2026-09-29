# Increment 22 acceptance: Bygdin (@perf, 2026-09-29): summary

**Verdict: ACCEPTED** against the design's acceptance, which asks for the
area to be within 2 % of NVE. The reduced catchment is **304.9096 km²
against NVE's 305.54 km², -0.206 %**.

There is **no performance baseline**: this is the first `catchment` run and
the first mesh of this domain. Its speed figures are the baseline for the
next run.

- **Run conditions**: AC power at the start and end of every block. The
  commit is `6757925`, with a Release `_core`. Method, tables and raw logs:
  `bygdin/README.md`.

**`rasputin catchment`** (seed 8.5425 61.3512, CORINE lake, DTM10):

- **Windows**: three (4.35 M, 8.15 M and 16.94 M nodes). The design
  expected two.
- **Result**: 3,049,095 nodes. The fine outline has 17,812 vertices and the
  reduced one 740, at 20 m. Both are 304.909550 km²; they differ by
  -3e-07 m².
- **Hausdorff distance, fine to reduced**: 19.92 m at 20 m.
- **Times, median of 3**: flood 4.68 s, trace 0.13 s, reduce 0.04 s, wall
  6.53 s.
- **Memory**: max RSS 989 MB (peak footprint 920 MB).

`--outline-tolerance` ladder (the area is equal at every step; all four
outlines are valid):

| tolerance | vertices | Hausdorff |
|---|---|---|
| 0 | 5,513 | 0 m |
| 10 | 1,382 | 9.99 m |
| 20 | 740 | 19.92 m |
| 50 | 295 | 49.85 m |

**Against NVE's polygon** (delfelt 1187, fetched without GDAL and saved as
`bygdin/nve_delfelt_1187.geojson`), comparing DEM nodes:

- 99.12 % of NVE's nodes are in ours;
- 99.33 % of ours are in NVE's;
- 2.05 km² are ours only, and 2.67 km² are NVE's only.

This is the design's prototype figure, reproduced by the built code.

**`rasputin mesh --domain`** (10 threads, median of 3):

| domain | tol | triangles | refine | total | max RSS | worst angle | max degree |
|---|---|---|---|---|---|---|---|
| reduced | 10 m | 52,198 | 0.034 s | 0.36 s | 646 MB | 0.356° | 12 |
| fine | 10 m | 72,292 | 0.043 s | 0.39 s | 647 MB | 0.655° | 13 |
| reduced | 1 m | 1,129,026 | 0.497 s | 0.93 s | 925 MB | 0.046° | 18 |
| fine | 1 m | 1,137,703 | 0.483 s | 0.92 s | 942 MB | 0.083° | 16 |
| reduced + CORINE | 10 m | 83,171 | 0.046 s | 1.97 s | 776 MB | 0.0041° | 16 |

- **Checks**: every mesh is within its tolerance, has 0 uncovered nodes and
  0 constrained-Delaunay violations.
- **The reduction** saves 27.8 % of the triangles at 10 m and 0.8 % at 1 m.
- **Land cover inside the catchment** (CORINE 2018):
  - sparse vegetation (333) 42 %;
  - bare rock (332) 21 %;
  - heath (322) 17 %;
  - water (512) 16 %;
  - glacier (335) 2.4 %.

**For the main session**:

- The design says NVE's polygon is "not committed", but this brief asked for
  it to be saved. It is committed, and the conflict is flagged in the
  README.
- No defects were found.
