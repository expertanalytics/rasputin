# Increment 16b-1/2 acceptance (@perf, 2026-09-28/29): summary

**Verdict (2026-09-28), superseded by the 2026-09-29 addendum (ACCEPTED): NOT ACCEPTED as designed, on one measure: the candidate read from
Ola's European GeoPackage.** It takes **17.5-19.9 s**, against the
design's admitted ~1 s (R5). Everything else passes. All runs were on
**battery** and are compared against battery. The branch is `5f3a522`; the
base is `origin/master` `14f5fe3`. Method, tables and raw data are in
`16b12-acceptance/README.md`.

- **The defect, with its cause measured.** `io/geopackage.py`
  `query_features` joins the feature table to the R-tree. SQLite plans that
  as a scan of all 2.4 M feature rows, with one R-tree probe per row
  (`EXPLAIN QUERY PLAN`: `SCAN t`). The query takes 17.5-18.7 s. The
  R-tree scan alone, which is what the design timed, takes 0.71 s. The same
  345 rows fetched by primary key through an `IN` subquery take 0.85 s. On
  top of that, `features clip` takes 8.4-8.6 s. It is thought to be the
  reprojection of the 14.7 M candidate vertices, not profiled. 13.8 M of
  those vertices are in rows that only the widening admits, and all of them
  are dropped. With the Norway file (25833), read and clip take 0.9 s and
  1.3 s. Reported, not fixed.
- **The quarter circle (M5), committed extract, `--features-map corine`:**
  M5 is reproduced.
  - Start vertices: 5 627, as in M5.
  - Triangles: 441 855 at 1 m (M5 441 863) and 54 840 at 10 m (equal).
  - Tolerance held, with 0 nodes uncovered.
  - `vtkPolyDataReader` reads the file with R10's fields.
  - 8 242 of 9 345 constraint edges carry `land_cover` and 3 905 carry
    `water`.
  - `node` 0.013 s, against M5's 0.33 s.
- **Ola's case:** a 287 km² ring in EPSG:4326 on the 4-tile corner of
  `6602_1`/`6602_2`/`6603_3`/`6603_4` in DTM10_UTM33_20260925, run with
  both GeoPackages.
  - **The two routes give bit-identical meshes:** 1 235 226 triangles at
    1 m and 63 020 at 10 m.
  - stderr: `87 features kept, 258 (Europe) / 175 (Norway) dropped outside,
    28 clipped, 0 empty skipped`, then `11316 input vertices, 5642 noded
    vertices`. There is no scan line.
  - At 1 m the whole run takes 30.3 s with the Europe file, 6.3 s with the
    Norway file, and 3.4 s without features.
- **The 1 m benchmark and thread sweep:** 4 pairs, no change.
  - The `_core` sha256 and both mesh sha256s are identical in all 8 runs.
  - Pooled over 42 cells, the change is a median of +0.27 % (range −3.3 %
    to +3.4 %), with no cell above +5 %.
  - bench.py flagged 12 quarter cells in pair 3. The same build moves up to
    +13.6 % between its own runs, and the fourth pair, run as a tie-break,
    was ACCEPTED.
  - Ceiling 2.5-2.8x, as the base.
- **CORINE baseline, M4's 48 km square, through the real CLI** (Europe file):
  - It holds against 16b-0: 6 753 084 triangles at 1 m (16b-0: 6 752 517),
    and 480 941 at 10 m (480 961), which is still **3.61x** the
    featureless mesh. That is 20c's input.
  - `node` is 0.15 s, because the CLI feeds the rings unmerged (16b-0's
    layout B: 0.14 s).
  - Refine takes 5.9 s at 1 m. Peak RSS is 3.1 GB, which includes writing
    the file and `--stats`; 16b-0's scripts did neither.
  - Feature input is 26 s of 48 s: the same defect.
- **Also found:** `tools/bench.py@4d3ec3a:605` fails when `--label` has a `/` and
  `--mesh-dir` is given, because the mesh's parent directory is never
  created. The workaround was to create it first.
- **Not measured:** a profile of `features clip`, and any AC run.

## Addendum, 2026-09-29: re-timed after `d58d693` (battery, 64 → 63 %)

`features read` from the Europe file, before → after `d58d693`'s
primary-key fetch, median of 3 runs:

- Catchment at 1 m: 17.90 → 1.77 s.
- Catchment at 10 m: 19.08 → 0.98 s.
- 48 km square at 1 m: 17.72 → 1.01 s.

The totals fall from 30.3 to 14.0 s, from 28.2 to 9.9 s and from 48.0 to
30.6 s. The three meshes are byte-identical to the earlier runs.
**Verdict: the read defect is fixed; 16b-1/2 is ACCEPTED.**
`features clip` (8.2-8.5 s) is a known open item, pending Ola's decision,
and is not counted against the verdict. Details are in the README's
addendum.
