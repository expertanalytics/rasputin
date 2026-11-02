# Increment 122 (lakes overrule land cover): section 8's checks before review

These checks compare master with the branch. They are **not** a timed
acceptance: the increment touches no refine or mesh code, so it needs no
`@perf` acceptance run (`docs/increments/122-lakes-overrule.md`, Status
line and section 8). The only timing figure here is the `land cover read`
phase, which section 8 predicts.

- **Master:** public `origin/master` `5921a878`. **Branch:**
  `worktree-lakes-overrule` `221e2cf7` (the full Python suite passes).
- **Machine:** Apple M1 Max, 10 cores, on AC power (`pmset -g batt`:
  "Now drawing from 'AC Power'", 100 %, charged). Nothing else ran
  alongside.
- **One `_core` for both sides.** `git diff --stat 5921a878 221e2cf7`
  changes no C++, so the branch's `_core` (rebuilt with
  `tools/rebuild_core.py`, sha256 `a353ac1efe06…`) was copied into a
  `git archive 5921a878 src_python` export. `scripts/run_master.py` takes
  the editable-install finder out of `sys.meta_path` and asserts that both
  `tin_engine` and `_core` load from the export.
- **How a mesh was compared** (`scripts/compare.py`, with the branch's
  `vtk_read.read_mesh` reading both sides): the vertex arrays are compared
  exactly, and so is their order; the triangles are compared as a set,
  each one's corners sorted; codes are compared triangle by triangle,
  matched through that set. "Marked" means the branch's
  `land_cover_filled` is 1 on that triangle.
- **Data**, all dated before the run: DTM10 `DTM10_UTM33_20260925`; the
  clean products `rasputin_data/clean/*_25833_20261101` (modified
  2026-11-01), called through Ola's `~/.config/rasputin/sources.toml` ids;
  NVE lakes `rasputin_scratch/temp-check-20261021/nve-lakes-25833.geojson`
  (2026-10-21, 201 polygons); the Bygdin domain
  `docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson`. The
  Numedalslågen and Lagan arguments are the ones in
  `docs/benchmarks/quick/cases.toml`. The Bygdin arguments are section 2's:
  `--tolerance 10 --edges all --threads 8 --prefer lake=<NVE>`.

## Verdicts

| # | Check (section 8) | Measured | Verdict |
|---|---|---|---|
| 1 | Numedalslågen: master's vertices and triangles; codes equal wherever master's is not 0 | 334,189 vertices and 653,037 triangles on both sides, identical (same order). Master has 0 triangles with code 0, and **every** code is equal. `land_cover_filled` 0 | PASS |
| 2 | Lagan: the same | 260,190 vertices and 463,342 triangles, identical. 0 triangles with code 0 on master; every code equal; `land_cover_filled` 0 | PASS |
| 3 | Bygdin, `clc2018-clean-25833 --prefer lake=` NVE: master's vertices and triangle set; codes equal except on the 145 island triangles; `land_cover_gaps` 0 | 30,367 vertices and 59,623 triangles, identical. 145 marked triangles, 0.1078 km² in 61 groups. Codes differ on 2 triangles only, both marked (master 322, branch 333). Unmarked triangles that differ: 0. `land_cover_gaps` 0 | PASS. Section 8 expected 3 changed from `_islands`; there are 2 |
| 4 | Bygdin, `ar50-innlandet-clean-25833 --prefer lake=` NVE: exit 2 before the mesh phase starts; 5 parts, 7.8634 km²; the largest 5.8199 km² around 136216 6830714 (±1 m); each part names `ar50-vestland-clean-25833` and no other id | Exit 2 after 5.6 s. No mesh, no `--stats` file and no record were written. The message lists 5 parts, 7.8634 km² in all, and every part names only `ar50-vestland-clean-25833`. The largest is 5.8199 km² **around 136214 6830715** | PASS on the exit, the parts, the areas and the ids. The place is **2.2 m off the predicted point**, outside ±1 m. This is a flaw in the prediction, not in the code (see below) |
| 5 | The same with `--land-cover-partial`: the gap line and the fill line, with 2,623 triangles, 8.006 km² and 92 groups, the largest 5.82 km² (±0.01); `land_cover_gap_m2` 7,863,437 (±1) | Exit 0. Gap line: "7.8634 km2 of the domain in no file (5 parts, the largest 5.8199 km2 around 136214 6830715 …)". Fill line: "2623 triangles (8.0061 km2) … 92 groups, the largest 5.8199 km2". Record: `land_cover_gaps` 5, `land_cover_gap_m2` 7863437.2, `land_cover_filled` 2623 | PASS |
| 6 | Innlandet plus `ar50-vestland-clean-25833`: exit 0 with `land_cover_gaps` 0 | Exit 0, `land_cover_gaps` 0, `land_cover_gap_m2` 0.0. 63,582 vertices and 125,973 triangles, the same as master. 394 marked islands; codes differ on 4 marked triangles (5051 and 5052 swapped, 2 each way), 0 unmarked | PASS |
| 7 | Numedalslågen: `land_cover_gaps` 0, and `land cover read` grows by under 1 s | `land_cover_gaps` 0. `land cover read`, median of 3 runs taken in turn with master after one warm-up: master 2.197 s (2.072, 2.197, 2.268), branch 1.974 s (2.402, 1.959, 1.974), so a change of −0.22 s, which is within the run-to-run spread. Wall clock, median: master 14.16 s, branch 14.15 s | PASS |

Also measured, though section 8 does not ask for it: the mesh from
Innlandet alone on master against the branch's `--land-cover-partial`
mesh. Vertices and triangles are identical. Master's 2,245 triangles with
code 0 are all marked on the branch and all take 5051. The other 4 changed
codes are marked islands, and no unmarked code changes
(`raw/compare-bygdin-innlandet-partial.txt`).

## The gap's place: why it is 2.2 m from the prediction

What runs: `coverage_gap` takes each part's place as the centre of its
largest inscribed circle, computed with
`shapely.maximum_inscribed_circle(part, GAP_NOISE / 10)`, a tolerance of
0.001 m (`land_cover_map.py`). Section 3.6 step 3 chose this tolerance
on purpose. Section 2's prototype used shapely's default tolerance, a
thousandth of the part's bounding box (several metres here). The
prediction of 136216 6830714 was copied from that prototype.

The probe (`scripts/gap_place_probe.py`, output in
`raw/gap-place-probe.txt`) runs the branch CLI, catches the domain and
the cover polygons that `coverage_gap` receives, and computes the largest
part both ways:

| Tolerance | Centre (EPSG:25833) | Distance to the gap's edge |
|---|---|---|
| 0.001 m (the code) | 136214.23 6830715.43 | 591.364 m |
| shapely's default (the prototype) | 136215.64 6830714.29 | 589.716 m |
| the predicted point 136216 6830714 | | 589.294 m |

The tiled computation and the untiled one give the same centre. So the
2.2 m gap comes from the tolerance, not from the tiles. The code's point
lies 2 m farther from the edge than the predicted one, so it is the better
answer to "the point of the part farthest from the data". The design text
should take the code's place: section 3.6's sample lines (136214 6830715,
and 141916 6830809 for the second part) and the figure in section 8.

## Method and raw data

- `scripts/all.sh` is the order the runs were made in, and
  `scripts/runs.sh` runs one of them (`--out .vtk --record --stats`).
  Both refer to the scratch folder of the session that made them. To
  rerun, set `S` to a scratch folder, export master into `$S/master`, and
  copy `_core` in as described above.
- `raw/<case>-<side>.json` holds each run's record, `.md` its `--stats`
  output, `.stderr.txt` its standard error (Python's start-up
  `RuntimeWarning` lines removed), and `.code` its exit status.
  The records and the stats name scratch-folder paths for `--out`.
- **The comparisons were rerun.** The first commit, 1ec9d178, held only
  tracebacks in its comparison files. Its loop ran under zsh, which does
  not split `$p` into words, so `compare.py` was handed file names with a
  space in them. The meshes had already been deleted, so
  `scripts/rerun.sh` (bash) ran each case once more, untimed: one master
  run and one branch run, then `compare.py` on each pair. Its outputs are
  `raw/compare-numedalslagen.txt`, `raw/compare-lagan.txt`,
  `raw/compare-bygdin-clc.txt`, `raw/compare-bygdin-innlandet-vestland.txt`
  and `raw/compare-bygdin-innlandet-partial.txt` (master with Innlandet
  alone against the branch with `--land-cover-partial`). Every figure in
  the verdict table matches them, and so do the rerun's records
  (`land_cover_gaps`, `land_cover_gap_m2`, `land_cover_filled`) and its
  refusal text. No verdict changed. The rerun's records, stats and stderr
  were not committed again. Every `.json`, `.md`, `.stderr.txt` and
  `.code` under `raw/` is from the first run (`scripts/all.sh`), the run
  that gave the timings.
- The meshes are not committed. To regenerate them, run
  `scripts/rerun.sh`.
