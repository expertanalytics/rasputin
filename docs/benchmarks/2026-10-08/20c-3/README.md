# 20c-3: short timed check and gate runs (2026-10-08)

Branch `worktree-soft-quality-3` at `745bcb0d` against master `483221d2`
(20c-2 merged). The design's "Speed judgment" check
(`docs/increments/20c-soft-quality.md`): Lagan and Numedalslagen once each
with the defaults, not a sweep, then the 20c-3 gate table's runs.

**Method.** Apple M1 Max, 10 cores, on AC power (charged), threads 10 (the
default), one run each, no warm-up, bounds checks on (libc++ fast) in both
trees. Branch: the worktree's own `.venv` and `_core` (`tin_engine.__file__`
in the worktree); master: the main checkout's `.venv`, rebuilt from master
`483221d2`. Arguments: the Lagan and Numedalslagen cases of
`docs/benchmarks/quick/cases.toml` on `worktree-perf-process` (tolerance 10 m,
CORINE 2018), plus `--out` and `--stats`. Times are the `**total**` row of
`--stats` (the `mesh` command body). Slivers are counted from the written
`.vtk` by `scripts/slivers.py` / `scripts/gate.py` (plan-view smallest angle
under 1°); their worst angles agree with `--stats`. Scripts: `scripts/`; raw
`--stats` reports: `stats/` (`branch-*`, `master-*` the defaults; `g-*` the
gate runs). Meshes stay out of the repository; rerun `scripts/run.sh`,
`gates.sh` and `off.sh` with `S` set to a scratch folder to regenerate them.

## Slivers with the defaults

| catchment | tree | triangles | slivers under 1° | worst angle |
|---|---|---|---|---|
| Lagan | master | 799 378 | 437 | 0.002384° |
| Lagan | branch | 799 408 | 186 | 0.1917° |
| Numedalslagen | master | 1 141 207 | 296 | 0.008623° |
| Numedalslagen | branch | 1 142 108 | 169 | 0.08145° |

Other quality, branch against master: largest height error unchanged
(Lagan 9.99993 m, Numedalslagen 9.99994 m); max vertex degree Lagan 14 (15),
Numedalslagen 19 (19). Delaunay check not run.

## Speed with the defaults (one run each)

| catchment | master total | branch total | change | `features clip` master | branch | of it, clean-up |
|---|---|---|---|---|---|---|
| Lagan | 17.64 s | 22.87 s | +5.24 s (+30 %) | 6.87 s | 12.13 s | 8.04 s |
| Numedalslagen | 6.81 s | 14.10 s | +7.29 s (+107 %) | 1.80 s | 8.84 s | 8.78 s |

**Hotspot (40 % or more of a run):** on the branch, `features clip` is 53 %
of Lagan's run and 63 % of Numedalslagen's; its sub-row `features clip:
clean-up` (the new land-cover stage) is 35 % and 62 %. The design estimated
the stage at 1.1 s per catchment; measured 8.0 s and 8.8 s. Which part of the
clean-up costs the time (repair, same-class merge, outline rule) is not
profiled; the gate runs below, which turn the same-class merge off, spend
3.0 to 6.3 s in it.

## Gate runs (branch, same-class merge off, every flag given)

| run | catchment | triangles | under 1° | share | worst | under 1° with a side < 10 cm | under 1°, centre within 20 m of outline | clean-up / total |
|---|---|---|---|---|---|---|---|---|
| tolerance 2 | Lagan | 742 931 | 151 | 0.0203 % | 0.006697° | 19 | 83 | 3.03 / 17.34 s |
| tolerance 2 + repair 5 cm | Lagan | 742 929 | 149 | 0.0201 % | 0.02865° | 17 | 83 | 3.11 / 17.60 s |
| tolerance 2 + outline 5 m | Lagan | 743 664 | 65 | 0.0087 % | 0.006697° | 2 | 4 | 6.16 / 20.60 s |
| + repair 5 cm | Lagan | 743 662 | 63 | 0.0085 % | 0.2320° | 0 | 4 | 6.25 / 20.65 s |
| tolerance 2 | Numedalslagen | 1 013 581 | 29 | 0.0029 % | 0.012446° | 4 | 25 | 2.71 / 7.24 s |
| tolerance 2 + repair 5 cm | Numedalslagen | 1 013 605 | 29 | 0.0029 % | 0.012446° | 4 | 25 | 2.95 / 7.46 s |
| tolerance 2 + outline 5 m | Numedalslagen | 1 016 427 | 15 | 0.0015 % | 0.08145° | 4 | 11 | 4.99 / 9.73 s |
| + repair 5 cm | Numedalslagen | 1 016 320 | 15 | 0.0015 % | 0.08145° | 4 | 11 | 5.29 / 9.74 s |

Against the gate table:

* Tolerance 2 alone: Lagan side < 10 cm 19 ≤ 25, share 0.0203 % ≤ 0.03 %,
  worst 0.006697° ≥ 0.005°; Numedalslagen share 0.0029 % ≤ 0.006 %, worst
  0.012446° ≥ 0.008°. **Met.**
* Outline rule at 5 m on top, Lagan: 65 ≤ 90, 4 near the outline ≤ 5,
  triangles +0.10 % ≤ +1 %, worst 1.000 × the tolerance-only figure. **Met.**
* Outline rule at 5 m on top, Numedalslagen: triangles +0.28 % ≤ +1 %
  (met), but **count 15 > 10, 11 near the outline > 5, worst 0.0814° < 0.1°:
  not met** (the design's figures were 4, 0 and 0.832°). "Near the outline"
  here is the sliver's centroid within 20 m of the domain outline; the design's
  own counting method was not checked against this one.
* Repair at 5 cm on top of each: every threshold of the run without it holds
  where that run held (Numedalslagen's outline-rule run fails the same three
  with or without it). Lagan's slit: corner 413 602.500 6 331 638.254 is not
  a mesh vertex (the nearest is the other corner, 1.0 cm away), so the two
  open corners are not both vertices. **Met.** The far-corner condition (no
  triangle with its smallest angle between two constraint edges at
  413 677.295 6 331 664.081) was **not checked**: that corner is a vertex in
  every run.
* All four off (`--features-repair 0 --features-outline-snap 0
  --features-tolerance 0 --no-features-merge-same-class`): points and cells
  bit for bit equal to master's default mesh on both catchments
  (`scripts/same.py`, from `POINTS` to the end of file). **Met.**
* Not run: the "no shared border inside the catchment left unmatched by the
  rule" check, and the "population-3 lines carry no constraint edge with the
  merge on" check.

## Re-time after rulings T1 and T2 (`35b7a873`, `8716f062`)

Apple M1 Max, AC power, one run each, defaults as above; master figures are
the run above at `16143f4c` (not re-run). Both commits change Python only
(`feature_input.py`), so the editable install's extension was not rebuilt.
Script: `scripts/t12.sh` (T1's `feature_input.py` checked out for its two runs,
then restored); raw reports: `stats/t12/`. Meshes stay out of the repository
(scratch `rasputin_scratch/perf-20c3-t12-20261008/`).

**T2 leaves the mesh unchanged:** `.vtk` SHA-256 on T1 and on T2 equal for both
catchments (Lagan `ef67e65c…`, Numedalslagen `1f11f291…`). **Met.**

| catchment | master total | T1 total | T2 total | T2 vs master | clean-up T1 | clean-up T2 | limit |
|---|---|---|---|---|---|---|---|
| Lagan | 17.64 s | 22.68 s | 20.59 s | +2.95 s (+17 %) | 7.98 s | 5.92 s | ≤ 6.5 s, met |
| Numedalslagen | 6.81 s | 13.67 s | 10.51 s | +3.70 s (+54 %) | 8.68 s | 5.48 s | ≤ 6.0 s, met |

**Phases at 40 % or more (T2):** `features clip` 48.5 % of Lagan's run;
`features clip` 52.9 % and its sub-row `features clip: clean-up` 52.2 % of
Numedalslagen's.

Slivers with the defaults, T2 (master as above):

| catchment | triangles | under 1° | worst angle | master under 1° | master worst |
|---|---|---|---|---|---|
| Lagan | 799 368 | 180 | 0.1917° | 437 | 0.002384° |
| Numedalslagen | 1 141 983 | 158 | 0.6302° | 296 | 0.008623° |

Max vertex degree Lagan 14, Numedalslagen 19; largest height error 9.99993 m
and 9.99994 m (unchanged). Delaunay check not run.

Outline-rule gate on T2 (`--no-features-merge-same-class --features-tolerance 2
--features-outline-snap 5 --features-repair 0`), against the tolerance-2 run
above (`g-tol2-*`, made before T1; T1 and T2 touch only the outline rule, so it
was not re-run):

| catchment | triangles (vs tolerance 2) | under 1° | centre within 20 m of outline | worst | verdict |
|---|---|---|---|---|---|
| Numedalslagen | 1 016 301 (+0.27 %) | 4 (≤ 10) | 0 (≤ 5) | 0.8322° (≥ 0.1°) | met |
| Lagan | 743 616 (+0.09 %) | 61 (≤ 90) | 0 (≤ 5) | 0.006697° (1.000 × tolerance-2) | met |

## Ruling R3: the gate's unrun checks (`aaee898b`)

Apple M1 Max, AC power, one session (04:17 to 04:23), the catchment arguments
of `scripts/gates.sh`. `aaee898b` changes Python only (`feature_input.py`,
`run_record.py` since `8716f062`), so the editable install was not rebuilt.
Scripts: `scripts/r3.sh` ((a), (b)), `scripts/r3de.sh` with
`scripts/r3drive.py` ((d), (e): wraps `snap_to_outline` and runs `rasputin
mesh` in-process), `scripts/farcorner.py` ((c)), `scripts/r3e_where.py` (the
(e) follow-up); `scripts/gate.py` and `scripts/corners.py` as before. Raw
output: `stats/r3/`. Meshes stay out of the repository (scratch
`rasputin_scratch/perf-20c3-r3-20261008/`; rerun the scripts to regenerate).

**(a) Defaults, mesh unchanged by R1: met.** `.vtk` SHA-256 Lagan
`ef67e65c…`, Numedalslagen `1f11f291…`, equal to T2's.

**(b) Repair with the outline rule** (`--no-features-merge-same-class
--features-tolerance 2 --features-outline-snap 5 --features-repair 0.05`),
against the tolerance-only meshes `g-tol2-*`:

| catchment | triangles (vs tolerance 2) | under 1° | centre within 20 m of outline | worst | verdict |
|---|---|---|---|---|---|
| Lagan | 743 614 (+0.09 %, ≤ +1 %) | 59 (≤ 90) | 0 (≤ 5) | 0.2320° (≥ 0.006362°) | met |
| Numedalslagen | 1 016 194 (+0.26 %, ≤ +1 %) | 4 (≤ 10) | 0 (≤ 5) | 0.8322° (≥ 0.1°) | met |

`corners.py` on Lagan: open corner 1 has no vertex within 1 mm (nearest
1.00 cm), open corner 2 has one; not both: **met**. Largest height error
9.99993 m and 9.99996 m. Max vertex degree and the Delaunay check not run.

**(c) The slit's far corner: met.** Both Lagan repair meshes (the new one
from (b) and `g-tol2rep-lagan.vtk`) have a vertex at 413 677.295
6 331 664.081; 6 triangles meet there, none with a smallest angle under 1°.

**(d) No shared border left unmatched: met.** `snap_to_outline`'s polygons
from the (b) runs, measured as `covfar.py` does: 0.000 m on both. The probe
fails when it should: one vertex of one shared border moved 0.5 m gives
208.170 m (Lagan) and 337.081 m (Numedalslagen).

**(e) Population 3 with the merge on: not met as worded.** Length of
`snap_to_outline`'s land-cover line segments in EPSG:3035 with both ends
within 1 m of N 3 811 923.31 or E 4 585 680.19, inside the domain, Lagan:

| run | length |
|---|---|
| defaults (merge on) | 26.194 m (expected 0) |
| `--no-features-merge-same-class` | 161 164.495 m (more than 0: the probe can fail) |

The 26.194 m is 13.097 m of border counted once for each of its two owners
(`stats/r3/e-where.txt`): four stretches of 0.69 to 9.46 m, each a border
between two polygons that the merge left separate, so of different classes
(merged polygons 10 and 12, and 10 and 11). Whether these are population-3
artefacts or land-cover borders that happen to lie within 1 m of the line has
not been measured.
