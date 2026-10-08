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
