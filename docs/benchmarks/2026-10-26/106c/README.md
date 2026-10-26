# 106c acceptance: `mesh --catchments` (@perf, 2026-10-26)

Head `worktree-demo-fixes-c` at `423b74cd` (106a, 106b and 106c on master
`f17b3e11`). Base: master `79535812`, from a detached scratch worktree, run
back to back with the head. The rule is the `@perf` paragraph of
`docs/increments/106-wednesday-cells.md`, section 4 (106c).

Apple M1 Max (8 P + 2 E cores), 32 GiB, macOS 27.0.1, Python 3.14.7,
shapely 2.2.0 (GEOS 3.14.1). Every run on AC (`pmset -g batt`: "Now drawing
from 'AC Power'", 100 %), read before and after. 106c has no C++; the
head's and master's Release `_core` are the same bytes (sha256 `a353ac1e…`),
installed into the worktree's venv before the catchment runs. No sanitizer
run: no C++ patch. Desktop applications (Teams, a browser) were running
throughout; nothing else of mine ran beside a timing run, except that the
first seconds of the Ljungan run overlapped a short check of Kultsjön's
mesh (its `catchments` phase was then re-timed alone, table 2).

## Verdict

**ACCEPTED** for 106c: the 1 m benchmark's meshes are byte-identical to
master's in the same time; both catchments give the predicted triangle
counts and keep every cell within its unit to the predicted share; the time
outside the land-cover step is within the prediction plus 10 %.

**Found beside it, not 106c's: a 4.3-fold slower land-cover step from 106a**
(same branch). Kultsjön, `--land-cover clc2018 --triangles 6000`, no
`--catchments`: "land cover fractions" 10.43 s on master, 45.20 s on the
head. Profiled (cProfile, `profile-land-cover-head.txt`): `shapely.make_valid`
in 106a's `_repaired` (`src_python/tin_engine/feature_input.py`) takes
33.19 s for the 43 land-cover polygons it repairs, 0.77 s each. Ljungan's
land-cover step is 54.5 s at 80,000 triangles. Not timed on master for
Ljungan.

## 1. The 1 m benchmark (`tools/bench.py`, `--hardening on`)

| run | commit | tile hash | quarter hash | tile t0 | tile t1 | tile t12 | tile t20 | tile ceiling | quarter t0 | quarter t1 | quarter t20 | quarter ceiling | verdict against base |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `106c-head-ac-r1` | 423b74cd | 11741a81 | a60fb597 | 0.2753 | 0.6558 | 0.2729 | 0.2710 | 2.42x | 0.1973 | 0.5388 | 0.1951 | 2.76x | ACCEPTED (`bench.py compare`) |
| `106c-base-ac-r1` | 79535812 | 11741a81 | a60fb597 | 0.2758 | 0.6509 | 0.2684 | 0.2716 | 2.40x | 0.2044 | 0.5655 | 0.2105 | 2.69x | (base) |
| `106c-head-ac-r2` | 423b74cd | 11741a81 | a60fb597 | 0.2800 | 0.6729 | 0.2839 | 0.2838 | 2.37x | 0.2038 | 0.5209 | 0.1925 | 2.71x | 3 of 42 counts over 5 % (tile t14-t16, +6.2 to +8.6 %) |

Refine medians in seconds, 5 repeats; ceiling is t1 over t20. The meshes
are the hashes of 2026-10-25 (104b). Quality as there: worst angle
0.6296° / 0.3955°, max degree 74 / 18, within tolerance, 0 Delaunay
violations of 692,056 / 641,792 edges. r2's three counts over 5 % come from
the same `_core` bytes as the base, so they are noise of the day.

Against 2026-10-25's stored AC run (`104b-head-ac-r2`, the same `_core`
bytes), `bench.py` prints REGRESSION on every count, +11 to +48 %; master
measured today is as slow (tile t1 0.6509 against 0.4770). The machine was
slower today for both, thought to be the desktop load; not profiled.

## 2. Kultsjön and Ljungan with `--catchments smhi-svar2022-subcatchments`

Command (section 6 of the increment, without `--reservoirs`, which is 106d):
`rasputin mesh --out OUT.gpkg --dem glo30 --out-crs EPSG:3006 --domain
~/rasputin-output/opal-se2-2026-10-21/catchments/NAME.geojson --catchments
smhi-svar2022-subcatchments --land-cover clc2018 --triangles N --stats
OUT.stats.json`, one run each. `--land-cover` is required by `--catchments`.

| catchment | `--triangles` | band | border vertices | triangles (predicted) | largest share of a cell outside its unit (predicted) | ids wrong | worst angle | max degree |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Kultsjön, 75 units | 6,000 | 101.19 m | 2,701 | **6,057** (6,057) | **8.06e-5** (8.1e-5) | 0 | 1.09° | 10 |
| Ljungan, 809 units | 80,000 | 76.10 m | 36,970 | **78,046** (78,046) | **3.00e-4** (3.0e-4) | 0 | 0.995° | 12 |
| Kultsjön, `--features-tolerance 200` | 4,000 | 200 m | | **3,972** (3,972) | | | 0.535° | 10 |

Every mesh within its tolerance (budget sentence in the logs). Kultsjön: no
cell over 1e-4 outside its unit; Ljungan: 3. Every unit's id appears;
5 to 726 cells per unit (Kultsjön), 2 to 946 (Ljungan). No "centroid in no
unit" line. The share is checked by `straddle_units.py` against the units
rebuilt by the CLI's own `_units` at the recorded band; "ids wrong" counts
cells whose `catchment_id` is not the unit holding most of the cell.

**Times** (seconds; the `--stats` phase table, and `units_phases.py` for the
split, run alone):

| | Kultsjön 6,000 | Ljungan 80,000 |
|---|---:|---:|
| reading and selecting | 0.14 | 0.46 |
| repairing (and the union's area) | 0.59 | 6.65 |
| simplifying | 0.57 | 6.08 |
| `catchments` phase in the run | 1.30 | 13.65 |
| everything else but land cover | 0.87 | 28.99 |
| run total less "land cover fractions" | **2.18** | **42.64** |
| prediction: mesh + read and repair + simplify, plus 10 % | (1.2 + 0.3 + 0.6) x 1.1 = 2.31 | (30.9 + 3.0 + 6.1) x 1.1 = 44.0 |
| land cover fractions (not in the prediction) | 45.20 | 54.52 |
| run total | 47.38 | 97.16 |

The three figures the design left unmeasured, now set from these runs:
the largest share outside a unit, 8.06e-5 and 3.00e-4; reading and
repairing, 0.73 s and 7.11 s (the design's 0.3 s and 3.0 s are under it,
the repair being the larger part); simplifying, 0.57 s and 6.08 s.

Land-cover step, Kultsjön 6,000 without `--catchments`:
`kultsjon-6000-master.stats.md` (master's package, 10.43 s of 11.45 s) and
`kultsjon-6000-nocatch.stats.md` (head, 45.23 s of 46.24 s). Same 5,720
triangles.

## Files

- `../106c-head-ac-r1/`, `../106c-base-ac-r1/`, `../106c-head-ac-r2/`:
  `bench.py` output (`run.json`, `raw.tsv`, `README.md`), bench.py blob
  `16ad2ff5`.
- `*.stats.md`, `*.log`: each catchment run's `--stats` output and stderr.
- `straddle_units.py`, `units_phases.py`: the checks above, with their
  commands in their docstrings.
- `profile-land-cover-head.txt`: the cProfile summary of the head's
  Kultsjön run without `--catchments`.

The meshes were written to `rasputin_scratch/perf-106c-out/` and removed;
rerun the commands above to regenerate them.
