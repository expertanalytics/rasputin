# Increment 111: the cost of repairing land-cover polygons up front

Measured by `@perf`, 2026-10-29, for code review round 2's fix in
`docs/increments/111-one-land-cover-map.md` (last `## Review` entry):
`fa3164aa` makes `land_cover_map.build` check and repair each file's polygons
before they are joined (`land_cover_map._fixed`). Questions: where the added
time goes, how master repaired, and the cheapest form that keeps the crash on
a self-crossing polygon fixed.

## Method

- Machine: Apple M1 Max (8 performance and 2 efficiency cores, 32 GB), macOS
  27.0.1, Python 3.14.7, Shapely 2.2.0, GEOS 3.14.1. **AC power** before and
  after both runs (`pmset -g batt`, first and last lines of the `.txt` files);
  nothing else timed ran beside them.
- Quick check (`tools/bench_quick.py run --save-baseline`; `tools/bench.py`
  and `tools/bench_quick.py` are the same on both trees), back to back: master
  `9c317c16` (a scratch clone, its own `--features` cases, Ola's
  `~/.config/rasputin/sources.toml`), then the branch at `fa3164aa` (its
  `--land-cover ... --edges all` cases, `XDG_CONFIG_HOME` pointing at a copy of
  that file with its 15 `rules =` lines renamed `classes =`).
  `RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data`. Release `_core`,
  hardening on. Medians of 5 (tile, quarter, velhas) or 3 after one warm-up.
  The branch's stored baseline file was restored after the run (`git checkout`);
  the records are kept here as `master-9c317c16.json` and `branch-fa3164aa.json`.
- Profile: `probe/prof_fixed.py` runs `rasputin mesh` under cProfile with
  `_fixed`'s parts timed. `probe/dump.py` saves the polygons and the region
  `_fixed` receives (to the scratchpad, 87 MB and 255 MB, not kept; rerun it
  with the case's arguments to regenerate them); `probe/alts.py`,
  `probe/perpoly.py`, `probe/reuse.py` and `probe/proposal.py` time other
  forms on them (medians of 3 or 5). Their outputs are the `.txt` files beside
  them.

## Totals, seconds (median [min-max])

| Case | master 9c317c16 | branch fa3164aa | Change | Same mesh | Max error / tolerance |
|---|---|---|---|---|---|
| tile | 0.44 [0.44-0.45] | 0.44 [0.43-0.44] | -1.7 % | yes | 1.000 / 1 |
| tile t=1 | 0.74 [0.73-0.75] | 0.73 [0.72-0.73] | -1.4 % | yes | 1.000 / 1 |
| quarter | 0.38 [0.37-0.38] | 0.37 [0.37-0.38] | -1.3 % | yes | 1.000 / 1 |
| numedalslagen | 9.10 [9.08-9.27] | 9.44 [9.42-9.50] | +3.7 % | yes | 10.000 / 10 |
| lagan | 15.00 [14.87-15.86] | 15.30 [15.21-15.33] | +2.0 % | yes | 10.000 / 10 |
| geilo-al-ramp | 0.34 [0.34-0.34] | 0.34 [0.34-0.34] | -0.5 % | yes | 19.984 / 20 |
| romsdal-slope | 3.71 [3.70-3.75] | 3.67 [3.64-3.69] | -1.1 % | yes | 10.000 / 10 |
| velhas | 0.13 [0.13-0.13] | 0.13 [0.13-0.13] | -0.8 % | yes | 0.500 / 0.5 |

`bench_quick.verdict` on the branch against this master record says NO CHANGE:
master's numedalslagen spread this time (9.08-9.27, one slow run) widens the
band past 3.7 %. Against the stored baseline (8.83 s), the branch reads
numedalslagen +6.9 % and lagan +5.3 %, both SLOWER; master reads neither.
This is the same slowdown the main session saw (+6.0 %), not noise.

## Where the time goes (medians, s)

| Case | master: features read + clip (clean-up) | branch: land cover read + edges (clean-up) | Added |
|---|---|---|---|
| numedalslagen | 0.24 + 6.13 (6.07) = 6.37 | 0.62 + 6.07 (6.06) = 6.69 | +0.32 |
| lagan | 0.98 + 6.61 (5.63) = 7.59 | 2.29 + 5.63 (5.62) = 7.92 | +0.33 |

The whole difference is `_fixed`, inside `land cover read`. Profiled on
numedalslagen (7,932 polygons, 5.42 M vertices, none invalid): `_fixed` 0.325 s,
of it `shapely.clip_by_rect` to the widened box 0.085 s and `shapely.is_valid`
of the clipped parts (1.27 M vertices) 0.242 s. Seven clipped polygons over
10 k vertices take half of the validity time (CORINE 312 forest, up to 103 k
vertices, 32 parts and 655 holes); on lagan one polygon (174 k vertices,
1,290 holes) takes 0.13 s of 0.25 s (`probe/perpoly.txt`).

## How master repaired

- **The `--features` (edges) path checked nothing.** `_Tally._take` and
  `_clean` (`feature_input.py` at `9c317c16`) send each polygon straight to
  `shapely.intersection(g, region)`, which throws on a bowtie (`probe/throws.txt`);
  code review round 2 reports the crash there. So master's quick-check cases
  pay nothing because they repair nothing.
- **The `--land-cover` shares path repaired everything read, whole.**
  `read_class_polygons` calls `_repaired` on every polygon, after the move to
  the DEM's CRS and with no clip first; `_repaired` starts with the whole
  polygon's `geometry.is_valid`, and only an invalid one is cut to the box of
  the mesh's vertices widened by 100 m and rebuilt. On the same polygons that
  whole check takes 3.46 s (numedalslagen) and 11.75 s (lagan), against the
  branch's 0.32 s and 0.33 s (`probe/alts.txt`). The bowtie probe master passed
  ran through this path.

Which invalid shapes make GEOS's intersection throw (Shapely 2.2.0, GEOS 3.14.1):
a bowtie, a hole crossing its shell, overlapping parts and a zero-width spike
throw; a hole outside its shell does not, and gives a wrong area. So repairing
only where an overlay throws would leave silent errors in the class shares.

## Other forms, timed on the dumped polygons (median s)

| Form | numedalslagen | lagan |
|---|---|---|
| branch check: clip to box + `is_valid` | 0.316 | 0.331 |
| the same, 8 threads, 64 equal-count chunks | 0.195 | 0.241 |
| `is_simple` of the clipped boundary | 0.230 (739 false alarms) | 0.214 (720) |
| whole `is_valid` (master's shares path) | 3.458 | 11.745 |
| whole `is_valid`, 8 threads | 3.158 | 4.880 |

The clean-up's first step, `intersection(g, region)` per polygon, costs
0.615 s on the raw polygons and 0.501 s on the box-clipped ones (numedalslagen;
lagan 0.500 and 0.383 s), with the same vertex counts and equal areas
(`probe/reuse.txt`).

Check plus that intersection, together (`probe/proposal.txt`):

| Form | numedalslagen | lagan |
|---|---|---|
| no check (as master's edges path) | 0.617 | 0.494 |
| branch: check, then intersection of the raw polygons | 0.931 | 0.818 |
| proposed: check in 8 threads over 32 groups of equal vertex count, the clipped polygons kept and intersected | 0.597 | 0.590 |

The proposed form costs -0.02 s (numedalslagen) and +0.10 s (lagan, 0.6 % of
the run) against no check at all. Its clean-up inputs are not bit-identical:
719 and 735 polygons start their rings elsewhere, two differ by under
2e-8 m², and vertex counts are equal; whether the mesh hash changes is for
the quick check on the commit that makes it.

## Check

```bash
.venv/bin/python -c "import sys; sys.path.insert(0, 'tools'); import bench_quick as q; \
r = lambda p: q.QuickRecord.model_validate_json(open(p).read()); \
d = 'docs/benchmarks/2026-10-29/quick-111-repair/'; \
print(q.verdict(r(d + 'branch-fa3164aa.json'), r(d + 'master-9c317c16.json')).status)"
```
