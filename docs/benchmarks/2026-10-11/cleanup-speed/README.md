# Increment 35 acceptance: the land-cover clean-up step, made faster (@perf, 2026-10-11)

Branch `worktree-cleanup-speed`, head `357808b1`, against its base, master
`3e6f87eb` (the merge of increment 32, `git merge-base master HEAD`). The
rule: `docs/increments/35-cleanup-speed.md`, section 7: the quick check
only (`tools/bench_quick.py run`), every mesh byte-identical to the base's,
Lagan `features clip` at most 0.65 of the base's median, Lagan's run at most
0.85, Numedalslagen `features clip` at most 1.0. No refine or mesh code is
touched, so no `bench.py` 1 m run and no thread sweep.

## Verdict: ACCEPTED

| gate | measured (medians of 3, AC) | limit | result |
|---|---|---|---|
| mesh SHA-256 unchanged, all 7 cases | equal on every case (`quick/delta.md`, column "mesh same") | equal | met |
| Lagan `features clip` | 11.156 s -> 6.342 s, ratio **0.569** | at most 0.65 | met |
| Lagan run total | 19.296 s -> 14.542 s, ratio **0.754** | at most 0.85 | met |
| Numedalslagen `features clip` | 6.092 s -> 6.000 s, ratio **0.985** | at most 1.0 | met |
| the other cases' totals | -0.7 % to +0.8 % | the case's band, at least 5 % | met |
| tool verdict (head against base) | `FASTER`, exit 0 (`quick/head-run.log`) | not `SLOWER`/`BROKEN` | met |

`max_error` is within tolerance on every case and the same on both sides
(records). Numedalslagen's clean-up stays near 6 s (5.94 s), as section 9's
question 2 expects.

The quick check's table, base `3e6f87eb` / head `357808b1`
(`scripts/quick_delta.py`, full output `quick/delta.md`):

| case | base total s | head total s | change | features clip | decode | refine |
|---|---:|---:|---:|---|---|---|
| tile | 0.417 | 0.416 | -0.3 % | - | 0.061 / 0.060 | 0.186 / 0.190 (+2.0 %) |
| tile t=1 | 0.704 | 0.699 | -0.7 % | - | 0.063 / 0.060 | 0.472 / 0.472 |
| quarter | 0.352 | 0.355 | +0.8 % | - | 0.044 / 0.044 | 0.171 / 0.172 |
| numedalslagen | 8.891 | 8.829 | -0.7 % | 6.092 / 6.000 (-1.5 %) | 0.921 / 0.921 | 0.660 / 0.659 |
| lagan | 19.296 | 14.542 | **-24.6 %** | 11.156 / 6.342 (**-43.1 %**) | 4.545 / 4.578 | 0.238 / 0.238 |
| geilo-al-ramp | 0.328 | 0.329 | +0.1 % | - | 0.139 / 0.138 | 0.100 / 0.099 |
| romsdal-slope | 3.582 | 3.578 | -0.1 % | - | 0.068 / 0.068 | 2.000 / 2.009 |

Spreads (max - min) / median: Lagan total 0.6 % base, 0.1 % head;
Lagan `features clip` 0.5 % / 0.1 %; Numedalslagen `features clip` 0.2 % /
0.2 %. The increment file's one-run figures (11.74 -> 6.74 s, 0.57) agree.

## Hotspots (`docs/benchmarks/quick/hotspots.toml`)

The ruled entry is Lagan `features clip` at 40.7 % (2026-10-07). Shares of
the case's total, head run:

| case | phase | base share | head share | tool line |
|---|---|---:|---:|---|
| lagan | features clip | 57.8 % | **43.6 %** | raised on the base run (grown 10 points or more over 40.7 %); not raised on head |
| lagan | features clip: clean-up | 36.7 % | 37.1 % | under 40 % |
| lagan | decode | 23.6 % | 31.5 % | under 40 % |
| numedalslagen | features clip | 68.5 % | 68.0 % | `HOTSPOT` (not in the file) |
| numedalslagen | features clip: clean-up | 67.8 % | 67.2 % | `HOTSPOT` (not in the file) |

Lagan's `features clip` share had grown to 57.8 % on master `3e6f87eb`
(increment 32 made the refine phase there much smaller, 0.65 s -> 0.24 s,
while the clip stayed); increment 35 brings it to 43.6 %, within 10 points of
the ruling, so the tool is quiet on head. Numedalslagen's `features clip` has
printed `HOTSPOT` on every quick check since the first (53 % on 2026-10-10,
`2026-10-10/slope-tol/quick/head-run.log`); it is unruled, and its share grew
with increment 32 for the same reason (refine 1.19 s -> 0.66 s). The other
lines (`tile` and `quarter` refine, `tile t=1` scan, `geilo-al-ramp` decode,
`romsdal-slope` refine) are the same unruled lines as the last two checks.
Both unruled sets go to Ola.

## Method

- `RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data`, this worktree's
  `.venv` (Python 3.14.7, numpy 2.5.4, shapely 2.2.0, GEOS 3.14.1, pyproj
  3.8.0, PROJ 9.8.1), Apple M1 Max, macOS 27.0, AC power before and after
  every run (`quick/power-*.txt`); nothing else ran.
- **First run** (`quick/head-run1-vs-39e3ae34.log`): head against the stored
  baseline `39e3ae34` (increment 34's). The tool said `NO BASELINE: machine
  numpy: 2.5.4 vs 2.5.3` (this venv's numpy is newer than the one that made
  the baseline) and judged nothing; its record is
  `quick/head-357808b1-run1.json`, and `quick/delta-vs-39e3ae34.md` compares
  the two by hand (`scripts/quick_delta.py base-39e3ae34.json
  head-357808b1-run1.json`, with the old baseline from `git show
  357808b1:docs/benchmarks/quick/baseline-ac.json`, not kept here): Numedalslagen's and Lagan's meshes differ from `39e3ae34`'s
  and their refine phases are 45 % and 64 % shorter. Increment 32 (the border
  simplifier, merged as `3e6f87eb` after that baseline was saved) changes the
  clean-up's output, so those differences are master's, not this branch's.
- **Back to back**, as `docs/increments/README.md` ("Acceptance") says for a
  `NO BASELINE`: `bench_quick.py run --tree <scratch worktree of 3e6f87eb>
  --save-baseline` (record `quick/base-3e6f87eb.json`, log
  `quick/base-run.log`, judged against `39e3ae34`: `SLOWER` on Lagan, which is
  the numpy-matched measure of master's own drift since increment 34), then
  `bench_quick.py run --save-baseline` at head (record
  `quick/head-357808b1.json`, log `quick/head-run.log`): `FASTER`, exit 0.
  Both the Release `_core` builds came from `bench.build` into each tree's
  `build-bench`; the Python package ran from each tree.
- **Baseline refreshed**: `docs/benchmarks/quick/baseline-ac.json` is now the
  head record (`357808b1`, numpy 2.5.4), per `perf-quick-check.md` section 3
  rule (a) (a `FASTER` accepted, at the branch's head). A checkout whose venv
  has numpy 2.5.3 will see `NO BASELINE: machine numpy` until it upgrades.
- Meshes went to the tool's temporary directory and are not kept; the runs
  regenerate them. Wall clock: 20:29 to 20:38, three quick checks.
