# Quick check, increment 111 (one land-cover map): branch against master

Measured by `@perf`, 2026-10-28, for
`docs/increments/111-one-land-cover-map.md` section 6 (*Gate*): no run that
does not use `--features` gets slower, and a run following every class border
(`--land-cover ... --edges all`, no `--shares`) costs what the `--features`
CORINE run costs on master.

## Method

- Machine: Apple M1 Max (8 performance and 2 efficiency cores, 32 GB), macOS
  27.0.1, Python 3.14.7, NumPy 2.5.4, Shapely 2.2.0, GEOS 3.14.1. **AC power**
  before and after both runs (`pmset -g batt`); nothing else timed ran beside
  them.
- Script: `tools/bench_quick.py` and `tools/bench.py` at `7974b15f`, the same
  on both trees. Each run builds Release `_core` into its own tree's
  `build-bench` (hardening on, the default) and times the cases in that tree's
  `docs/benchmarks/quick/cases.toml`: master's own `--features` commands on
  master, the `--land-cover ... --edges all` commands on the branch. Inputs
  are the same files (fingerprints equal in both records).
- Runs back to back: master `9c317c16` (a detached scratch worktree), then
  branch head `c42645f8`. Medians of 5 (tile, quarter, velhas) or of 3 after
  one warm-up (the rest); threads 0 (the CLI default) except `tile t=1`.
- `RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data`. On the branch,
  `XDG_CONFIG_HOME` pointed at a copy of Ola's `~/.config/rasputin/sources.toml`
  with its 15 `rules =` lines renamed `classes =`, as the design's section 2
  has the main session do after the merge; the unedited file is refused by
  the branch, and every `rasputin mesh` with it (the first branch attempt
  stopped on the `tile` case with `source.ar50-innlandet.rules: Extra inputs
  are not permitted`).
- Records: `master-9c317c16.json` and `branch-c42645f8.json` (the
  `--save-baseline` output, kept here; the stored baseline was not changed).
  Printed output: `master-9c317c16.txt`, `branch-c42645f8.txt`. Meshes were
  written to a temporary folder and not kept; rerun the two commands to
  regenerate them.

## Totals, seconds (median [min-max])

| Case | master 9c317c16 | branch c42645f8 | Change | Same mesh | Max error / tolerance |
|---|---|---|---|---|---|
| tile | 0.44 [0.43-0.45] | 0.44 [0.43-0.65] | -0.5 % | yes | 1.000 / 1 |
| tile t=1 | 0.73 [0.73-0.73] | 0.73 [0.73-0.74] | +0.5 % | yes | 1.000 / 1 |
| quarter | 0.38 [0.37-0.38] | 0.37 [0.37-0.38] | -1.4 % | yes | 1.000 / 1 |
| numedalslagen | 9.06 [9.06-9.08] | 9.10 [9.08-9.12] | +0.4 % | yes | 10.000 / 10 |
| lagan | 14.88 [14.87-14.89] | 14.82 [14.81-14.85] | -0.4 % | yes | 10.000 / 10 |
| geilo-al-ramp | 0.34 [0.34-0.35] | 0.34 [0.34-0.34] | -0.4 % | yes | 19.984 / 20 |
| romsdal-slope | 3.64 [3.64-3.70] | 3.68 [3.67-3.75] | +1.0 % | yes | 10.000 / 10 |
| velhas | 0.13 [0.13-0.13] | 0.13 [0.13-0.14] | +2.1 % | yes | 0.500 / 0.5 |

Branch judged against the master record with `bench_quick.verdict`: **NO
CHANGE** (no measure outside its band, 5 % or the spread). "Same mesh" is the
quick check's mesh hash, from `POINTS` on.

## The two land-cover cases, by phase (medians, s)

| Case | master: features read / features clip / its clean-up | branch: land cover read / edges / its clean-up |
|---|---|---|
| numedalslagen | 0.21 / 6.12 / 6.06 | 0.29 / 6.04 / 6.04 |
| lagan | 0.92 / 6.57 / 5.63 | 1.85 / 5.62 / 5.62 |

On Lagan the vertex move left the clip for the read (section 8, *Phases*):
read plus clip is 7.49 s on master and 7.47 s on the branch.

## Against the stored baseline (`quick/baseline-apple-m1-max-ac.json`, tree d886847d)

Both runs print the same one line: `SLOWER: numedalslagen decode 0.90 ->
0.95 s` (+5.6 % on master, +5.8 % on the branch, band 5.0 %). Master shows it
too, so it is not this branch's; the DEM decode is not code this increment
touches. Hotspot lines are the same on both, with names changed: the
branch's `numedalslagen edges: clean-up 66 %` replaces master's
`features clip: clean-up 67 %` (neither is ruled; the ruled `edges` at 68 %
is not raised, at 66 %).

## Check

```bash
.venv/bin/python -c "import sys; sys.path.insert(0, 'tools'); import bench_quick as q; \
r = lambda p: q.QuickRecord.model_validate_json(open(p).read()); \
d = 'docs/benchmarks/2026-10-28/quick-111/'; \
print(q.verdict(r(d + 'branch-c42645f8.json'), r(d + 'master-9c317c16.json')).status)"
```
