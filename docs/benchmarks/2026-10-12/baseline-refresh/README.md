# Quick-check baseline refresh on master `c317c7e4` and `969c49e2`, and one look at Lagan `features read` (@perf, 2026-10-12)

Asked by the main session: refresh `docs/benchmarks/quick/baseline-ac.json`
(stored on increment 35's head `357808b1`, 2026-10-11, master `3e6f87eb` plus
increment 35) on master `c317c7e4`, after #230, #231 (36), #232 (26a), #233,
#234 and #235 (26b) merged; and look once at Lagan's `features read`, which
read +12 to +17 % against that baseline on the main session's three idle runs
of 2026-10-12 (`docs/increments/26-land-cover.md`, section 14). The refresh
rule: `docs/increments/perf-quick-check.md`, section 3, *The baseline*, (b).

## Result in one paragraph

Two quick checks on master this afternoon (15:19 and 15:30 local, AC) read
**every** case 6 to 15 % slower than the stored baseline, including
single-thread C++ refine on byte-identical meshes (`tile t=1 refine` +15.5 %)
and a `_core.so` byte-identical to the one the main session's idle run used
at 11:15 (SHA-256 `55bb3993…`, both built from `c317c7e4`). The machine was
in interactive use (load average 3 to 4.5, WindowServer at 30 % CPU; no video
call), and the same `.so` was within band on every phase but one at 12:10 UTC,
so these runs measure the machine, not master. **The record saved by run 2 is
kept here as `baseline-ac-candidate-loaded-1530.json`, not installed as the
baseline**: a baseline 6 to 15 % slow on unchanged code would hide a
regression of that size until the next refresh. Lagan `features read` is
**not a code change**: the baseline's own commit and `.so`, run back to back
with master, reads 1.107 s against master's 1.103 s (medians of 3), both
against 0.922 s stored; the profile shows the same call counts and times on
both trees. See *Lagan `features read`* below. The idle refresh: *Run 3* found no idle minute in
the afternoon; **the evening's *Run 4* did it** (20:29 local, idle, AC) and
installed the new baseline, and *Run 5* is increment 38's quick check against
it: `NO CHANGE`, every mesh identical.

## Method

- Tree: `/Users/skavhaug/projects/rasputin/.claude/worktrees/baseline-1012`
  at `c317c7e4` (master), its own `build-bench` (Release, `-O3 -DNDEBUG`,
  hardening ON, AppleClang 21.0.0) and its own `.venv`, pinned to the stored
  baseline's machine record: Python 3.14.7, numpy 2.5.4, shapely 2.2.0 on
  GEOS 3.14.1, pyproj 3.8.0 on PROJ 9.8.1 (`uv pip freeze` of the
  `mapbiomas-26b` worktree's venv, `tin_engine` itself through `bench.py`'s
  `pkg` tree). The main checkout's `.venv` has numpy 2.5.3 and would give
  `NO BASELINE: machine numpy`.
- Command: `RASPUTIN_DATA=/Users/skavhaug/projects/rasputin_data
  .venv/bin/python tools/bench_quick.py run [--save-baseline]`, run from the
  worktree, nothing of ours beside it. Script version: `tools/bench_quick.py`
  and `tools/bench.py` at `c317c7e4`.
- Machine: Apple M1 Max (8 P + 2 E cores, 32 GB), macOS 27.0. Power: AC,
  battery 100 %, `pmset -g batt` before and after each run (in the logs).
  `pmset -g therm`: no thermal or performance warning recorded.
- Foreign load, checked before each run (`ps`, `top`, `uptime`): no video
  call; WindowServer 29 to 32 % CPU, loginwindow 13 %, iTerm2 11 %, Safari
  8 %, Slack idle; load averages 2.7 / 4.3 / 5.6 before run 1, 3.0 / 3.4 /
  4.3 before run 2, 4.5 after it. The machine was being used.
- Old tree for the profile: a scratch worktree of `357808b1` (the stored
  baseline's commit; its `io/`, `feature_input.py`, `crs.py` and `cli.py`
  are identical to master `3e6f87eb`'s, `git diff --stat 3e6f87eb 357808b1 --
  src_python/tin_engine/io src_python/tin_engine/cli.py` is empty), built the
  same way; its `.so` SHA-256 `51872f0a…` equals the `cleanup-speed`
  worktree's, which made the stored baseline. Removed after the run.

## Run 1 (15:19, judged against the stored baseline): `SLOWER`, exit 1

`quick-run-1.log`. Fourteen `SLOWER` lines, no `MESH CHANGED`, no `BROKEN`:

| line | stored | run 1 | change | band |
|---|---:|---:|---:|---:|
| tile t=1 total | 0.70 | 0.78 | +10.9 % | 5.0 % |
| tile t=1 refine | 0.47 | 0.52 | +11.2 % | 5.0 % |
| numedalslagen total | 8.83 | 9.48 | +7.3 % | 5.0 % |
| numedalslagen decode | 0.92 | 1.00 | +8.5 % | 5.0 % |
| numedalslagen features clip | 6.00 | 6.38 | +6.3 % | 5.0 % |
| numedalslagen features clip: clean-up | 5.94 | 6.31 | +6.3 % | 5.0 % |
| lagan total | 14.54 | 15.49 | +6.5 % | 5.0 % |
| lagan features read | 0.92 | 1.10 | +19.8 % | 5.5 % |
| lagan features clip | 6.34 | 6.69 | +5.5 % | 5.0 % |
| lagan features clip: clean-up | 5.39 | 5.68 | +5.4 % | 5.0 % |
| romsdal-slope total | 3.58 | 3.92 | +9.5 % | 5.0 % |
| romsdal-slope refine | 2.01 | 2.24 | +11.3 % | 5.0 % |
| romsdal-slope refine: scan (parallel) | 0.56 | 0.63 | +11.5 % | 5.0 % |
| romsdal-slope refine: split + flip (serial) | 1.18 | 1.33 | +12.5 % | 5.0 % |

## Run 2 (15:30, `--save-baseline`): `SLOWER`, exit 1; record not installed

`quick-run-2-save.log`; the record it wrote is
`baseline-ac-candidate-loaded-1530.json` (commit `c317c7e4`, power `ac`).
The full table against the stored baseline, every phase at 5 % or more of its
case's total, is `raw/delta-stored-vs-loaded-1530.md`; the totals:

| case | stored s | run 2 s | change | mesh SHA-256 |
|---|---:|---:|---:|---|
| tile | 0.416 | 0.453 | +8.9 % | same |
| tile t=1 | 0.699 | 0.789 | +12.8 % | same |
| quarter | 0.355 | 0.388 | +9.3 % | same |
| numedalslagen | 8.829 | 9.343 | +5.8 % | same |
| lagan | 14.542 | 15.483 | +6.5 % | same |
| geilo-al-ramp | 0.329 | 0.364 | +10.9 % | same |
| romsdal-slope | 3.578 | 3.876 | +8.3 % | same |

**Phases that moved more than their band, and why.** Every one. None ties to
a merge: the refine phases are C++ in a `.so` identical to the one within
band at 11:15 today, on meshes identical to the stored ones (`tile t=1
refine` +15.5 %, `romsdal-slope refine: split + flip (serial)` +12.9 %, both
serial C++); `decode` is tifffile and numpy, unchanged since the baseline
(`tile decode` +12 %); `features clip` is GEOS, unchanged (+5.7 % and +6.2 %).
A merge that slowed one phase would not slow every phase by about the same
share. The cause is the machine's state this afternoon (measured: the load
averages above; the same binary within band at 12:10 UTC per increment 26,
section 14). What in the machine's state, I did not measure.

The seven `HOTSPOT` lines (tile and quarter `refine`, `tile t=1 refine: scan
(parallel)`, Numedalslågen `features clip: clean-up`, `geilo-al-ramp decode`,
`romsdal-slope refine`) have been printed on every quick check since
2026-10-10 (`docs/benchmarks/2026-10-10/slope-tol/quick/head-run.log`);
nothing new. Numedalslågen `features clip` is ruled in `hotspots.toml`
(2026-10-12, 68 %) and prints nothing; its `clean-up` sub-row at 67 % is the
same work under the same ruling.

## Lagan `features read`: machine, not code (one look, 12 minutes)

What the phase does: `cli._open_features` times `feature_input.open_features`
and subtracts the clip (`src_python/tin_engine/cli.py@c317c7e4:1854`): the
CORINE GeoPackage
(8.9 GB, EPSG:3035) is queried through its rtree for the Lagan box
(`io/geopackage.query_features`, 7867 calls), each WKB blob decoded to
shapely (`decode_geometry`), reprojected to EPSG:3006 through pyproj
(`crs.reprojector`, `_transformer` builds the `Transformer`), and tested
against the domain. Of the files the task named, `io/repository.py` and
`io/models.py` changed since the baseline (`create_sqlite`,
`CachedObject.fetched`, `ClassWindow`, `ClassFractions`, `SrsEntry`), but
none of the new code is on this path: `open_geopackage` and `query_features`
are unchanged, and `crs.py`'s only change is the added `srs_entry` (36).

Measured, both trees back to back, alternating, same venv and inputs,
`tools/bench.py _child` without the profiler (`raw/lagan-plain-runs.jsonl`):

| tree | features read s (3 runs) | median | total s median |
|---|---|---:|---:|
| old `357808b1` (the baseline's own commit and `.so`) | 1.107, 1.133, 1.093 | **1.107** | 15.491 |
| new `c317c7e4` (master) | 1.131, 1.103, 1.102 | **1.103** | 15.541 |
| stored baseline (2026-10-11 evening) | | 0.922 | 14.542 |

The baseline's own code reads 1.10 s today. So the +14 to +20 % is not a
change in `io/repository.py`, `io/models.py` or anything else between the
two commits; it is the machine today against the machine on 2026-10-11
evening. The input fingerprints match (the tool would have said
`NO BASELINE: input` otherwise), so it is not the data.

Under cProfile (`raw/lagan-profile-read-path.txt`, two runs per tree):
`open_features` cumulative 8.52 / 8.51 s old, 8.53 / 8.50 s new;
`read_source` 1.07 / 1.09 s old, 1.11 / 1.07 s new; `query_features` 7867
calls, 1.05 / 1.07 s old, 1.09 / 1.05 s new; `decode_geometry` 7866 calls,
0.20 s on both. Same call counts, same times. The only sub-step that leaned
one way in both pairs is `crs._transformer` (pyproj `Transformer.from_crs`,
PROJ database work, about 60 primitive calls): 0.46 / 0.59 s old, 0.73 /
0.75 s new, under the profiler; it did not show in the unprofiled phase
above, and the function is the same text on both trees. Not a hotspot: the
phase is 7 % of the run (1.09 of 15.48 s), under the 40 % rule of section 5;
no `hotspots.toml` line is due. (Bookkeeping: the profiled runs passed
`--stats --binary`, and `--stats` took `--binary` as its report path, so they
wrote a text mesh and a report file named `--binary`, removed; all four used
the same arguments, so the old-new comparison holds. The unprofiled runs
above were `--binary`.)

## Run 3 (idle wait, then `--save-baseline`)

A waiter (`raw/wait_idle_then_save.sh`) sampled the 1-minute load average
every 30 s and was to run `--save-baseline` once it read under 1.5 twice in
a row, within 45 minutes. See `quick-run-3-idle-save.log`.

**It never ran.** From 15:34 to 16:19 the 1-minute load average stayed between
1.50 and 5.11 (46 samples, in the log), never under 1.5 on two samples in a
row; the waiter gave up at the cap, exit 4. **The stored baseline
(`357808b1`, 2026-10-11) stayed installed this afternoon; the evening's run 4
replaced it (below).**
The refresh is one command on an idle machine (about 4 minutes, nothing else
timed): from a checkout whose venv matches the stored machine record (numpy
2.5.4), `RASPUTIN_DATA=<data> .venv/bin/python tools/bench_quick.py run
--save-baseline` on master. Its judge lines against the stored baseline should
then show the refine phases within band; if they again read +6 to +15 % on an
idle machine, that is a finding about the machine since 2026-10-11, not about
master, and needs a look at the machine, not the code.

## Evening: run 4, the idle refresh, and run 5, increment 38's quick check

### Method

- The Swedish catchment batch finished at 18:05 UTC; Ola away. `origin/master`
  `f09b9863` merged into this branch as `969c49e2` (what master gained since
  the afternoon's `c317c7e4`: #236, increment 37, the MCP server, Python
  only; #237 the generality audit and #238 a retrospective, both prose).
  `_core` rebuilt Release into this tree's `build-bench`: SHA-256
  `55bb3993…`, byte-identical to the afternoon's and to the main session's
  11:15 build, so master's C++ is unchanged since `c317c7e4`.
- Same venv, machine record, command and data as the afternoon (above);
  script version `tools/bench_quick.py` and `tools/bench.py` at `969c49e2`
  (unchanged since `c317c7e4`). The waiter
  `raw/evening_idle_save_then_38.sh` sampled the 1-minute load every 30 s and
  started run 4 once it read under 2.0 twice in a row, then run 5 straight
  after; nothing of ours beside either. The 30 samples are in
  `quick-run-4-idle-save.log`: 28 at 2.03 to 3.37 from 20:15 to 20:28, then
  the two that fired, 1.96 and 1.76 at 20:29. `ps` before each run: no video call, no
  batch; iTerm2's `pidinfo` helper at 75 % of one core (it polls the agent
  panel), WindowServer no longer in the top six. Power: AC, battery 100 %,
  `pmset -g batt` before and after each run in the logs.
- Increment 38's tree: `/Users/skavhaug/projects/rasputin/.claude/worktrees/threads`,
  branch `worktree-threads` at `8f9289c4` (`--threads` and `--memory-budget`
  as settings; Python only). Its merge base with master is `5d2e706d`, master
  before #238 (prose), and `git diff --name-only 5d2e706d 8f9289c4` lists no
  C++ or CMake file; its `build-bench/_core…so`, built by the tool, is the
  same `55bb3993…`. Run through this tree's venv and `tools/bench_quick.py
  run --tree <threads>`, so it was judged against the baseline run 4 saved.

### Run 4 (20:29, `--save-baseline`, idle, AC): baseline installed

`quick-run-4-idle-save.log`, exit 2: `baseline saved`, then **`NO BASELINE:
input $RASPUTIN_DATA/sweden_glo30_cache changed`** — the Lagan case's GLO-30
cache directory gained tiles when the Swedish batch's fetch filled it, so its
fingerprint moved from `dir:3deaa80f…` to `dir:f00da080…`; the tool refuses
to judge across an input change, and the new baseline carries the new
fingerprint. The other nine inputs and the machine record are unchanged. The
installed record is `docs/benchmarks/quick/baseline-ac.json` (commit
`969c49e2`, power `ac`), copied here as `baseline-ac-new-969c49e2.json`.

Judged against the stored `357808b1` record as the tool would have, had the
fingerprint matched (`verdict` over the two records with the old inputs
substituted; the full table is `raw/delta-old-357808b1-vs-new-969c49e2.md`):
**one line beyond its band, no `MESH CHANGED`, no `BROKEN`**; every case's
mesh SHA-256 and `max_error` equal the stored ones.

| line | stored | run 4 | change | band |
|---|---:|---:|---:|---:|
| lagan features read | 0.92 | 1.03 | +12.1 % | 5.5 % |

The totals and the C++ phases:

| case | stored s | run 4 s | change | refine (change) | mesh |
|---|---:|---:|---:|---|---|
| tile | 0.416 | 0.430 | +3.3 % | +1.2 % | same |
| tile t=1 | 0.699 | 0.711 | +1.7 % | +0.2 % | same |
| quarter | 0.355 | 0.366 | +3.1 % | -0.4 % | same |
| numedalslagen | 8.829 | 8.894 | +0.7 % | +0.4 % | same |
| lagan | 14.542 | 14.736 | +1.3 % | (under 5 % of the run) | same |
| geilo-al-ramp | 0.329 | 0.344 | +4.7 % | +1.5 % | same |
| romsdal-slope | 3.578 | 3.624 | +1.3 % | +0.1 % | same |

So the afternoon's +6 to +15 % on every phase was the machine, as the
afternoon section says: idle, the same `.so` reads within 1.5 % of the stored
baseline on every refine phase. No merge since `357808b1` moved a phase.
`lagan features read` stays +12 % on an idle machine (1.030 to 1.071 s over
the three timed runs, so not a cache warming within the run); it is not code
(the baseline's own commit read the same 1.10 s back to back with master this
afternoon, *Lagan `features read`* above), and what in the machine's state
since 2026-10-11 makes the CORINE GeoPackage query slower by 0.11 s was not
measured. The small `decode` phases (+7 to +16 %, 0.007 to 0.010 s each) are
under the tool's 0.05 s floor. The seven `HOTSPOT` lines are the ones printed
on every quick check since 2026-10-10; nothing new.

### Run 5 (20:31, increment 38 `8f9289c4` against the refreshed baseline): `NO CHANGE`

`quick-run-5-increment-38.log`, exit 0. The tool's judge lines are the
seven standing `HOTSPOT` lines and `NO CHANGE`: no `SLOWER`, no `FASTER`,
no `MESH CHANGED`, no `BROKEN`, no `NOT JUDGED`. So every case's mesh
SHA-256 equals the baseline's (the tool prints `MESH CHANGED` for any case
whose hash differs), every `max_error` is within its tolerance, and no phase
at 5 % or more of its case's total moved beyond its band (5 %, or the
baseline's spread) and 0.05 s. The tool keeps no record without
`--save-baseline`, so run 5's per-phase figures are not stored; its verdict
is the evidence. Consistent with the design's claim of bit-identical meshes
and no timing change.

## Verdict

- Refresh: **done in the evening (run 4)**: `docs/benchmarks/quick/baseline-ac.json`
  is now master `f09b9863`'s code on `969c49e2`, idle, AC, with the Lagan
  GLO-30 cache's new fingerprint. Against the old `357808b1` baseline: no
  mesh changed, `max_error` within tolerance on every case, refine phases
  within 1.5 %; one line beyond its band, `lagan features read` +12.1 %
  (0.11 s), machine-side (see above), cause within the machine not measured.
  The afternoon's runs 1 and 2 measured a loaded machine, not a merge.
- Increment 38 (`8f9289c4`): **NO CHANGE** against the refreshed baseline;
  meshes identical on every case; no phase outside its band.
- Lagan `features read`: **machine, not code**; 7 % of the run, not a hotspot.

## Files

- `quick-run-1.log`, `quick-run-2-save.log`, `quick-run-3-idle-save.log`:
  the tool's lines with `pmset` and `uptime` before and after.
- `baseline-ac-old-357808b1.json`: the stored baseline as found.
- `baseline-ac-candidate-loaded-1530.json`: run 2's record, not installed.
- `raw/delta-stored-vs-loaded-1530.md`: stored vs run 2, every phase at 5 %
  or more of its case.
- `raw/lagan-plain-runs.jsonl`, `raw/lagan-profile-read-path.txt`,
  `raw/prof_mesh.py`: the features-read look.
- Evening: `quick-run-4-idle-save.log` (the load samples, then the refresh)
  and `quick-run-5-increment-38.log`; `baseline-ac-new-969c49e2.json`, the
  record installed as `docs/benchmarks/quick/baseline-ac.json`;
  `raw/delta-old-357808b1-vs-new-969c49e2.md`, old vs new, every phase at
  5 % or more of its case; `raw/evening_idle_save_then_38.sh`, the waiter
  that ran both.
- Meshes: each run's mesh went to a `rasputin-quick-*` temporary folder and
  the session scratchpad; regenerate with the `cases.toml` arguments.

## Review

Round 1, 2026-10-12, @reviewer (`c31e55e7..fe06568c`, 0 production lines): CHANGES REQUESTED, one wording, applied in this recording commit: /Users/skavhaug/projects/rasputin/.claude/worktrees/baseline-1012/docs/benchmarks/2026-10-12/baseline-refresh/README.md@fe06568c:171-172 said the old baseline stays installed; run 4 replaced it. Suggestion applied: the title names `969c49e2`. Times before the evening section are local (CEST, UTC+2). Verified: the refresh rule is section 3 (b); runs 1 and 2 match their logs; both `_core` hashes `55bb3993`; the waiter's 30 samples (14-minute wait); run 4 saved, then NO BASELINE for the GLO-30 cache; installed baseline equal to the new copy, machine record equal, all seven meshes and max errors equal; run 5 (increment 38) NO CHANGE; shipped tests pass; check_citations exit 0; merges clean onto master.
