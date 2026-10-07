# Increment pq: the 15-minute performance check

**Status:** designed by `@architect` on `85a3e6dd`, 2026-10-07; read by
`@perf` the same night and corrected (section 8); next, Ola reads it.
Tooling only: `tools/bench.py` and a new `tools/bench_quick.py`. One PR. No
refine or mesh code, so no acceptance run of its own. Its first baseline is
taken after 30d (the outline-buffer fix, branch `worktree-buffer-speed-2`)
merges. Deferred to a later PR: the timing lock (section 6) and the `profile`
subcommand (section 7).

## 1. What Ola asked

Ola, 2026-10-07: "Make architect and perf design a process to measure
performance using maximum 15 minutes. Every now and then, more extensive tests
should be run. But wasting hours on things that could be checked in minutes is
not really very impressive." And: "I want architect and perf to design a system
that is lean and mainly catches regressions. Also, if anything takes a very
large portion of the run time, this hotspot should be recorded and discussed
with me."

Why. The outline-buffer question (`docs/benchmarks/2026-10-07/30d-buffer/README.md`,
branch `worktree-buffer-speed-2` at `d9f9c025`) took 89 minutes and needed one
profile of one run. `tools/bench.py` judges only `refine_s`, so it cannot see a
slower `decode` or `features read` at all. 15f-3's acceptance found `app_s`
+100 % with `refine_s` within noise (`docs/benchmarks/2026-10-04/15f-3-acceptance.md`).

## 2. Prior art: legacy and literature

*Literature.* Kalibera and Jones, "Rigorous benchmarking in reasonable time",
ISMM 2013 (doi:10.1145/2464157.2464160): choose repetitions from measured
variance, and report effect sizes. Followed in spirit: the noise band below is
each case's own measured spread, with a fixed floor; no confidence interval,
because 3 to 5 repeats do not support one. Daly et al., "The use of change
point detection to identify software performance regressions in a continuous
integration system", ICPE 2020 (doi:10.1145/3358960.3375791): a series of
runs per commit and change-point detection on it. Not used: we have one
machine and a run per increment, not a dense series; a stored baseline and a
band is the method that fits. No novelty is claimed.

*Legacy.* Nothing. `git grep -n -i -E "benchmark|cProfile|timeit|perf_counter|time\.time" legacy-archive -- legacy`
returned one line, `legacy-archive:legacy/rasputin/avalanche.py:26`, a
`datetime.timedelta` (a false match); the legacy tree has no timing code.

## 3. The quick check

`python tools/bench_quick.py run [--tree DIR] [--budget 900] [--save-baseline]`.

**When it runs** (Ola's ruling on question 3, section 10). There is no blanket
trigger. Each increment's design says in one line whether the change could
change speed, and why. It could when it changes what work `rasputin mesh` does:
a step turned on by default, how files are read or written, refine, decode or
feature reading. A new flag off by default, a reworded error or help text is
"no effect on speed". On yes, the main session runs the quick check before the
push; on no, nothing runs. The monthly extensive run (section 4) catches what
that judgment misses.

**Cases**, fixed, in `docs/benchmarks/quick/cases.toml`, run in this order:

| case | inputs | threads | runs | per run (`@perf`, 2026-10-07) |
|---|---|---|---|---|
| `tile` | bench.py's 1 m tile (`tests/fixtures/dem_archive/7908_3_10m_z33.tif`), `--tolerance 1` | default and 1 | 5 | 0.92 s; 1.21 s on 1 thread |
| `quarter` | same DEM, `docs/benchmarks/2026-09-26/quarter.geojson` | default | 5 | 0.80 s |
| `numedalslagen` | DTM10, NVE outline, `corine2018_dtm10_utm33.gpkg`, `--tolerance 10` | default | 1 warm-up + 3 | 6.97 s |
| `lagan` | GLO-30 cache, SMHI outline, the European CORINE GeoPackage (`U2018_CLC2018_V2020_20u1.gpkg`, 8.9 GB), `--out-crs EPSG:3006`, `--tolerance 10` | default | 1 warm-up + 3 | 16.97 s after 30d |

The catchments are two of 30d's cases, one per feature path that a real run
takes: a pre-cut national GeoPackage in the DEM's CRS, and the continental
GeoPackage queried and reprojected. `lagan` replaces the earlier `ljungan`
(GeoJSON window): the window is a hand-cut stopgap, and Ljungan's 53 % `decode`
hotspot is gone once 30d merges; `lagan` carries the hotspot that remains
(`features clip`, 40.7 % of 16.97 s,
`docs/benchmarks/2026-10-07/30d-buffer/built.md` on 30d's branch). Neither
case needs the network (meshing is offline, `sources.py:6`). Every run is
`--binary`. Paths are relative to `$RASPUTIN_DATA` (default
`../rasputin_data`); before the first baseline, the Numedalslågen outline is
copied there from
`rasputin_scratch/norway/numedalslagen/numedalslagen_outline_nve.geojson`
(a results folder).

**What is measured.** Per run: the process wall time, and each `--stats` phase
row with the clock's total. The child (`bench.py _child`) captures these by
replacing `cli.PhaseClock` with a subclass that keeps its instance, and adds
`phases` and `total_s` to its `BENCH` line (old records still parse). Also
`max_error` and the mesh's SHA-256 from the last run of each case.

**The baseline** is one file per power state, `docs/benchmarks/quick/baseline-ac.json`
and `baseline-battery.json`, a run of the check on master: commit, machine,
input fingerprint (for a file: path, size, mtime, and content SHA-256 under
100 MB; for a directory, such as DTM10's 762 files or the GLO-30 cache: the
SHA-256 of its sorted relative file names, sizes and mtimes, recursive, no
content read),
and per case and measure the median, min and max. A machine, power state or
input mismatch is `NO BASELINE: <field>`, never a cross comparison. It is
refreshed with `--save-baseline` (a) in the PR whose check showed a change that
was accepted (`FASTER`, or a `SLOWER` Ola accepted), measured at that branch's
head, (b) after a `NO BASELINE`, on a master checkout, and (c) by the monthly
extensive run.

**Regression.** Per case, for the total and for each phase at 5 % or more of the
baseline's total: `SLOWER` when the new median exceeds the baseline's by more
than the band *and* by more than 0.05 s. The band is the larger of 5 % (Ola's
2026-09-27 ruling for bench.py) and the baseline's own (max − min) / median;
the new run's spread does not widen it, so one slow outlier in the check cannot
hide a slowdown. Symmetric for `FASTER`. A `max_error` over the tolerance is
`BROKEN`. A changed mesh hash is reported, not judged.
*Scale:* phases of 0.1 s and up, runs of 1 to 20 s. `@perf`'s check on 30d's
master runs (3 each, 7 to 103 s): totals spread 0.5 to 4.4 %, so 5 % holds for
them; small phases spread more (Ljungan `decode` 5.7 to 8.2 %, Numedalslågen
`decode` 12.8 %, its `features read` 37.5 % from one fast run). A phase with a
wide baseline spread stays loosely judged until the next baseline; its
slowdown is still caught by the total's 5 % band once it moves the total by
5 %. bench.py's `refine_s` spread was 0.3 to 5.9 %.

**Report.** One line per finding, then the verdict:

```
SLOWER: lagan features clip 6.89 -> 8.10 s (+17.6 %, band 5.0 %)
HOTSPOT: lagan features clip 41 % of the run (not in hotspots.toml)
NO CHANGE | FASTER | SLOWER | BROKEN | NO BASELINE: <field> | OUT OF TIME: <cases not measured>
```

Exit 0 (no change, faster), 1 (slower or broken), 2 (no baseline), 3 (bad
input), 4 (out of time). A hotspot does not change the exit.
Evidence goes to the scratchpad; it is committed under `docs/benchmarks/<date>/`
only with a finding or a baseline refresh.

**The 15 minutes is enforced in the tool.** The deadline is the command's start
plus `--budget` (default and maximum 900 s; more is refused). Every subprocess
(CMake and each child) gets `timeout = deadline − now − 10 s`, through a new
`timeout` argument of bench.py's `Runner.run`; on expiry the child is killed.
Before each case, if the baseline's median times its runs exceeds what is left,
the case is skipped rather than started. Either way the verdict is
`OUT OF TIME` naming the cases not measured, and what was measured is still
judged. Estimate (`@perf`, 2026-10-07): the Release build of `_core` 18 s
fresh, 0 s when nothing changed; the cases about 2 min (tile and quarter 15 s,
Numedalslågen 28 s, Lagan 68 s, plus start-up). About 2.5 minutes in all, so
the 15 minutes is a guard, not a target.

## 4. Extensive runs

The existing acceptance (`bench.py run`: the 1 m benchmark, the thread sweep,
the catchment runs an increment names). They run (a) before the push of an
increment that touches refine or mesh code, after its quick check, as today;
(b) on master on the first unattended night of each month, on AC, which also
refreshes the baseline; (c) on Ola's ask. They never answer "is it slower"
(the quick check does) or "why is it slow" (one profile does, section 7).
Trimming the sweep is left to `@perf`, with evidence, outside this increment.

## 5. Hotspots go to Ola

A **hotspot** is a `--stats` phase row at 40 % or more of a case's total, from
the quick check, or a call whose own time is 25 % or more of the run, from a
profile (section 7). *Scale:* runs of 1 to 100 s; at 40 %, `lagan` after 30d
(`features clip`, 40.7 % of 16.97 s) is flagged, and so were 30d's master
Lagan (`decode` 49 %, `features read` 41 % of 72.7 s) and Ljungan (`decode`
53 %); Numedalslågen's largest phase (`features clip`, 26 % of 6.97 s) is not.

**One place:** `docs/benchmarks/quick/hotspots.toml`, one entry per ruled
hotspot: case, phase or call, share, date, Ola's ruling word for word. The tool
prints `HOTSPOT:` for one not in the file, or grown 10 points or more beyond its
ruled share. `@perf` (or the main session) turns each into an `ASK OLA:` line;
after Ola rules, the entry is added in the next commit that touches the file.
The first quick run will raise `lagan`'s `features clip` (already reported by
30d's `@perf`), and may raise the refine phases on `tile` (not measured); one
ruling each then silences them.

## 6. The timing lock: deferred

Not in this PR (`@perf`: it is not needed to catch regressions, and it is about
30 lines). Meanwhile the present rule stays in force: nothing runs beside a
timing run. `tools/brief.py` refuses to brief anyone beside `@perf`
(`tools/brief.py@85a3e6dd:162-163`), and the main session runs the quick check
itself only while no subagent is running. The later PR, if Ola wants `@perf`'s
write-up to overlap other work, is the design that stood here at `3c885768`: a
quiet check, a lock file held only while children are timed, and `brief.py`
refusing spawns while the lock is live.

## 7. Diagnosis: one profile first

The rule is in this PR; the tool is deferred (it saves typing, not hours). A
"why is this slow" brief gets one run of the case under cProfile, by hand:
`python -m cProfile -s tottime "$(command -v rasputin)" mesh <the case's
arguments from cases.toml>`, the top 15 calls by own time, and each as a share
of the run, a hotspot at 25 %. cProfile sees a `_core` or GEOS call as one
entry, which is what such a question needs; its overhead inflates
Python-heavy calls, so shares are reported, not times. The brief gets this,
its report, and stops. Candidate fixes, control catchments and sweeps are each
a new ask. A later PR may add `bench_quick.py profile <case>` (about 25 lines).

## 8. Size, seams, and what `@perf` confirmed

Counted as CLAUDE.md §2 counts: `bench.py` +20 (the `timeout` argument and
`Completed.timed_out` 10, phases in the child 10); `bench_quick.py` about 200
(models 30, cases, baseline and fingerprint I/O 35, the loop with the deadline
45, verdict 45, hotspots 20, Typer 25). About 220 in all.
`tools/bench_quick.py` joins mypy's `files`.

Seams for `tests/python/test_bench_quick.py`: `verdict(new, base)` and
`hotspots(record, ruled)` pure over records built in the test; the directory
fingerprint over a `tmp_path` tree (a renamed, resized or touched file changes
it); the deadline with a fake `Runner` whose child returns `timed_out`; one
real child (skipped without `_core`) to prove the `PhaseClock` replacement
still sees `mesh`'s clock (`src_python/tin_engine/cli.py@85a3e6dd:787`).

**Confirmed by `@perf`, 2026-10-07**, with these corrections folded in above:
(1) the per-run times, now in section 3's table; (2) the `_core` Release build,
18 s fresh and 0 s warm, not 1 to 4 min; (3) the 5 % / 0.05 s band holds for
whole runs but not for small phases, and the new run's spread could hide a
slowdown of up to about 35 % in one, so the band now uses the baseline's
spread only; (4) the catchments need no network. Also from `@perf`: capturing
phases by replacing `cli.PhaseClock` works (parsing the `--stats` table would
tie the tool to its layout, so it is not used); directories needed a
fingerprint rule (section 3); Ljungan's hotspot goes stale with 30d, so the
baseline is taken after 30d merges and `lagan` replaces `ljungan`.

## 9. Rule text this needs (proposed, not applied)

**`docs/increments/README.md`, *Acceptance*, a new first paragraph (+137 words):**
"**The quick check, when speed could change.** Each design says in one line
whether the change could change speed, and why. It could when it changes what
work `rasputin mesh` does: a step on by default, how files are read or written,
refine, decode or feature reading. A new flag off by default, a reworded error
or help text has no effect. On yes, the main session runs
`python tools/bench_quick.py run` before the push: fixed cases against the
stored master baseline, stopped at 15 minutes
(`docs/increments/perf-quick-check.md`). On no, nothing runs. It answers "is it
slower"; the run below never does. The run below also runs on master on the
first unattended night of each month, which catches what the one-line judgment
misses, and on Ola's ask. A `HOTSPOT` line goes to Ola as an `ASK OLA:` line."

**`.claude/agents/perf.md`, §1, two bullets (+44 words):**
"* **The quick check** (`tools/bench_quick.py`): 15 minutes at most; its cases,
baselines and ruled hotspots are in `docs/benchmarks/quick/`.
* **Hotspots go to Ola.** A `HOTSPOT` line, from the quick check or a profile,
is an `ASK OLA:` line in your handback, regression or not."

**`.claude/agents/perf.md`, §2, one bullet (+34 words):**
"* **"Why is it slow" gets one profile first**: one run of the case under
`python -m cProfile -s tottime`, and its report. Candidate fixes, control
catchments and sweeps are each a new ask."

**`CLAUDE.md` §3, *Step order*, appended (+34 words):**
"The main session runs the quick check itself, with no subagent running, and
spawns `@perf` only on a verdict other than `NO CHANGE`, on a `HOTSPOT`, for
the extensive run, or for a diagnosis."

Net: +249 words across three files; none removed. The `description:` line
of `perf.md`'s front matter stays as it is.

## 10. Questions for Ola (all ruled 2026-10-07)

Asked in the design's handback at `a7a15f1b`; recorded here as asked.

1. Should the main session run the quick check itself, and spawn `@perf` only
   when something is found? Default: yes.
2. Are 40 % of a run (for a phase) and 25 % (for a single call in a profile)
   the right hotspot thresholds? Default: yes.
3. Should every increment whose change touches code that `rasputin mesh` runs
   get the quick check before its push, not only refine and mesh changes?
   Default: yes.
   **Ola, 2026-10-07: "I'm not sure I support q3." "It's very dogmatic."
   "Aay you change the CLI - why would we want this?"** The main session then
   proposed the one-line judgment now in section 3 (*When it runs*), and
   **Ola: "Ok, this is good".** *Ruled: no blanket trigger; each design's
   one-line speed judgment decides.*
4. Should the full extensive run happen on master on the first unattended
   night of each month? Default: yes.

**Ola, 2026-10-07, on the main session's list "3. The measuring process: the main session runs the quick check itself; the hotspot thresholds are 40 % for a phase and 25 % for a single call; the full run happens once a month on an unattended night": "3: yes".** *Ruled: questions 1, 2 and 4 take their defaults. Recorded by the main session.*
