# Increment pq: the 15-minute performance check

**Status:** designed by `@architect` on `85a3e6dd`, 2026-10-07; to be read by
`@perf` before Ola sees it (section 8 lists what `@perf` should confirm).
Tooling only: `tools/bench.py`, a new `tools/bench_quick.py`, `tools/brief.py`.
One PR. No refine or mesh code, so no acceptance run of its own.

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

| case | inputs | threads | runs | est. per run |
|---|---|---|---|---|
| `tile` | bench.py's 1 m tile (`tests/fixtures/dem_archive/7908_3_10m_z33.tif`), `--tolerance 1` | default and 1 | 5 | ~0.9 s |
| `quarter` | same DEM, `docs/benchmarks/2026-09-26/quarter.geojson` | default | 5 | ~0.9 s |
| `numedalslagen` | DTM10, NVE outline, `corine2018_dtm10_utm33.gpkg`, `--tolerance 10` | default | 1 warm-up + 3 | ~7 s |
| `ljungan` | GLO-30 cache, SMHI outline, CORINE GeoJSON window, `--out-crs EPSG:3006`, `--tolerance 10` | default | 1 warm-up + 3 | ~9.3 s |

The two catchments are 30d's controls; `ljungan` has Lagan's staircase outline
at an eighth of its run time, so today's buffer hotspot shows on it (`decode`
4.94 of 9.29 s, 53 %). Every run is `--binary`. Paths are relative to
`$RASPUTIN_DATA` (default `../rasputin_data`); the Numedalslågen outline is
copied there from `rasputin_scratch`, which is a results folder.

**What is measured.** Per run: the process wall time, and each `--stats` phase
row with the clock's total. The child (`bench.py _child`) captures these by
replacing `cli.PhaseClock` with a subclass that keeps its instance, and adds
`phases` and `total_s` to its `BENCH` line (old records still parse). Also
`max_error` and the mesh's SHA-256 from the last run of each case.

**The baseline** is one file per power state, `docs/benchmarks/quick/baseline-ac.json`
and `baseline-battery.json`, a run of the check on master: commit, machine,
input fingerprint (path, size, mtime; content SHA-256 for files under 100 MB),
and per case and measure the median, min and max. A machine, power state or
input mismatch is `NO BASELINE: <field>`, never a cross comparison. It is
refreshed with `--save-baseline` (a) in the PR whose check showed a change that
was accepted (`FASTER`, or a `SLOWER` Ola accepted), measured at that branch's
head, (b) after a `NO BASELINE`, on a master checkout, and (c) by the monthly
extensive run.

**Regression.** Per case, for the total and for each phase at 5 % or more of the
baseline's total: `SLOWER` when the new median exceeds the baseline's by more
than the band *and* by more than 0.05 s. The band is the largest of 5 % (Ola's
2026-09-27 ruling for bench.py), the baseline's own (max − min) / median, and
the new run's. Symmetric for `FASTER`. A `max_error` over the tolerance is
`BROKEN`. A changed mesh hash is reported, not judged.
*Scale:* phases of 0.1 s and up, runs of 1 to 10 s. Checked against 30d's
master table (3 runs each, 7 to 103 s): the totals' spread is 0.5 to 4.4 %,
the phases' 0.5 to 15 % with one outlier at 37 % (Numedalslågen's
`features read`, 0.26 to 0.41 s), which that phase's band then absorbs;
bench.py's `refine_s` spread was 0.3 to 5.9 %.

**Report.** One line per finding, then the verdict:

```
SLOWER: ljungan decode 4.94 -> 6.10 s (+23.5 %, band 5.0 %)
HOTSPOT: ljungan decode 53 % of the run (not in hotspots.toml)
NO CHANGE | FASTER | SLOWER | BROKEN | NO BASELINE: <field> | OUT OF TIME: <cases not measured>
```

Exit 0 (no change, faster), 1 (slower or broken), 2 (no baseline), 3 (bad input
or not quiet, section 6), 4 (out of time). A hotspot does not change the exit.
Evidence goes to the scratchpad; it is committed under `docs/benchmarks/<date>/`
only with a finding or a baseline refresh.

**The 15 minutes is enforced in the tool.** The deadline is the command's start
plus `--budget` (default and maximum 900 s; more is refused). Every subprocess
(CMake and each child) gets `timeout = deadline − now − 10 s`, through a new
`timeout` argument of bench.py's `Runner.run`; on expiry the child is killed.
Before each case, if the baseline's median times its runs exceeds what is left,
the case is skipped rather than started. Either way the verdict is
`OUT OF TIME` naming the cases not measured, and what was measured is still
judged. Estimate: the Release build of `_core` 1 to 4 min (unmeasured), the
cases about 2 min.

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
profile (section 7). *Scale:* runs of 1 to 100 s; at 40 %, 30d's Lagan
(`decode` 49 %, `features read` 41 % of 72.7 s) and `ljungan` (53 %) are
flagged, Numedalslågen's largest phase (`features clip`, 26 % of 6.97 s) is not.

**One place:** `docs/benchmarks/quick/hotspots.toml`, one entry per ruled
hotspot: case, phase or call, share, date, Ola's ruling word for word. The tool
prints `HOTSPOT:` for one not in the file, or grown 10 points or more beyond its
ruled share. `@perf` (or the main session) turns each into an `ASK OLA:` line;
after Ola rules, the entry is added in the next commit that touches the file.
The first quick run may well raise the refine phases on `tile` (not
measured); one "expected" ruling then silences them.

## 6. The lock covers the timing only

`bench.py`'s `timing(label, deadline)` context, used by `run` and by the quick
check around the measured children (not the build, not the write-up):

1. **Quiet check.** It refuses (exit 3, naming the files) while any note file
   other than `session.md` and `perf-*.md` is in the main checkout's
   `.claude/current-task/`: a subagent is running.
2. **Lock.** It writes `<git-common-dir>/timing.lock` (pid, label, start,
   deadline) and removes it on exit. Live means the pid is alive and the
   deadline, if any, not passed; a dead lock is ignored with a warning.
3. **`tools/brief.py`** refuses any spawn while the lock is live, naming the
   label and start. Its rule "nothing runs beside `@perf`"
   (`tools/brief.py@85a3e6dd:162-163`) is dropped: `@perf` may be briefed
   beside others, and simply cannot time until they are done.

## 7. Diagnosis: one profile first

`python tools/bench_quick.py profile <case> [--budget 900]`: one run of the
case under `cProfile` (the child under `python -m cProfile`), the top 15 calls
by own time, and each as a share of the run, `HOTSPOT` at 25 %. cProfile sees
a `_core` or GEOS call as one entry, which is what such a question needs; its
overhead inflates Python-heavy calls, so shares are reported, not times. A
"why is this slow" brief gets this, its report, and stops. Candidate fixes,
control catchments and sweeps are each a new ask.

## 8. Size, seams, and what `@perf` should confirm

Counted as CLAUDE.md §2 counts: `bench.py` +45 (the `timeout` argument and
`Completed.timed_out` 10, phases in the child 10, `timing` 25);
`bench_quick.py` about 215 (models 30, cases and baseline I/O 25, the loop with
the deadline 45, verdict 45, hotspots 20, profile 25, Typer 25); `brief.py`
about +5. About 265 in all. `tools/bench_quick.py` joins mypy's `files`.

Seams for `tests/python/test_bench_quick.py`: `verdict(new, base)` and
`hotspots(record, ruled)` pure over records built in the test; the deadline
with a fake `Runner` whose child returns `timed_out`; `timing` against a
`tmp_path` git directory; one real child (skipped without `_core`) to prove the
`PhaseClock` replacement still sees `mesh`'s clock
(`src_python/tin_engine/cli.py@85a3e6dd:787`).

`@perf`, please confirm or correct: (1) the per-run times in section 3's table
(mine come from `docs/benchmarks/2026-10-04/23bm-*/run.json`, 210 samples in
3.1 min, and 30d's master table); (2) the Release build time of `_core` into a
fresh and a warm `build-bench`; (3) that the 5 % / 0.05 s band holds on 3 runs
of the catchments; (4) that `ljungan` reads only the GLO-30 cache, with no
network.

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

**`.claude/agents/perf.md`, §2, one bullet (+43 words):**
"* **"Why is it slow" gets one profile first**: `bench_quick.py profile <case>`,
and its report. Candidate fixes, control catchments and sweeps are each a new
ask. The lock is the timing: `bench.py` holds it while it measures, and your
write-up runs beside others."

**`CLAUDE.md` §3, *Step order*, appended (+30 words):**
"The main session runs the quick check itself and spawns `@perf` only on a
verdict other than `NO CHANGE`, on a `HOTSPOT`, for the extensive run, or for a
diagnosis."

Net: +254 words across three files; none removed. The `description:` line
of `perf.md`'s front matter stays as it is.

## 10. Questions for Ola (question 3 ruled 2026-10-07; 1, 2 and 4 open)

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
