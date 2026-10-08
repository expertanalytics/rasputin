# Increment pq: the 15-minute performance check

**Status:** designed by `@architect` on `85a3e6dd`, 2026-10-07; read by
`@perf` the same night and corrected (section 8); Ola's rulings in section 10.
Red step `fa3e03ef` (`@tester`), green step `b3439629` (`@developer`), the two
data files and the as-built update in `807130ec` (`@architect`). Section 11
rules on the green step's assumptions; its one change, the file fingerprint,
is done: red `166f2528`, green `db1ce78e`. Code review round 1 asked for
changes (the Review section); master merged in as `8c05d65e`. 242 net
production lines (`python3 tools/count_loc.py 483221d2 8c05d65e`) against the
estimate of about 220. Code review round 2 approved; merged as #217
(`64481aae`). Follow-up (section 12): record the shapely, GEOS, pyproj and
PROJ versions, designed 2026-10-08; next `@tester`'s red step.
Tooling only: `tools/bench.py` and a new `tools/bench_quick.py`. One PR. No
refine or mesh code, so no acceptance run of its own. Its first baseline is
taken on master after this PR merges (30d, the outline-buffer fix, merged as
#215). Deferred to a later PR: the timing lock (section 6) and the `profile`
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
`--binary`. Paths start with `$TREE` (the tree measured) or `$RASPUTIN_DATA`
(default `../rasputin_data` beside the checkout the tool runs from; set it when
running from a worktree). The Numedalslågen outline was copied there from
`rasputin_scratch/norway/numedalslagen/numedalslagen_outline_nve.geojson`
(a results folder) on 2026-10-07; `cmp` of the two files finds no difference.

**What is measured.** Per run: the process wall time, and each `--stats` phase
row with the clock's total. The child (`bench.py _child`) captures these by
replacing `cli.PhaseClock` with a subclass that keeps its instance, and adds
`phases` and `total_s` to its `BENCH` line (old records still parse). Also
`max_error` and the mesh's SHA-256 from the last run of each case.

**The baseline** is one file per power state, `docs/benchmarks/quick/baseline-ac.json`
and `baseline-battery.json`, a run of the check on master: commit, machine,
input fingerprint, keyed by the input as written in `cases.toml` (for a file
under 100 MB: size and content SHA-256; at 100 MB or more: size and mtime; for
a directory, such as DTM10's 762 files or the GLO-30 cache: the
SHA-256 of its sorted relative file names, sizes and mtimes, recursive, no
content read),
and per case and measure the median, min and max. A machine, power state or
input mismatch is `NO BASELINE: <field>`, never a cross comparison. The
baselines, `cases.toml` and `hotspots.toml` are read from, and
`--save-baseline` writes to, the `docs/benchmarks/quick` of the checkout the
tool runs from, even with `--tree` naming another. It is
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
plus `--budget` (default and maximum 900 s; more is refused). CMake and each
child get `timeout = deadline − now − 10 s`, through a new
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
45, verdict 45, hotspots 20, Typer 25). About 220 in all; as built, 242
(`bench.py` +25, `bench_quick.py` +217).
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

## 11. As built: the green step's assumptions, ruled

`@architect`, on `b3439629`. The data files: `docs/benchmarks/quick/cases.toml`
is `@developer`'s draft unchanged. It agrees with section 3's table and with
the arguments of 30d's runs (`docs/benchmarks/2026-10-07/30d-buffer/scripts/run_mesh.sh`,
cases `gpkg33` and `gpkg`), and every input it names exists.
`docs/benchmarks/quick/hotspots.toml` holds Ola's first ruling: `lagan`'s
`features clip`, 40.7 %, recorded as expected for now (2026-10-07). The phase
name is the `--stats` row `mesh` writes (`src_python/tin_engine/cli.py@b3439629:1466`).

1. **`$TREE` and `$RASPUTIN_DATA` in `cases.toml`**, expanded when the tool runs.
   *Kept.* The default data folder is `rasputin_data` beside the checkout the
   tool runs from, which is wrong in a worktree under `.claude/worktrees/`.
   The run then stops with exit 3, naming the missing path. From a worktree,
   set `RASPUTIN_DATA`.
2. **Fingerprints keyed by the input as written, not the expanded path.**
   *Kept*: a master baseline and a branch worktree expand `$TREE` differently.
   **Changed alongside it:** a file under 100 MB is fingerprinted by size and
   content only, not mtime. Git sets a file's mtime when it checks the file
   out, so the same fixture has a different mtime in each checkout. With mtime
   in the fingerprint, a check run with `--tree` on another checkout would
   always print `NO BASELINE: input $TREE/... changed`. (Checked on this
   machine: `7908_3_10m_z33.tif` and `quarter.geojson` have the same size and
   SHA-256 in the main checkout and in this worktree, but different mtimes.)
   A file of 100 MB or more, and a directory, keep size and mtime, because
   their content is not read. *For `@tester`:* a file under 100 MB that is
   touched keeps its fingerprint (this replaces
   `test_a_file_fingerprint_sees_a_touch`); a file of 100 MB or more that is
   touched gets a new one.
3. **Machine mismatch compares every field of `bench.Machine`**, including the
   macOS, Python and numpy versions; `bench.py` compares only the CPU and the
   core counts. *Kept*: the quick check times whole runs, and decode and
   feature reading run in Python. A version change gives `NO BASELINE`, and
   rule (b) of section 3 takes a new master baseline in about 2.5 minutes. A
   known gap: the shapely, pyproj and GEOS versions are not recorded, and
   GEOS does `features clip`'s work. Closed by the follow-up in section 12.
4. **The verdict:** BROKEN over SLOWER over FASTER over NO CHANGE; `OUT OF
   TIME` replaces any of them, with exit 4, and the lines of what was measured
   are still printed. *Kept*, as section 3 says. A known gap, from code review
   round 1: `NO BASELINE` is returned before `max_error` is checked
   (`tools/bench_quick.py@db1ce78e:121-123`), so a run with no matching
   baseline, such as the first `--save-baseline`, cannot report `BROKEN`. Left
   as is, because the fix is code and this round is prose; a later PR moves
   the check ahead. `bench.py`'s acceptance still judges its own cases'
   tolerance (`tools/bench.py@db1ce78e:352`).
5. **Two extra lines,** `MESH CHANGED` (the mesh hash moved, reported, not
   judged) and `NOT JUDGED` (a case missing from the baseline). *Kept.*
6. **`--save-baseline` judges first, then writes**; it writes nothing after
   `OUT OF TIME` or when the power state changed during the run, and the first
   save, with no baseline yet, exits 2. *Kept.*
7. **The skip rule** uses the baseline's median total times the case's runs,
   without its warm-up or start-up. *Kept*, as section 3 says. The estimate
   can be short by the warm-up; a case that then overruns is killed at the
   deadline and named in `OUT OF TIME`, the same verdict a skip gives.
8. **A ruled hotspot is raised again at 10 points or more** above its ruled
   share. *Kept*, as section 5 says.
9. **Each run's mesh goes to a `tempfile.mkdtemp` folder that is not
   removed.** *Kept.* It holds one mesh, overwritten case by case (about the
   size of Lagan's binary mesh), in the per-user temporary folder.

## 12. Follow-up: the GEOS version

A second PR, after #217 merged pq (`64481aae`). Designed by `@architect`,
2026-10-08.

**The ruling.** The main session asked, in its list after pq's code review
round 2 (as given to `@architect`, with its elisions): "3. Record the GEOS
version in the quick check ... The fix: about 2 lines, plus a test.
@architect's default: no, not in this PR ... My recommendation: yes, because
the upgrade lands this week. Your call." **Ola, 2026-10-08: "On 1, try
increasing to 1m as well. I'm pretty sure 5cm is still only noise. Yes to
rest."** ("On 1" answers another question of that list, not this
increment's.) *Ruled: yes, as its own PR.*

**What is recorded.** Four versions, added to `bench.py`'s `Machine`:
`shapely` (`shapely.__version__`), `geos` (`shapely.geos_version_string`, the
GEOS library shapely runs on), `pyproj` (`pyproj.__version__`) and `proj`
(`pyproj.proj_version_str`). GEOS does `features clip`'s work and the buffer
of the outline; PROJ does the reprojection that `lagan` runs
(`--out-crs EPSG:3006`, the continental GeoPackage reprojected), so a PROJ
upgrade can move `lagan`'s `features read` as a GEOS upgrade can move its
`features clip`. The two Python package versions are kept beside them
because the wheels bundle the libraries: an upgrade normally arrives as a
new wheel, and the package version is what `pip` and the lock file show.
*Object identity:* `_machine` runs in the parent, and the children run with
the same `sys.executable` (`tools/bench_quick.py@64481aae:227`); `--pkg`
redirects only `tin_engine`, so the parent's shapely and pyproj are the ones
the children import. (Checked on this machine's main venv: shapely 2.1.2,
GEOS 3.13.1, pyproj 3.8.0, PROJ 9.8.1.)

**Where.** The four fields are `str | None = None` on `Machine`. The default
is needed: `bench.py` loads committed `run.json` records under
`docs/benchmarks/` (for example `docs/benchmarks/2026-09-27/21a/run.json`),
which lack them. `bench.py`'s own comparison key stays the CPU and the core
counts, so its acceptance is unchanged. `bench_quick.py` changes not at all:
it already compares every `Machine` field (`tools/bench_quick.py@64481aae:110-113`),
so a version that differs from the baseline's gives
`NO BASELINE: machine geos: <new> vs <old>`, exit 2, as numpy does today.
`bench.py`'s report names the four in its *Method* line beside numpy.

*Consequence, for the upgrade itself:* the quick check does not judge a run
across a GEOS or PROJ upgrade; it says `NO BASELINE`, and rule (b) of
section 3 takes a new master baseline. To measure what the upgrade did to
speed, run the quick check on master with `--save-baseline` before
upgrading, upgrade, and compare the two baseline files' medians by hand, or
run `bench.py`'s acceptance, whose key ignores the versions.

**The one red test** (`@tester`, `tests/python/test_bench_quick.py`): build
a `Machine` through `bench._machine` with a fake `Runner` for `sysctl` and
`sw_vers`; assert its `geos` and `proj` equal `shapely.geos_version_string`
and `pyproj.proj_version_str`; then judge a record carrying that machine
against a baseline whose machine differs only in `geos`, and assert
`NO BASELINE` naming `geos`, exit 2. Red today: `Machine` has no `geos`
field. Existing tests that build `Machine` from a dict keep passing, since
the new fields default to `None`.

**Speed.** No effect on speed: it adds four version strings to the run
record and changes nothing `rasputin mesh` does, so the quick check does not
run for this PR.

**Size**, counted as CLAUDE.md section 2 counts, all in `tools/bench.py`:
the four fields 4, the two imports 2, two more arguments in `_machine`'s
packed `Machine(...)` call (already under `# fmt: skip`) 2, the *Method*
line 1. About 9 net production lines.

## Review

pq code review round 1, 2026-10-08, @reviewer (`85a3e6dd..db1ce78e`, 242 net production lines against about 220): CHANGES REQUESTED. (1) The pq row in /Users/skavhaug/projects/rasputin/ROADMAP.md@db1ce78e:72 is false as built: it still says Ljungan, "the timing lock covers only the timing", "about 265" lines and "designed; @perf to read". The as-built count in /Users/skavhaug/projects/rasputin/docs/increments/perf-quick-check.md@db1ce78e:202-203 is also stale (241 and +216; it is now 242 and +217), and so is the Status at lines 1-10, which still sends the fingerprint change back to @tester. (2) /Users/skavhaug/projects/rasputin/pyproject.toml@db1ce78e:137-138 conflicts with master `483221d2` (`git merge-tree`). Merge master in and keep both sides' mypy `files` entries. Code, the 194 tests, mypy, ruff and the governance gates are green; CI has not run (no PR yet). The four `# fmt: skip` regions (three new, three extended) each keep one call on two lines instead of five to eight; accepted for compactness, not size. Suggestions: `verdict` returns NO BASELINE before checking `max_error`, so a first `--save-baseline` cannot report BROKEN; section 3's "every subprocess" gets the deadline should say "CMake and each child"; `--save-baseline` writes to this checkout's `docs/benchmarks/quick` even with `--tree`, which should be stated.

pq code review round 2, 2026-10-08, @reviewer (`db1ce78e..69f18499`, 242 net production lines against about 220, unchanged from round 1): APPROVED. (1) Fixed. /Users/skavhaug/projects/rasputin/ROADMAP.md@69f18499:74 now matches the build: it names Lagan, says the timing lock and the `profile` subcommand are deferred, gives 242 lines and says "built, not pushed". The Status at /Users/skavhaug/projects/rasputin/docs/increments/perf-quick-check.md@69f18499:1-13 and the count at lines 206-207 (242: bench.py +25, bench_quick.py +217) are also correct, and both counts were re-run. (2) Fixed. Merge `8c05d65e` keeps both sides' mypy `files` entries at /Users/skavhaug/projects/rasputin/pyproject.toml@69f18499:137-139. `git diff 483221d2 69f18499 --stat` shows the same nine files as `85a3e6dd..db1ce78e`, and the tools, tests and `docs/benchmarks/quick` are unchanged since `db1ce78e`. `git merge-tree --write-tree origin/master HEAD` is clean. The two bench suites pass (194 tests), and mypy, ruff check, ruff format, check_citations and check_prohibited_deps are all green. CI has not run (no PR). The three round 1 suggestions are taken up in the prose and checked against the code. The NO BASELINE known gap in section 11 item 4 cites /Users/skavhaug/projects/rasputin/tools/bench_quick.py@db1ce78e:121-123 and /Users/skavhaug/projects/rasputin/tools/bench.py@db1ce78e:352, and both citations point at the right lines.
