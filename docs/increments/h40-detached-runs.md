# h40 — Detached runs: heavy model runs survive the agent

A heavy perf or land-surface-model (lsm) run goes out as a **detached job** —
a process in its own session, writing artifacts to disk, that lives
independently of the agent that launched it. The persona launches it and later
**re-attaches** by reading the finished artifacts. No agent ever holds a live
heavy run on its monitored stream.

Pipe 3.2 (`detached-runs`). Fast lane by the pipe, but this is a process-shaping
change touching how `@perf` works, so it is written as a full increment file
first (design only; no production code and no push in this step). The pipe calls
it "150 lines"; the estimate below is smaller because the on-disk artifact
format already exists in `bench.py`.

Note on numbering: the brief and the ROADMAP/pipe name three different ids for
this ("pipe 3.13", "3.2", "3.3 reap"). The live pipe (`work_pipes.py show`) is
authoritative: **this is 3.2 detached-runs; reap is 3.3; timed-stages is 3.13**.
The brief's "3.13" is a typo for the item it then describes correctly.

---

## Prior art: legacy and literature

**Legacy.** `git grep -n -iE 'detach|re-?attach|checkpoint' legacy-archive -- legacy`
returns nothing: the legacy archive is the old CGAL meshing code, which had no
agent harness and no long detached runs. Nothing is carried across.

**In-tree prior art (the real base).** Three tools already hold the pieces, and
this increment composes them rather than inventing a new store:

* `tools/bench.py` already writes a **run directory** per run:
  `docs/benchmarks/<date>/<label>/` with `run.json` (a Pydantic record,
  `RunRecord`), `raw.tsv`, per-sample raw files (`<sample>.json`,
  `<sample>.code`, `<sample>.stderr.txt`) and a `README.md`. A run that
  finishes writes `run.json`; a run that did not never gets one
  (`bench.py:521`, `write_evidence`). That "the final JSON is the done-marker"
  convention is reused verbatim.
* `tools/reap.py` lists and (`--kill`) stops every process whose working folder
  is inside a removed or merged worktree, or that has run over 6 h
  (`reap.py:57`, `select`). It keys on `cwd`, so a detached job started inside a
  worktree is already visible to it — **provided the job's cwd is inside that
  worktree** (see §5, boundary with reap).
* `tools/session_state.py` (the recap) already prints "Running Bash-tool jobs"
  (`background_jobs`, keyed on the Claude-Code shell wrapper) and "Stray
  processes" (via `reap.listing(reap.found())`). A detached job is **not** a
  Bash-tool wrapper, so it is invisible to `background_jobs` today; the recap
  gains one small reader of the detached-run manifests (§5).

**Literature.** This is standard Unix job detachment (`setsid(2)` /
`start_new_session`, a backgrounded process in its own session group, immune to
the launcher's SIGHUP and to a controlling-terminal / stream death) plus the
long-known "write to a temp name, atomic rename on completion" pattern for the
done-marker (POSIX `rename(2)` atomicity within a filesystem). No novelty claim;
nothing here is intended for publication. The one thing we must not assume is
cross-filesystem atomicity — the temp and final names must share a directory
(§2).

**Ideas file.** `../rasputin_scratch/ideas/ideas.md` and
`storage-state-facts-2026-10-20.md` were read. No existing idea entry covers
detached runs; reap is already a ROADMAP/pipe item (3.3) and is **not** designed
here. This increment overlaps reap at the cwd/liveness boundary only (§5).

---

## The problem, probed

2026-11-02: a `@perf` model run launched as a subagent was lost when the laptop
suspended ~20 min mid-run. Wall-clock elapsed during sleep, so Anthropic's
harness stream watchdog saw >600 s of no stream progress and killed the
subagent. The watchdog is not ours to change, and Ola ruled out raising its
timeout ("I believe fixing problems, not band-aiding them!").

The mechanism (not blame): the heavy run's liveness was coupled to the agent's
**monitored stream**. Any event that silences the stream — suspend, agent kill,
context loss, restart — destroys the run. Ola's framing: this is a "bleeding
effect across agents" — a live run held in one agent's context is lost or
corrupted when that agent is killed, restarted or suspended, and the next agent
inherits broken live state instead of a clean artifact.

Probe that the fix works (scratch, 2026-11-02,
`scratchpad/probe_detach.py`): a child launched with
`subprocess.Popen(..., start_new_session=True)` lands in its **own session and
process group** (`os.getpgid(child) == child != launcher pgid`), kept writing
its log, and finished on its own **after the launcher exited** without being
waited on. So a detached job's liveness does not depend on the launcher's
process, group, controlling terminal or stream. A suspend pauses and resumes the
job's own wall clock with no watchdog watching it; nothing to trip.

---

## The invariant

> **A heavy run never lives on an agent's monitored stream.** A run that is
> *heavy* — a `@perf` acceptance/profile run, or an lsm model run, expected to
> exceed a few minutes of compute — is launched as a **detached job** and
> re-attached from disk. The launching persona hands back as soon as the job is
> launched and its manifest is written; it does not hold the live run.

What is *heavy* (the rule that decides detach):
* any `@perf` timed acceptance run or scaling sweep (`bench.py run`);
* any lsm model run (per-frame / per-triangle field output), e.g. the north-face
  melt instrumented run (q146);
* any run a persona's own judgement expects to exceed **~3 minutes** of compute.

What stays on-stream (not detached): quick checks, probes, builds, the test
suite, `bench_quick.py`, anything a `Monitor`/`run_in_background` Bash job
already covers and finishes well inside the watchdog window. A short Bash
`run_in_background` job is still a monitored harness job and is fine; the
invariant is about *heavy* runs, not about all background work.

This is one direction of Ola's "reduce bleeding effects": the only thing that
crosses between agents is a clean on-disk artifact, read by the next agent,
never inherited live process state.

---

## The mechanism

A new standard-library tool, **`tools/detach.py`**, with two verbs.

### `detach.py start`

```
python3 tools/detach.py start --run-dir <DIR> --name <NAME> -- <cmd> [args...]
```

* Creates `<DIR>` if absent. `<DIR>` is a run directory under
  `docs/benchmarks/<date>/<label>/` (perf) or
  `docs/benchmarks/<date>/<label>/` likewise for an lsm run — the same tree
  `@perf` already commits, so detachment adds no new storage location.
* Launches `<cmd>` with `subprocess.Popen(cmd, start_new_session=True,
  stdout=<DIR>/<NAME>.stdout, stderr=<DIR>/<NAME>.stderr, stdin=DEVNULL,
  cwd=<launcher's repo/worktree root>)`. `start_new_session=True` is the whole
  detachment (probed above). `stdin=DEVNULL` so the job can never block on input
  (the 2026-10-30 reap case, `reap.py:17`).
* Writes a **manifest** `<DIR>/<NAME>.detached.json` atomically (write to
  `<NAME>.detached.json.tmp` in the *same* directory, then `os.replace`):

  ```json
  {
    "name": "<NAME>",
    "pid": 20660,
    "pgid": 20660,
    "cmd": ["python3", "tools/bench.py", "run", "..."],
    "cwd": "/abs/worktree/root",
    "started": "2026-11-02T10:46:00Z",
    "host": "<platform.node()>",
    "stdout": "<NAME>.stdout",
    "stderr": "<NAME>.stderr",
    "done_marker": "<NAME>.done.json",
    "state": "running"
  }
  ```

* Prints the manifest path and exits **0 immediately** — it never waits on the
  job. The persona hands back here.

The launched `<cmd>` owns its own progress and completion artifacts (the detach
tool does not interpose on the job's stdout stream). The **completion contract**
is: the job writes its final artifact *last*, atomically, under a name the
manifest's `done_marker` field names. For a `bench.py run`, `run.json` is the
natural done-marker (it is written only on success, `bench.py:521`); the
manifest's `done_marker` is set to `run.json`. For an lsm run, the run writes
per-frame fields as it goes (partial, see below) and a final `summary.json`
last; `done_marker` is `summary.json`.

### Partial vs complete — the on-disk contract

* **Complete:** the `done_marker` file exists and parses. A re-spawned persona
  that sees it reads the finished result and does **not** redo the run.
* **Partial:** the manifest exists, the `done_marker` does not. The persona
  checks liveness (§3) to tell *running* from *crashed*.
* **Checkpoint / per-frame:** a long lsm run writes each frame to its own file
  (`frame_00001.npy`, …) as it computes it, and never rewrites an earlier frame.
  So a partial run's completed frames are intact and a resumed run can skip
  them. This is the job's own responsibility, not `detach.py`'s; `detach.py`
  only owns the manifest and the launch. The q146 run — a single-config lsm run
  writing per-frame fields — fits this exactly: one `start`, per-frame files,
  final `summary.json`.
* **Atomicity:** every marker and every frame file is written to a `.tmp` name
  in the same directory and `os.replace`d, so a kill mid-write never leaves a
  half-file a reader will accept. (Shared-directory requirement from the
  literature note.)

### `detach.py status`

```
python3 tools/detach.py status --run-dir <DIR> [--name <NAME>]
```

Reads each `*.detached.json` manifest under `<DIR>` and prints, per job, one of:

* `done <NAME> <DIR>/<done_marker>` — done-marker present and parses;
* `running <NAME> pid <PID>` — no done-marker, and `os.kill(pid, 0)` on the
  manifest's `pid` succeeds **and** the running pid's start matches (host guard,
  see below);
* `failed <NAME> (no done-marker, process gone)` — no done-marker and the pid is
  gone (crashed / killed / suspended-host-powered-off);
* `unknown <NAME> (host <H>, not this machine)` — manifest `host` ≠ this host, so
  liveness can't be judged here.

`status` is a **pure read of disk + one `kill(pid,0)` probe**; it starts no
heavy work and is safe to run from any persona or the main session.

A pid-reuse guard: on Linux compare `/proc/<pid>` start time; on macOS compare
`ps -o lstart= -p <pid>` against `started`. A mismatch is reported `failed`
(the original job is gone; some unrelated process now holds the pid). This keeps
`running` from ever being a false positive.

---

## The re-attach contract

How a persona or the main session discovers and reads a detached run:

1. The run's `--run-dir` is deterministic from the increment/task (perf:
   `docs/benchmarks/<date>/<label>/`), so a re-spawned persona knows where to
   look from its brief — it does not need inherited live state.
2. `python3 tools/detach.py status --run-dir <DIR>` tells it running / done /
   failed / unknown.
3. On **done**, the persona reads the done-marker and the raw artifacts exactly
   as `@perf` reads a finished `bench.py` run today (`run.json`, `raw.tsv`,
   per-sample files). Nothing new to read.
4. On **running**, the persona hands back and is re-spawned later; it never
   waits in-context (that would re-couple liveness to its stream).
5. On **failed**, the persona reports the failure (stderr file named) and, for a
   checkpointed lsm run, may relaunch with the completed frames left in place.

The main session's recap surfaces the same four states so a cold start sees
them (§5).

---

## Boundary with existing tools (flag overlaps; do not design them)

* **Recap / `session_state.py` (not redesigned here).** The recap lists
  Bash-tool wrapper jobs and stray processes but **not** detached jobs, because a
  detached job is neither a Claude-Code shell wrapper nor (while its worktree is
  live and it is under 6 h) a stray. Overlap: this increment adds **one small
  reader** to the recap — "Detached runs:" printing
  `detach.status()` over the known run-dir roots (default
  `docs/benchmarks/<newest date>/`), in the same costs-one-line-on-fault style
  the recap already uses for reap and pipes (`session_state.py:343`). That is
  the only `session_state.py` change and it is read-only.
* **reap (3.3, not designed here).** reap keys on `cwd` and on the worktree being
  removed/merged or the age > 6 h (`reap.py:57`). A detached job launched with
  `cwd = <worktree root>` is therefore already reapable. **Overlap / hazard:** a
  detached job is *meant* to outlive the agent, so reap must not kill a
  still-wanted detached job just because the launching agent is gone. The
  mechanism keeps them disjoint: reap kills on *worktree removed/merged* or
  *over 6 h*, never on "agent gone". A detached run's worktree should not be
  removed while its done-marker is absent (next bullet). reap's existing rules
  are sufficient; **no change to reap is designed here** — flagged for the reap
  builder: its listing should show the detached manifest's name when a reaped
  pid matches a manifest, so a human killing stray processes sees it was a
  detached run.
* **Worktree removal.** A detached job must not be orphaned when its worktree is
  cleaned. **Overlap / hazard:** `detach.py start` writes the run dir under
  `docs/benchmarks/`, which lives in the worktree. Two safe options, decided by
  Ola (open question Q1): (a) run dirs go under the **main checkout's**
  `docs/benchmarks/` (found via `--path-format=absolute --git-common-dir`, as
  `session_state.main_checkout` already does), so removing a worktree never
  touches a run; or (b) worktree removal refuses while any manifest under it
  shows `running`/partial (a `detach.status` check in the removal tool). The
  design does **not** build worktree-removal changes; it names the hazard and the
  default (Q1).

---

## How `@perf`'s file and the rule docs change (rule-text budget)

The pattern lives in **`tools/detach.py`'s docstring** (not counted by
`rule_sizes.py`, which counts `.claude/agents/*.md`, `.claude/skills/**/*.md` and
the core rule files — `rule_sizes.py:22-40`). So the governed rule text pays for
itself:

* `.claude/agents/perf.md` §2 "How a run is made" today says "Alone: timed
  stages run first…". Replace the first sentence of that bullet with: *"Heavy:
  a timed acceptance, scaling or lsm run goes out detached
  (`python3 tools/detach.py start`), and you hand back; re-attach with
  `detach.py status` and read the run dir."* This **removes** the now-redundant
  "then hand back" wording and the implied in-context wait, so the net word
  change in `perf.md` is near zero (add the detach sentence, cut the stream-wait
  words it replaces). The reviewer confirms the count with
  `python3 tools/rule_sizes.py`.
* No change to `CLAUDE.md`, `REQUIRED-READING.md`, `README.md` or
  `common.md` — the invariant is a `@perf`-and-lsm rule, enforced by the tool
  and the agent file, not a project-wide one. If the reviewer finds the
  invariant needs a home outside `perf.md` (e.g. an lsm persona that does not
  yet exist), that is a word add that must be paid for by a cut elsewhere, and
  it is flagged, not assumed.

---

## Testability — without an actual 20-min suspend

The whole point is that liveness no longer depends on the stream, so the tests
plant the *events* the stream-coupling used to lose, cheaply:

* **T1 detachment.** `detach.py start` on a `sleep` command: assert the child's
  pgid ≠ the test's pgid and the manifest `pid`/`pgid` match the child. Fixture:
  a 2 s `/bin/sh -c 'sleep 2; ...'`; measured in the probe above, the launcher
  returns immediately (<50 ms).
* **T2 survives launcher exit.** Start a job that writes a marker after a short
  sleep; let the Python launcher process **exit** (the stand-in for the agent
  dying / being killed); assert the marker appears. This reproduces the lost-run
  condition (agent gone, job alive) with no suspend. Measured in the probe: the
  detached child finished and wrote `DONE` after the launcher exited.
* **T3 suspend stand-in.** SIGSTOP the job, assert `status` still reports
  `running` (pid alive), SIGCONT it, assert it completes and `status` flips to
  `done`. SIGSTOP is exactly what a laptop suspend does to the process; this is
  the 20-min suspend in milliseconds.
* **T4 partial vs complete.** With no done-marker present, `status` reports
  `running`/`failed` by pid liveness; with the done-marker written, `done`. Kill
  a job mid-run (no done-marker, pid gone) → `failed`.
* **T5 atomic write.** Write a frame/marker through the `.tmp`+`os.replace` path,
  interrupt before replace (leave the `.tmp`), assert `status` does **not** read
  the `.tmp` as done.
* **T6 pid-reuse guard.** Point a manifest at this test's own pid (a process that
  is not the original job) with a mismatched `started`; assert `status` reports
  `failed`, not `running`.
* **T7 recap reader.** Give `session_state`'s new detached reader a run dir with
  one done and one failed manifest; assert the two lines, and that a malformed
  manifest costs one line, not a crash (the recap's standing fault policy).

Tests: `tests/python/test_detach.py` (T1–T6) and an addition to
`tests/python/test_session_state.py` (T7). Standard library only; no C++, no
`@perf` acceptance (this touches neither refine nor mesh code).

---

## Size, exclusions, degeneracy

* **Estimate.** `tools/detach.py` ~90 net production lines (two verbs, manifest
  I/O, pid-liveness with the host/start guard). The `session_state.py` recap
  reader ~12 lines. The `perf.md` edit is rule text, budget-neutral. Total under
  **~110 net production lines**, well inside 700.
* **Exclusions.** Not built here: reap changes, worktree-removal changes, any
  resume/replay logic for an lsm run (the job owns its own checkpoint skip —
  `detach.py` only distinguishes partial from complete). No daemon, no
  scheduler, no cross-host tracking beyond the `host` field reporting `unknown`.
* **Degeneracy.** Missing run dir → `status` prints nothing (or "none"), exit 0.
  Manifest present, `done_marker` and process both gone → `failed`. Manifest on
  another host → `unknown`. Corrupt/`.tmp` done-marker → treated as not-done.
  Two jobs same name in one dir → second `start` refuses (manifest already
  exists and shows running); overwriting is an explicit `--force`.

---

## Open questions for Ola (each with a default)

* **Q1 — where run dirs live vs worktree removal.** Default: run dirs go under
  the **main checkout's** `docs/benchmarks/<date>/<label>/` (option a), so a
  worktree can be removed without orphaning or losing a detached run, and the
  next agent reads the artifact from a stable path. Recommended over blocking
  worktree removal, which couples cleanup to job state. (Affects one line in
  `detach.py start`.)
* **Q2 — the heavy threshold.** Default: "~3 minutes of compute, or any
  `@perf` timed run or lsm run" as the rule for what must detach, left to the
  persona's judgement rather than a hard gate. Recommended: judgement, since a
  hard timer would itself need to watch the stream. Say if you want a fixed
  number instead.
* **Q3 — does the invariant bind beyond `@perf` and lsm now.** Default: it binds
  `@perf` and lsm runs only (the files that have them), and widens when an lsm
  persona exists. Recommended, to keep the rule-text budget flat today.
