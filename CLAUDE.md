# CLAUDE.md - Core Agent & Orchestration Guide

This file drives the automation loops for Claude Code in this repository. The
main session dispatches them (§3).

## 1. Personas
Each task goes to its persona, one file per persona in `.claude/agents/`; the
file's `description` says what the persona does.

## 2. Core Constraints & Technical Mandates
* **Strict Size Limit:** Under **700 net lines of production code per pull
  request**: lines added minus lines removed, both under the same exclusions. A
  line counts unless it is blank, a comment, a docstring, or the body of a raw
  literal. Files under `tests/` or `docs/`, and files that are not Python, C++,
  CMake or shell, are not counted; a file renamed across that boundary counts in
  full. Added and removed lines are the `+` and `-` ranges of the PR's
  `git diff -U0` hunks against its base, judged at the new and the old
  revision; `python3 tools/count_loc.py <base> <head>` counts them. Lines count
  as written: packing code by hand under `# fmt: skip` /
  `# fmt: off` is allowed, provided the packed lines stay readable and the
  review says why each new region is packed. This is the only statement of the
  rule; everywhere else points here.
* **Rule text pays for itself:** a change that adds words to the files
  `tools/rule_sizes.py` counts cuts as many in the same change, unless Ola
  explicitly allows the growth.
* **Prohibited Dependencies:** Never introduce `CGAL`, `GDAL`, `OGR`, `Fiona`,
  `Rasterio` (it wraps GDAL), or external `date` libraries.
  Enforced by `tools/check_prohibited_deps.py` over imports, includes, declared
  dependencies and build directives.
* **Core Stack:** Modern C++ (C++20 Concepts, Pybind11, `std::chrono`) + Async Python 3.12+ (Pydantic V2, Typer, Shapely, PyProj, NumPy, tifffile).
* **I/O Boundary:** File decoding is Python's. The C++ core never opens a file, sees a path, or links a codec, and CRS never crosses into it. See the `raster` section of `project_structure.md`.

## 3. Test-Driven Development (TDD) Protocol

Before acting on this repository — in the main session as well as in any persona
— read `.claude/REQUIRED-READING.md` and `docs/increments/README.md`. A session
that starts cold or resumes after a context loss runs the recovery steps at the
top of `.claude/REQUIRED-READING.md` **before** its first spawn, commit or edit;
that is the only statement of the rule and of where in-flight state lives.
Every code alteration follows the step order below, dispatched by the main session.

### The main session dispatches
The main session, started with no agent name, is the dispatcher: it spawns
the personas. Its rules:
* **Step order.** Design (`@architect`), failing tests (`@tester`), code
  (`@developer`), test run, review (`@reviewer`), and `@perf`'s acceptance run
  when refine or mesh code is touched (`docs/increments/README.md`).
  A performance fix is timed by `@perf` before review. Failing
  tests go back to `@developer`.
* **Resume while small.** Resume a persona only while its transcript is under
  about 150 k tokens; past that, spawn fresh on the increment file.
* **Worktrees.** Make each new worktree with `tools/new_worktree.py`.
* **Finish before starting.** At most three items are in flight: product,
  harness and process work, retrospectives included. An item counts from
  its first spawn until its PR merges or Ola parks or drops it, except while
  it waits only on Ola. A new item starts only when fewer than three count.
  If Ola asks for one beyond that, warn him with the count ("this makes
  4 of 3") and start it only after he confirms. A free writer goes first to
  the item nearest to merging, then to product work before harness and process work.
  The next item is chosen with Ola from his private ideas file.
* **No asking between internal steps, but stop at the remote.** Do not ask
  Ola between steps (tests to code, code to review); loop until `@tester` and
  `@reviewer` are satisfied. `.claude/REQUIRED-READING.md` says which acts
  need a fresh yes.
* **Conflicts before a push question.** Before asking Ola to push a branch,
  run `git fetch` and `git merge-tree --write-tree origin/master <branch>`,
  and the same against each branch queued or asked for ahead of it. When
  master conflicts, a writer merges it in first (`tools/merge_master.py`);
  when a branch ahead conflicts, the question names the files.
* **Prose edits.** Wordings and review records: `docs/increments/README.md`
  step 4. With both writer slots busy, the main session may merge
  append-only blocks on a prose branch; it composes no text in a persona's
  file and writes no code.
* **Recap each round.** At every new round or increment, run
  `python3 tools/session_state.py` and open with its recap.
* **Milestone updates.** After each milestone, one line to Ola on where the
  round stands.
* **Report only finished, verified results, briefly.** Answer Ola's question
  first, then stop. No file that does not exist yet, no number from a run
  still in progress, no cause not checked. No unasked images, no undefined
  jargon. Two or more questions or acts awaiting Ola are numbered in one
  sequence, each explained, with its default.
* **Lessons.** A persona reports a lesson in its handback; pass it to
  `@orchestrator`, which records it in `docs/retrospectives/`.
* **When to spawn `@orchestrator`.** After each increment merges (a check of
  that increment, and its lessons recorded); the morning after each
  unattended night (idle time, guard refusals and false positives, work done
  out of role); and for a research round, weekly or when Ola asks.
  `tools/duties.py` shows which are owed.

## 4. Operational Commands

### Python Layer (Async & CLI)
```bash
pytest                 # Run all modern tests (uses pytest-asyncio)
pytest tests/python/   # Run target Python testing suite
```

### C++ Core Layer (C++20 Concepts)
```bash
cmake -S . -B build && cmake --build build -j   # -j alone: nproc is Linux-only
ctest --test-dir build                          # all registered suites; see tests/cpp/CMakeLists.txt
```
The compiler gate: the C++ targets build with `-Wall -Wextra -Wpedantic -Werror`
(`CMakeLists.txt`, `tests/cpp/CMakeLists.txt`), so any warning fails the build.

### Static gates (Python)
```bash
mypy                   # strict, over src_python/tin_engine
ruff check .
ruff format --check .  # docs/ and .claude/ are excluded
```

### Governance gates
```bash
python tools/check_prohibited_deps.py   # section 2, checked against real imports
python tools/check_detria_boundary.py   # detria.hpp stays in one TU, zero headers
python3 tools/check_citations.py        # cited lines resolve; lists the at-risk ones
python3 tools/check_code_history.py     # no history labels in Python docstrings or comments
```

### CI
The gates above also run in GitHub Actions (`.github/workflows/main.yaml`), and
CI is authoritative: local green does not mean the branch is green, so
**verify check status before declaring anything merge-ready**:
```bash
gh pr checks <pr>      # must be green; required for merge on master
```
`CI result` is green only when every gating job passed, or was skipped on a
prose-only change. A PR that changes
only prose (Markdown at the root or under `docs/`, less the files the suites
read; `tools/ci_changes.py` decides) skips the C++, sanitizer and Python jobs
and runs only the governance gates, so its green says nothing about code.
`master` merges through a merge queue: `gh pr merge <pr>` enqueues the PR (the
queue's method, merge commit, applies), and the queue tests it on top of the
PRs ahead and merges it. No update-branch loop. A method flag such as `--merge`
only prints a warning (the queue's method wins); `-d`/`--delete-branch` is
refused. Never pass `--admin`: it bypasses the queue.

