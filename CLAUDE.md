# CLAUDE.md - Core Agent & Orchestration Guide

This file drives the automation loops for Claude Code in this repository. The
main session dispatches them (§3).

## 1. Multi-Agent Ecosystem
This project is governed by specialized sub-agents. Always defer tasks to the correct persona in `.claude/agents/`:
* `@orchestrator`: Workflow watcher. Flags deviations from these rules, researches agentic design and proposes improvements, and owns `docs/retrospectives/`.
* `@architect`: Enforces declarative structures, component boundaries, and interface decoupling. Checks the literature and any novelty claim before a design.
* `@tester`: Owns the test suites. Enforces the coverage floor of `testing.md` and adversarial geometry fuzzing.
* `@developer`: Writes clean, high-performance C++20 and async Python code.
* `@reviewer`: Final gatekeeper. Audits CI status, LOC, red-step scaffolding, and prose claims against code.
* `@perf`: Performance owner. Runs the benchmark and scaling acceptance for refine- and mesh-touching increments (`tools/bench.py`), profiles the serial phase, and keeps the evidence in `docs/benchmarks/<date>/`.

## 2. Core Constraints & Technical Mandates
* **Strict Size Limit:** Under **700 lines of production code per pull request**,
  where a line counts unless it is blank, a comment, a docstring, or the body of
  a raw literal; tests excluded. The exclusions exist so the ceiling does not
  penalise the comment density this project asks for, and blank lines add no
  reading. Lines count as written: packing code by hand under `# fmt: skip` /
  `# fmt: off` is allowed, provided the packed lines stay readable and the
  review says why each new region is packed. This is the only statement of the
  rule; everywhere else points here.
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
Every code alteration must execute this strict pipeline, dispatched by the main session:
1. `@architect` defines interfaces and types.
2. `@tester` writes failing unit/async test cases *first* (including happy path and edge cases).
3. `@developer` writes the minimal code needed to pass the active tests.
4. `@tester` and `@reviewer` validate results and type consistency before merge readiness.

Step 4 is not conditional on the branch containing code, and the push that
would publish it is not yours to make unasked. `.claude/REQUIRED-READING.md`
rules on both — the approval boundary and when the assessment fires.

### The main session dispatches
The main session, started with no agent name, is the dispatcher: it spawns
the personas. Its rules:
* **Step order.** Design (`@architect`), failing tests (`@tester`), code
  (`@developer`), test run, review (`@reviewer`), and `@perf`'s acceptance run
  when refine or mesh code is touched (`docs/increments/README.md`).
  A performance fix is timed by `@perf` before review. Failing
  tests go back to `@developer`.
* **No asking between internal steps, but stop at the remote.** Do not ask
  Ola between steps (tests to code, code to review); loop until `@tester` and
  `@reviewer` are satisfied. `.claude/REQUIRED-READING.md` says which acts
  need a fresh yes.
* **Recap each round.** At every new round or increment, run
  `python3 tools/session_state.py` and open with its recap.
* **Milestone updates.** After each milestone, give Ola a one-line log of
  where the round stands (e.g. "`@tester` has 8 failing tests; on to
  `@developer`").
* **Report only finished, verified results, briefly.** Answer Ola's question
  first, then stop. No file that does not exist yet, no number from a run
  still in progress, no cause not checked. No unasked images, no undefined
  jargon.
* **Prohibited dependencies.** Hold every persona to §2's list.
* **Lessons.** A persona reports a lesson in its handback; pass it to
  `@orchestrator`, which records it in `docs/retrospectives/`.
* **When to spawn `@orchestrator`.** After each increment merges (a check of
  that increment, and its lessons recorded); the morning after each
  unattended night (idle time, guard refusals and false positives, work done
  out of role); and for a research round, weekly or when Ola asks.

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
```

### CI
The gates above also run in GitHub Actions (`.github/workflows/main.yaml`), and
CI is authoritative: local green does not mean the branch is green, so
**verify check status before declaring anything merge-ready**:
```bash
gh pr checks <pr>      # must be green; required for merge on master
```
`master` merges through a merge queue: `gh pr merge <pr>` enqueues the PR (the
queue's method, merge commit, applies), and the queue tests it on top of the
PRs ahead and merges it. No update-branch loop. A method flag such as `--merge`
only prints a warning (the queue's method wins); `-d`/`--delete-branch` is
refused. Never pass `--admin`: it bypasses the queue.

