# Harness h17: CI time spent in tests

Status: design by `@architect` on master `44fa7f5`; design approved
(review round 3, below). PR 1 is in: red tests `a655aa0`, harness marker
`cc8c309`, oracle `b47c283`/`9d6105a`, rulings `f50a9c7` (§9), collection
hook `04e0090`/`bc69744`, CMake and NaN comments `6db86e3`, workflow
`00aa0eb`, fixes `f92c70e`. Code review round 1 asked for two prose fixes,
both taken. Next `@reviewer` round 2 (re-read the two lines), then push as
one PR on Ola's yes; acceptance (§4e) is read from the PR run. §8 question
1 runs on its default (the harness job on Python 3.12); questions 2
and 3 bind PR 2 only and stay open. Takes the test audit's Q1 and Q2 (`docs/increments/test-audit.md`,
R2, R3, R11) and h15's PR 2 (`docs/increments/h15-ci-speed.md`, §5 A and
§7), as two PRs (§5). Ola, 2026-10-05: "defaults on all, CI speed first".

**In one paragraph.** In the reference run (37314235492, the queue run of
PR #185, after h15's Python split) the wait for a green `CI result` was
9.0 min, set by the asan+ubsan job (8.7 min: Build 194 s, Test 317 s). The
slowest Python job (3.12) came next, ending 5.7 min into the run, 4.8 min
of it in the main suite. The
thread sanitizer job ran to 11.9 min but gates nothing. PR 1 replaces one
quadratic test oracle, moves the tests of the harness tools into a job of
their own, and trims the thread sanitizer list: estimated 9.0 to about
6.8 min on a slow runner. PR 2 then splits the asan+ubsan tests over three
jobs: about 4.9 min. No test is dropped from any gating job; every check a
run makes today is still made, on the same build.

## 1. Prior art: legacy and literature

*Literature.* The new oracle (§4) is uniform-grid bucketing for point
location: points are filed by the grid cell holding them and each query
(here a triangle's bounding box) visits only the cells it covers. This is
the standard grid method (Akman, Franklin, Kankanhalli, Narayanaswami,
"Geometric computing and uniform grid technique", *Computer-Aided Design*
21(7), 1989, 410-420). The departure: the grid here is the DEM's own
lattice, which the frame already uses, so no cell size is chosen. Nothing
new is claimed. The CI side (test sharding by stride, pytest markers) is
the mechanism h15 §1 already cites.

*Legacy.* Nothing to carry: the legacy tree had no CI workflow and no
bucketed oracle. The grep:

```
$ git grep -l -i -E "fsanitize|pytest.mark|bucket" legacy-archive -- legacy
legacy-archive:legacy/tests/test_polygons.py
$ git ls-tree -r --name-only legacy-archive -- legacy | grep -iE "\.github|workflow|\.yml|\.yaml"
(nothing)
```

The one hit is a `@pytest.mark.skip # Segfaults` on a polygon test.

## 2. The reference run, read job by job

From `gh run view 37314235492 --json jobs` (start and end per job and
step) and the job logs. Minutes from the run's first job (13:06:14Z).

| job | ends at | inside |
|---|---|---|
| C++ sanitizers (asan+ubsan) | 8.9 | Build 194 s, Test 317 s; 994 tests; per-test sum 1179 s |
| Python 3.12 | 5.7 | Install 37 s, Tests 287 s (4996 passed, 119 skipped, coverage 98.62 %) |
| Python 3.13 / 3.14 | 3.9 / 5.2 | |
| Python extras, 3.12 / 3.13 / 3.14 | 3.4 / 3.8 / 3.5 | |
| C++ core ubuntu / macOS / unchecked | 3.2 / 2.1 / 3.3 | |
| **CI result** | **9.0** | needs every job above, not TSan |
| C++ thread sanitizer | 11.9 | Build 126 s, Test 566 s, 20 suites one after another |

The asan+ubsan Build of 194 s puts this run on h15 §3b's slow runner kind
(Build 188-204 s; the fast kind 115-134 s). Per-test, under asan+ubsan:
ES9 317 s, RP3 93 s, `legalise_around ... matches today's` 86 s, the
`row_spans` T1 sweeps 22-50 s, RP3's control 30 s.

Two corrections to the test audit's estimates follow, and are recorded in
its Rulings section:

- **TSan does not set the wait.** R2 shortens the run and frees a runner;
  it moves `CI result` by nothing.
- **The asan+ubsan Test step has two bounds**, as h15 §3b found: the
  longest test (ES9) and the total work over the four test slots
  (1179 / 4 = 295 s). R3 removes the first, not the second.

## 3. Estimates (to be replaced by PR 1's run)

Method as h15 §5: each job's end is its start offset plus the steps it
keeps; the Test step is the larger of the longest test and the per-test
sum over 4.

- **asan+ubsan after R3.** ES9, RP3 and RP3's control (440 s together)
  fall to an estimated 50 s (the oracle was most of their time: ES9 runs
  3.3 s in Release). Per-test sum about 790 s, so Test about 200 s and the
  job about 6.6 min on a slow runner, about 3.8 on a fast one.
- **Python 3.12 after R11.** The harness tests were 63 % of the local
  suite's time (test audit R11); at that share the main suite drops from
  287 s to about 105 s and the job to about 2.6 min. The new harness job:
  install about 10 s, tests about 180 s, about 3.3 min.
- **TSan after R2 and R3.** Eight thread-free suites gone (§4c; the first
  four took 253 s, the other four's times not read), the remaining 12
  four at a time: Test at most about 70 s, job about 3.4 min. Off the
  critical path either way.
- **PR 1, `CI result`:** about 6.8 min (slow runner), about 4.0 (fast).
- **PR 2, three shards:** each shard Build 194 s plus about 75 s (790 /
  12 plus 10 % for an uneven split, but at least the longest test, 86 s):
  about 4.9 min slow, 3.7 fast. Below that the Python extras job (3.4-3.8)
  and the C++ core jobs (about 3.2) are the floor.

## 4. PR 1: the oracle, the harness job, the TSan list

Production code (`CLAUDE.md` §2 count): **none**. Workflow, `pyproject.toml`
and test files only, none counted. About 30 workflow lines, about 120 test
lines added and 60 removed.

### 4a. One J2 tolerance oracle, bucketed (R3) — `@tester`

Today `sources_over`
(`tests/cpp/property/prop_refinement_edge_strip.cpp@44fa7f5:645-677`) and
`violations` (`tests/cpp/property/prop_refinement_refine_points.cpp@44fa7f5:142-173`)
are one oracle written twice: every output triangle against every check
point, with a `std::set` lookup per pair. Both move to one header,
`tests/cpp/support/j2_oracle.hpp`, namespace `j2_oracle`, with no Catch2
include:

```cpp
struct Findings { std::size_t over = 0; std::size_t not_ccw = 0; };

// Frame: (col, -row), lattice units, as both callers already build it.
// Returns the (point, triangle) pairs over tolerance, exactly as today.
Findings violations(std::span<const Point2> vertices,          // output, framed
                    std::span<const double> z,                  // output z (or a planted copy)
                    std::span<const std::uint8_t> valid,
                    std::span<const std::array<std::uint32_t, 3>> triangles,
                    std::span<const Point2> points,             // check points, framed
                    std::span<const float> point_z,
                    std::span<const std::uint8_t> skip,         // 1 = a start vertex (J2 excludes it)
                    std::size_t cols, std::size_t rows, double tol);
```

Rules the header keeps, so the verdict is today's:

1. `skip` is computed by the caller once per point, from world coordinates,
   as today's set lookup does, not once per (triangle, point) pair.
2. Points are filed in `cols x rows` buckets by `(floor(col), floor(row))`,
   `row = -y`. A point outside `[0, cols-1] x [0, rows-1]` throws
   `std::invalid_argument`: no point is silently left out.
3. A triangle visits the buckets from `floor` of its box's low corner to
   `floor` of its high corner, inclusive. A point in the closed triangle is
   a convex combination of its corners, so each coordinate lies in the
   box exactly (min and max of doubles are exact), and `floor` is monotone:
   its bucket is visited. The pairs the buckets skip are pairs the
   orientation test rejects, so the set of pairs counted is today's.
4. The orientation test (`DefaultKernel::orient2d`, closed), the plane
   expression in today's operand order, `zmax = max(1, max |z|)` over the
   `z` given, and `tol + 1e-9 * zmax` stay character for character.
5. `not_ccw` counts triangles not counter-clockwise, over all triangles,
   valid or not (RP3's `REQUIRE` becomes `REQUIRE(f.not_ccw == 0)`).
6. A triangle's bucket range is clamped to the grid, `[0, cols-1] x
   [0, rows-1]`, before it is walked. Rule 2 puts every point inside that
   range, so the clamp drops no bucket that holds a point; it only keeps a
   triangle that reaches past the last column or row from indexing outside
   the bucket array.

Cases (in `prop_refinement_refine_points.cpp`, so no new CMake target):

- RP3's control (`tests/cpp/property/prop_refinement_refine_points.cpp@44fa7f5:259-270`,
  z shifted by twice the tolerance) still reports violations, through the
  new header.
- New, a hand-built case with a known count: a 3 x 3 node grid, two
  triangles, check points on a vertex, on an edge, on a cell boundary
  inside a triangle, on the last column and last row, and one marked
  `skip`; expected `over` worked out by hand for a planted plane.
- New: a point outside the grid throws.

Before deleting the two local copies, `@tester` runs RP3's control
through both the old copy and the new header and records in the handback
that the two return equal `over` counts at each of its three tolerances
(0.5, 2.0, 8.0). It does the same for one ES9 generator combination
(`Terrain::Smooth`, seed 1, tolerance 0.5) through ES9's old copy,
`sources_over`, and the new header: equal counts on ES9's output (0) and
on that output with z shifted by twice the tolerance (more than 0, so the
comparison can fail). Only then: ES9 and RP3 call the header; the two local
copies go. The ES9 body and its
16 generator combinations are not split (h15's four-case split is not
needed once ES9 is fast, if Ola says yes to §8 question 3; §5 says when it
comes back).

### 4b. The harness tests in their own job (R11) — `@tester`, then `@developer`

A test file is a **harness test** when it tests a file under `tools/` or
`.claude/` and imports nothing from `tin_engine`. At `44fa7f5` that is the
16 files `test_away`, `test_brief`, `test_check_citations`,
`test_check_citations_pinning`, `test_ci_changes`, `test_count_loc`,
`test_guard_governance`, `test_guard_push`, `test_guard_spawn`,
`test_guard_targets`, `test_guard_unattended`, `test_harness_mode`,
`test_rule_sizes`, `test_session_state`, `test_settings_wiring` and
`test_spawn_rule_text` (found by: references `tools/`, `.claude` or `hooks`,
and zero `tin_engine` imports; `@tester` confirms the list). `test_bench`
and `test_scratch_copy` import `tin_engine` and stay where they are.

- `@tester`: `pytestmark = pytest.mark.harness` in each file; the marker
  registered under `markers` in `pyproject.toml`'s pytest section
  (`--strict-markers` requires it; test configuration, so `@tester`'s).
- `@developer`, `main.yaml`: the `python` job's command becomes
  `pytest -m "not harness"`; a new job `harness` (name `Python harness
  tools`), gated on `changes` like `python`, Python **3.12** (question 1),
  checkout with `fetch-depth: 0` (`test_count_loc` recounts recorded PRs),
  installs `pytest`, `pytest-asyncio` and `pytest-cov` at the `dev`
  extra's bounds and **not** the package, and runs
  `pytest --no-cov -m harness`. It joins `CI result`'s `needs` and its
  shell list (h11's T4 tests already fail otherwise).
- `@tester`, before `@developer` (added by §9 A): pytest imports every
  test file before `-m` deselects, so without the package the harness
  job fails at collection on every product test file (86 errors in
  `@tester`'s trial). A `pytest_ignore_collect` hook in
  `tests/python/conftest.py` closes this. When the mark expression is
  exactly `harness` (`config.getoption("markexpr") == "harness"`), it
  returns `True` for a `test_*.py` file whose text does not contain
  `mark.harness`, and `None` (pytest decides) for everything else; any
  other expression, or none, leaves collection as today. Text, not the
  module's AST, on purpose: the text test over-approximates "carries the
  harness mark", so its only error is to collect a file it need not, which
  fails loudly on import. A file can carry the mark without the text only
  through `getattr` or a hook that adds marks; the suite uses neither
  (`grep -n "add_marker\|pytest_collection_modifyitems" tests/python/*.py`
  finds none). The test, in a new harness-marked file: under `-m harness`
  a fixture directory with a marked file and an unmarked file that imports
  a module that does not exist passes, and the unmarked file is never
  imported (red at `9d6105a`: collection error); a marked file importing
  that module still errors (the loud mode stays); under `-m "not harness"`
  the unmarked file is imported. Test configuration under `tests/`, so
  `@tester`'s.

Why this passes h15 §4's test: the two marker expressions are complements,
so every collected test runs in exactly one of the two jobs, in every run;
the hook ignores only files that hold no harness-marked test, which
`-m harness` would deselect anyway, so it removes nothing the harness job
would run;
coverage is unchanged because harness tests never import `tin_engine`.
What it gives up, by Ola's ruling: harness tests run on one Python version,
not three. Two failure modes stay loud: a harness test that needs
`tin_engine` fails on import in the harness job (the package is not
installed); a new harness file without the marker runs in the product legs
(slower, not unsafe). `pytest` run locally with no `-m` still runs both.

### 4c. The thread sanitizer list (R2) — `@developer`

Drop `test_mesh_row_spans`, `test_mesh_lawson_stack`,
`test_mesh_lattice_split` and `prop_refinement_scan_equivalence` from both
the Build targets and the Test list: none can start a thread. Checked by a
walk of each suite's transitive `#include`s for `std::thread`,
`std::jthread`, `std::async`, `<thread>`, `<future>` or `pthread`: zero files
for the four, while the controls `prop_refinement_refine` and
`test_refinement_chunks` each find `parallel_util/chunks.hpp` (script run
against `44fa7f5`; `src/` has no thread code). Drop also
`test_refinement_scan`, `test_mesh_lawson`, `test_refinement_scan_offnode`
and `test_mesh_quality` (§9 B): the same walk, rerun at `9d6105a` over all
20 suites, finds zero files for these four and for the first four, and
`parallel_util/chunks.hpp` for each of the other 12; the only file under
`include/` or `src/` that matches is `chunks.hpp`. Run the remaining 12 suites
`getconf _NPROCESSORS_ONLN` at a time, each suite's output held and printed
whole, a failing suite named, and the step failing if any suite fails
(`TSAN_OPTIONS=halt_on_error=1` stays). `@tester` corrects the comment at
`tests/cpp/CMakeLists.txt@44fa7f5:285` ("All three run in the TSan job"),
and the one at `tests/cpp/CMakeLists.txt@44fa7f5:217-219`, which reads as
if all three suites start threads: `test_mesh_lattice_split` links
`Threads` but starts none (the walk above).

Constant: **four at a time**, the runner's CPU count. Assumes 12 suites,
the longest about 60 s after R3 (refine_points 61 s today); not checked for
memory, since four TSan processes on a 16 GB runner are an assumption until
PR 1's run.

### 4d. Red tests — `@tester`, before `@developer`

In `tests/python/test_ci_changes.py` (it parses the workflow already):

- H1: a job `harness` exists; its pytest step carries `-m harness` and
  `--no-cov`; no step installs the package (`pip install -e` or `.`); its
  checkout has `fetch-depth: 0`. Red at master.
- H2: the `python` job's pytest line carries `-m "not harness"`, the exact
  complement of H1's expression. Red at master.
- H3: the TSan Build target list and Test suite list are the same set, and
  each name is an `add_terrain*_test` target in `tests/cpp/CMakeLists.txt`.
  Green at master; it keeps a later edit from building a suite it does not
  run, or the reverse.

That the TSan list drops exactly the eight, and that the parallel loop
reports failures, is `@reviewer`'s check, by `git diff` and by reading the
step, not a test that freezes the command (h15 §7's rule).

### 4e. Acceptance, on PR 1's own run

Read with `gh run view <id> --json jobs` (job and step times), and the job
logs (`gh run view <id> --log --job <job id>`: pytest's summary line,
ctest's `Test #N ... sec` lines and its `out of N`, each TSan suite's
`All tests passed`). The runner kind is read from the asan+ubsan Build
time (slow: 180 s or more; fast: 140 s or less); compare with this run's
figures on a slow runner, h15 §3b's on a fast one. A Build time between
the two (141-179 s) names no kind: judge it against the slow bounds, and
say in the acceptance note that the kind was not identified.

- **Same checks:** per Python version, main-suite passed and skipped plus
  the harness job's equal the base commit's run (4996 and 119 on 3.12 at
  `44fa7f5`) plus the tests the branch adds, counted from its diff;
  coverage on 3.12 unchanged (98.62 % at `44fa7f5`); ctest's `out of N`
  equals the base's 994 plus the cases the branch adds; 12 TSan suites,
  each `All tests passed`.
- **ES9 and RP3** each run shorter, under asan+ubsan, than
  `legalise_around ... matches today's` in the same run (86 s today), so
  neither is the longest test any more.
- **Times, from the run's first job:** the asan+ubsan Test step at most
  240 s (slow) or 130 s (fast); Python 3.12 at most 3.5 min; the harness
  job at most 4.5 min; `CI result` at most 7.5 min (slow) or 5.0 (fast).

Expected on PR 1's run, from code review round 1's local counts: on
Python 3.12, main suite plus harness job 5126 (4996 + 119 at `44fa7f5`,
plus the 11 tests the branch adds); ctest `out of 996` (994 plus 2); 12
TSan suites.

If ES9 misses its bound, `@perf` profiles it under ASan before PR 2, and
PR 2 takes h15's four-case ES9 split back.

## 5. PR 2: the asan+ubsan tests in three jobs (h15 §5 A)

h15 §5 A and §7 PR 2 as designed, with two changes:

1. **No ES9 split** unless PR 1's acceptance says so (§4e), if Ola says
   yes to §8 question 3.
2. **The order with h12**, if Ola says yes to §8 question 2. h15 §6 rule 3 put PR 2 after h12's PR A (the
   queue skip, C2), because PR 2 renames the sanitizer job and C2 looks for
   it by name. h12's PR A has not merged (branch `worktree-h12-design`, not
   an ancestor of master at `44fa7f5`; its review fixes are still queued).
   The rule no longer needs to hold: h12's own test
   `test_t7_names_the_tool_looks_for_are_job_names` reads the job names from
   the workflow and fails once the job is renamed, so whichever of the two
   merges second carries the name change, loudly. If PR 2 is first, h12
   then requires all three shard names, `success`, in one check suite with
   `CI result` (h15 §5 A), not one shard's name. If h12 is first, PR 2
   changes `tools/ci_changes.py` as h15 §7 says (about 15 lines,
   `@developer`, with `@tester`'s T7 and T8 changes first). PR 1 adds the
   `harness` job, which h12 must also skip on a tested tree; same rule.
   The loud check there is h12's `test_t4_every_other_job_is_gated_on_changes`
   (branch `worktree-h12-design`, `tests/python/test_ci_changes.py`): every
   job gated on `changes` must carry h12's skip condition, so a `harness`
   job without it fails that test. h12's list of jobs that skip
   (`CODE_JOBS` in the same file) was written before h15's PR 1 added
   `python-extras`, and needs `python-extras` and `harness` added by
   whichever merges second.

`@tester` first: a test that the sanitizer job's matrix lists shards 1 to
k and its ctest line runs `-I ${{ matrix.shard }},,k` with that k and
`--no-tests=error` (the union of the strides is every test). Then
`@developer`: the matrix. About 10 workflow lines.

Constant: **k = 3 shards**, as h15 §5 A, now on a per-test sum of about
790 s instead of 1179 s. Re-derive from PR 1's measured sum: the smallest
k for which sum / (4 k) x 1.1 is at most the longest test or 90 s,
whichever is larger; 3 at 790 s. Checked only at 994 tests.

**Acceptance, on PR 2's run:** the shards' `out of N` add up to PR 1's
count, plus 3 if PR 2 takes the ES9 split; each shard's Test step at
most 2 min; `CI result` at most 5.5 min (slow) or 4.5 (fast).

## 6. Who does what

| step | persona | PR | files |
|---|---|---|---|
| J2 oracle, ES9 and RP3 rewired, oracle cases, CMake comment | `@tester` | 1 | `tests/cpp/support/j2_oracle.hpp`, the two property files, `tests/cpp/CMakeLists.txt` |
| harness marker in 16 files, marker registered | `@tester` | 1 | `tests/python/test_*.py`, `pyproject.toml` |
| collection hook and its test (§4b, §9 A); CMake comments for §9 B; NaN comment (§9) | `@tester` | 1 | `tests/python/conftest.py`, a new `tests/python/test_*.py`, `tests/cpp/CMakeLists.txt`, `tests/cpp/support/j2_oracle.hpp` |
| H1-H3 | `@tester` | 1 | `tests/python/test_ci_changes.py` |
| harness job, `-m "not harness"`, TSan list and parallel loop | `@developer` | 1 | `.github/workflows/main.yaml` |
| acceptance read from the PR run | `@reviewer` (with the main session's `gh` output) | 1 | none |
| shard test | `@tester` | 2 | `tests/python/test_ci_changes.py` |
| matrix; `tools/ci_changes.py` only if h12 merged first | `@developer` | 2 | `main.yaml`, maybe `tools/ci_changes.py` |

No `@perf` run: neither PR touches refine or mesh code.

## 7. Exclusions

- No flag, build type or environment change in any gating job (h15 §8).
- No change to which checks are required.
- The other test audit findings (R1, R4-R10) are later PRs; R4's
  `legalise_around` copy (86 s) is the next bound after PR 2.
- ccache stays h15's PR 3, decided after PR 2's times.

## 8. Questions for Ola

1. **Which Python runs the harness tests?** Default: **3.12**, the floor
   the governance job and mypy already use. Your Mac runs the hooks on
   3.14, so a 3.14-only break in a hook would show there first, not in CI.
2. **PR 2 before h12's queue-skip PR?** Default: **yes**, whichever lands
   second takes the shard names; h12's own test makes that impossible to
   miss. A yes changes your h15 ruling 2 of 2026-10-05, which put the
   shards after h12's queue-skip change (h15, Ola's rulings).
3. **PR 2 drops h15's four-way split of ES9 (the 317 s edge-strip test)
   unless PR 1's measured run shows it is still needed.** Default: **yes**.
   A yes also changes your h15 ruling 2, which accepted the split as part
   of PR 2 (h15, Ola's rulings).

## 9. PR 1's test side: findings and rulings

`@tester`'s commits `a655aa0` (H1-H3, red), `cc8c309` (the marker on the
16 files, registered in `pyproject.toml`), `b47c283` (the shared oracle,
`tests/cpp/support/j2_oracle.hpp`) and `9d6105a` (the clamp case). The
equal-counts step of §4a held in all five cases (`b47c283`'s message: RP3's
control 11422 / 8537 / 2193 pairs at 0.5 / 2 / 8, old copy and header
alike; ES9 0 and 0 on its output, 8620 and 8620 with z shifted). Locally
under asan+ubsan, by `@tester`: ES9 52.6 s to 3.3 s, RP3 15.1 s to 0.66 s.

Rulings, made on defaults while Ola is away:

- **A. The harness job cannot collect without the package.** Ruled: the
  `pytest_ignore_collect` hook, specified in §4b, written by `@tester`
  with its test. Not chosen: installing the package in the harness job
  (it would cost the job the install time, about 37 s on the reference
  run, and lose the loud failure of a harness test that imports
  `tin_engine`); moving the harness tests to a directory of their own
  (the cleanest boundary, but it moves 16 files and their fixture modules,
  splits `conftest.py`'s prose check over two directories, and breaks
  every unpinned citation of those paths; worth it only if the marker
  route proves fragile). `@tester`'s trial of the hook, in a venv holding
  only the three test tools: 1078 passed, 12 skipped in 143 s; the 12
  skips came from a scratch copy that was not a git work tree.
- **B. Four more thread-free TSan suites.** Ruled: PR 1 drops them too
  (§4c, with the walk's result). TSan finds data races; a suite that
  starts no thread has none to find, and asan+ubsan still runs all eight
  on every PR. No rule asks invariant-critical suites to run under TSan
  (`grep -n -i tsan testing.md`). `@tester` corrects the CMake comments that
  then say otherwise: `tests/cpp/CMakeLists.txt@9d6105a:217-221` (it names
  only `test_refinement_chunks` as threaded, which stays true, but says
  `test_mesh_lattice_split` alone leaves the job) and
  `tests/cpp/CMakeLists.txt@9d6105a:261` ("Both run in the TSan job";
  `test_mesh_quality` no longer does). H3 is unchanged: it compares the
  two lists, whatever they hold.
- **C. A citation the marker block moved.**
  `docs/increments/29-nve-reference-catchments.md` cited line 117 of
  `pyproject.toml`, unpinned; the marker block moved the quoted line
  (`exclude = ["docs/benchmarks"]`) to 122. Pinned to `0130889f`, the
  commit that wrote the citing line, by `@architect`, with a note in the
  record's lead-in. Of the other 13 citations `check_citations.py` lists
  as at risk on this branch, six cite lines above every hunk this branch
  makes in their file, so it does not move them; the other seven
  (`test_brief.py`, `test_away.py` and `prop_refinement_edge_strip.cpp`
  lines) were already off their quoted lines at `fe9fec1`, before this
  branch, and most name their commit in the prose. This branch changes
  what none of them says.

Departures `@tester` took, recorded and accepted:

- The marker registration in `pyproject.toml` is `@tester`'s, as §4b and §6
  assign; not a departure.
- `j2_oracle::violations` throws `std::invalid_argument` on more than rule
  2's off-grid point: a NaN check point, mismatched span sizes, an empty
  grid, and a triangle naming a vertex past the end. Each turns a silent
  wrong read into a loud failure; none changes a count on valid input.
- A NaN corner. The header's comment
  (`tests/cpp/support/j2_oracle.hpp@9d6105a:49-51`) says a triangle with a
  NaN corner visits every bucket. That holds only on an axis where
  `std::min` or `std::max` returns the NaN, which depends on the corner's
  position (the leftmost element that no other is less than: `min({NaN, 1,
  2})` is NaN, `min({1, NaN, 2})` is 1). The verdict is today's all the
  same, for another reason: with a NaN corner the plane is NaN, and
  `abs(NaN - z) > slack` is false, so such a triangle counts no pair over
  in either the old copies or the header, and `not_ccw` walks every
  triangle regardless. `@tester` rewrites the comment to give that reason.

## Review

### Round 1: `@reviewer`, design, `e7bbec4` on `44fa7f5`

`@reviewer`'s record, word for word (one citation pinned to `44fa7f5` by
`@architect`):

> Design review round 1 (`@reviewer`, `e7bbec4`): CHANGES REQUESTED, five [now] items: h15 paragraph contradicts its new status line; §4e runner-kind gap of 141-179 s and the 5.5 vs 5.7 min figure; test-audit Rulings record no questions or defaults; record old = new oracle counts on RP3's control before deleting the copies; citation 258 to 259. Checked against run 37314235492: CI result waits only on non-TSan jobs, 9.0 min, ES9 317 s, sum 1179 s; four TSan suites reach no thread code; the oracle argument holds; the 16 harness files pass without `tin_engine` (1083 passed).
>
> Later: j2_oracle clamp bucket range to grid; tests/cpp/CMakeLists.txt@44fa7f5:217-218 "All three link Threads" false for test_mesh_lattice_split; cite h12 test_t4 as the loud check for the harness job, note CODE_JOBS predates python-extras; ~17 jobs may hit GitHub's concurrent-job limit (unchecked).

Fixed in the round-2 commit: the five [now] items (h15's order sentence
and §6 rule 3 marked superseded by §5 here; §4e's 141-179 s rule and one
figure, 5.7 min, for the Python 3.12 job's end; the test audit's Rulings
now carry the questions and defaults as put to Ola; §4a's equal-counts
record; the citation), and three of the later ones: the clamp (§4a rule 6),
the CMake comment (§4c), h12's T4 test and `CODE_JOBS` (§5). Left: the
concurrent-job limit, unchecked; PR 1's run would show it (a job queued, not started,
at the run's start).

### Round 2: `@reviewer`, design, `e7bbec4..edaf4e0`

`@reviewer`'s record, word for word:

> Design review round 2 (`@reviewer`, `e7bbec4..edaf4e0`): CHANGES REQUESTED, two [now] items: docs/increments/test-audit.md@edaf4e0:370 misquotes the third question; quote the chat's sentence verbatim ("And should the CI-speed PRs (sanitizer lists, the fast check, the harness job) go first, before the rest of the audit? I'd recommend it, given what CI time does to your flow."); docs/increments/h17-ci-test-time.md@edaf4e0:276 drops the ES9 split that Ola's h15 ruling 2 (docs/increments/h15-ci-speed.md@edaf4e0:728) accepted, but no question puts this to Ola (docs/increments/h17-ci-test-time.md@edaf4e0:339-342) and docs/increments/h15-ci-speed.md@edaf4e0:3 states it unconditionally; put it to Ola (default yes) and make h15's clause conditional. Round 1's five [now] items and the three later ones taken: fixed and true. Round 1's "All three link Threads false" was wrong (tests/cpp/CMakeLists.txt@44fa7f5:225-227); the design's "links Threads but starts none" is right. check_citations clean.

Fixed in the round-3 commit: the test audit's third question now quotes
the chat's sentence, checked against this session's transcript; §8
question 3 puts the dropped ES9 split to Ola (default yes), and h15's
status line, §4a and §5 item 1 make it conditional on that yes. Both
suggestions taken: §4a's equal-counts step also runs one ES9 combination
through `sources_over` and the header, and §5 item 2 is conditional on §8
question 2.

### Round 3: `@reviewer`, design, `edaf4e0..e09f466`

`@reviewer`'s record, word for word:

> Design review round 3 (`@reviewer`, `edaf4e0..e09f466`): APPROVED. Round 2's two [now] items fixed and true: test-audit's third question matches the chat verbatim; §8 question 3 puts the ES9 split to Ola, and h15's status line, §4a and §5 item 1 are conditional on it. §4a's ES9 check holds against tests/cpp/property/prop_refinement_edge_strip.cpp@e09f466:682-705. h15 rulings 1-5 are consistent with h17. Fold in: docs/increments/h17-ci-test-time.md@e09f466:349 "(h15, Review)" should read "(h15, Ola's rulings)"; docs/increments/h17-ci-test-time.md@e09f466:315 needs "plus 3 if PR 2 takes the ES9 split". check_citations clean. Ready for `@tester` on PR 1 once Ola answers §8, or on its defaults.

Fixed in the round-4 commit: both fold-ins (§8 question 2 now points at
h15's "Ola's rulings" section; §5's acceptance count adds 3 when PR 2
takes the ES9 split, one test becoming four). Design approved; the next
step is `@tester` on PR 1.

### Code review round 1: `@reviewer`, `44fa7f5..00aa0eb`

`@reviewer`'s record, word for word:

> Code review round 1 (`@reviewer`, `44fa7f5..00aa0eb`): CHANGES REQUESTED, two [now] items, both stale prose: tests/python/test_ci_changes.py@00aa0eb:703 "drops exactly the four thread-free suites" is eight since ruling B; docs/increments/h17-ci-test-time.md@00aa0eb:3-8 status line names steps already done. Checked and true: 0 production LOC; j2_oracle against §4a rules 1-6, the equal-counts record (b47c283, quoted in §9), both old copies gone, hand counts 9/4/2/10/11 worked by hand; ctest 996/996 in build and build-san from local logs (not rebuilt by the reviewer); marker on exactly the 16 files plus the hook test; split 1094 + 4031 = 5125 = 5114 + 11 new; harness command 1094 passed without the package (Python 3.13); full suite 5008 passed, 118 skipped; harness job, python job's -m "not harness", TSan 12 = 12 with exactly the eight dropped, parallel loop's failure reporting simulated; gates clean; 21 at-risk citations re-read. Push as one PR. After the two fixes: ready to push on Ola's yes; acceptance (§4e) is read from the PR run.

Fixed: item 1 (the test comment now says eight suites) in `f92c70e`
(`@tester`); item 2 (the status line) in this commit, which also records
the expected acceptance counts in §4e.
