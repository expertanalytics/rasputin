# Audit: the test suites, judged by what they do for control of the code

Status: audit, read-only, by `@tester`, against master `44fa7f5`. Every
`file:line` is at `44fa7f5`. Companion to the Python audit
(`docs/increments/python-audit.md` on branch `worktree-python-audit`, pinned
to `12dace7`); its X1 (CLI harness), X2 (private `cli` names) and X3
(import-firewall tests) are cited, not redone.

## How it was measured

- CI: run 37314235492 (the merge-queue run of PR 185, 2026-10-05), 11 min 57 s
  wall. Job and step times from `gh run view --json jobs`; per-suite TSan
  seconds from the job log's line timestamps; per-case ASan seconds from
  ctest's lines in the job log.
- Python: whole suite, main checkout's venv, run from the test-audit worktree
  with coverage on as in CI: 5216 passed, 14 skipped, 350 s locally, line
  coverage 98.6 %. Per-file sums from a second run with `--junitxml`.
- C++ Release: the last full ctest log on disk
  (`.claude/worktrees/23c-1/build/Testing/Temporary/LastTest.log`, 983 cases):
  13.1 s of test time in total. No C++ was built.

## Verdict

1. The suites are not out of control as a mass; they are out of control at a
   few seams. 75,000 test lines for 18,700 production lines is a lot, but the
   cost Ola feels comes from three places: copies of old implementations kept
   as oracles, tests that reach through `cli.py`'s privates and namespace, and
   two oracles that are quadratic under the sanitizers.
2. The CI wall time is set by tests, not compilers. In Python, 63 % of the
   suite's time is tests of the harness tools, run on all three legs. In C++, the TSan job (705 s) and
   the ASan job (522 s) are the two longest jobs, and in both the Test step is
   most of it. One test case (ES9) is the whole ASan test step; four
   TSan suites that never start a thread are 253 s of the TSan step.
3. Python line coverage is 98.6 %: the Python gap is not lines but the lack of
   byte-for-byte goldens on the paths the Python audit's PRs A and H will
   move. C++ has no line coverage measurement at all.
4. Exact-wording pins are not a problem: across the product tests there are
   19 asserted strings of 40 characters or more, and refusal tests mostly pin
   short key phrases.

## Findings, ranked (control of the production code first, then size and CI time)

Each: evidence, shape, test lines saved, CI seconds saved, what it does for a
refactor. "Measured" numbers come from the runs above; "estimate" numbers do
not and say why.

### R1. Tests reach through `cli.py` into its privates and its namespace (extends X2)

The Python audit's X2 counts six private `cli` names. The whole suite pins
more, and in two ways X2 does not list:

- Binding-layer tests build their start mesh through `cli` privates:
  `cli._engine` and `cli._constraint_arrays` at
  `tests/python/test_core_edge_strip.py@44fa7f5:63`,
  `tests/python/test_core_refine_points.py@44fa7f5:41`, and `cli._engine` also in
  `tests/python/feature_fixtures.py@44fa7f5:136`, `test_core_refine.py`,
  `test_refine_golden.py` (5 files in all). The lowest layer's tests depend on
  the top layer's private helpers, so the F8 move of `_engine` to
  `start_mesh.py` breaks the binding suites, not just the CLI ones.
- Collaborators are patched by their name inside `cli`'s namespace:
  `monkeypatch.setattr(cli, "refine", ...)` at 9 sites in 6 files
  (`tests/python/test_cli_constraint_feet.py@44fa7f5:81, 223, 245`, `tests/python/test_cli_start_quality.py@44fa7f5:74,
  216`, `tests/python/test_cli_mesh_stats.py@44fa7f5:317`, `tests/python/test_cli_mesh_edge_strip.py@44fa7f5:181`,
  `tests/python/test_cli_mesh_plain_output.py@44fa7f5:574`, `tests/python/test_refine_golden.py@44fa7f5:134`);
  `cli.open_dem` (`tests/python/test_cli_mesh_geographic.py@44fa7f5:358`,
  `tests/python/test_cli_mesh_domain_crs.py@44fa7f5:326`); `cli.final_check.run`
  (`tests/python/test_cli_mesh_geographic.py@44fa7f5:359`); `cli._dem_mesh`
  (`tests/python/test_cli_mesh_domain_crs.py@44fa7f5:299, 305`); `cli.hardening`
  (`tests/python/test_hardening.py@44fa7f5:73`). Each pins "which module calls `refine`". When F8
  moves the call into `mesh_run.py`, the patch lands on a name nobody calls
  any more, and the spy sees nothing: the test fails if it asserts the spy
  was called, and passes vacuously if it does not.
- Remaining X2 names, re-counted at `44fa7f5`: `cli._triangulated` 4
  (`tests/python/test_cli_draw.py@44fa7f5:652, 805`, `tests/python/test_cli_mesh.py@44fa7f5:96`,
  `tests/python/test_cli_mesh_vtk.py@44fa7f5:54`), `cli._chain_masks` 3, `cli._undirected` 1
  (`tests/python/test_cli_mesh_vtk.py@44fa7f5:175`); also outside `cli`:
  `dem_input._domain_plan` (`tests/python/test_cli_fetch.py@44fa7f5:48`), `DemTile._adopt`
  (`tests/python/test_mosaic.py@44fa7f5:1630-1662`, `tests/python/test_target_grid.py@44fa7f5:534-540`).

Shape: one test-side seam module, `tests/python/seams.py`, that is the only
file naming these: `start_mesh(xy, chains)` and `constraint_arrays(run)`
(today calling `cli._engine` / `cli._constraint_arrays`), and a
`spy_on_refine(monkeypatch)` fixture (today patching `cli.refine`). The ~25
sites call the seam. When F8 or E moves the code, one file changes, in the
same PR. Lines saved: about 40 (the nine spy bodies are near-copies). CI: 0.
Refactor: turns every python-audit PR that touches `cli.py` from a ~10-file
test edit into a one-file edit. Do it first.

### R2. CI: TSan runs thread-free suites, serially

TSan job: Build 126 s, Test 566 s, the longest job of the run. Per suite
(log timestamps): `prop_refinement_edge_strip` 165, `test_mesh_row_spans`
149, `test_mesh_lawson_stack` 68, `prop_refinement_refine_points` 61,
`prop_refinement_frozen` 46, `prop_refinement_refine` 24,
`test_mesh_lattice_split` 20, `prop_refinement_scan_equivalence` 16, the
other twelve under 4 s each.

`test_mesh_row_spans`, `test_mesh_lawson_stack`, `test_mesh_lattice_split`
and `prop_refinement_scan_equivalence` never start a thread: the only thread
code is `include/terrain/parallel_util/chunks.hpp`, reached only through
`refine.hpp` and `refine_points.hpp`, and none of the four includes either
(checked by grep of their includes and of the headers they include). TSan
cannot see a race in a one-thread process, so those 253 s buy nothing that
the ASan job (which runs them too) does not. The list is in
`.github/workflows/main.yaml` (the TSan Build and Test steps), and the loop
runs the binaries one after another.

Shape: drop the four from the TSan list; run the rest with `xargs -P
"$(getconf _NPROCESSORS_ONLN)"`, as the ASan job already runs ctest in
parallel. Test lines: 0 (workflow only; not a `tests/` change). CI: measured
-253 s serial; with the parallel run and R3, the TSan Test step is bounded by
its slowest remaining suite, estimate 60-70 s (the refine_points suite at 61 s
today, less after R3). The TSan job goes from 705 s to about 200 s.
Refactor: none directly; it shortens every loop.

### R3. CI: one quadratic oracle is the whole ASan test step

ASan job: Build 194 s, Test 317 s, and ctest's own total is 317.16 s because
one case, ES9 (`tests/cpp/property/prop_refinement_edge_strip.cpp@44fa7f5:679`), ran
316.8 s. In Release it runs 3.3 s. Next: RP3 92.6 s, `legalise_around ...
matches today's` 86.2 s, the row_spans T1 sweeps 22-50 s each.

ES9's J2 oracle `sources_over` (`tests/cpp/property/prop_refinement_edge_strip.cpp@44fa7f5:645-677`)
loops every output triangle against every source point (2,048 on the 33 x 33
grid), and does a `std::set<pair<double,double>>` lookup per (triangle,
point) pair before the orientation test, over 16 generator configurations.
RP3's `violations` (`tests/cpp/property/prop_refinement_refine_points.cpp@44fa7f5:142-173`) is the same
oracle, written a second time, with a bounding-box test, but also after the
set lookup. Under Debug + ASan the set lookups dominate (mechanism, not
profiled: no C++ was built).

Shape: one `j2_violations` in `tests/cpp/support/` (the strip oracle header
or a new `oracles.hpp`, see R5): map each point to the frame and decide
"start vertex" once (P work), bucket points by lattice cell, and visit only
the cells a triangle's box covers. Same verdict, same frame. Lines: about
-30 (two copies to one). CI: estimate ASan Test step from 317 s to about
90 s (then bounded by `legalise_around`, R4), TSan `edge_strip` and
`refine_points` suites down by most of their 226 s. Verifiable by the
`prop_refinement_edge_strip` time in the next CI log. Refactor: one oracle,
not two, for 15c's J2. `@tester` writes it; the oracle's ability to fail
(RP3's planted-z case, `tests/cpp/property/prop_refinement_refine_points.cpp@44fa7f5:259`) must still
pass against the new oracle.

### R4. Copies of old implementations kept as oracles after the change landed

Several suites pin "today's" behaviour by carrying a verbatim copy of an
earlier version of the production code:

- `tests/cpp/unit/test_mesh_lawson_stack.cpp` (281 lines): `legalise_around`
  as at `09b005b`, copied, compared insertion by insertion over 8 seeds x 3
  frames x 2 modes x 300 attempts (32.5 million assertions in 3 cases; TSan
  68 s, ASan 86 + 45 s).
- `tests/cpp/unit/test_refinement_active.cpp@44fa7f5:36` (109 lines): refine's sort +
  unique at `09b005b`.
- `tests/cpp/unit/test_mesh_lattice_incircle.cpp@44fa7f5:704-934`: `must_flip` as at
  `93a8066`, plus its grid and insertion helpers, for two cases.
- `tests/cpp/unit/test_raster.cpp@44fa7f5:384-417, 497`: the four-corner bilinear
  expression before increment 27.
- `tests/cpp/property/prop_noding_verifier_sweep.cpp@44fa7f5:837`: six recorded
  digests of `node<K>` at `6b4fcb9`.

Each was the right tool for its increment (prove a rewrite changed nothing).
Kept afterwards, each is a second implementation that must not be edited,
pins one algorithm's exact choices, and fails on any deliberate improvement,
with the reason in a comment naming a commit the code has long moved past.
The property oracles (Delaunay, tolerance, thread determinism, rule D) and
the golden digests at the public boundary (R6) carry the same protection
without pinning internals.

Not in this list: `tests/cpp/support/verifier_oracle.hpp` (the brute-force
guarantee check stays as the oracle by ruling R8 of increment 16b) and
`tests/cpp/support/scan_oracle.hpp` (the box walk is a simpler independent
algorithm, not a copy of the same one). Keep both.

Shape: retire the copies; keep the cases that test the new contract itself
(for `lawson_stack`, "no flip to make", `tests/cpp/unit/test_mesh_lawson_stack.cpp@44fa7f5:232`, and
"discards a non-empty stack", `:251`, rewritten to compare a stale-loaded
stack with an empty one, both production calls, over a few seeds instead of
the copy over 8 seeds x 300 attempts; together about 60 lines). The raster S3 copy and the noder digests cost
no measurable time; keep them, or move the digests into R6's goldens file. Each retirement needs
`@architect`'s yes, because the increment files name these tests. Lines:
about -470 (lawson_stack -240, active -80, must_flip copy and its case
-150). CI: ASan -86 s of the critical case after R3, TSan -68 s. Refactor:
the meshing kernel (`lawson.hpp`, `refine.hpp`'s active set) can then
change under the property oracles alone.

### R5. The Delaunay oracle is written eight times in `tests/cpp`

`delaunay_violations` at `tests/cpp/unit/test_mesh_lawson.cpp@44fa7f5:125`,
`tests/cpp/unit/test_mesh_quality_void.cpp@44fa7f5:127`, `tests/cpp/unit/test_mesh_quality.cpp@44fa7f5:280`,
`tests/cpp/support/strip_oracle.hpp@44fa7f5:317`; `delaunay_oracle` at
`tests/cpp/property/prop_refinement_refine_points.cpp@44fa7f5:175`,
`tests/cpp/property/prop_refinement_constraint_feet.cpp@44fa7f5:329`; `check_delaunay` at
`tests/cpp/property/prop_refinement_refine.cpp@44fa7f5:134`; inline at `tests/cpp/property/prop_refinement_quality.cpp@44fa7f5:197`.
Tolerance oracles: `tests/cpp/property/prop_refinement_refine.cpp@44fa7f5:165`,
`tests/cpp/property/prop_refinement_constraint_feet.cpp@44fa7f5:235`, plus the two J2 copies of R3. 17
to 38 lines each, about 210 in all. Rule D (`.claude/agents/tester.md` §3D)
names one reference implementation, but eight copies exist, in two frames
(lattice frame and world frame).

Shape: `tests/cpp/support/oracles.hpp` with `delaunay_violations(mesh,
frame)` for a `LatticeMesh` and for an output mesh, and
`tolerance_violations(dem, out, tol)`; the eight call it. Lines: about
-150. CI: compile time only, small. Refactor: a change to the output types
(`Outcome`, `PointRefineOutcome`) touches one oracle file, not eight; and
rule D's "decided in the producer's frame" is written once.

### R6. Gap: no byte-for-byte goldens on the paths the Python refactors move

Goldens exist: refine on the real tile and the quarter circle
(`tests/python/test_refine_golden.py@44fa7f5:52-55, 176-179, 190`), the stride-2 VTK
(`tests/python/test_cli_mesh_refine.py@44fa7f5:75`), feet on/off (`tests/python/test_cli_constraint_feet.py@44fa7f5:54-56`),
the catchment bowl (`tests/python/test_catchment.py@44fa7f5:418-419`). One value is recorded twice
(`a8e8720d...` at `tests/python/test_cli_constraint_feet.py@44fa7f5:54` and
`tests/python/test_refine_golden.py@44fa7f5:177`). None covers the written files of `rasputin mesh`
with features (rivers, lakes, land cover), with `--domain` in another CRS, or
on a mosaic, which is the code the Python audit's PR A (lattice arithmetic,
"meshes byte-identical") and PR H (`mesh_run.py`) move. Their gate today is
`@perf`'s 1 m benchmark, run by hand, on one input.

Shape: one `tests/python/goldens.py` holding every recorded digest and the
digest function, and about six new CLI goldens over the written `.vtk`/`.ply`
of committed fixtures (features, domain-CRS, mosaic, land cover, edge strip,
start quality), recorded at master before PR A. Lines: about +120. CI:
about +6 s (each run is 1 s or less on the committed fixtures). Refactor:
this is the safety net for PRs A and H; with it, "byte-identical" is a CI
check, not a manual run.

### R7. Gap: no C++ line coverage, and `bindings/core.cpp` is measured by nothing

`--cov=tin_engine` measures Python only. The C++ headers (8,800 lines) and
`bindings/core.cpp` (1,422 lines) have no coverage number. By include count,
the thinnest direct tests are `refinement/strip_scan.hpp` (166 lines; one unit
file tests one function, `test_refinement_coincidence_radius.cpp`),
`refinement/seam.hpp` (102; one property file), `hydrology/flood.hpp` (66;
only through `accumulate` and `upstream`), `noding/broad_phase.hpp` (84; one
property file). That says where to look, not that they are uncovered.

Shape: a local `llvm-cov` report (a `tools/` script or a CMake preset, not a
gate), run once before any C++ refactor PR to list unexecuted lines. Lines:
about +40 in tools. CI: 0 if local. Refactor: the C++ side gets the evidence
Python already has.

### R8. Tests that pin module layout, prose or the stub's text

- Stub text: `TestStubs` in `tests/python/test_core_cdt.py@44fa7f5:981`, `tests/python/test_core_noding.py@44fa7f5:925`,
  `tests/python/test_core_edge_strip.py@44fa7f5:602`, `tests/python/test_core_refine_points.py@44fa7f5:385`, and
  `tests/python/test_core_seam.py@44fa7f5:334`, `tests/python/test_core_frozen.py@44fa7f5:240` read `_core.pyi` as text
  and look for hand-picked declarations (about 90 lines). `python -m
  mypy.stubtest tin_engine._core` checks the whole stub against the module;
  run here it reports 88 errors, all pybind11 shapes (metaclass, `__init__`
  arguments, enum `__members__`/`__int__`/`__index__`, enum values typed as
  the enum), which a short allowlist absorbs. Replace the text tests by one
  stubtest step.
- Red-step scaffolding still in the suite: `test_the_module_exists`
  (`tests/python/test_viz_scene.py@44fa7f5:388`, `tests/python/test_viz_svg.py@44fa7f5:558`),
  `test_viz_package_directory_exists` (`tests/python/test_viz_protocols.py@44fa7f5:121`),
  `test_stub_file_exists` (`tests/python/test_core_cdt.py@44fa7f5:986`),
  `test_the_writer_has_moved_out_of_cli` (`tests/python/test_io_geojson.py@44fa7f5:77`).
  They were the red step's first failure; green, they test that a file exists.
- Prose pins: `test_the_package_docstring_says_so`
  (`tests/python/test_io_repository.py@44fa7f5:408`), `test_features_never_mentions_the_extension_at_all`
  (`tests/python/test_features.py@44fa7f5:605`, greps source text, X3 notes it),
  `test_describes_docstring_no_longer_claims_one_enumeration`
  (`tests/python/test_core_noding.py@44fa7f5:916`).
- Call-site pins: `test_adopt_is_called_only_from_the_mosaic_and_the_resampler`
  (`tests/python/test_mosaic.py@44fa7f5:1670`) and `test_exactly_one_from_crs_site_and_it_is_in_crs_py`
  (`tests/python/test_crs.py@44fa7f5:180`) name the modules allowed to call a function. These are
  layering rules; they belong in X3's `test_layering.py` table, where PR A's
  moves update one row.

Lines: about -150 (stub text 90, scaffolding and prose 60). CI: negligible.
Refactor: stops failing on moves and docstring edits that change no behaviour.

### R9. `raising=False` patches can go silently vacuous

`monkeypatch.setattr(..., raising=False)` on a name that does not exist adds
it and patches nothing. Today: `tests/python/test_edge_strip.py@44fa7f5:75` sets
`edge_strip._core`, which `edge_strip.py` never reads (it imports the names,
`src_python/tin_engine/edge_strip.py@44fa7f5:14`); `tests/python/test_hardening.py@44fa7f5:72-73`,
`tests/python/test_catchment.py@44fa7f5:286`, `tests/python/test_decompose_partition.py@44fa7f5:259`,
`tests/python/test_feature_input.py@44fa7f5:537`, `tests/python/test_target_grid.py@44fa7f5:443, 475, 491, 506` (the
last four target a name that exists, so `raising=False` only removes the
alarm a rename would ring). Shape: drop `raising=False` where the name exists;
patch through R1's seam where it is a "whichever spelling" patch. Lines: 0.
Refactor: a rename then fails loudly instead of leaving a test that checks
nothing.

### R10. Smaller duplication (Python)

- The GIL "ticker" is copied (`tests/python/test_core_cdt.py@44fa7f5:751`,
  `tests/python/test_core_noding.py@44fa7f5:721`), and its self-test runs twice (2.1 + 2.2 s);
  four more GIL tests time a call by hand (`tests/python/test_core_refine_points.py@44fa7f5:348`,
  `tests/python/test_core_edge_strip.py@44fa7f5:511-595`). One `gil_probe.py` and one parametrised
  "every long call releases the GIL" test: about -90 lines, -2 s, and one
  place to adjust the thresholds if a runner is slow.
- The real-tile quarter circle is meshed at 1 m by four suites, each asserting
  `achieved <= tolerance` plus one thing (`tests/python/test_cli_mesh_domain.py@44fa7f5:499`,
  `tests/python/test_cli_constraint_feet.py@44fa7f5:284`, `tests/python/test_cli_start_quality.py@44fa7f5:388`,
  `tests/python/test_refine_golden.py@44fa7f5:146`). A module-scoped run shared by the asserts would
  save about 5 s, in the codecs step only. Low priority.
- Helper copies outside X1's CLI harness: `repo` (14 definitions, 7 distinct
  bodies, 98 lines), `mesh` (9), `valley` (3). About -100 lines; fold into
  X1's PR.
- Lazy `importlib.import_module("tin_engine...")` in 28 test files (red-step
  scaffolding, so the suite collects before the module exists). Harmless at
  run time, but `grep "from tin_engine.x import"` and an IDE rename miss them.
  Convert to plain imports when a file is next touched; no PR of its own.

### R11. CI: tests of the harness tools are most of the Python suite's time, three times over

Per file, from the JUnit run (329 s of test time in all, 5,230 cases): the
tests of `tools/` and `.claude/hooks/` (`test_check_citations_pinning` 43 s,
`test_guard_push` 26, `test_brief` 25, `test_guard_targets` 20,
`test_guard_governance` 18, `test_guard_spawn` 15, `test_count_loc` 13,
`test_session_state` 12, and ten smaller) take **206 s, 63 %**, over 1,234
cases at a median 0.2-0.3 s each (each starts a subprocess or a git
repository). All `test_cli_*` take 53 s, all `test_core_*` 8 s. They are
8,099 non-blank lines, 18 % of `tests/python`. They do not exercise
`tin_engine`'s behaviour (two import it: `test_bench.py`, and `test_scratch_copy.py` to locate the installed package; those two stay in the product legs or install the package), add nothing to
`--cov=tin_engine`, and do not depend on the Python version, yet each of the
three Python legs runs them (Tests step 287 s on 3.12 in the reference run).

Shape: a `harness` marker (or the path list), excluded from the three product
legs, and run once in its own job beside them (install `pytest` only). Lines:
0 in tests (a marker line in `pyproject.toml` and the workflow). CI: estimate
-200 s on each Python leg; the 3.12 leg (332 s) drops to about 140 s, and
the new job (about 220 s) runs in parallel. Refactor: none; it is the
largest single cut in Python CI time.

## CI wall time if R2, R3, R4 and R11 land (estimate)

Reference run: 11 min 57 s, set by TSan (705 s), then ASan (522 s), then
Python 3.12 (332 s). After: TSan about 200 s (R2, R3), ASan about 194 s build
+ 50-60 s test (R3, R4; then bounded by the row_spans T1 sweeps at 22-50 s),
Python legs about 140-180 s, C++ ubuntu 186 s. The run would be set by ASan
at about 4.5 min. The estimates are bounded by the measured per-case times
above; the real figure is the next CI run's.

## Proposed PR order, after the Python audit's

The Python audit's order is T2, T1, B, A, F, C, D, E, G, H. These fit around
it; every one is tests or CI only (0 production lines).

| # | Branch | Takes | Test lines | CI seconds (wall) | Waits for | Before |
|---|---|---|---|---|---|---|
| Q1 | `audit-ci-sanitizers` | R2, R3 | about -30 | about -500 (TSan) and -230 (ASan) | nothing | anything else; a workflow change and one C++ oracle |
| Q2 | `audit-ci-harness-split` | R11 | 0 | about -200 per Python leg | nothing | can go with Q1 |
| Q3 | `audit-test-seams` | R1, R9 | about -40 | 0 | nothing | python-audit E, G, H (and fold into T1 if T1 is not yet open) |
| Q4 | `audit-goldens` | R6 | about +120 | about +6 | nothing | python-audit A and H |
| Q5 | `audit-retire-copies` | R4, R5 | about -620 | about -90 ASan, -68 TSan | `@architect`'s yes per retired test | none; anytime after Q1 |
| Q6 | `audit-stub-and-layout` | R8 (stubtest replaces stub text; call-site rules into X3's table) | about -150 | negligible | python-audit T2 (`test_layering.py`) | none |
| Q7 | `audit-cpp-coverage` | R7 | +40 in tools | 0 (local) | nothing | any C++ refactor |
| - | folded | R10 into python-audit T1 | about -190 | about -7 | T1 | - |

Totals: about -910 test lines net (about -1,030 removed, +120 goldens; +40 in tools), CI wall
from about 12 min to about 4.5 min (estimate), and two new safety nets for
refactors (Q4 goldens before PR A and H; Q7 C++ coverage before any C++
refactor).

## Rulings

The questions, as put to Ola in the main session's chat on 2026-10-05:

- "Once an increment has merged, retire tests that compare against a copy
  of older code, and let the property checks and recorded outputs guard
  it? Default: yes, with @architect confirming each one."
- "Run the harness-tool tests once per CI run in their own job? Default:
  yes."
- "should the CI-speed PRs go first? (recommended)"

2026-10-05, Ola, verbatim: "defaults on all, CI speed first". So:

1. **Tests that compare against a copy of older code** (R4) are retired once
   the increment that needed them has merged. The property oracles and the
   recorded digests carry the protection instead, and `@architect` confirms
   each retirement, test by test, because the increment files name them.
2. **Tests of the harness tools** (R11) run once per CI run, in a job of
   their own, not on every Python leg.
3. CI speed comes first: Q1 and Q2 (R2, R3, R11) are the next work. Their
   design is `docs/increments/h17-ci-test-time.md`, which also takes over
   h15's PR 2 (the asan+ubsan tests split over three jobs).

One correction, found while designing h17 (`gh run view 37314235492 --json
jobs`): the thread sanitizer job gates nothing (`CI result` does not need
it), so in the reference run `CI result` finished at 9.0 min, set by the
asan+ubsan job (8.7 min), not at the run's 11 min 57 s. R2 shortens the
run, not the wait for a green check. And the asan+ubsan Test step is
bounded by its total work over four test slots as well as by ES9 (h15
§3b), so R3 alone brings it to about 200 s, not 90 s. h17 §3 has the
figures.
