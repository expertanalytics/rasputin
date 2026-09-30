# Rule-file history (moved out on 2026-09-30)

Rule files (`CLAUDE.md`, `.claude/agents/*`, `.claude/skills/*`,
`docs/increments/README.md`, `docs/PRINCIPLES.md`, `testing.md`) carry rules and
the commands that check them, not the history of how each rule was learned.
Ola, 2026-09-30: the setup is to be reused on other projects, where the history
means nothing. This file keeps that history, word for word where it was a
passage, under the file and section it came from. `.claude/REQUIRED-READING.md`
has its own record, `2026-09-27-required-reading-incidents.md`.

## `CLAUDE.md`

**§2, Strict Size Limit.** The exclusions (blank lines, comments, docstrings,
raw-literal bodies) were Ola's ruling on 2026-09-28; before that, 20b counted
blank lines. Hand-packing under `# fmt: skip` / `# fmt: off` was allowed by
Ola on 2026-09-27, on `tools/bench.py`.

**§4, CI.** "A whole branch once merged with CI red on every commit because
the workflow was never checked" is the reason for "verify check status before
declaring anything merge-ready".

## `.claude/agents/architect.md`

**§4, Literature and Novelty.** This is rule 1 of
`docs/retrospectives/2026-09-27-increments-14-to-20b.md`. The examples behind
it: "Increments 14, 14b and 18 rebuilt Garland and Heckbert 1995 without
reading it; increment 20 used Chew/Ruppert refinement but snapped to DEM nodes,
which drops its guarantees, and did not consider off-centres (Üngör)." And "The
DDT flip rule was designed against the literature just cited
(`docs/research/data-dependent-triangulation.md`)." The rule about unchecked
novelty claims exists because Ola wants the option to publish kept open.

## `.claude/agents/developer.md`

- A test that looks wrong is reported as a specification disagreement: "that
  has produced real design fixes in every increment so far".
- `-ffp-contract=off`: "three tests once passed only because clang contracted a
  determinant to a single fma, and the ubuntu CI leg would have gone red."

## `.claude/agents/orchestrator.md`

- Step 5 (`@reviewer`) "was skipped across four consecutive merges".
- The recap is rule 5 and "report only what is finished" is rule 4 of the
  2026-09-27 retrospective (`2026-09-27-increments-14-to-20b.md`, which has the
  incidents). "Ola works across days and is not always fully present; the
  recap is what lets Ola pick the thread back up."
- "Guard the Context" named CGAL, GDAL and "legacy `lib/date`"; it now points at
  `CLAUDE.md` §2, which lists them.

## `.claude/agents/perf.md`

- "Parallelism was a founding goal and went unmeasured until after increment
  20b; your job is that it is never unmeasured again."
- The serial-phase profile, 2026-09-27
  (`docs/benchmarks/2026-09-27/serial-profile/README.md`): the serial part is
  about a third of single-thread refine, mostly Lawson legalisation, and the
  scan itself stops speeding up near 5× from load imbalance.
- The scaling ceiling: refine speeds up at most about 2.0-2.2× from 1 to 20
  threads, flat from about 7-8. The 2026-09-26 sweep
  (`docs/benchmarks/2026-09-26/scaling/`, medians per thread count) gives 2.2×
  on AC and 2.1× on battery; 2026-09-27 gives 2.02-2.03× on battery.
- Power state: "Ola develops while travelling". The 2026-09-26 sweep used 5
  runs.
- Plants: five crash dialogs showed on Ola's screen on 2026-09-27 from one 21c
  plant (`docs/benchmarks/2026-09-27/21c/data/a0eo/plant_eo_stalecommit/outcome.txt`).
- Evidence: the 2026-09-26 evidence sat only in `/tmp` until it was moved.

## `.claude/agents/reviewer.md`

- Sections 2-4 (readability, typing, style) were removed on 2026-09-27 because
  the gates cover them.
- CI: "a branch was reviewed twice, passed every local gate both times, and
  merged with CI failing on every commit — the workflow still referenced a
  build system the branch had deleted, and no review pass had looked at
  `.github/` at all."
- The three checks went unchecked across four merged PRs before they were
  written down.
- Scaffolding: "Three such comments survived three merges in
  `tests/cpp/CMakeLists.txt`, each describing headers that by then existed."
- Prose claims: "Eleven false or stale statements accumulated across six files:
  a README advertising a program that no longer exists, three documents naming
  a property-testing framework that four test files explicitly decline to use,
  two naming three ctest targets when there are sixteen."
- LOC: "Increment 3's design stated 'no split' *and* specified the seam to use if
  the implementation overran, naming the likely cause; the overrun happened in
  exactly that place, and nobody re-measured, so the contingency never fired."

## `.claude/agents/tester.md`

- Coverage: "the Python surface it covers is 37 lines against ~1,200 of
  unmeasured C++" (a figure from when it was written).
- Section B, generic security hardening, was removed on 2026-09-27.
- §3D is rule 3 of `2026-09-27-increments-14-to-20b.md`. "The Lawson bug fixed
  in `c23583b` was present from 14b to 20b because no refinement suite carried
  the first oracle."

## `docs/increments/README.md`

- **Why on disk.** "During increment 1 the architect's design was settled in a
  conversation, that conversation was lost, and the next round had to restate
  the entire design from scratch — twice. Worse, one restatement was wrong: it
  claimed `testing.md` already carried a `[planned]` gate that had only ever
  been *recommended*, and `@developer` caught it by grepping the history rather
  than by trusting the brief."
- **Step 1, Literature.** This is rule 1 of
  `2026-09-27-increments-14-to-20b.md`: "increments 14, 14b and 18 rebuilt
  Garland and Heckbert 1995 without reading it. It applies to increment files
  written after 2026-09-27."
- **Step 1, Legacy.** "the post-CGAL increments replace what CGAL *did*, and
  `legacy/triangulate_dem.h` only ever *called* it", which is why "nothing"
  was the commonest answer. "Increment files 01-05 predate this half."
- **Amendments.** "increment 3 pinned three behaviours that way, and increment
  2 retuned three constants *after* implementation when they turned out to
  depend on FMA contraction."
- **ROADMAP row.** "Increments 5b, 5c, 5d and 7 all shipped or were designed
  while that table said 5b was 'not started', and the whole `raster/` module is
  in the tree with no row and no record. The cause was structural: nothing in
  this protocol referred to that file, so the one document describing the
  state of the work was the one document no step maintained."
- **Merge commits.** "a governance audit found every production file in this
  repo had previously landed in the same commit as its test."
- **Acceptance.** This is rule 2 of `2026-09-27-increments-14-to-20b.md`. The
  power-state rule exists because Ola develops while travelling. "The
  2026-09-26 logs are not in `bench.py`'s format and are never a baseline."
- **Mutation testing** "is what has caught the real defects so far — a
  normalization inversion in 1a, and a coverage gap in 1b where 1436 of 1437
  assertions passed under a mutant."
- **Model tiering.** "the one attempt at tiering stalled for 600s and produced
  nothing."
- **Parallel suites.** "most are not independent — increment 2's `ring.hpp`
  returns `Box2` and `Segment2`, so its suite could not start before theirs."
- **Documentation defects.** "Increments 2 and 3 each carried a 'documentation
  debt this increment should clear' section. Between them they cleared
  nothing, grew from three entries to five, and listed 5 of the 11 defects an
  audit later found — reading as exhaustive while being less than half."

## `.claude/REQUIRED-READING.md`

- "The harness": "(auto mode has let unapproved pushes through)" becomes the
  risk stated plainly. The incident is recorded in
  `2026-09-27-required-reading-incidents.md`, "The harness: auto mode let
  pushes through".

## `testing.md`

- `noding`: "The list below replaces the five bullets that stood here until
  increment 5b. Three of them were false, and `docs/increments/05-noder.md`
  gives the measurement for each. A fourth — the feature-property bullet — was
  already corrected upstream by increment 7, from a one-bit `is_river` OR to a
  union over property sets; 5b keeps that correction and sharpens where the
  merge happens and what the oracle may be built from."
- The Hausdorff bullet: "The sum-of-lengths form that stood here was false
  twice over."

## Not moved, for the generic-layer design

`docs/PRINCIPLES.md` keeps history as structured fields rather than
narrative: each principle has an **Origin** pointer (commits, retrospectives)
and a **Status: Last exercised** field, and its retirement rule reads that
field. Stripping them would change how the file governs itself, so it is left
for `@architect`'s generic-layer design.
