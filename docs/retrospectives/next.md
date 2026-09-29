# Next retrospective: agenda

When: after auto-catchment is ready to use (Ola, 2026-09-28). Run by
`@orchestrator`. Items are added as they come up; the retrospective itself
gets its own dated file here, and this file is then emptied.

## Agents taking on each other's work (Ola, 2026-09-28): the main focus

Ola, 2026-09-29: "when we do our retrospective, we should have a specific
focus on the role bleed issue." Names for it: role drift or role bleed in
practice; "disobey role specification" (failure mode 1.2) in the MAST
taxonomy of multi-agent failures (Cemri et al. 2025, arXiv 2503.13657),
which traces most such failures to weak role definitions and missing checks.
Ola's analogy (2026-09-28): over-smoothing in GNNs, where repeated message
passing makes every node look alike. Here every hand-off carries the whole
context, and each persona picks up a bit of the others' jobs until the roles
blur. The GNN remedies map across: skip connections (restate each persona's
role and limits in every brief) and a bounded reach for messages (hard limits
on what each persona may touch). Inside transformers the same effect is
called rank collapse or over-smoothing, which skip connections also counter
(Dong, Cordonnier and Loukas 2021; from memory, not checked).

Ola: "One concern I have is that the agents 'leak' responsibilities to one
another." All personas are the same model with different briefs, and each
fills a gap it sees rather than handing it back. Seen on 2026-09-28:

- `@tester`, in the 16b-1/2 red step, wrote about 560 lines of throwaway
  production modules under `src_python/` to check its own tests.
- `@tester` wrote into the design document (its "Pinned by the red suite"
  section); the main session also edited it, to correct how Ola's rulings
  were attributed.
- `@architect` recorded two text corrections as Ola's rulings, following
  the wording of the main session's brief.

Causes: the boundaries are written in prose, not enforced (`@tester` can
write to `src_python/`); some briefs invite the drift ("show the test fails
against a wrong implementation"); an agent that finishes looks better than
one that hands back.

To decide:

1. **Hard limits per persona** (Ola: "add the hard limits to the
   retrospective list"). Enforced, not asked: e.g. a hook that refuses
   `@tester` writes under `src_python/` and `@developer` edits under
   `tests/`, and who may edit `docs/increments/`. Needs Ola's yes, since it
   touches `.claude/` hooks or settings.
2. Briefs that state what the agent must not do, and where to hand work
   back.
3. `@reviewer` audits who changed which files, and whether each change was
   that persona's to make.

Found by `@reviewer`'s first role-boundary audit (16b-1/2, 2026-09-28):
`@tester` commits e99c8ea and 3990449 added about 60 lines to the design
document (following the 16b-0 precedent); the main session's 413e91b edited
the design document too. `@developer` and `@architect` stayed in their
areas.

Second audit (16b-1/2 re-review): `@architect` fixed files outside
`docs/increments/` at `@reviewer`'s request (`project_structure.md`, a
benchmark README, another increment's record). No persona's remit names
those files; the hard limits need a stated owner for each area. Also: nine
test suites already on master still say, in the present tense, how they go
red; left for a follow-up.

2026-09-29: `@tester`, checking its query-plan test, temporarily edited
`src_python/tin_engine/io/geopackage.py` and restored it with
`git checkout` (reported, not committed). Brief-level limits did not stop
it; a hard limit would have.

2026-09-29: the main session's brief asked `@perf` to fix a `tools/bench.py`
bug without a red test. `@perf` declined, citing its own rule (a failing
test in `tests/python/test_bench.py` first). The brief crossed the line, the
persona held it. The fix (create the mesh's parent directory for a `--label`
with `/`) waits for a `@tester` red step.

2026-09-29, night: `@architect` built an uncommitted Python prototype of
the 22 design to check its area against NVE (304.91 against 305.54 km²).
Useful evidence, but it is implementation work in the design step.

2026-09-29, night: `@developer` changed the window rule during the 22 PR 1
green step (the design's rule never stopped early) without stopping to hand
back; the design amendment and the pinning test came after the code. The
fix was right and disclosed, but the order was code first.

2026-09-29, night: `@developer` added user-visible behaviour (lake refusals
reported under `--lakes`) in a review-fix commit with no failing test first;
`@reviewer` caught it and `@tester` pinned it afterwards.
