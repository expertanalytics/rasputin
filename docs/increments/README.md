# Increment records

One file per full-lane increment, `docs/increments/NN-name.md`, read by
every persona working on it.

## The loop

1. `@architect` writes the design: types, invariants, exclusions, degeneracy
   policy, LOC estimate; no production code. Under about 150 lines: one page
   (equations, sources, prediction, test list), its effect probed and
   reported to Ola before it is built. First comes **Prior art: legacy and
   literature**.

   *Literature:* the published method it builds on, cited; each departure
   with its reason, named when it drops the method's guarantees; for a
   novelty claim, what was searched and found (an unchecked one is a defect:
   publication is kept possible).

   *Legacy:* `git grep -n <pattern> legacy-archive -- legacy`, its files,
   and what is carried across (its intent stated) or why nothing is.
   `@reviewer` reviews the design once, after Ola's answers and before
   `@tester` starts, naming every hole; code review checks them fixed; only
   a design wrong at its core is reviewed again.
2. `@tester` commits a failing suite, **red**, before the implementation.
3. `@developer` makes it green; the green commit touches **no test file**. A
   ruling or retuned constant after red is `@tester`'s amendment, in its own
   commit with the reason.
4. `@reviewer` audits before the push (verdicts: `reviewer.md`). The main
   session copies each round's verdict line, commit range and LOC count into
   the increment file's `## Review` section, applies any wordings of
   `APPROVED WITH WORDINGS` exactly in that commit, and commits before the
   push; that commit needs no second round.

The main session writes the as-built record (pins, kill record, counts) from
the handbacks and updates `ROADMAP.md`'s row in the PR. No stacked branches;
one increment file per PR.

## Two lanes

**The full lane** is the loop above, for the meshing core, the code whose
mistakes give wrong numbers without a crash, and the land-surface physics.

**The fast lane** is everything else: one `@developer`
(`tools/brief.py --lane fast`) writes the test and the fix on one branch,
CI runs, and it merges; no increment file, `@reviewer` or retrospective. CI
refuses a full-lane path without an increment file (`tools/check_lane.py`
lists the paths). After three red CI runs, stop and ask Ola through
`tools/questions.py`.

## Acceptance

A diff touching `include/terrain/refinement/`, `include/terrain/mesh/` or
what drives them needs `@perf`'s acceptance run (`perf.md` §1).

## Cost constraints

- **Reference, do not restate:** increment files are the source of truth;
  fix a wrong fact there, not in a prompt.
- **Mutation testing** only for a suite the increment file names
  invariant-critical, with its targets.
- **A documentation defect found in an increment is fixed in its PR, or not
  recorded.**
