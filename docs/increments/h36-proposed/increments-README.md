# Increment records

One file per increment, `docs/increments/NN-name.md`, written by `@architect`
before `@tester` starts and read by every persona working on it.

## The loop

1. `@architect` writes the design: types, invariants, exclusions, degeneracy
   policy, LOC estimate; no production code. Under about 150 lines: one page
   (equations, sources, prediction, test list), its effect probed and
   reported to Ola before it is built. First comes **Prior art: legacy and
   literature**.

   *Literature.* The published method it builds on, cited, and each
   departure with its reason, saying so when one drops the method's
   guarantees. A novelty claim says what was searched and found; an unchecked
   one is a defect, as publication is kept possible.

   *Legacy.* `git grep -n <pattern> legacy-archive -- legacy`, the files it
   returned, and what is carried across or why nothing is; `@architect`
   reads that source and states its intent. `@reviewer` reviews the design
   before `@tester` starts, and again when an update adds a step.
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

## Small changes: the short path

Under 50 net lines (`CLAUDE.md` §2), no C++ file and nothing *Acceptance*
covers: the brief says short path, and `@reviewer` may send it back. The
increment file holds the change, its tests, what is left out, the estimate
and *Prior art*. One `@reviewer` round judges design, red and green together
and rules on the red step's pins; a pin that changes the design goes to
`@architect` before green. Red still comes before green.

## Acceptance: an increment that touches refine or mesh code

A diff touching `include/terrain/refinement/`, `include/terrain/mesh/` or
what drives them needs `@perf`'s acceptance run (`.claude/agents/perf.md`,
§1).

## Cost constraints

- **Reference, do not restate:** increment files are the source of truth;
  fix a wrong fact there, not in a prompt.
- **Mutation testing** only for a suite the increment file names
  invariant-critical, with its targets.
- **A documentation defect found in an increment is fixed in its PR, or not
  recorded.**
