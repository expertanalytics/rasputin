# Increment records

One file per increment of the post-CGAL core, written by `@architect` before
`@tester` writes anything, and read by every persona that works on that
increment. An audit of shipped work may live here; its status line says so.

## The loop

Per `CLAUDE.md` §3, with the artifact each step produces:

1. `@architect` writes `docs/increments/NN-name.md`. Design only — types,
   invariants, exclusions, degeneracy policy, LOC estimate. No production code.
   The file carries a **Prior art: legacy and literature** section, written
   before the design, not after it.

   *Literature.* Name the published method the increment builds on, with a
   citation, and say what differs from it. Design with the literature, not
   against it: every departure says why, and a departure that drops the
   method's guarantees is named as such. If the increment claims something
   new, say what was searched and what was found; a novelty claim is not made
   without that check. The option to publish is kept open, so an unchecked
   claim is a defect, not a detail.

   *Legacy.* What the `legacy-archive` tag holds on the subject, with the
   command and the files it returned
   (`git grep -n <pattern> legacy-archive -- legacy`), and what is carried
   across or why nothing is. Where something is, `@architect` reads that source
   and states its intent before `@tester` starts.
   `@reviewer` reviews the design before `@tester` starts (on the short path,
   below, with the code), and again when a design update adds a step.
2. `@tester` reads it and writes a failing suite. No production code. The
   suite is committed **red**, before the implementation exists.
3. `@developer` reads both and makes it green. The green commit touches **no
   test file** — that is what makes the trace mean anything.
4. `@reviewer` audits before merge (CI is authoritative, `CLAUDE.md` §4). An increment that
   touches refine or mesh code also needs `@perf`'s acceptance run (below).
   **The review leaves a trace in the increment file.** `@reviewer` is
   read-only, so the main session itself, not an agent it spawns, copies the
   handback's verdict, the commit range
   it reviewed and its LOC count, verbatim, into a `## Review` section of
   `docs/increments/NN-name.md` (one entry per review round), and commits it
   before the push. Check: `grep -n '^## Review' docs/increments/NN-name.md` on the branch.
   `@reviewer`'s verdict is `APPROVED`, `APPROVED WITH WORDINGS` or `CHANGES REQUESTED`, and its handback ends with exactly one line `Review: <verdict> <sha>`,
   and no other line starts with `Review:`: its
   verdict and the full sha of the head it reviewed; a hook records it, and
   the away-mode push grant reads only that record. When `@reviewer`'s only
   findings are wordings in prose files (Markdown at the root or under
   `docs/`) that it names line by line with the new text, its verdict is
   `APPROVED WITH WORDINGS`. The main session applies them exactly, in the
   recording commit, and a branch whose only change after the review is that
   commit counts as reviewed: no second round. A wording in a code file, or
   any other finding, is `CHANGES REQUESTED`.

The main session writes the as-built record (pins, kill record, counts) into
the increment file from `@tester`'s and `@developer`'s handbacks or the
commits they name; no separate
`@architect` step.

A ruling or a platform-retuned constant after the red commit is `@tester`'s
amendment, in its own commit with the reason. `@developer` never edits tests,
and no test change hides in an implementation commit.

**The merge updates `ROADMAP.md`'s row for that increment**, in the same PR,
before the merge rather than after it.

**Increment PRs merge with a merge commit, never a squash**, so the red commit
stays ahead of the green one in history.

## Small changes: the short path

An increment estimated under 50 net lines (`CLAUDE.md` §2) that touches no
C++ file and nothing the *Acceptance* section covers takes the short path;
the brief says so, and `@reviewer` may send it back to the full path.

- Its increment file is short: what changes, the tests that pin it, what is
  left out, the estimate, and the *Prior art* section (step 1).
- One `@reviewer` round judges the design, the red step and the green step
  together, and rules on the red step's pins. A pin that changes the design
  still goes to `@architect` before green.
- The main session commits the verdict, the Status line and the ROADMAP
  row together, before the push.

Never dropped: the red commit ahead of the green one, a green commit with no
test file, and the approval boundary (`.claude/REQUIRED-READING.md`).

## Acceptance: an increment that touches refine or mesh code

This is the one statement of this rule. An increment whose
diff touches refine or mesh code (`include/terrain/refinement/`,
`include/terrain/mesh/`, and what drives them) is accepted only with:

- **the 1 m benchmark and a thread-scaling sweep**, run from a checked-in
  script, `tools/bench.py`, and compared with the previous increment's run;
- **the power state recorded with each run** (`pmset -g batt` on macOS;
  `/sys/class/power_supply` on Linux). A battery run
  is compared only against a battery baseline, and an AC run only against an
  AC one, because the power state changes timings;
- the evidence committed under `docs/benchmarks/<date>/`.

`@perf` owns all of it (`.claude/agents/perf.md`); the tool's design is
`docs/benchmarks/bench-py.md`. When `bench.py run` finds no comparable stored
run (it prints `NO BASELINE` and exits 2), measure the previous increment's
merge commit with `--tree` back to back with the new one: the older run is the
newer one's baseline. Only runs in `bench.py`'s format are baselines.

## Cost constraints

**Reference, do not restate.** The increment files are the source of truth,
and `.claude/REQUIRED-READING.md` rules on what is currently being asked.
Prompts point at them. If a fact is wrong there, fix the file rather than
correcting it in a prompt — the correction is otherwise lost with the
transcript.

**Mutation testing is required only for the invariant-critical suite of an
increment**, named as such in its increment file with the mutation targets its
kill record must cover.

**A documentation defect found during an increment is fixed in that increment's
PR, or it is not recorded.** No debt ledger.
