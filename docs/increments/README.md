# Increment records

One file per increment of the post-CGAL core, written by `@architect` before
`@tester` writes anything, and read by every persona that works on that
increment. The directory also holds the occasional audit — a question about
work already shipped, answered with evidence and dated to the commit that
answered it, rather than a design for work not yet done. Such a file says so
in its own status line; `kernel-sufficiency-audit.md` is the first.

## Why these exist on disk rather than in a conversation

Agent transcripts do not survive a session. During increment 1 the architect's
design was settled in a conversation, that conversation was lost, and the next
round had to restate the entire design from scratch — twice. Worse, one
restatement was wrong: it claimed `testing.md` already carried a `[planned]`
gate that had only ever been *recommended*, and `@developer` caught it by
grepping the history rather than by trusting the brief.

A design that lives only in a prompt is a design that gets re-derived, and
re-derivation is where the errors enter. These files are the source of truth.
A prompt should say "read `docs/increments/02-core-geometry.md`", not contain
it.

The same argument applies to *what is currently being asked*, one level up from
what was designed; `.claude/REQUIRED-READING.md` rules on that, and on what a
cold session must do before it acts.

## The loop

Per `CLAUDE.md` §3, with the artifact each step produces:

1. `@architect` writes `docs/increments/NN-name.md`. Design only — types,
   invariants, exclusions, degeneracy policy, LOC estimate. No production code.
   The file carries a **Prior art: legacy and literature** section, written
   before the design, not after it.

   *Literature.* Name the published method the increment builds on, with a
   citation, and say what differs from it. If the increment claims something
   new, say what was searched and what was found; a novelty claim is not made
   without that check. This is the retrospective's rule 1
   (`docs/retrospectives/2026-09-27-increments-14-to-20b.md`): increments 14,
   14b and 18 rebuilt Garland and Heckbert 1995 without reading it. It applies
   to increment files written after 2026-09-27.

   *Legacy.* What the legacy tree
   holds on this increment's subject, and either what is being carried across or
   why nothing is. "Nothing" is a legitimate answer and the commonest one — the
   post-CGAL increments replace what CGAL *did*, and `legacy/triangulate_dem.h`
   only ever *called* it — but it is an answer, with the grep that reached it,
   not an omission — and the section pastes the command **with the file list it
   returned**, because a cited grep and an unrun one look identical on the page.
   The legacy tree no longer ships in the working tree; it is kept at the
   `legacy-archive` tag, so the grep runs against the tag
   (`git grep -n <pattern> legacy-archive -- legacy`).
   Increment files 01-05 predate this half.
   Where the answer is not "nothing", `@architect` reads the archived source
   from the `legacy-archive` tag and reports intent before `@tester` is spawned
   — before, because a suite written against re-derived intent pins the
   re-derivation, and a domain constant guessed wrong is then guarded by a test
   that agrees with it.
2. `@tester` reads it and writes a failing suite. No production code. The
   suite is committed **red**, before the implementation exists.
3. `@developer` reads both and makes it green. The green commit touches **no
   test file** — that is what makes the trace mean anything.
4. `@reviewer` audits before merge. CI is authoritative. An increment that
   touches refine or mesh code also needs `@perf`'s acceptance run (below).

Steps 2 and 3 are not strictly once each. A ruling can land after the red
commit, and the suite that encodes it is still `@tester`'s to write — increment
3 pinned three behaviours that way, and increment 2 retuned three constants
*after* implementation when they turned out to depend on FMA contraction. Such
amendments land as their own commit with the reason in the message, never
folded into the green one. The rule is not "tests are frozen after red"; it is
"`@developer` does not edit tests, and no test change hides inside an
implementation commit".

**The merge updates `ROADMAP.md`'s row for that increment**, in the same PR,
before the merge rather than after it. The table is an index of what shipped and
it is the only file a newcomer reads to orient themselves.

This is a step rather than an expectation because the expectation failed.
Increments 5b, 5c, 5d and 7 all shipped or were designed while that table said
5b was "not started", and the whole `raster/` module is in the tree with no row
and no record. The cause was structural: nothing in this protocol referred to
that file, so the one document describing the state of the work was the one
document no step maintained. A ledger nobody settles is worse than none, and the
same argument that retires a documentation-debt section applies to an index.

**Increment PRs merge with a merge commit, never a squash.** The whole protocol
rests on the red commit staying ahead of the green one in history, and one
squash destroys that evidence silently and irreversibly.

The red commit stays ahead of the green one in history. That trace is the only
thing that makes the test-first claim verifiable after the fact; a governance
audit found every production file in this repo had previously landed in the
same commit as its test.

## Acceptance: an increment that touches refine or mesh code

This is the one statement of the retrospective's rule 2. An increment whose
diff touches refine or mesh code (`include/terrain/refinement/`,
`include/terrain/mesh/`, and what drives them) is accepted only with:

- **the 1 m benchmark and a thread-scaling sweep**, run from a checked-in
  script, `tools/bench.py`, and compared with the previous increment's run;
- **the power state recorded with each run** (`pmset -g batt`). A battery run
  is compared only against a battery baseline, and an AC run only against an
  AC one, because Ola develops while travelling;
- the evidence committed under `docs/benchmarks/<date>/`.

`@perf` owns all of it (`.claude/agents/perf.md`); the tool's design is
`docs/benchmarks/bench-py.md`. When `bench.py run` finds no comparable stored
run (it prints `NO BASELINE` and exits 2), measure the previous increment's
merge commit with `--tree` back to back with the new one: the older run is the
newer one's baseline. The 2026-09-26 logs are not in `bench.py`'s format and
are never a baseline.

## Cost constraints

These exist because a round that costs three times what it needs to is a round
that gets skipped next time.

**Reference, do not restate.** Prompts point at the increment file. If a fact
is wrong there, fix the file rather than correcting it in a prompt — the
correction is otherwise lost with the transcript.

**Mutation testing is required only for the invariant-critical suite of an
increment**, named as such in its increment file. Writing a throwaway
implementation, compiling the suite against it, and killing deliberate mutants
is what has caught the real defects so far — a normalization inversion in 1a,
and a coverage gap in 1b where 1436 of 1437 assertions passed under a mutant.
It is also the slowest part of a round. Spend it where the topology decisions
are, not on value types whose failure mode is a typo.

**Template cross products are opt-in, not default.** Running every algorithm
as a `TEMPLATE_TEST_CASE` over every kernel × every model multiplies each
compile in the edit-validate loop. Use the full cross product only where an
instantiation proves something a later increment depends on — the increment
file says which.

**Match the model to the work — but the applicable surface is small.** Use a
smaller model only for a value-type header with no kernel parameter and no
exactness claim; otherwise do not spend the round setting it up. Predicate,
kernel and topology work is not delegable downward, and the one attempt at
tiering stalled for 600s and produced nothing.

**Independent suites run as parallel agents.** Two suites that do not share a
header do not need to share a round. Note that most are not independent —
increment 2's `ring.hpp` returns `Box2` and `Segment2`, so its suite could not
start before theirs.

**A documentation defect found during an increment is fixed in that increment's
PR, or it is not recorded.** Increments 2 and 3 each carried a "documentation
debt this increment should clear" section. Between them they cleared nothing,
grew from three entries to five, and listed 5 of the 11 defects an audit later
found — reading as exhaustive while being less than half. A ledger nobody
settles converts a one-line fix into a permanent entry and gives false
assurance that the rest is clean. Fix it, or leave it to be found.
