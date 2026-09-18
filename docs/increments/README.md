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
   The file carries a **Prior art in `legacy/`** section: what the legacy tree
   holds on this increment's subject, and either what is being carried across or
   why nothing is. "Nothing" is a legitimate answer and the commonest one — the
   post-CGAL increments replace what CGAL *did*, and `legacy/triangulate_dem.h`
   only ever *called* it — but it is an answer, with the grep that reached it,
   not an omission — and the section pastes the command **with the file list it
   returned**, because a cited grep and an unrun one look identical on the page.
   Applies to increment files written from this commit onward; 01-05 predate it.
   Where the answer is not "nothing", `@migration-expert` reads
   the legacy source and reports intent before `@tester` is spawned — before,
   because a suite written against re-derived intent pins the re-derivation, and
   a domain constant guessed wrong is then guarded by a test that agrees with it.
2. `@tester` reads it and writes a failing suite. No production code. The
   suite is committed **red**, before the implementation exists.
3. `@developer` reads both and makes it green. The green commit touches **no
   test file** — that is what makes the trace mean anything.
4. `@reviewer` audits before merge. CI is authoritative.

Steps 2 and 3 are not strictly once each. A ruling can land after the red
commit, and the suite that encodes it is still `@tester`'s to write — increment
3 pinned three behaviours that way, and increment 2 retuned three constants
*after* implementation when they turned out to depend on FMA contraction. Such
amendments land as their own commit with the reason in the message, never
folded into the green one. The rule is not "tests are frozen after red"; it is
"`@developer` does not edit tests, and no test change hides inside an
implementation commit".

**Increment PRs merge with a merge commit, never a squash.** The whole protocol
rests on the red commit staying ahead of the green one in history, and one
squash destroys that evidence silently and irreversibly.

The red commit stays ahead of the green one in history. That trace is the only
thing that makes the test-first claim verifiable after the fact; a governance
audit found every production file in this repo had previously landed in the
same commit as its test.

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
