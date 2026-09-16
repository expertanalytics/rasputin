# Increment records

One file per increment of the post-CGAL core, written by `@architect` before
`@tester` writes anything, and read by every persona that works on that
increment.

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

## The loop

Per `CLAUDE.md` §3, with the artifact each step produces:

1. `@architect` writes `docs/increments/NN-name.md`. Design only — types,
   invariants, exclusions, degeneracy policy, LOC estimate. No production code.
2. `@tester` reads it and writes a failing suite. No production code. The
   suite is committed **red**, before the implementation exists.
3. `@developer` reads both and makes it green, touching no test.
4. `@reviewer` audits before merge. CI is authoritative.

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

**Match the model to the work.** Exact-predicate and topology reasoning earns
the larger model. Mechanical value types and their tests do not.

**Independent suites run as parallel agents.** Two suites that do not share a
header do not need to share a round.
