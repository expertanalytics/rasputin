# Increment records

One file per increment of the post-CGAL core, written by `@architect` before
`@tester` writes anything, and read by every persona that works on that
increment. The directory also holds the occasional audit — a question about
work already shipped, answered with evidence and dated to the commit that
answered it, rather than a design for work not yet done. Such a file says so
in its own status line; `kernel-sufficiency-audit.md` is the first.

## Why these exist on disk rather than in a conversation

Agent transcripts do not survive a session. A design that lives only in a
prompt is a design that gets re-derived, and re-derivation is where the errors
enter. These files are the source of truth.
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
   citation, and say what differs from it. Design with the literature, not
   against it: every departure says why, and a departure that drops the
   method's guarantees is named as such. If the increment claims something
   new (for example an exact sup-norm guarantee on a DEM with constraints),
   say what was searched and what was found; a novelty claim is not made
   without that check. The option to publish is kept open, so an unchecked
   claim is a defect, not a detail.

   *Legacy.* What the legacy tree holds on this increment's subject, and either
   what is being carried across or why nothing is. "Nothing" is a legitimate
   answer, but it is an answer, with the grep that reached it, not an omission
   — and the section pastes the command **with the file list it returned**,
   because a cited grep and an unrun one look identical on the page.
   The legacy tree no longer ships in the working tree; it is kept at the
   `legacy-archive` tag, so the grep runs against the tag
   (`git grep -n <pattern> legacy-archive -- legacy`).
   Where the answer is not "nothing", `@architect` reads the archived source
   from the `legacy-archive` tag and reports intent before `@tester` is spawned
   — before, because a suite written against re-derived intent pins the
   re-derivation, and a domain constant guessed wrong is then guarded by a test
   that agrees with it.
2. `@tester` reads it and writes a failing suite. No production code. The
   suite is committed **red**, before the implementation exists.
3. `@developer` reads both and makes it green. The green commit touches **no
   test file** — that is what makes the trace mean anything.
4. `@reviewer` audits before merge (CI is authoritative, `CLAUDE.md` §4). An increment that
   touches refine or mesh code also needs `@perf`'s acceptance run (below).
   **The review leaves a trace in the increment file.** `@reviewer` is
   read-only, so its spawner copies the handback's verdict, the commit range
   it reviewed and its LOC count, verbatim, into a `## Review` section of
   `docs/increments/NN-name.md` (one entry per review round), and commits it
   before the push. Check: `grep -n '^## Review' docs/increments/NN-name.md` on the branch.

Steps 2 and 3 are not strictly once each. A ruling can land after the red
commit, and the suite that encodes it is still `@tester`'s to write; so can a
retuned constant, when implementation shows it depends on the platform. Such
amendments land as their own commit with the reason in the message, never
folded into the green one. The rule is not "tests are frozen after red"; it is
"`@developer` does not edit tests, and no test change hides inside an
implementation commit".

**The merge updates `ROADMAP.md`'s row for that increment**, in the same PR,
before the merge rather than after it. The table is an index of what shipped and
it is the only file a newcomer reads to orient themselves. It is a step rather
than an expectation because an index that no step maintains goes stale: a
ledger nobody settles is worse than none.

**Increment PRs merge with a merge commit, never a squash.** The whole protocol
rests on the red commit staying ahead of the green one in history — the only
thing that makes the test-first claim verifiable after the fact — and one
squash destroys that evidence silently and irreversibly.

## Acceptance: an increment that touches refine or mesh code

This is the one statement of this rule. An increment whose
diff touches refine or mesh code (`include/terrain/refinement/`,
`include/terrain/mesh/`, and what drives them) is accepted only with:

- **the 1 m benchmark and a thread-scaling sweep**, run from a checked-in
  script, `tools/bench.py`, and compared with the previous increment's run;
- **the power state recorded with each run** (`pmset -g batt`). A battery run
  is compared only against a battery baseline, and an AC run only against an
  AC one, because the power state changes timings;
- the evidence committed under `docs/benchmarks/<date>/`.

`@perf` owns all of it (`.claude/agents/perf.md`); the tool's design is
`docs/benchmarks/bench-py.md`. When `bench.py run` finds no comparable stored
run (it prints `NO BASELINE` and exits 2), measure the previous increment's
merge commit with `--tree` back to back with the new one: the older run is the
newer one's baseline. Only runs in `bench.py`'s format are baselines.

## Cost constraints

These exist because a round that costs three times what it needs to is a round
that gets skipped next time.

**Reference, do not restate.** Prompts point at the increment file. If a fact
is wrong there, fix the file rather than correcting it in a prompt — the
correction is otherwise lost with the transcript.

**Mutation testing is required only for the invariant-critical suite of an
increment**, named as such in its increment file. Writing a throwaway
implementation, compiling the suite against it, and killing deliberate mutants
catches defects a passing suite hides. It is also the slowest part of a round.
Spend it where the topology decisions are, not on value types whose failure
mode is a typo.

**Template cross products are opt-in, not default.** Running every algorithm
as a `TEMPLATE_TEST_CASE` over every kernel × every model multiplies each
compile in the edit-validate loop. Use the full cross product only where an
instantiation proves something a later increment depends on — the increment
file says which.

**Match the model to the work — but the applicable surface is small.** Use a
smaller model only for a value-type header with no kernel parameter and no
exactness claim; otherwise do not spend the round setting it up. Predicate,
kernel and topology work is not delegable downward.

**A documentation defect found during an increment is fixed in that increment's
PR, or it is not recorded.** A "documentation debt" section is a ledger nobody
settles: it converts a one-line fix into a permanent entry and gives false
assurance that the rest is clean. Fix it, or leave it to be found.
