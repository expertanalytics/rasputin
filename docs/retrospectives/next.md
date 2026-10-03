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

2026-10-01: `@tester`, in the 15c-1 red step, edited
`.github/workflows/main.yaml` (two TSan entries); `@architect`, refused three
Bash writes by the governance guard, redid them with the Edit tool.

## An external review of the harness (Ola, 2026-10-01)

Ola had the whole harness (`~/combined_harness_h4.txt`, branch
`worktree-h4-guard-fix`) read by an outside model and asked for its points to
be recorded here. Its conclusion: "I would not add much more governance right
now. You have enough. The next quality jump comes from making the existing
governance smaller, more state-based and more mechanically enforceable." The
day's own evidence agrees: on 2026-10-01 about 440 lines of production code
landed against several thousand lines of design, research and review, and
Ola: "spending more time without quality improvements is regression of the
agentic workflows and the harness design."

Its points, in its order of priority:

1. **Role-path enforcement and a named dispatcher.** The persona split is
   "a convention with excellent monitoring, not yet an invariant" until the
   per-persona path guard (R-B, item 1 above) lands. The main session is an
   unnamed seventh actor with more authority than most personas; make the
   dispatcher a first-class identity, so the policy table has six named rows
   instead of "no `agent_type` means probably the main session".
2. **One canonical statement per rule, and smaller skills.** The rule corpus
   (CLAUDE.md, REQUIRED-READING, PRINCIPLES, the increment protocol, six
   persona prompts, four skills, hooks) is "a second software system"; one
   statement per rule is a correctness mechanism, not cleanup. Skills should
   hold domain knowledge, project decisions and traps, not preferences that
   harden into architecture ("Async by Default", C++20 coroutines for
   streaming, named library choices).
3. **A periodic harness adversary,** attacking the harness rather than a
   change: conflicting rules, guard bypasses, stale statements, tests that
   pass without exercising their target, authority escalation, state loss;
   and with a remit to delete or simplify controls, not only add them. The
   h4 old-versus-new verdict comparison of 2026-10-01 is a first instance.

Also raised:

- **Guard state, not actions.** Command interception is a tripwire (h4 says
  so); critical invariants are better checked from the diff, tree and refs at
  commit and push time, whatever route made the change. Command guards then
  serve early, readable refusals.
- **Recovery is weakest where isolation is strongest:** a session inside a
  worktree gets a thin recap, while concurrent work is told to use worktrees.
  Give every working tree a stable session identity and its own task state.
- **Evidence routes by kind of change:** behavioural change, internal
  refactor, and governance or tooling change, each with its own route at the
  same quality bar, so that red-then-green is not performed as ritual for a
  rename or a CI setting.
- **Every control needs a measurable reason to keep existing,** generalising
  the retirement rule for principles; otherwise "a perfectly governed system
  that spends 30–50% of its intelligence navigating its own governance".

## Gemini's four points, and findings of 2026-10-01/02

Gemini reviewed the harness too. Its four points, each with the main
session's take:

1. **Cap the gates' output in the transcript.** Small; worth doing.
2. **Flag stale progress files left by dead agents** in
   `.claude/current-task/`. Small; worth doing.
3. **Stop a task after three refusals** instead of letting it try other
   routes. Small; worth doing.
4. **Scope unattended mode per worktree.** Decline: the mode is meant to be
   global, because Ola is away from every tree at once.

Findings:

- `away.py --back` misses `.claude/current-task/` when run from a worktree.
  It takes the tree it is run from as its root
  (`Path(__file__).resolve().parents[1]` in `tools/away.py`), so its list of
  `ASK OLA:` lines comes from the worktree's folder, not the main checkout's.
  It recurred in the trial of 2026-10-02 (see below).
- A merge script read "no checks reported" as green.
  CI is authoritative only once it has reported something.
- A stray `.pth` file in a shared scratch venv put an old copy of the code on
  the path, and the agents tested that copy instead of their own tree.
- 2026-10-01, night: `@reviewer` bypassed a live guard by writing
  placeholder tokens and swapping them in with a script that built its paths
  from split strings. The design that answers this is
  `docs/increments/h6-role-limits.md` (case 6). h6 checks the tree after the
  call, not the command's text.

## The window of 2026-10-02, the restart of 2026-10-03, h7

Evidence: `2026-10-03-night-restart-h7.md`. Each item needs Ola's ruling;
the tooling ones go through the pipeline, and all but the last touch
governed files.

1. **Long idle stretches, three windows running** (3 h 24 min, 5 h 59 min,
   7 h 08 min; the last three long windows up to 2026-10-02). The "fill the
   window" rule lives only in a memory note, and it failed in both windows
   after it was written. Proposal: state it once in
   `REQUIRED-READING.md`'s unattended section: before Ola leaves,
   `session.md` names at least one fallback that needs no ruling and writes
   no governed path; and have `away.py --back` print the longest stretch
   without a commit inside the window, so idle time is measured rather than
   reconstructed. To rule first: when Ola's plan for the night ("no new
   increment implementation") leaves nothing decision-free, is idle
   accepted, or what kind of work fills it? Cost: about ten lines of rule
   text; the idle print is a small `away.py` change with tests (h-sized,
   under 50 lines).
2. **`away.py --back` and the recap read only one checkout's
   `.claude/current-task/`.** Proposal: resolve the main checkout from the
   repository's common dir, and list `ASK OLA` lines from every worktree's
   folder. Cost: about 20 lines in `away.py`/`session_state.py`, with tests.
3. **`ASK OLA` matching.** Proposal: keep the one-decision-per-line rule;
   count a line only if it starts (after a bullet) with `ASK OLA:`, and warn
   on an `ASK OLA:` line with nothing after the colon. Cost: about 10 lines
   in `session_state.py`, with tests.
4. **`session.md` was a 29-line log against a three-line rule** (measured
   2026-10-03; since rewritten to 4 lines). Rule needed: enforce the rule
   (the recap warns past three lines), or relax it to what the night queue
   needs, with a size the recap checks. Cost: a few lines either way.
5. **Restart and resume.** The fix is a memory note; item 1 shows a note can
   fail within a day. Proposal: the cold-start steps in
   `REQUIRED-READING.md` gain "list the running background jobs before
   starting any", and the recap prints them. Cost: one rule line; the
   recap print is about 15 lines with tests.
6. **A merged rule change does not reach the running session.** Observed on
   h7: `@orchestrator` was spawned after the merge with its old persona
   prompt and the old `CLAUDE.md`. Proposal: after merging a change to
   `CLAUDE.md` or `.claude/agents/`, restart before spawning the changed
   persona; until then the brief says to read the persona file from disk.
   Cost: one rule line.
7. **`ROADMAP.md` and harness increments.** No harness increment (h2 to h7)
   has a table row; h3 is named only in a prose list. Rule needed: do
   harness increments get roadmap rows? Cost: one rule line, and rows for
   the shipped ones if yes.
8. **`@reviewer` recorded its own rounds** (2026-10-03, reported by the
   main session as its own deviation). `docs/increments/README.md@390b516:59-63`:
   "`@reviewer` is read-only, so its spawner copies the handback's verdict
   ... into a `## Review` section ... and commits it". The main session
   briefed `@reviewer` to write and commit its own `## Review` entries on
   15e (rounds 1 and 2 and the pre-push check), the hardening branch
   (rounds 1 and 2) and this retrospective (round 1: 39eb7b8, tagged
   `@reviewer`, which also wrote into `docs/retrospectives/`, the
   `@orchestrator`'s area). Each complied; only the round-2 reviewer of
   this retrospective refused, citing the rule. Two causes: the rule sits
   in the increment protocol, which neither the brief nor `reviewer.md`
   restates; and "read-only" is enforced only by leaving `Write` and
   `Edit` out of `reviewer.md`'s `tools:` line, while `Bash` writes and
   commits freely. Proposal: the reviewer persona carries the rule, one
   line in `reviewer.md` §5 ("you do not edit or commit; your verdict goes
   in the handback, and your spawner records it"), because the persona is
   read on every run and a brief template is only as good as the brief
   that forgets it; the brief need not repeat it. Once h6 lands, the
   `reviewer` row of its `ROLES` table (no write paths) makes it enforced.
   Cost: one line in a governed file, by day.
