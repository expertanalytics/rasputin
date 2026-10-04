# Next retrospective: agenda

When: after auto-catchment is ready to use (Ola, 2026-09-28). Run by
`@orchestrator`. Items are added as they come up; the retrospective itself
gets its own dated file here, and this file is then emptied.

## Taking recurring decisions off the main session (Ola, 2026-10-03)

Ola: "I need even more control. My gut feeling is to offload the recurring
decisions, or make them more mechanical." Evidence and plan:
`docs/research/2026-10-03-dispatcher-control.md`. It counts 28 main-session errors from
27 September to 3 October, 22 of them of a kind a script, template or hook
could have stopped, and proposes five stages, each for Ola to rule: (1)
briefs assembled from files by a script, with a hook that refuses a spawn
without one (about 140 lines); (2) one script that writes `session.md` and
one that says what step comes next on each branch (about 230 to 380); (3) a
checked-in merge runner with a lock (about 120); (4) hooks on the main
session's own turns: concurrency cap, idle guard, unexplained codes, rule
changes (about 125, probes first); (5) handbacks routed, reviews recorded,
`ROADMAP.md` status generated. Each stage deletes the prose rule it
replaces.

**Ruled 2026-10-04 (Ola):** "1: yes to stage 1, then stages 2 and 3", and
"A, B, C: yes": stage 4 runs its probes first, then Ola decides; the plan
moves to `docs/research/` (done); the recap's "Next on ROADMAP.md" rows are
capped at about 200 characters (first sentence plus status). What follows:
stage 1 is the next harness increment (h9), then stages 2 and 3, each
through the pipeline; all write governed paths (`tools/`, `.claude/`), so
by day. The row cap is a change to `tools/session_state.py` (governed).

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

Status: items 1-6 were ruled by Ola on 2026-10-03 (about 15:42-15:45 UTC);
each ruling is recorded under its item, and the proposal text is kept.
Items 7 and 8 were ruled on 2026-10-04, below. The rulings are to be implemented as one harness
increment, by day, through the pipeline: they touch governed files
(`tools/away.py`, `tools/session_state.py`, `.claude/REQUIRED-READING.md`),
which are not written in unattended mode.

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

   **Ruled 2026-10-03 (Ola, about 15:42-15:45 UTC):** yes, both
   parts: the rule in `REQUIRED-READING.md` (before Ola leaves, `session.md`
   names a fallback that needs no ruling and writes no governed path) and
   the longest-stretch-without-a-commit print in `away.py --back`. On the
   question to rule first: idle is accepted only when the fallback list is
   genuinely empty.

2. **`away.py --back` and the recap read only one checkout's
   `.claude/current-task/`.** Proposal: resolve the main checkout from the
   repository's common dir, and list `ASK OLA` lines from every worktree's
   folder. Cost: about 20 lines in `away.py`/`session_state.py`, with tests.

   **Ruled 2026-10-03 (Ola, about 15:42-15:45 UTC):** yes, as proposed.

3. **`ASK OLA` matching.** Proposal: keep the one-decision-per-line rule;
   count a line only if it starts (after a bullet) with `ASK OLA:`, and warn
   on an `ASK OLA:` line with nothing after the colon. Cost: about 10 lines
   in `session_state.py`, with tests.

   **Ruled 2026-10-03 (Ola, about 15:42-15:45 UTC):** yes, as proposed: a line
   counts only if it starts (after a bullet) with `ASK OLA:`, and an empty
   one warns.

4. **`session.md` was a 29-line log against a three-line rule** (measured
   2026-10-03; since rewritten to 4 lines). Rule needed: enforce the rule
   (the recap warns past three lines), or relax it to what the night queue
   needs, with a size the recap checks. Cost: a few lines either way.

   **Ruled 2026-10-03 (Ola, about 15:42-15:45 UTC):** option (a), enforce.
   `session.md` holds exactly one `NOW` line, one `QUEUE` line, and one
   `ASK OLA:` line per open decision, and nothing else: no rulings and no
   history, which go to the increment files, `ROADMAP.md` or this file.
   The recap warns on any other kind of line.

5. **Restart and resume.** The fix is a memory note; item 1 shows a note can
   fail within a day. Proposal: the cold-start steps in
   `REQUIRED-READING.md` gain "list the running background jobs before
   starting any", and the recap prints them. Cost: one rule line; the
   recap print is about 15 lines with tests.

   **Ruled 2026-10-03 (Ola, about 15:42-15:45 UTC):** yes: list the running
   background jobs before starting any, and the recap prints them.

6. **A merged rule change does not reach the running session.** Observed on
   h7: `@orchestrator` was spawned after the merge with its old persona
   prompt and the old `CLAUDE.md`. Proposal: after merging a change to
   `CLAUDE.md` or `.claude/agents/`, restart before spawning the changed
   persona; until then the brief says to read the persona file from disk.
   Cost: one rule line.

   **Ruled 2026-10-03 (Ola, about 15:42-15:45 UTC):** yes, as proposed.

7. **`ROADMAP.md` and harness increments.** No harness increment (h2 to h7)
   has a table row; h3 is named only in a prose list. Rule needed: do
   harness increments get roadmap rows? Cost: one rule line, and rows for
   the shipped ones if yes.

   **Ruled 2026-10-04 (Ola):** "On next, yes to all escept 15." No rows for
   harness increments; one pointer line under `ROADMAP.md` to
   `docs/increments/h*`. `@architect`, `ROADMAP.md`, not governed.
8. **`@reviewer` recorded its own rounds** (2026-10-03, reported by the
   main session as its own deviation). `docs/increments/README.md@390b516:59-63`:
   "`@reviewer` is read-only, so its spawner copies the handback's verdict
   ... into a `## Review` section ... and commits it". The main session
   briefed `@reviewer` to write and commit its own `## Review` entries, and
   seven runs complied: on 15e 3e01580, 78df916 and 70f8403 (rounds 1 and 2
   and the pre-push note), on release hardening 22b0106 and 11db227, and on
   this retrospective 39eb7b8 and e2230b1 (rounds 1 and 2, which also wrote
   into `docs/retrospectives/`, `@orchestrator`'s area). One refused: the
   reviewer of the citation-fix branch `worktree-agent-aa31dcdcb718cfc28`
   ("I did not record the round or commit", citing the rule above).
   Two causes: the rule sits in the increment protocol, which neither the
   brief nor `reviewer.md` restates; and "read-only" is enforced only by
   leaving `Write` and `Edit` out of `reviewer.md`'s `tools:` line, while
   `Bash` writes and commits freely. Proposal: the reviewer persona carries
   the rule, one line in `reviewer.md` §5 ("you do not edit or commit; your
   verdict goes in the handback, and your spawner records it"), because the
   persona is read on every run and a brief template is only as good as the
   brief that forgets it; the brief need not repeat it. Cost: one line in a
   governed file, by day.

   h6 would not have caught these, once it is implemented (its design
   merged in #145; no `ROLES` table is in `tools/` yet). Its `reviewer` row
   allows no writes (`docs/increments/h6-role-limits.md@390b516:70`), but
   §3.3 finds a Bash write by comparing the tree's dirty set
   (`git status --porcelain` plus `git hash-object`) before and after the
   call, not `HEAD`. All seven reviewers wrote and committed in one Bash
   call, so the paths were clean before and after, and h6 reports nothing;
   it sees a Bash write only while that write is uncommitted. §5 names the
   neighbouring gap ("a change made and undone within one Bash call") but
   not this one. Proposal for h6: the `pre-tool` snapshot also records
   `HEAD`, and `post-tool` attributes the paths of any new commit
   (`git diff --name-only <old>..<new>`) to the persona, judged against its
   row like a dirty path; and §5's list gains "a write committed in the
   same call" until then. Cost: about 10 lines and one test in h6.

   **Ruled 2026-10-04 (Ola):** yes (item 7's quote). One line in
   `reviewer.md`: you do not edit or commit; your spawner records the
   verdict. `@architect`, governed. The h6 `HEAD` proposal was not in the
   recommendation put to Ola and stays open.

## The merges of 2026-10-03: 15e, 15f-1, 15f-2, hardening, citations

Evidence: `2026-10-03-day-merges-15e-15f.md` (section numbers below are
that file's). Numbering continues from the section above. Each item needs
Ola's ruling; all but 15 touch governed files or tools.

Status: ruled by Ola on 2026-10-04: "On next, yes to all escept 15." Each
ruling is under its item.

9. **No `@orchestrator` check after five merges, until Ola asked** (2a). The
   session had loaded `CLAUDE.md` before #146 added "When to spawn
   `@orchestrator`", and all 33 of its subagents got that older copy too:
   item 6 again, with a cost this time. Proposals: (a) item 6's rule, with
   `/compact` as the cheap route: Claude Code re-reads the project
   `CLAUDE.md` from disk after `/compact` (code.claude.com/docs/en/memory,
   read 2026-10-03). Cost: one rule line. (b) Make the trigger mechanical:
   the recap prints the PRs merged since the newest commit under
   `docs/retrospectives/` ("3 merges with no `@orchestrator` check"). Cost:
   about 20 lines in `tools/session_state.py`, with tests.

   **Ruled 2026-10-04:** yes to (b): the recap prints the merges since the
   last `@orchestrator` check, about 20 lines in `tools/session_state.py`
   (governed), through the pipeline. Stage 2's `pipeline.py` would carry
   the same line; build it once.
10. **The required mutation tests were left out of 15e and 15f-2** (2b).
    The README requires them for a suite the increment file names
    invariant-critical; both files named one; the main session's briefs said
    "no mutation rounds (lean brief)", going past the lean-briefs note's own
    exception. To rule: does the README's requirement stand? If yes,
    `reviewer.md` §5 gains one line ("a named invariant-critical suite has
    its mutants run, with the kill record in a handback"), and the repair is
    a `@tester` task on master: run 15f-2's three named mutants (ES2, ES3,
    ES5) and 15e's on the store's gather/scatter. About 30 to 60 minutes, no
    production change, no governed path, so it is a good unattended
    fallback. If no, the README line becomes "optional". Cost: one line
    either way.

    **Ruled 2026-10-04:** the requirement stands for named
    invariant-critical suites. `@tester` runs 15e's and 15f-2's missing
    mutants on master as unattended filler (no governed path); `@architect`
    adds the line to `reviewer.md` §5 (governed).
11. **Absolute tolerances get a scale check** (3a, and ES16 in 3c). L12 and
    L14 used a fixed 1e-10 that stops working once lattice coordinates reach
    about 10⁶; ES16 used a fixed 1e-9 where the true bound is slope × offset.
    Proposal: one line in `architect.md` ("each absolute constant in a design
    states the scale it assumes and the largest input it was checked at, in
    a table like L16's") and the same line in `tester.md` for numeric bounds
    in assertions. Source: Dawson, "Comparing Floating Point Numbers, 2012
    Edition" (read 2026-10-03): a fixed epsilon fails once values grow;
    compare relative to magnitude or in ulps. Cost: two lines in governed
    files.

    **Ruled 2026-10-04:** yes. `@architect` adds one line each to
    `architect.md` and `tester.md`: every absolute constant states the scale
    it assumes (governed).
12. **Fused multiply-adds** (3b). `.claude/agents/developer.md@586fbc1:30` already makes
    `@developer` build with `-ffp-contract=off` too, and it caught ES13
    before CI. `tester.md` has no such line. Proposal: copy that line into
    `tester.md` (cost: one line), unless Ola rules the open question in
    `session.md` to turn contraction off project-wide (one CMake line;
    stored mesh hashes on the Mac change, and `@perf` would measure what
    arm64 loses without fused multiply-adds).

    **Ruled 2026-10-04:** dropped. Ola, same day: "turn off fused
    multiply-add" (increment 28), so contraction is off project-wide and the
    `tester.md` line is moot. Also ruled, for the record: the 15b sentence
    becomes "except the edge strip's points, which agree to rounding".
13. **A flagged test assumption went unconfirmed before green** (3c).
    `@tester`'s red handback flagged the ES15 choice about node (15, 8) as
    "The ruling did not say this". The main session passed it to Ola as
    information and asked no one to confirm it, and it was the assertion
    `@developer` later found wrong. ES16's absolute bound was not flagged.
    Proposal: when a red handback lists choices the ruling did not make,
    the spawner sends that list to `@architect` to confirm or correct before
    green, and `@tester` lists each such expected value with its derivation.
    This would have caught ES15. It would have caught ES16 only if its bound
    had been listed. Cost: one short `@architect` turn per red step that has
    such choices, and one line in the main session's dispatch rules. The
    alternative is no new rule: green caught both, and the fix was one
    commit.

    **Ruled 2026-10-04:** yes, as a line in h9's brief templates (stage 1),
    not a rule of its own.
14. **`ROADMAP.md`: conflicts, stale status, no owner** (3d, 2e, 2f). Two
    hand-resolved conflicts today, and master still says 15e is "awaiting
    the push and CI". Options: (a) the row's status says only designed, in
    progress or shipped, with the PR number, and the details stay in the
    increment file's status line; (b) a small tool builds the status column
    from the increment files' status lines, with a gate that the table is
    current, the pattern `towncrier` uses for changelogs (each change writes
    its own fragment, so nobody edits the shared file;
    towncrier.readthedocs.io, read 2026-10-03), about 60 lines with tests;
    (c) keep hand merges and only fix the stale text. Any of them also needs
    an owner for `ROADMAP.md` and for `.github/workflows/` (edited today by
    `@architect`, `@perf`, the main session and `@developer`); they belong on
    the hard-limits list above. Whatever the choice, master's two stale
    lines (`ROADMAP.md` row 15, `docs/increments/15e-memory-fixes.md@586fbc1:8`) need a docs fix.

    **Ruled 2026-10-04:** (a) now: the status column says designed, in
    progress or shipped with the PR number; `@architect`, `ROADMAP.md`, not
    governed, with the two stale lines. (b) later, in stage 5, owner
    `@architect`.
15. **Concurrency** (2d). Three agents ran at once several times, and three
    edited one worktree in parallel by design at 14:30, against the "agents
    in pairs" note. Nothing broke. To rule: do read-only `@reviewer` runs
    count toward the two, and may personas share a worktree when their files
    do not overlap? Cost: an update to the memory note.

    **Ruled 2026-10-04 (Ola):** "Read only agents should be allowed even
    when two agents are working, as long as they don't require the same
    files, ie that the files the read-only agent reads could change." Read-
    only agents do not count toward the two, provided nothing they read is
    being changed. The main session has updated its memory note; no file
    edit. Shared worktrees were not ruled.
16. **Numbers passed to Ola** (3e). 15f-1's overrun went out as 45 % (gross)
    instead of 36 % (net), copied from a handback. Proposal: `@developer`
    reports lines as §2 counts them, net, and the main session quotes line
    counts from `@reviewer`'s record only. Cost: one line in `developer.md`.

    **Ruled 2026-10-04:** the main session quotes line counts from
    `@reviewer` only. A dispatch rule: a line in h9's templates, or one in
    `CLAUDE.md` §3 (`@architect`, governed).
17. **A fallback's paths are checked by hand** (section 4). The first
    fallback queued today (h5's green step) writes three governed files;
    it was caught only because the morning check had just looked at h5.
    Proposal, to go with item 1 of "The window of 2026-10-02, the restart
    of 2026-10-03, h7" above (state the fill-the-window rule): a fallback
    line in `session.md` names the paths it will write, and the recap flags
    any governed one. Cost: about 15 lines in `tools/session_state.py`, with
    tests. Also: `session.md` was 10 lines at 15:26 UTC against the
    three-line rule (item 4 again).

    **Ruled 2026-10-04:** folded into stage 2's state tool, not a separate
    change.

## Recovery round of 2026-10-04: proposals waiting on Ola

Evidence and wording: `docs/retrospectives/2026-10-04-recovery-round.md`.
Rule proposals P1 (refusal tests assert the reason; whole-suite run against a
fix the design names), P2a (re-read the citation at-risk list after every
merge), P3a ("no checks reported" is not green; read `mergeStateStatus`), P5a
(a recorded publishing yes lapses with the session that heard it); tool
proposals P2b (`check_citations.py` labels at-risk lines unchanged, moved or
rewritten), P3b (the recap lists open PRs with merge state), P4 (`brief.py`
names the test command a worktree can run); one cut, C1 (about 90 words of
`REQUIRED-READING.md`). Questions Q1 to Q3 there. **Ruled 2026-10-04:** Q1 and Q3 yes (the whole-suite run in P1; P5a). Q2: no review-round cap, and `orchestrator.md`'s "past two rounds" line is dropped. The other proposals still wait on Ola.
