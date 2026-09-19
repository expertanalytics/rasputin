# Required reading

## Before you act in a cold or resumed session

"Ok, let's continue" is not an instruction. Neither is `git log`. Before the
first subagent spawn, the first commit, or the first edit, in this order:

1. `python3 tools/session_state.py` — prints `.claude/current-task/` (the
   session's ask first, each subagent's after it, as context) and the
   last human turns of the predecessor session, **including prompts the user
   queued and the harness absorbed mid-turn**. An absorbed prompt never appears
   as a normal turn, which is how a session on 2026-09-17 lost the instruction
   "can we verify that we haven't built a kernel that is insufficient for our
   purposes" and spent a full round designing an unrequested increment instead.
2. Judge each surfaced turn against the tree: a turn with no answering commit,
   file or PR is still pending, and an earlier pending turn outranks your
   reconstruction of "what comes next".
3. If anything is still ambiguous after both, **ask the user one question**. A
   question costs a line; a guessed round costs ~150k tokens.

Recovery is sufficient when every surfaced turn is either answered on disk or
explicitly cancelled by the user. `git log` plus `docs/increments/` is not
sufficient: it shows what was finished, never what was asked.

## While you act: in-flight state goes on disk

`docs/increments/README.md` puts designs on disk because transcripts do not
survive. The same applies one level up, to *what is currently being asked*:

- **The current ask lives in `.claude/current-task/`** — a directory, untracked
  and gitignored, because it is per-working-tree and worthless after the fact.
  Three lines or fewer per file (the ask, the persona it went to, the file it
  will produce). Work here runs in parallel, so there is one file per writer:

  - `session.md` is the **main session's** record, and only the main session
    writes it. It is what a cold session reads first.
  - every other file is one subagent's, named `<persona>-<HHMMSS>.md`
    (`tester-142530.md`). **The spawner assigns the path in the prompt**; a
    subagent writes that path and no other. Assigning it is what keeps two
    parallel `@tester`s from colliding, and it is why "do not touch
    `session.md`" no longer has to be retyped into every brief.

  This became one file per writer after a single shared file failed within an
  hour of being legislated: `@tester` wrote `.claude/current-task.md` and
  destroyed the session's record verbatim, leaving three lines of subagent ask
  where the only on-disk note of the round had been. A mechanism built to
  survive a context loss would have caused one.

  It is still not a ledger, and several files make a ledger easier to get wrong,
  so the deletions are assigned rather than left to good intentions:
  `session.md` is overwritten in place and deleted when the round lands; a
  subagent's file is deleted by **its spawner**, on reading the handback — never
  by the subagent, which may be exactly the thing that died, and the rule
  recurses, so a nested subagent's file is its own spawner's and never the
  session's; and a session **sweeps `.claude/current-task/` of every file
  except those of agents it is currently waiting on, directly or transitively,
  after step 1 above and never before it.** A file from a dead agent therefore
  survives at most one round. That last sweep is the backstop, because it is the
  only deletion that still happens when the agent owing one is gone.

  The sweep test is an exception for what is live, and the two simpler spellings
  both failed. "Every file I did not just spawn" is provenance: a session already
  mid-round with two live subagents has, at the instant it starts a parallel
  round, *just* spawned none of them, so that criterion deletes its own live
  agents' files. "Every file I am no longer waiting on" is its dual and fails the
  other way — an orphan from a dead predecessor session is one you were never
  waiting on, so it is never one you are *no longer* waiting on, and nothing ever
  licenses deleting it. *Transitively* is what settles the nested case: you wait
  on your subagent, not on its subagent, so without it a live nested file is
  swept and a dead one lives forever, and the text supports whichever reading
  the reader arrives with.

  Sweep after the read, never before, because the file a premature sweep
  destroys is precisely the dead step's — the one the backstop exists for.

  **One session per working tree.** `session.md` is a reserved singular name and
  nothing assigns it, so two sessions in one tree collide on it exactly as
  subagents once did. The subagent fix does not apply, because the thing that
  hands out a subagent's path is its spawner and no one hands out a session's.
  Concurrent sessions need separate worktrees; `git worktree add` is the cheap
  answer and `.claude/current-task/` is per-tree by construction.

  Nothing mechanically stops a persona writing `session.md`; the control is that
  it is handed a different path and never has cause to guess one. A `PreToolUse`
  hook is the only real enforcement, and it needs approval — see *The harness*.
- **A subagent whose product is a file creates that file first and writes
  incrementally.** The kernel audit above died having emitted only a progress
  line; had it created its target file on arrival, the round would have been
  resumable instead of lost. A persona's prompt names the output path.

## Before you write, judge or plan code

Invoke the Skill tool for `modern-cxx`, `computational-geometry`,
`python-development` and `geospatial-data-formats` — whichever touch the task —
before writing, judging or planning any code.

This is stated as a step rather than left to the `skills:` frontmatter key
because that key does not reliably preload them and subagents do not inherit
skills from the caller. Both are re-testable in about a minute, which is why the
instruction is an imperative rather than a declaration.

A session that skips it re-derives decisions the project has already written
down. Each `SKILL.md` is the settled answer to a question that otherwise comes
back as an open one — `geospatial-data-formats/SKILL.md` rules on where GeoTIFF
decoding sits, `computational-geometry/SKILL.md` on degeneracy policy.

When you write a justification into this project's governance, name the command,
file or commit a sceptical reader would use to check it — and **run it before
you write the claim down**. Naming a check is not performing one. `debd8d1`'s
message stated that `tools/session_state.py` "cannot identify which transcript
is the current session"; it was reading `CLAUDE_SESSION_ID` where the harness
exports `CLAUDE_CODE_SESSION_ID`, so a one-word bug shipped as an inherent
limitation, from reading the source instead of running it. One second of
execution refutes it. `0fa05e3` is the correction.

The probe set is a claim too. Derive it from the code as fixed, not from the bug
as found: a fix that widens what the code accepts widens the inputs that can
break it, and re-running the old probes tests the old code. Measured — five
broken inputs verified a fix that had just widened a glob from `*.md` to every
file, all five were valid UTF-8 because the bug had been, and the first non-UTF-8
file killed the tool (`fdbd532`).

And make the probe able to fail. If a pass looks the same as a probe that never
ran, it measured nothing. A reviewer timed this tool's FIFO case with
`signal.alarm`, but `TimeoutError` is an `OSError`, so the blocking `open()`
raised the alarm *inside* the `except OSError` it was meant to expose: the probe
reported "no hang" in output produced entirely by the hang. Time a hang from
outside the process, and prefer a probe whose pass and its own absence look
different — `b0bf129` is the correction.

This tree has a standing instance of that hazard: **`pytest` does not rebuild
the C++ extension.** scikit-build-core's editable auto-rebuild needs `ninja`,
which is not installed, so it silently no-ops, and a run after a change to
`bindings/core.cpp` exercises the previously installed `.so`. A green handback
from a session that never rebuilt is byte-identical to one that did. Rebuild and
reinstall explicitly before every `pytest` that is meant to measure C++:

```bash
cmake --build build-pyext -j --target _core
cp build-pyext/_core.cpython-*-darwin.so .venv/lib/python3.*/site-packages/tin_engine/
```

`pip install -e . --no-build-isolation` is not the route — `scikit_build_core`
is absent from the venv. This paragraph was verified the way it asks you to
verify: increment 7's fixed `bindings/core.cpp` was replaced on disk with the
pre-fix version that four tests were written to refute, and `pytest` reported
94 passed. A broken X, and Y did not notice.

When there is nothing to run — a comment, a design invariant, a claim of the
form "X is verified by Y" — one question catches the same defect by inspection,
in a line, with no build: **is the claim about the same object the code
evaluates?** In every occurrence recorded below it was not. A §2 gate
matched path-shaped keys (`boost/geometry`, `date/date.h`) against bare CMake
tokens, so `find_package(Boost COMPONENTS geometry)` and `date::date` passed
until `1807e73`. Three increment-2 tests bounded the compiler's FMA choice
rather than `FastKernel`'s error (`dd67a68`). Increment 5's oracles used exact
incidence where the producer used hot-pixel proximity, a "second" candidate
generator that was the same broad phase, and the dedup's own records in place of
the input (`docs/increments/05-noder.md`, guarantees 14 and 15). A comment credited the
NaN fixture with killing the `!= 0.0` mutant of a function guarded by
`> 0.0 && <= max()`, where the second conjunct rejects NaN under either
spelling, so that fixture kills nothing (`7ece838`). That sentence was itself
committed with the two objects swapped and caught by this branch's own pre-push
pass — the check applied to the paragraph that defines it. Each repair was made in the same
mode as the defect — reasoning about the claim instead of running it — which is
how the repairs kept seeding the next occurrence.

One corollary, learned the expensive way on increment 7's review rounds, which
found this defect repeatedly and **found no defect in behaviour at all**. Several
instances were in the text repairing the previous one. The common property: each
was a sentence that had to be **re-derived whenever anything around it changed**
— a commit hash naming "the current anchor", a count of how many commits had
moved a figure, a list of which ones. Correcting such a sentence is a move that
reproduces itself, because the commit applying the correction is itself a commit
and can falsify the line it just wrote.

So when a claim's truth depends on the state of the tree, **write the rule and
the command that resolves it, not the resolved value**. `07-edge-properties.md`
now says "measured at the last commit that touches a production file" and prints
`git log --oneline -1 -- include/ bindings/ src_python/`, where it used to name a
hash; the enumeration beside it was deleted rather than corrected again. Those
were the first corrections in that sequence that did not falsify themselves on
landing. A resolved value is a citation with an expiry date nobody can see —
and note that no counter appears in this paragraph, which is the rule applied to
itself: "how many rounds" is exactly such a value.

So **do not write "X is verified by Y" until Y has been run against a broken
X**, and where Y cannot be run, apply the object-identity question instead.
Verify a gate the way this paragraph's first example was verified before being
written down: append the prohibited directive to `CMakeLists.txt`, run
`python3 tools/check_prohibited_deps.py`, confirm it fails naming both, restore
the file. Note the asymmetry the list exposes — `tools/*.py` is executable code
with no suite behind it, unlike the C++ it gates.

This binds claims about behaviour — what a script does, what a suite covers,
what a file says — not design opinions, which have no command to run. If you
cannot, state the mechanism rather than the incident: a mechanism someone can
re-run is stronger evidence than an event they can only be told about.

Also read `docs/increments/README.md` — the increment protocol, and the cost
constraints on a round — plus the increment file for whatever you are working
on. Designs live on disk precisely so they are not re-derived from a prompt.

## Before you publish: the approval boundary, and the assessment

One line settles both questions, because publishing is the only irreversible
thing an agent here does: **the working tree is yours, the remote is the
user's.**

### What runs on an instruction's momentum

No fresh approval, however many steps it takes: editing any file, spawning any
persona, running builds, tests and gates, writing under `.claude/current-task/`,
and `git commit`. All of it is undone by `git reset` and none of it is visible
to anyone else. An instruction carries through as much of this as satisfying it
requires — "re-verify with `@reviewer`" does authorise acting on what the
reviewer finds, including new edits and another persona.

A rule that made each of these ask would be the wrong answer. The autonomy is
the point.

### What needs a fresh yes, every time

`git push` — including a branch's first — `gh pr create`, `gh pr merge`, any
force-push or history rewrite, and any edit to `.claude/settings*.json` or the
permission system. Editing `CLAUDE.md`, `.claude/**` or `docs/` is ordinary
tree work; it becomes a publication at the same moment everything else does.

Push, not commit, is the line because pushing a branch that has an open PR
re-triggers CI, re-notifies a live reviewer, and silently changes what the user
is being asked to merge while they are holding the merge. A local commit does
none of that.

### The standing grant, and its limit

An instruction that names a publishing act is the yes for **one** occurrence of
that act. "Push it and re-verify with `@reviewer`" authorised one push. What
followed on 2026-09-17 was a reviewer's CHANGES REQUESTED, an `@architect`
spawn, edits, a probe, two commits and a **second** push to a branch targeting
master, with no further user input. Everything before the second push was
momentum. The second push published a different tree under a spent grant.

The test is mechanical: *have I already performed this named act once, and has
the tree changed since?* If yes, ask. An ask costs one line.

### The assessment fires at that same boundary

`CLAUDE.md` §3 step 4 reads as the last step of a code increment, and so a
branch whose product is prose or tooling has been entering no loop at all.
That is the gap. **`@reviewer` runs once on any branch before its first push,
whether or not the branch contains production code, and again before a push
that follows a round of findings.**

The trigger is the push boundary — not per-round, not per-commit, not on
everything. It is the same event as the approval boundary above, which is why
neither needs bookkeeping of its own: when you stop to ask, that is when the
pass has already happened.

On a prose or tooling branch the pass has a different scope from a code one:

- every claim the branch asserts, checked by running the thing it asserts
  about, per the rule above;
- `python3 tools/check_citations.py`, with its at-risk list **re-read as
  quotations**. Resolving a line number is not resolving the quotation — a
  sweep on this branch resolved all 46 citations in the kernel audit and was
  still wrong three times, because three pointed at text the same branch had
  rewritten to say the opposite;
- for any rule the branch adds, the incident it prevents, named with a commit.

One pass per branch is affordable on `docs/increments/README.md`'s cost terms —
it is the same pass the code path already pays — and it is the only step here
that has repeatedly found real defects. On PR #68 alone `@reviewer` caught
the false limitation in `debd8d1`, the unswept citations in `760dbd9`, and two
of the three self-confirming invariants in increment 5. The session driving
those rounds produced every one of them and saw none.

## The harness

A `SessionStart` hook that runs `tools/session_state.py` and prints its output
is worth wiring, because this check fails exactly where a human is least likely
to run a command by hand. Propose it for approval separately; do not add it to
`.claude/settings.json` on your own initiative.

A `PreToolUse` hook denying `Write`/`Edit` on `.claude/current-task/session.md`
from a subagent is the only thing that would make the session's file structurally
safe rather than conventionally safe; propose it the same way.

The same applies to a `PreToolUse` hook on `Bash(git push*)`, and it is the one
mechanism that would make the approval boundary self-enforcing. Do not assume
the permission system is that backstop: `git push` appears nowhere in
`.claude/settings.local.json`'s allow-list, yet two unapproved pushes went
through on 2026-09-17 under the auto mode configured in
`~/.claude/settings.json`. Until a hook is approved, the rule above is the
agent's own to keep, with nothing underneath it.
