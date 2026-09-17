# Required reading

## Before you act in a cold or resumed session

"Ok, let's continue" is not an instruction. Neither is `git log`. Before the
first subagent spawn, the first commit, or the first edit, in this order:

1. `python tools/session_state.py` — prints `.claude/current-task.md` and the
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

- **The current ask lives in `.claude/current-task.md`** — untracked and
  gitignored, because it is per-working-tree and worthless after the fact.
  Write it when a user request starts a round, in three lines or fewer (the
  ask, the persona it went to, the file it will produce); delete it when the
  round lands. One file, overwritten, never a ledger.
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
persona, running builds, tests and gates, writing `.claude/current-task.md`,
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
that has repeatedly found real defects. On this branch alone `@reviewer` caught
the false limitation in `debd8d1`, the unswept citations in `760dbd9`, and two
of the three self-confirming invariants in increment 5. The session driving
those rounds produced every one of them and saw none.

## The harness

A `SessionStart` hook that runs `tools/session_state.py` and prints its output
is worth wiring, because this check fails exactly where a human is least likely
to run a command by hand. Propose it for approval separately; do not add it to
`.claude/settings.json` on your own initiative.

The same applies to a `PreToolUse` hook on `Bash(git push*)`, and it is the one
mechanism that would make the approval boundary self-enforcing. Do not assume
the permission system is that backstop: `git push` appears nowhere in
`.claude/settings.local.json`'s allow-list, yet two unapproved pushes went
through on 2026-09-17 under the auto mode configured in
`~/.claude/settings.json`. Until a hook is approved, the rule above is the
agent's own to keep, with nothing underneath it.
