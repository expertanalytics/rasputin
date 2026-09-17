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
file or commit a sceptical reader would use to check it. If you cannot, state
the mechanism rather than the incident: a mechanism someone can re-run is
stronger evidence than an event they can only be told about.

Also read `docs/increments/README.md` — the increment protocol, and the cost
constraints on a round — plus the increment file for whatever you are working
on. Designs live on disk precisely so they are not re-derived from a prompt.

## The harness

A `SessionStart` hook that runs `tools/session_state.py` and prints its output
is worth wiring, because this check fails exactly where a human is least likely
to run a command by hand. Propose it for approval separately; do not add it to
`.claude/settings.json` on your own initiative.
