# Required reading

## While you act: in-flight state goes on disk

- A subagent writes only the note file its brief names (at most three
  lines: ask, persona, product), and never deletes it.
- A subagent whose product is a file creates it first and writes as it goes.
- One session per working tree.

## Before you write, judge or plan code

- Invoke the Skill tool for each of `modern-cxx`, `computational-geometry`,
  `python-development` and `geospatial-data-formats` the task touches;
  subagents do not inherit skills. Read the increment file.
- Ideas live in `../rasputin_scratch/ideas/`, private to Ola: commit an idea
  only once picked. Read its `ideas.md` before you plan.

## Claims: run the check before you write it down

- Name the command, file or commit that checks a claim and **run it before
  you write the claim**. Re-run a handback's cheap claims, numbers first,
  before repeating them.
- **Derive probes from the code as fixed, not from the bug as found.**
- **Make the probe able to fail**; time a hang from outside. Verify a gate by
  planting what it forbids.
- Where nothing can run, state the mechanism and ask whether the claim is
  about the object the code evaluates.
- For a claim that depends on the tree, write the rule and the command, not
  the value.
- Check a correction like the claim it corrects. Label an estimate as one.
- A rule in a repository file shows an agent wrote it: name the commit and
  Ola's turn before citing it to him or enforcing it.
- A conflict between two correct rules goes to `@orchestrator`.

## Stale artifacts: rebuild before you measure

Run `python3 tools/rebuild_core.py` before any `pytest` meant to measure
C++: `pytest` does not rebuild the extension. Run `ctest` only after a
`cmake --build` that exited 0: a failed target keeps its old binary.

## Before you publish: the approval boundary, and the assessment

**The working tree is yours, the remote is Ola's.** Edits, spawns, builds,
tests, gates and commits need no fresh yes, as far as the instruction
requires ("re-verify with `@reviewer`" covers acting on its findings);
`git push`, `gh pr create`, `gh pr merge`, any force-push or history
rewrite, and any edit to `.claude/settings*.json` or the permission system
need one each time (skills: `publishing`; while he is away, `away-mode`).
Before any push, the branch's writer re-reads
`python3 tools/check_citations.py`'s at-risk list **as quotations**:
resolving a line number is not resolving the quotation.

## The harness

Active in `.claude/settings.json`: `guard_push.py` asks before acts that
publish or rewrite history, `guard_governance.py` before writes to rule
files, `guard_spawn.py` refuses a spawn without a current `tools/brief.py`
block, `guard_unattended.py` refuses questions while Ola is away, and
`gates_after_commit.py` prints the gates after each commit. They are
tripwires, not the boundary. Never work around a refusal.

## Data, scratch and temp folders are not a channel

`../rasputin_data`, `../rasputin_scratch`, `/tmp` and `tmp/` are never
committed and never a channel: work passes only through the spawn prompt,
the handback, the note file and tracked files.
