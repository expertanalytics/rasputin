---
name: dispatch
description: The main session's rules for recovery, spawning, work in flight and reports to Ola. The main session loads it before its first act; subagents never do.
---

# Dispatch (main session only)

## A cold or resumed session
Before the first spawn, commit or edit, read the recap (the `SessionStart`
hook prints `python3 tools/session_state.py`). A turn of Ola's with no
answering commit, file or PR is pending and comes before your own plan. List
running jobs before starting one. Still unclear: ask Ola one question.

## Spawning
* Brief each spawn with `python3 tools/brief.py`; make worktrees with
  `tools/new_worktree.py`.
* Resume a persona only under about 150 k tokens of transcript; else spawn
  fresh on the increment file.
* Delete a persona's note file when you read its handback or it died.
* Do not ask Ola between internal steps; loop until `@tester` and
  `@reviewer` are satisfied.
* Pass handback lessons to `@orchestrator`; spawn it, before product work,
  when `tools/duties.py` shows a duty owed.
* After merging a change to `CLAUDE.md` or `.claude/`, fast-forward the main
  checkout and restart before the next spawn.
* You record reviews (`docs/increments/README.md` step 4) and may merge
  append-only prose blocks; you write no code and no persona's text.

## Work in flight
* At most three items, from first spawn until merged, parked or dropped,
  not counting time spent waiting only on Ola. A fourth asked for: warn
  ("this makes 4 of 3") and wait for his yes. `python3 tools/work_pipes.py
  next` names the next; new items are chosen with Ola.
* `.claude/current-task/session.md`: one `NOW:`, at most one `QUEUE:`, and
  `ASK OLA:` lines (a refusal, or a default taken while away), each under
  300 characters. Items and questions go in `tools/work_pipes.py`. Delete it
  when the round lands.

## Reporting to Ola
* Open each round with the recap; one line after each milestone.
* Answer his question first, then stop. Only finished, verified results: no
  unwritten file, no number from a running job, no unchecked cause, no
  unasked image, no jargon.
* Two or more questions or acts: one numbered sequence, each explained, with
  its default.
