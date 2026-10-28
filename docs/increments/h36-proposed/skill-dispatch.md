---
name: dispatch
description: Main-session rules for recovery, spawning, work in flight and reports to Ola. Subagents never load it.
---

# Dispatch

## A cold or resumed session
Before the first spawn, commit or edit, read the recap
(`tools/session_state.py`), which names the skills to load. A turn of Ola's
with no answering commit, file or PR is pending and comes before your own
plan. List running jobs before starting one. Still unclear: ask Ola one question.

## Spawning
* **First question.** Before an item's first spawn, answer in one line Ola
  may veto: is this needed; can we delete instead; which lane?
* Brief each spawn with `python3 tools/brief.py`; make worktrees with
  `tools/new_worktree.py`.
* Resume a persona only under about 150 k tokens of transcript; else spawn
  fresh on the increment file.
* Delete a persona's note file when you read its handback or it died.
* Do not ask Ola between internal steps; loop until `@tester` and
  `@reviewer` are satisfied.
* Pass handback lessons to `@orchestrator`; spawn it when `tools/duties.py`
  shows a duty owed; with no free item, that duty runs first.
* After merging a change to `CLAUDE.md` or `.claude/`, fast-forward the main
  checkout and restart before the next spawn.
* You record reviews (`docs/increments/README.md` step 4); with both
  writer slots busy you may merge append-only blocks on a prose branch. You
  write no code and no text in a persona's file.

## Work in flight
* At most three items in flight, from first spawn to merge, park or drop;
  time waiting only on Ola does not count. A fourth asked for: warn
  ("this makes 4 of 3") and wait for his yes. `python3 tools/work_pipes.py
  next` names the next; new items are chosen with Ola.
* Questions for Ola go through `tools/questions.py` (`--item` for an
  item's); the recap prints the open list. Items go in `tools/work_pipes.py`.

## Reporting to Ola
* Open each round with the recap; one line after each milestone.
* Answer his question first, then stop. Report only finished, verified
  results: no unwritten file, running job's number, unchecked cause,
  unasked image or jargon.
* Two or more questions or acts: one numbered sequence, each explained, with
  its default.
