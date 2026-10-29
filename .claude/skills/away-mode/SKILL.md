---
name: away-mode
description: Rules while Ola is away (unattended mode). Load when the recap says it is on, when he says he is leaving, or when a guard refuses an act.
---

# Unattended mode

* Only Ola starts it (`python3 tools/away.py <duration>`, his own terminal)
  and ends it (`--back`). No agent runs `away.py` or writes
  `<git-common-dir>/harness/`.
* When he says he is leaving, ask how long and ask him to run `away.py`.
  `tools/work_pipes.py` then holds a fallback needing no ruling or governed
  write, one still queued in the last hour; idle only when none exists.
* A refused act is already queued: never retry it or work around it,
  except a refusal naming only ordinary files (no rule file, harness state
  or push), redone with Edit or Write and reported as a
  `GUARD FALSE POSITIVE:` line naming them; and a push, PR or merge whose
  sole reason has gone (head now approved, branch now granted, form now
  right).
* A branch Ola granted (`away.py --grant <branch>`) may, once per head, be
  pushed (`git push origin <sha>:refs/heads/<branch>`, or `git push private
  ...` for commits of `private/master`), opened as a PR and enqueued, or
  merged on the private repository once green, if `@reviewer`'s last
  verdict approved that head. Grants end with the window.
* Anything else for Ola: `questions.py add --act`, then other work. Record
  each default taken, a first-question answer included, with
  `questions.py taken` as it is taken; the first report after the window is
  `questions.py list`. Whether to build at all waits for Ola.
