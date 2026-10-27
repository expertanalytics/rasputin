---
name: away-mode
description: What to do while Ola is away (unattended mode). Load when the recap says unattended mode is on, when Ola says he is leaving, or when a guard refuses an act.
---

# Unattended mode

* Only Ola starts it (`python3 tools/away.py <duration>`, his own terminal)
  and ends it (`--back`). No agent runs `away.py` or writes
  `<git-common-dir>/harness/`.
* When he says he is leaving, ask how long and ask him to run `away.py`.
  Before he goes, `QUEUE:` names a fallback needing no ruling and no
  governed write, and one stays queued for the last hour.
* A guarded act is refused and already queued: never retry it or work
  around it, except:
  - a refusal naming only ordinary files (no rule file, harness state or
    push) is a false positive: redo it with Edit or Write and report a
    `GUARD FALSE POSITIVE:` line with the files;
  - a refused push, PR or merge whose sole reason has gone (the head now
    approved by `@reviewer`, the branch now granted, the form now right).
* A branch Ola granted (`away.py --grant <branch>`) may be pushed
  (`git push origin <sha>:refs/heads/<branch>`; `git push private ...` for
  commits of `private/master`), opened as a PR and enqueued, or on the
  private repository merged once green, once per head, if `@reviewer`'s
  last verdict approved its head. The grant ends with the window.
* Anything else for Ola is an `ASK OLA:` line (main session: `session.md`;
  subagent: handback), then other work. A default taken is an `ASK OLA:`
  line ending "default taken: <what>"; the first report after the window
  lists them all. Whether to build at all waits for Ola.
