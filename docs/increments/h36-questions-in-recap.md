# Harness h36: the question tool in the recap, in place of ASK OLA lines

Status: **designed**, 2026-10-27, on `worktree-question-tool` off master
`7c2aa7c0`; built on its own branch after `h35-question-tool.md` merges. Short
path (`docs/increments/README.md`): about 40 counted lines added and 45
removed. Every code file it changes is governed, so the green step runs with
Ola at the guard's prompt, by day.

Source: as `h35-question-tool.md`; this is the half that stops the old routes.

## Prior art: legacy and literature

As `h35-question-tool.md`: no literature or legacy applies, no novelty is
claimed, no ideas-file entry is taken up.

## What changes

* `tools/session_state.py`: the recap's *Waiting on Ola* block prints
  `questions.list_lines()`, which is read-only (it stores no number mapping
  and marks no `taken` line shown), instead of reading `ASK OLA:` lines.
  `session.md` holds `NOW:` and `QUEUE:` only; an `ASK OLA:` line is warned as
  an unknown line. `tools/questions.py` joins the governed list
  (`.claude/hooks/guard_governance.py`), since a hook then imports it.
* `tools/work_pipes.py`: `ask` and `answer` are cut, with no alias; `show`
  marks an item waiting when the store holds an open question naming it, and
  prints no question block of its own.
* `tools/harness_mode.py` and `.claude/hooks/guard_unattended.py`: the refusal
  text tells the main session to run `questions add --act` (a subagent: put
  the question in its handback) instead of writing an `ASK OLA:` line.
* `tools/away.py --back` prints `questions.list_lines()`.
* Rule text, net words cut: `CLAUDE.md`'s "Two or more questions or acts
  awaiting Ola are numbered in one sequence, each explained, with its default"
  becomes "Questions to Ola are the output of `python3 tools/questions.py
  list`"; `.claude/REQUIRED-READING.md`'s `ASK OLA:` rules shrink to that
  pointer. Counted with `tools/rule_sizes.py`.

## Tests

1. The recap prints the open questions as `list_lines` does, and the store is
   byte-identical after it.
2. An `ASK OLA:` line in `session.md` is warned as not a `NOW` or `QUEUE` line.
3. `work_pipes ask` and `answer` are gone (the parser refuses them); `show`
   marks an item with an open question as waiting.
4. The unattended refusal text names `questions add --act`.

## Left out

Moving answered questions into a tracked file: answers about private work
must not reach the public repository, and the increment files already carry
Ola's rulings.
