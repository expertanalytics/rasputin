# Harness h23: duties owed, computed from git and the away log

Status: **design**, 2026-10-17, on `worktree-h23-duties` off master
`81ea2f73`. One PR. Every production file is governed (`tools/`,
`.claude/settings.json`), so the green step and the settings change wait for
Ola, by day.

Ola, 2026-10-17: "Go, is there a way to get you to work more regularly on
"duties"?"

Implements proposal P1 of the 2026-10-17 overdue-check retrospective
(`docs/retrospectives/2026-10-17-overdue-check-night-merges-audit.md`,
section 0 and section 5 P1, on `worktree-retro-1017` at `0e1616a2`, PR #282).
Its finding: every trigger for an `@orchestrator` run is the main session's
choice, so a lapse leaves no trace; 10 merges and three away windows passed
unchecked.

Acceptance: the harness suite (`pytest -m harness`) and the three live
checks of section 7, recorded in this file. No `@perf` run (no refine or mesh
code).

Words used here. **Duty**: an `@orchestrator` run the main session owes
under `CLAUDE.md` §3, *When to spawn `@orchestrator`*. **Away window**: a
stretch of unattended mode, from `tools/away.py <duration>` to `--back` or
its end time. **Check commit**: a commit whose subject ends
` (@orchestrator)`, the suffix every persona commit carries (the
`tools/brief.py` block asks for it). **Retrospective time**: the first time
a check commit added a given `.md` file under `docs/retrospectives/`, on any
local branch or on `origin/master`. **Owed**: the event that calls for a
duty is later than the last retrospective time.

## 1. Prior art: legacy and literature

Tooling; no novelty claimed.

*Literature* (each read 2026-10-17):

- Healthchecks.io, *Docs*: "Healthchecks.io works as a dead man's switch for
  processes that need to run continuously or on a regular, known schedule.
  … It raises an alert as soon as a ping does not arrive on time." This
  design is that switch with git as the ping log: a retrospective commit is
  the ping, and a merge or a finished night opens a deadline. Departure: no
  server and no schedule of its own; the deadline is "before the next
  event", not a clock, because the rule says "after each".
- Prometheus operator runbook, *Watchdog*: "This is an alert meant to ensure
  that the entire alerting pipeline is functional. … This alert is always
  firing." Taken up in the status line only: it always prints something
  ("Duties: none owed"), so a blank status line means the tool broke, not
  that all is well. The per-prompt hook stays silent when nothing is owed,
  because it costs context on every turn.
- Claude Code, *Hooks reference*, `UserPromptSubmit`
  (code.claude.com/docs/en/hooks): it fires on typed prompts and also on "A
  background subagent reporting back to the session that started it" and
  scheduled tasks; on exit 0 "Claude Code adds stdout it treats as plain
  text to Claude's context"; "Neither channel produces a visible transcript
  entry" (a system reminder Claude reads, which Ola does not see); exit 2
  "Blocks the prompt, so it never reaches Claude"; default timeout 30 s, and
  a hook that reaches it "is canceled and its output … is discarded"; output
  is capped at 10,000 characters. The page does not say whether the hook
  fires inside a subagent (section 8, pin 2).
- Claude Code, *Status line* (code.claude.com/docs/en/statusline): set in
  user or project settings as `statusLine` with `type: "command"`; runs at
  session start and "when a new assistant message arrives", debounced at
  300 ms; `refreshInterval` re-runs it every N seconds, for when "the main
  session is idle"; "Slow scripts block the status line from updating";
  stdin carries JSON with `workspace.project_dir`. Ola has no status line
  today (`grep -n statusLine ~/.claude/settings.json` finds none).

*Legacy.* Nothing: the legacy tree is the old meshing product, with no
harness.

```
$ git grep -n -i -l 'duties\|overdue\|retrospective\|statusline' legacy-archive -- legacy
$ (no files)
```

*Ideas file.* Overlaps two entries and takes up neither: *Recap: items in
flight and open ideas* and *Lead time from design to merge, per increment*
(both add lines to the recap; this adds one more, so the recap's 10,000
character cap is shared).

## 2. The duties, and how each is computed

No state of its own: every input is already on disk, and the tool writes
nothing. Four git calls and one file read, in all.

- **R, the last retrospective time.** One call:
  `git log --branches origin/master --diff-filter=A --format='C %cI %s' --name-only -- 'docs/retrospectives/*.md'`.
  Only check commits count. Not every new file there is a check: the
  retrospective agenda `next.md` (`1b60375f`, by the main session for Ola)
  and the rule-file sweep (`403e2c9e`, which added
  `2026-09-30-rule-file-history.md` beside rule edits) would each have
  cleared every duty. Since 2026-09-30 every commit that adds a file there
  carries the suffix (`git log --branches origin/master --diff-filter=A
  --since=2026-09-30 --format=%s -- 'docs/retrospectives/*.md' | grep -v
  '(@orchestrator)$'` prints nothing); the older ones without it are all
  earlier than any R that matters. `origin/master` beside the local
  branches, because a check merged and its local branch deleted is still a
  check. For each path, keep its earliest add time; R is the latest of
  those. Earliest, because a retrospective carried from an
  unpushed branch onto a newer one (`a62b3edb`, `f35f8969` carrying
  `b3e05434`) is re-added with a new date and would otherwise look like a
  fresh check. Any branch, not only master: a retrospective written and not
  yet pushed is a duty done; the unpushed branch is its own line (D4). On
  this tree the rule counts 24 of the 32 files there as checks,
  and R is the same as before (`0e1616a2`).
- **The away log**: `windows.jsonl` in `harness_mode.state_dir(repo)`, read
  only. The tool never names that path in a shell line, so the guards do
  not see it.

The four duties, in output order:

| # | Duty | Owed when | Rule it checks |
|---|---|---|---|
| D1 | post-merge check | a merge on `origin/master` (first parent) is later than R and is not a check-recording merge (below) | "After each increment merges" |
| D2 | morning check | an away window of at least `NIGHT` (4 h) ended after R and before now | "the morning after each unattended night" |
| D3 | research round | no file under `docs/retrospectives/` whose name contains `-research` was first added in the last `WEEK` (7 days) | "for a research round, weekly" |
| D4 | push retrospective branches | a local branch adds a retrospective file that `origin/master` does not hold | the five branches the retrospective found unmerged; recap and status line only |

D1's query: `git log --first-parent --merges --diff-merges=first-parent
--name-only --format='M %cI %s' --since=<R> origin/master` (`master` if
`origin/master` is missing). A **check-recording merge** is left out, or
merging a check would make the next one owed: its diff against the first
parent touches at least one path under `docs/retrospectives/` and nothing
outside `docs/retrospectives/`, `ROADMAP.md` and `docs/increments/`. That
covers #228 (ROADMAP and a retrospective) and #268 (ROADMAP, the increment
file of 46 and a retrospective). The retrospective path is required: #208
and #237, audits that touch only ROADMAP and increment files, are work and
owe a check. The PR number comes from
the subject `Merge pull request #N`. It reads `origin/master` as last
fetched and never fetches: no network on a prompt. The main session fetches
before every push question (`CLAUDE.md` §3), so the lag is short.

D2's window ends, from the events in order: a `back` event ends the window
it names (its `until`); an `enter` with no `back` before the next event ends
at its own `until`, or at the next `enter`'s `since` if that is earlier
(log lines 39-40 and 46-47 are such cases). A window whose end is after now
is still running and owes nothing yet. Any retrospective time after the end
clears it, a merge check included: simple, and in practice the morning
retrospective is the next one written.

D3's marker is new: `.claude/agents/orchestrator.md` §2 gains "A research
round's retrospective file name contains `-research`." Until the first one
lands, D3 reads "none on record", which is true (no research round since
about 2026-10-12, per the retrospective).

D4's query, two calls:
`git log --branches --not origin/master --no-merges --diff-filter=A --source --name-only --format='S %S' -- 'docs/retrospectives/*.md'`,
minus the paths in `git ls-tree --name-only origin/master docs/retrospectives/`.
A branch Ola drops is deleted, and the line goes. Only Ola can push, so D4
stays out of the per-prompt hook, where the main session could do nothing
with it; it shows in the recap and the status line, and the hook skips its
two calls.

**Constants.** `NIGHT = 4 h`: the away log holds 43 ended windows,
2026-10-01 to 2026-10-17 (83 lines, 41 `back` events and two windows that
ran out), lasting 15 min to 9.6 h. The 15 nights run 6.6 to 8.9 h; 23 day
windows are under 4 h (the longest 3.9 h), and five day windows are longer
(4.2, 4.4, 6.3, 7.4 and 9.6 h). Those five owe a check too: idle time and
guard refusals pile up the same way by day, and the retrospective counted
the 7.4 h one of 2026-10-16 as unchecked. `WEEK = 7 days`, from the rule's
"weekly". Both checked against that log only (durations read from it by
hand).

**Not a duty here: the recap.** "Recap each round" leaves no trace in git or
the away log, so it cannot be computed without state. The `SessionStart`
hook already prints it on start, resume, `/clear` and compaction.

## 3. The output: one line

Nothing owed: no output (hook and recap); `Duties: none owed` (status line).
Owed: one line, the owed items joined by `; `, in D1 to D4 order; the hook
leaves D4 out, so D4 alone prints nothing there:

```
Duties owed: post-merge check (1 merge since the last retrospective, 2026-10-17 09:40 UTC: #280); morning check (1 away window ended since, latest 2026-10-17 06:45 UTC); research round (none on record); push 6 retrospective branches (worktree-retro-1004c, worktree-retro-1005, +4)
```

Item forms, exactly:

- D1: `post-merge check (<k> merge[s] since the last retrospective, <R>: <#a, #b, #c>[, +<n>])`
  (the newest three PR numbers; a merge with no number shows its short sha);
  with no retrospective at all, `post-merge check (no retrospective on record)`.
- D2: `morning check (<k> away window[s] ended since, latest <end>)`.
- D3: `research round (last <date>)` or `research round (none on record)`.
- D4: `push <k> retrospective branch[es] (<b1>, <b2>[, +<n>])`, sorted.
- Times: `YYYY-MM-DD HH:MM UTC`; dates `YYYY-MM-DD`.
- On any fault: `Duties: unknown (<ExceptionName>: <message, 120 chars>)`,
  in all three places, because silence would read as "none owed".

Run against this tree on 2026-10-17 (a throwaway prototype in the session
scratchpad, reading git only): 32 retrospective files, 24 of them checks
(rerun with the check-commit rule), R =
2026-10-17 09:40 UTC (`0e1616a2`), one merge since (#280), no research
file, six unmerged branches (`worktree-retro-1004c`, `-1005`, `-1005b`,
`-1008`, `-1015`, `-1017`). The away window that ended 06:45 UTC is before R,
so D2 is not owed; the line above shows its form only.

## 4. Where it appears

1. **The recap**, `tools/session_state.py`: one line after `Last landed:`,
   through `import duties` inside a `try` that costs one line on a fault,
   as `harness_mode` and `rule_sizes` already do there.
2. **A `UserPromptSubmit` hook**, the main session's context on every turn
   and on every background handback. A new entry in `.claude/settings.json`:
   `{"type": "command", "command": "python3 \"$CLAUDE_PROJECT_DIR/tools/duties.py\"", "timeout": 5}`.
   The tool always exits 0: exit 2 would block Ola's prompt. The 5 s timeout
   bounds a hung git; a timeout drops the line and the prompt goes through.
   With no arguments the tool runs as the hook: D1 to D3, two git calls.
3. **Ola's status line**: `"statusLine": {"type": "command", "command":
   "python3 \"$CLAUDE_PROJECT_DIR/tools/duties.py\" --status",
   "refreshInterval": 600}` in `.claude/settings.json`. `--status` reads the
   stdin JSON and runs against `workspace.project_dir`. The 600 s refresh
   keeps it current overnight, when no assistant message triggers it.
   Cost to Ola, from the status-line page (read 2026-10-17): "With a custom
   status line configured, Claude Code stops showing most of the footer's
   keyboard hints, including `esc to interrupt`, the `? for shortcuts`
   fallback, and the `hold space to speak` voice dictation hint." The
   built-in footer badges stay. A question for Ola; the default is to add
   it. If he says no, this item and pin 1 go, and the recap and hook
   remain.

Both settings entries need Ola's yes (`.claude/REQUIRED-READING.md`, the
approval boundary), and the hook list there gains one sentence for each.

**Time cost per prompt**: the prototype took 64-65 ms inside Python over
three runs and 90 ms wall clock with interpreter start (`/usr/bin/time -p`;
bare `python3 -c pass` is 20 ms), at 211 local branches, 32 retrospective
files, 83 away-log lines, running all four git calls; the hook runs two of
them. Context cost while something is owed: one line of
about 300 characters (about 75 tokens) per turn.

## 5. If the main session ignores the line

The hook only informs; nothing can make the main session spawn. Three things
keep the lapse visible anyway:

1. **The status line is Ola's, not the model's.** It sits in his terminal
   whatever the main session does, and the counts grow with each merge and
   night ("3 merges since"), so a lapse is seen within one turn of Ola
   looking, not three days later.
2. **It is computed, not remembered.** A new session, a compaction or a
   resume gets the same line from the recap; no note can go stale.
3. **The watcher measures the lapse.** `.claude/agents/orchestrator.md` §1
   gains: "Each run starts with `python3 tools/duties.py` and records, for
   each duty listed, how long it was owed." Ignoring the line then becomes a
   finding with a duration, in the next retrospective.

## 6. Components

```
git refs ──┐                    ┌─> session_state.py (recap line)
           ├─> gather(repo) ──> owed(inputs, now) ──> render() ─┼─> stdout (UserPromptSubmit hook)
windows.jsonl ┘   (I/O shell)      (pure)            (pure)     └─> --status (status line)
```

`tools/duties.py`, mypy-strict (added to `pyproject.toml`'s `files`):

- `Window(start, end)`, `Merge(when, label, paths)`, `Inputs(retro: datetime | None,
  research: datetime | None, merges, windows, unmerged)`: frozen dataclasses.
- `window_ends(events: list[dict[str, Any]], now) -> list[Window]`: pure;
  a malformed line is skipped.
- `owed(inputs: Inputs, now) -> list[str]`: pure; the item strings of
  section 3.
- `render(items, place: Literal["hook", "recap", "status"]) -> str | None`:
  the hook drops D4; only the status line prints `Duties: none owed`.
- `gather(repo: Path, unmerged: bool) -> Inputs`: the four git calls (two
  when `unmerged` is false: R and D1) and the log read.
- `main(argv, stdin) -> int`: always 0; any exception becomes the
  `Duties: unknown` line.

Estimate: about 85 net production lines (`duties.py` about 80,
`session_state.py` 5). Over the 50-line short path, so the full loop. The
settings JSON is not counted.

Prose in the same PR: `CLAUDE.md` §3 *When to spawn `@orchestrator`* gains
"`tools/duties.py` shows which are owed (recap, every prompt, Ola's status
line)"; the orchestrator lines of section 2 and section 5;
`.claude/REQUIRED-READING.md` *The harness* names the hook and the status
line.

## 7. Tests

`tests/python/test_duties.py`, marked `harness`. Not invariant-critical: no
mutation round.

1. `window_ends`: a `back` pair gives its `until`; an `enter` with no `back`
   ends at the next `enter`'s `since` when earlier than its `until`, else at
   its `until`; a running window (end after now) is dropped; a malformed
   line is skipped. Fixture: six synthetic events copied in shape from log
   lines 39-40, 46-47, 81-83.
2. D2 threshold: a window of 3 h 59 min ending after R owes nothing; 4 h 0
   min owes; a 9 h window ending one minute before R owes nothing.
3. D1: a merge one second after R owes; one before R does not; check-recording
   merges do not, with paths `docs/retrospectives/x.md` alone and with
   `ROADMAP.md`, `docs/increments/46-x.md` and `docs/retrospectives/x.md`
   (the shape of #268); a merge with only `ROADMAP.md` and
   `docs/increments/46-x.md` owes (the shape of #237); five merges list the
   newest three numbers and `+2`.
4. D3: no research file gives `none on record`; first added 6 days 23 h ago
   owes nothing; 7 days 1 h ago gives `last <date>`.
5. `render`: the golden line of section 3 for a fixed `Inputs` and now, in
   `recap` and `status`; in `hook` the same line without its D4 item, and
   `None` when D4 is the only item; nothing owed gives `None`, or
   `Duties: none owed` in `status`.
6. `gather` on a git fixture in `tmp_path` (fixed `GIT_COMMITTER_DATE`s,
   `refs/remotes/origin/master` set with `update-ref`): a retrospective
   added on branch A at t1 and the same path re-added on branch B at t2 > t1
   gives R = t1; a code merge on master after t1 is counted, a
   retrospective-only merge is not; A and B are listed in D4 until the path
   is on `origin/master`. Non-check files: a later commit at t3 > t1 on
   branch C, subject `Retrospective agenda (Ola)`, adding
   `docs/retrospectives/next.md`, leaves R = t1 (the shape of `1b60375f`);
   the same commit with subject ending ` (@orchestrator)` gives R = t3. A
   check commit only on `origin/master` (no local branch holds it) at
   t4 > t3 gives R = t4. With `unmerged=False`, the D4 list is empty. About
   11 commits.
7. `main`: in a directory that is not a repository it prints the
   `Duties: unknown` line and returns 0; it never returns 2.
8. The recap prints the line after `Last landed:`, and a `duties` that
   raises costs one line, not the recap.

Live checks, by day, recorded here: (a) the hook's line in a fresh session's
debug log (`claude --debug`) after one prompt, and its absence in a
subagent's; (b) the status line in Ola's terminal, started from a worktree
directory; (c) `time python3 tools/duties.py` on the main checkout, under
the 5 s timeout with margin.

## 8. Pinned or assumed beyond the sources

1. `$CLAUDE_PROJECT_DIR` is set for the status line command. The status-line
   page documents it for hooks only; live check (b) settles it. Fallback:
   the command `python3 "$(git rev-parse --show-toplevel)/tools/duties.py" --status`.
2. `UserPromptSubmit` does not fire inside subagents (their turns are not
   prompts). Not stated in the docs; live check (a). If it does, the line
   reaches subagents too: noise, not harm.
3. Commit dates are trusted: a commit with a forged or skewed committer date
   moves R.
4. A check is a commit by `@orchestrator`, known by its subject suffix. A
   check written by another persona, or with the suffix left off, clears
   nothing; an `@orchestrator` commit that adds a non-check `.md` there
   clears duties. The first has not happened since 2026-09-30 (the command
   in section 2); the second is not checked, and is the price of keeping no
   state.

## Review

- **Design review, round 1, 2026-10-17.** Commit `d72cd744`, 0 counted lines, @reviewer: CHANGES REQUESTED. B1: any new file under `docs/retrospectives/` counts as a check, and non-check files exist there (`next.md`, the rule-file sweep); count only commits whose subject ends `(@orchestrator)`, with a test case. Held: an unmerged check clears D1/D2 because D4 lists the branch; hook claims match the docs; 0.10 s measured; "none owed" is not noise; the 3-day gap would have been flagged from #268 on. Suggestions: D4 out of the per-turn hook; check-recording PRs that touch only ROADMAP and an increment file should not owe a check; read `origin/master` too; four git calls, not three.
