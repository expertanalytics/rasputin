# The unattended night of 2026-10-02, the restart of 2026-10-03, and h7

Checked by `@orchestrator` on 2026-10-03 at master 390b516, its first run in
the role h7 gave it. Times are UTC. Proposals that need Ola's ruling are in
`next.md`; this file holds the evidence.

## 1. The night: Ola away 12:02 to 21:39

The window, from `<git-common-dir>/harness/windows.jsonl`: entered
2026-10-02T12:02:25, ended by `away.py --back` at 21:38:58 (9 h 36 min).

Work done: 15c-2 (geographic DEMs) to approval, `rasputin fetch` without
install metadata, the ANADEM accuracy note, the per-phase memory run on
sub-basin 761, the basin as nine level-3 pieces, the basin memory note, and
the h6 design. All reviewed, none pushed (Ola pushed them the next morning as
#139 to #145). No guard refused anything: there is no queue file for the
night (only `queue-2026-10-01.jsonl` exists).

### 1a. Seven hours idle

The last commit on any ref is 0b6305b, at 14:31:06. Nothing was committed on
any ref after it until Ola came back:

```
$ git log --all --format='%h %cI %s' --since='2026-10-02T14:32:00Z' --until='2026-10-02T21:39:00Z'
(no output)
```

So 7 h 08 min of a 9 h 36 min window was idle. This is the second night
running: on 2026-10-01 the session sat idle 3 h 25 min for the same reason
("the rest needs Ola's rulings"). Ola then gave the rule "fill the
unattended window" (stage a deep queue, fallbacks that need no rulings, a
blocked item only blocks itself), and it was kept as a memory note, not as a
rule file. One day later the session record said "idle after that BY
DESIGN: every remaining item needs Ola's ruling", and the night's queue in
`session.md` had no fallback items. A memory note did not hold through the
next night.

There was also a conflict the session settled on its own: Ola's plan for the
night was "no new increment implementation tonight", and filling the window
asks for work. Which one wins was not lifted to Ola.

**Correction to the brief.** The brief named h5's green step (red tests and
rulings on branch `worktree-h5-state-check`) as a decision-free fallback that
was missed. It was not decision-free at night. h5's design puts its main file
at `.claude/hooks/state_check.py`, moves the governed set out of
`.claude/hooks/guard_governance.py`, and changes `tools/away.py`. All three
are governed paths (`GOVERNED` and `GOVERNED_PREFIXES` in
`guard_governance.py`), so in unattended mode each write is refused and
queued, and a refused write to a rule file is never redone by any route
(`.claude/REQUIRED-READING.md`, *The harness*). A green step there would have
produced queued refusals, not code. The finding stands without it: the queue
named no fallback at all. Whether a decision-free item existed that night is
not settled here.

### 1b. `away.py --back` run from a worktree shows no `ASK OLA` lines

`tools/away.py` takes its root as `Path(__file__).resolve().parents[1]`, the
checkout the script was run from, and reads `root/.claude/current-task/`.
Probe from this worktree:

```
$ python3 -c "from session_state import pending_decisions; ..."   # in tools/
pending_decisions(<this worktree>/.claude/current-task)  ->  []
```

The main checkout's `session.md` held five open decisions at the time. This
was already recorded in `next.md` (findings of 2026-10-01/02) and not fixed;
it then bit in the trial it was found in.

### 1c. Decisions written as bullets under an `ASK OLA:` heading are missed

`pending_decisions` in `tools/session_state.py` lists any line containing
"ask ola", in any case. Run on a copy of `session.md` as it stands:

```
'session.md: ASK OLA:'
'session.md: 21:39 UTC Ola back (no decisions tonight). Trial findings ...
 "ASK OLA" are listed ...'
```

The five decisions under the heading (basin route, h6 questions 1 and 2, h5
by day, the pushes, the side notes) are not listed. What is listed is the
empty heading and a line that only talks about the marker. Two faults, one
on each side:

- the writer broke the format `.claude/REQUIRED-READING.md` sets (a decision
  is "a line containing `ASK OLA:`");
- the tool takes any mention for a decision and says nothing about an
  `ASK OLA:` line with no text after it.

### 1d. `session.md` has become a log

`.claude/REQUIRED-READING.md`: "Three lines or fewer per file", and
`session.md` "is overwritten in place and deleted when the round lands". It is
29 lines and 5,007 characters (`wc -lc`), appended to since the evening of
2026-10-02 across several rounds. The `SessionStart` recap is capped at
10,000 characters, so this one file takes half of it. Finding 1c follows from
the same habit: a log grows headings, and headings hide decisions.

## 2. The restart, 2026-10-03

Ola asked whether restarting the CLI for an update was safe. The main session
said yes, and that the background merge scripts "die on restart". They do
not: they run under the Claude Code daemon and kept going. After the restart
the main session started a second merge chain without checking `ps`, and the
two raced on #142 and #143. GitHub merged #142 once and, as the main session
reports, refused one early merge of #143; a refused merge leaves no event in
the PR timeline, so that part is not checked here. The merges landed in order
(`gh pr list --state merged`): #142 08:58:05, #143 09:10:01, #144 09:20:46.
No harm.

Two rules were broken, both in `.claude/REQUIRED-READING.md`:

- *Claims*: a claim about behaviour (what a restart does to running jobs)
  was given to Ola without running the check that settles it;
- *Before you act in a resumed session*: the resumed session acted (a new
  chain) before establishing what was still in flight.

Fixed so far: the `session.md` note now says background scripts survive a
restart and to check `ps` before resuming, and a memory note says the same.
Both are notes, not rules; section 1a shows a note can fail within a day.

## 3. h7, as merged (#146)

| Step | Evidence |
|---|---|
| design, `@architect` | 0d15560 |
| review round 1, CHANGES REQUESTED (five points) | `## Review` in `docs/increments/h7-orchestrator-role.md` |
| fixes, `@architect` | b6af28e |
| review round 2, APPROVED | same section; recorded in 43b6fe8 |
| push and PR, Ola ("merge once green") | #146 |
| merge commit, not squash | 390b516, 09:43:15 |

Docs only, 0 production lines: five files, all prose (`git diff --stat
7810cf8 390b516`). The review ran before the push, stayed within two rounds,
and left its trace. The route was right for a rules-only change: no red or
green step to perform.

Lessons:

- **The rewrite did not reach the run it was written for.** This run was
  spawned after the merge, yet the persona prompt it was given is the old
  one ("Master Agent & Project Orchestrator", the TDD loop), and the
  `CLAUDE.md` in its context has no "The main session dispatches" section.
  Both files on disk at 390b516 have the new text. Observed, not checked
  against Claude Code's documentation: persona files and `CLAUDE.md` are
  read when the session starts, so a session that merges a change to them
  keeps the old text, and so do the agents it spawns. This run followed the
  new role only because the brief said to read the persona file from disk.
  The main session is in the same position: the dispatcher rules h7 moved
  into `CLAUDE.md` are not in its context until it restarts.
- **`ROADMAP.md` has no row for h7**, nor for h4 or h6; h3 has one.
  `docs/increments/README.md` says the merge updates the increment's row. It
  is not settled whether harness increments belong in the roadmap.

## Review

**Round 1**, `@reviewer`, on 15ff2b1 (base master 390b516):
**CHANGES REQUESTED**. 0 production lines (two prose files). Not pushed, so
no CI. `python3 tools/check_citations.py`: all resolve, none at risk.
Checked and holding: the window (12:02:25 to 21:38:58), the last commit
0b6305b at 14:31:06 and the 7 h 08 min; no refusals (`--back` printed
"Queued while away (0)"); `away.py`'s root and `pending_decisions`' match;
h5's three files are governed; the h7 table (0d15560, b6af28e, 43b6fe8,
390b516 with two parents); the merge times; the old persona prompt and the
29 lines and 5,007 characters (both in this run's transcript); the quoted
session lines. Proposals are framed for Ola's ruling. Blocking:

1. "Second night running" is false by the file's own method. The window
   2026-10-01T22:43:29 to `--back` at 2026-10-02T07:06:21 had no commit on
   any ref from 2331110 (01:06:58) until after the window: 5 h 59 min idle.
   That makes three windows running. Fix §1a and `next.md` item 1.
2. `session.md` is described in the present tense ("It is 29 lines", "as it
   stands", `next.md` item 4 "is a 29-line log"). It has since been
   rewritten to 4 lines. State it as measured on 2026-10-03, in the past
   tense.
3. "h3 has one" (§3 and `next.md` item 7): h3 has no row in `ROADMAP.md`'s
   table. It is item 1 of the "Order of work from 2026-09-30" list. No
   harness increment (h2 to h7) has a table row.
4. §1b says the fault "bit in the trial" but names no evidence (persona §1:
   each finding names its file, commit or transcript line). The evidence
   exists: main-session transcript `85c14e7c`, 2026-10-02T21:38:58Z,
   `--back` run with cwd `.claude/worktrees/basin-memory`, printing "ASK OLA
   lines in .claude/current-task/: (none)". Cite it.

Non-blocking: the "night" ran 14:02 to 23:39 local time. None of the seven
proposals states what it would cost to adopt. h5's new `tools/governed.py`
and `tools/git_hooks.py` are not yet governed paths, so "queued refusals,
not code" is a little too strong (the main file is governed, so the finding
holds).
