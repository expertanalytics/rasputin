# 2026-10-26 team review: context load, flow, and the watcher's blind spot

Times are UTC. Sources: today's 110 subagent transcripts under
`a8ebd8a5-7c29-4760-8e5d-791db2aa83cc/subagents/` (05:03 to 16:27), the main
transcript (line numbers below), `git log`, `gh run list`, and the two
unpushed checks on `worktree-retro-merge-365` (`b76f7285`, `e000d3d1`), whose
flow figures this file uses rather than repeats. Token counts are characters
divided by four. **This file is public**; it names flaws, not culprits.

`tools/duties.py` at 16:26: post-merge check owed for #367 (merged 15:36,
50 min) and #364 (16:01, 25 min); both are covered by `e000d3d1`, not yet on
master. **The brief's own example does not hold**: no agent read the
Wednesday increment file whole except three designers; see §1.

## 1. Context load: what each agent was given against what it used

**Fixed load.** Every subagent's first call carries 24 to 26 k tokens
(median 25.5 k, all personas alike) before it reads anything: the harness
prompt and tool schemas (about 14 k, not ours), the skill list (about 4.6 k),
`CLAUDE.md` (2 k), the memory index (3.5 k; 27 of its 87 lines are
land-surface, inflow or weather) and the persona file.

**Reads, all 110 agents (tokens of tool results).**

| what | agents | tokens | used by a subagent |
|---|---|---|---|
| reading code, `grep`/`sed` over the tree | most | 954 k | the work |
| `git` / `gh` output | most | 227 k | the work |
| `.claude/REQUIRED-READING.md` | 90 | 232 k | about a third: *Before you write*, *Claims*, *Stale artifacts*, *Data*, the note-file lines (about 570 of 1,759 words); the rest is the main session's (recovery, `session.md`, publishing, the harness) |
| increment files | 55 | 213 k | the agent's own fix: for a one-fix agent about 130 of `106-wednesday-cells.md`'s lines |
| persona file read again from disk | 25 | 51 k | none: it is already the system prompt |
| `docs/increments/README.md` | 20 | 51 k | the step it is on |

**Sample of 13** (peak context; rule and increment reads in k tokens):

| persona, task | peak | REQUIRED-READING | increment | work reads |
|---|---|---|---|---|
| developer, 106d reservoirs | 85 | 4.2 | 4.9 | 14.4 |
| tester, 106d red | 94 | 3.0 | 3.4 | 20.5 |
| reviewer, 106d code | 41 | 1.1 | 5.3 | 0 |
| reviewer, merge check fix 4 | 33 | 0 | 0 | 2.9 |
| architect, merge master into fix 4 | 38 | 1.0 | 0 | 2.7 |
| architect, design 106e | 139 | 3.8 | 10.8 | 27.9 |
| tester, 106e red | 154 | 3.2 | 5.2 | 35.0 |
| developer, 106e green | 117 | 3.3 | 5.4 | 33.0 |
| perf, 106e acceptance | 92 | 4.2 | 4.1 | 14.4 |
| developer, one serial CI run | 39 | 3.4 | 0 | 3.6 |
| reviewer, parallel runner | 33 | 0 | 0 | 2.2 |
| developer, lsm slice 1 | 80 | 3.6 | 0 | 15.7 |
| orchestrator, check of #362 | 90 | 2.3 | 2.7 | 16.4 |

**Unneeded, per agent:** about 6.5 k tokens of rules and memory (the
main-session parts of `REQUIRED-READING.md`, `CLAUDE.md` §3 and the memory
index, the persona re-read), plus about 3 to 7 k of other fixes' sections
for an agent on one of the five fixes (`adc140a73`, 106d green, read 9.1 k
of the file; its own sections are about 2 k). That is 20 to 25 % of a small
agent's context, under 10 % of a large one's. Rule files carry no
land-surface text; the memory index does.

**The bigger load is the dispatcher.** The main session went from 107 k
tokens at 05:02 to 912 k at 16:02 in 947 turns, no compaction, while
dispatching four pipes: of the 110 spawns, 56 Wednesday fixes, 22 harness
(pipes, heartbeat, parallel tests, CI), 14 land-surface, 9 retrospectives.
Whether turn latency grew with it is not measured here.

## 2. Flow: where the five fixes' time went

From `b76f7285` §1 and `e000d3d1` §8. Red step to last approval:
40 min (106a), 1 h 23 (106b+c), 55 min (106d), 1 h 11 (106e): **4 h 09 in
all**. Approval to merge: 2 h 25, 3 h 37, 3 h 36, 3 h 30: **13 h 08**.

| cause | measured | evidence |
|---|---|---|
| waiting on the branch below (stack) | 61 + about 130 + 69 min | `b76f7285` §1, `e000d3d1` §8 |
| base merges, each a serial full suite | 6 agent runs of 22 to 30 min | agents `a743c457`, `ab967ed2`, `a60f75c5`, `a8af331b`, `abdd313f`, `a23bf60f` |
| Ola's push yes | 48 min (106a, 09:02 to 09:50), 56 min (106b+c, 12:01 to 12:57) | transcript 36137, 38749 to 39027 |
| macOS-only digest tests | 2 red CI runs, 1 review round, 1 PR out of the queue | `86fd5208`, `31101d5b`, `54b228a4` |
| CI and queue per pass | 30 to 37 min before #367, about 17 after | `gh run list` |
| shared design before any red step | 06:27 to 07:40, fixes 106a to 106d | agents `a05612fb`, `a36f5213` |

106a's 52-minute rework (09:57 to 10:49) was a real slowdown found by
`@perf`; that time bought something. The short path was used (106a, b, d)
and worked: the build was not the slow part.

## 3. Why my checks named this late

- **My watch list is about rules, not time.** `orchestrator.md` §1 lists
  skipped steps, role bleed, guards worked around, idle away windows and
  false "done". Nothing asks where wall time went between approval and merge,
  or what agents carry. `lead_time.py` splits lead time only into build and
  wait.
- **I run after a merge, one PR at a time.** The first check (12:00) saw only
  #362, whose wait was its own rework; the stack showed only across PRs, and
  those merged 14:03 to 16:01. No duty runs while work is in flight.
- **Briefs hand me the findings.** The 12:00 brief listed what to cover
  ("the away window ... no idle stretch over 13 min") and "Keep it small";
  the check followed the list.

So the stack, the serial suites and the waits went unnamed until Ola asked at
13:24, and were measured at 15:24 (`b76f7285` §4).

## 4. Proposals, ranked by measured gain (Ola decides)

1. **No stacked branches.** Each fix branches from master; a fix that needs
   another waits for its merge or shares its PR. Gain: up to 4 h 20 min of
   waiting today, and the merge of a base into its stack (`3204b29d`, which
   also carried the macOS digest). Which fixes truly needed the one below is
   unchecked. Cost: about +25 words in
   `docs/increments/README.md`.
2. **Ola's push yes for named branches while he is present**, as `away.py
   --grant` does while he is away. Gain: 1 h 44 min today. Cost: about +20
   words in `REQUIRED-READING.md`; it moves the approval boundary.
3. **One increment file per PR**, with a shared design kept short in its
   own file. Gain: 5 of today's 6 base merges on the fixes conflicted, each
   only in that file (`b1eab4dd`, `d945ecc6`, `3204b29d`, `7f106eb1`,
   `a991d3d7`), and 3 to 7 k tokens per later-fix agent. Cost:
   about +8 words.
4. **A draft PR at the red step**, so CI runs from the first commit (CI
   runs only on pull requests, `.github/workflows/main.yaml`). Gain: 106c's
   macOS digest would have shown at 09:03, not 12:57. Cost: a push yes per
   branch, about +15 words.
5. **A subagent reading list.** Split `REQUIRED-READING.md` into the
   subagent third and the main-session rest; `brief.py` asks for the
   persona file only after a rule merge; land-surface memory lines move to a
   file the land-surface pipe points at (Ola's file). Gain: about 6.5 k
   tokens per agent, 0.7 M a day; minutes unmeasured. Cost: net 0 words
   (moved), one pointer line.
6. **One dispatcher per pipe** (mesh, land-surface, harness), each its own
   session and worktree; the three-item cap read across them from
   `tools/work_pipes.py`. Gain: the 107 k to 912 k growth splits; minutes
   unmeasured. Cost: about +40 words, and Ola talks to several sessions.
7. **Flow in my watch list, and a mid-day trigger**: each check splits the
   wait (base, suites, Ola, CI, queue), and the recap warns when an approved
   branch is unmerged after an hour. Cost: about +30 words in
   `orchestrator.md`, paid by proposal 8.
8. **Cut**: `CLAUDE.md` §3 *When to spawn `@orchestrator`* restates what
   `tools/duties.py` computes; becomes "when `tools/duties.py` shows a duty
   owed". About -35 words. Plus the open `reviewer.md` cut (-17,
   `2026-10-26-merge-357.md`).

Already taken: the parallel suite (#367): 22 min to 3 min locally, CI 16 to
18 min down to 9.

## 5. Rule text

`python3 tools/rule_sizes.py`: 9,494 words, unchanged since
`2026-10-26-merge-362.md`. Proposals 1 to 7 add about 140 and move about
1,200; proposal 8 and the reviewer cut take out 52. Anecdote still in a rule
file: the `cbrt` example, `tester.md` line 26 (`b76f7285` proposal 1).
