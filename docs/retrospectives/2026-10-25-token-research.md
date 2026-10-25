# 2026-10-25 research round: what fills the main session's context, and a plan

Asked for by Ola after he rejected `/clear` by hand (proposal 1 of
`2026-10-25-merge-352.md` §4, now withdrawn). Times UTC. Source: the main
transcript `a8ebd8a5-7c29-4760-8e5d-791db2aa83cc.jsonl`, from the
auto-compaction at 09:56:55 to 20:46 (593 model calls). **Public file**:
private-package work is counted, never described.

## 0. Duties owed

`python3 tools/duties.py` at the start printed nothing: none owed.

## 1. Measured: what fills the main session's context

Method: each call's context (input + cache read + cache write, from its
`usage` block) minus the previous call's is the growth; it is shared among the
items that arrived in between by their size (characters / 4; an image counted
as 1,600 tokens, the documented ceiling, unchecked against today's page).
Growth 09:57-20:46: **783 k tokens** (65 k after the compaction to 851 k).

| What | Share of growth | How |
|---|---|---|
| **Brief blocks** | **34 %** (263 k) | `brief.py` output written to a temp file and read back (median block 2,894 characters), then retyped as the `Agent` prompt (median task wording only 996 characters). Each of 84 spawns carries the block twice, once as output tokens |
| Handbacks, launch notices | 20 % (157 k) | 84 handbacks, median 3,556 characters (about 550 of them the harness's frame); each background launch notice about 420 tokens |
| Reading files, diffs, logs, `gh`, conflict checks | 17 % (134 k) | 92 shell reads/greps, 32 `git show/diff/log`, 62 `gh`, 14 merge-tree/merge |
| Main session's own prose and thinking | 14 % (111 k) | answers to Ola, recaps |
| Harness reminders and hook output | 9 % (71 k) | per-turn reminders, duty line, guard texts |
| Plots shown to Ola | 2 % (13 k) | 10 images |
| Wording edits by hand | 1 % | 2 `Edit`, 16 `Write` calls (mostly `session.md`) |
| Not attributed | 3 % (about 26 k) | |

Cost follows size: every call re-reads the whole context from cache, so the
weighted cost of a call grows with everything that came before it (Claude Code
docs, *Why usage climbs in a long session*). Since the session runs a 1 M
window, auto-compaction fires only at about 967 k: three times this week
(22 Oct 17:45, 24 Oct 03:45, 25 Oct 09:56, `compact_boundary` entries), each
about 50 s.

**The state layer has overflowed too.** `session.md` is 84 lines, 21,969
characters: 46 `ASK OLA:`, 21 queue lines, 7 ruling lines, 24 lines over 300
characters (the rule: one `NOW:`, one `QUEUE:`, `ASK OLA:` lines, no rulings,
300 characters each). 51 subagent note files linger (the spawner deletes each
on the handback). The recap is **41,932 characters against its cap of
10,000** (`python3 tools/session_state.py | wc -m`), and it is reprinted at
every compaction.

## 2. Sources (each read for this round)

- Anthropic, *Effective context engineering for AI agents* (29 Sep 2025):
  compaction; "structured note-taking ... persisted to memory outside of the
  context window"; a subagent "returns only a condensed, distilled summary of
  its work (often 1,000-2,000 tokens)"; "tool result clearing" as the lightest
  compaction; agents keep "lightweight identifiers (file paths ...)".
- Anthropic, *Effective harnesses for long-running agents* (26 Nov 2025): each
  session starts fresh and catches up from `claude-progress.txt` and git; the
  model changes JSON status fields more safely than Markdown.
- Anthropic, *How we built our multi-agent research system*: the lead saves
  its plan to memory before truncation; "implement artifact systems" so
  subagents' output does not pass through the lead (the "game of telephone").
- Claude Code docs, *Model configuration*, *Set the auto-compact window*:
  `autoCompactWindow` in a settings file, 100 k to 1 M; default for 1 M models
  about 967 k. *Explore the context window*: after compaction, `CLAUDE.md`
  and memory are re-injected and `SessionStart` hooks matching `compact` run.
  *Prompt caching*: a warm-cache compaction "costs a fraction of what the
  context size suggests".
- Chroma, *Context rot* (14 Jul 2025, 18 models): "LLMs do not maintain
  consistent performance across input lengths". Long context costs quality,
  not only tokens.

## 3. The plan (one page)

Ola's frame: layers with separate concerns. **The dispatcher holds state and
decisions; tools hold mechanics; agents hold content; files are the
interfaces.** Today the main session crosses layers in three places: it
carries the brief's mechanics through its own context (34 %), it takes in
agents' content in full (20 %) and reads diffs and logs to check it (17 %),
and its state file has become a history (§1).

1. **Compact automatically at 300 k** (settings, no rule words).
   `"autoCompactWindow": 300000` in `.claude/settings.json`. This is the
   scalable `/clear`: no one types anything, the recap hook already reloads
   state from disk after each compaction, and in-flight work is on disk by
   rule. Replayed on today's calls: the main session's weighted cost from
   09:57 falls from 30.2 M to 13.5 M (**-55 %**), with 3 compactions instead
   of 0 (about 2.5 min). Cost: Ola's yes on a settings edit; one check that
   the project scope is honoured (`/autocompact` reports the window in force).
   Unchecked: a summary may drop a detail; the recap and `session.md` are the
   guard against that, so step 2 comes with it.
2. **State back to state** (a tool, no new rule). The pipes' first slice
   (`docs/increments/h31-pipes.md`) moves the queue out of `session.md`;
   rulings go to the increment file or here, never `session.md`. A
   `SubagentStop` hook deletes the finished agent's note file (mechanics, not
   the dispatcher's memory). Target: the recap under its 10,000 characters.
3. **Briefs by file** (pipes' second slice). `brief.py` writes the block to a
   file beside the note file and prints one line with its path and hash;
   `guard_spawn.py` checks the file; the `Agent` prompt is that line plus the
   task (at most 150 words). The persona reads its brief itself. Saves about
   1,450 tokens of dispatcher context per spawn; replayed, -27 % on its own.
4. **Handbacks to 150 words**, `Result` first; details stay in the commit
   message. About -10 % of growth.
5. **One status tool for merge-readiness.** `tools/ready.py <branch>` prints
   one line: conflicts with master and with queued branches, CI, last
   recorded review. It replaces the dispatcher's merge-tree, `gh` and log
   reads; content claims are `@reviewer`'s, not the dispatcher's.
6. **Shorter dispatcher sessions run by the harness: not now.** Steps 1-2
   already give a fresh context through files, inside Claude Code's own
   compaction. A driver of `claude -p` sessions would cost the prompt Ola
   talks to and a new permission path. Revisit only if, a week after steps
   1-5, the main session is still above 25 % of the day's weighted cost
   (`2026-10-25-merge-352.md` §3.2 method).

Order: 1 and 2 first (no design round; 2 is in the first slice), then 3-5 as
the pipes' second slice, one design page.

## 4. Design growth and wording rounds

Findings in `2026-10-25-merge-352.md` §3.1 and §3.3: every design doubled
between first commit and approval, each review answered by adding text; four
red/green pairs on output text after the first green took 2 h 55 min of the
triangle budget's 5 h 56 min; a 52-minute tester round reconciled pinned
predicted figures. Fewer rounds also means fewer spawns, and each spawn costs
the dispatcher about 3 k tokens (block twice, task, launch notice, handback) of context it then re-reads (§1).

## 5. Rule inventory for this round

Read in full: `CLAUDE.md`, `.claude/REQUIRED-READING.md`,
`docs/increments/README.md`, `docs/PRINCIPLES.md`, `.claude/briefs/common.md`,
the six persona files, `.claude/settings.json`, `guard_spawn.py`,
`record_review.py`, `tools/brief.py`. The four skills, the other hooks and
tools: the 10-21 inventory (`2026-10-21-research.md` §1) stands, plus the
changes since `e5c12a95` (`git diff e5c12a95 HEAD -- .claude/hooks tools
.claude/skills`), read; none touches these rows.

| # | Rule | Owner | Restated or pointed at | Conflict |
|---|---|---|---|---|
| T1 | In-flight state on disk; `session.md` format | `REQUIRED-READING.md` | `session_state.py` warns; h31 §1 | Broken in practice: 84 lines (§1) |
| T2 | Recap capped at 10,000 characters | `REQUIRED-READING.md` | none | 41,932 now |
| T3 | Spawner deletes the note file | `REQUIRED-READING.md` | none | 51 linger; no tool |
| T4 | Resume a persona only under about 150 k | `CLAUDE.md` §3 | memory note | Nothing bounds the main session |
| T5 | Spawn carries the unchanged block | `brief.py`, `guard_spawn.py` | `REQUIRED-READING.md` *The harness* | Forces the block through the dispatcher's context |
| T6 | Handback at most 300 words | `common.md` | memory note | none |
| T7 | Verify a handback: re-run its cheap claims | `PRINCIPLES.md` D2 | `reviewer.md` "verify recorded evidence" | **Conflict** with layering: puts reading in the dispatcher |
| T8 | Conflicts checked before a push question | `CLAUDE.md` §3 | none | Done by hand, no tool |
| T9 | One-page design under 150 lines | README step 1 | none | Designs double on review |
| T10 | A refusal test asserts the reason | `tester.md` §4 | none | Pins whole sentences |

| Proposal | Rows | File | Net words |
|---|---|---|---|
| Plan 1 | T4 | `.claude/settings.json` only | 0 |
| Plan 2 | T1-T3 | hook and pipes, no rule text | 0 |
| Plan 3 | T5 | `REQUIRED-READING.md`: "without an unchanged, current block from" becomes "without a current brief file from" | about -2 |
| Plan 4 | T6 | `common.md`: "300 words" becomes "150 words, `Result` first" | +2 |
| Plan 5 | T7, T8 | D2's rule becomes "Re-run the cheap claims in any report before repeating them to Ola; the main session uses a tool's status line, never diffs." | about +8 |
| merge-352 P2 | T9 | README step 1: "A design covers only the next PR, about one page per 150 lines ...; a review round's answer replaces text, it does not add it." | +14 |
| merge-352 P3 | T10 | `tester.md`: "asserts the reason's key words, never a whole sentence or a figure the design only predicts" | +7 |
| merge-352 P4 | | README step 1: cut "The option to publish is kept open ..." | -17 |
| merge-352 P5 | | cut the two "only statement of the rule" sentences | -20 |
| **Net** | | | **about -8** |

Rule text now: `python3 tools/rule_sizes.py` **9,496** words, 0 since
`88994fac`. Anecdotes in rule files: none new.

## Ola, verbatim

"Instead of me pretending to understand what's the way forward here, we have
an @orchestrator that should be very much on top of this?" and "You want _me_
to /clear your context? That's not scalable." (2026-10-25, quoted in the
brief, checked by `brief.py`). His steer on the plan, relayed by the main
session: "I my mind, layering and separation of concerns is the way forward,
but I'm not an LLM expert :)".
