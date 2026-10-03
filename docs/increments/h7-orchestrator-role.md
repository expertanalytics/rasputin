# Harness h7: what @orchestrator does now

Status: **design**, @architect, 2026-10-03. Rules and docs only, no code. One
PR. Implements Ola's rulings of 2026-10-03, which also answer questions 1 and
2 of `docs/increments/h6-role-limits.md`.

## 1. Prior art

Tooling; no novelty claimed. The split is the usual one between the process
that runs the work and the process that watches it: an operator and an
auditor. The failure it answers is role bleed, "disobey role specification"
(failure mode 1.2 in Cemri et al. 2025, arXiv 2503.13657), already cited in
`docs/retrospectives/next.md`. *Legacy*: nothing.

```
$ git grep -l -iE 'orchestrator|retrospective' legacy-archive -- legacy
(no output, exit 1)
```

## 2. What changes

**Who dispatches.** The main session, started with no agent name, is the
dispatcher ("the main loop"). It runs the step order, briefs the personas,
recaps each round and reports to Ola. Those rules move out of
`.claude/agents/orchestrator.md` into `CLAUDE.md` §3, under "The main session
dispatches". `CLAUDE.md` is loaded into every session, so the dispatcher
always has them in context. They live there and nowhere else.

**What @orchestrator does.** Three jobs, all in
`.claude/agents/orchestrator.md`:

1. *Watch the workflow.* Read the transcripts' traces (commits, increment
   files, handbacks quoted in them, the harness queue, `gh pr checks`) and
   flag deviations: a skipped step, a persona doing another persona's work, a
   guard worked around, a review going past two rounds, an unattended window
   left idle, a result reported that was never checked.
2. *Research agentic design.* Papers, Anthropic guidance, other harnesses.
   Each proposal carries its source and the incident in this repo it would
   have prevented.
3. *Own the retrospectives.* `docs/retrospectives/`, `next.md` included.
   @orchestrator is the only writer there. Other personas put a lesson in
   their handback; the dispatcher passes it on; @orchestrator records it.

**Its write limit.** @orchestrator writes only under `docs/retrospectives/`.
A change to a rule file (`CLAUDE.md`, `.claude/**`, `REQUIRED-READING.md`,
persona files, hooks) it proposes in its report; Ola decides, and the
dispatcher briefs the persona that owns the file.

The persona file states what the role does, and keeps a limit only where it
is a boundary the agent could plausibly cross and that can be enforced; the
write limit is the one such limit.

## 3. When it runs

After each increment merges, the morning after each unattended night, and
for a research round weekly or when Ola asks. The schedule is stated once,
in `CLAUDE.md` §3, *The main session dispatches*, where the dispatcher reads
it.

A check produces a short report to the dispatcher and the entries it writes
in `docs/retrospectives/`.

## 4. How the write limit is enforced

**Now:** rule text only. `.claude/agents/orchestrator.md` states the limit,
and `guard_governance.py` already asks (by day) or refuses (unattended)
before any write to a rule file, whoever makes it. Nothing yet stops a write
to, say, `src_python/`; the frontmatter grants `Write` and `Edit` because
recording a retrospective needs them, and Claude Code's `tools:` key cannot
scope a tool to a path.

**Later, with h6:** the `ROLES` table in `state_check.py` gets the row
`orchestrator` → `docs/retrospectives/`. From then on a write elsewhere is
denied (Edit, Write) or reported (Bash), like any other persona's. The
`architect` and `dispatcher` rows, which held all of `docs/`, now exclude
`docs/retrospectives/`, so @orchestrator is the only writer there; h6's test
list gains test 14 for both, the dispatcher's write being an ask by day and a
queued deny at night (h6 §3.1). The `ROLES` lookup gains exclusions inside an
allowed prefix for this. The h6 design file is updated in this PR.

## 5. Files

| File | Change |
|---|---|
| `.claude/agents/orchestrator.md` | rewritten for the three jobs; description and tools |
| `CLAUDE.md` | header, §1 roster line, §3 pipeline no longer "via @orchestrator"; new "The main session dispatches" with the rules moved from the persona (the milestone-update rule included, as its own bullet), the lessons chain and the @orchestrator schedule |
| `docs/PRINCIPLES.md` | "write it in the log" becomes "report it in your handback", pointing to `CLAUDE.md` §3 |
| `docs/increments/h6-role-limits.md` | questions 1 and 2 ruled; `orchestrator` row; `docs/retrospectives/` taken out of the `architect` and `dispatcher` rows; test 14; exclusions in the `ROLES` lookup, estimate ~125 → ~130 |
| this file | the design |

No persona file or `REQUIRED-READING.md` told a persona to write a
retrospective (`git grep -n -i -E 'retrospective|next\.md' -- .claude CLAUDE.md`
finds only citations), so none needs redirecting. The one rule that did is
`docs/PRINCIPLES.md`'s "When you learn something, write it in the log".

Estimate: 0 production lines.

## 6. Rulings

1. **Principle D3** (`docs/PRINCIPLES.md`), "a compromise between two
   correct principles goes to `@orchestrator`": ruled by Ola, 2026-10-03.
   D3 stays unchanged; @orchestrator proposes, Ola rules.

## Review

**Round 1**, @reviewer, on 0d15560: **CHANGES REQUESTED**. Blocking:
the milestone-update rule was dropped in the move to `CLAUDE.md`; h6's
test 14 lacked the dispatcher case that §4 claimed; the lessons rule was
stated in full in three places. Non-blocking: put the @orchestrator schedule
in the dispatcher rules, once; say that h6's `ROLES` lookup supports an
exclusion. All five fixed in b6af28e, and D3 marked ruled.

**Round 2**, @reviewer, on b6af28e: **APPROVED**. 0 production lines.
Remaining non-blocking note: when h6 lands, its persona-file pointer to the
`ROLES` table replaces the write-limit line in `.claude/agents/orchestrator.md`.
