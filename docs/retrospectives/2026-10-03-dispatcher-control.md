# Why the main session makes most of the mistakes, and how to take decisions off it

Research by `@orchestrator`, 2026-10-03, at master 586fbc1. Ola asked for it:
"I have registered that a lot of the problems in the harness is related to
the main loop (you). Why is that?" and then "I need even more control. My gut
feeling is to offload the recurring decisions, or make them more mechanical."

The brief asked for this file under `docs/research/`. `@orchestrator` writes
only under `docs/retrospectives/` (its persona file, section 4), so it is
here. If Ola wants it under `docs/research/`, moving it is one `git mv`.

Words used below:

- **the main session**: the Claude Code session Ola talks to. It starts the
  personas (`@architect`, `@tester` and so on) and is called the
  *dispatcher* here.
- **a brief**: the text the main session writes when it starts a persona.
- **a handback**: the report a persona returns when it is done.
- **`session.md`**: the main session's note of what is running, what is
  queued and what waits on Ola (`.claude/current-task/session.md`).
- **the recap**: what `tools/session_state.py` prints when a session starts.
- **a hook**: a script Claude Code runs at a fixed point (before a tool call,
  when a persona starts, when a turn ends) that can add text, or refuse.

Transcripts are cited as *session, line*: the session's first eight
characters and the line in its `.jsonl` file under
`~/.claude/projects/-Users-skavhaug-projects-rasputin/`. Times are UTC.

## 1. The short answer

The main session's own five reasons are mostly right; one is half wrong, and
two are missing.

- **Right:** nobody checks it; it does the glue between roles; its rules are
  prose while the personas are partly fenced by hooks.
- **Half wrong:** "personas reread their rules, I carry mine from memory".
  Persona files are re-read: Claude Code watches `.claude/agents/` and "the
  next delegation uses the updated definition, with no restart needed"
  (code.claude.com/docs/en/sub-agents). This run got the new
  `orchestrator.md`. But `CLAUDE.md` is not re-read, by the main session
  *or by the personas it starts*: this run, started at 16:04, carries the
  `CLAUDE.md` from before 09:43 ("`@orchestrator`: Master project driver")
  although the file on disk changed at 09:43. So every persona works from
  the main session's stale copy of `CLAUDE.md`; it is not only the main
  session's problem. The docs say `/compact` re-reads the project
  `CLAUDE.md` from disk (code.claude.com/docs/en/memory, "Instructions seem
  lost after `/compact`").
- **Missing, 1: volume.** By 16:04 today the main session had started 38
  personas (this run included), sent 38 follow-up messages to running ones,
  run 187 shell commands, and answered 39 turns from Ola. 76 of those shell
  commands name `session.md`, most of them hand-written `sed` substitutions
  followed by a `grep -c` to see whether they matched (af318b82, lines 575
  to 811). From 27 to 29 September one session ran 97 shell commands
  containing `git commit` and 92 naming `ROADMAP.md` (38487caf). Every persona starts with an
  empty context; the main session runs for hours or days (85c14e7c ran from
  30 September to 3 October, through two compactions). Anthropic's own
  guidance: "As the number of tokens in the context window increases, the
  model's ability to accurately recall information from that context
  decreases" (anthropic.com/engineering/effective-context-engineering-for-ai-agents,
  2025-09-29).
- **Missing, 2: its rules live in three places, one of them private.** The
  dispatcher's rules are in `CLAUDE.md` section 3, in
  `.claude/REQUIRED-READING.md`, and in about ten memory notes (agents in
  pairs, lean briefs, no implicit labels, fill the window, caffeinate,
  call him Ola, recap each round, ...). The memory notes are loaded once at
  session start, are not in the repository, and no reviewer reads them. Two
  of this week's errors cite one: the red-step briefs for 15e and 15f-2
  dropped the required mutation tests as "lean brief", going past the rule
  file and past the note's own exception (section 2, row 5).

Claude Code's documentation states the general point plainly: "Claude treats
them [`CLAUDE.md` and memory] as context, not enforced configuration. To
block an action regardless of what Claude decides, use a PreToolUse hook
instead" (code.claude.com/docs/en/memory).

## 2. The catalogue: this week's main-session errors

From 27 September to 3 October. Each row says which recurring decision it
was, and what could have stopped it: a **script** (a tool the main session
runs), a **template** (fixed text the brief must contain), a **hook**
(Claude Code runs it, so the main session cannot forget it), a **persona**,
or **nothing mechanical**.

### Briefing

| # | What happened | Evidence | Could have been stopped by |
|---|---|---|---|
| 1 | Ola's sentence about multicolouring was paraphrased in a brief; the design then printed the paraphrase as his verbatim words | 38487caf 1732 ("my brief to `@architect` paraphrased your ... sentence") | Template: Ola's words are pasted from the transcript, never retyped. `tools/session_state.py` already extracts his turns |
| 2 | `@perf` was asked to fix a `tools/bench.py` bug with no failing test first; `@perf` refused | 38487caf 9213; `next.md`, 2026-09-29 | Template: the brief cannot remove a step the persona's own file requires |
| 3 | The 16b design's mutant checks were never briefed; the main session first told Ola they had not been asked for, then corrected itself | 38487caf 7690 | Script: the brief generator copies the increment file's required checks into the brief |
| 4 | Seven `@reviewer` runs were told to write and commit their own review records, against `docs/increments/README.md` ("its spawner copies the handback's verdict") | af318b82 600: "record your round in the increment file as the rules require (commit on the branch; do not push)"; `next.md` item 8 | Template, plus the reviewer persona carrying the rule (already proposed, `next.md` item 8) |
| 5 | The mutation tests the increment files require were left out of 15e and 15f-2 ("Lean: no throwaway or mutation rounds") | af318b82 484 and 1842; `2026-10-03-day-merges-15e-15f.md` section 2b | Script: as row 3. Dropping a required step needs a flag that writes an `ASK OLA:` line instead |
| 6 | An `@architect` was told to "pick a free increment id"; it picked 24, already taken on an unmerged branch | af318b82 3902, then 4129 ("Id clash") | Script: next free id across all branches, not only master |
| 7 | This run was told to write `docs/research/`, outside `@orchestrator`'s write limit | af318b82 4428 | Template: the brief carries the persona's write limit from its file |
| 8 | A resumed `@tester` started in another agent's worktree: the main session had `cd`'d into `agent-ae452511...` at 15:54 and resumed the tester at 15:59 from there | af318b82 4254, 4334, 4405 ("my environment note pointed at `agent-ae452511f0527a114`") | Hook: refuse a spawn or resume while the main session's working directory is not the main checkout (hook input carries `cwd`) |

### Relaying and answering Ola

| # | What happened | Evidence | Could have been stopped by |
|---|---|---|---|
| 9 | Internal codes sent to Ola, four times in six days: "I don't understand what R5 is ... You have started with some strange abbrevs" (28 Sep 17:18), "Can you stop referring to things as A, B, R2?" (28 Sep 20:18), "We're slowly developing tribal language here again" (1 Oct 16:30), "What's L14? I'm human and you throwing these codes at me is silly" (today 15:46). A memory note against it was rewritten again today, after the fourth | 38487caf 7724, 8221; 85c14e7c 3250, 4767 (admitted); af318b82 4064 | Hook, partly: a check when the turn ends that lists code-like tokens (letter plus digits) not explained in the same message. A tripwire, with false positives to measure |
| 10 | 15f-1's size overrun quoted as 45 % (gross) instead of 36 % (net), copied from a handback | af318b82 1745; day retrospective section 3e | Script: one line counter that implements `CLAUDE.md` section 2's rule; numbers quoted from it or from `@reviewer`'s record |
| 11 | `@tester` flagged a test assumption the ruling had not made; it went to Ola as "worth knowing" and to nobody to confirm; it was the assertion later found wrong | day retrospective section 3c | Template plus routing: a handback has a fixed "Assumptions beyond the ruling" heading, and anything under it is queued for `@architect` |
| 12 | A cause relayed before it was checked: "`@perf` was debugging a bug in its simulation", wrong (the crashes were a planted bug) | 38487caf 3689 | Nothing mechanical; the claims rule in `REQUIRED-READING.md` |
| 13 | "About one line of code at no real cost", wrong | 38487caf 8405 | Nothing mechanical |
| 14 | "No access to the artifact upload tool", wrong | 7954d3a8 1889 | Nothing mechanical |
| 15 | Told Ola he could start unattended mode from the session with `!`; it cannot (`ENXIO`) | 85c14e7c 9003, 9008 | Nothing mechanical; the claims rule |
| 16 | "Restarting is safe; the merge scripts die on restart": they do not, and two merge chains raced | af318b82 172; `2026-10-03-night-restart-h7.md` section 2 | Script: a merge runner that holds a lock and refuses a second chain |
| 17 | Read "later today" as "next" and reordered the day | 85c14e7c 1480 | Nothing mechanical |
| 18 | Read "stick to the plan" as a new order; Ola: "The plan was 21d first" | 38487caf 3458 to 3496 | Script, partly: a queue file in a fixed order, changed only by a named command, makes a reordering visible |
| 19 | Drafted a feedback report to Anthropic unasked, quoting Ola, then put the first draft's content into the second | 85c14e7c 4843, 4903 | Nothing mechanical; a judgement error |

### Queueing, concurrency, merging, state, restarts

| # | What happened | Evidence | Could have been stopped by |
|---|---|---|---|
| 20 | Three long unattended windows with long idle stretches: 3 h 24 min, 5 h 59 min, 7 h 08 min | `2026-10-03-night-restart-h7.md` section 1a; 85c14e7c 1345, 10473 (both admitted) | Hook: while unattended mode is on, refuse to end a turn when nothing is running and the queue has an item that needs no ruling. Already ruled: `next.md` item 1 |
| 21 | The first fallback queued for today's window (a green step) writes three governed files, which unattended mode refuses | day retrospective section 4 | Script: a queue line names the paths it will write; the recap flags governed ones (`next.md` item 17) |
| 22 | Three agents at once several times, and three in one worktree, against the "agents in pairs" note | day retrospective section 2d | Hook: count running personas (recorded when each starts and stops) and refuse a third |
| 23 | Merge scripts rewritten by hand in the scratchpad each session (`merge.sh`, `merge2.sh`, `merge7.sh`); one read "no checks reported" as green | af318b82 2034; `next.md`, findings of 2026-10-01/02 | Script: one checked-in, tested merge runner |
| 24 | `ROADMAP.md` conflicts resolved by hand twice; stale "awaiting the push" text merged to master | day retrospective sections 2e, 3d | Script: the status column built from increment files (`next.md` item 14) |
| 25 | `session.md` grew into a 29-line log; five decisions written as bullets under one `ASK OLA:` heading were invisible to the recap | `2026-10-03-night-restart-h7.md` sections 1c, 1d | Script: a writer for `session.md` in the format Ola ruled today (NOW, QUEUE, one `ASK OLA:` line per decision) |
| 26 | No `@orchestrator` check after five merges, until Ola asked | af318b82 3601, 3606; day retrospective section 2a | Script: the recap prints merges since the last retrospective (`next.md` item 9b) |
| 27 | A rule merged at 09:43 never reached the session or its 33 personas | day retrospective section 2a; this run's own context | Hook: after a merge or pull that changes `CLAUDE.md`, say "run `/compact` before the next spawn" |

### Role

| # | What happened | Evidence | Could have been stopped by |
|---|---|---|---|
| 28 | About to write `@orchestrator`'s reflection document itself; Ola: "You're not the @orchestrator" | 7954d3a8 1678, 1685; `2026-09-29-orchestrator-and-hooks-audit.md` section 1.2 | Hook: the role-limits design merged today (`docs/increments/h6-role-limits.md`, not yet built) gives the main session a row that does not include `docs/retrospectives/` |

### What the catalogue says

28 errors in seven days. By kind:

| Recurring decision | Errors | Rows |
|---|---|---|
| Writing a brief | 8 | 1 to 8 |
| Relaying and answering Ola | 11 | 9 to 19 |
| Queueing and idle time | 2 | 20, 21 |
| Concurrency | 1 (repeated) | 22 |
| Merging and `ROADMAP.md` | 3 | 16, 23, 24 |
| State files and triggers | 2 | 25, 26 |
| Picking up rule changes | 1 | 27 |
| Role | 1 | 28 |

(Row 16 is counted under both relaying and merging.)

- **22 of 28 could have been stopped or caught by a script, template or
  hook** (two of them, rows 9 and 18, only partly). All eight briefing
  errors are in that group, and they cost the most: two increments merged
  without their required mutation tests, seven review records written by the
  wrong persona, and an increment number collision.
- **6 cannot be made mechanical** (rows 12 to 15, 17 and 19): four claims
  made without trying the thing first, one misreading of Ola, and one
  judgement error. They need the claims rule in `REQUIRED-READING.md`, and
  Ola's eye.
- **The glue work** (rows 16, 20 to 27) is all mechanical in kind. None of it
  needs judgement; it needs bookkeeping that never forgets.

## 3. The recurring decisions, and who should own each

"Owner" is the most mechanical thing that can do the job. The main session
keeps only what needs judgement.

| Decision | Today | Proposed owner | What stays with the main session |
|---|---|---|---|
| **What step comes next on each branch** | Remembered, written into `session.md` by hand | Script `tools/pipeline.py`: reads git and each increment file (design, red and green commits, `## Review` rounds, `@perf` acceptance, required mutation suite, `ROADMAP.md` row) and prints one line per branch: "15f-3: green step next" | Choosing between branches when several are ready, inside the order Ola set |
| **Writing a brief** | Hand-written, 900 to 2,800 characters each, paraphrasing rules from memory | Script `tools/brief.py <persona> <worktree> <step>` prints the fixed part: files to read from disk (including `CLAUDE.md`), the worktree and isolation, the output file under `.claude/current-task/`, the persona's write limit, the steps the increment file requires, the handback headings, and one sentence: "if this brief contradicts your persona file or the increment file, they win; say so in your handback". A hook refuses a spawn whose brief lacks the generated part | The task itself, in a few sentences, and Ola's words pasted, not retyped |
| **Who may run now** | A memory note ("pairs") | Hook: records each persona's start and stop, refuses a third, refuses a second persona in a worktree that is busy, refuses anything beside `@perf` while it times | Asking Ola for an exception |
| **Increment numbers** | The designing persona guesses | `tools/pipeline.py next-id`: the next id free across all branches | none |
| **Merging** | A scratchpad script rewritten each session | Checked-in `tools/merge_queue.py`: one chain at a time (lock), "no checks reported" is not green, update branch, merge commit only, stops on a conflict and says which file | Getting Ola's word for each PR; resolving a conflict that is not mechanical |
| **`ROADMAP.md` status** | Hand-edited in every PR | Built from the increment files' status lines (`next.md` item 14, option b) | none |
| **Recording a review** | The main session appends by hand (sometimes the reviewer did, wrongly) | Script `tools/record_review.py` appends the round from the handback's fixed fields (verdict, range, line count, blocking items) | none |
| **`session.md`** | 76 hand edits today | Script `tools/session.py now / queue / ask / answered`, writing the format Ola ruled today; the recap warns on anything else | Wording the ASK for Ola |
| **When to run `@orchestrator`** | Remembered; missed five times today | The recap prints "N merges since the last retrospective" (`next.md` item 9b) | none |
| **Routing what a handback raises** | Read and relayed by the main session | Handbacks end in fixed headings: Assumptions beyond the ruling, Questions for Ola, Lessons. A script (later a hook when a persona stops) turns them into queue lines: assumptions to `@architect`, questions to `ASK OLA:` lines, lessons to `@orchestrator` | Translating the questions into plain words for Ola |
| **Filling an unattended window** | A memory note | Hook (ruled in principle, `next.md` item 1): while unattended mode is on, a turn may not end with nothing running and a decision-free item queued | Building the queue before Ola leaves |
| **Picking up a rule change** | Not done | Hook after a merge or pull that changes `CLAUDE.md`: "run `/compact` before the next spawn" | Running it |
| **Line counts and figures for Ola** | Copied from handbacks | Script: one line counter for `CLAUDE.md` section 2's rule | Choosing what Ola needs to hear |
| **Plain language to Ola** | A memory note | Hook when a turn ends: list unexplained codes (warn first, count false positives) | Writing plainly |
| **Settling a conflict between two rules** | The main session, silently (lean briefs against mutation testing) | Ola. The brief script makes it visible: leaving out a required step needs `--skip <step> --why`, which writes an `ASK OLA:` line | none |
| **Pushes, merges, rulings, priorities** | Ola | Ola (unchanged) | Asking, one decision per line |

**What stays judgement in the main session, and why.** Three things:
turning Ola's words into a task (it needs the conversation); explaining
results to Ola in plain words (it needs to know what he knows); and handling
the unexpected (a crash, a persona that died, a finding that changes the
plan). Everything else in the table is bookkeeping, and bookkeeping is what
a long context does worst and a script does best.

## 4. What others do

Each source was read on 2026-10-03 unless marked.

1. **Use fixed code paths where the task is known.** Anthropic, "Building
   effective agents" (Erik Schluntz and Barry Zhang, 2024-12-19): "Workflows
   are systems where LLMs and tools are orchestrated through predefined code
   paths. Agents are systems where LLMs dynamically direct their own
   processes"; "Workflows offer predictability and consistency for
   well-defined tasks". The orchestrator-workers pattern there "suits complex
   tasks where subtasks can't be predicted in advance". Our pipeline (design,
   red, green, review, `@perf`) *is* predicted in advance. The same article:
   "You can add programmatic checks ('gate' ...) on any intermediate steps".
   Applies to: `tools/pipeline.py`, the brief script.
2. **Fixed pipelines beat free agents on software tasks.** Xia, Deng, Dunn
   and Zhang, "Agentless" (arXiv 2407.01489, 2024): "a simplistic three-phase
   process of localization, repair, and patch validation, without letting
   the LLM decide future actions" outperformed the open-source agents on
   SWE-bench Lite at $0.70 per issue. Wu et al., "StateFlow" (arXiv
   2403.11322, 2024) separates "process grounding" (which state the task is
   in, driven by a state machine) from "sub-task solving" (the LLM's work in
   that state), and reports 13 % and 28 % higher success than ReAct at a
   fifth and a third of the cost. Applies to: the split in section 3,
   where the state machine says what comes next and the personas do the work.
3. **Practice says the same.** HumanLayer's "12-Factor Agents"
   (github.com/humanlayer/12-factor-agents): "Most of the products out there
   billing themselves as 'AI Agents' are not all that agentic. A lot of them
   are mostly deterministic code, with LLM steps sprinkled in at just the
   right points." Its factors include "Own your control flow", "Unify
   execution state and business state" and "Small, Focused Agents" (titles
   only checked; the bodies were not read).
4. **Delegation needs a fixed shape.** Anthropic, "How we built our
   multi-agent research system" (2025-06-13): "Each subagent needs an
   objective, an output format, guidance on the tools and sources to use, and
   clear task boundaries"; "Without detailed task descriptions, agents
   duplicate work, leave gaps, or fail to find necessary information"; and
   "We combine the adaptability of AI agents built on Claude with
   deterministic safeguards like retry logic and regular checkpoints".
   Applies to: the brief script and the handback headings.
5. **State goes in files, and structured files survive better than prose.**
   Anthropic, "Effective harnesses for long-running agents" (2025-11-26):
   agents keep a progress file and a feature list, edited "only by changing
   the status of a passes field", because "the model is less likely to
   inappropriately change or overwrite JSON files compared to Markdown
   files"; and work goes "one feature at a time". Anthropic, "Effective
   context engineering for AI agents" (2025-09-29): "structured note-taking,
   or agentic memory", and sub-agents with "clean context windows" returning
   condensed summaries. Applies to: `session.md` written by a script, not by
   `sed`.
6. **Most multi-agent failures are structural.** Cemri et al., "Why Do
   Multi-Agent LLM Systems Fail?" (arXiv 2503.13657, 2025): 14 failure modes
   in three groups (system design, inter-agent misalignment, task
   verification), and the failures "require more sophisticated solutions"
   than better prompts. Already cited in `next.md` for role bleed.
7. **Claude Code's own mechanisms** (code.claude.com/docs):
   - `CLAUDE.md` and memory are "context, not enforced configuration. To
     block an action regardless of what Claude decides, use a PreToolUse
     hook" (memory page).
   - Hook matchers use the exact tool names, and the tool that starts a
     persona is `Agent`; resuming one is `SendMessage` (tools reference). So a
     `PreToolUse` hook can check, rewrite (`updatedInput`) or refuse a brief.
   - `SubagentStart` fires when a persona starts and can add text to its
     context; it cannot refuse. `Stop` and `SubagentStop` can refuse the end
     of a turn with a reason (hooks page).
   - Persona files reload without a restart (sub-agents page); `CLAUDE.md`
     is re-read after `/compact` (memory page).
   - `claude --agent <name>` runs the main session itself as a named agent,
     whose system prompt replaces the default one; hooks then see its
     `agent_type` (sub-agents and hooks pages). This would make the
     dispatcher a named role with its own file, as the outside review of
     2026-10-01 asked (`next.md`, "An external review", point 1). Whether
     that file also reloads without a restart for the main session is not
     documented; unchecked.
   - Not checked: whether `SubagentStop` fires for personas run in the
     background, and whether a background persona's completion counts as a
     prompt for `UserPromptSubmit`. Both need a probe before a hook relies
     on them.

The common thread: the coordinator keeps the steps that need judgement, and
everything with a known answer moves into code that runs every time.

## 5. A staged plan

Ola rules each stage. Line counts are estimates of tool code, tests not
counted. Everything that touches `.claude/settings.json` or a hook needs
Ola's yes and goes through the normal design, red, green and review steps as
a harness increment.

Already ruled today and assumed below: `next.md` items 1 (idle is accepted
only if the fallback list is truly empty), 2, 3, 5 and 6, and item 4 as
option (a): `session.md` holds NOW, QUEUE, and one `ASK OLA:` line per open
decision, nothing else.

**Stage 1: briefs come from files (the first step).** The largest group of
costly errors, rows 1 to 8.

- `tools/brief.py`, about 100 lines: prints the fixed part of a brief from a
  short template per persona plus the increment file (required steps, the
  invariant-critical suite if named, review rounds so far), the worktree,
  the output file, and the persona's write limit.
- A `PreToolUse` hook on `Agent` and `SendMessage`, about 40 lines: refuses
  a spawn whose prompt lacks the generated block, and refuses any spawn or
  resume while the main session's working directory is not the main
  checkout (row 8).
- Six templates of about 15 lines of prose each, under `.claude/`.
- In return, delete `CLAUDE.md` section 3's "Briefs" bullet and the "lean
  briefs" memory note: the templates replace both.
- Cost: one harness increment, about 140 lines. Would have stopped rows 2
  to 8, and row 1 if the template says Ola's words are pasted.

**Stage 2: one place for state, and "what comes next".** Rows 6, 21, 25, 26.

- `tools/session.py`, about 80 lines: `now`, `queue add|done`, `ask`,
  `answered`; writes the format Ola ruled. This is how `next.md` items 3 and
  4 get implemented, not an extra.
- `tools/pipeline.py`, about 150 lines: one line per open branch saying
  which step is next and which required step is missing; merges since the
  last retrospective; the next free increment id; queue items that write
  governed paths. The recap prints it.
- Cost: one or two harness increments, about 230 lines.

**Stage 3: merges by a checked-in runner.** Rows 16, 23, and the restart
race.

- `tools/merge_queue.py`, about 120 lines: a lock so only one chain runs;
  "no checks reported" is not green; update the branch, wait, merge commit
  only; stop on a conflict and name the file. `next.md` item 5 (list the
  running background jobs in the recap) reads its lock.
- Cost: one harness increment.

**Stage 4: hooks on the dispatcher's own turns.** Rows 9, 20, 22, 27. Probes
first (does `SubagentStop` fire for background personas; how often the code
check fires wrongly on a day's messages).

- Concurrency: `SubagentStart` and `SubagentStop` record who runs; the
  stage 1 hook refuses a third persona, a second one in a busy worktree, or
  anything beside a `@perf` timing run. About 50 lines.
- When a turn ends (`Stop`): in unattended mode, refuse to stop with
  nothing running and a decision-free item queued (`next.md` item 1); at any
  time, list unexplained codes in the message, as a warning only for the
  first week. About 60 lines.
- After a merge or pull that changes `CLAUDE.md`: "run `/compact` before the
  next spawn". About 15 lines.
- Cost: one harness increment, about 125 lines, and one week of measured
  false positives before the code check may refuse.

**Stage 5: handbacks routed, reviews recorded, `ROADMAP.md` generated.**
Rows 10, 11, 24.

- Fixed handback headings in every template (stage 1), and
  `tools/record_review.py` (about 50 lines) to append a review round.
- `ROADMAP.md`'s status column built from the increment files (`next.md`
  item 14, option b, about 60 lines).
- A line counter for `CLAUDE.md` section 2 (about 60 lines), so a size
  figure has one source.
- Later, once stage 4 shows `SubagentStop` works for background personas:
  route the handback headings into the queue automatically.

**Stage 6, not now: a driver outside the chat.** A Python state machine that
runs the personas headless and leaves the chat session only to talk to
Ola. It is what the sources in section 4 point to at the end of the road,
but it is weeks of work and would replace much of what exists. Decide after
stages 1 to 4 have run for two weeks and the count in section 2 shows what
is left.

**Not re-proposed: a dispatcher persona file.** Ola ruled today that the
main session is the dispatcher by having no agent name, with no
`.claude/agents/dispatcher.md` (`docs/increments/h6-role-limits.md`,
section 7, question 1). One fact for a later look, not a request: persona
files reload without a restart and `CLAUDE.md` does not, so the
dispatcher's rules would stay fresh in a persona file, if a main session
started with `claude --agent` reloads its file too. That last part is
undocumented and unchecked.

**How to tell whether it works.** `@orchestrator`'s check after each merge
and each unattended night counts main-session errors by the kinds in
section 2. This week's baseline is 28, of which 22 were mechanical in kind.
After stages 1 to 3 the mechanical kinds should be near zero; if they are
not, the tools are not being used, and stage 4's hooks are the next step.

**What the plan does not add.** The outside review of 2026-10-01 warned
against more governance; each stage above replaces prose with a tool and
deletes the prose (`CLAUDE.md`'s brief rules, `session.md` prose in
`REQUIRED-READING.md`, and the memory notes that hold dispatcher rules:
pairs, lean briefs, fill the window). The rule count should go down, not up.

## 6. Ola's rulings needed

1. Stage 1 (brief script and spawn hook) as the next harness increment?
2. Stages 2 and 3 after it, in that order?
3. Stage 4: probe first, then decide?
4. Should this file move to `docs/research/`?
