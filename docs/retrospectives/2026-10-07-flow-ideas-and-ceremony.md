# Flow: three items in flight, a short path for small changes, and an ideas file

Proposal by `@orchestrator`, 2026-10-07, written on `5c412383`. **No rule file
is changed by this document.** Ola approves the wording; then the main session
briefs the persona that owns each file. Times are UTC (Ola's clock is +2 h).

What Ola said, verbatim (main session transcript
`6a989a29-4492-430f-99b0-9fe3058f931a.jsonl`, lines 23844 to 24072):

- 11:42: "I think we see regression in effieciency of the development
  process now. The number of things to do grow faster than we can finish
  them."
- 11:45, answering the main session's two questions (push guard PR B; finish
  before starting, at most 3 items in flight, product first): "1: yes 2: yes,
  this is a good idea. We still need a place to record ideas and future work."
- 11:49: "One line per idea might be quite small. Or could one line contain an
  extensive idea description? We also should have agents explicitly knowing
  about the file during planning. Is it the orchestrator or you (main) that
  need it in your descriptions?"
- 11:52: "Shall I give @orchestrator this shape to write up, together with the
  lighter-ceremony rule? Yes"

So part 1 (the cap) is **ruled**; this document only words it. Parts 2 and 3
are **proposals**.

Labels used below: **31** = increment 31, binary `.vtk` and `.ply` by default
(branch `worktree-binary-default`). **h18** = harness increment 18, the
worktree and master-merge tools (`worktree-harness-tools`). **h19** = harness
increment 19, the after-commit gate that follows the commit into its worktree
(`worktree-commit-gate`, called T1 in chat). **Guard PR B** = the second PR of
harness increment 16, closing hidden routes past the push and governance guards
(`worktree-h16b`, now #210). **20c-2** = the second PR of increment 20c, the
soft quality rule and line split (`worktree-soft-quality-2`). **30d** = growing
staircase outlines in pieces (perf audit, #208). **Pins** = choices the red
step's tests make that the design left open (`.claude/agents/tester.md@5c412383:27-29`).

## 0. Evidence from today

### What was spent, and what landed

From the subagent transcripts under `.../6a989a29.../subagents/` (first and
last timestamp of each run, worktree from its brief block), between Ola's
return (05:34) and this run's start (11:52):

- **54 spawns, 611 agent-minutes.** By kind: 13 `@reviewer`, 11 `@tester`,
  8 `@developer`, 4 `@architect` designs, **15 `@architect` runs that
  recorded a verdict, ruled pins, fixed findings or wrote the as-built
  section**, 2 `@orchestrator`, 1 `@perf`.
- **Merged today: #194, #205, #206, #207, #209.** All five were started
  before today (#194 opened 2026-10-06 08:50; #205 and #207 are 2026-10-06
  increments; #206 and #209 are retrospectives). **No item started today had
  merged by 11:52**, and the last merge was #209 at 08:37: 3 h 15 min with no
  merge while 4-5 items were open.
- Items started today: 30d's measurement (06:46), 31 (07:33), h18 (07:51),
  20c-2's red step (08:24), h19 (09:08).

| Item | Net lines (`count_loc.py`) | Design file, lines | Spawns | Design / code review rounds | Record, pin and fix spawns | Agent-min | State at 11:52 |
|---|---|---|---|---|---|---|---|
| 31 | 4 (est. about 1) | 271 | 10 | 2 / 1 | 4 (`@architect` 3, `@tester` citation fix 1) | 65 | code approved 08:42; not pushed |
| h19 | 67 (est. about 45) | 434 | 12 | 2 / 2 | 3 | 95 | approved 11:50; recording running |
| h18 | 385 (est. about 275) | 870 | 10 | 2 / 1 | 3 | 102 | 3 blockers; fix round queued |
| Guard PR B | 116 | (h16 file) | 23 since 2026-10-05, 7 today | "round 11" by the file's own count | 7 `@architect` overall | 421 overall, 82 today | pushed as #210 at 11:46 |
| 20c-2 | 130 | (20c file, 6 design rounds 2026-10-06/07) | 5 | 0 / 0 so far | 1 pin ruling | 151 | mutation round done; `@perf` and review not started |

Commands: `python3 tools/count_loc.py $(git merge-base master <branch>) <branch>`;
`git show <branch>:docs/increments/<file> | wc -l`.

The harness items (h18, h19, guard PR B) took 279 of the 611 agent-minutes,
46 %. The product items took 312.

### What the 4-line change cost (31)

Seven commits for 4 net production lines (`git log ad91b5dd~1..b7bc1210`):
design `ad91b5dd`; design round 1 answered `09cbee96`; red `e3fe391a`; design
round 2 recorded and four pins ruled `777f66d1`; a test citation fixed
`9008811c`; green `29bc5ad3`; as-built `b7bc1210`. Ten spawns, 65
agent-minutes, 69 minutes from design to code approval (07:33 to 08:42).

Then it **waited more than three hours for a one-line recording commit.** At
08:42 the main session said it was "ready to push once its approval is
recorded"; Ola said yes to the push at 08:59 (transcript line 22271). The
main session's 08:59 queue put the recording behind guard PR B's code *and*
h19's new design ("2. The T1 design. 3. Recording binary-default's
approval", line 22265), because a recording `@architect` counts as one of the
two writers. At 11:52 it was still unpushed.

What design review round 1 caught was real: two golden tests hash the text
output (reviewer handback, "Design review 31 round 2"). On a shorter path the
green step's whole-suite run would have failed on them instead, and the fix
would have gone back to `@tester`: later, not lost.

### Findings (role section 1)

1. **A new start ranked ahead of a finish.** The 08:59 queue (line 22265)
   put h19's design ahead of 31's recording, 17 minutes after Ola's 08:40
   "Nono, we work in order, no need to rush." The main session named the
   cause itself at 11:42: "I overcorrected."
2. **Recording is the spawner's job, and was delegated as a writer.**
   `docs/increments/README.md@5c412383:59-63` and `.claude/agents/reviewer.md@5c412383:46` say the spawner (the
   main session) copies the verdict. Today five `@architect` runs did only
   that ("20c-1 record cr1...", "20c-1 record survivor kills", "20c-1 record
   review round 2", "Guard PR B record r11", "h19 recording commit"), and the
   recording was then queued for a writer slot. Ola approved a tool for it this
   morning (P6 of `2026-10-07-night-20c-1.md`), not yet built.
3. **The design review round is in no rule file.** `grep -rn -i 'design
   review' CLAUDE.md .claude/ docs/increments/README.md` finds nothing. It is
   the main session's practice. Today it ran twice on each of 31, h18 and h19.
4. **A brief dropped a step the increment file requires.** h18's design
   (§5) asks `@reviewer` for one real run of `new_worktree.py`; the brief
   forbade running the tools, so round 1 skipped it (h18 code review round 1
   handback; the main session owned it at 11:37). `.claude/briefs/common.md@5c412383:5` says "a brief
   cannot drop a step".
5. **A claim to Ola overstated.** At 11:42 the main session said h18, h19 and
   guard PR B "took most of today's agent time". Measured: 46 %, the largest
   share, not most.
6. **Approved work forgotten.** #198 (rename `cli_driver.py`, approved by
   `@reviewer`) has been open since 2026-10-06 09:01 and is not enqueued.
   22 worktree branches have commits not on `origin/master`, 16 of them last
   touched before today (`git rev-list --count origin/master..<branch>` over
   `git worktree list`; not checked one by one, some may be superseded).
7. **`Monitor` is in every brief and in no persona's tool list.**
   `.claude/briefs/common.md@5c412383:12` says "Wait for a background run with the
   Monitor tool, never `sleep`"; the `tools:` line of all six persona files
   (`.claude/agents/*.md:4`) omits it. h19's round-2 reviewer waited with an
   `until ...; do sleep 5` loop instead (its Lessons). The morning check's P12
   covered `@architect` only.
8. **This run started over the cap.** The cap was adopted at 11:45 with five
   items open; at 11:46 the main session said this write-up "starts as soon as
   one item lands"; it started at 11:52, on Ola's direct yes, with none landed
   (#210 was pushed, not merged). Not a breach if Ola's direct ask overrides
   the cap; the rule below says so (question Q1).
9. **#209 was pushed with no `@reviewer` run.** `.claude/REQUIRED-READING.md@5c412383:126-128`
   requires one on every branch before its first push, prose included. No
   spawn has `worktree=.../retro-1007b` in its brief except the morning check
   itself (07:33-07:50); the push is at transcript line 22052 (08:35), on
   Ola's yes. The cap's pressure runs the other way, which is why the short
   path in part 2 keeps the review in its "never dropped" list.
10. **`session.md` carried a ruling.** At 08:22 the main session appended a
   `RULED (Ola 2026-10-07 ...)` line (transcript line 21859);
   `.claude/REQUIRED-READING.md@5c412383:34-37` allows only `NOW:`, `QUEUE:` and `ASK OLA:`
   lines, "no rulings, no history". At 11:34 it reported clearing a stale
   question from the same file.

## 1. The cap: finish before starting (ruled; wording only)

**Where it lives.** `CLAUDE.md` §3, *The main session dispatches*, because it
binds the main session only, and every dispatch rule is there. Not in
`REQUIRED-READING.md`, which every persona reads. The count is shown by the
recap tool (below).

**Diff, `CLAUDE.md` after line 59** (the end of *Step order*):

```diff
   A performance fix is timed by `@perf` before review. Failing
   tests go back to `@developer`.
+* **Finish before starting.** At most three items are in flight. An item
+  counts from its first spawn until its PR merges or Ola parks or drops it.
+  A new item starts only when one lands, or when Ola names it and is told
+  the count. A free writer goes first to the item nearest to merging, then
+  to product work before harness and process work. The next item is chosen
+  with Ola from `docs/ideas.md`.
```

**The recap.** `tools/session_state.py` prints "In flight:" for the current
branch only (`in_flight`, lines 143-151). Proposal: one line, "In flight: N
of 3", listing each worktree branch with an open PR or with commits not on
`origin/master` and touched in the last 7 days, each with its PR number. About
25 lines with tests, in a governed file; its own small item, after one of the
current five lands.

Cost: about 70 words of rule text; the recap line above.

## 2. A short path for small changes (proposal)

**Source.** Anthropic, "Building effective agents" (read 2026-10-07):
"we recommend finding the simplest solution possible, and only increasing
complexity when needed", and "Agentic systems often trade latency and cost
for better task performance, and you should consider when this tradeoff makes
sense." Little's law (lead time = work in progress / throughput) is the
standard argument behind work-in-progress caps in Kanban; cited from memory,
unchecked.

**The threshold.** Estimate under 50 net lines (`CLAUDE.md` §2's count) and
no C++ under `include/` (so no predicate, kernel, refine or mesh code, and no
`@perf` run). Today 31 qualifies; h19 (67) and h18 (385) do not. The brief
names the path; `@reviewer` may send an item back to the full path.

**What changes for 31-sized items.** From 10 spawns to 5: design, red, green,
one review, one recording. Saved on 31: two design review rounds and their fix
run, the pin ruling run, the as-built run (about 20 agent-minutes), and, the
larger part, the three-hour wait for a recording slot.

**Diff, `docs/increments/README.md`, new section after line 82** (before
`## Acceptance`):

```diff
+## Small changes: the short path
+
+An increment estimated under 50 net lines (`CLAUDE.md` §2) that touches no
+C++ under `include/` takes the short path; the brief says so, and `@reviewer`
+may send it back to the full path.
+
+- Its increment file is short: what changes, the tests that pin it, what is
+  left out, the estimate, and the *Prior art* section (step 1).
+- One `@reviewer` round judges the design, the red step and the green step
+  together, and rules on the red step's pins. A pin that changes the design
+  still goes to `@architect` before green.
+- The verdict, the Status line and the ROADMAP row land in one commit,
+  before the push.
+
+Never dropped: the red commit ahead of the green one, a green commit with no
+test file, `@reviewer` before every push, Ola's yes for each push, and the
+prompt on a governed file.
```

**Diff, `.claude/agents/tester.md@5c412383:27-29`:**

```diff
 * **Choices beyond the design:** list every choice your tests pin that the
   increment file leaves open, under the handback heading "Pinned or assumed
-  beyond the design"; `@architect` confirms or rules on each before green.
+  beyond the design"; `@architect` confirms or rules on each before green
+  (on the short path, only a pin that changes the design: `docs/increments/README.md`).
```

**Diff, `.claude/agents/reviewer.md` after line 44** (check 5):

```diff
+6. **On the short path** (`docs/increments/README.md`), the one round also
+   judges the design and rules on the red step's pins.
```

**Conflicts lifted, each a question for Ola** (Q2 to Q4 below): who writes the
recording commit; whether a recording counts as a writer; whether the full
path's design review is written down.

## 3. The ideas file, `docs/ideas.md` (proposal)

### The shape (as Ola approved it at 11:49-11:52)

```markdown
# Ideas and future work

Nothing here is in flight. An idea becomes a ROADMAP row when it starts.
Newest last. A long idea keeps its own file; its entry points to it.

## <short name>

- **What:** the idea in plain words, as long as it needs.
- **Why:** the problem it solves, or the measurement behind it.
- **From:** who, when, and a pointer (Ola's message, a review, a retrospective).
- **Size:** a rough guess, if known.
- **Status:** idea | approved by Ola, <date> | started as ROADMAP row N | dropped: <reason>
```

"approved by Ola" is added to Ola's three statuses: under the cap, approved
work waits, and today's approved-but-unstarted work (the morning check's
proposals, 30d, 20c-3) needs a place (Q5).

### Who reads it and who adds to it

| Who | Reads | Adds |
|---|---|---|
| Main session | when a slot frees, to choose the next item with Ola | Ola's asks and ideas, the moment he voices them; ideas from `@tester`, `@developer` and `@perf` handbacks; a reviewer's untaken suggestions |
| `@architect` | before every design | ideas it finds that are out of scope |
| `@reviewer` | (read-only) | its untaken suggestions, through its spawner |
| `@orchestrator` | weekly, for stale or finished entries | process ideas from retrospectives |
| `@tester`, `@developer`, `@perf` | not needed while working | in the handback, under *Ideas* |

### How they learn it: diffs

**`.claude/REQUIRED-READING.md`, after line 62:**

```diff
 Also read `docs/increments/README.md` (the protocol and a round's cost
 constraints) and the increment file for whatever you are working on.
+
+**Ideas and future work live in `docs/ideas.md`.** Read it before you plan.
+An idea outside your task goes in your handback under *Ideas*; the main
+session adds it.
```

**`.claude/briefs/common.md@5c412383:17`:** "Hand back under: **Result; Pinned or
assumed beyond the design; Questions for Ola; Lessons; Ideas; ASK OLA and
GUARD FALSE POSITIVE lines**" (adds "Ideas").

**`.claude/agents/architect.md`, after line 40** (section 5; the file has no
planning section, and section 5 is where its working rules are):

```diff
+* **Ideas:** before a design, read `docs/ideas.md`; name each entry the
+  design takes up or overlaps, and add an entry for each idea it leaves out.
```

**`.claude/agents/reviewer.md@5c412383:52`:**

```diff
-4. **Suggestions:** Non-blocking, and only where no gate would catch it.
+4. **Suggestions:** Non-blocking, and only where no gate would catch it. One
+   the branch does not take goes into `docs/ideas.md`, added by your spawner.
```

**`.claude/agents/orchestrator.md`, section 3 and section 4:**

```diff
 Text (`python3 tools/rule_sizes.py`) and proposes a cut.
+Once a week, check `docs/ideas.md` for entries gone stale or done unmarked.
 ...
-- **Write only under `docs/retrospectives/`.**
+- **Write only under `docs/retrospectives/`, and add entries to `docs/ideas.md`.**
```

**`CLAUDE.md` §3:** covered by the last sentence of the cap bullet in part 1.

**The recap:** "Open ideas: N" beside "In flight: N of 3", in the same small
`tools/session_state.py` item as part 1.

### Seed entries for the first version

Written here because `docs/ideas.md` is outside this persona's write limit;
`@architect` creates the file after Ola's yes.

#### Morning check proposals P1-P12 and the tester.md cut
- **What:** the twelve proposals of the night retrospective: work by day in parallel first, `@perf` timing at quiet times, a prose-only neighbour measured beside `@perf`, two writers kept, `brief.py` cutting at a word boundary, review records copied by a tool (P6), the mutation round placed before review, `session.md` format enforced, CI's GCC locally, worktree checks not scratch, three guard false positives, `Monitor` for `@architect`; and a cut of about 120 words of `tester.md`.
- **Why:** this morning's idle time and the brief refusals.
- **From:** `@orchestrator`, 2026-10-07; `docs/retrospectives/2026-10-07-night-20c-1.md@6ae60089:373-384`.
- **Size:** P6 about 40 lines; the rest mostly rule text.
- **Status:** approved by Ola, 2026-10-07 ("yes to the morning check proposals", 08:22).

#### 30d: grow staircase outlines in pieces
- **What:** grow a raster-traced outline in short overlapping pieces and unite them, instead of one GEOS buffer of the whole staircase; test that the region matches GEOS's within 1e-6 m.
- **Why:** 18.5 s of Lagan's run after F1a, 3.7 s in pieces; minutes or worse for outlines with steps finer than the grid.
- **From:** perf audit (#208, `docs/increments/perf-audit.md`, finding F1c); Ola chose pieces over node distance tests, 2026-10-07 06:48 ("yes to both").
- **Size:** about 25 lines.
- **Status:** approved by Ola, 2026-10-07; measured (`@perf`, worktree `buffer-speed`); no ROADMAP row.

#### 20c-3: input clean-up, coarsening, and closing gaps in the input
- **What:** 20c-3 as designed, plus closing slits between land-cover polygons that should share an edge, such as the 1 cm by 79 m slit in Lagan's CORINE data.
- **Why:** that slit alone forces a 0.0071° angle. Ola, 08:19: "I consider the gap in the CORINE-data an error. The edges should have been shared."
- **From:** Ola, 2026-10-07 08:19 and 08:22 ("build 20c-2 and 20c-3"); `docs/increments/20c-soft-quality.md`.
- **Size:** in the 20c file.
- **Status:** approved by Ola, 2026-10-07; part of ROADMAP row 20c; design update waits on 20c-2's measurements.

#### Cap on runner and shell words per line in the push guard
- **What:** deny a line with more runner and shell words than a small cap, with a plain reason, instead of the crash-worded deny at the recursion limit (about 990 words).
- **Why:** bounds the guard's time, which grows with shell words times line length.
- **From:** h16 review round 10, suggestion S1; `docs/increments/h16-harness-fixes.md` §6 on `worktree-h16b` (#210), line 1368; Ola, 10:07, "Ok, go for defaults" (no cap in PR B, a later harness increment).
- **Size:** small.
- **Status:** idea.

#### Make the after-commit hook formatter-clean
- **What:** apply the three hunks `ruff format --diff` gives for `.claude/hooks/gates_after_commit.py`.
- **Why:** `.claude/` is outside the format gate, so nothing else will.
- **From:** h19 code review round 1, S4.
- **Size:** three hunks; a governed file, so Ola's prompt.
- **Status:** idea.

#### h18: explicit capture flag, and the merge-head check at the stop
- **What:** (1) `step` in `tools/new_worktree.py` captures output only when its last argument starts with `--`; make it a keyword argument. (2) In `merge_master.py start`, compare `MERGE_HEAD` with the recorded hash before printing the conflict stop.
- **Why:** (1) a later trailing flag would silently swallow output; (2) the race refuses before the persona resolves anything.
- **From:** h18 code review round 1, suggestions "explicit-capture-flag" and "merge-head-check-at-stop".
- **Size:** about 2 lines each.
- **Status:** idea; h18's fix round may take them.

#### Report a missing shell_scan outside a work tree
- **What:** when `tools/shell_scan.py` fails to import and the command ran outside any work tree, the after-commit hook should say so.
- **Why:** today the import failure shows only on a commit made inside a tree.
- **From:** h19 code review round 2, suggestion 1 (reviewer's default: leave it).
- **Size:** a line or two in a governed file.
- **Status:** idea.

#### Monitor in every persona's tool list
- **What:** add `Monitor` to the `tools:` line of all six persona files, or make `.claude/briefs/common.md@5c412383:12` name the fallback.
- **Why:** every brief says to wait with `Monitor`, never `sleep`; no persona has it, and h19's reviewer fell back to a `sleep` loop.
- **From:** h19 code review round 2, Lessons; extends the morning check's P12 (architect only).
- **Size:** six one-word edits, governed files.
- **Status:** idea.

#### Governed-file asks answered without a prompt in auto mode
- **What:** make the governance guard's "ask" reach Ola even in auto permission mode: an explicit ask rule in the settings, or sessions that edit governed files outside auto mode.
- **Why:** when `@developer` edited h19's hook, the guard returned "ask" and no prompt reached Ola. Inferred from the session's mode, not checked.
- **From:** h19 code review round 2; main session, 11:50.
- **Size:** a settings edit; Ola's yes needed.
- **Status:** idea (a question for Ola).

#### Recap: items in flight and open ideas
- **What:** "In flight: N of 3" across worktrees, with PR numbers, and "Open ideas: N" in `tools/session_state.py`.
- **Why:** the recap shows only the current branch; #198 sat approved and unenqueued for over a day.
- **From:** this proposal, part 1.
- **Size:** about 25 lines with tests, governed.
- **Status:** idea (follows Ola's yes to this proposal).

#### Clear out stale worktrees
- **What:** list the worktrees whose branch is merged or abandoned, and remove them with Ola's yes.
- **Why:** 111 worktrees; 22 branches with commits not on master, 16 untouched since before today.
- **From:** this proposal, finding 6.
- **Size:** a listing script, or one `git worktree remove` pass.
- **Status:** idea.

## 4. Rule text: size and a matching cut

`python3 tools/rule_sizes.py` at `5c412383`: **10,564 words**, unchanged since
the night retrospective (`6ae60089`). Parts 1 to 3 add about 340 words (measured over the diff blocks above): the cap
70, the short path 150, `tester.md` 15, `reviewer.md` 40, `REQUIRED-READING.md`
30, `architect.md` 28, `orchestrator.md` 25, `common.md` 1.

**Cut C1, about 300 words, so the net is about +40:**

1. `.claude/REQUIRED-READING.md@5c412383:175-188` (149 words, the `SessionStart`
   paragraph) becomes:
   "`SessionStart` runs `tools/session_state.py` on every source; if the recap
   is missing, run it by hand (step 1). Its output is capped at 10,000
   characters (`python3 tools/session_state.py | wc -m`), so keep
   `session.md` and the note files short." (45 words; saves about 105.) What
   goes, to a retrospective: the one-time check for a recap in the first
   persona spawned with the hook live, and how the tool finds the main
   checkout, which its own code states (`main_checkout`, lines 154-157).
2. `.claude/REQUIRED-READING.md@5c412383:137-151` (149 words): the list of what
   `guard_push.py` asks about is the hook's own docstring. Replace lines
   137-143 with "`guard_push.py` asks before any act that publishes, rewrites
   history or writes refs, remotes or config (its docstring lists them);"
   (saves about 45).
3. `docs/increments/README.md@5c412383:10-20` (87 words, "Why these exist on disk") is
   reasoning, not a rule, and *Reference, do not restate* (lines 108-110)
   already states the rule. Move it to a retrospective (saves about 85).
4. `docs/increments/README.md@5c412383:75-77` and `:80-82`: drop the reasons after the
   rules ("It is a step rather than an expectation because ..."; "The whole
   protocol rests on ...") (saves about 65).

Cuts shift line numbers: run `python3 tools/check_citations.py` on the branch
that makes them, and re-read its at-risk list.

## Questions for Ola

Each has a default; the default is what the diffs above assume.

- **Q1. Does your direct ask override the cap?** Default: yes; the main
  session tells you the count when it starts a fourth item.
- **Q2. Who writes the recording commit?** The rules already give it to the
  main session (`docs/increments/README.md@5c412383:59-63`, `.claude/agents/reviewer.md@5c412383:46`), but today it went to
  `@architect` five times, after role-bleed worries about the main session
  editing design files (`next.md`, 2026-09-28 item). Default: the main session
  runs the review-copy tool (morning check P6) once built; until then it
  copies the verdict line by hand, verbatim, and touches nothing else.
- **Q3. Does a docs-only recording commit count as one of the two writers?**
  31 waited three hours for a slot. Default: no, when no other agent is writing
  in that worktree. This changes your "agents in pairs" rule, so it is yours.
- **Q4. Write the full path's design review into the rules?** It is practice,
  not rule (finding 3), so "the short path has no design review" points at
  nothing. Default: yes, one line in `docs/increments/README.md` step 1:
  "`@reviewer` reviews the design before `@tester` starts."
- **Q5. Add "approved by Ola" as a status in the ideas file?** Default: yes.
- **Q6. What happens to ROADMAP rows that are approved but not started**
  (5d, 16d, the rows "approved by Ola 2026-10-06; not designed", the parked
  backlog)? Default: they stay where they are, since moving them shifts the
  ROADMAP lines that about 45 citations quote (h18's design review counted 46
  at-risk quotations of lines 50 and 54); from now on a new unstarted item goes
  to `docs/ideas.md` only.
- **Q7. Does an open PR waiting only on you (#198) count toward the three?**
  Default: yes, until you enqueue it or park it; the recap line makes it
  visible.
- **Q8. Do retrospectives and process write-ups count toward the three?**
  Default: yes, since they use writer slots and review rounds (#206 took
  three review rounds on the night of 2026-10-06).
