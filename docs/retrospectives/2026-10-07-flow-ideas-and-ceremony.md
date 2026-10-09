# Flow: three items in flight, a short path for small changes, and an ideas file

**Status: approved in review round 2; next the push, on Ola's yes.** Ruled by Ola on
2026-10-07 (11:45, 12:40, 13:20 and 13:36, below), except the recap line of
part 1, which waits on Ola. Proposal by `@orchestrator`, written on
`5c412383` (`c61c1aa0`), amended to Ola's rulings the same day (`db51001c`)
and to review round 1. **No rule file is changed by this document**; the main
session briefs the persona that owns each file. Times are UTC (Ola's clock is
+2 h). Transcript line numbers count from 1, as `grep -n` and `sed -n` do.

Ola's rulings, verbatim (the main session showed him the questions at the end
renumbered; the mapping is the main session's):

- 12:40: "1: Yes, but after approving a warning 2: Yes, you do it. 3: No 4:
  Yes 5: Yes". His 1 to 5 are Q1, Q2, Q3, Q7 and Q8 below.
- 13:20: "yes, yes, yes. 1, 2, 3." His 1, 2, 3 are Q4, Q5 and Q6, each with
  its default.
- 13:36: "Ok, yes to both then." (transcript line 24753), answering the
  main session's two questions at 13:34 (line 24670): Q9, the short path
  itself, and Q10, the reading of his "3: No". Both with their defaults.

What each ruling changes is under *Rulings* at the end; the diffs in parts 1
to 3 already follow them.

What Ola said, verbatim (main session transcript
`6a989a29-4492-430f-99b0-9fe3058f931a.jsonl`, lines 23845 to 24073):

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

So part 1 (the cap) was approved at 11:45 and this document words it; its
details are Q1, Q7 and Q8 (12:40). The 11:52 yes was to *write up* the short
path, not to the rule: part 2 was approved at 13:36 (Q9), its details being
Q2 to Q4. Part 3's shape was given at 11:49-11:52 and its details ruled at
13:20 (Q5, Q6). The recap line (part 1) was not put to Ola: the cap bullet
does not name it, so it waits on Ola.

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
return (05:34) and this proposal's start (11:52), this proposal's own run left
out:

- **53 spawns, 609 agent-minutes.** By kind: 13 `@reviewer`, 11 `@tester`,
  8 `@developer`, 4 `@architect` designs, **15 `@architect` runs that
  recorded a verdict, ruled pins, fixed findings or wrote the as-built
  section**, 1 `@orchestrator` (the morning check), 1 `@perf`.
- **Merged today: #194, #205, #206, #207, #209.** Four were started before
  today (#194 opened 2026-10-06 08:50; #205 and #207 are 2026-10-06
  increments; #206 is a 2026-10-06 retrospective). **#209, the morning
  check, was started today** (spawned 07:33, its one commit `6ae60089` at
  07:49) and merged at 08:37 with no `@reviewer` run (finding 9). **No other
  item started today had merged by 11:52**, and #209's was the last merge:
  3 h 15 min with no merge while 4-5 items were open.
- Items started today: 30d's measurement (06:46), the morning check (#209,
  07:33), 31 (07:33), h18 (07:51), 20c-2's red step (08:24), h19 (09:08).

| Item | Net lines (`count_loc.py`) | Design file, lines | Spawns | Design / code review rounds | Record, pin and fix spawns | Agent-min | State at 11:52 |
|---|---|---|---|---|---|---|---|
| 31 | 4 (est. about 1) | 271 | 10 | 2 / 1 | 4 (`@architect` 3, `@tester` citation fix 1) | 65 | code approved 08:42; not pushed |
| h19 | 67 (est. about 45) | 434 | 12 | 2 / 2 | 3 | 95 | approved 11:50; recording running |
| h18 | 385 (est. about 275) | 870 | 10 | 2 / 1 | 3 | 102 | 3 blockers; fix round queued |
| Guard PR B | 116 | (h16 file) | 23 since 2026-10-05, 7 today | "round 11" by the file's own count | 7 `@architect` overall | 421 overall, 82 today | pushed as #210 at 11:46 |
| 20c-2 | 130 | (20c file, 6 design rounds 2026-10-06/07) | 5 | 0 / 0 so far | 1 pin ruling | 151 | mutation round done; `@perf` and review not started |

Commands: `python3 tools/count_loc.py $(git merge-base master <branch>) <branch>`;
`git show <branch>:docs/increments/<file> | wc -l`.

After 11:52 (`gh pr view <n> --json mergedAt`): h19 merged as #211 at 12:19.
31's verdict was recorded by the main session itself, on Ola's 12:40 ruling
(`69e3551b`, "Recorded by the main session, per Ola 2026-10-07"), and 31
merged as #212 at 13:04, 25 minutes after the ruling. #198 and #210 are still
open.

The harness items (h18, h19, guard PR B) took 279 of the 609 agent-minutes,
46 %. The product items (31, 20c-1's last runs before its merge as #207, 20c-2, the perf audit
and 30d's measurement) took 312, and the morning check most of the rest
(16.5).

### What the 4-line change cost (31)

Seven commits for 4 net production lines (`git log ad91b5dd~1..b7bc1210`):
design `ad91b5dd`; design round 1 answered `09cbee96`; red `e3fe391a`; design
round 2 recorded and four pins ruled `777f66d1`; a test citation fixed
`9008811c`; green `29bc5ad3`; as-built `b7bc1210`. Ten spawns, 65
agent-minutes, 69 minutes from design to code approval (07:33 to 08:42).

Then it **waited more than three hours for a one-line recording commit.** At
08:42 the main session said it was "ready to push once its approval is
recorded"; Ola said yes to the push at 08:59 (transcript line 22272). The
main session's 08:59 queue put the recording behind guard PR B's code *and*
h19's new design ("2. The T1 design. 3. Recording binary-default's
approval", line 22266), because a recording `@architect` counts as one of the
two writers. At 11:52 it was still unpushed.

What design review round 1 caught was real: two golden tests hash the text
output (reviewer handback, "Design review 31 round 2"). On a shorter path the
green step's whole-suite run would have failed on them instead, and the fix
would have gone back to `@tester`: later, not lost.

### Findings (role section 1)

1. **A new start ranked ahead of a finish.** The 08:59 queue (line 22266)
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
   (#210 was pushed, not merged). Ola ruled at 12:40 that his direct ask
   overrides the cap only after a warning he approves (Q1); the rule below
   says so. This run had no such warning.
9. **#209 was pushed with no `@reviewer` run.** `.claude/REQUIRED-READING.md@5c412383:126-128`
   requires one on every branch before its first push, prose included. No
   spawn has `worktree=.../retro-1007b` in its brief except the morning check
   itself (07:33-07:50); the push is at transcript line 22053 (08:35), on
   Ola's yes. The cap's pressure runs the other way, which is why the short
   path in part 2 keeps the review in its "never dropped" list.
10. **`session.md` carried a ruling.** At 08:22 the main session appended a
   `RULED (Ola 2026-10-07 ...)` line (transcript line 21860);
   `.claude/REQUIRED-READING.md@5c412383:34-37` allows only `NOW:`, `QUEUE:` and `ASK OLA:`
   lines, "no rulings, no history". At 11:34 it reported clearing a stale
   question from the same file.

## 1. The cap: finish before starting (ruled)

**Where it lives.** `CLAUDE.md` §3, *The main session dispatches*, because it
binds the main session only, and every dispatch rule is there. Not in
`REQUIRED-READING.md`, which every persona reads. The recap tool would show
the count (below; waiting on Ola).

**Diff, `CLAUDE.md` after line 59** (the end of *Step order*):

```diff
   A performance fix is timed by `@perf` before review. Failing
   tests go back to `@developer`.
+* **Finish before starting.** At most three items are in flight: product,
+  harness and process work, retrospectives included. An item counts from
+  its first spawn until its PR merges or Ola parks or drops it; a PR waiting
+  only on Ola still counts. A new item starts only when one lands. If Ola
+  asks for one beyond that, warn him with the count ("this makes 4 of 3")
+  and start it only after he confirms. A free writer goes first to the item
+  nearest to merging, then to product work before harness and process work.
+  The next item is chosen with Ola from `docs/ideas.md`.
```

Rulings in it: Ola's direct ask overrides the cap only after a warning he
approves (Q1); a PR waiting only on Ola counts (Q7); retrospectives and
process write-ups count (Q8).

**The recap.** `tools/session_state.py` prints "In flight:" for the current
branch only (`in_flight`, lines 143-151). Proposal: one line, "In flight: N
of 3", listing each worktree branch with an open PR or with commits not on
`origin/master` and touched in the last 7 days, each with its PR number. About
25 lines with tests, in a governed file; its own small item, counted toward
the three. **Waiting on Ola:** the cap bullet does not name the recap line,
so his yes to the cap does not cover it, and it was not put to him on its
own.

Cost: about 105 words of rule text; the recap line above, if approved.

## 2. A short path for small changes, and who records (ruled 13:36)

**Source.** Anthropic, "Building effective agents" (read 2026-10-07):
"we recommend finding the simplest solution possible, and only increasing
complexity when needed", and "Agentic systems often trade latency and cost
for better task performance, and you should consider when this tradeoff makes
sense." Little's law (lead time = work in progress / throughput) is the
standard argument behind work-in-progress caps in Kanban; cited from memory,
unchecked.

**The threshold.** Estimate under 50 net lines (`CLAUDE.md` §2's count),
touching no C++ file and nothing the *Acceptance* section of
`docs/increments/README.md` covers (`@5c412383:86-88`: refine or mesh code
"and what drives them", which is what needs a `@perf` run). C++ lives outside
`include/` too (`src/predicates/detria_exact.cpp`,
`src/cdt/detria_backend.cpp`, `bindings/core.cpp`), so "no C++ under
`include/`", the first wording, left a hole. Today 31 qualifies: its
production change is Python only (`src_python/tin_engine/cli.py` and the
`.ply` and `.vtk` writers, `git diff --stat ad91b5dd~1 b7bc1210`), the output
format, not the refine or mesh path; h19 (67) and h18 (385) do not. The brief
names the path; `@reviewer` may send an item back to the full path.

**What changes for 31-sized items.** From 10 spawns to 4: design, red, green,
one review; the main session writes the recording commit (Q2). Saved on 31:
two design review rounds and their fix run, the pin ruling run, the as-built
run (about 20 agent-minutes), and, the larger part, the three-hour wait for a
recording slot.

**Diff, `docs/increments/README.md`, new section after line 82** (before
`## Acceptance`):

```diff
+## Small changes: the short path
+
+An increment estimated under 50 net lines (`CLAUDE.md` §2) that touches no
+C++ file and nothing the *Acceptance* section covers takes the short path;
+the brief says so, and `@reviewer` may send it back to the full path.
+
+- Its increment file is short: what changes, the tests that pin it, what is
+  left out, the estimate, and the *Prior art* section (step 1).
+- One `@reviewer` round judges the design, the red step and the green step
+  together, and rules on the red step's pins. A pin that changes the design
+  still goes to `@architect` before green.
+- The main session commits the verdict, the Status line and the ROADMAP
+  row together, before the push.
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

**`.claude/agents/reviewer.md`: no edit.** The first version added a check 6
saying the same as the short-path section's second bullet; `@reviewer` reads
`docs/increments/README.md` already (`.claude/REQUIRED-READING.md@5c412383:61-62`)
and the brief names the path, so the line was a second copy. Dropped in this
amendment to keep the budget (part 4).

**The full path's design review, written down (Q4).** Diff,
`docs/increments/README.md@5c412383:50-52`, the end of step 1 (it runs from
line 26 to line 52; the file indents step text by three spaces):

```diff
    — before, because a suite written against re-derived intent pins the
    re-derivation, and a domain constant guessed wrong is then guarded by a test
    that agrees with it.
+   `@reviewer` reviews the design before `@tester` starts.
```

The short path's "one `@reviewer` round judges the design ..." is then the
named exception to this line.

**Who records the verdict (Q2, Q3).** The rule already gives it to the
spawner; Ola ruled that the main session does it itself (Q2), and that an
agent spawned only to record still counts as one of the two writers (Q3 at
12:40, its reading confirmed as Q10 at 13:36).
Diff, `docs/increments/README.md@5c412383:59-63` (step 4; three-space indent,
as in the file):

```diff
    **The review leaves a trace in the increment file.** `@reviewer` is
-   read-only, so its spawner copies the handback's verdict, the commit range
+   read-only, so the main session itself, not an agent it spawns, copies the
+   handback's verdict, the commit range
    it reviewed and its LOC count, verbatim, into a `## Review` section of
```

`.claude/agents/reviewer.md@5c412383:46` ("your spawner records it") stays
true and needs no edit. `.claude/REQUIRED-READING.md` gets no line: it
already sends every reader to `docs/increments/README.md`
(`.claude/REQUIRED-READING.md@5c412383:61-62`), so step 4 stays the one
statement of the recording rule, and a second copy would drift from it.
Until the review-copy tool (the morning check's P6) is built, the main
session copies the verdict line by hand, verbatim; it did so for 31 in
`69e3551b`.

## 3. The ideas file, `docs/ideas.md` (ruled)

### The shape (Ola approved it at 11:49-11:52; Q5 and Q6 ruled at 13:20)

```markdown
# Ideas and future work

Nothing here is in flight. An idea becomes a ROADMAP row when it starts.
Newest last. A long idea keeps its own file; its entry points to it.
ROADMAP rows approved before this file existed stay in ROADMAP; every new
idea, approved or not, comes here only.

## <short name>

- **What:** the idea in plain words, as long as it needs.
- **Why:** the problem it solves, or the measurement behind it.
- **From:** who, when, and a pointer (Ola's message, a review, a retrospective).
- **Size:** a rough guess, if known.
- **Status:** idea | approved by Ola, <date> | started as ROADMAP row N | dropped: <reason>
```

"approved by Ola" is added to Ola's three statuses (Q5, ruled yes): under
the cap, approved work waits, and today's approved-but-unstarted work (the
morning check's proposals, 30d, 20c-3) needs a place. The ROADMAP sentence in
the header is Q6, ruled with its default: the old approved-but-unstarted rows
(5d, 16d, the rows "approved by Ola 2026-10-06; not designed", the parked
backlog) stay where they are, because moving them would shift ROADMAP lines
that about 45 citations quote.

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
 text (`python3 tools/rule_sizes.py`) and proposes a cut.
+Once a week, check `docs/ideas.md` for entries gone stale or done unmarked.
 ...
-- **Write only under `docs/retrospectives/`.**
+- **Write only under `docs/retrospectives/`, and add entries to `docs/ideas.md`.**
```

**`CLAUDE.md` §3:** covered by the last sentence of the cap bullet in part 1.

**The recap:** "Open ideas: N" beside "In flight: N of 3", in the same small
`tools/session_state.py` item as part 1; waiting on Ola, as that item is.

### Seed entries for the first version

Written here because `docs/ideas.md` is outside this persona's write limit;
`@architect` creates the file after Ola's yes.

#### Morning check proposals P1-P12 and the tester.md cut
- **What:** the twelve proposals of the night retrospective: P1 parallel work first and `@perf` timing at quiet times; P2 staging the night before Ola leaves; P3 a quiet neighbour measured beside `@perf`; P4 two writers kept; P5 `brief.py` cutting at a word boundary; P6 review records copied by a tool; P7 the mutation round placed in the step list; P8 `session.md`'s format enforced; P9 CI's GCC locally; P10 worktree checks, scratch exempt; P11 three guard false positives; P12 `Monitor` for `@architect`; and a cut of about 120 words of `tester.md`.
- **Why:** this morning's idle time and the brief refusals.
- **From:** `@orchestrator`, 2026-10-07; `docs/retrospectives/2026-10-07-night-20c-1.md@6ae60089:373-384`.
- **Size:** P6 about 40 lines; the rest mostly rule text.
- **Status:** approved by Ola, 2026-10-07 ("yes to the morning check proposals", 08:22).

#### 30d: grow staircase outlines in pieces
- **What:** grow a raster-traced outline in short overlapping pieces and unite them, instead of one GEOS buffer of the whole staircase; test that the region matches GEOS's within 1e-6 m.
- **Why:** 18.5 s of Lagan's run after F1a, 3.7 s in pieces (one run each, from the audit); minutes or worse for outlines with steps finer than the grid.
- **From:** perf audit (#208, `docs/increments/perf-audit.md`, finding F1c); Ola chose pieces over node distance tests, 2026-10-07 06:48 ("yes to both").
- **Size:** about 25 lines.
- **Status:** approved by Ola, 2026-10-07; not measured: only the performance review's one-run figures (`docs/increments/perf-audit.md@1da4a144:182`, PR #208); the `@perf` run (worktree `buffer-speed`, started 06:46, silent after 06:51, interrupted 07:19, no handback) was lost to sleep; no ROADMAP row.

#### 20c-3: input clean-up, coarsening, and closing gaps in the input
- **What:** 20c-3 as designed, plus closing slits between land-cover polygons that should share an edge, such as the 1 cm by 79 m slit in Lagan's CORINE data.
- **Why:** that slit alone forces a 0.0071° angle. Ola, 08:19: "I consider the cap [gap] in the CORINE-data an error. The edges should have been shared."
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
- **What:** apply the two hunks `ruff format --diff` gives (on master `89e35797`, with the repository's ruff settings) for `.claude/hooks/gates_after_commit.py`.
- **Why:** `.claude/` is outside the format gate, so nothing else will.
- **From:** h19 code review round 1, S4.
- **Size:** two hunks; a governed file, so Ola's prompt.
- **Status:** idea.

#### h18: explicit capture flag, and the merge-head check at the stop
- **What:** (1) `step` in `tools/new_worktree.py` captures output only when its last argument starts with `--`; make it a keyword argument. (2) In `merge_master.py start`, compare `MERGE_HEAD` with the recorded hash before printing the conflict stop.
- **Why:** (1) a later trailing flag would silently swallow output; (2) the race refuses before the persona resolves anything.
- **From:** h18 code review round 1, suggestions "explicit-capture-flag" and "merge-head-check-at-stop".
- **Size:** about 2 lines each.
- **Status:** idea; h18's fix round may take them.

#### After-commit hook: keep the cause when two faults meet
- **What:** when `tools/shell_scan.py` cannot be imported and the session's directory is in no work tree, the after-commit hook drops the note that says why; keep it. The hook still prints NOT CHECKED naming the directory and exits 2; only the cause is lost.
- **Why:** h19's design (§8 pin 7) asks for the cause; the fix for review suggestion S2 drops it in this one case (`committed_trees` returns the notes only when every entry is a `Path`).
- **From:** h19 code review round 2, suggestion 1, left as noted (`docs/increments/h19-commit-gate-tree.md@2f84131e:421-426`); h19 merged as #211.
- **Size:** a line or two in a governed file.
- **Status:** idea.

#### Push guard: PowerShell flags beyond -c and -Command
- **What:** the guards read a PowerShell (`pwsh`) command string given with `-c` or `-Command`; make them read the other flags that pass one too. Which flags, unchecked here: #210's design file names them.
- **Why:** a command string passed by another flag is a route past the push and governance guards.
- **From:** the main session's brief for this amendment, 2026-10-07 ("being fixed in #210 now"); #210 is guard PR B, branch `worktree-h16b`.
- **Size:** small; in `tools/shell_scan.py`, a governed file.
- **Status:** started in #210 (guard PR B, harness increment 16), in progress.

#### #198: rename `tests/python/cli_driver.py` to `cli_helpers.py`
- **What:** the rename, as in the open PR #198.
- **Why:** Ola ruled on 2026-10-05 that the old name was cryptic (PR body).
- **From:** opened 2026-10-06 09:01; approved by `@reviewer` (PR body: `5c34390`, and the master merge `39f61b2`).
- **Size:** a rename, no production lines.
- **Status:** approved by `@reviewer`, waiting on Ola's yes to enqueue; counts toward the three in flight (Q7).

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
- **Status:** waiting on Ola: not put to him on its own, and the cap bullet he approved does not name it. Its own item, counted toward the three, if approved.

#### Clear out stale worktrees
- **What:** list the worktrees whose branch is merged or abandoned, and remove them with Ola's yes.
- **Why:** 111 worktrees; 22 branches with commits not on master, 16 untouched since before today.
- **From:** this proposal, finding 6.
- **Size:** a listing script, or one `git worktree remove` pass.
- **Status:** idea.

## 4. Rule text: size and a matching cut

`python3 tools/rule_sizes.py` at `5c412383`: **10,564 words**, unchanged since
the night retrospective (`6ae60089`). Parts 1 to 3, as ruled, add about **370
words** (the `+` lines less the `-` lines of the diff blocks above): the cap
105, the short path 155, `tester.md` 12, the design-review line 7, the
recording line 7, `REQUIRED-READING.md` 28, `architect.md` 26, `reviewer.md`
13, `orchestrator.md` 17, `common.md` 1. The first version counted about 340;
the rulings added the warning step and the two counting rules to the cap, the
design-review line and the recording line, and dropped `reviewer.md`'s check 6.
Review round 1's fixes added 5 to the short path and took 35 off cuts 2 and 3.

**Cut C1, about 285 words, so the net is about +85:**

1. `.claude/REQUIRED-READING.md@5c412383:175-188` (149 words, the `SessionStart`
   paragraph) becomes:
   "`SessionStart` runs `tools/session_state.py` on every source; if the recap
   is missing, run it by hand (step 1). Its output is capped at 10,000
   characters (`python3 tools/session_state.py | wc -m`), so keep
   `session.md` and the note files short." (45 words; saves about 105.) What
   goes, to a retrospective: the one-time check for a recap in the first
   persona spawned with the hook live, and how the tool finds the main
   checkout, which its own code states (`main_checkout`, `tools/session_state.py@5c412383:153-156`).
2. `.claude/REQUIRED-READING.md@5c412383:137-143` (62 words): the list of
   what `guard_push.py` asks about is in the hook's own code
   (`.claude/hooks/guard_push.py@5c412383:44-45` and `:107-152`; its
   docstring does not name rebase, `reset --hard`, filter-branch,
   `commit --amend`, `gh release` or `gh repo`). Lines 137-143 become, with
   line 144 onward unchanged:
   "Active in `.claude/settings.json`: `guard_push.py` asks before any act
   that publishes, rewrites history or writes refs, remotes or git config,
   and before `--no-verify` (its code lists them; not `gh pr close` or
   `gh pr comment`); `guard_governance.py` asks" (36 words; saves about 26).
3. `docs/increments/README.md@5c412383:10-20` ("Why these exist on disk") is
   mostly reasoning, and *Reference, do not restate* (lines 108-110) restates
   lines 15-16. Two parts are stated nowhere else and stay: line 14, "These
   files are the source of truth" (`git grep` over the rule files), and lines
   18-20, the README's only pointer to `.claude/REQUIRED-READING.md` for what
   is currently being asked. So lines 10-20 go, and lines 108-110 become:
   "**Reference, do not restate.** The increment files are the source of
   truth, and `.claude/REQUIRED-READING.md` rules on what is currently being
   asked. Prompts point at them. If a fact is wrong there, fix the file rather
   than correcting it in a prompt — the correction is otherwise lost with the
   transcript." (122 words become 50; saves about 70.)
4. `docs/increments/README.md@5c412383:75-77` and `:80-82`: drop the reasons after the
   rules ("It is a step rather than an expectation because ..."; "The whole
   protocol rests on ...") (saves about 65).
5. `.claude/REQUIRED-READING.md@5c412383:191-194`: drop "Proposed and not
   approved: ..." to the end of the paragraph, a status list of two hooks
   that belongs in a retrospective; the audit it cites keeps them (saves 18).

Cuts shift line numbers: run `python3 tools/check_citations.py` on the branch
that makes them, and re-read its at-risk list.

## Rulings

Ola's answers to the questions of the first version (`c61c1aa0`, Q1 to Q8)
and to the two asked after review round 1 (Q9, Q10), with what each changed
above.

- **Q1. Does Ola's direct ask override the cap?** "Yes, but after approving a
  warning." The main session warns with the count ("this makes 4 of 3") and
  starts the item only after Ola confirms. In the cap bullet (part 1).
- **Q2. Who writes the recording commit?** "Yes, you do it": the main session
  itself, not an agent it spawns. In step 4 of `docs/increments/README.md`
  (part 2); first done for 31 in `69e3551b`.
- **Q3. Does a docs-only recording commit by an agent count as one of the two
  writers?** "No". The main session read it as "no to the proposed
  exception: it still counts", and asked again at 13:34 (Q10); Ola
  confirmed that reading at 13:36. Moot in practice once the main session
  records (Q2); Ola's "agents in pairs" rule is unchanged.
- **Q4. Write the full path's design review into the rules?** Yes: one line
  at the end of step 1 of `docs/increments/README.md` (part 2).
- **Q5. "approved by Ola" as a status in the ideas file?** Yes (part 3).
- **Q6. ROADMAP rows approved but not started?** The default: they stay in
  ROADMAP; a new idea goes only to `docs/ideas.md` (the file's header, part 3).
- **Q7. Does a PR waiting only on Ola count toward the three?** Yes (the cap
  bullet). #198 is such a PR today.
- **Q8. Do retrospectives and process write-ups count toward the three?**
  Yes (the cap bullet).
- **Q9. Does Ola approve the short path itself?** Asked at 13:34 after review
  round 1: "Changes under 50 lines that touch no C++ and nothing `@perf`
  times get one review round covering design, tests and code together, and
  I record the verdict. Failing tests first, the review, your yes for every
  push and the hook prompts all stay." Default yes. "Ok, yes to both then."
  (13:36). Part 2's heading and the Status line say so.
- **Q10. Did "3: No" mean that a docs-only recording by an agent still
  counts as one of the two writers?** Default "Yes, it still counts". Yes, in
  the same 13:36 message (Q3 above).
- **Not ruled:** the recap line of part 1 (`tools/session_state.py`, "In
  flight: N of 3" and "Open ideas: N"). It waits on Ola.

## To apply

The main session briefs the owner of each file: `CLAUDE.md` (the cap),
`docs/increments/README.md` (the short path, the design-review line, the
recording line, cuts 3 and 4), `.claude/agents/tester.md`,
`.claude/agents/architect.md`, `.claude/agents/orchestrator.md`,
`.claude/briefs/common.md`, `.claude/REQUIRED-READING.md` (the ideas line,
cuts 1, 2 and 5), and the new `docs/ideas.md` with the seed entries of part 3.
Every one of these is a governed file, so each edit meets Ola's prompt. The
recap line (`tools/session_state.py`) is not ruled; if Ola approves it, it is
its own small item, counted toward the three.

## Review

Review round 1 (@reviewer, 2026-10-07, 5c412383..db51001c, docs only): CHANGES REQUESTED; 5 blocking: 30d-not-measured, 209-started-today, short-path-cpp-hole, short-path-not-ruled, cut2-docstring; 11 suggestions: line-numbers-0-based, q3-reading, window-and-product-split, cap-word-count, ruff-hunks, morning-check-list, corine-quote, citation-offsets, step1-placement, orchestrator-context, cut3-keep-source-of-truth.

Review round 2 (@reviewer, 2026-10-07, db51001c..d482ff8c, docs only, 0 net production lines): APPROVED; 0 blocking; 2 suggestions: harness-minutes-rounding, cut2-forge-writes.

Recorded by the main session (Ola, 2026-10-07: the main session records review verdicts). Both suggestions left as noted: harness-minutes-rounding is cosmetic (279.8 minutes); cut2-forge-writes (whether cut 2 names `gh api`/`curl` writes to the forge) is for whoever applies cut 2.

Review round 1 of the applied rule edits (@reviewer, 2026-10-08, 4cd7e050..de723bbc, rule text and docs only, 0 net production lines; rule text +63 words against the +85 accepted): CHANGES REQUESTED; 2 blocking: cut2-unknown-under-runner (the new REQUIRED-READING harness paragraph says the push guard asks before an unknown git/gh command also when another program runs it; `caffeinate git myalias` and `find -exec git myalias` pass with no prompt, /Users/skavhaug/projects/rasputin/.claude/hooks/guard_push.py@4cd7e050:318), design-review-cases (the new README step 1 line asks for timing "on a case with land cover and one without"; the approved proposal 2 says "both catchments with default flags"); 2 suggestions: ideas-path-absolute, ci-after-push.

Review round 2 of the applied rule edits (@reviewer, 2026-10-08, e090a5fe..d57a4750, rule text and docs only, 0 net production lines; rule text +84 words against the +85 accepted): CHANGES REQUESTED; round 1's cut2-unknown-under-runner and design-review-cases fixed; 1 blocking: inventory-docstrings-only (the new research-round paragraph in /Users/skavhaug/projects/rasputin/.claude/worktrees/rules-213/.claude/agents/orchestrator.md@d57a4750:36-40 reads hooks and tools only by their docstrings, but Ola approved reading every harness file in full; the push guard's docstring at /Users/skavhaug/projects/rasputin/.claude/worktrees/rules-213/.claude/hooks/guard_push.py@d57a4750:30 claims more than the code enforces); 2 suggestions: name-the-catchments, docstring-g7-overclaim.

Review round 3 of the applied rule edits (@reviewer, 2026-10-08, 529db4f3..1171ca7b, and the branch as a whole against base 4cd7e050; rule text and docs only, 0 net production lines; rule text +85 words against the +85 accepted): APPROVED; round 2's inventory-docstrings-only blocker fixed (/Users/skavhaug/projects/rasputin/.claude/worktrees/rules-213/.claude/agents/orchestrator.md@1171ca7b:36-38 reads the rule files, `.claude/settings.json` and every hook and tool in full, code included); 0 blocking; 1 suggestion: quick-check-cases-wording (reviewer.md check 6 says "the quick check's catchments" while cases.toml holds two catchments and two tiles); CI has not run (no PR yet) and must be green before the PR is enqueued.

Docstring fix review, rounds 1 and 2 (@reviewer, 2026-10-09, short path; guard_push.py G7 sentence, closes round 2's suggestion "docstring-g7-overclaim"; 742165dc..f56abd7e, module docstring only, 0 net production lines). Round 1 CHANGES REQUESTED: the first wording (348248b3) said any runner skips G2, but shell_scan.py strips env, timeout, xargs and other wrappers, so G2 still applies to them. Round 2 APPROVED (f56abd7e): matches runs() and main() in /Users/skavhaug/projects/rasputin/.claude/worktrees/guard-doc/.claude/hooks/guard_push.py@f56abd7e and the wrapper list at /Users/skavhaug/projects/rasputin/.claude/worktrees/guard-doc/tools/shell_scan.py@f56abd7e:49; 16 probes agree. Suggestion: wrappers are stripped only at the head of a command (caffeinate env git myalias passes). CI must be green before enqueue.
