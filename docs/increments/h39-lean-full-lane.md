# Harness h39: a leaner full lane

Status: **approved**, 2026-10-29 (Ola, section 5), designed off master
`43f1b2e9`, applied after merging master `e0aafd73`. Rule text only: no
tool, no production line, no path `tools/check_lane.py` names. The rule
texts in section 4 are applied as written.

## Ola, verbatim

2026-10-29 06:55 UTC: "yes, we need to lean up the development process!"

2026-10-29 10:32 UTC, to section 5's five questions: "yes to all five"

## 0. The plain picture

Increment 111 (one land-cover map, branch `worktree-land-cover-edges`) has
run 17 hours from its first design commit and is still in code review. Its
record shows two loops that went round more than they needed to: the
design was reviewed six times, and code review found the same kind of bad
input (a polygon whose ring crosses itself, a "bowtie") four times, one
code path at a time. This design cuts the first loop to one round and
makes the second one name every path the first time.

## 1. Prior art: legacy and literature

*Literature:* no published method is used and nothing is claimed new.

*Legacy:* `git grep -n -i "design review\|code review" legacy-archive --
legacy` finds nothing (0 lines); the legacy code has no process text.
In this repository: h34 (`docs/increments/h34-lean-process.md`, merged as
#389) moved small changes into the fast lane and dropped the "short path"
(one combined review round for a small full-lane change). The full lane's
own loop was left as it was; this increment changes that loop.

## 2. Where 111's time went

Source: the branch's commit times (`git log --format='%h %ad %p | %s'
master..worktree-land-cover-edges`) and its `## Review` section, to code
review round 4's record `edaa4622`. A review round's time is from the
commit it reviews to the commit that records it.

| Stage | Commit times | Elapsed |
|---|---|---|
| Design, 6 review rounds | 14:33 to 18:49 | 4 h 16 min |
| of which: review rounds and their fixes | 14:33-15:03, 17:53-18:49 | 1 h 26 min |
| of which: round 2 to the revised design (Ola's answers, the revision) | 15:03-17:53 | 2 h 50 min |
| Red and 30 test-migration commits | 18:49 to 20:45 | 1 h 56 min |
| Green | to 21:37 | 52 min |
| Six test amendments | to 21:53 | 16 min |
| `@perf` quick check (the design's own section-6 gate) | to 23:23 | 1 h 30 min |
| Code review rounds 1 to 4 and their fixes | 23:23-01:10, 07:01-07:59 | 2 h 45 min |
| Overnight wait on Ola's re-pin ruling | 01:10-07:01 | 5 h 51 min |

Findings, by what they cost:

1. **Design rounds 2 to 6.** Reviewer time 3, 5, 4, 16 and 4 min: 32 min.
   Rounds 3 to 6 ran because the design changed after Ola's answers
   (README step 1: "again when an update adds a step"). Their findings: a
   false premise in a question for Ola (round 3), several land-cover files
   per run dropped without saying (round 4), the wrong overlap predicate
   and test 10's figures computed from it (round 5). In 107 to 110 the
   rounds after the first found wordings or one pin each, and took minutes
   (108's rounds 2 and 3: 11:26 to 11:32; 109's two rounds: 16:48 to
   16:58).
2. **The bowtie, four code-review rounds.** The design moved the polygon
   validity check out of the read for speed (section 4 of 111, read
   probe: 3.8 s per read, 90 % of it `is_valid`). Master's reader repaired
   every polygon on read, so every later step met valid polygons; on the
   branch, the overlap join, the `--prefer` cut and the edges clean-up meet
   unchecked ones. Code review found that class in round 1 (an overlapping
   pair), round 2 (a polygon in no pair), round 3 (outside the read box,
   in the overlap join) and round 4 (outside the read box, in the
   `--prefer` cut, `preferred.py:70`). Each fix covered the path named.
   Rounds 2 to 4 took 1 h 56 min of active time, and round 5 is still to
   come. Round 3 also found `THREADS = 8` ignoring `--threads` and the
   re-pin evidence missing; the rest of rounds 2 to 4 was this class.
3. **The test interface came after red started.** Section 8 of 111 ("the
   Python interface and wordings the test migration needs") was committed
   at 19:38, after migration commits from 19:06. Migration commits stop
   at 19:27 and resume at 19:52. The rule says "before `@developer`
   starts", so it was not late by the rule.

Checked and not proposed:

- **Master merges.** The branch has none:
  `git rev-list --merges --count 75ab2b4f..worktree-land-cover-edges` is
  0, and its reflog holds only its creation from `origin/master` and
  commits (70 now, not 67). Nothing to cut.
- **Review-record commits** (10 on 111). Each lands within minutes of the
  handback it records; no time to save shown.
- **The `@perf` quick check.** 111's own section 6 set the speed gate it
  ran; the README's acceptance rule did not ask for it. Speed is a stated
  goal, so it stays the design's call.

## 3. Design

Three changes, each paid for in words (section 4):

- **One design review round** (finding 1). `@reviewer` reviews the design
  once, before `@tester` starts, after Ola's answers are in it. Its
  findings are fixed by `@architect` in the design; code review checks
  them as part of the behaviour it already checks. A design found wrong at
  its core (the reviewer says so in its blocking issues) is redesigned and
  reviewed again. A change after the review goes through the tester's
  amendment or the code review, as a pin after red does today. Kept as
  they are: red before green, the green commit touching no test file,
  code review of the behaviour, `@perf`'s acceptance run.
- **A blocking finding names its input class and every path that input
  reaches** (finding 2). In 111 round 1 that would have read: "a
  self-crossing polygon reaches the overlap join, the `--prefer` cut and
  the edges clean-up, inside and outside the read box"; `@tester` then
  plants it on each.
- **The interface tests call is in the design before `@tester` starts**
  (finding 3), not only before `@developer`.

Savings, each an estimate from 111's commit times:

| Change | Saving in 111 | Elsewhere |
|---|---|---|
| One design round | about 30 min of review time, and 5 record commits; the fixes of rounds 3 to 5 (27 min) move to red or code review rather than vanish | one or two rounds of a few minutes per increment (107 to 110: 5 rounds) |
| Input class and every path | about 1.5 to 2 h of code-review rounds 2 to 4, and round 5 | any crash class found in review |
| Interface before red | up to 25 min of a stalled migration (19:27 to 19:52) | increments that migrate many suites |

The largest single cost in 111, the 5 h 51 min overnight wait, was for a
ruling of Ola's and is not a process step.

Risk of the first change: 111's round-4 and round-5 findings would have
reached red or code review instead. Round 5's wrong figures would have
failed at red, where `@tester` measures them; round 4's dropped behaviour
is what code review compares against master. Round 3's false premise was
in a question for Ola, which the main session checks before relaying it.

## 4. Proposed rule texts

Net words: 0 (`len(text.split())`, what `tools/rule_sizes.py` counts, on
each old and new text below; scratch script, not kept).

| # | File | Change | Words |
|---|---|---|---|
| R1 | `docs/increments/README.md` step 1 | one design round | +9 |
| R2 | `.claude/agents/reviewer.md` handback item 3 | input class, every path | +11 |
| R3 | `.claude/agents/architect.md` §5 Blueprint | before `@tester` | +3 |
| C1 | `.claude/agents/reviewer.md`: "The ceiling" section | deleted: CI already refuses over 700 lines (`.github/workflows/main.yaml:251`, `count_loc.py --max 700`) | -18 |
| C2 | `.claude/agents/reviewer.md` §5 item 6 | deleted: "once" in R1 says it | -5 |

**R1.** Old:

> `@reviewer` reviews the design before `@tester` starts, and again when an
> update adds a step.

New:

> `@reviewer` reviews the design once, before `@tester` starts; code review
> checks its findings fixed; only a design wrong at its core is reviewed
> again.

**R2.** Old: `3. **Blocking Issues.**` New:

> 3. **Blocking Issues,** each naming the failing input's class and every
>    path it reaches.

**R3.** Old:

> * **Blueprint:** a complex task gets the data flow and interface boundaries
>   before `@developer` starts.

New:

> * **Blueprint:** a complex task gets the data flow and the interface its
>   tests call before `@tester` starts.

**C1.** Delete:

> ## 1. The ceiling
> Reject a PR over `CLAUDE.md` §2's ceiling (`tools/count_loc.py`) and say how
> to split it.

**C2.** Delete from item 6: "All holes in one pass."

No tool changes: the verdict words stay `APPROVED`, `APPROVED WITH
WORDINGS` and `CHANGES REQUESTED` (`tools/push_grant.py:44` parses them).

## 5. Questions for Ola

1. Review a full-lane design once, with its findings fixed in the design
   and checked at code review, and a second review only when the design is
   wrong at its core? Default: yes.
2. Hold that one review until your answers to the design's questions are
   in it? Default: yes.
3. Have the reviewer name, for each blocking bug, the kind of input that
   breaks and every place in the code it reaches, so one fix covers them
   all? Default: yes.
4. Put the functions the tests will call in the design before the tester
   starts? Default: yes.
5. Delete the reviewer's size-limit section, since CI already refuses a
   pull request over 700 lines? Default: yes.

Answers (2026-10-29 10:32 UTC, "yes to all five"): 1 yes, 2 yes, 3 yes,
4 yes, 5 yes.

## 6. Size, tests, left out

0 counted lines. No test: the change is rule text. Left out: the review
record's format, master merges and the `@perf` quick check (section 2,
"Checked and not proposed").

Ideas file: overlaps "Lead time from design to merge, per increment"
(approved 2026-10-08, not started); section 2's table is that figure for
one increment, by hand. "Repair the land cover once per area of interest"
would remove finding 2's input class at its source; not taken up here.

## Review

- Code review round 1: CHANGES REQUESTED, e0aafd73..ee4dac92, 0 (rule text only, net words 0; R1, R2, R3, C1, C2 match section 4; citations, lane, rule-text suites pass, 1282 passed; reviewer.md kept as §5 is fine; blocking: Ola's yes to question 2, the one review waits for his answers, is in no rule text, so R1 must say it; C2's deletion of "All holes in one pass." rests on "once", which counts rounds and does not ask for every hole, so fold "naming every hole" into R1, about +6 words to pay; section 4 cites `tools/push_grant.py:44`, the verdict pattern is at line 56).
