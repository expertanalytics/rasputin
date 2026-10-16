# 2026-10-04: the check after 23c-1 (#164), and the day's lessons

Recorded by `@orchestrator` on 2026-10-04, on master at e56a2bf. Part 1 is
the post-merge check that `CLAUDE.md` §3 asks for after each merge. Part 2
records the twelve lessons the main session passed on from today's
handbacks, each with its source. Part 3 lists the proposals for Ola. Part 4
measures the rule text.

Times are UTC unless marked "local" (UTC+2). Commit times printed by git
are local. Transcripts are under
`~/.claude/projects/-Users-skavhaug-projects-rasputin/`: the night session
is `af318b82…` (af318b82 below), the midday one `965f70c2…`, and the one
that spawned this run `38de7f93…`, a fork of `965f70c2` made at 11:47:39.
"Line n" is a line of that JSONL file.

## 1. 23c-1 against the rules

23c-1 is the first half of increment 23c: the partition rule
(`decompose.py`), the piece vocabulary (`features.PIECE_VOCABULARY`), the
mesh index model (`io/mesh_index.py`) and `_core.indexed_mesh`. It merged
as #164 at e56a2bf, with a merge commit whose parents are 72d4608 and
03fc584.

### What went by the rules

- **The order of steps held.** The red commits came first: 58fb1f1, 6dd4ffe
  and 714b0b3 (`@tester`, cherry-picked from the 23c red step at 02:43,
  af318b82 line 10413). Then `@architect` settled the red step's 15 pins
  before any code (6ca9a8c). The green step followed (3ffae13, `@developer`).
  `@developer`'s three points went to `@architect` (6ce8e57), then to a new
  red commit (b9f154a) and a green one (1426245). Round 1's findings went
  through the same order again: c4ec906 and 9633282, red 7e5fe0a and
  e5881a5, green 0679950. This is retro item 13 of 2026-10-03, done before it
  was a rule: the red step's choices went to `@architect` before the green
  step.
- **No green commit touches a test file.** 3ffae13, 1426245 and 0679950 touch
  only `src_python/`, `bindings/`, `project_structure.md` and one citation in
  `docs/increments/15f-edge-strip.md` (`git show --name-only`).
- **The line count fits.** `@reviewer` counted 156 net lines against an
  estimate of about 215. My own rough count of the merge diff gives 161: a
  script that drops blank lines, comments, docstrings and raw-literal bodies,
  not the reviewer's method. The difference is in `bindings/core.cpp` (35
  against 30) and `decompose.py` (55 against 53). Either count is far under
  700.
- **No `@perf` run was needed, and that holds.** The merge diff changes
  nothing under `include/` (`git diff --stat e56a2bf^1 e56a2bf`).
- **CI was green at the merged head.** At 03fc584, all nine jobs pass
  (`gh pr checks 164`). The merge commit keeps the red-before-green history.
- **Each review round was recorded in the increment file.** Rounds 1 and 2
  were recorded by the main session (af318b82 lines 10656 and 10806), round
  3 by `@architect` on the main session's brief (42c5f53;
  `965f70c2` line 1098). That is within the README's rule, which leaves the
  copying to the spawner.
- **Every push had Ola's yes.** At 04:44:18 Ola wrote "push all five and
  merge on green" (af318b82 line 11404), and 027e8f9 was pushed at 04:45:59.
  At 11:36:43 he queued "Push #164" (`965f70c2` line 1209, absorbed mid-turn),
  and 42c5f53 was pushed at 11:36:50. `guard_push.py` asked both times
  (line 1242), and asked again before the merge job's
  `gh pr update-branch` and `gh pr merge` (line 1253).
- **`check_citations.py` exits 0 on master at e56a2bf.**

### Findings

**F1. The merge left ROADMAP row 23 and the increment's status line
stale.** `docs/increments/README.md@e56a2bf:73` says "The merge updates
`ROADMAP.md`'s row for that increment, in the same PR, before the merge".
On master at e56a2bf, ROADMAP row 23 still says 23b is "in review", though
23b merged as #162 at 4f56551. It also says 23c-1 is "built at 156 lines and
in review". The status line of `docs/increments/23-basin-scale.md` says
"23c-1 (about 215 estimated, 156 built) in review". So two merges in a row
skipped the step. This is lesson 4 below, now on master.

**F2. The round-3 brief said the branch was not pushed.** The brief at
11:11:28 (`38de7f93` line 758) said: "The branch is not pushed yet, so CI
is not available". But #164 had existed since 04:46:26 with 027e8f9 at its
head, and that head's CI run (37177966246) went red at 04:59:17 on GCC 13.
The red CI then sat for six hours with no one assigned to it. The fix,
4cd28dc, was made at 06:02 for 23b and reached this branch only through
9c791d6. The reviewer caught the false brief by checking the PR. This is
lesson 5 below.

**F3. The pushed and merged heads differ from the approved one, and no one
looked at them.** Round 3 asked for no further round "unless the pushed
head differs". 42c5f53 adds only the round-3 record. 03fc584 is GitHub's
merge of master into the branch: its committer is "GitHub", and it was
made by `gh pr update-branch` in the background job of `965f70c2` line 1251.
It brings master's #170 and #171, which are rule files and
retrospectives. Neither head got a review. The risk was low: both changes
are docs only, CI covered 03fc584, and round 3 had named remote master as
"merges cleanly". No proposal: once the merge queue of h10 is in use, the
merge of master is no longer a separate step on the branch.

**F4. Mac sleep, not CI, set the time from push to merge.** 42c5f53 was
pushed at 11:36:56. CI was green at 11:58:42 (`gh run list`). The local
polling job noticed only at 13:44, when the lid opened (`pmset -g log`:
"Wake … lid" at 15:44:16 local). After `update-branch`, CI ran from 13:45
to 14:11. The merge waited for the next lid wake, at 14:19:54, and landed at
14:23:36. In total: 2 h 47 min, of which 48 min was CI. Ola had switched on
"Allow auto-merge" at 11:41 (`965f70c2` line 1286). GitHub's auto-merge
"merges a pull request automatically after all required reviews and status
checks pass" (docs.github.com, *Automatically merging a pull request*, read
2026-10-04), with no laptop needed. h10 covers this; this is the evidence
for lesson 12 below.

**F5. Commits without a persona tag.** All the night's step commits on #164
end with `(@persona)`. The ones after h9's brief template went live
(acec5d0, 10:06) do not: 3b9408f (`@tester`), b2f6176 (`@architect`),
42c5f53 (`@architect`, recording round 3). Elsewhere today, so do the red
and green commits of increment 29's PR 1 (8d87f77, e43020b) and of h10
(fad1843, 391578c). These are the commits the red-before-green record rests
on. `.claude/briefs/common.md` asks for the Co-Authored-By trailer, not the
tag. No rule file requires the tag, so this is a gap, not a breach.
Proposal P8.

**F6. The main session resolved docs conflicts itself, then handed the
next ones to `@developer`.** At 04:45 the main session merged 23b into
23c-1 and fixed the ROADMAP conflict with its own script (027e8f9;
af318b82 line 11484). At 10:43 it gave the next merge, 9c791d6, to
`@developer`, and allowed it to write the conflicting docs files as
"mechanical" (`38de7f93` line 417). Since 2026-10-04's ruling 14(a),
`ROADMAP.md` belongs to `@architect`. Neither route is clearly wrong, but
the hard-limits list in `next.md` needs to say whose job a merge conflict
in a docs file is. Added there as evidence.

**F7. A persona edited a doc owned by another.** `@developer`'s green
commit 3ffae13 moved a citation in `docs/increments/15f-edge-strip.md`,
because its own change had moved the cited binding line, and it said so in
the commit message. Same pattern as F6, same place in `next.md`.

### The main session's deviations

- **D1 (rule broken). Personas were spawned after a persona file changed,
  with no restart.** #171 changed `orchestrator.md`, `reviewer.md` and
  `tester.md` and merged at 11:15. The session forked at 11:47:39
  (`SessionStart:fork`, `38de7f93` line 1046). After that it spawned
  `@reviewer` (14:48), `@tester` (14:46) and this `@orchestrator` (14:46),
  all without a restart. The rule:
  "After merging a change to `CLAUDE.md` or `.claude/agents/`, restart
  before spawning the changed persona" (`.claude/REQUIRED-READING.md@e56a2bf:47`).
  The proof that the fork did not reload the files is this run's own
  system prompt. It still has §1's "a review that went past two rounds"
  line, which 06a513b removed. It lacks §3's `rule_sizes.py` sentence,
  which f940bfd added. `tools/brief.py`'s "read the persona file from
  disk; the files win" was the backstop, and it worked: the personas
  noticed the difference (lesson 9). Whether a fork counts as a restart is
  unchecked; this case says it does not.
- **D2 (rule broken). A brief asked `@developer` to edit a test file**
  (lesson 10). The PR-1 green brief for increment 29 at 14:37:51
  (`38de7f93` line 1811) said: "Remove the REMOVE-AT-GREEN guard (that is a
  tests/ file outside your write limit: if the guard refuses, leave it and
  say so)". The brief knew the step was outside the role, and left it to a
  guard to refuse. `@developer` declined by its own rule
  (`.claude/agents/developer.md@e56a2bf:20`). `@tester` then made the
  removal as its own commit, edef966.
- **D3. A false statement of fact in a brief** (F2).
- **D4 (rule broken). The merge did not update the ROADMAP row** (F1). This
  is the main session's step, because it performs the merge.

## 2. Lessons reported on 2026-10-04

Each lesson gives its source, what I checked, and whether a rule already
covers it.

**L1. A chain must be "taut", not just 8-connected** (`@reviewer`,
increment 29 design round 3; `@architect`: simulate such claims before
writing them). The round-3 entry in
`docs/increments/29-nve-reference-catchments.md` (worktree
`nve-catchments`) explains it. The flood names a node's downstream
neighbour when it first queues the node. So if any two non-consecutive
chain nodes are 8-neighbours, the earlier node can drain into the later
one and skip the nodes between. The entry simulates it with the flood of
`upstream.hpp`: the path (2,4), (3,4), (3,5) gives `flow_to` (2,4)→(3,5).
In round 4, the rule was checked on 2,000 random chains: 1,738 raw chains
failed, and no taut one did. Covered? The general rule is
`REQUIRED-READING.md`'s "run the check before you write it down", and
`architect.md` checks the literature, not behaviour. The specific lesson is
domain knowledge. It belongs in the `computational-geometry` skill as a
trap, not in a rule (P6).

**L2. Rule out a simulated probe's own artefacts before its failures
count** (`@reviewer`, increment 29 round 4). Round 4 says: "the remaining
failures all start at a lower node beside the chain, which `drains`
reports". Those were failures of the probe's setup, not of the rule. It
pairs with the existing "make the probe able to fail": a probe must also be
able to pass for the right reason. Proposal P7 is one sentence.

**L3. A ruling Ola gives in chat goes into the design file as a ruling,
with his words** (`@reviewer`, increment 29 round 4). Round 4 records the
figure-script ruling with Ola's words ("map questions: 1: No, 2: Yes.").
Its wording lesson was about where the file gave the architect's argument
for the same conclusion in place of the ruling. Covered in part: the brief
template now carries "Ola, verbatim" (`.claude/briefs/common.md`). Nothing
says the design file records a ruling as quoted. Same family as the
role-bleed item of 2026-09-28 in `next.md` (`@architect` recording text
corrections as Ola's rulings). P5.

**L4. A status line copied into the ROADMAP row drifts** (`@reviewer`).
The cases: increment 29 round 2's blocker 7, round 5's "`ROADMAP.md:54`
still says 'awaiting round 4'", 23b round 6 item 3 and round 8, and now F1
on master. Ola has ruled 14(a): the row says only designed, in progress or
shipped, with the PR number. He ruled 14(b), a generated column, for stage
5. The ruled row format would have stopped most of these, but master's row
23 is still in the old style. P1.

**L5. Briefs should carry the PR's remote head** (`@reviewer`, 23c-1 round
3). See F2. `tools/brief.py` builds the block already, so the remote head
can be read there (`gh pr view <pr> --json headRefOid,statusCheckRollup`)
instead of written by hand. P2.

**L6. Pin line citations into rule files as `path@sha:n`** (`@architect`,
h10). A bare citation into a rule file made `check_citations.py` exit 1, and
the retrospectives already use the pinned form for this reason. The rule
files change often, so a bare number is false within a day. P3.

**L7. The push guard's list missed `gh pr update-branch`, and its text
fallback is untested** (`@architect`, h10; `@reviewer`, h10 round 2). Ola
said yes at 11:20:08 ("yes, add update-branch to the guard"). The red step
is fad1843 and the green step 391578c. h10 round 2's mutant check: removing
`update-branch` from the parsed verb tuple
(`.claude/hooks/guard_push.py@35673c8:151`) is caught; removing it from the
text pattern (`:40`) survives, "because no row reaches the text fallback".
The verb list now has three copies: two in the guard and one in prose
(`.claude/REQUIRED-READING.md@e56a2bf:136-142`, which still lacks
`update-branch`). P4 and cut C2.

**L8. The main checkout's `.venv` has an old `_core.so`, and pybind11 reads
`PYTHON_EXECUTABLE`** (`@developer`, h10 and the 23c-1 merge). Checked: the
file `.venv/lib/python3.14/site-packages/tin_engine/_core.cpython-314-darwin.so`
is dated 29 September, while `bindings/core.cpp` last changed today. The
pybind11 docs (*Build systems*, read 2026-10-04) say "an exact Python
installation can be specified with `PYTHON_EXECUTABLE`" in the classic mode,
which is the default with `find_package(pybind11 CONFIG)` (`CMakeLists.txt`
line 87). They also say pybind11 switches to FindPython where CMake 3.27 has
removed the old mechanism; I did not check which mode this machine uses.
Covered: `REQUIRED-READING.md`'s "Stale artifacts" section has the rebuild
commands. Not covered: a brief that asks for "the full pytest" in a run that
may not build C++. P2 adds that to the brief block.

**L9. The persona files on disk differ from the copies in the agents'
system prompts** (`@reviewer`, `@tester`). Confirmed and explained under
D1.

**L10. Removing a test-file scaffold goes to `@tester` as its own commit**
(`@developer`, increment 29 PR 1). Covered by
`docs/increments/README.md` step 3 and `developer.md` line 20. The main
session broke it (D2). No new rule. The brief template could say it: P2.

**L11. Mac sleep stalled `@tester` twice, and committing the draft first
saved the work** (`@tester`). `pmset -g log` shows two lid-closed sleeps on
battery: from 15:30:05 to 15:44:16 local, and from 16:15:29 to 16:19:54
local. 0624bbf ("red WIP, draft accumulate suite from the stalled attempt",
16:21 local) is the draft that survived. Covered:
`REQUIRED-READING.md`, "A subagent whose product is a file creates that
file first and writes incrementally". This is that rule working.

**L12. Merging three PRs one after another under "require up to date" re-ran
CI in full each time** (main session, with Ola). That led to h10, the merge
queue. F4 gives #164's timing: of 2 h 47 min, 48 min was CI and the rest was
Mac sleep. In the queue, GitHub does the merge itself. Covered by h10.

## 3. Proposals for Ola

Each is one change, with its evidence above. All but P6 touch governed
files, so they go through the pipeline by day.

- **P1 (stale rows; L4, F1, D4).** First, fix master's row 23 and the 23
  status line now (`@architect`, not governed; a few minutes). Then a check:
  for each increment whose PR `gh` reports as merged, its ROADMAP row and
  `Status:` line must not say "in review". It goes in the merge step or the
  recap, and replaces 14(b) until stage 5. About 30 lines, with tests, in
  `tools/session_state.py`. The source of the pattern is towncrier, already
  cited in 14(b). The other choice is to wait for stage 5 and accept stale
  rows until then.
- **P2 (briefs; L5, L8, L10, D2).** Two lines in
  `.claude/briefs/common.md`, filled in by `tools/brief.py`. First, the PR's
  remote head and the state of its checks, when the branch has a PR.
  Second, the C++-build allowance stated next to any test command, as
  P4 of 2026-10-04 already proposes. For a green brief, also the line
  "scaffold removal is `@tester`'s, as its own commit". Cost: about 15 lines
  and tests in `tools/brief.py`.
- **P3 (citations; L6).** `reviewer.md` §5 and `architect.md`: a line
  citation into `CLAUDE.md`, `.claude/**` or `docs/increments/README.md` is
  written as `path@sha:n`. One line each. A tool alternative:
  `check_citations.py` warns on a bare citation into a rule file, about 10
  lines.
- **P4 (guard lists; L7).** In `guard_push.py`, build the text pattern
  from the same verb tuple as the parsed check, so there is one list. Add
  test rows that force the text fallback (a line `shell_scan` cannot read),
  one for each verb. Then add a test that runs `gh pr --help`, lists the
  subcommands, and fails on any writing subcommand that is neither guarded
  nor named as exempt (`close` and `comment` are exempt today). This finds
  the next `update-branch`. About 25 lines and tests (`@tester`, then
  `@developer`).
- **P5 (rulings; L3).** One line in `architect.md`: "A ruling Ola gives in
  chat goes into the design file under its ruling heading, with his words
  quoted and the time; your reasoning, if any, follows it, labelled as
  yours."
- **P6 (taut chains; L1).** A trap entry in the `computational-geometry`
  skill. It is not a rule file. "A priority flood sets a node's receiver
  when it is first queued, so a chain meant to drain node by node must be
  taut: no two non-consecutive nodes equal or 8-neighbours." About 40
  words.
- **P7 (probes; L2).** One sentence added to `REQUIRED-READING.md`,
  "Claims", after "Make the probe able to fail": "and rule out the probe's
  own setup before you count its failures".
- **P8 (persona tags; F5).** One line in `.claude/briefs/common.md`:
  "End each commit subject with `(@$persona)`". The red/green record and
  this check depend on it.

Not proposed: a rule for a pushed head that differs from the approved one
only by a review record or GitHub's merge of master (F3). The merge queue
removes the second case.

Open, for Ola to rule: who runs the mutants of a suite named
invariant-critical? `reviewer.md` §5's fourth check says the mutants "had"
been run, "with the kill record in a handback". `reviewer.md` also says
"You do not edit or commit". But the review brief for increment 29's PR 1 at
14:48:12 (`38de7f93` line 2015) asked `@reviewer` for "the mutation round".
Running mutants means editing production files for a while. The review
entries of 23a-1 and 23a-2 also report mutation passes; I did not check
who ran them. Default: `@tester` runs the mutants as its own
task, and `@reviewer` checks the kill record.

## 4. Rule text: size and a cut

`python3 tools/rule_sizes.py` at e56a2bf: **10,156 words**, up 266 since
the recovery-round retrospective (dd48544). Biggest changes:
`.claude/briefs/common.md` is new (150), `tester.md` +70,
`REQUIRED-READING.md` +44, `reviewer.md` +35. `CLAUDE.md` is down 46.

The recovery round proposed cut C1 (about 90 words of
`REQUIRED-READING.md`); it has not been made, and the file grew instead.

**C2.** In `REQUIRED-READING.md`'s harness section
(`.claude/REQUIRED-READING.md@e56a2bf:136-142`), replace the list of what
`guard_push.py` asks before (about 70 words) with one sentence: "asks before
anything that writes to the remote, rewrites history or changes git
configuration; `tests/python/test_guard_push.py` lists each case." About
50 words saved. It also removes the third copy of a list that is already
out of date (L7). The guard's docstring does not hold a full list today
(read at e56a2bf), so the tests are the place to look.
