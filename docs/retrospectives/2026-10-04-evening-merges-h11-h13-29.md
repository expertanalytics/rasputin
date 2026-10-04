# 2026-10-04 evening: the check after h11 (#176), h13 (#177) and increment 29 (#173, #178)

Recorded by `@orchestrator` on the night of 2026-10-04, on master at
529613a, during Ola's unattended night. Part 1 is the post-merge check that
`CLAUDE.md` §3 asks for. Part 2 records the ten lessons the main session
passed on. Part 3 lists proposals for Ola, Part 4 measures the rule text.

What the labels mean: **h11** is the harness increment that skips the code
jobs in CI on a pull request (PR) that changes only prose, and adds one
required check, `CI result`; **h13** limits the number of parallel build
jobs on CI, which made the macOS build slow; **increment 29** builds
catchments for NVE's reference gauging stations (NVE is the Norwegian Water
Resources and Energy Directorate), and its **PR 1** (#173) is flow
accumulation in C++, its **PR 3** (#178) the station and river readers.

All times are UTC. Git prints commit times in local time (UTC+2); they are
converted here. Merge times are GitHub's `mergedAt` (`gh pr view <n> --json
mergedAt`), the moment the merge queue landed the PR; the queue's merge
commit carries an earlier committer time, when the queue built it (#173's
9fb1fe8 says 17:32, #175's d926644 says 18:30). Transcripts are under
`~/.claude/projects/-Users-skavhaug-projects-rasputin/`: `86af816a…` is the
main session from 15:08 to 18:56, `806b4380…` the one after the restart at
18:56 that spawned this run. "Line n" is a line of that JSONL file.
Subagent timings come from
`806b4380…/subagents/agent-*.jsonl` (first and last timestamp).

## 1. The merges against the rules

### 1.1 #173, increment 29 PR 1 (merged 18:02)

The design rounds and the red and green commits were checked in
`2026-10-04-23c-1-and-day-lessons.md`. What came after: code review rounds
1 to 3 (15:19 to 15:40, `86af816a` lines 276-718), each recorded by
`@architect` on the main session's brief (36c2370, 93474dc), the push on
Ola's "yes, push and open PRs" (17:04:20) and the enqueue on "merge #173 and
#174 when CI is green" (17:17:53; `gh pr merge` at 17:30:35, line 1140).
Step order holds. Four `@architect` commits (81392f5, 0b86517, 725c8cd,
5dbfec0) carry no `(@persona)` tag; the tag rule merged later, at 18:54
(#175), so that is not a deviation. #173 was the first PR through the merge
queue, and the queue's commit is the one that landed (h11 §5 step 0, checked
by `@reviewer` in h11 round 1).

### 1.2 #176, h11 (merged 20:22)

Order: design 2545c07 (17:27), Ola's three rulings at 17:38:45 recorded by
`@architect` (0e116a8), red 21d675e and 8db5bfc (`@tester`), a design
amendment on who owns the test hook and on the install trap (9e5203e), a red
amendment 8986083, green 8b5c4ee (`@developer`, 18:11), rules text b2e0a7b
(`@architect`), review rounds 1 and 2, push at 19:28:58 on "push h11", PR at
19:38:15 on "Yes, open PR". Holds.

- **The enqueue raced its own question.** At 19:38:29 the main session asked
  "Shall I put it in the merge queue?" and ran `gh pr merge 176` at 19:39:26
  with no answer in chat (line 818). `guard_push.py` asked
  (`permissionDecision: ask`), and the result came back 27 s later, so Ola
  approved the prompt. The permission prompt is a valid yes; the chat
  question was left unanswered and became noise. Not a deviation; worth
  one line in Part 3 (proposal P5).
- **Root file written by the main session.** 2faa3c4 (`testing.md`, 19:01,
  line 294, a `python3` heredoc) carries no tag and no persona wrote it.
  `testing.md` is in no persona's write limit (lesson 1, below).
- **Status line stale on master.** h11's status line says "the branch has
  no CI run yet"; the PR ran and merged. Nobody updates a design file's
  status after the merge.

### 1.3 #177, h13 (merged 20:54)

Order: diagnosis and design 7c4b178 (`@architect`, 20:07), Ola's rulings
("1: all three. 2: yes.", 20:10:18) recorded by the main session (f12a0e6,
line 1495), workflow edit a77879f (`@developer`), review round 1 (prose
only), fixes 14d6adc (`@architect`), round 2 approved, master merged in by
the main session (97bde1d), push, PR and enqueue on Ola's "yes" (20:26:51)
to the four acts named in one question (line 1855). No red step, accepted
by `@reviewer` on the h10 precedent (the effect is visible only on CI).

The result, checked with `gh run view` on the macOS job's "Build" step:
140 s on the PR run (97bde1d) and 84 s in the queue (accd52a). #178 had
three runs: the first PR run (3dff7bc, run 37232282246, started before h13
merged) took 1,189 s; the second PR run (84c0a4d) 143 s and the queue run
(529613a) 141 s. h13 §2a's seven earlier runs took 19.7 to 25.1 min
(checked by `@reviewer` in round 1).

- **Step owed.** h13 §6 says the CI numbers are recorded as a further
  review round after the PR's run. Not done; the main session told Ola it
  would "fold that into tomorrow's bookkeeping" (line 2480). The status
  line on master still says "awaiting the push".

### 1.4 #178, increment 29 PR 3 (merged 21:25)

Red 213a2ef (18:15), Ola's rulings recorded by `@architect` (b519bd9), red
amendments, green aae91bb (`@developer`), then six code review rounds.
What each round blocked on, and who wrote the text it blocked on:

| Round | Reviewer time | Blockers, and whose text |
|---|---|---|
| 1 | 9.8 min | NVE's lake number 0 read as a lake (design and code: a real defect); status line and `ROADMAP.md:54` (`@architect`); a citation the green step moved (`@developer`'s change); a file-table row (`@architect`); `project_structure.md` entries (nobody's file) |
| 2 | 7.8 min | a test citation moved by the main session's ae78666; status line (`@architect`); red-step paragraphs now false (`@tester`) |
| 3 | 1.6 min | status line says "two citations", round 2 asked for one (`@architect`, c7d427a) |
| 4 | 0.7 min | approved |
| 5 | 6.0 min | (after PR CI failed and master came in) the ruling record says only the 3.14 leg failed, all three did; "Data use" claims a test the ruling removed (main session, 9607b98) |
| 6 | 1.5 min | approved |

Lesson 6 says four of the six rounds were on prose the main session or
`@architect` wrote. The record says: rounds 2, 3 and 5 blocked on prose
alone; rounds 4 and 6 were one-line re-checks; round 1 had one real defect
among five blockers. Of the prose errors, `@architect` wrote the round-3
one and the main session the round-5 one. Measured cost of the prose-only
tail: reviewer rounds 2 to 6 took 17.6 min of agent time and about 3.5
million cache-read tokens, against round 1's 9.8 min and 5.3 million; the
fix agents for round 2 add 2.7 min. The two rounds that exist only because
a recorded sentence was wrong (4 and 6) took 0.7 and 1.5 min. The rounds
are cheap because each is scoped to the delta; the waste is in where the
errors come from, not in the review.

**The status line is the most frequent single blocker.** Across all 14
review records in `29-nve-reference-catchments.md`, the status line or
`ROADMAP.md:54` is a blocker in five (design round 2, PR 1 round 1, PR 3
rounds 1, 2 and 3) and a "required" or suggestion in two more. The status
line is now one sentence of about 170 words that narrates every round, and
it is stale on master ("awaiting the push and CI"), as is row 29.

- **The main session made a design call.** When the h11 hook failed #178's
  CI, the main session's brief to `@tester` (line 2157, 20:34:47) said "the
  design requires that credit, so the read stays" and ordered the route
  that adds `NOTICE.md` to the files the suites may read. That is
  `@architect`'s question. Ola questioned it ("Having a code test checking
  the NVE credit seems quite odd", 20:40:38) and ruled the test out ("yes,
  remove the NOTICE.md check. So this was good! We found a logical flaw!",
  20:42:00). The same brief is where "failed on Python 3.14" started; it went
  into `@tester`'s 4e877ae message and then into the main session's ruling
  record, which round 5 caught.
- **The main session reverted a test commit.** 09a44c6 reverts `@tester`'s
  4e877ae in `tests/python/test_ci_changes.py` (line 2245). Small and on
  Ola's ruling, but a test file changed by the main session; a one-line
  `@tester` brief would have done it, and `@tester` was spawned anyway
  minutes later for the deletion.
- **Re-enqueued without a fresh merge yes.** The enqueue at 20:29:23 (on
  Ola's "yes") left GitHub's auto-merge on. After the first PR run failed,
  the branch changed (test removed, master merged, two records) and the main
  session pushed 84c0a4d on "push #178 when landed". Auto-merge was still on,
  so GitHub enqueued and merged the new head with no new `gh pr merge`.
  The main session said so in advance ("Your earlier enqueue is still in
  place, so it goes into the merge queue by itself", 20:52:53) and Ola did
  not object, so this is disclosed, not hidden. But
  `.claude/REQUIRED-READING.md`'s test ("have I already performed this named
  act once, and has the tree changed since?") does not see it, because no
  act was performed: GitHub keeps auto-merge across pushes by anyone with
  write access ("Auto-merge is disabled if someone without write permissions
  pushes new changes to the head branch …", GitHub Docs, *Automatically merging
  a pull request*, read 2026-10-04). Proposal P5.

### 1.5 The main session's own work this evening

Recorded from `806b4380` (lines given), beyond what is above:

- **Prose written by the main session.** `project_structure.md` (ae78666,
  line 467: "Text drafted by @architect … applied by the main session"),
  `testing.md` (2faa3c4), a citation in `15f-edge-strip.md` (line 777),
  Ola's rulings in h13 (f12a0e6) and in 29 (9607b98), review records and
  status lines for h11, h13 and 29 PR 3, and two prose fixes "in the commit
  that records this round" (29 PR 3 rounds 3 and 5). `reviewer.md` says the
  verdict's "spawner records it", so recording is the main session's; the
  fixes and the rulings are not. In `86af816a` the same records were
  `@architect`'s (spawns "Record 29 PR 1 review round 3" and others); after
  the restart the main session did them itself. Every one of these edits
  went through `python3` heredocs or `sed -i`, none through Edit or Write.
  None was a refused write, and none of the files is governed, so no guard
  was worked around; but the route keeps these edits out of any per-tool
  check a later hook might add.
- **Two briefs reconstructed by hand, both refused by `guard_spawn.py`.**
  At 19:47:36 (line 947) the brief's "Review rounds recorded" line ended
  "Verdict: CHANGES REQUESTED"; `tools/brief.py` cuts that line at 200
  characters and had printed "Verdict: CHANGES". At 20:00:05 (line 1357)
  the main session had printed only the first line and the tail of the
  block (`sed -n '1p;/Ola, verbatim/,$p'`, line 1342) and rebuilt the middle,
  with a different note-file name. Both refusals said "brief: edited; paste
  the output of tools/brief.py unchanged", and both retries pasted a fresh
  block unchanged. The guard worked as designed; the retry is the route it
  asks for, not a workaround. The first refusal comes partly from the tool:
  a cut that drops "REQUESTED" invites completion (proposal P4).
- **`session.md` misuse, small.** One `ASK OLA:` line ("queued after:
  @orchestrator (29 PR1, rules, h11); h12; cite by phrase; …") is a queue,
  not a decision waiting on Ola.
- **No idle and no concurrency breach seen.** Writers never exceeded two
  (checked spawn by spawn against handbacks). The unattended window began at
  about 21:27; `@tester` (29 PR 2 red) and `@architect` (h14 design, done at
  21:38:54) started at 21:28. The night itself is for the morning check.

### 1.6 Not covered here

#172 (h10, the merge queue), #174 (ROADMAP row 23) and #175 (the rules
branch, which carried the previous retrospective) merged since
`2026-10-04-23c-1-and-day-lessons.md` and were not in this brief. Their
records were not re-read.

## 2. Lessons reported this evening

Each lesson as the main session passed it on, then what the record shows
and where it leads.

**L1. Root Markdown files are outside every write limit.** `tools/brief.py`'s
`WRITES` gives `@architect` `docs/` (less retrospectives), `ROADMAP.md`,
`CLAUDE.md` and `.claude/**/*.md`, and nobody the other root files:
`README.md`, `INSTALL.md`, `NOTICE.md`, `testing.md`, `project_structure.md`,
`auto_catchments.md`, `parallel_refinement.md` (`ls *.md` at 529613a). So the
main session wrote `project_structure.md` (ae78666) and `testing.md`
(2faa3c4) itself. Before h9 made the limits explicit (#169, merged 10:06,
merge acec5d0), `@developer` (23d4dad) and `@perf` (acdfce0, c6e0049) wrote
root files in their own commits. After h9, one write went outside a limit:
aae91bb, 29 PR 3's green step (`@developer`, 18:39), added 13 lines to
`NOTICE.md` (NVE's credit). The main session briefed it (`86af816a` line
2376, 18:28:40) with the write limit "src_python/, include/, src/,
bindings/, tools/, .claude/hooks/, .github/, CMakeLists.txt,
pyproject.toml", and the brief does not name `NOTICE.md`; but the design
asked for the credit (the "What is committed" ruling, and the file table's
`NOTICE.md` row marked "docs") and `@tester`'s red 213a2ef tested it
(`test_notice_md_credits_nve_under_nlod`). So the green step could pass only
by writing outside its limit, and no brief, guard or review round flagged
it. The gap starts in the table, and after h9 it showed in conduct: a
design and a red test that require a file no persona owns, and a write
outside the limit that nobody noticed. Proposal P1.

**L2. A design rule over an external field needs a count of its sentinel
as well as null.** NVE's river service sends `vatnlnr` (lake number) = 0 for
"no lake" on 956,447 features; the design's "set" meant "not null", so a
blank-type river with 0 became a lake (29 PR 3 round 1). The design's own
measurement ("set on 496 river features") counted the zeros and so looked
like evidence. Proposal P7 (a trap line in the `geospatial-data-formats`
skill, not a rule file).

**L3. Rebuild `_core` after merging master, before a full suite.**
Incident: 851a497 merged master (with #173's `accumulate`) into the h11
branch; `@reviewer`'s h11 round 1 (handback, `806b4380` line 507, 19:16)
found the worktree's `.venv` extension stale, so `test_core_accumulate.py`
failed there, and rebuilt `_core` in a scratch build to measure.
`.claude/REQUIRED-READING.md`, "Stale artifacts", already says to rebuild
"before every `pytest` meant to measure C++". After a master merge every
full run measures C++, because master's C++ changed under the venv's
extension; the reader of a Python-only branch does not think of it as
"meant to measure C++". The main session's briefs now say it in words
(line 2157: "master brought C++ changes since the venv's build: rebuild per
REQUIRED-READING"). Proposal P8 puts it in `tools/brief.py`.

**L4. Briefs should say where "green" ran.** Folded into P8, and into P2
of `2026-10-04-23c-1-and-day-lessons.md` (the PR's remote head and check
state in the brief), still waiting on Ola.

**L5. A status line or a brief summary is a claim to check.** Two cases
tonight: "two citations" for one (`@architect`, c7d427a, round 3), and
"failed on Python 3.14" for all three legs (main session's brief, line
2157, copied into 4e877ae and 9607b98, round 5). The second travelled
through three hands before review caught it. The rule exists ("Claims: run
the check before you write it down"); what is missing is a habit for the
dispatcher's own summaries, and a status line that does not restate counts
at all (P3).

**L6. Six review rounds on #178.** Measured in §1.4: the prose-only tail
cost 17.6 min of reviewer time and about 3.5 million cache-read tokens; the
two rounds caused only by a wrong recorded sentence took 0.7 and 1.5 min.
No review-round cap is written down, and Ola declined to add one this
morning (`2026-10-04-recovery-round.md`, Q2); the numbers do not argue for
one.
The cut is at the source: status lines that narrate rounds (five blockers
in 14 rounds), and the main session writing prose it then has reviewed.
Proposals P2 and P3.

**L7. h11's prose-read hook caught a design flaw on its first PR.** The
hook (`tests/python/conftest.py`) failed #178 because
`test_notice_md_credits_nve_under_nlod` read `NOTICE.md`: a unit test that
checks the wording of a prose file, which a prose-only PR could then change
without the suite running. Ola ruled the test out (20:42:00). The guard did
its job; the dispatcher's first response (exempt the file) did not, and the
hook's message steers that way: "Add it to NOT_PROSE in tools/ci_changes.py,
or stop reading it." Proposal P6.

**L8. `pytest -q | tail` can hide the hook's report.** The hook turns the
session red at session end; a pipe into `tail` reports `tail`'s exit
status, and a short tail can cut the hook's lines. Check `pytest`'s own
exit status (`set -o pipefail`, or `${PIPESTATUS[0]}`). No incident: it
is a precaution `@tester` reported in its ea489e5 handback (`806b4380` line
2392, 20:51); no run in the record had the report hidden. Into P8.

**L9. The git-archive scratch copy fails 7 or 8 `test_settings_wiring`
tests.** Reproduced at 529613a: in a `git archive HEAD | tar -x` copy,
`tests/python/test_settings_wiring.py` gives 7 failed, 20 passed, all in
`test_hook_is_executable_in_the_checkout`, which asks git's index for each
hook's mode (`git ls-files -s`); the copy is not a git work tree. `tester.md`
(merged in #175 today) requires exactly that kind of copy for mutants, and
`@reviewer` ran a full suite in one (29 PR 3 round 2, recorded in c7d427a:
"8 failed, all eight `test_settings_wiring`'s executable-hook check"). A suite that is red by
construction in the sanctioned copy makes every full run there need a
footnote. Proposal P9.

**L10. Three CI builds per merge, now two.** With the merge queue (h10) and
the old `push: master` trigger, #175 built three times: PR run 18:07 to
18:29, queue run 18:30 to 18:53, push run 18:54 to 19:20 (`gh run list`).
h11 removed the push trigger; no push run follows #176, #177 or #178. Ola's
words at 19:52:54: "No the thing is that it does not merge after the merge
queue. It still build one more effing time." The duplicate that remains,
PR run plus queue run (#178: 14 min plus 12 min), is h12's question.

**Waiting, as Ola saw it.** From PR open to merge: #176 44 min (PR run 20.7,
queue 22.3), #177 26 min, #178 56 min including the failed first run. Ola's
turns at 19:57:09 ("#176 still not done? How long now?"), 20:06:06 ("I won't
do anyting more until #176 is merged and active") and 21:25:31 ("Is #178
about done soon?") are the cost in his attention. h13 took the slowest job
from about 22 min to about 2; the sanitizer jobs, at about 10 min, now set
the floor (h14 is being designed tonight).

## 3. Proposals for Ola

Each is one change, with its evidence above, its source and its cost. All
but P7 touch governed or rule-stating files, so they go through the
pipeline by day.

- **P1 (root files have an owner; L1).** Every tracked path has exactly one
  writer. Add "root `*.md` files other than `CLAUDE.md`" to `@architect`'s
  entry in `tools/brief.py`'s `WRITES`, and a test that every file
  `git ls-files` lists matches at least one persona's limit or a short
  named list of Ola-only paths (`.claude/settings*.json`, `LICENSE`). The
  test finds the next gap before a main session fills it. Source: GitHub's
  CODEOWNERS, whose catch-all line gives "the default owners for everything
  in the repo", and where "the last matching pattern takes the most
  precedence" (GitHub Docs, *About code owners*, read 2026-10-04). Cost:
  one entry and about 20 lines of test (`tools/brief.py` is governed).
  Question for Ola below on `testing.md`.
- **P2 (the dispatcher records verbatim, and only that; L5, L6, §1.5).**
  `reviewer.md` says the verdict's "spawner records it". Narrow it: the
  spawner pastes the reviewer's verdict paragraph unchanged under
  `## Review`; any other edit to the increment file, the status line
  included, goes to `@architect`; a ruling of Ola's goes in quoted, with the
  time, by `@architect` (this repeats P5 of the 23c-1 retrospective, which
  waits on Ola). The two wrong sentences that cost rounds 4 and 6 were the
  dispatcher's and `@architect`'s own summaries, not the reviewer's words.
  Stage 5 of the dispatcher-control plan ("reviews recorded",
  `docs/research/2026-10-03-dispatcher-control.md`) would make it a tool;
  this is the one-line version now. Source: Anthropic, *Building effective
  agents* ("Change the arguments so that it is harder to make mistakes",
  read 2026-10-04): a paste cannot misquote. Cost: one sentence in
  `reviewer.md`.
- **P3 (a status line is a state, not a log; L5, L6, §1.2-1.4).** The
  status line says one of design, red, green, in review, approved, merged
  (with the PR number), and points at `## Review` for the rest; no counts,
  no round narrative. After the merge, the merge step sets it to "merged as
  #N" and the ROADMAP row the same, which also covers P1 of the 23c-1
  retrospective. Evidence: five status-line blockers in 14 review rounds of
  increment 29; three merged increments (h11, h13, 29) whose status lines on
  master are stale tonight; increment 29's line is about 170 words. Cost:
  one sentence in `docs/increments/README.md`; a check in
  `tools/session_state.py` (merged PR, status not "merged") of about 30
  lines. A first fix of tonight's three stale lines is bookkeeping and
  needs no ruling.
- **P4 (`brief.py` prints the verdict, not 200 characters; §1.5).** The
  "Review rounds recorded" line prints the first 200 characters of the last
  review paragraph, which on long records ends mid-verdict ("Verdict:
  CHANGES"). Print the round heading and the verdict word instead. The
  hand-completed "REQUESTED" was one of tonight's two `guard_spawn.py`
  refusals (`806b4380` line 947; a refused spawn leaves no commit, so the
  transcript line is the only record). Cost: about 5 lines and a test (`tools/brief.py` is governed).
- **P5 (a push to a PR with auto-merge on is also a merge; §1.4).** One
  sentence in `.claude/REQUIRED-READING.md`, approval boundary: "A push to a
  PR whose auto-merge is on also enqueues it: ask for both, or run
  `gh pr merge --disable-auto` first." And in `guard_push.py`'s ask reason
  on `git push`, name the PR's auto-merge state when it is on (about 15
  lines and tests). Also: ask in chat or let the permission prompt ask, not
  both at once (#176, 19:38 to 19:39). Source: GitHub Docs, *Automatically
  merging a pull request*: auto-merge is disabled only when "someone without
  write permissions pushes new changes" (read 2026-10-04).
- **P6 (the hook's message names the right first question; L7).** In
  `tests/python/conftest.py` (`@tester`'s file), reword the message: "A test
  that checks the wording of a prose file belongs in review, not in the
  suite: remove the read. Add the file to NOT_PROSE only if the build or a
  tool really reads it." Cost: two lines. The main session's first route
  on #178 followed today's wording.
- **P7 (sentinel values; L2).** A trap line in the
  `geospatial-data-formats` skill (not a rule file): "When a design says a
  service field is 'set', count null, empty and the service's sentinel
  (often 0 or -1) separately, and say which one 'set' excludes." About 35
  words. Incident: the defect landed in aae91bb and was fixed in c4e50ff
  (red b59c648).
- **P8 (the brief's test line; L3, L4, L8).** `tools/brief.py` already
  states the C++-build allowance. Add one line wherever a full suite is
  asked for: "After merging master, rebuild `_core` first; read `pytest`'s
  own exit status, not a pipe's; say whether green ran locally or on CI."
  Merges with P4 of the recovery-round retrospective and P2 of the 23c-1
  one, both still waiting. Cost: about 10 lines and a test. Incident for
  L3: 851a497, the master merge that left the h11 worktree's `_core` stale
  (`806b4380` line 507). L8 has no incident; it is a precaution.
- **P9 (the sanctioned copy runs green; L9).** In
  `test_hook_is_executable_in_the_checkout`, skip with a stated reason when
  the tree is not a git work tree, as the index check cannot run there; the
  on-disk executable check could stay, but `git archive` keeps modes, so
  the remaining failures are the index lookups. Owner `@tester`. Cost:
  about 4 lines. The other choice is for `tester.md` to name the suites
  that do not run in a copy; that adds rule text, so it is not the default.
  Incident: c7d427a records round 2's "8 failed" run in the copy that
  `tester.md`'s copy rule (merged in d926644, #175) requires.

Not proposed: a rule against the main session editing with `python3`
heredocs or `sed -i`. The files were ordinary and no guard was refused;
P2 removes most of the occasions.

### Questions for Ola, each with a default

1. **Who owns `testing.md`?** It states the coverage floor that `@tester`
   enforces, but it is prose about the suites. Default: `@architect`, like
   the other root Markdown files (P1), with `@tester` named in the review.
2. **May the main session apply text a persona drafted, when the file is in
   no persona's limit?** It did tonight (ae78666). Default: no; P1 closes the
   gap instead.
3. **Do P2 and P3 come before stage 5, or wait for it?** Default: now, as
   the one-sentence versions; stage 5 replaces them.

## 4. Rule text: size and a cut

`python3 tools/rule_sizes.py` at 529613a: **10,506 words**, up 350 since
the 23c-1 retrospective (4e8a740). The growth: `CLAUDE.md` +155 (h10's
merge-queue paragraph and h11's `CI result` paragraph), `tester.md` +94
(mutants in a scratch copy, from #175), `reviewer.md` +47,
`.claude/briefs/common.md` +31, `docs/increments/README.md` +9,
`REQUIRED-READING.md` +7, `PRINCIPLES.md` +7. The cuts proposed before,
C1 (about 90 words of `REQUIRED-READING.md`) and C2 (about 50 words, the
`guard_push.py` list), have not been made.

**C3.** In `CLAUDE.md` §4, h11's paragraph has a sentence and a half written
for the time between the merge and Ola's settings change: "it is meant to be
the one required check on `master`, and `gh api …` shows whether it is yet.
While the per-job checks are still required, a prose-only PR waits on matrix
checks that never report." That time is over: `gh api
repos/expertanalytics/rasputin/branches/master/protection/required_status_checks`
returns `CI result` alone (app 15368, GitHub Actions; read tonight, and
confirmed by the main session at 20:25:38). Replace the 37 words with "it is
the one required check on `master`" (8 words). About 29 words saved. The
command can stay in h11 §5, where the record is.

Taken together, C1, C2 and C3 would remove about 170 words, half of what was
added since the last retrospective.

## Review

**`@reviewer`, round 1, 2026-10-04.** Range `529613a..43fbd5f`. Verdict:
CHANGES REQUESTED. About 60 claims checked, five wrong: L1 put aae91bb
before h9, but h9 merged at 10:06 (acec5d0 is its ancestor), so its
`NOTICE.md` lines are a write outside `@developer`'s limit after h9; the
merge times of #173 and #175 were the queue commits' times, not GitHub's
`mergedAt` (18:02 and 18:54); #178 had three CI runs, and the first macOS
build (3dff7bc) took 1,189 s; rounds 4 and 6 took 0.7 and 1.5 min, not
"about 2 min each"; P4, P7, P8 and P9 named no incident commit. Two
suggestions: the review-round cap is "not written down, and Ola declined to
add one", and an ellipsis where the GitHub quotation stops. All fixed in the
commit after 43fbd5f.

**Round 2, 2026-10-04 (copied from `@reviewer`'s handback).** Range `43fbd5f..9774a68`. Verdict: APPROVED. All five fixes checked by running what each asserts; both round-1 suggestions taken; `check_citations.py` exits 0. Non-blocking: line 507 of `806b4380` is 19:17:10, not 19:16.
