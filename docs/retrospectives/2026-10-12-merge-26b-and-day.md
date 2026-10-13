# 2026-10-12: increment 26b merged, the ruling record and the morning retrospective merged, increment 37 pushed and the generality audit waiting on a push as of 14:43Z; the day's lessons

`@orchestrator`, one post-merge check for #235 (increment 26b, exact
land-cover class fractions per triangle, `docs/increments/26-land-cover.md`
section 14, merge commit `c317c7e4` created 11:30:05Z, merged 11:39:51Z),
with one line each for #233 (the hotspot ruling record, `ac9ee6a4`) and #234
(the morning retrospective, `28a96ddb`), and a check of two finished items,
unmerged at this file's first commit (`f4ebb3e5`, 14:43Z; the tree moved
while it was written, see the Review section): increment 37 (the MCP
server, `worktree-mcp-server` at `33061c0f`, review-approved, pushed as
#236 at 14:38:03Z, queued 14:47:09Z) and the generality audit
(`worktree-generality-audit` at `1c405ade`, review-complete; Ola's product
review, not re-reviewed here). Sources: `git log`, the increment files, the
GitHub timeline of each PR (`gh api repos/.../issues/<n>/timeline`),
`tools/count_loc.py`, the harness queue and window files, the main session's
transcript of today (`6a989a29`), `docs/benchmarks/2026-10-12/`. Nothing was
rebuilt or re-run. **Every time below is UTC from git (`%ci`) or GitHub;
the main session's own figures for today are local time (CEST, UTC+2)
labelled UTC**, see D4.

## 1. Did each item follow the pipeline

| Item | Design, review | Red | Green | Mutation | Code review | `@perf` | LOC (`count_loc.py`, merge parents) | Deviations |
|---|---|---|---|---|---|---|---|---|
| 26b (#235) | 26's design of 2026-10-03, two rounds then; no freshness pass (section 14 says why) | `1820c1d4` 07:57Z (`@tester`, 92 cases, 68 red) | `5b8bb3e0` 08:23Z (`@developer`), test-side fixes `f165282b` 09:01Z (`@tester`) | none (section 9 names none in 26b) | 2 rounds, APPROVED `68d4b75b` 09:54Z | `98eea44f` ACCEPTED, `docs/benchmarks/2026-10-12/26b-fractions/` | **621 net** (estimate 420) | D1, D4, D5, D6 |
| ruling record (#233) | prose, `50f19dae` 04:53Z (main session) | - | - | - | 2 rounds, `7bba17f8` 07:57Z | - | 0 | D1 |
| retrospective (#234) | prose, `a749eab7` 07:38Z (`@orchestrator`) | - | - | - | 2 rounds, `2e35f70e` 08:03Z | - | 0 | D1 |
| 37 (#236, unmerged at 14:43Z) | `b990e38c` 08:34Z (`@architect`), 3 rounds, APPROVED `ba81cde1` 09:07Z | `ee5b5d52` 10:13Z (62 tests, all red) | `851592ba` 10:54Z | none (section 9) | 3 rounds, APPROVED `33061c0f` 12:15Z | not needed (no refine or mesh code) | **328 net** (`count_loc.py c317c7e4 33061c0f`; estimate 260, worst case 360) | D1, D2, D7 |
| audit (unmerged) | `6d468ba6` 12:23Z (`@architect`) | - | - | - | 1 round, `1c405ade` 12:44Z | - | 0 | D1 |

**26b in detail.** Order held: red before green, the red commit touches
only `tests/` (six files, 1,711 lines), the green commit no test file (ten
production files), the test-side fixes are `@tester`'s own commit with the
reasons in its message (`git show --stat` on each); `@perf` ran before code
review (09:37Z before 09:48Z); both review records are one paragraph with
full-path pinned citations; the as-built record is the main session's
(`4149b50d`), and names what it decided on defaults. Red/green separation
holds by the letter and the spirit: the four cases that stayed red after
green were test-side (a rounded-input expectation, a ragged `.ply` reader,
a missing directory, four layering rows), fixed by `@tester`, not folded
into the green commit.

- *LOC 621 against 420.* `overlap.hpp` 231 counted lines against the
  estimate's 150, and `cli.py` 129 against 55 (the seven refusals, each
  tested). Review round 1 rejected the first explanation (57 comment lines
  the count keeps) and the record now carries the counts the reviewer
  measured (318 lines: 22 blank, 65 comment-only). Under 700; the estimate
  was 48 % low, the second underestimate of the land-cover series (26a: 163
  against 170, fine; 36: 378 against 305). Nothing to rule; the reviewer's
  lesson ("a size explanation needs the grep that supports it") stands.
- *`@perf`'s evidence.* `docs/benchmarks/2026-10-12/26b-fractions/README.md`
  gives the tree, the `.so` hash, the machine, the power log and the
  foreign load (a Teams call, load average 5.3 to 7.6); gate 1 (exactness
  on the Velhas piece, all 19 classes inside the reference band, 0.77 s per
  million triangles at 10 m) is a measurement with its raw files; gate 2
  (no change to existing commands) is judged back to back against the base
  `5f28292c` because the stored baseline could not be used under the load.
  The main session's idle re-run (`15b77a67`) found Lagan `features read`
  +12 to +17 % on master too; `@perf`'s afternoon branch
  (`worktree-baseline-1012`, `c31e55e7`) traced it to the machine, not the
  code (the baseline's own commit and `.so` read 1.107 s against master's
  1.103 s), and declined to install a refreshed baseline measured on a
  loaded machine. The right refusal: a baseline 6 to 15 % slow would hide a
  regression of that size.
- *Out of role.* None in the personas' commits: every one carries its tag
  and touches its files. The main session's are D1 and D5.

**#233 in one line.** Prose: the Numedalslågen clip ruled expected at 68 %
in `hotspots.toml` and three ROADMAP rows to merged; two review rounds;
merged from the queue 10:46:35Z, 8 min 51 after enqueue.

**#234 in one line.** Prose: the morning retrospective; two review rounds,
the second closed on the reviewer's wording (D1); enqueued 10:30:09Z,
merged 10:31:30Z (prose-only, governance gates only).

**37, finished and pushed.** Three design rounds, three code rounds,
master merged in (`030fd32b`) with the one broken test fixed by `@tester`
(`d8a1729d`, the parameter count dropped for the identity set); every
gate green per the round-3 record; 328 net lines. Approved 12:15Z, pushed
as #236 at 14:38:03Z (2 h 23 on the push), queued 14:47:09Z. Section 3
has its lead time.

**The audit exists**: `docs/increments/generality-audit.md` on
`worktree-generality-audit`, one review round, three wordings and four
acronym expansions applied, two cites given their path (`1c405ade`).
Ola's to read.

## 2. Deviations judged (the main session's, as the brief asked)

- **D1, seven review rounds closed on the reviewer's exact wording with no
  further round**: ruling record round 2 (`7bba17f8`), retrospective round 2
  (`21acc6e1`), 37 design rounds 2 and 3 (`37f811c1`, `ba81cde1`), 26b code
  round 2 (`68d4b75b`), 37 code round 2 (`2a871c50`), the audit round 1
  (`8ced7ef2`). Each subject names the fix; several of the records say "as the
  recording-commit practice allows". This is P1 of
  `2026-10-12-merges-35-36-docs.md` applied on its default, seven times in
  one day, with no rule text yet, so three review records today restate the
  permission themselves. Judged consistent with P1 as proposed; the text
  should land (P1 below), and it should draw one line the morning's wording
  did not: 37 design round 2 went past wording. The main session wrote a
  design sentence (the child inherits the server's working directory) and a
  test-case specification (T6's relative-`out` case) into `@architect`'s
  file, verbatim from the reviewer. The sentence states an assumption the
  code depends on, which is design content; `@tester` then wrote the test
  itself in the red step (`ee5b5d52`), so no test was written out of role.
  Low risk, disclosed, and cheaper than a spawn. The morning's P1
  (`docs/retrospectives/2026-10-12-merges-35-36-docs.md@c317c7e4:250-253`)
  says "wording fixes"; P1 below rewords it to "text the reviewer gives
  verbatim" and leaves a decision, a type or a new invariant to the owner.
- **D2, one `@orchestrator` check for 26b and 37 together.** `CLAUDE.md`
  §3 spawns the check after each merge; 37 is unmerged, so this check reads
  it as finished work and the merge itself (queue, CI, merge commit) is not
  checked. One line in the next retrospective covers it. Not a deviation.
- **D3, the main session ran the quick check itself, three times.**
  `docs/increments/perf-quick-check.md@c317c7e4:262` and Ola's ruling of
  2026-10-07 at `:288` name it as the main session's job. Not a deviation.
- **D4, local time labelled UTC in records.** Section 14 of the increment
  file says the idle re-run was at "~12:10 UTC"; the commit that records it
  (`15b77a67`) is at 12:22:49 +0200, so the run was about 10:10Z. This
  check's brief says the red commit was at 07:19 UTC (git: 07:57Z; the
  `@tester` spawn was 07:22Z), 37's design at 07:58 UTC (git: 08:34Z; the
  `@architect` spawn 08:04Z) and its approval at 14:25 UTC (git:
  12:15Z), and `session.md`'s `NOW:` says "~15:50 UTC" on a file written
  before this check's 14:24:39Z spawn (transcript). A lead-time table built from those
  figures would be wrong by up to two hours. The mechanism: the session
  reads a wall clock in local time and writes "UTC". Proposal P7.
- **D5, an unpushed prose branch merged into the product branch.**
  `18a233ac` merges `worktree-hotspot-ruling` (#233, then review-complete
  and unpushed) into `worktree-mapbiomas-26b` so the Status line and ROADMAP
  row could read consistently (26b code review round 2 says why). The
  merge commit carries git's default subject, no "by the main session" and
  no persona tag, the one commit on 26b's own branch that names no actor;
  and #235's
  timeline lists #233's eight commits as its own. It coupled two PRs: had
  #233 been sent back, 26b would have carried its content. Expedient and
  harmless today; recorded so the next such merge says who made it and why
  in its subject.
- **D6, the dequeue, and the account of it.** The GitHub timeline of #235:
  `gh pr merge 235` ran at 10:29:51Z (transcript); GitHub enabled auto-merge
  at 10:36:10Z; the queue added #235 at 10:37:50Z, six seconds after #233 (10:37:44Z), and **removed it at
  10:38:01Z**, eleven seconds later and 8 min 34 before #233 merged
  (10:46:35Z). The main session noticed at 11:27Z, asked ("the queue
  dropped it after #233 merged ahead of it"), Ola said "yes, re-enqueue" at
  11:29:49Z, the queue took it at 11:30:05Z and merged it at 11:39:51Z. Two
  findings. (a) The account is wrong in its order: the drop came before
  #233's merge, not after it, so "#233 merged ahead" was an inference, not a
  reading of the timeline; lesson 6 of today (a merge-order claim is checked
  with the record, not a guess) applies to the main session's own
  sentence. The timeline does not say why the bot removed it. (b) **52
  minutes out of the queue, 49 of them unnoticed**: the session was on 37's
  fix chain and no process watched the queue. The same happened to #232
  last night (the main session's transcript at 06:00Z today: "dropped out
  of the queue with a conflict against #231"). Proposal P8.
- **D7, 37's four review fixes ran in series** (`@tester` 11:11Z,
  `@developer` 11:25Z, `@architect` 11:32Z; 38 min from the round-1 record
  `dbc77f3b` 11:04:32Z to the round-2 record `2a871c50` 11:42:48Z).
  `tools/brief.py` refuses two writers in one
  worktree by design; the parallel route is two branches and a merge, with
  a second build. On a 38-minute chain the merge would eat the saving.
  Judged right; no change.
- **Checked, no finding.** Ola's one word "push" at 10:26:12Z answered the
  main session's question that named the three branches, each as its own
  PR and enqueue; the three pushes, three `gh pr create` and three
  `gh pr merge` at 10:28Z to 10:29Z are the named acts, one occurrence
  each. The re-enqueue had its own yes. `@developer`'s `git commit --amend`
  at 08:33Z was refused (`queue-2026-10-12.jsonl`, "rewrites history"),
  correct, and not retried: the verification record went into section 14
  instead. No persona commit is out of role.

**Guards.** Three refusals today, all in the away window (07:43:12Z to
10:25:54Z, back on `away.py --back`). The amend above. Two `guard_governance`
refusals of reviewers' read-only Python heredocs: 07:47Z (the morning
retrospective's reviewer calling `guard_governance.governed()` on governed
path strings) and 09:42Z (26b's reviewer running `tools/check_citations.py`
through `subprocess` inside `python3 - <<'EOF'`). The hook's text rule
(`.claude/hooks/guard_governance.py@c317c7e4:91-95`) counts `python3 -` and a
heredoc as writing constructs, so any governed path in the body is a write;
`tools/check_citations.py` is governed by the `tools/check_` prefix (`:81`).
By `.claude/REQUIRED-READING.md`'s definition these are not false positives
(the named files are governed), and 26b's round-1 record labels the second
one "Guard false positive" all the same. The reviewer lost nothing (both
reran the check as a plain command). The guard cannot be taught to read
Python; the personas can be told the shape that passes. Proposal P9.

**Idle.** The away window had work running throughout: 23 commits from
07:44Z to 10:22Z (`git log --all` between the window's bounds), the longest
gap 29 min (1,746 s, between 37's design approval `0a03b224` and `@perf`'s
26b acceptance `98eea44f`); the window ended
with three items waiting on a push and 37 in its green step. No idle to
report. The ideas file was checked this morning; not repeated.

## 3. Lead time

| Item | First spawn | First commit | Approved | Pushed | Merged | Design or red to merge | Waiting |
|---|---|---|---|---|---|---|---|
| 26b (#235) | 07:22Z (`@tester`) | red 07:57Z | 09:54Z | 10:29Z | 11:39:51Z | **3 h 42** from red (4 h 17 from the spawn) | 35 min on Ola (away until 10:25Z); 64 min from auto-merge enabled (10:36:10Z) to merged, 52 of them out of the queue |
| ruling (#233) | - | 04:53Z | 07:57Z | 10:28Z | 10:46:35Z | 5 h 53 | 2 h 31 on the push, 9 min queue |
| retrospective (#234) | 07:38Z | 07:38Z | 08:03Z | 10:29Z | 10:31:30Z | 2 h 53 | 2 h 26 on the push, 1 min queue |
| 37 (#236) | 08:04Z (`@architect`) | design 08:34Z | 12:15Z | 14:38:03Z | unmerged at 14:43Z (queued 14:47:09Z; the Review section has the rest) | 3 h 41 design to approval | 2 h 23 on the push |

26b is the second item in two days with the shape to keep (35 was the
first): built, measured and reviewed inside Ola's away window, pushed four
minutes after his return. Its waiting was the queue, not Ola. 37 moved from
design to approval in 3 h 41 with no wait on Ola at all, then waited 2 h 23
for the push (14:38Z); the mechanism is the one the morning named (the away-mode push
grant, P2 of the morning, not built), only by day this time: Ola was
present and the session asked once, in a list of ten `ASK OLA:` lines.

## 4. Lessons handed back, recorded and judged

Disposition: *rule* (a line in a rule file, proposed below), *checklist*
(a line in a brief or persona checklist), *lesson* (recorded here only),
*drop*.

1. A red suite computes expected values from the inputs as rounded, not as
   intended (26b C++ case 6: 0.5 against a triangle whose area as doubles is
   4.9999997e-19; fixed to the shoelace area of the rounded corners,
   `f165282b`). The object-identity question of `.claude/REQUIRED-READING.md`
   already covers it. **Lesson.**
2. The red step carries the layering rows and the reader support for any
   file structure it asserts (the `.ply` list reader and four
   `test_layering.py` rows came in `f165282b`, after green). **Checklist**,
   `@tester` brief wording; no rule text.
3. An amend of an unpushed commit is refused in away mode, so the record
   goes into the as-built section. Already the rule (a history rewrite needs
   a yes) and already the practice. **Drop.**
4. Check the 1-minute load average, not only `ps`, before a timing run: a
   video call shifted every phase 10 to 30 % (`26b-fractions/README.md`),
   and interactive use alone shifted serial C++ 15 %
   (`baseline-refresh/README.md`, run 1 and 2). **Rule**, one line in
   `perf.md` (P10).
5. A percentage in a comment names the README line it sums from. **Lesson.**
6. A merge-order claim is checked with the commit's parents, not merge
   timestamps; today's D6 adds: or with the forge's timeline, never from
   memory. **Lesson.**
7. A cite written without its directory path breaks the citation gate, and a
   review record quoting it becomes an at-risk entry; four times today, each
   caught by a reviewer; the main session now runs the worktree's own
   `check_citations.py` after each recording commit. The hook,
   `gates_after_commit.py`, gates the worktree the commit line's `cd` or
   `git -C` names (`.claude/hooks/gates_after_commit.py@c317c7e4:77-87` and
   `:111-131`), and it did its job on this file: the unpinned
   `PRINCIPLES.md` cite in section 6 was caught after the first commit
   `f4ebb3e5` and fixed in `1d3fb3b6`, not an amend. What is still open is
   narrower: P2 of `2026-10-11-merges-32-33-34.md` (approved, unbuilt) has
   `step()` check that a joined relative path is a directory before the
   hook gates it. **Carried proposal**, nothing new.
8. The citation gate's at-risk tally moves with every file master later
   touches; a tally names its command and tree. Done in #234. **Lesson.**
9. An assumption named only in a handback is not in the design (37 design
   round 2). **Lesson**; the increment README's "the increment files are the
   source of truth" says it.
10. Help text is a rendered claim; check it as the default renderer shows
    it (Rich ate `[mcp]`, 37 code round 2). **Lesson.**
11. A design claim about existing CLI wording is checked by running the
    command before a test pins it (37's T6 "naming `--out`"). The
    claims rule already binds it. **Lesson.**
12. A design that names a README section needs the as-built record to say
    it landed (37 code round 1). **Lesson.**
13. A pinned snapshot of another module's shape (a CLI's option count) in a
    derivation suite is a resolved value where a rule belongs (`d8a1729d`).
    `.claude/REQUIRED-READING.md` already says "write the rule and the
    command, not the resolved value". **Lesson.**
14. A brief saying "both red" is checked against the code; a case pinning
    current behaviour is a regression guard, and the handback says which.
    **Checklist**, `@tester` brief wording.
15. Small-model summaries of SDK documentation were wrong twice; introspect
    the installed package (37 design). **Lesson**, for `@architect`'s
    practice: a `WebFetch` summary is a lead, the package is the source.
16. `docs/` is ruff-linted; a probe under it needs a file-level `noqa` with
    its reason. **Lesson.**
17. A venv with an editable install imports its own checkout even with
    `sys.path` set; assert `tin_engine.__file__` in every probe.
    **Lesson** (the object-identity question again).
18. The full Python suite takes 12 minutes, past the 10-minute foreground
    limit; run it in the background from the start. **Lesson**; it belongs
    in the brief of any persona that runs the full suite, not in a rule.
19. `--stats` takes an optional path, so `--stats --binary` writes a file
    named `--binary`. A product defect, not a lesson: **to Ola's ideas
    file** (handback).
20. The two heredoc refusals, judged above under Guards. **Checklist** (P9).
21. `tools/brief.py` refuses two writers in one worktree, so 37's fixes ran
    in series (D7). **Lesson**; the series was the cheaper route.
22. The queue dropped #235 and the re-enqueue needed a fresh yes (D6).
    **Rule or tool** (P8): the yes is right (one yes, one enqueue); the 49
    unnoticed minutes are the defect.
23. A claim about an unmerged item carries the time it was true: this file
    said 37 was "waiting on a push" while #236 had been pushed five minutes
    before its first commit (14:38:03Z against `f4ebb3e5` at 14:43Z), and
    both #236 and #237 merged before its review. The tree moves while a
    retrospective is written. The reviewer's lesson; the claims rule's
    "write the rule and the command, not the resolved value" covers it,
    with a date where a resolved value is unavoidable. **Lesson.**

## 5. The morning's proposals: which took effect today

| Proposal (`2026-10-12-merges-35-36-docs.md`, section 8) | State on `c317c7e4` |
|---|---|
| P1, the recording commit applies the reviewer's named wording | applied in practice seven times (D1); no text in `docs/increments/README.md` (`git log --since=2026-10-12 -- docs/increments/README.md`: nothing) |
| P2, the away-mode push grant | not built; no commit to `.claude/hooks/` or `tools/away.py` today; 37 waits on a push by day too |
| P3, `check_citations.py` prints drifted cites, counts the rest | not built; no commit to `tools/check_citations.py` today |
| P4, a `mesh_sha256` change refreshes the baseline in the same PR | not applied (`perf.md` unchanged); the refresh itself was attempted and held (`c31e55e7`) |
| P5, the -155 word cut | not applied (`rule_sizes.py`: 0 change) |
| P6, carry 2026-10-11 B1 and P2 | carried; B1 still unapplied (this check's brief says "none named"); P2 unbuilt (lesson 7) |

The mechanism is the morning's: product work filled the three slots all
day (26b, 37, the audit), and `CLAUDE.md` sends a free writer to product
first. Eight approved harness items now wait (P1 to P5 above,
2026-10-11 B1 and P2, 2026-10-07 P8).

## 6. Rule text: size, anecdotes, repeats

`python3 tools/rule_sizes.py` on `c317c7e4`: **10,935 words, +0** since this
morning (`a749eab7`); against the 2026-10-03 baseline (`cc5c85b4`, 9,773
words) **+1,162**, +944 like for like. Not one word moved in a day with
seven review rounds that each restated a permission the rules do not hold.

**Line-level pass.**

- *Carried, unapplied:* the morning's three cuts (-155: the pipeline stated
  twice in `CLAUDE.md`, the branch-protection history there, the
  `guard_governance.py` test cases restated in `REQUIRED-READING.md`).
- *New, a worked example in a rule:* `.claude/REQUIRED-READING.md@c317c7e4:79-81`
  spells the gate-planting check as a recipe for one gate ("append the
  prohibited directive to `CMakeLists.txt`, run `check_prohibited_deps.py`,
  confirm it fails naming it, restore the file"). The rule is the first
  clause ("verify a gate by planting what it forbids"); the recipe is one
  case of it and belongs with the tool (`CLAUDE.md` §4 already names it),
  or nowhere. About -30 words.
- *Restatement outside the rules, counted for P1:* "as the recording-commit
  practice allows" (`git grep -l "recording-commit practice" <tree> -- docs/`)
  stands in two files on `c317c7e4` (increment 35's file from yesterday
  and the morning retrospective), plus increment 37's file on
  `worktree-mcp-server` at `33061c0f` and the audit on
  `worktree-generality-audit` at `1c405ade`: four in all, and all four on
  master once #236 and #237 merge. Each is 20 words a
  one-clause rule would make unnecessary.
- *Checked, not an anecdote:* `docs/PRINCIPLES.md@c317c7e4:71` and `:118` name
  retrospective files as pointers to reasoning, which is what a pointer is
  for.

**Proposed cut this round (default: yes):** the new -30 with the carried
-155, about -185 words, to about 10,750.

## 7. Proposals (Ola decides; defaults in brackets)

1. **P1 (carried, reworded)** [yes]: `docs/increments/README.md` step 4
   gains "the recording commit applies text the reviewer gives verbatim,
   each named in the subject; a decision, a type or a new invariant goes to
   the owner", `@architect`, about +30 words. D1.
2. **P2, P3, P4, P6 (carried)** [yes, unchanged]: the push grant, the
   drifted-cite printer, the baseline-refresh line, 2026-10-11 B1 and P2.
   Section 5 says none moved.
3. **P5 (carried) plus the new -30** [yes]: -185 words, `@architect`.
4. **P7, times carry their zone from the command** [yes]: the recap and
   every record take their time from `date -u` or git's `%ci`, never from
   a wall clock plus the word "UTC"; one line in `.claude/briefs/common.md`
   (+15 words) and `session_state.py` prints `date -u` in the recap header.
   D4.
5. **P8, watch the queue after an enqueue** [yes]: `gh pr merge` is followed
   by a background loop (`gh pr view <n> --json state,autoMergeRequest`
   every 60 s) that prints one line when the PR leaves the queue unmerged;
   a shell loop or an ungoverned `tools/queue_watch.py`, about 30 lines,
   `@developer`. D6: 49 min today, #232 yesterday.
6. **P9, governed paths never go inside a heredoc or `python3 -c`** [yes]:
   one line in `.claude/briefs/common.md`: run the tool as a plain command
   and put any script in the scratchpad. +25 words. Guards, lesson 20.
7. **P10, the load average before a timing run** [yes]: one line in
   `perf.md`: `uptime`'s 1-minute load average is recorded with `pmset`,
   and a run under a load above the baseline's is judged back to back or
   not at all. +25 words, `@perf`. Lesson 4.

## Review

Round 1, 2026-10-12, `@reviewer` (`c317c7e4..1d3fb3b6`): CHANGES REQUESTED,
six wordings, applied in this round's recording commit (the one after
`1d3fb3b6` on `post-merge-1012b`). The reviewer verified every timestamp,
LOC figure, rule size, disposition and citation in the file; the six are
pinned at `1d3fb3b6`, and each was re-checked here before it was applied
(`gh api repos/expertanalytics/rasputin/issues/235/timeline`, `gh pr view
236`, `git log -1 --format=%ci`, `git show` of the hook and the two cited
files).

1. `/Users/skavhaug/projects/rasputin/.claude/worktrees/post-merge-1012b/docs/retrospectives/2026-10-12-merge-26b-and-day.md@1d3fb3b6:135`
   said `gh pr merge` enabled auto-merge at 10:29:51Z; the timeline has
   `auto_merge_enabled` at 10:36:10Z, so D6 now gives both times, and the
   lead-time row at `:191` reads 64 min from auto-merge enabled, not 70.
2. `/Users/skavhaug/projects/rasputin/.claude/worktrees/post-merge-1012b/docs/retrospectives/2026-10-12-merge-26b-and-day.md@1d3fb3b6:151`
   and `:154`: `dbc77f3b` 11:04:32Z to `2a871c50` 11:42:48Z is 38 minutes,
   not 48.
3. `/Users/skavhaug/projects/rasputin/.claude/worktrees/post-merge-1012b/docs/retrospectives/2026-10-12-merge-26b-and-day.md@1d3fb3b6:1`,
   `:8`, `:77-78`, `:194`, `:199`: 37 was pushed as #236 at 14:38:03Z
   (queued 14:47:09Z), before this file's first commit; "waiting on a push"
   dropped for 37, "pushed 14:38Z, 2 h 23 on the push" written, the claims
   dated (lesson 23).
4. `/Users/skavhaug/projects/rasputin/.claude/worktrees/post-merge-1012b/docs/retrospectives/2026-10-12-merge-26b-and-day.md@1d3fb3b6:234-236`
   said the hook judges the main checkout's tree; it gates the worktree the
   commit line's `cd` names and caught this file's own cite, and the open
   2026-10-11 P2 is about `step()` checking a joined relative path.
   Reworded.
5. `/Users/skavhaug/projects/rasputin/.claude/worktrees/post-merge-1012b/docs/retrospectives/2026-10-12-merge-26b-and-day.md@1d3fb3b6:26`:
   "section 10" is section 9 of increment 37's file.
6. `/Users/skavhaug/projects/rasputin/.claude/worktrees/post-merge-1012b/docs/retrospectives/2026-10-12-merge-26b-and-day.md@1d3fb3b6:103-105`:
   the "verbatim text, type, invariant" phrasing is this file's P1, not the
   morning's, which says "wording fixes"; now "P1 below", with the
   morning's cite.

Two suggestions taken: the spawn time at `:120` (14:24:39Z, transcript) and
the grep at `:316`, which the reviewer counted at three files on master;
re-run here it is two on `c317c7e4` and four on `origin/master` at
`5d2e706d`, after #236 (merged 14:56:10Z) and #237 (the audit, pushed
14:53:40Z, merged 14:57:17Z) landed, so the file now names the trees. Those
two merges also date this file: its title and sections 1 and 3 are as of
14:43Z, and the next check covers both merges (D2).

One finding from running the gate on this record. `python3
tools/check_citations.py` was green before and after the six cites above
were planted wrong (`@1d3fb3b6:9135` in a 353-line revision): its `PATH`
pattern starts with `(?!/)` and its comment says a token starting with `/`
is "never a citation" (`tools/check_citations.py@c317c7e4:62-64`). A
relative cite planted wrong the same way (the hook cite in lesson 7 with
its line changed to 9077) was caught at once ("has 227 lines"). So the full-path
pinned citations the review records are asked to carry are not checked by
the gate at all; "every cite resolves" is true only of the relative ones.
On `c317c7e4`, 251 such cites in 20 files under `docs/` (`git grep -o
'/Users/skavhaug/projects/rasputin/[^ \`]*@[0-9a-f]\{7,\}:[0-9][0-9-]*'
c317c7e4 -- docs/`) have never been checked. Proposal P11 in the handback.

Round 2, 2026-10-12, @reviewer (`1d3fb3b6..49dc4ddd`): CHANGES REQUESTED, one line, applied in this recording commit as specified: /Users/skavhaug/projects/rasputin/.claude/worktrees/post-merge-1012b/docs/retrospectives/2026-10-12-merge-26b-and-day.md@49dc4ddd:418 said "a 410-line file"; the file at `1d3fb3b6` has 353 lines, so "a 353-line revision". The six round-1 wordings verified against their sources; the grep counts hold; the 251 absolute cites in 20 files confirmed, and planting a broken absolute and a broken relative cite in a scratch copy reproduced the skip (absolute silent, relative caught). Ideas for the owner at the next touch: cite `tools/check_citations.py@c317c7e4:62-65` (the pattern is line 65); "2 h 23" is 2 h 22 to the second; the pipefail lesson belongs in the lesson list; the fix should strip the worktree prefix, since the review records require full paths, and expect the 251 newly checked cites to yield a tally. Pushed on the reviewer's wording without a third round (the recording-commit practice, asked of Ola 2026-10-12, default yes).
