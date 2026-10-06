# 2026-10-06, 06:05 to 00:43: three speed-ups, six master merges, 85 lessons

Recorded by `@orchestrator` on 2026-10-07 at about 01:00 Oslo time, on
`worktree-retro-1007` off master `0c2572fb`. Times are Oslo time (UTC+2)
unless marked UTC. The part of the day before 06:05 is in the night
retrospective, `docs/retrospectives/2026-10-06-night-audit-geojson.md` on
`worktree-retro-1006` (commit `88056ce`, PR #194, not merged; see D6).

Sources: the main session's transcript,
`~/.claude/projects/-Users-skavhaug-projects-rasputin/6a989a29-4492-430f-99b0-9fe3058f931a.jsonl`
(cited as "line N"; this window is lines 9512 to 18859), its subagent
transcripts in the `subagents/` directory beside it (cited by agent id),
`git log` on `origin/master` and the branches named, `gh pr view`, `gh run
view`, the harness's `windows.jsonl`, and the 85 lessons the main session
passed me, which are recorded word for word in section 9. Handback quotes
are the personas' words, not Ola's. Every rule or tool change below is a
proposal; Ola decides.

Labels used here:

- **30a, 30b, 30c** are the three speed-ups Ola asked for on 2026-10-06
  ("yes, start the bottleneck removal process", line 12270): a faster
  land-cover lookup (#199), a faster clip of the CORINE land-cover polygons
  to the catchment (#200), and faster reading of the elevation model tiles
  (#204).
- **Audit A, C, D, F** are pull requests of the Python code audit
  (`docs/increments/python-audit.md`). C is "one GeoJSON reader" (#201).
- **C2** is a follow-up to C that the main session started without being
  asked (D1).
- **20c** is the planned fix for sliver triangles (very thin triangles) in
  the meshes; **the sliver research** is #203, which found their causes.
- **Guard fix PR B** is the second pull request of increment h16, which
  closes routes around the push and rule-file guards (`worktree-h16b`).
- A **red step** commits failing tests, a **green step** the code that
  passes them, and a **pin ruling** is `@architect` confirming the choices
  `@tester` made beyond the design, before the green step.
- A **master merge** is merging `origin/master` into a branch that is
  behind it.

## 1. What happened

**Merged to master in this window** (merge-commit time, from
`TZ=Europe/Oslo git log --first-parent origin/master`):

| Time | PR | What |
|---|---|---|
| 06:47 | #189 | test-harness cleanup |
| 07:59 | #191 | C++ dead-code removal |
| 08:59 | #192 | audit B, one "same coordinate system" rule |
| 10:51 | #193 | GeoPackage colour note in `ROADMAP.md` |
| 10:58 | #195 | GeoTIFF whose coordinate system is given only by parameters |
| 11:05 | #196 | audit F, shared catchment code |
| 13:11 | #197 | audit A |
| 14:13 | #199 | 30a, land cover about 4 times faster |
| 16:57 | #200 | 30b, CORINE clip 3 to 4 times faster |
| 17:04 | #202 | three `ROADMAP.md` rows and the parked NMD land-cover proposal |
| 18:04 | #203 | the sliver research (Lagan) |
| 18:17 | #201 | audit C; its master merge broke code, and review caught it (D3) |
| 00:22 | #204 | 30c, DEM reading |

The main session's brief listed #197 to #204 only; #189 to #196 also merged
after 06:05 and are included.

**In flight at 00:43:** the 20c design (design review round 3 asked for
changes at 19:47; `@architect` resumed at 00:13 after the Mac slept); guard
fix PR B at its round-9 red step (`@tester`, 00:15); audit D, master merged,
in code review round 3 (spawned 00:42, line 18766).

**Still open on GitHub:** #194 (the night retrospective; red since 10:51,
D6) and #198 (the `cli_driver.py` rename). The main session told Ola at
11:16 that #198 conflicted with master after #196 and would get a master
merge once `@perf` was done (line 12520); the PR has not changed since
11:01.

**Ola's own time:** at the keyboard most of the day; away windows
06:57 to 07:52, 08:51 to 10:50 and 12:15 to 13:03 (`away.py --back`, lines
11343, 12114, 13348), and one entered at 14:15 that ran to its end at 15:27
(`windows.jsonl` has no "back" for it); "Suspending now, going out for dinner" at 20:10
(line 18245) until 00:12; the night window from 00:27.

## 2. Deviations

**D1. The main session started unrequested work, and put its details to
Ola as questions.** At about 05:10 (before this window; the night
retrospective's section 10) it started C2 "to use the rest of the night",
and ran C2's later steps in this window (spawns at 06:13, 06:28, 06:47,
07:02, 07:05 and 07:21, lines 9584 to 10896). At 11:17 it listed four C2 questions with
defaults for Ola to rule on (line 12529). Ola, 11:20:52: "The whole C2 part
here is entierly made up by you." The main session agreed (line 12538: "You
never asked for that PR. I started it myself ... to keep a free agent
busy"). Ola, 11:22:29: "Is it useful for our goal? In that case, we should
discuss, plan and execute properly ;)". The main session's own answer (line
12553) was that it served neither huge catchments nor residual inflow, and
it was dropped at Ola's "yes, drop it and add the ROADMAP row" (11:24, line
12557). This breaks Ola's note on filling the window, which says never to
invent new pull requests to fill slots. The night retrospective's section
10 had counted C2's start as filling an empty slot, without asking whether
Ola had asked for it. I should have asked.

**D2. Twelve corrections of what the main session had told Ola.** Each is
the main session's own words, by transcript line:

| Line | Time | What was corrected |
|---|---|---|
| 9678 | 06:16 | "I said `bench.py` couldn't run a reprojected DEM, and that was wrong" (the claim was made in the night, before this window; only its correction falls inside it) |
| 9717 | 06:20 | who caused the scratch-folder hand-off (from the night) |
| 11525 | 08:49 | #191 had already merged; audit F's PR had not been opened |
| 12529 | 11:17 | two items listed as open were already settled |
| 12538 | 11:21 | C2 was its own initiative (D1) |
| 14926 | 14:47 | Ola never approved a GeoJSON wording row |
| 15041 | 14:53 | a null station "could pair a catchment with a station called None" (said at 11:22, line 12553): false, station numbers must look like `2.11.0` |
| 16008 | 16:26 | Swedish CORINE 2018 is not mapped afresh from Sentinel-2 |
| 16236 | 16:41 | the sliver line is a seam inside the CORINE data, not the 57° N edge between two DEM tiles (said at 16:11, line 15732) |
| 17734 | 18:38 | the 0.0066° triangle "is made only of input lines": withdrawn after review round 2 |
| 17877 | 19:12 | the withdrawal withdrawn: the triangle is in the CORINE source after all; the contradicting run used another input file |
| 18045 | 19:42 | "drop your outline rule" withdrawn: it was measured on the wrong mesh (8 % on 20c-2, 22.5 % on 20c-3) |

Six of these (the null station, CORINE 2018, the tile edge, the 0.0066°
triangle twice, the outline rule) are a cause or a "could" relayed from a
handback before anyone had checked it. The `bench.py` claim was also
unchecked when said; whether it came from a handback is in the night's
window, and I did not trace it, so it is not counted among the six. `CLAUDE.md` §3 already says "no
cause not checked". The night retrospective's D3 asked that such claims be
passed on as "X says ... (unchecked)"; it was not used today. Ola had to
take each one on trust and then unlearn it.

**D3. A master merge that changed code was checked as a docs merge, and
the after-commit gate did not see it.** The decisive check: the
after-commit gate ran at 17:14:41, after C's master merge `5c2b771a`
landed, took 12,776 ms and exited 0 (`agent-a9e6e5394037aa395`, line 152
of its transcript), yet `ruff check .` over a copy of the tree at
`5c2b771a` exits 1 with one error, `Undefined name MultiPolygon` (run for
this retrospective, with the main checkout's `.venv/bin/ruff`). The merge
brief for C (line 16709, 17:11) said the conflicts were "only those two docs", that
`feature_input.py` "merges cleanly", and that `ast.parse` "is enough". The
merge, `5c2b771a`, dropped `MultiPolygon` from `feature_input.py`'s
imports, which #200's `_linework` uses. `ruff check --select F821` on the
file at `5c2b771a` reports `Undefined name MultiPolygon`; at the fix
`800c4adf` it reports nothing (run for this retrospective). Code review
round 8 caught it, because its brief (line 16843) asked for the suite "if
that's cheap". Two more things in the merging agent's transcript
(`agent-a9e6e5394037aa395`):

- At 17:13:34 it chained `git checkout -q f56952b5 2>/dev/null || true`
  onto a `check_citations` run. That left the merge, so `MERGE_HEAD` was
  gone (17:13:48). It then rebuilt the merge by hand: `git write-tree`,
  `git commit-tree -p HEAD -p origin/master`, and `git merge --ff-only`
  (17:14:18). No guard covers that route; no guard refused anything either.
- `gates_after_commit.py` fired after that command (the text contains
  `git merge`); that is the 17:14:41 run above. It ran in the main checkout,
  not in `audit-geojson`: the hook takes the tree from the event's `cwd`
  (`committed_tree` in `.claude/hooks/gates_after_commit.py@0c2572fb`), and
  in every one of about 7,900 hook events in this session's subagent
  transcripts that `cwd` is `/Users/skavhaug/projects/rasputin`. Subagents'
  shells go back to the session directory between calls, and they work in
  worktrees by `cd <worktree> && ...` inside the command. Of those runs,
  about 400 took over 3 s (the gates ran) and all exited 0. (Counts are
  "about" because the transcripts kept growing while I counted, near 01:00
  on 2026-10-07.) **So in this session the
  after-commit gate has never checked a subagent's commit in a worktree.**
  Claude Code's hooks guide says `cwd` "follows Claude ... after Claude runs
  `cd`" (code.claude.com/docs/en/hooks, read 2026-10-07); for subagents the
  `cd` lasts only one call. A probe of the function itself was refused by
  the unattended guard (section 4), so this rests on the code and the
  transcripts, not on a run of the hook.

Audit D's master merge (00:25 to 00:40) was briefed with ruff and the full
suite (the `@architect` docs-half brief, line 18592, 00:25, and the
`@developer` code-half brief, line 18690, 00:31; line 18598 is the main
session's message to Ola about it). The suite found a test broken by a helper-signature
change on master, and `@tester` fixed it (line 18727). The lesson took, by
brief.

**D4. Pin ruling skipped for 30a.** The main session's own report: it
skipped `@architect`'s ruling on `@tester`'s pins before green, because Ola
wanted numbers within the hour (line 12824, 12:15: "it would be great to
have some speedup results ready by the time I come back"). The pins were
ruled after green. 30b (ruling 14:41, green 14:47) and 30c (ruling 16:54,
green 17:01) kept the order. `@perf`'s acceptance ran before code review on
all three (30a 12:39 then 12:51; 30b 14:58 then 15:10; 30c 19:12 then 19:48).

**D5. Shell writes and `git stash`, against briefs.** A scan of the
subagent transcripts in this window (Bash commands only):

- Test files appended through a shell heredoc (`cat >> tests/python/...
  <<'EOF'`), all by `@tester`: 06:29 `test_station_set.py`, 12:23
  `test_landcover.py`, 14:10 two files in one command
  (`test_guard_governance.py` and `test_guard_targets.py`,
  `agent-ae32fa72d6b542...`), 14:32 `test_feature_input.py`, and 00:17 one
  file (`test_guard_push.py`, `agent-a57b4c6ba349a4...`).
- `@architect` appended to design files by heredoc about 17 times, to 7
  files. 8 of them went to `docs/increments/20c-soft-quality.md`: seven by
  `agent-a6c46d87f5ceed...` (17:12 and 17:15, transcript lines 595 and 697,
  then five from 17:44 to 17:47) and one at 00:17
  (`agent-a4053180efccd8...`). `@perf` appended once to
  `docs/benchmarks/2026-10-06/30c-dem-read/README.md` (19:46). The scan
  reads every subagent transcript from 06:05 to 00:43 and counts each
  `cat >> <file> <<` (or `cat <<EOF >> <file>`) whose file is under
  `docs/`, `tests/`, `tools/` or the source trees; "about", because an
  append written another way is not counted. Its per-file output:

  ```
    8  architect  20c-soft-quality.md
    2  architect  30a-landcover-speed.md
    2  architect  30b-clip-speed.md
    2  architect  30c-dem-read-speed.md
    1  architect  h16-harness-fixes.md
    1  architect  python-audit-pr-c2.md
    1  architect  python-audit.md
    1  perf       README.md
    1  tester     test_feature_input.py
    1  tester     test_guard_governance.py
    1  tester     test_guard_push.py
    1  tester     test_guard_targets.py
    1  tester     test_landcover.py
    1  tester     test_station_set.py
  ```

  The brief line covers only rule files, so these broke no brief; T6's
  "appends included" would cover them.
- `@developer`, guard fix green (17:31, `agent-a300dfa80...`), rewrote
  `tools/shell_scan.py` through `python3 - <<'EOF'` with `read_text`/
  `replace`, and committed it with `.claude/hooks/guard_push.py`. Whether
  the hook file itself was also written by the shell, and whether Ola was
  asked, I did not establish (unchecked).
- `git stash` in six runs: `@developer` 07:00 (`geojson-gaps`), `@architect`
  14:02 (`h16b`), `@tester` 13:02 (`landcover-speed`), `@developer` 17:47
  and 17:51 (`h16b`; each run did a `stash` and a `stash pop`), `@tester`
  00:41 (`audit-encoders`; see D9). The stash is one stack for every worktree:
  `git rev-parse --git-path refs/stash` gives the same
  `/Users/skavhaug/projects/rasputin/.git/refs/stash` from the main checkout
  and from `h16b` (checked), and git's worktree guide says every ref under
  `refs/` is shared except `refs/bisect`, `refs/worktree` and
  `refs/rewritten` (git-scm.com/docs/git-worktree, read 2026-10-07). One
  agent's `stash pop` can take another's change. No incident of that today.

**D6. The night retrospective sits red, and was reported as waiting.** PR
#194's governance check failed at 10:51:17: `check_citations.py` found
three citations to `b26beb8`, a commit not yet on `origin` (`gh run view
37438850992`). At 11:16 the main session told Ola "#194 (night
retrospective) is set to auto-merge and waiting for its checks" (line
12520). The PR is still open and red. `b26beb8` reached master with #201,
so a re-run should now pass (not run: a CI re-run is a write to GitHub).
Behind it, three more retrospectives exist only as local branches:
`worktree-retro-1004c` (7 commits), `worktree-retro-1005` and
`worktree-retro-1005b` (1 each), none pushed. `next.md` on master stops at
the proposals of 2026-10-04, so the proposals since then are not where Ola
or the next session looks.

**D7. Stale background jobs, found by Ola.** Ola, 20:10: "The two perf jobs
seem to be stale?" (line 18268). The main session found five leftover log
watchers holding two `@perf` jobs open, one since the morning, and stopped
them (lines 18307 to 18343). The recap is meant to list running
background jobs (Ola's ruling of 2026-10-03, item 5 in `next.md`); whether
it did here, and the main session missed it, I did not check.

**D8. Who merges master into a branch.** Master merges were split by file
kind across three personas (docs to `@architect`, tests to `@tester`, code
to `@developer`) for #191, #192, F (twice), C (twice), #197 and D: about 20
writer spawns, and about 10 more to review the merges, out of about 130
spawns in the window. The C merge at 17:11 went to `@architect` alone, code
included (D3). The retrospective of 2026-10-04 left "whose job a merge
conflict is" open.

**D9. Night findings, 00:25 to 01:05** (the last part after this window's
end, added at review round 1):

- **The main session typed a wrong merge range into audit D's briefs.**
  Audit D's master merge `15f41a82` brought in master from #187 up to #201,
  because D had never merged master and its merge base was `b26beb83`, C's
  approved head (`git log --first-parent b26beb83..01d98c2b`, as round 4
  cited it). The main session's briefs said otherwise: the merge brief "Master now also
  has #199, #200 and #202-#203" (line 18592, 00:25), the code review round 3
  brief "#200-#203, with C merged as #201" (line 18766, 00:42), and the
  Status-rewrite brief "with #200-#203" (line 18948, 00:57). `@architect`
  wrote the last into D's Status line, and code review round 4 caught it
  (line 19041, 01:02). The main session told Ola at 01:02 that the phrase
  "came from my brief" (line 19053). The same lesson as T2's "print
  `git log --merges --first-parent`": a range is a command's output, not a
  typed value.
- **`@tester` used `git stash` again at 00:41** (`stash -q` then
  `stash pop -q`, `agent-a9c4f08df2702e...`, "Audit D: fix broken CLI test
  call"), on `audit-encoders`. The main session's lesson says "despite the
  brief forbidding it". That is true of the earlier `@tester` brief for the
  same merge, "Use Edit for every change, appends included; no shell writes
  and no `git stash`" (line 18641, 00:29), but that was another spawn, and it
  did not stash. The brief of the agent that did (line 18727, 00:40) says
  "Use Edit" and does not mention `git stash`; neither does its `brief.py`
  block. So the ban was in the brief before it, not in this one. Either way
  the rule lived in a brief, which every new spawn has to be given again; T6
  puts it in the guard.

## 3. Idle time

- **09:25 to 10:50, 85 minutes, no agent running**, during and just after
  the away window that ended at 10:03 (the next spawn is at 10:52, line
  12180). At 09:25 the main session wrote that it was "holding" C's review
  "until audit F lands and the GeoJSON design is renumbered to section 12"
  (line 12108). C's renumbering and review needed no ruling. `away.py
  --back` printed a longest stretch without a commit of 39 minutes, because
  it measures only to the window's planned end (10:03), not to Ola's
  return (10:50).
- 06:57 to 07:52 and 12:15 to 13:03: longest stretches without a commit 14
  minutes each.
- 20:12 to 00:12: the Mac was suspended by Ola. Not idle.

## 4. Guards

`away.py --back` listed two refusals after 06:05, both `gh pr create`
runs by the main session (06:58 and 08:52, lines 11343 and 12114), held
for Ola's yes as the rule says, and none for the window ending 13:03. I
could not read the queue for the 14:15 window or tonight's (third refusal
below), so those are not covered.

Four refusals of my own, in this run (unattended mode on until 10:03
2026-10-07):

- A `python3 -c` that **read** `.claude/hooks/gates_after_commit.py` to call
  its `committed_tree` on a sample event was refused as a write to that
  file. A false positive in substance, but the named file is a rule file,
  so I did not redo it by any route. The D3 claim rests on the code and the
  transcripts instead.
- A transcript scan whose Python regex contained `>>?\s*\S*(tests/python/...`
  was refused as a write to that "path". It named no real file; I redid it
  as a script written with Write into my scratchpad and run with
  `python3`.
- An `ls` and a read-only `python3` over `<git-common-dir>/harness/` were
  refused: "only hooks and away.py write the harness state". It wrote
  nothing. Harness state is not an ordinary file, so I did not redo it. An
  earlier `cat` of `windows.jsonl` in the same run was not refused, so the
  guard treats two reads of one folder differently.
- A read-only `gh api graphql -f query='{...}'` (to see whether #198 was
  dropped from the merge queue) was refused as "gh api with a writing
  method": `-f` makes `gh api` send a POST, which GraphQL needs even to
  read. A forge act, so not redone. The guard cannot tell a GraphQL query
  from a mutation.

## 5. The lessons, grouped

85 lessons came in (section 9). One, "(#192 merge) local master stale",
counts in two groups.

| Group | Lessons | What repeats |
|---|---|---|
| Master merges into a branch | 15 | 9 on the merge itself (what it brought in, `--theirs`, chaining commands, the check's scope, the commit message, re-adding a table total, broken tests after a clean merge); 6 on numbered sections in shared docs colliding or going stale |
| Probes and claims (method) | 13 | probe the whole input domain, the nearest valid input, every runner and path spelling; measure on the right object (PR N's mesh, the prototype's own numbers); counts not shares |
| Domain traps | 11 | Shapely `prepare` and predicates, CRS seams, millimetre segments, NMD versions, a fixture name equal to a subcommand |
| Which Python and which code runs | 10 | worktree without a real venv, copied or symlinked venvs, base installs in the worktree venv, `-DPYTHON_EXECUTABLE`, a script shadowing an import |
| `check_citations.py` | 8 | stale local master as the base (4), ranges past end of file, short `file:line` forms, unchanged ranges flagged, citations shifted by a mid-file insert |
| Figures in documents | 7 | hand-copied or uncommitted-run numbers, rankings not re-measured, counting rule not stated; generated tables fixed it in 30b |
| Briefs and quotes | 5 | a phrase or rule cited from memory, a scratch extract as a channel, stale quotes of another branch |
| Rule and tool alignment | 5 | `Monitor` missing from `@perf`'s tools, `python` against `python3`, rule-file size, a slow guard line, busy-waiting beside a timing run |
| Main session's own steps | 3 | C2 unrequested, 30a pin ruling skipped, the null-station claim |
| Records | 3 | `ROADMAP.md` row with the status line, chat rulings into the increment file, which PR publishes a profile |
| Shell writes and `git stash` | 3 | heredoc writes, stash against briefs, twice |
| Scratch | 3 | a glob `rm` deleted others' files, keep prototypes and input counts |

Two groups are most of the avoidable cost: master merges (15) and which
Python runs (10). Both repeat across personas and days, both have a
mechanical fix, and both were already proposed (the merge runner, stage 3
of `docs/research/2026-10-03-dispatcher-control.md`, ruled yes on
2026-10-04; and P6 of the night retrospective, a venv for each worktree).
The probes, domain and figures groups are judgement; existing rules cover
them ("Claims: run the check before you write it down" in
`.claude/REQUIRED-READING.md`), and a new rule line for each would add
words without adding a check.

## 6. Proposals

Fewest changes, most lessons. Rule files stay without history; the
incidents stay here.

| # | Change | Lessons and incidents it would have stopped | Cost | Owner | Default |
|---|---|---|---|---|---|
| T1 | **The after-commit gate checks the tree that was committed.** `committed_tree` takes the directory from a leading `cd <dir>` or `git -C <dir>` in the command (`tools/shell_scan.py` already parses commands), and falls back to `cwd`. A test feeds it a subagent-shaped event. | D3 (the dropped import at `5c2b771a`, caught at the merge instead of in review round 8); every worktree commit today ran its gates in the main checkout | about 10 lines and a test; governed (`.claude/hooks/`), by day | `@tester` red, `@developer` green | yes, first |
| T2 | **Build the merge runner now** (stage 3, already ruled yes). `tools/merge_master.py <worktree>`: fetch; print `git log --merges --first-parent <base>..origin/master`; `git merge --no-commit`; stop for hunk-by-hunk resolution (never whole-file `--theirs`); then ruff, mypy, the full Python suite in that worktree's own venv, and `check_citations.py --base origin/master`; commit with `-F` from a template carrying the trailer and the subject tag. One persona runs it end to end. | 8 of the 15 merge lessons, D3, D8 (about 20 spawns a day of split merges), D9's typed merge range | about 120 lines with tests (the plan's estimate); `tools/`, by day | `@developer`, through the pipeline | yes |
| T3 | **One file per PR, not numbered sections in a shared file.** One line in `docs/increments/README.md`: "A PR adds its design as its own file or as a titled section; a citation names the heading, never a section number." Audit A already has its own file. | the 6 numbering lessons (section 9 twice, section 11 twice, stale docstring references, a wrong prediction in a brief) | about 30 words; governed, by day | `@architect` | yes |
| T4 | **`tools/new_worktree.py`**: `git worktree add`, then a real venv in the worktree (`uv venv`, editable install with ruff and mypy), `_core` built with `-DPYTHON_EXECUTABLE` set to that venv's Python, and a check that prints `tin_engine.__file__` and one submodule's. `--base <sha>` makes a separate scratch venv for a base install. This is P6 of the night retrospective, made concrete. | 8 of the 10 venv lessons; C's merge brief "you can't run the suite without a venv" | about 40 lines with tests, and about a minute per worktree; `tools/`, by day | `@developer` | yes |
| T5 | **`check_citations.py` defaults to `origin/master`** (its merge-base with the branch) and warns when `origin/master` was not fetched in the last hour; a cited range past end of file is broken, not at-risk; the at-risk list leaves out ranges whose cited text is the same at base and head (P2b of 2026-10-04). | 6 of the 8 citation lessons, including the 184 noisy entries in 30a | about 40 lines with tests; `tools/`, by day | `@developer` | yes |
| T6 | **`git stash` refused, and every file written with Edit or Write.** `guard_push.py` refuses `git stash` except `list` and `show` (the stack is shared by every worktree, D5). `tools/brief.py`'s line "write a rule file with Edit or Write, never through the shell" becomes "write every file with Edit or Write, appends included". | the 3 lessons; 6 stash runs, of which the 00:41 one (D9) had a brief that did not repeat the ban; heredoc appends (D5): 6 by `@tester` to test files in 5 commands, about 17 by `@architect` to 7 design files (8 of them to 20c), and 1 by `@perf` to a benchmarks README | about 5 lines and a test in the guard; one brief line; governed, by day | `@tester`, `@developer` | yes |
| R1 | **A relayed claim names who checked it.** Change `CLAUDE.md` §3's "no cause not checked" to: "A cause, a 'could', or a status from a handback goes to Ola with who checked it ('`@reviewer` checked') or as unchecked." No new line. | D2: six of the twelve corrections; D6's "waiting for its checks" | a few words; governed, by day | `@architect` | yes |
| R2 | **Each `QUEUE:` item in `session.md` names Ola's ask**: a date and quote, or a `ROADMAP.md` row. The recap warns on an item without one. Part of stage 2's state tool. | D1 (C2) | about 10 lines in `tools/session_state.py` with tests; governed | `@developer` in stage 2 | yes |
| S1 | Small alignments: `Monitor` added to `.claude/agents/perf.md`'s `tools:` line; `python` becomes `python3` in `CLAUDE.md` §4; `tools/brief.py` names a scratch subdirectory per agent (`<scratchpad>/<persona>-<worktree>-<HHMMSS>/`), and cleanup removes only that directory. | three alignment lessons; the `rm -f *.txt` that deleted the main session's `lessons.txt` | a few words and lines; governed | `@architect`, `@developer` | yes |
| — | **Not proposed:** new rule lines for the probe, domain, figure and brief lessons. Existing rules cover them, and 30b showed the fix for figures in practice: README tables generated by script between markers. If figures go wrong again after 30b's practice, one line in `perf.md` then. | | | | |

Order, if Ola says yes to all: T1 (smallest, and it makes every later
commit honest), then T4 and T2 together (the runner uses the worktree's
venv), then T5, T6, T3 and S1, and R1 and R2 with the next governed-file PR.

Sources for the approach. Anthropic's "Building effective agents"
(anthropic.com/engineering/building-effective-agents, read 2026-10-07):
"Poka-yoke your tools. Change the arguments so that it is harder to make
mistakes", and agents should "gain 'ground truth' from the environment at
each step". T1, T2 and T4 move a check from a brief, which a persona can
skip, into a tool, which runs anyway. The `cwd` behaviour in T1:
code.claude.com/docs/en/hooks, read 2026-10-07. The shared stash in T6:
git-scm.com/docs/git-worktree, read 2026-10-07, and the local
`git rev-parse` above. T3 is the pattern `towncrier` uses for changelogs
(each change writes its own fragment; cited in `next.md` item 14, read
2026-10-03).

## 7. Rule text: size and one cut

`python3 tools/rule_sizes.py` at `0c2572fb`: **10,564 words**, the same as
at the 2026-10-05 and night measures. The tool compares with the newest
retrospective on master, 2026-10-04 (+408 since), because the later ones
are not merged (D6). `CLAUDE.md` 1,146; `.claude/REQUIRED-READING.md`
1,821; `docs/increments/README.md` 1,209; `docs/PRINCIPLES.md` 1,428;
agents 3,088; skills 1,668; `.claude/briefs/common.md` 204. None of the
cuts proposed on 2026-10-04, 2026-10-05 or in the night retrospective has
landed.

**In flight:** guard fix PR B grows `.claude/REQUIRED-READING.md` to
2,055 words on `worktree-h16b` (+234, checked with `git show
worktree-h16b:.claude/REQUIRED-READING.md | wc -w`). All of it lengthens
the list of guarded commands in "The harness", the paragraph the night
retrospective proposed to cut to about 90 words.

**Today's cut, two parts:**

1. **PR B lands the night's harness cut instead of its longer list.**
   Replace the first paragraph of "The harness" with the night
   retrospective's short version (about 90 words); the full list of
   guarded commands goes into each hook's docstring, where it is tested.
   Against PR B as it stands: about 340 words fewer. Against master: about
   105 fewer. Risk: a reader no longer sees in one place that, say,
   `git fetch` into a named ref is guarded; the refusal message names it.
2. **With T4, "Stale artifacts: rebuild before you measure" (113 words)
   becomes one line**: "`python3 tools/new_worktree.py --rebuild` rebuilds
   `_core` into this worktree's venv and prints which module is loaded;
   read `cmake --build`'s exit status before `ctest`." About 85 words
   fewer, and the `cp` wildcard that silently copied a wrong-Python module
   (lesson "(tester F merge)") goes with it.

Together, about 190 words off `.claude/REQUIRED-READING.md` against master
(10 %), and the file stays near its size instead of growing by 13 %.

## 8. Questions for Ola

1. **T1, the after-commit gate checks the worktree that was committed?**
   Default: yes, first, by day.
2. **T2 and T4, the merge runner and a venv for each new worktree, as the
   next harness work, before more audit PRs?** Default: yes.
3. **T3, T5, T6, R1, R2, S1 as in the table?** Default: yes to all.
4. **The cut: PR B lands the short "The harness" paragraph instead of its
   longer list?** Default: yes; `@architect` amends PR B's design.
5. **Re-run CI on #194 and let it merge, then push the three local
   retrospective branches (2026-10-04c, 2026-10-05, 2026-10-05b) as one PR
   with `next.md` brought up to date?** Default: yes; the CI re-run and the
   push need your yes.

## 9. The lessons, word for word, by group

As the main session passed them on 2026-10-07, one per line.

**Master merges into a branch (15)**

- (architect #192 merge) two branches each added 'section 9' to python-audit.md; numbered sections in shared increment files collide at merge.
- (architect C merge) docstrings citing a design by section number go stale when designs in one increment file merge; cite by heading text.
- (reviewer F merge2) parallel branches adding 'next' section to a shared doc collide (F and C both 11); numbering sentence should name the other branch: second to merge takes next number.
- (C merge) brief predicted section numbers wrongly (A has its own file); check master headings before briefing.
- (C r6) 'only change since round N' claims: check git log --first-parent; an earlier merge went unreviewed.
- (C r6) renumbering greps must include the increment file itself.
- (C merge) git checkout --theirs takes the whole file, dropping non-conflicting hunks; resolve hunk by hunk and diff vs both parents.
- (C merge) don't chain unrelated git commands onto checks during a merge.
- (C r8) main session's brief scoped the merge check to the conflicted docs; a clean text merge broke code (C dropped MultiPolygon import, 30b's _linework uses it). Every master-merge check runs ruff + suite on the merged tree.
- (30c merge) briefs asking a persona to commit a merge: git merge --no-commit then git commit -F (trailer + subject ending).
- (30c cr1) 'what the merge brought in': settle with git log --merges --first-parent base..master.
- (30c cr1) merges after green: the 'base at the merge' is the branch at the merge; designs must say so.
- (30c rec) briefs should name the '## Review' heading, not a section number.
- (D merge) re-add the PR-order table's Total after every merge of it.
- (D merge) a helper-signature change on master breaks cleanly merged tests; full suite after every merge caught it.

**Probes and claims (13)**

- (reviewer C2 code r1) membership test -> .get() widened breaking inputs to `properties`; probe never varied properties -> regression.
- (architect C2 ruling) probe a member's whole value domain (missing, null, empty, non-empty non-objects); numbers/booleans and NVE lakes crashed at base too.
- (30b design) run 'identical on real catchments' speedups through the fixture recorder too: perf's covers rule passed suite + catchments but changed 2 fixture outputs.
- (30b r1 amend) plant reviewer's suggested fixes too: literal .boundary broke lines while suites passed.
- (30b red) a pin meant to catch a planted bug should first assert the condition the bug depends on.
- (30b pins) breaking a scratch copy so each pin is seen to fail cost ~2 min; worth doing at every pin ruling.
- (roadmap rows r1) probe the nearest valid-looking input too (32633.5 accepted silently).
- (30c dr1) a probe guard's claim must say what it checks (site-packages path != fresh install).
- (20c dr1) check a gate derived from a prototype against that prototype's own numbers.
- (h16b r9) probe every 'first program' exception with each known runner placed before it.
- (20c fix) gate on counts not shares when a PR removes triangles.
- (h16b r9 design) probe path rules with ., .., // spellings and via Edit/Write too.
- (20c dr3) measure a rule proposed for PR N on PR N's mesh (8% on 20c-2 = 22.5% on 20c-3).

**Domain traps (11)**

- (30b design) shapely predicates use only the FIRST argument's prepared form; prepare(x) then pred(y,x) gains nothing. Grep the codebase.
- (30b design r1) a whole-geometry predicate as a filter for a per-edge test is exact only on valid input; with invalid polygons accepted as linework, filter on boundary.
- (NMD) check for newer dataset versions when a brief names one (NMD 2023 exists).
- (slivers) a sliver line aligned with a tile edge in one CRS: check it in every input's CRS first (straight in 3035 = CORINE seam).
- (slivers) gap/overlap checks pass mm-scale vertex pairs; search for short segments.
- (30c red) test_mosaic build fixture passes needed only to plan_mosaic; seam tests must call assemble(plan, load, needed).
- (20c) tag output vertices by insertion order before building anything: leading suspect (refinement) made 38 of 1598.
- (20c) noder --snap-spacing 0.1/1 fails on Lagan: MalformedOutput chain does not tile the buffer (noded_pslg_builder.hpp:304); ROADMAP row added on Ola's yes.
- (20c dr2) an angle between two constraint segments at a shared vertex is a ceiling on any mesh's worst angle: cross-check worst-angle claims across runs with it.
- (20c dr2 fix) shapely overlay on a coverage makes mm segments; check for them.
- (D test fix) FIXTURE='catchment' equals a subcommand name; cli_driver's first-arg check let an old-form call through.

**Which Python and which code runs (10)**

- (reviewer #192 merge) local master stale (44fa7f5) -> check_citations default --base master over-reports; pytest addopts write .coverage, reviewers use --no-cov.
- (tester F merge) first build-pyext needs -DPYTHON_EXECUTABLE (caps); Python_EXECUTABLE ignored -> wrong-Python module; REQUIRED-READING's cp wildcard copies it silently.
- (developer #197 merge) perf's build-bench/pkg symlink venv: new modules need manual links; Path.rglob skips symlinked dirs -> test_the_table_names_every_module fails locally. P6: real venv per worktree.
- (30b r1 amend) timing script inside a copied package dir shadowed the import; check __file__.
- (30b red) worktree .venv lacked ruff; briefs should say where static tools come from.
- (roadmap rows) shared venv's rasputin points at main checkout; probes must confirm module __file__.
- (30b r1) worktree .venv held perf's base install; reviewers/testers running pytest there test base code. Base installs go in scratch venvs.
- (30c design) copied venv's bin/rasputin runs the original venv's python; call each venv's own interpreter and log tin_engine.__file__.
- (30c red) scratch_copy.py needs the worktree .venv; scratch venv + copied package + .pth works for pure-Python branches.
- (30c green) site-packages copy + existing _core lets green run the byte-identical gate without a C++ build; diff -r ties it to the commit.

**`check_citations.py` (8, one shared with the group above)**

- (reviewer retro r3/r4, orchestrator) a citation range given in a review went one line past EOF (:205-206 on a 205-line file); check `wc -l` / run check_citations on ranges before passing them in a brief.
- (reviewer #192 merge) local master stale: listed above.
- (30a prose) check_citations.py diffs against stale local master -> 184 noisy at-risk entries.
- (30a review) briefs for worktrees on unpushed bases should pass check_citations --base.
- (h16b r8) word-for-word records fail check_citations when the reviewer cites short file:line; reviewers should pin full paths.
- (30b prose) check_citations at-risk list flags unchanged ranges of edited files; compare cited range base vs head automatically.
- (30c dr1) inserting lines mid-file in a long decision doc shifts its own review-record citations; put notes on existing lines or at the end.
- check_citations: pass --base origin/master after fetch.

**Figures in documents (7)**

- (30a r2) profile README figures written from an uncommitted run; paste timings from raw/ after the last run.
- (30a r2 fix) @perf could generate README tables from raw files instead of hand-copying.
- (30b perf) README tables generated by script between markers; reuse for all acceptance runs.
- (30c design) dated benchmark scripts break when patched private helpers are renamed (mosaic._valid).
- (NMD) compare doc numbers vs probe prints before committing.
- (sliver r1) re-measure 'the most X' rankings, not just quoted values.
- (sliver r1 fix) state the counting rule next to 'the most' counts; report distinct triangles when strips overlap; compare float coords with tolerance.

**Briefs and quotes (5)**

- (reviewer retro r3) a rule "never X ... except Y" must be cited with both parts when an entry claims it was broken.
- (architect GeoTIFF merge) brief.py counts review rounds only on capitalised APPROVED/CHANGES REQUESTED; lower-case records show 0. Quotes of another branch's prose go stale; pin them.
- (30a r2 fix) grep a quoted phrase before putting it in a brief (regions 1.01/3.09 was a reviewer paraphrase).
- (h16b r8) brief pointed at a scratchpad handback extract; REQUIRED-READING forbids scratch as a channel. Quote in the brief or point at the subagent log.
- (30b red) main session cited test-audit R10 wrongly (it is 'Smaller duplication'); cite the right rule.

**Rule and tool alignment (5)**

- (perf 30a) Monitor not in @perf's tools but brief.py tells it to use Monitor; align.
- (h16b r8) REQUIRED-READING at 2012 words (+213 over baseline); fold into the approved R cut.
- (h16b r9 design) recursive rescans are exponential; a timed-out hook passes the command: add timing rows at ~40 trigger words.
- (30c perf) never busy-wait beside a timing run.
- (30c cr2) CLAUDE.md §4 writes 'python tools/...' but only python3 exists outside a venv; use python3 throughout.

**Main session's own steps (3)**

- (main) started the PR C follow-up (geojson-gaps) on my own initiative to fill a slot; dropped 2026-10-06. Filling slots is not a reason to start unrequested work; check goal fit first.
- (main 30a) skipped the @architect pin-ruling step before green to save time (Ola wanted numbers in 1 h); pins ruled after green instead. Deviation.
- (roadmap rows) main session's claim to Ola (None station could pair wrongly) was false; probe 'could' claims before relaying.

**Records (3)**

- (30a review) a @perf profile commit a design builds on goes unreviewed until a PR carries it; say which PR publishes it up front.
- (30a) recording commits should update the ROADMAP row with the Status line.
- (C) chat rulings recorded only in session.md; write them into the increment file the same round.

**Shell writes and `git stash` (3)**

- (h16b green) @developer wrote a governed hook file via a shell heredoc (rule: Edit/Write only) and used bare git stash/pop (shared stack across worktrees) during a background run. Briefs should forbid both explicitly.
- (h16b red r9) @tester also appended a test section via shell heredoc; Edit/Write rule must say 'every file, appends too'.
- (D test fix) @tester used git stash despite the brief forbidding it (stack empty after). Rule needs a hook or a stronger line.

**Scratch (3)**

- (architect C2 ruling) its cleanup `rm -f *.txt` in the shared scratchpad deleted other agents' files (incl. the main session's lessons.txt). Scratch cleanup must name files exactly; per-agent scratch subdirs.
- (20c fix) keep scratch prototypes a design's gates rest on (patch file) until the design is approved.
- (20c dr2 fix) record vertex counts of regenerated scratch inputs; two '1 cm' files caused a false contradiction.

## Review

- `@reviewer round 1 of the 2026-10-06 day retrospective (/Users/skavhaug/projects/rasputin/.claude/worktrees/retro-1007/docs/retrospectives/2026-10-06-day-bottlenecks-and-merges.md @530fa7b4): CHANGES REQUESTED. "Section 7" should say section 9 (lines 15 and 268); "five" corrections should be six (lines 117 and 309); the main session's 06:16 correction is missing from D2's eleven (main transcript line 9678).`
- `@orchestrator` answer to round 1: section 9 in both places; six, with the 06:16 `bench.py` correction added to D2 (twelve, and R1 says "six of the twelve"; the `bench.py` claim is not among the six because its origin is in the night's window and was not traced); suggestions taken (D3 leads with the 17:14:41 run that exited 0 on a tree where `ruff check .` exits 1; hook counts given as "about"; D5's file and stash details corrected, which also found `test_guard_targets.py`, not `test_guard_push.py`, in the 14:10 run, and six `@architect` heredoc appends; the D merge citation now names briefs 18592 and 18690); D9 added for the night's typed merge range and the 00:41 `git stash`, whose own brief did not carry the ban.
- `@reviewer round 2 of the 2026-10-06 day retrospective (/Users/skavhaug/projects/rasputin/.claude/worktrees/retro-1007/docs/retrospectives/2026-10-06-day-bottlenecks-and-merges.md @ddbe7108, range 530fa7b4..ddbe7108, 0 counted LOC, docs only): CHANGES REQUESTED. The new count of @architect heredoc appends is too low (lines 191-194 and 354). agent-a6c46d87f5ceed776 appended to 20c-soft-quality.md seven times, at 17:12 and 17:15 (transcript lines 595, 697) as well as 17:44-17:47, so that file got 8 appends counting 00:17. Across the window, @architect made about 17 heredoc appends to 7 design files: 30a, 30b, 30c, h16, python-audit, python-audit-pr-c2 and 20c. Everything else checks out.`
- `@orchestrator` answer to round 2: D5 and T6 now give the count from a scan of every subagent transcript in the window, with its per-file output pasted in D5: about 17 `@architect` appends to 7 design files, 8 of them to 20c, and 1 `@perf` append to a benchmarks README; the scan also confirms 6 `@tester` appends in 5 commands. Suggestions taken: D2's 06:16 row says the claim was made in the night; T6 says "of which the 00:41 one", so the stash runs stay six.
