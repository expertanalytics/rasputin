# 2026-10-17 overdue check: the night of 16 to 17 October, the merges #269 to #278, a code audit for machine-written oddities, and why my checks stopped

Recorded by `@orchestrator` on 2026-10-17 at about 09:30 UTC, on
`worktree-retro-1017` off master `e5b9b044`. **Times are UTC**; Oslo is
UTC+2. Every rule or harness change below is a proposal; Ola decides.

Ola, 2026-10-17 08:47: "I havent really heard from orchestator lately. Also,
reviewer. Are developing proper code, without AI-spesific oddities, and are
we keeping in control by checks and balances?" Later, through the main
session: why did `@orchestrator` never propose a check on its own cadence?

Sources: `<git-common-dir>/harness/windows.jsonl`, `queue-2026-10-16.jsonl`,
`queue-2026-10-17.jsonl`, `grants-used.jsonl` and `reviews.jsonl` (read
only); `.claude/current-task/session.md` as of about 09:05; the main session's
transcript `~/.claude/projects/-Users-skavhaug-projects-rasputin/dd80fd6d-….jsonl`,
cited as `T:<line>`; `TZ=UTC git log --all` for the window; `gh pr list`;
`python3 tools/count_loc.py`, `tools/check_code_history.py`,
`tools/rule_sizes.py`; `git diff -U0 7d1338ed e5b9b044` (all code merged
since 2026-10-14) for the audit.

Words used here. **Away window**: the time Ola is away and `tools/away.py`
has unattended mode on. **Grant**: a branch Ola named with `away.py --grant`,
which may be pushed and merged while he is away. **Red step**: `@tester`
commits failing tests; **green step**: `@developer` makes them pass.
**Recording commit**: the main session's commit that writes a review's
verdict into the increment file. **Governed file**: a file the governance
guard asks about (rule files, hooks, `tools/`). Increment numbers: 57 vector
datasets in the source catalogue (#275), 58 "class map" renamed "legend"
(#277), 60 `DemInput.label` deleted (#276), 61 lake-inflow routing
(private), 62 observed melt-out (private, parked), 63 `rasputin fetch` for
lakes and land cover, 64 `plugin_api` exports `accumulate` (#279), 65 the
refined mesh is Delaunay at its written coordinates (#278), 66 no two
vertices written at one point, h22 harness fixes (private push grants,
after-commit false reds, scratch copies of the private package), 23c-2 the
basin-pieces branch.

## 0. Why the checks stopped (the finding about my own design)

**What happened.** The last retrospective on master is #268, merged
2026-10-15 05:39. The next one, `worktree-retro-1015` (`b3e05434`, the night
of increment 50), was committed at 2026-10-14 23:06 and never pushed; it is
carried on this branch unchanged. Since then: 10 pull requests merged
(#269 to #278), plus private #10; three away windows ended
(`windows.jsonl`: back 2026-10-16 03:51, 22:35, 2026-10-17 06:45); no
`@orchestrator` run until this one. The main session told Ola at 08:48
(`T:27115`): "That's on me."

**Not only this week.** Five retrospective branches never reached master:
`worktree-retro-1004c` (3 files), `-1005` (the cut list), `-1005b`,
`-1008` (with `lead_time.py`) and `-1015`
(`git rev-list --count origin/master..<branch> -- docs/retrospectives`).
Nothing reports them.

**The design fault.** Every trigger for my runs is the main session's choice
(`CLAUDE.md` §3, *When to spawn `@orchestrator`*), and so is pushing what I
write. The watcher's schedule is kept by the one it watches. When the main
session is busy, or behind on context, the check that would say so is the
one that is skipped. No research round has run since about 2026-10-12,
though the rule says weekly.

**Why my retrospectives missed it.**

1. Each run is scoped by its trigger: one night, one increment. It reads that
   window's commits, queue and transcript. The gap between two runs lies
   outside every window I am given.
2. A run that does not happen leaves nothing to read. My duty list
   (`.claude/agents/orchestrator.md` §1) flags acts that left a trace: a step
   skipped, work out of role, a guard worked around, idle, a false report. A
   missing duty is none of these.
3. `tools/rule_sizes.py` measures against the last retrospective **on
   master**, so an unmerged one is invisible to it, and so to me.
4. My proposals so far were guards on what agents do, never on what the
   dispatcher fails to do. `grep -ril 'overdue\|duties' docs/retrospectives`
   finds no earlier mention.

**What the literature says** (each source read today):

- Claude Code's hooks guide: hooks give "deterministic control: certain
  actions always happen rather than relying on the LLM to choose to run
  them." `UserPromptSubmit` "always fires on every occurrence", also "on
  turns Claude Code starts on its own", and its stdout or
  `additionalContext` goes into Claude's context.
  ([hooks guide](https://code.claude.com/docs/en/hooks-guide))
- The Prometheus operator's `Watchdog` alert: "This is an alert meant to
  ensure that the entire alerting pipeline is functional. … If not firing
  then it should alert external systems that this alerting system is no
  longer working." Silence is the signal, and something outside the watched
  system must notice the silence.
  ([runbook](https://runbooks.prometheus-operator.dev/runbooks/general/watchdog/))
- MAST (Cemri et al. 2025, arXiv 2503.13657) groups multi-agent failures
  into "system design issues", "inter-agent misalignment" and "task
  verification". A verifier that runs only when the verified agent calls it
  is a design issue, not a verification one. (The abstract page only; the
  14 modes are not named there.)
- Anthropic, *Managed agents*: "Harnesses encode assumptions that go stale as
  models improve." The parts should be interfaces that "could fail or be
  replaced independently".
  ([post](https://anthropic.com/engineering/managed-agents))

All four say the same thing: a check that only runs when the checked process
calls it fails silently, exactly when that process is overloaded. The fix is
a trigger that neither the dispatcher nor the watcher owns, with silence
reported out loud. See P1.

## 1. The night of 16 to 17 October

Away mode on at 23:24:03, `--back` at 06:45:45 (granted until 11:24;
grants: `worktree-legend`, `worktree-dem-label`). Ola said "Suspend now
until I get home" at 22:44 (`T:20440`), then asked at 23:22 "Do you have
enough to do over the night?" and gave the grant command
(`T:20804`-`T:20838`). The suspension was lifted by that exchange, but its
`SUSPENDED` line is still in `session.md`.

### 1.1 Idle time

57 commits in the window, across all branches. Gaps of 20 minutes or more:

| Gap | Length | What ran (handback that ended it) |
|---|---|---|
| 06:07 to 06:41 | 34 min | `@developer`, 66 green (06:41:39, `T:25986`) |
| 04:59 to 05:24 | 25 min | `@developer` h22 (refused, 05:03), then `@architect` 66 design (05:24) |
| 05:43 to 06:07 | 24 min | `@tester`, 66 red (06:08) |
| 00:28 to 00:52 | 24 min | `@tester`, 63 red (00:54) |
| 03:57 to 04:21 | 24 min | `@architect`, h22 design (04:21) |
| 01:06 to 01:29 | 23 min | `@developer`, 64 green (01:28) |
| 23:57 to 00:19 | 22 min | `@architect`, 63 design (00:19) |

**No idle.** An agent was at work through every gap; the window ended early
because Ola came back. This is the second night in a row with none.

### 1.2 The nine queued refusals

`queue-2026-10-17.jsonl` holds nine entries, all in the window.

| # | Time | Who | Act | Verdict |
|---|---|---|---|---|
| 1 | 23:25:54 | main session | `gh pr create --head worktree-legend … --body "…"` | **False positive by meaning.** The grant covered `worktree-legend`, but `tools/push_grant.py@df923398:42-44` accepts `--body-file`, not `--body`. The PR waited 7 h 31 min, to 06:56 (#277). Told to Ola at 23:27 (`T:21096`); no `GUARD FALSE POSITIVE` line. |
| 2 | 23:49:13 | main session | `git pull --ff-only origin master` in the main checkout, to reinstall `.venv` | Correct by the rule (writes a ref). Not retried: `git reflog master` shows the next move at 07:01, attended. |
| 3 | 01:38:32 | `@developer` (63) | a heredoc that runs `tools/count_loc.py` against a temporary tree commit | **False positive**: it reads the tool; the guard saw the path in text it could not parse. |
| 4 | 03:42:00 | main session | a `bash -c` loop of `git merge-tree` over branch pairs | **False positive** (a permission prompt, refused while away). Recorded. The check `CLAUDE.md` §3 *Conflicts before a push question* asks for could not run as one line; one plain `git merge-tree` per pair would have. |
| 5 | 04:24:58 | `@reviewer` (h22) | `git push private 0123…` fed to `guard_push.py` as a probe | Not a false positive: a probe that wrote a fake request into Ola's queue. Recorded, with the lesson. |
| 6 | 04:27:46 | `@reviewer` (h22) | `python3 -c` reading `.claude/settings.json` | **False positive** (a read). Recorded. |
| 7, 8 | 04:39:14 | `@architect` (h22) | `git config` in a throwaway repository in the scratchpad, after a `cd` into the worktree | **False positive by meaning**, one command, both guards. The scratchpad exemption needs a single plain command. Recorded once. |
| 9 | 05:02:46 | `@developer` (h22) | a Python script rewriting `tools/push_grant.py` and `tools/away.py` | **Correct**: governed files while away. Not redone in the window (files changed 07:53, after Ola's yes at 06:55). Two faults upstream: the main session queued a green step whose files are all governed, and the persona edited them by script, not with Edit. |

Five of the nine commands (six entries) were false positives, two were
correct refusals, one was a probe. Today two more of mine:
`away.py --back` (the brief asked for it, but `.claude/REQUIRED-READING.md`
*Unattended mode* says no agent runs `away.py`: the guard was right and the
brief wrong), and a read-only `python3 -c` over `reviews.jsonl`, refused as
a harness write (false positive; read with Read instead).

### 1.3 Work out of role, and other deviations

1. **The morning report misstated a fact.** "Nothing was pushed: none of
   the pushes were covered by the grant" (06:46, `T:26079`). In fact
   `grants-used.jsonl` shows `worktree-legend` pushed at 23:25:45, and
   `worktree-dem-label` pushed, opened and enqueued at 23:26. It merged as
   #276 at 23:48.
2. **"Yes, we're within the rules"** (08:36, `T:26889`), then at 08:48 "four
   items are in flight (63, 64, 66, h22) against the cap of three"
   (`T:27115`). By the rule's text the count was higher: 23c-2 also had
   spawns and no merge, and was not parked. Items blocked behind a PR that
   waits on Ola are ambiguous under the rule as written.
3. **Session note out of form.** `.claude/REQUIRED-READING.md` allows one
   `NOW:`, one `QUEUE:` and `ASK OLA:` lines, "no rulings, no history".
   `session.md` has 69 lines: 4 `NOW`, 5 `QUEUE`, 17 `LESSON`, and rulings,
   findings and decisions besides (13.6 kB). The recap is 19,030 characters
   against its 10,000 cap (`python3 tools/session_state.py | wc -m`). The
   lessons piled up there because I was not run to take them.
4. **Stale note files.** `perf-refine-world-053047.md` and
   `reviewer-fetch-vectors-050656.md` are still in `.claude/current-task/`.
   Both personas handed back before 04:00. The spawner deletes the file on
   reading the handback.
5. **The main session edited design files outside a recording commit**:
   `e23131f1` (61's flow-direction sentence) and `aa171207` (51 P2's
   `project_structure.md` and ROADMAP record). Both belong to `@architect`.
   Small, but the same drift as before.
6. **A merge brief forbade a C++ build** that master's 65 change needed
   (07:49, `T:26655`). Respawned at once, with building allowed.
7. **The merge watcher misread the queue.** Its loop exited at 07:09 with
   "277 left the queue (CLEAN)" (`T:26471`) because `mergeQueueEntry` was
   empty. #277 merged at 07:19 and #278 at 07:23. The main session found out
   at 07:47, when `gh pr merge 277` answered "already merged" (`T:26526`).
8. **Leftover `tail -f` waiters**: 2 stopped at 03:57 (`T:24624`) and 8 at
   08:36 (`T:26889`), between 5 hours and 2 days old. The cause is in the
   files: `.claude/briefs/common.md:12` tells every persona to "Wait for a
   background run with the Monitor tool", and no persona's `tools:` line
   names Monitor (`grep -n '^tools:' .claude/agents/*.md`). The private
   ideas file has had the fix since h19: *Monitor in every persona's tool
   list*.

### 1.4 Lessons passed on, recorded here

From `session.md`, as the personas and the main session wrote them:

- A worktree needs its own venv for tests that spawn subprocesses. The
  shared `.venv` is an editable install of master: 14 false "Unknown tool"
  failures.
- #273 failed on Python 3.12 and 3.13 only: `Path.exists` raises
  `ENAMETOOLONG` there, not on 3.14. Judge pytest by its exit code: #272
  failed CI on the prose-read check with all tests passing.
- `tools/scratch_copy.py` breaks for rasputin-lsm (now h22).
  `push_grant.py` had no private-remote form (now h22).
- Agents must not `git stash` or check out in a worktree while its suite
  runs; the stash stack is shared.
- `merge_master.py --continue` runs the full suite (over 10 min): brief it
  as a background run. A `@developer` started three such runs on one
  worktree. Idea: a lock in `merge_master.py`.
- Merge-tree's conflict list misses a file that merges clean but keeps a
  renamed name. A rename's merge brief should let `@tester` make mechanical
  renames in the branch's own tests.
- Before a merge brief that forbids a build, check
  `git diff --stat <base> origin/master -- bindings src include`.
- Mutation rounds: rebuild each mutant and confirm the build. A mutant that
  removes a variable's only use fails to compile under `-Werror` and can
  pass for a kill. A survivor in `ctest` may still be killed by pytest.
- A header edit moves line-number citations below it: grep
  `'<file>:[0-9]'` before committing.
- `pytest -m harness` alone always fails the 85 % coverage floor: read the
  pass count, or pass `--no-cov`.
- Worktree venvs have no pip: rebuild the extension with
  `uv pip install --reinstall`.
- `bench_quick` prints no times without a stored baseline for the exact OS
  version. Idea: print the medians anyway.
- Probing `guard_push.py` with a push command while away queues a refusal:
  probe with a test fixture.
- A green step whose files are all governed cannot run while away. Check the
  governed list before queuing it for a night.
- Under the away grant: get an APPROVED verdict on the exact head, commit
  nothing after it, then push with the bare
  `git push origin <sha>:refs/heads/<branch>`.
- After each handback, check for leftover `tail -f` and stop finished
  agents.

## 2. The merges since 2026-10-15

`git log --merges --first-parent origin/master`, with `count_loc.py <merge>^1
<merge>` (net counted lines) and `reviews.jsonl`.

| PR | Increment | Net lines | Check |
|---|---|---|---|
| #269 | h21 PR A, rule text | 0 (prose) | Four design rounds, three sent back, then approved. Rule text by the main session at Ola's prompt. In order. |
| #270 | 51 P1, the open-core plugin seam | +157 | Design (two rounds), red, green, test run, review (two rounds, one CI-only). In order. |
| #271 | 51 P2, snow and radiation leave the core | −2,317 | Red (absence test), green, review. The docs commit `aa171207` is by the main session: @architect's file. |
| #272 | 54, the open core reads as one tool | −80 | In order. The CI fix `3cc2b848` (+6 lines in `tools/ci_changes.py`) came after the last review, with no recording commit. Small, but code after review is unreviewed. |
| #273 | 55, local datasets in the catalogue | +242 | Two review rounds, CI fix recorded. In order. |
| #274 | 56, `rasputin draw` deleted | −720 | Review sent it back once. In order. |
| #275 | 57, vector sources and the MCP output root | +340 | Two design–red–green rounds (Ola's mid-build ruling), three review rounds. In order. |
| #276 | 60, `DemInput.label` deleted | −3 | Pushed and merged overnight under the grant on the reviewed head plus a master merge. In order. |
| #277 | 58, legend rename | +2 | In order. The overnight PR refusal is in 1.2. |
| #278 | 65, Delaunay as written | +107 | Design (two rounds), red, green, mutation round, review (two rounds, a stale header comment caught), `@perf` acceptance (hashes unchanged, time within 5 %). In order; refine code, so `@perf` was required and ran. |

No PR came near the 700-line limit. **`@reviewer` is working:** 47 recorded
verdicts since 2026-10-15, 18 of them `CHANGES REQUESTED` (38 %). It caught
real faults: a header comment that no longer matched the code in 65, two
fail-open paths in h22's private push check, and a snap rule in 66 that can
move a point past its edge. Ola does not see this because the verdicts reach
him only through the main session's summaries (P7).

## 3. Code audit: machine-written oddities in code merged since 2026-10-14

Scope: `git diff -U0 7d1338ed e5b9b044`: 3,132 added lines outside `tests/`,
8,769 in `tests/`. Searched by pattern, then read: `source_report.py` (all
387 lines), `refine_points.hpp` (65's diff), `cli_kit.py`, `push_grant.py`
(parts), `record_review.py`, and `test_cli_help.py`, `test_legend_name.py`
and `test_plugin_api.py`.

### 3.1 What is fine

- **Comments say why, not what.** Of 91 added comment lines, a sample of 30
  had none that restates its line.
- **The mannered register is old, not new.** Phrases like "the whole
  point" or "blesses exactly the disagreement" occur 9 times in
  `src_python`, `include` and `tools`. All predate 2026-10-14 (e.g.
  `src_python/tin_engine/features.py@df923398:139-151`, from 2026-09-18); none were added since.
- **The new Python modules are clean under the history gate.**
  `source_report.py` and `toml_data.py` give no hits. Docstrings are 12 %
  of the documented functions' lines.
- **No defensive code that nothing reaches.** Three broad `except
  Exception` were added, each with a reason: a broken plugin is reported
  (`src_python/tin_engine/plugins.py@df923398:52`, `src_python/tin_engine/sources.py@df923398:326`), and a hook never blocks
  (`.claude/hooks/record_review.py@e5b9b044:43`). The `assert x is not None` lines narrow types
  for mypy.
- **Tests test behaviour.** One private import (`_core`, the extension).
  No test reads source by AST. Test names added since 2026-10-14 carry no
  design labels. The older suites hold 227 such names (`test_t5_…`) in 27
  files.

### 3.2 What to cut (counts, places, proposed cut)

| # | Pattern | Where | Count | Cut | A gate? |
|---|---|---|---|---|---|
| A | History in docstrings, **outside the gate's scope** | 14 modules not in `[tool.check_code_history].paths`: `feature_input.py` (42 hits), `stats.py` (16), `run_record.py` (13), `chains.py` (9), `palettes.py` (8, "since increment 41" at `src_python/tin_engine/palettes.py@61755c04:3-7`) … | 116 hits; 4 on lines touched since 2026-10-14 (`src_python/tin_engine/feature_input.py@61755c04:71`, `:566` "R5:", `src_python/tin_engine/legends.py@61755c04:3`, `src_python/tin_engine/palettes.py@61755c04:4`) | Make the gate cover all of `src_python/tin_engine` minus an exclusion list. A new module now escapes it by default. Sweep the 116 (comments only, no counted lines). | Yes: `check_code_history.py`, run with these paths, finds all 116. |
| B | History labels in **C++ comments** | `include/`, `bindings/`: "65:" on new fields at `include/terrain/refinement/refine_points.hpp@61755c04:133-135` and `:388`; "step 6" at `:28`, `:465`; pointers into `docs/increments/` in 29 files | about 156 comment lines (regex: `//.*(\b[0-9]{2}[a-z]?:\|\b[RNDQB][0-9]+\|step [0-9]\|docs/increments/)`) | Extend the gate to `//` comments in `.hpp` and `.cpp`, then sweep. | Yes, with the extension. |
| C | **A legacy alias kept and pinned by a test** | `src_python/tin_engine/cli.py@61755c04:66-68` imports `cli_kit` helpers under their old names (`_destination`, `_resources`, `_words`); `tests/python/test_plugin_api.py@61755c04:49-52` asserts the old names exist | 3 aliases, 1 parametrised test | Use the new names in `cli.py`; delete the test. Ola's rule of 2026-10-16: no aliases. | No lint rule sees a rename alias. |
| D | **Absence tests**: tests that pin that deleted things stay deleted | `tests/python/test_cli_help.py@df923398:224`, `:264`, `:271`, `:284`; `tests/python/test_open_core.py@61755c04:66`; `tests/python/test_mcp_vector_sources.py@61755c04:325`; `tests/python/test_reference.py@61755c04:546`; all of `test_legend_name.py` (169 lines: a repository scanner for old spellings, with 6 self-tests) | about 12 tests plus the 169-line file | Keep the ones that test behaviour (old options are refused: `tests/python/test_legend_name.py@61755c04:151-169`). Delete the scanners and the "not importable" checks, or give absence tests a lifetime. Code kept for code deleted. | None. A rule for `@tester`. |
| E | **A word built from a variable**, beside a helper that does it | `src_python/tin_engine/stats.py@61755c04:176-177` builds `"hole{plural}"`; `run_record.plural` exists and is exported | 1 | Use `plural(holes, "hole", "holes")`. | None. |
| F | **Duplicated helpers** | a `git`/`_git` runner in 9 files: `tools/check_code_history.py`, `check_citations.py`, `merge_master.py`, `session_state.py`, `brief.py`, `count_loc.py`, `rule_sizes.py`, `push_grant.py`, `.claude/hooks/gates_after_commit.py` | 9 | One helper in `tools/harness_mode.py`. Low priority: tools, not product. | None in the toolchain. |
| G | **Positional booleans and `None`** | `src_python/tin_engine/source_report.py@61755c04:257` `Facts("features", None, suffix, shapes, got.attributes, True, entry.legend)`, `:264`, `:329`; `src_python/tin_engine/source_report.py@61755c04:175` multiplies a tuple by a boolean | 15 hits across `src_python`, `tools`, hooks | Keyword arguments; an `if` for `:175`. | **Yes: ruff `FBT003`** (`ruff check --select FBT003` finds 15 today). |
| H | **A 78-line lambda with coded values** | `include/terrain/refinement/refine_points.hpp@61755c04:397-474` (`as_written`): names `wv`, `ctris`, `cbits`, `cmasks`; a `kind` vector of 0/1/2 explained in a comment | 1 | A named function in `detail::` with an enum for the three kinds. The clearest machine-look spot in the C++. | No function-size check in the C++ toolchain. |

Two patterns from the brief did not show: names that narrate history (none
in the new code, beyond the C++ "65:" labels in B), and long docstrings.
Only 9 functions have a docstring of 8 lines or more that is at least half
their code: `src_python/tin_engine/cli.py@61755c04:430` `_chain_masks`, `src_python/tin_engine/features.py@61755c04:110`, `:139`,
`src_python/tin_engine/feature_input.py@61755c04:115`, `src_python/tin_engine/io/models.py@61755c04:306`, `.claude/hooks/guard_push.py@e5b9b044:243`,
`src_python/tin_engine/cli_kit.py@61755c04:33`, `src_python/tin_engine/cli.py@61755c04:411`, `:388`. Most are older; trim them as files
are touched.

**Answer to Ola, short:** the new code is mostly sound. It explains why,
tests behaviour, and the reviewer sends back 38 %. The oddities are
systematic, not random: history labels where no gate looks (C++ comments,
14 modules outside the Python gate), a habit of keeping tests and aliases
for deleted code, and dense packing where the line limit pushes. A, B and G
can be gated; C, D, E and H need a reviewer's eye, so P6 adds them to
`@reviewer`'s checklist.

## 4. Rule text

`python3 tools/rule_sizes.py`: **11,550 words, +35** since
`2026-10-14-merge-46.md` (`CLAUDE.md` +55 for *Conflicts before a push
question*; `docs/increments/README.md` +19; `reviewer.md` −39). The −347
words carried from that retrospective (its P9) are still not applied.

**Proposed cut:** apply the carried −347. In addition:

- `.claude/REQUIRED-READING.md` *The harness*: the guard-by-guard account
  restates the hooks' code, including examples (`cp json.py tools`, "with or
  without the trailing `/`", `tools/./json.py`). Point to the hook files
  and keep the one rule each enforces: about −90 words. These are examples
  in a rule file, so they move to the hooks' docstrings.
- `.claude/briefs/common.md:12`, the Monitor sentence: replace it under P2.
  About 0 net words.

## 5. Proposals

- **P1. Duties computed by a hook, owned by neither the main session nor
  me** [proposal]. A pure read, `tools/duties.py` (about 60 lines with
  tests, governed), prints at most three lines, and only when something is
  overdue:
  - "Retrospective owed: last on master \<date\> (#\<n\>), \<k\> merges since";
  - "Night unchecked: away window ended \<time\>, no retrospective since"
    (from `windows.jsonl` `back` events against the retrospectives' commit
    times);
  - "Research round owed: last \<date\>" (weekly);
  - "Unmerged retrospective branches: \<names\>".

  Call it from `session_state.py` (the recap), from a `UserPromptSubmit`
  hook (into the main session's context on every turn, per the hooks
  guide), and, for Ola, from the status line, the one channel neither agent
  controls. Evidence: section 0, five unmerged retrospectives, three
  unchecked windows. Cost: one hook entry and one status-line entry in
  `.claude/settings.json` (Ola's yes); deletes nothing. The prose trigger
  list in `CLAUDE.md` §3 stays as the rule the hook checks.
- **P2. A waiter that exits.** Change `.claude/briefs/common.md:12` to say:
  "Run a long command with Bash `run_in_background` and wait for its exit
  notification; never `tail -f` or `sleep`." Alternatively add `Monitor` to
  the `tools:` lines (the ideas-file entry). Whether subagents get Monitor
  is unchecked. The recap also lists `tail -f` processes older than an
  hour. Evidence: 1.3 item 8.
- **P3. `brief.py` decides whether C++ may be built.** It runs `git diff
  --quiet <merge-base> origin/master -- bindings include src CMakeLists.txt`
  and writes "you may build C++" from the result; the main session no longer
  writes it. About 10 lines. Evidence: 1.3 item 6.
- **P4. A checked-in merge waiter.** It polls `gh pr view --json state`
  until MERGED or CLOSED, and reads an empty queue entry as unknown, never
  as "left". This is stage 3 of the dispatcher-control plan Ola approved on
  2026-10-04, not yet built. Evidence: 1.3 item 7.
- **P5. The grant accepts `gh pr create --body`** under the same no-expansion
  check as `--title`: 1 line and a test in `push_grant.py`. Evidence: 1.2
  row 1, a 7.5-hour wait.
- **P6. Code cuts A to H** (section 3), in one cleanup increment through the
  pipeline: A and B extend `check_code_history.py` first (tests first), G
  turns on `FBT003` in `pyproject.toml`, and C to H go to `@reviewer`'s
  checklist as "aliases, absence tests, words from variables, packed
  lambdas". It removes code; the 700-line limit counts the sweep as
  deletions or not at all (comments).
- **P7. Ola sees the reviewer's verdicts directly**: the recap prints the
  day's verdict counts from `reviews.jsonl` (sent back / approved, with the
  last three findings' one-liners). About 15 lines in `session_state.py`.
  Evidence: Ola's question.
- **P8. `session.md` held to its form.** The recap prints how many lines
  break the form, and its own length against the 10,000-character cap,
  first rather than as a warning at the end. Lessons go to the handback and
  then to me (P1 makes sure I am run). Evidence: 1.3 item 3.
