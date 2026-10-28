# Harness h36: cut the rule payload by half

Status: **designed**, 2026-10-28, round 3 (amended after its review, §7's sweep): origin/master `14a17091` merged
in, and every proposed text rebased onto today's version of the file it
replaces (§7). About 37 counted lines (§3), no C++ file, no refine or mesh
code, and no path `tools/check_lane.py` names as full-lane. Kept in the full
lane by the main session (design and review rounds already run). Most of the
PR is rule text, and six of the changed tools are governed
(`guard_governance.py` asks before each write), so both the main session's
prose step and `@developer`'s green step run with Ola present.

## Ola, verbatim

From the outside review he adopted on 2026-10-27, pasted with "we should do
this asap":

> Slash Rule Payload by 50%: Audit CLAUDE.md, REQUIRED-READING.md, and the
> persona briefs. Remove historical incident narratives (e.g., commit
> references explaining why Boost was banned). Retain only the active
> imperative command. Shift secondary guidance into dynamic skills that load
> strictly on-demand.

2026-10-28 13:33 UTC, answering "Should the rule cut also shorten the guide
to increments and the four domain skills, delete the principles file, and add
a small rebuild tool, beyond the three files the outside review named?":

> yes to the rule cut, all four files

Read here as a yes to all four parts of the wider scope: the increments
README, the four domain skills, deleting `docs/PRINCIPLES.md`, and
`tools/rebuild_core.py` (an open question in the round-3 handback, with that
reading as its default).

**Default the main session took while Ola was away** (not his ruling; he can
reverse it): the lean-process rulings, now merged (§5), are carried in
shortened; Ola sees each shorter wording at the final review.

The main session also reported that Ola likes tools that do one job and would
build harness tools where they save context (its report, not his words).

## Prior art: legacy and literature

*Literature.* Anthropic's Agent Skills documentation
(<https://platform.claude.com/docs/en/agents-and-tools/agent-skills/overview>,
read 2026-10-27) describes progressive disclosure: a skill's name and
description are always in context (about 100 tokens each), its body loads only
when triggered, and a script's code never enters context, only its output.
This design follows it: rules only some agents need become skills, and rules a
script can carry out become scripts. No novelty is claimed.

*Legacy.* `git grep -n -i -E 'CLAUDE.md|skill|required reading' legacy-archive -- legacy`
returns nothing; nothing is carried across.

*Ideas file.* Two entries overlap: "Morning check proposals P1-P12 and the
tester.md cut" (the cut is taken up here: `tester.md` 799 to 374 words) and
the `docs/PRINCIPLES.md` review (stale entries, D1 against two writers),
which the file's deletion settles.

*History in the counted files.* No incident narrative is left on master: the
example the review gives (a commit explaining a banned library) is not there
(`grep -n -i boost` over the counted files returns nothing). What history
remains is pointers (PRINCIPLES.md's *Origin* lines, perf.md's pointer to
`24-release-hardening.md`, REQUIRED-READING's pointer to its incidents file)
and label examples; all go, since `docs/retrospectives/` already holds the
incidents. The cut is mostly restatement, description of what tools already
do, and generic advice any model follows unprompted.

## 1. What happens to each file

Full new texts: `docs/increments/h36-proposed/` (file names there are not
the real names, so that no copy is loaded as a rule or guarded: Claude Code
loads a `CLAUDE.md` in any subdirectory it reads).

| File | Stays as an imperative | Moves to a skill | Replaced by a tool | Dropped |
|---|---|---|---|---|
| `CLAUDE.md` | §1 to §4, each cut to its commands; step order with both lanes | dispatch rules of §3, the first question included → `dispatch`; CI and merge queue → `publishing` | size-limit definition → `tools/count_loc.py` docstring; "merge commit, never squash" → the merge queue's method; `@orchestrator` timing → `tools/duties.py`; gate list → `gates_after_commit.py`, except `check_code_history.py`, which the hook does not run (`gates_after_commit.py` `gates()`), so its command stays in §4's block | method-flag and `-d` behaviour of `gh pr merge` |
| `REQUIRED-READING.md` | note file, skills, ideas, claims, approval boundary (what needs no fresh yes, and what needs one) and the at-risk citation re-read, harness hook list, data folders | cold start, `questions.py` and `work_pipes.py` → `dispatch`; grants, full-lane reviewer-before-push → `publishing`; unattended mode, fast-lane branches while away → `away-mode` | rebuild recipe → `tools/rebuild_core.py`; hook descriptions → the hooks' own refusal texts | recap-size rule; `touch` after a `cp` restore (mutants now run in scratch copies); "editing rule files is ordinary tree work" |
| `docs/increments/README.md` | the loop, the two lanes, acceptance pointer, cost constraints | verdict format → `reviewer.md`; acceptance details → `perf.md` | merge commit → the merge queue | |
| `docs/PRINCIPLES.md` | deleted; A6, B3, B4, D2, D3 become five lines in REQUIRED-READING's *Claims*; the six files that name it change (§3) | | | *Origin* and *Status* lines |
| `.claude/briefs/common.md` | all; tested phrases kept word for word | | | "never by a name pattern" |
| `architect.md` | §1 to §5 condensed | | | generic pillar prose |
| `developer.md` | all, condensed, the fast-lane paragraph included | | | |
| `orchestrator.md` | all, condensed | | | |
| `perf.md` | all, condensed; gains the acceptance run | | `--no-build` rebuild → `rebuild_core.py` | py-spy hint |
| `reviewer.md` | all, condensed; gains the verdict format | | | |
| `tester.md` | §1, §3A, §3C, §3D, §4 | | coverage floor → `--cov-fail-under` in the pytest config | Catch2 and tier prose |
| four domain skills | condensed, section numbers kept | | | marketing lines, generic advice |

Agent and skill `description` fields are cut to one sentence where long:
they are loaded into every session.

## 2. New skills and when each loads

| Skill | Holds | Loads |
|---|---|---|
| `dispatch` | recovery, spawning, work in flight, reports to Ola | by the main session before its first act (`CLAUDE.md` says so, and the recap names it after every start, `/clear` and compaction); never by a subagent |
| `publishing` | one yes per act, reviewer before push, conflict check, merge queue, CI | before asking Ola to push or merge, before any push, PR or merge, before calling a branch merge-ready |
| `away-mode` | unattended mode, refusals, grants, defaults while away | when the recap names it (unattended mode on), when Ola says he is leaving, when a guard refuses an act |

The approval boundary itself stays in REQUIRED-READING: every agent must know
the remote is Ola's without loading anything, and so does the writer's
re-read of the at-risk citations before a push. `publishing` says the review
comes first on a full-lane branch.

**Why the recap names the skills.** `CLAUDE.md` stays in context through a
compaction, but a skill body loaded earlier does not, and a line already
acted on is easy to pass over. The `SessionStart` hook runs
`tools/session_state.py` on startup, resume, `/clear` and compaction
(`.claude/settings.json`), so a line in its output reaches the main session
at exactly the moments the skill is lost. The line says "main session" so
that a subagent running the script reads it as not its own.

## 3. Tools

**New: `tools/rebuild_core.py`**, standard library, not governed. Rebuilds
the extension of the worktree that holds it and installs it where its venv
imports it.

- Refuses (exit 2, one `rebuild_core:` line) when `build-pyext/CMakeCache.txt`
  is missing, naming `tools/new_worktree.py` as what configures it.
- Runs `cmake --build <root>/build-pyext -j --target _core`; a non-zero exit
  is returned as its own, with one line, and nothing is copied.
- Asks `<root>/.venv/bin/python` for the `tin_engine._core` file it imports
  and for its extension suffix; copies `build-pyext/_core<suffix>` over it;
  exit 2 if either is missing.
- Prints one line: the destination and the first 12 hex digits of its
  SHA-256, equal to the built file's.

**Changed: `tools/count_loc.py`**, docstring only (governed): it gains the
file rule now in `CLAUDE.md` §2 (Python, C++, CMake and shell files outside
`tests/` and `docs/` count; a file renamed across that boundary counts in
full), so the tool is the one statement.

**Changed: `tools/session_state.py`** (governed). Right after
`== recap ==`, one line, before `Last landed:`:

- unattended flag on, or unreadable (the guards then refuse as when on):
  `Main session: load the dispatch and away-mode skills before your first act.`
- otherwise (off, ended, or state unknown):
  `Main session: load the dispatch skill before your first act.`

The mode is taken from the `read_mode` result `print_recap` already reads,
not from the header text. Both positions other tests pin hold: the header
stays above `== recap ==` (`test_session_state.py` `_header`, line 170) and
the duties line stays right after `Last landed:` (`test_duties.py`
`_after_last_landed`).

**Changed: the six files that name `docs/PRINCIPLES.md` or the recipe**
(found with `git grep -n PRINCIPLES -- tools .claude/hooks tests` and
`grep -n recipe tools/new_worktree.py`):

| File | Change | Counted lines |
|---|---|---|
| `.claude/hooks/guard_governance.py@b2207408:60` | drop the entry from `GOVERNED`; its docstring's line 14 mention stays (history, uncounted) | -1 |
| `tools/rule_sizes.py@14a17091:25` (governed) | drop it from `RULE_FILES` | -1 |
| `tools/ci_changes.py@14a17091:38` (governed) | drop it from `NOT_PROSE` and from the comment above | -1 |
| `tools/check_citations.py@b2207408:96` (governed) | drop it from `GOVERNED` | -1 |
| `tools/new_worktree.py@14a17091:11`, line 125 | the docstring names `tools/rebuild_core.py` instead of REQUIRED-READING's recipe; the printed `rebuild` line becomes `python3 <wt>/tools/rebuild_core.py   (configured, not built)` | 0 |
| `tests/python/test_session_state.py@14a17091:731` | stop copying the file | tests |

Four more test files pin the file name and change with it (§6).

**Changed: `docs/increments/h2-citation-pinning.md` line 109**: its example, a
section citation of A6 in the principles file, is cut, since the file goes (the only citation the move
breaks; checked by a script over every ID and heading citation into the
counted files, run against the proposed texts). Citations of PRINCIPLES.md
pinned to a commit (`docs/PRINCIPLES.md@94a157b2:...`) still resolve.

**Counted lines:** `rebuild_core.py` about 35, `session_state.py` about +6,
the four list entries -4, `new_worktree.py` 0: about **37**.

## 4. Words, by `tools/rule_sizes.py`'s count

"Now" is the merged tree at `8356c3a3` (`python3 tools/rule_sizes.py`);
"proposed" is `wc -w` over `docs/increments/h36-proposed/`, the same count
(`rule_sizes.words` is `len(text.split())`).

| File | Now | Proposed | Change |
|---|---|---|---|
| `CLAUDE.md` | 1175 | 286 | -889 |
| `.claude/REQUIRED-READING.md` | 1752 | 502 | -1250 |
| `docs/increments/README.md` | 990 | 412 | -578 |
| `docs/PRINCIPLES.md` | 437 | 0 | -437 |
| `.claude/briefs/common.md` | 244 | 227 | -17 |
| `architect.md` | 495 | 218 | -277 |
| `developer.md` | 251 | 174 | -77 |
| `orchestrator.md` | 504 | 310 | -194 |
| `perf.md` | 623 | 334 | -289 |
| `reviewer.md` | 454 | 328 | -126 |
| `tester.md` | 799 | 374 | -425 |
| skill `computational-geometry` | 346 | 176 | -170 |
| skill `geospatial-data-formats` | 249 | 109 | -140 |
| skill `modern-cxx` | 308 | 131 | -177 |
| skill `python-development` | 634 | 254 | -380 |
| skill `dispatch` (new) | 0 | 365 | +365 |
| skill `publishing` (new) | 0 | 158 | +158 |
| skill `away-mode` (new) | 0 | 272 | +272 |
| **Total** | **9261** | **4630** | **-4631 (50.0 % of now; half is 4630.5)** |

The total counts the new skills in full, though they load only on demand.
What each agent loads before its first act (`CLAUDE.md`, REQUIRED-READING,
the increments README, plus the brief and its persona file, or `dispatch`
for the main session):

| Agent | Now | Proposed |
|---|---|---|
| main session | 3917 | 1565 |
| `@architect` | 4656 | 1645 |
| `@developer` | 4412 | 1601 |
| `@orchestrator` | 4665 | 1737 |
| `@perf` | 4784 | 1761 |
| `@reviewer` | 4615 | 1755 |
| `@tester` | 4960 | 1801 |

## 5. The lean-process change (h34), merged first

h34 merged as #389, with its follow-ups #390 to #393 (question tool, session
file retired, lane follow-ups). Its rulings are now in the files themselves,
so each proposed text was rebased onto them (§7) rather than carried over
from h34's design. h34's `tests/python/test_lean_process_rule_text.py` pins
its texts word for word; this PR's red step deletes it, since most of those
texts no longer exist in that form. Its words are in "now" above, so §4's
total already includes them; no estimate remains.

## 6. Tests and steps

`tests/python/test_rebuild_core.py`, on a temporary tree holding a copy of
the tool, a fake `cmake` on `PATH` and a fake `.venv/bin/python`:

1. No `build-pyext/CMakeCache.txt`: exit 2, one line naming
   `tools/new_worktree.py`; the fake `cmake` never runs.
2. `cmake` exits 1: exit 1, one line naming the exit status; the destination
   is unchanged byte for byte.
3. `cmake` exits 0 and writes `_core<suffix>`: exit 0; the destination's
   bytes equal the built file's; the printed digest is their SHA-256's first
   12 hex digits.
4. The built file for the venv's suffix is missing (another suffix present):
   exit 2, nothing copied.

`tests/python/test_session_state.py`, with the existing `harness_repo`
fixture and `harness_fixtures` flag helpers (each recap is a few dozen lines;
the asserted line is the one right after `== recap ==`):

5. Flag off (no flag file): `Main session: load the dispatch skill before
   your first act.`
6. `flag_on`: `Main session: load the dispatch and away-mode skills before
   your first act.`
7. `flag_broken(repo, "not-json")`: as 6.
8. `flag_expired`: as 5.

Existing tests changed because `docs/PRINCIPLES.md` goes (the red step; each
keeps its intent):

- `tests/python/test_session_state.py@14a17091:731`: the copied list loses the file.
- `test_rule_sizes.py`: the `RULE_FILES` test (line 63) loses it; the
  `STARTING` fixture (line 131, 10 words) loses it, so its row assert (184)
  goes and the totals become 97 (lines 191-192, still `+7` on 90) and 97
  (line 229); the order test (line 205) loses it; line 293's anchor becomes
  `docs/increments/README.md`. Line 110 tests the pure `table()` with a
  literal dict and is unchanged.
- `tests/python/test_ci_changes.py@14a17091:85`: `READ_BY_CODE` loses it.
- `tests/python/test_check_citations_pinning.py@14a17091:570`: `GOVERNED` loses it.
- `tests/python/test_new_worktree.py@14a17091:328`: the `rebuild` line names `tools/rebuild_core.py`
  (replacing the `cmake --build` and `--target _core` asserts); the
  `build-pyext` assert goes.

`tests/python/test_lean_process_rule_text.py` is deleted in the red step
(§5).

Three suites read the real rule files and pass unchanged against the
proposed texts (checked by a script over `docs/increments/h36-proposed/`
that applies each suite's own matching: its regexes, `collapsed()` and
phrase lists; the suites run against the applied files in CI):
`tests/python/test_spawn_rule_text.py` (the step-order phrase "a performance
fix is timed by `@perf` before review;", the harness section naming
`guard_spawn.py` and `tools/brief.py`, the tester and reviewer phrases),
`tests/python/test_settings_wiring.py` (*The harness* has one paragraph
opening "Active in `.claude/settings.json`" that names `guard_unattended.py`,
`guard_push.py`, `guard_governance.py` and `gates_after_commit.py`), and
`tests/python/test_brief.py` through `brief_fixtures.PHRASES` (the handback
headings, "DEFAULTS TAKEN" included). No invariant-critical suite.

Steps: red (`@tester`), green (`@developer`, with Ola present for the
governed files: the tool, the recap line, the four list entries,
`new_worktree.py`, the count_loc docstring), then the main session, with Ola
present, copies the proposed texts into place, creates the three skills,
deletes `docs/PRINCIPLES.md` and edits h2's example; one `@reviewer` round
over all of it.

Left out: a `rule_sizes.py` column for the per-agent payload of §4 (counted
here by `wc -w`); the history paragraph in `guard_governance.py`'s
docstring.

## 7. Round 3: the proposed texts rebased onto 2026-10-28's rules

The rule files changed between `b2207408` and `14a17091`
(`git diff b2207408 14a17091 -- <file>`). Each change, and where the
proposed text now holds it:

| Today's rule | Proposed text |
|---|---|
| First question before an item's first spawn (`CLAUDE.md`) | `dispatch`, *Spawning* |
| Step order with full and fast lane (`CLAUDE.md`) | `CLAUDE.md`, word for word |
| Lessons go to the lessons log; a duty runs first with no free item (`CLAUDE.md`) | `dispatch` |
| Merge queue text without the method flag (`CLAUDE.md`) | `publishing` |
| Questions through `tools/questions.py`, items in `work_pipes.py`; `session.md` gone (REQUIRED-READING) | `dispatch`, *Work in flight*; the `session.md` bullet is gone |
| `@reviewer` before push only on a full-lane branch (REQUIRED-READING) | `publishing` |
| Writer re-reads the at-risk citations as quotations before any push (REQUIRED-READING) | REQUIRED-READING, approval boundary |
| "On an instruction's momentum, with no fresh approval": edits, spawns, builds, tests, gates, note files and commits, as far as the instruction requires (REQUIRED-READING) | REQUIRED-READING, approval boundary, first sentence after "the remote is Ola's" |
| Fast-lane branch for origin needs no grant or review while away; `questions.py add --act`, `taken`, `list`; fallback held in `work_pipes.py` | `away-mode` |
| Two lanes replace the short path; three red CI runs, ask Ola (README) | increments README, *Two lanes* |
| Developer description and fast-lane paragraph | `developer.md` |
| Tester: no short-path exception to *Choices beyond the design* | `tester.md` |
| Orchestrator: description, lessons log and weekly retrospective, no `-research` name, writes the log | `orchestrator.md` |
| Reviewer: semantic mandate, no CI precondition, build-workflow check, *Follow-ups* | `reviewer.md` |
| Handback heading "DEFAULTS TAKEN" (brief) | `briefs-common.md` |

**Dropped or rewritten from rounds 1 and 2,** because a rule merged today
replaced them:

- `reviewer.md` no longer gains the three prose-branch checks: h34 removed
  them (its R15); adding them back would have reverted Ola's ruling.
- `dispatch`'s `session.md` bullet and `away-mode`'s `ASK OLA:` and `QUEUE:`
  lines: replaced by `questions.py` and `work_pipes.py` (#385, #391).
- The README's short path, and every "short path" in this file: the short
  path no longer exists.
- §5's carry-over plan and its +60-word estimate: h34 is merged, so its
  words are counted in "now".

**Restored**, dropped by rounds 1 and 2 without a reason given in §1: the
brief's scratch-copy write grant (`25f01c49`) and "never `sleep` or a
polling loop" (`1d778553`); `dispatch`'s condition "with both writer slots
busy" on merging append-only blocks.

**Defects of round 2 found while rebasing**, now fixed: the proposed
REQUIRED-READING's *The harness* lacked the "Active in
`.claude/settings.json`" paragraph and `gates_after_commit.py`, which
`test_settings_wiring.py` pins (already so at `b2207408`); the step-order
phrase and the brief's "ASK OLA" heading would have failed today's
`test_spawn_rule_text.py` and `brief_fixtures.PHRASES`.

**Cut further to stay at half:** eleven `description` fields (the four
domain skills, `architect.md`, `orchestrator.md`, `perf.md`, `reviewer.md`,
`tester.md`, `publishing`, `away-mode`) cut to one short sentence;
`.claude/agents/tester.md@14a17091:40` (*Tiers*, a description of the suites), the duplicate
baseline sentence in `perf.md`, `CLAUDE.md`'s gate-hook sentence (now only
in REQUIRED-READING's hook list) and REQUIRED-READING's opening line and
*Ideas* handback clause (the brief's headings carry it) are gone.

**Sweep of today's REQUIRED-READING, paragraph by paragraph** (after the
round-3 review found the momentum paragraph unplaced). Each paragraph is now
kept, moved, replaced by a tool, or dropped as below:

- Kept, condensed: the lost-session case of deleting a note file and "or
  explicitly cancelled by Ola" (`dispatch`); "state the mechanism" where
  nothing can run (*Claims*); "idle is accepted only when no such item
  exists" (`away-mode`); "the permission system is not the push backstop:
  auto mode can let an unapproved push through" (REQUIRED-READING, *The
  harness*: "tripwires, not the boundary").
- Kept as a permission, not a restatement (design round 4): "an idea enters
  git only when picked for building" is REQUIRED-READING's "commit an idea
  only once picked" (+2 words, paid by `CLAUDE.md`'s "also" and "the
  growth"); "commit nothing from it" alone would forbid an increment file
  built from a picked idea.
- Fixed: `away-mode` required the head's sha under `## Review` in the
  increment file; today's rule is the verdict a hook records from
  `@reviewer`'s handback (`.claude/hooks/record_review.py`), so it says that.
- Dropped as restatement: "including a branch's first" (every push needs
  the yes); "binds claims about behaviour, not design opinions" (an opinion
  has no command to run); "A result another persona needs is named by path
  in the handback" (each agent's system prompt asks for absolute paths);
  "prompts Ola queued mid-turn" (the recap prints them); the `skills:`
  preload note (reasoning); "`pip install -e` is not the
  route" (part of the rebuild recipe the tool replaces).
- Replaced by a tool: the harness lines on the scratchpad pass, the duties
  hook, the status line and "the guard judges each redo again" (the hooks'
  own texts); the end of a window at the flag's time or `!` (`away.py`).
- To pay for the added words, these restatements go: `dispatch`'s "the
  recap prints the open list" and "Items go in `tools/work_pipes.py`" (named
  the line before); `perf.md`'s "Report measured figures only" (its
  description and §4 say it); `publishing`'s "A fast-lane branch has no
  review" (the README's fast lane says it); `away-mode`'s "A subagent puts
  both in its handback" (the brief's headings); `CLAUDE.md`'s "CI runs it
  with `--max 700`" (CI's workflow states it).

**Gate fixes after the copy into place** (both missed by the design
reviews, which ran the rule-text suites but not the governance gates on the
proposed texts): `tools/check_prohibited_deps.py` parses `CLAUDE.md`'s list
only in the form "**Prohibited Dependencies:** Never introduce <list>.", so
the bullet takes that form (+2 words, paid by "The C++ targets build" becoming
"C++ builds"); `tools/check_citations.py` found four broken citations, the
two here now pinned to `14a17091`.

## Review

- Design round 1: CHANGES REQUESTED, b2207408..439e8f81, 0 net production lines (reviewer rules unreachable in a skill; six files naming PRINCIPLES.md; a false hook claim; nothing loads dispatch; scope beyond the quoted words).
- Design round 2: CHANGES REQUESTED, ..626d5e83, 0 (all five fixed; only Ola's ruling on the wider scope remains).
- Design round 3: APPROVED WITH WORDINGS, 626d5e83..20f0338d, 0 (today's rules all kept where section 7 says; dropped items gone from master; rule-text suites pass on the proposed texts, proven by a planted break; 9261 against 4629 words, 50.0 %). Wording applied in this commit: `h36-proposed/skill-dispatch.md:17` "* Brief each spawn and new-step resume: `tools/brief.py`; make worktrees with" (keeps "a resume that starts a new step carries a fresh block"). Follow-up for @architect: section 1 does not list REQUIRED-READING's "On an instruction's momentum" paragraph, and no proposed text keeps it.
- Design round 4: APPROVED, ab1add36..fa53da14, 0 (each paragraph's fate in section 7 true against today's REQUIRED-READING; the six cuts each said by a line that stays; rule-text suites 153 passed, 12 skipped on the proposed texts, a planted break failed them; 4630 words). Follow-ups for the next writer: list the harness sentence "auto mode can let an unapproved push through" in the sweep (kept in "tripwires, not the boundary"); reclassify "an idea enters git only when picked" as a permission, not a restatement.
- Code review: APPROVED, b1e39351..7142b6e9, 45 (every gate exits 0; full Python suite 11271 passed, 18 skipped, 5 expected failures; all 17 live rule files byte-identical to their h36-proposed/ copies; green commit touches no test file; seven skill names match their folders; 4630 words; merges cleanly with origin/master 22d09fe1). Follow-ups, filed as pipe items: delete docs/increments/h36-proposed/ with its citation pin, the section 1/4/6 pointers and the ROADMAP row (45 lines); `publishing` asks a new review after any verdict other than APPROVED while README step 4 lets a wordings commit pass, a disagreement older than h36.
