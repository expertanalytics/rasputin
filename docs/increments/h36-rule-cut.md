# Harness h36: cut the rule payload by half

Status: **designed**, 2026-10-27, off master `b2207408`; round 2 after
`@reviewer`'s CHANGES REQUESTED on `439e8f81`. Short path under today's
rules: about 37 counted lines (§3), no C++ file, no refine or mesh code. Most
of the PR is rule text, and six of the changed tools are governed
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

**Defaults the main session took while Ola is away** (defaults, not his
rulings; he can reverse either):

- Scope: every file `tools/rule_sizes.py` counts, so `docs/increments/README.md`,
  `docs/PRINCIPLES.md` and the four domain skills are in, beside the files
  the review names; the new rebuild tool is in too.
- The lean-process rulings (§5) are carried in shortened; Ola sees each
  shorter wording at the final review.

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
tester.md cut" (the cut is taken up here: `tester.md` 811 to 409 words) and
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
| `CLAUDE.md` | §1 to §4, each cut to its commands; step order | dispatch rules of §3 → `dispatch`; CI and merge queue → `publishing` | size-limit definition → `tools/count_loc.py` docstring; "merge commit, never squash" → the merge queue's method; `@orchestrator` timing → `tools/duties.py`; gate list → `gates_after_commit.py`, except `check_code_history.py`, which the hook does not run (`gates_after_commit.py` `gates()`), so its command stays in §4's block | method-flag and `-d` behaviour of `gh pr merge` |
| `REQUIRED-READING.md` | note file, skills, ideas, claims, approval boundary (four lines), harness (four lines), data folders | cold start, `session.md` → `dispatch`; grants, reviewer-before-push → `publishing`; unattended mode → `away-mode` | rebuild recipe → `tools/rebuild_core.py`; hook descriptions → the hooks' own refusal texts | recap-size rule; `touch` after a `cp` restore (mutants now run in scratch copies); "editing rule files is ordinary tree work" |
| `docs/increments/README.md` | the loop, short path, acceptance pointer, cost constraints | verdict format → `reviewer.md`; acceptance details → `perf.md` | merge commit → the merge queue | "never dropped" list beyond red before green (each is stated where it applies) |
| `docs/PRINCIPLES.md` | deleted; A6, B3, B4, D2, D3 become five lines in REQUIRED-READING's *Claims*; the six files that name it change (§3) | | | *Origin* and *Status* lines |
| `.claude/briefs/common.md` | all; tested phrases kept word for word | | | "never by a name pattern" |
| `architect.md` | §1 to §5 condensed | | | generic pillar prose |
| `developer.md` | all, condensed | | | |
| `orchestrator.md` | all, condensed | | | |
| `perf.md` | all, condensed; gains the acceptance run | | `--no-build` rebuild → `rebuild_core.py` | py-spy hint |
| `reviewer.md` | all, condensed; gains the verdict format and REQUIRED-READING's three prose-branch checks (run what each claim is about; re-read the at-risk citations as quotations; name each new rule's incident) | | | |
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
the remote is Ola's without loading anything. The three prose-branch review
checks live in `reviewer.md`, the one file the agent that runs them loads;
`publishing` only says the review comes first.

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
stays above `== recap ==` (`test_session_state.py` `_header`, line 326) and
the duties line stays right after `Last landed:` (`test_duties.py`
`_after_last_landed`).

**Changed: the six files that name `docs/PRINCIPLES.md` or the recipe**
(found with `git grep -n PRINCIPLES -- tools .claude/hooks tests` and
`grep -n recipe tools/new_worktree.py`):

| File | Change | Counted lines |
|---|---|---|
| `.claude/hooks/guard_governance.py@b2207408:60` | drop the entry from `GOVERNED`; its docstring's line 14 mention stays (history, uncounted) | -1 |
| `tools/rule_sizes.py:25` (governed) | drop it from `RULE_FILES` | -1 |
| `tools/ci_changes.py:38` (governed) | drop it from `NOT_PROSE` and from the comment above | -1 |
| `tools/check_citations.py@b2207408:96` (governed) | drop it from `GOVERNED` | -1 |
| `tools/new_worktree.py:11, 125` | the docstring names `tools/rebuild_core.py` instead of REQUIRED-READING's recipe; the printed `rebuild` line becomes `python3 <wt>/tools/rebuild_core.py   (configured, not built)` | 0 |
| `tests/python/test_session_state.py:923` | stop copying the file | tests |

Four more test files pin the file name and change with it (§6).

**Changed: `docs/increments/h2-citation-pinning.md` line 109**: its example
`PRINCIPLES.md §A6` is cut, since the file goes (the only citation the move
breaks; checked by a script over every ID and heading citation into the
counted files, run against the proposed texts). Citations of PRINCIPLES.md
pinned to a commit (`docs/PRINCIPLES.md@94a157b2:...`) still resolve.

**Counted lines:** `rebuild_core.py` about 35, `session_state.py` about +6,
the four list entries -4, `new_worktree.py` 0: about **37**, under the short
path's 50.

## 4. Words, by `tools/rule_sizes.py`'s count

| File | Now | Proposed | Change |
|---|---|---|---|
| `CLAUDE.md` | 1155 | 308 | -847 |
| `.claude/REQUIRED-READING.md` | 1857 | 471 | -1386 |
| `docs/increments/README.md` | 1060 | 432 | -628 |
| `docs/PRINCIPLES.md` | 437 | 0 | -437 |
| `.claude/briefs/common.md` | 244 | 206 | -38 |
| `architect.md` | 495 | 219 | -276 |
| `developer.md` | 217 | 142 | -75 |
| `orchestrator.md` | 482 | 260 | -222 |
| `perf.md` | 623 | 359 | -264 |
| `reviewer.md` | 565 | 343 | -222 |
| `tester.md` | 811 | 409 | -402 |
| skill `computational-geometry` | 346 | 197 | -149 |
| skill `geospatial-data-formats` | 249 | 126 | -123 |
| skill `modern-cxx` | 308 | 153 | -155 |
| skill `python-development` | 634 | 277 | -357 |
| skill `dispatch` (new) | 0 | 342 | +342 |
| skill `publishing` (new) | 0 | 161 | +161 |
| skill `away-mode` (new) | 0 | 276 | +276 |
| **Total** | **9483** | **4681** | **-4802 (49.4 % of now)** |

Counted with `wc -w` over `docs/increments/h36-proposed/` (round 2: the
reviewer gains the three prose-branch checks, `publishing` loses them,
`CLAUDE.md` names the hook's gates and the code-history command, `dispatch`
is trimmed; the total is unchanged). The total counts the new skills in
full, though they load only on demand.
What each agent loads before its first act (`CLAUDE.md`, REQUIRED-READING,
the increments README, plus the brief and its persona file, or `dispatch`
for the main session):

| Agent | Now | Proposed |
|---|---|---|
| main session | 4072 | 1553 |
| `@architect` | 4811 | 1636 |
| `@developer` | 4533 | 1559 |
| `@orchestrator` | 4798 | 1677 |
| `@perf` | 4939 | 1776 |
| `@reviewer` | 4881 | 1760 |
| `@tester` | 5127 | 1826 |

## 5. With the lean-process change (h34), which lands first by default

h34 (branch `worktree-lean-process`, waiting for Ola) edits the same files
with 13 old-to-new texts. This PR then merges master, and the main session
carries each h34 ruling into the proposed text, condensed:

- R1 (lessons go to the lessons log): `dispatch` already says only "pass
  handback lessons to `@orchestrator`"; no change.
- R2 to R6 (`reviewer.md`): the opening sentence becomes h34's semantic
  mandate; "CI first" goes; check 3 becomes "a change to the build is
  exercised by the workflow"; "Suggestions" becomes "Follow-ups, filed by
  your spawner as pipe items".
- R7 to R10 (`orchestrator.md`): h34's description; the `-research` file
  name goes; §3 becomes h34's lessons log and weekly retrospective; §4 adds
  "and the log".
- R11 to R13 (README): the short path gains "no rule text, no new module or
  public type", the one-round rule, and its pointer from step 4.
- h34's `tests/python/test_lean_process_rule_text.py` pins h34's texts word
  for word; this PR's red step deletes it, since those texts no longer exist.

Estimated effect: about +60 words, so about 4,741, half of now, with no
margin. An estimate: if the applied total (`python3 tools/rule_sizes.py`)
comes out above 4,741, the main session cuts the excess from the carried
wordings before the final review, where Ola sees them (default above).
With h34 first, this PR no longer takes the short path (h34 bars rule text
from it).

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

- `test_session_state.py:923`: the copied list loses the file.
- `test_rule_sizes.py`: the `RULE_FILES` test (line 63) loses it; the
  `STARTING` fixture (line 131, 10 words) loses it, so its row assert (184)
  goes and the totals become 97 (lines 191-192, still `+7` on 90) and 97
  (line 229); the order test (line 205) loses it; line 293's anchor becomes
  `docs/increments/README.md`. Line 110 tests the pure `table()` with a
  literal dict and is unchanged.
- `test_ci_changes.py:85`: `READ_BY_CODE` loses it.
- `test_check_citations_pinning.py:570`: `GOVERNED` loses it.
- `test_new_worktree.py:328`: the `rebuild` line names `tools/rebuild_core.py`
  (replacing the `cmake --build` and `--target _core` asserts); the
  `build-pyext` assert goes.

`tests/python/test_spawn_rule_text.py` and `tests/python/test_brief.py`
pass unchanged: the phrases they pin are kept (checked with `grep` against
the proposed files; the suites run against the applied files in CI). No
invariant-critical suite.

Steps: red (`@tester`), green (`@developer`, with Ola present for the
governed files: the tool, the recap line, the four list entries,
`new_worktree.py`, the count_loc docstring), then the main session, with Ola
present, copies the proposed texts into place, creates the three skills,
deletes `docs/PRINCIPLES.md` and edits h2's example; one `@reviewer` round
over all of it.

Left out: a `rule_sizes.py` column for the per-agent payload of §4 (counted
here by `wc -w`); the history paragraph in `guard_governance.py`'s
docstring.

## Review

- Design round 1: CHANGES REQUESTED, b2207408..439e8f81, 0 net production lines (reviewer rules unreachable in a skill; six files naming PRINCIPLES.md; a false hook claim; nothing loads dispatch; scope beyond the quoted words).
- Design round 2: CHANGES REQUESTED, ..626d5e83, 0 (all five fixed; only Ola's ruling on the wider scope remains).
