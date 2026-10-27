# Harness h36: cut the rule payload by half

Status: **designed**, 2026-10-27, off master `b2207408`. Short path under
today's rules: about 35 counted lines (`tools/rebuild_core.py`), no C++ file,
no refine or mesh code. Most of the PR is rule text, which `guard_governance.py`
guards, so the main session applies it with Ola present.

Source: an outside review Ola adopted on 2026-10-27 and asked for "asap" (the
main session's report of his words). The point taken up: "Slash Rule Payload
by 50%: Audit CLAUDE.md, REQUIRED-READING.md, and the persona briefs. Remove
historical incident narratives (e.g., commit references explaining why Boost
was banned). Retain only the active imperative command. Shift secondary
guidance into dynamic skills that load strictly on-demand." Ola also said he
likes tools that do one job and would build harness tools where they save
context.

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
| `CLAUDE.md` | §1 to §4, each cut to its commands; step order | dispatch rules of §3 → `dispatch`; CI and merge queue → `publishing` | size-limit definition → `tools/count_loc.py` docstring; "merge commit, never squash" → the merge queue's method; `@orchestrator` timing → `tools/duties.py`; gate list → `gates_after_commit.py` | method-flag and `-d` behaviour of `gh pr merge` |
| `REQUIRED-READING.md` | note file, skills, ideas, claims, approval boundary (four lines), harness (four lines), data folders | cold start, `session.md` → `dispatch`; grants, reviewer-before-push → `publishing`; unattended mode → `away-mode` | rebuild recipe → `tools/rebuild_core.py`; hook descriptions → the hooks' own refusal texts | recap-size rule; `touch` after a `cp` restore (mutants now run in scratch copies); "editing rule files is ordinary tree work" |
| `docs/increments/README.md` | the loop, short path, acceptance pointer, cost constraints | verdict format → `reviewer.md`; acceptance details → `perf.md` | merge commit → the merge queue | "never dropped" list beyond red before green (each is stated where it applies) |
| `docs/PRINCIPLES.md` | deleted; A6, B3, B4, D2, D3 become five lines in REQUIRED-READING's *Claims* | | | *Origin* and *Status* lines |
| `.claude/briefs/common.md` | all; tested phrases kept word for word | | | "never by a name pattern" |
| `architect.md` | §1 to §5 condensed | | | generic pillar prose |
| `developer.md` | all, condensed | | | |
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
| `dispatch` | recovery, spawning, work in flight, reports to Ola | by the main session before its first act (`CLAUDE.md` says so); never by a subagent |
| `publishing` | one yes per act, reviewer before push, conflict check, merge queue, CI | before asking Ola to push or merge, before any push, PR or merge, before calling a branch merge-ready |
| `away-mode` | unattended mode, refusals, grants, defaults while away | when the recap says unattended mode is on, when Ola says he is leaving, when a guard refuses an act |

The approval boundary itself stays in REQUIRED-READING: every agent must know
the remote is Ola's without loading anything.

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

**Changed: `docs/increments/h2-citation-pinning.md` line 109**: its example
`PRINCIPLES.md §A6` is cut, since the file goes (the only citation the move
breaks; checked by a script over every ID and heading citation into the
counted files, run against the proposed texts).

## 4. Words, by `tools/rule_sizes.py`'s count

| File | Now | Proposed | Change |
|---|---|---|---|
| `CLAUDE.md` | 1155 | 303 | -852 |
| `.claude/REQUIRED-READING.md` | 1857 | 471 | -1386 |
| `docs/increments/README.md` | 1060 | 432 | -628 |
| `docs/PRINCIPLES.md` | 437 | 0 | -437 |
| `.claude/briefs/common.md` | 244 | 206 | -38 |
| `architect.md` | 495 | 219 | -276 |
| `developer.md` | 217 | 142 | -75 |
| `orchestrator.md` | 482 | 260 | -222 |
| `perf.md` | 623 | 359 | -264 |
| `reviewer.md` | 565 | 309 | -256 |
| `tester.md` | 811 | 409 | -402 |
| skill `computational-geometry` | 346 | 197 | -149 |
| skill `geospatial-data-formats` | 249 | 126 | -123 |
| skill `modern-cxx` | 308 | 153 | -155 |
| skill `python-development` | 634 | 277 | -357 |
| skill `dispatch` (new) | 0 | 354 | +354 |
| skill `publishing` (new) | 0 | 188 | +188 |
| skill `away-mode` (new) | 0 | 276 | +276 |
| **Total** | **9483** | **4681** | **-4802 (49.4 % of now)** |

The total counts the new skills in full, though they load only on demand.
What each agent loads before its first act (`CLAUDE.md`, REQUIRED-READING,
the increments README, plus the brief and its persona file, or `dispatch`
for the main session):

| Agent | Now | Proposed |
|---|---|---|
| main session | 4072 | 1560 |
| `@architect` | 4811 | 1631 |
| `@developer` | 4533 | 1554 |
| `@orchestrator` | 4798 | 1672 |
| `@perf` | 4939 | 1771 |
| `@reviewer` | 4881 | 1721 |
| `@tester` | 5127 | 1821 |

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

Estimated effect: about +60 words, total still under 4,741 (half of now).
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

`tests/python/test_spawn_rule_text.py` and `tests/python/test_brief.py`
pass unchanged: the phrases they pin are kept (checked with `grep` against
the proposed files; the suites run against the applied files in CI). No
invariant-critical suite.

Steps: red (`@tester`), green (`@developer`: the tool and the count_loc
docstring), then the main session, with Ola present, copies the proposed
texts into place, creates the three skills, deletes `docs/PRINCIPLES.md`
and edits h2's example; one `@reviewer` round over all of it.

Left out: a `rule_sizes.py` column for the per-agent payload of §4 (counted
here by a scratch script); changing `new_worktree.py`'s printed rebuild
recipe to name the new tool.
