# Harness h34: lessons log, semantic review, one review round

Status: **designed**, 2026-10-27, off master `662022b0`. Short path under
today's rules (`docs/increments/README.md`): about 15 counted lines in
`tools/duties.py` and `tools/brief.py`, no C++ file, no refine or mesh code.
`tools/brief.py` is governed, so the green step runs with Ola at the guard's
prompt, by day. The rule texts of §5 are applied in the same PR by the main
session with Ola's yes, as h33's were.

Source: an outside review of the workflow that Ola pasted on 2026-10-27 and
chose to adopt, prefaced by "we should do this asap". Its points 1, 3 and 4,
verbatim:

> "Abolish Automated Retrospective PRs: Restrict @orchestrator from filing
> standalone PRs or editing governance docs on every cycle. Have it append
> short Markdown records to a single uncommitted scratchpad or aggregate
> lessons on a strict weekly cycle. Cut the retrospective overhead
> immediately."
>
> "Refocus @reviewer on Semantics: Strip procedural checklist items from
> @reviewer. Let count_loc.py, CI, and static hooks enforce LOC limits,
> mutant targets, and formatting automatically. Instruct @reviewer to
> evaluate only two things: domain model validity and edge-case handling."
>
> "Relax the Micro-Increment Churn: Five review rounds for an 88-line PR
> point to conversational thrashing over trivialities. Raise the review
> threshold for trivial/short-path changes to let minor structural fixes
> land in a single round."

Defaults the main session took while Ola was away; they are not his ruling,
and he may reverse them: the research round moves into the weekly
retrospective (§1); non-blocking findings become pipe items on every path,
not only the short one (§2); the reviewer keeps a check whether the change's
claims are true, and keeps the procedural checks no tool enforces yet
(§2), until a tool does.

## Prior art: legacy and literature

*Literature.* None applies: these are repository process rules, not a method
with a published form. No novelty is claimed.

*Legacy.* `git grep -n -i -E 'retrospective|reviewer|lessons' legacy-archive -- legacy`
returns nothing; nothing is carried across.

*Ideas file.* No entry is taken up or overlapped. One new idea comes out of
§3 (a CI gate for the size ceiling), in the handback.

## 1. No retrospective PR per merge or per night

`@orchestrator` still runs when `tools/duties.py` shows a duty owed, but its
post-merge and morning checks file no branch and no PR. Each appends one
entry to the **lessons log**, `../rasputin_scratch/lessons.md`: the scratch
folder beside the main checkout, outside git, shared by every worktree and
kept when one is removed. An entry is a heading
`## <YYYY-MM-DDTHH:MMZ> merge #N` or `## <YYYY-MM-DDTHH:MMZ> night` (UTC), then
at most five lines, each a finding or lesson with its evidence. A finding
that needs a fix goes in that check's handback as a follow-up, and the main
session files it as a pipe item then, not after the weekly summary. Once a
week one PR, a dated file under `docs/retrospectives/`, groups the entries
since the last one by theme and holds the research round, the rule-text
measurement and its proposed cut; it carries no inflow, KGE or
precipitation-correction result (Ola's rule: those stay private). The
research round is no longer its own file.

The log is `@orchestrator`'s own data, written and read by it and read by a
tool; it is addressed to no other agent, so the "not a channel" rule of
`.claude/REQUIRED-READING.md` needs no exception. The path is in scratch, not
under `.git/`: Claude Code may ask before an edit under `.git/`, and an
unattended night refuses every prompt. Ola adds the log's row to
`rasputin_scratch/INDEX.md` (`tools/index.py` lists it as "Not yet described").

**`tools/duties.py`** (the hook line, the status line and the recap all read
`duties.line`, so all three follow):

- New `lessons_path(repo) -> Path`: the parent of
  `git rev-parse --path-format=absolute --git-common-dir`, then
  `../rasputin_scratch/lessons.md` beside it. Same path from every worktree.
- New pure `logged(text: str, now: datetime) -> datetime | None`: the latest
  stamp among lines matching
  `^## (\d{4}-\d{2}-\d{2}T\d{2}:\d{2}Z) (merge #\d+|night)\b`, in UTC, and
  not later than `now`; a later stamp (a typo, a wrong clock) is ignored, so
  it cannot silence the duties. Other lines are ignored; no match gives
  `None`. A missing file reads as empty. The entry is stamped from
  `date -u +%Y-%m-%dT%H:%MZ` at the start of the check, so a merge that lands
  while the check runs is still owed.
- `Inputs` gains `logged`; `research` goes.
- The **baseline** is the later of `logged` and `retro` (the latest
  retrospective first added by a check commit, as today), so the switch-over
  owes nothing for merges a retrospective already covered.
- Post-merge: merges after the baseline, less recording merges as today,
  give `post-merge check (N merge(s) since the last log entry, <minute>:
  <labels>)`; with no baseline, `post-merge check (no log entry on record)`.
- Morning: away windows of four hours or more that ended after the baseline,
  same text as today.
- The research item becomes `weekly retrospective (last YYYY-MM-DD)` when
  `retro` is seven days old or more, `weekly retrospective (none on record)`
  when there is none. Nothing else asks for a research round.
- The push item (local branches adding a retrospective master lacks) is
  unchanged; it now names the weekly branch.
- The `git log --since` of the merge query uses the baseline.

**`tools/brief.py`**: `WRITES["orchestrator"]` becomes
`("docs/retrospectives/", "the lessons log ../rasputin_scratch/lessons.md")`.

## 2. `@reviewer` judges semantics

Its mandate becomes: is the domain model valid (types, invariants, units,
degenerate inputs), are edge cases handled, is every claim the change makes
true. Each procedural check, with its enforcer:

| Check in `reviewer.md` today | Enforcer | In this design |
|---|---|---|
| CI status, including a prose-only PR's skipped jobs | Merge queue: `CI result` is master's one required check (`gh api repos/expertanalytics/rasputin/branches/master/protection/required_status_checks` returns `["CI result"]`) | removed |
| Readability, typing, style | ruff, ruff format, mypy in CI | removed (the sentence that waived them) |
| Citations resolve | `tools/check_citations.py`, governance job of CI | not in `reviewer.md`; re-reading the at-risk quotations stays in `.claude/REQUIRED-READING.md`, since no tool can |
| Verdict recorded | `.claude/hooks/record_review.py` and the push grant | the `Review:` line stays: it is the hook's input |
| Size ceiling | none: `tools/count_loc.py` measures, but CI only runs its tests (`grep -n count_loc .github/workflows/main.yaml`: comments only) | **kept** until a CI step refuses a count over 700 |
| LOC reconciled against the estimate | none | removed: no result depends on it; the ceiling, kept, is the part that does |
| Red-step scaffolding gone | none | kept |
| Mutant kill record covers the named targets | none | kept |
| `@perf` acceptance run recorded | none | kept |
| The workflow exercises a changed build | none | kept, one line |

The design-review checks (speed estimate, repair method and tolerance, kills
against fixture figures, scope against Ola's words) are judgements of the
design's content and stay.

Non-blocking findings become follow-up items in a pipe
(`tools/work_pipes.py add`, filed by the spawner), on every path, in place of
the private ideas file.

## 3. One review round on the short path

The short path (`docs/increments/README.md`) is today: under 50 net counted
lines, no C++ file, nothing the *Acceptance* section covers. It gains: no rule
text, and no new module or public type (the "no new concept" of the order).

`@reviewer` measures the real diff with `tools/count_loc.py` and checks
every condition above on it; if one fails, the change goes to the full path.
Such a change lands after its one combined round. Only three findings block,
and only they give `CHANGES REQUESTED`: a wrong result, a broken invariant, a
failing test. The round after their fix checks that fix alone and raises no
new non-blocking finding. Every other finding is a follow-up pipe item. The
away-mode push grant is unchanged: it still wants an approving verdict on the
head.

## 4. Tests that pin it (`tests/python/test_duties.py`, `test_brief.py`)

Fixture times are UTC; `NOW = 2026-10-27T18:00Z`.

1. `logged`: text with headings `## 2026-10-26T09:00Z merge #380`,
   `## 2026-10-27T12:00Z night`, `## 2026-10-27T13:00 merge #381` (no `Z`)
   and a body line `## not a heading` gives `2026-10-27T12:00Z` with
   `now = NOW`; `""` gives `None`.
2. A future stamp: adding `## 2026-10-28T09:00Z night` (15 hours after
   `NOW`) still gives `2026-10-27T12:00Z`; text holding only that heading
   gives `None`.
3. Post-merge after the log: `retro = 2026-10-20T10:00Z`,
   `logged = 2026-10-27T12:00Z`, code merges `#381` at 11:00 and `#382` at
   13:00 on 2026-10-27 owe exactly
   `post-merge check (1 merge since the last log entry, 2026-10-27 12:00 UTC: #382)`.
4. Baseline from the retrospective: `logged = None`, `retro =
   2026-10-27T12:00Z`, same merges: the same single item for `#382`; both
   `None`: `post-merge check (no log entry on record)`.
5. Morning: a five-hour window ending 2026-10-27T11:00Z owes nothing with
   `logged = 12:00Z`; one ending 13:00Z owes the morning item.
6. Weekly: `retro = 2026-10-19T10:00Z` gives
   `weekly retrospective (last 2026-10-19)`; `2026-10-21T10:00Z` gives no
   weekly item; `retro = None` gives `weekly retrospective (none on record)`.
   No item contains `research`.
7. `lessons_path` from a repository and from a worktree of it is the same
   path, `<parent of the main checkout>/rasputin_scratch/lessons.md`.
8. The orchestrator brief's write limit names `../rasputin_scratch/lessons.md`.

The existing tests that use `research=` or the "since the last retrospective"
text are amended by `@tester` in the red commit. No invariant-critical suite.

Left out: a CI gate for the size ceiling (an idea for Ola); deleting entries
from the log; any change to `tools/session_state.py`.

## 5. Rule texts, old to new

Applied verbatim by the main session with Ola's yes. Each old text occurs
once in its file at `662022b0`. Word change, by `tools/rule_sizes.py`'s count
on copies with these replacements applied: see §6.

#### R1 `CLAUDE.md`

````text
* **Lessons.** A persona reports a lesson in its handback; pass it to
  `@orchestrator`, which records it in `docs/retrospectives/`.
````

````text
* **Lessons.** A persona reports a lesson in its handback; pass it to
  `@orchestrator`, which appends it to the lessons log.
````

#### R2 `.claude/agents/reviewer.md`

````text
description: Final quality gatekeeper. Audits CI status, the LOC ceiling, red-step scaffolding, prose claims against code and the mutant kill record of invariant-critical suites, returning APPROVED or CHANGES REQUESTED. Read-only by design. Use before pushing anything.
````

````text
description: Final quality gatekeeper. Judges whether the domain model is valid, edge cases are handled and the change's claims are true, plus the procedural checks no tool makes yet. Read-only by design. Use before pushing anything.
````

#### R3 `.claude/agents/reviewer.md`

````text
Your mandate is the checks no gate makes: CI status, the size ceiling, leftover red-step scaffolding, prose claims that the change made false, and the mutant kill record of a suite named invariant-critical. Readability, typing and style are the gates' (`CLAUDE.md` §4); do not spend the round on them.
````

````text
Your mandate is semantics: is the domain model valid (types, invariants, units, degenerate inputs), are edge cases handled, are the change's claims true. Tools enforce CI status, style (the required `CI result` check) and your verdict's record (a hook reads your `Review:` line). The checks below are yours until a tool makes them.
````

#### R4 `.claude/agents/reviewer.md`

````text
### Precondition: CI status
Before any verdict, check what CI says (`CLAUDE.md` §4: CI is authoritative):
```bash
gh pr checks <pr>              # or: gh run list --branch <branch> --limit 1
```
**Red CI is an automatic `CHANGES REQUESTED`**, and so is a workflow that does
not exercise the current build: read `.github/` when the change touches the
build. On a prose-only PR (`CLAUDE.md` §4) the
C++, sanitizer and Python checks show as skipped; that is the design, not a
red, and `CI result` is the check that decides.

### The checks the gates cannot make
````

````text
### The checks no tool makes yet
````

#### R5 `.claude/agents/reviewer.md`

````text
3. **Actual LOC is reconciled against the increment doc's estimate.** Measure it;
   if the design named a split seam for an overrun, check whether it should fire.
````

````text
3. **A change to the build is exercised by the workflow** (read `.github/`).
````

#### R6 `.claude/agents/reviewer.md`

````text
4. **Suggestions:** Non-blocking, and only where no gate would catch it. One
   the branch does not take goes into the private ideas file, added by your spawner.
````

````text
4. **Follow-ups:** non-blocking findings; your spawner files each as a pipe
   item (`tools/work_pipes.py add`).
````

#### R7 `.claude/agents/orchestrator.md`

````text
description: Workflow watcher and retrospective owner. Checks how the agents worked and flags deviations from the rules, researches agentic design and proposes improvements with evidence, and records lessons in docs/retrospectives/, the only place it writes. Use after an increment merges, the morning after an unattended night, or for a research round.
````

````text
description: Workflow watcher and retrospective owner. Checks how the agents worked, appends deviations and lessons to the lessons log, and once a week files them with a research round as one retrospective. Use when tools/duties.py shows a duty owed.
````

#### R8 `.claude/agents/orchestrator.md`

````text
it would cost to adopt. Read the source before you cite it; mark anything
from memory as unchecked. A research round's retrospective file name
contains `-research`.
````

````text
it would cost to adopt. Read the source before you cite it; mark anything
from memory as unchecked.
````

#### R9 `.claude/agents/orchestrator.md`

````text
## 3. Own the retrospectives

`docs/retrospectives/`, `next.md` included, is yours. Record the lessons the
main session passes you (`CLAUDE.md` §3, *The main session dispatches*).
Quotes, dates and incidents go here. Each retrospective measures the rule
text (`python3 tools/rule_sizes.py`) and proposes a cut, restatements
included, and flags anecdotes in rule files (a case, incident or
example), moved to their owner or out.
Once a week, check the private ideas file for entries gone stale or done unmarked.
````

````text
## 3. The lessons log and the weekly retrospective

After a merge or a night, append to `../rasputin_scratch/lessons.md` (never
committed) `## <stamp> merge #N` (or `night`), stamped at the start by
`date -u +%Y-%m-%dT%H:%MZ`, then at most five lines of findings (§1) and
lessons passed to you. No PR. A finding needing a fix goes in your handback
too, for a pipe item.

Weekly, one PR: a dated file in `docs/retrospectives/` (yours, `next.md`
included) grouping the log by theme, with the research round (§2), no
inflow, KGE or pcorr result, the rule text measured
(`tools/rule_sizes.py`), a cut proposed and anecdotes in rule files flagged.
Check the private ideas file for stale or done entries.
````

#### R10 `.claude/agents/orchestrator.md`

````text
- **Write only under `docs/retrospectives/`.** New ideas, and stale or
````

````text
- **Write only under `docs/retrospectives/` and the log.** New ideas, and stale or
````

#### R11 `docs/increments/README.md`

````text
An increment estimated under 50 net lines (`CLAUDE.md` §2) that touches no
C++ file and nothing the *Acceptance* section covers takes the short path;
````

````text
An increment estimated under 50 net lines (`CLAUDE.md` §2) that touches no
C++ file, rule text or *Acceptance* code, nor adds a module or public type, takes the short path;
````

#### R12 `docs/increments/README.md`

````text
  still goes to `@architect` before green.
````

````text
  still goes to `@architect` before green.
- `@reviewer` checks these conditions on the real diff (`count_loc.py`);
  one failing means the full path. Else it lands after that round: only a
  wrong result, a broken invariant or a failing test blocks, and the next
  round checks that fix alone; other findings are pipe items.
````

#### R13 `docs/increments/README.md`

````text
   any other finding, is `CHANGES REQUESTED`.
````

````text
   any other finding, is `CHANGES REQUESTED` (short path: below).
````

## 6. Rule text size

Counted by `tools/rule_sizes.py`'s `current()` on copies of its files with
R1 to R13 applied (each old text found exactly once first):

| File | Words now | After | Change |
|---|---|---|---|
| `CLAUDE.md` | 1155 | 1157 | +2 |
| `docs/increments/README.md` | 1060 | 1116 | +56 |
| `.claude/agents/orchestrator.md` | 482 | 503 | +21 |
| `.claude/agents/reviewer.md` | 565 | 454 | -111 |
| All rule files | 9483 | 9451 | **-32** |
