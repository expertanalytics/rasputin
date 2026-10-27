"""The rule texts the lean-process increment changes, as its design says they will read.

Spec: `docs/increments/h34-lean-process.md` §5, R1 to R13 (old text, then new),
each pair copied here with its whitespace collapsed. Each check reads the real
file of the checkout: the new text is there and the old one is gone (unless
the new text contains it, as R12's does). Whitespace is collapsed on both
sides, so line wrapping does not count.

The rule texts are applied by the main session with Ola's yes, after the green
step of `tools/duties.py` and `tools/brief.py`; until then every case fails,
so each is a strict expected failure. Whoever applies R1 to R13 removes the
`xfail` mark in the same commit: a strict `xfail` that passes fails the run.
"""

from __future__ import annotations

import pytest

from harness_fixtures import REAL

pytestmark = pytest.mark.harness

#: Applied with Ola's yes, not by the green step (the design's status line).
PENDING = "R1 to R13 are applied by the main session with Ola's yes"


def collapsed(text: str) -> str:
    return " ".join(text.split())


#: §5's pairs: (file, old text, new text), whitespace collapsed.
RULE_TEXTS = [
    pytest.param(
        "CLAUDE.md",
        (
            "* **Lessons.** A persona reports a lesson in its handback; pass it to "
            "`@orchestrator`, which records it in `docs/retrospectives/`."
        ),
        (
            "* **Lessons.** A persona reports a lesson in its handback; pass it to "
            "`@orchestrator`, which appends it to the lessons log."
        ),
        id="R1",
    ),
    pytest.param(
        ".claude/agents/reviewer.md",
        (
            "description: Final quality gatekeeper. Audits CI status, the LOC ceiling, "
            "red-step scaffolding, prose claims against code and the mutant kill record of "
            "invariant-critical suites, returning APPROVED or CHANGES REQUESTED. Read-only by "
            "design. Use before pushing anything."
        ),
        (
            "description: Final quality gatekeeper. Judges whether the domain model is valid, "
            "edge cases are handled and the change's claims are true, plus the procedural "
            "checks no tool makes yet. Read-only by design. Use before pushing anything."
        ),
        id="R2",
    ),
    pytest.param(
        ".claude/agents/reviewer.md",
        (
            "Your mandate is the checks no gate makes: CI status, the size ceiling, leftover "
            "red-step scaffolding, prose claims that the change made false, and the mutant "
            "kill record of a suite named invariant-critical. Readability, typing and style "
            "are the gates' (`CLAUDE.md` §4); do not spend the round on them."
        ),
        (
            "Your mandate is semantics: is the domain model valid (types, invariants, units, "
            "degenerate inputs), are edge cases handled, are the change's claims true. Tools "
            "enforce CI status, style (the required `CI result` check) and your verdict's "
            "record (a hook reads your `Review:` line). The checks below are yours until a "
            "tool makes them."
        ),
        id="R3",
    ),
    pytest.param(
        ".claude/agents/reviewer.md",
        (
            "### Precondition: CI status Before any verdict, check what CI says (`CLAUDE.md` "
            "§4: CI is authoritative): ```bash gh pr checks <pr> # or: gh run list --branch "
            "<branch> --limit 1 ``` **Red CI is an automatic `CHANGES REQUESTED`**, and so is "
            "a workflow that does not exercise the current build: read `.github/` when the "
            "change touches the build. On a prose-only PR (`CLAUDE.md` §4) the C++, sanitizer "
            "and Python checks show as skipped; that is the design, not a red, and `CI "
            "result` is the check that decides. ### The checks the gates cannot make"
        ),
        ("### The checks no tool makes yet"),
        id="R4",
    ),
    pytest.param(
        ".claude/agents/reviewer.md",
        (
            "3. **Actual LOC is reconciled against the increment doc's estimate.** Measure "
            "it; if the design named a split seam for an overrun, check whether it should "
            "fire."
        ),
        ("3. **A change to the build is exercised by the workflow** (read `.github/`)."),
        id="R5",
    ),
    pytest.param(
        ".claude/agents/reviewer.md",
        (
            "4. **Suggestions:** Non-blocking, and only where no gate would catch it. One the "
            "branch does not take goes into the private ideas file, added by your spawner."
        ),
        (
            "4. **Follow-ups:** non-blocking findings; your spawner files each as a pipe item "
            "(`tools/work_pipes.py add`)."
        ),
        id="R6",
    ),
    pytest.param(
        ".claude/agents/orchestrator.md",
        (
            "description: Workflow watcher and retrospective owner. Checks how the agents "
            "worked and flags deviations from the rules, researches agentic design and "
            "proposes improvements with evidence, and records lessons in "
            "docs/retrospectives/, the only place it writes. Use after an increment merges, "
            "the morning after an unattended night, or for a research round."
        ),
        (
            "description: Workflow watcher and retrospective owner. Checks how the agents "
            "worked, appends deviations and lessons to the lessons log, and once a week files "
            "them with a research round as one retrospective. Use when tools/duties.py shows "
            "a duty owed."
        ),
        id="R7",
    ),
    pytest.param(
        ".claude/agents/orchestrator.md",
        (
            "it would cost to adopt. Read the source before you cite it; mark anything from "
            "memory as unchecked. A research round's retrospective file name contains "
            "`-research`."
        ),
        (
            "it would cost to adopt. Read the source before you cite it; mark anything from "
            "memory as unchecked."
        ),
        id="R8",
    ),
    pytest.param(
        ".claude/agents/orchestrator.md",
        (
            "## 3. Own the retrospectives `docs/retrospectives/`, `next.md` included, is "
            "yours. Record the lessons the main session passes you (`CLAUDE.md` §3, *The main "
            "session dispatches*). Quotes, dates and incidents go here. Each retrospective "
            "measures the rule text (`python3 tools/rule_sizes.py`) and proposes a cut, "
            "restatements included, and flags anecdotes in rule files (a case, incident or "
            "example), moved to their owner or out. Once a week, check the private ideas file "
            "for entries gone stale or done unmarked."
        ),
        (
            "## 3. The lessons log and the weekly retrospective After a merge or a night, "
            "append to `../rasputin_scratch/lessons.md` (never committed) `## <stamp> merge "
            "#N` (or `night`), stamped at the start by `date -u +%Y-%m-%dT%H:%MZ`, then at "
            "most five lines of findings (§1) and lessons passed to you. No PR. A finding "
            "needing a fix goes in your handback too, for a pipe item. Weekly, one PR: a "
            "dated file in `docs/retrospectives/` (yours, `next.md` included) grouping the "
            "log by theme, with the research round (§2), no inflow, KGE or pcorr result, the "
            "rule text measured (`tools/rule_sizes.py`), a cut proposed and anecdotes in rule "
            "files flagged. Check the private ideas file for stale or done entries."
        ),
        id="R9",
    ),
    pytest.param(
        ".claude/agents/orchestrator.md",
        ("- **Write only under `docs/retrospectives/`.** New ideas, and stale or"),
        ("- **Write only under `docs/retrospectives/` and the log.** New ideas, and stale or"),
        id="R10",
    ),
    pytest.param(
        "docs/increments/README.md",
        (
            "An increment estimated under 50 net lines (`CLAUDE.md` §2) that touches no C++ "
            "file and nothing the *Acceptance* section covers takes the short path;"
        ),
        (
            "An increment estimated under 50 net lines (`CLAUDE.md` §2) that touches no C++ "
            "file, rule text or *Acceptance* code, nor adds a module or public type, takes "
            "the short path;"
        ),
        id="R11",
    ),
    pytest.param(
        "docs/increments/README.md",
        ("still goes to `@architect` before green."),
        (
            "still goes to `@architect` before green. - `@reviewer` checks these conditions "
            "on the real diff (`count_loc.py`); one failing means the full path. Else it "
            "lands after that round: only a wrong result, a broken invariant or a failing "
            "test blocks, and the next round checks that fix alone; other findings are pipe "
            "items."
        ),
        id="R12",
    ),
    pytest.param(
        "docs/increments/README.md",
        ("any other finding, is `CHANGES REQUESTED`."),
        ("any other finding, is `CHANGES REQUESTED` (short path: below)."),
        id="R13",
    ),
]


@pytest.mark.xfail(strict=True, reason=PENDING)
@pytest.mark.parametrize(("relative", "old", "new"), RULE_TEXTS)
def test_the_rule_text_reads_as_the_design_says(relative: str, old: str, new: str) -> None:
    text = collapsed((REAL / relative).read_text())
    assert text.count(new) == 1, f"{relative}: the new text is not there exactly once"
    if old not in new:
        assert old not in text, f"{relative}: the old text is still there"


def test_each_old_text_of_the_design_is_in_its_file_once_or_already_replaced() -> None:
    """The pairs are copied right: each old or new text is found in its file.

    Before R1 to R13 are applied every old text occurs once (the design's own
    claim, at master `662022b0`); after, every new one does. A pair found
    neither way was copied wrong, and the expected failure above would hide it.
    """
    for case in RULE_TEXTS:
        relative, old, new = case.values
        text = collapsed((REAL / str(relative)).read_text())
        assert text.count(str(old)) == 1 or text.count(str(new)) == 1, case.id
