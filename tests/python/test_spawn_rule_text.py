"""The rule text h9 changes (test 26).

Spec: `docs/increments/h9-spawn-briefs.md` §3.3 and §4. Each check reads the
real file of the checkout. The persona-file lines are written against branch
`worktree-agent-lines` (8b8eb94), which is not merged; these checks pin only
h9's own phrases, not the text around them.
"""

from __future__ import annotations

import re

from harness_fixtures import REAL

CLAUDE = REAL / "CLAUDE.md"
REQUIRED_READING = REAL / ".claude" / "REQUIRED-READING.md"
TESTER = REAL / ".claude" / "agents" / "tester.md"
REVIEWER = REAL / ".claude" / "agents" / "reviewer.md"


def collapsed(text: str) -> str:
    return " ".join(text.split())


def bullet(text: str, opening: str) -> str:
    """The Markdown bullet that starts with `opening`, its wrapped lines joined."""
    match = re.search(rf"^\* {re.escape(opening)}.*?(?=^\* |^#|\Z)", text, re.MULTILINE | re.DOTALL)
    assert match is not None, f"no bullet starts {opening!r}"
    return collapsed(match.group(0))


def harness_section() -> str:
    text = REQUIRED_READING.read_text()
    section = re.search(r"^## The harness\n(.*?)(?=^## )", text, re.MULTILINE | re.DOTALL)
    assert section is not None, "REQUIRED-READING.md has no '## The harness' section"
    return section.group(1)


def test_claude_md_has_no_briefs_bullet() -> None:
    """§4: the "Briefs" bullet is deleted; the block and the persona files hold it."""
    lines = CLAUDE.read_text().splitlines()
    assert not [line for line in lines if line.startswith("* **Briefs.**")]


def test_claude_md_step_order_times_a_performance_fix_before_review() -> None:
    step_order = bullet(CLAUDE.read_text(), "**Step order.**")
    assert "A performance fix is timed by `@perf` before review." in step_order


def test_required_reading_names_the_hook_and_the_script_in_the_harness() -> None:
    section = harness_section()
    assert "guard_spawn.py" in section
    assert "tools/brief.py" in section


def test_required_reading_no_longer_says_the_brief_says_to_read_from_disk() -> None:
    assert "until then its brief says" not in collapsed(REQUIRED_READING.read_text())


def test_tester_md_has_the_choices_beyond_the_design_bullet() -> None:
    assert "Choices beyond the design" in TESTER.read_text()


def test_reviewer_md_reports_the_commit_range() -> None:
    assert "The commit range reviewed" in REVIEWER.read_text()
