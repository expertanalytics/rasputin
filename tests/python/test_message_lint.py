"""`.claude/hooks/message_lint.py`: a main-session stop whose message is not self-contained.

The rules are checked on the module's `findings`; the hook itself is run from a
copy in a temporary repository (`harness_fixtures.make_repo`), so its hash
ledger lands in `<tmp>/.git/harness` and the real one is never touched.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from harness_fixtures import REAL, SUBAGENT, make_repo, run_script, state_dir

pytestmark = pytest.mark.harness

HOOK = ".claude/hooks/message_lint.py"


def load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("message_lint", REAL / HOOK)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lint = load()


def rules(text: str) -> list[tuple[str, str]]:
    found: list[tuple[str, str]] = lint.findings(text)
    return found


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("See tools/brief.py now.", [("tools/brief.py", "path")]),
        ("Read `docs/increments/h3.md`.", [("docs/increments/h3.md", "path")]),
        ("Edit .claude/hooks/x.py.", [(".claude/hooks/x.py", "path")]),
        ("The bug is at brief.py:12.", [("brief.py:12", "path")]),
        ("Open increment 150 next.", [("increment 150", "number")]),
        ("Design 180 says so.", [("Design 180", "number")]),
        ("See item 2.58.", [("item 2.58", "number")]),
        ("Answer q421 first.", [("q421", "number")]),
        ("Rerun the 150 tool.", [("the 150 tool", "number")]),
        ("Run A14 won.", [("A14", "label")]),
        ("HC15 and F91 lost.", [("HC15", "label"), ("F91", "label")]),
        ("P beat PL and L.", [("P", "label"), ("PL", "label"), ("L", "label")]),
    ],
)
def test_each_rule_flags(text: str, expected: list[tuple[str, str]]) -> None:
    assert rules(text) == expected


@pytest.mark.parametrize(
    "text",
    [
        "See /Users/ola/projects/rasputin/tools/brief.py and /private/tmp/x/docs/a.md.",
        "At /Users/ola/rasputin/tools/brief.py:12 and https://github.com/o/r/tools/a.py.",
        "The mesh refinement tool (increment 150) and the boundary question (q421).",
        "The cold-start run (A14) and the [F91] probe.",
        "P95 latency, an L2 norm, C++20, Sentinel-2, CO2, an AR50 map, Q3.",
        "Pressure P = rho g h, and a joke :P",
        "Use `A14`, `q421` and `increment 150` as typed.",
        "```\ncd /x && pytest tests/python/t.py -k A14 q421\n```\nDone.",
        "The 2016 run held.",
        "",
    ],
)
def test_clean_text_passes(text: str) -> None:
    assert rules(text) == []


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    root = make_repo(tmp_path / "repo")
    shutil.copy2(REAL / HOOK, root / HOOK)
    return root


def transcript(path: Path, *texts: str) -> Path:
    """A JSONL transcript: a user line, then one assistant message split over `texts`."""
    lines: list[dict[str, Any]] = [
        {"type": "user", "message": {"role": "user", "content": "hi"}},
        {"type": "assistant", "message": {"id": "m0", "content": [{"type": "text", "text": "A1"}]}},
        {"type": "user", "message": {"role": "user", "content": [{"type": "tool_result"}]}},
    ]
    for text in texts:
        block = {"type": "text", "text": text}
        lines.append({"type": "assistant", "message": {"id": "m1", "content": [block]}})
    lines.append({"type": "system"})
    path.write_text("not json\n" + "\n".join(json.dumps(line) for line in lines) + "\n")
    return path


def stop(repo: Path, text: str, **extra: Any) -> dict[str, Any]:
    path = transcript(repo / "t.jsonl", "Intro.", text)
    return {
        "hook_event_name": "Stop",
        "session_id": "s1",
        "transcript_path": str(path),
        "stop_hook_active": False,
        **extra,
    }


def decision(result: Any) -> dict[str, Any] | None:
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout) if result.stdout.strip() else None


def test_violation_blocks_naming_fragment_and_rule(repo: Path) -> None:
    out = decision(run_script(repo, HOOK, stop(repo, "See tools/brief.py and q421.")))
    assert out is not None and out["decision"] == "block"
    assert "`tools/brief.py`: a relative path" in out["reason"]
    assert "`q421`: work named by a bare number" in out["reason"]


def test_only_the_last_message_counts(repo: Path) -> None:
    """The earlier message `m0` holds a label; the last one is clean."""
    assert decision(run_script(repo, HOOK, stop(repo, "All clean."))) is None


def test_last_assistant_message_field_wins(repo: Path) -> None:
    event = stop(repo, "Clean.", last_assistant_message="Run A14 won.")
    out = decision(run_script(repo, HOOK, event))
    assert out is not None and "`A14`" in out["reason"]


def test_same_message_is_blocked_once(repo: Path) -> None:
    event = stop(repo, "Run A14 won.")
    assert decision(run_script(repo, HOOK, event)) is not None
    assert decision(run_script(repo, HOOK, {**event, "stop_hook_active": True})) is None
    assert (state_dir(repo) / "message_lint.jsonl").exists()


def test_a_streak_of_blocks_ends(repo: Path) -> None:
    for n in range(lint.STREAK):
        assert decision(run_script(repo, HOOK, stop(repo, f"Run A{n + 10} won."))) is not None
    event = stop(repo, "Run A99 won.", stop_hook_active=True)
    assert decision(run_script(repo, HOOK, event)) is None


def test_subagent_passes_untouched(repo: Path) -> None:
    event = stop(repo, "See tools/brief.py, q421 and A14.", **SUBAGENT)
    assert decision(run_script(repo, HOOK, event)) is None
    assert not (state_dir(repo) / "message_lint.jsonl").exists()


@pytest.mark.parametrize(
    "stdin",
    [
        "not json",
        {"hook_event_name": "Stop", "transcript_path": "/nonexistent/t.jsonl"},
        {"hook_event_name": "Stop", "transcript_path": 7, "last_assistant_message": 7},
    ],
)
def test_crash_or_nothing_to_read_lets_the_stop_through(repo: Path, stdin: Any) -> None:
    result = run_script(repo, HOOK, stdin)
    assert decision(result) is None


def test_crash_reports_on_stderr(repo: Path) -> None:
    event = {"hook_event_name": "Stop", "transcript_path": "/nonexistent/t.jsonl"}
    assert "message_lint failed" in run_script(repo, HOOK, event).stderr


def test_hook_is_executable() -> None:
    assert (REAL / HOOK).stat().st_mode & 0o111
