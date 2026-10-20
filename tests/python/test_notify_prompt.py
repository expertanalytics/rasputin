"""`.claude/hooks/notify_prompt.py`: a Mac banner for a waiting permission prompt (h28).

The spec is `docs/increments/h28-away-prompt.md` §3 and §6, tests 12 to 14.
`notify(event, run, platform)` is called with a recording `run`, so no banner
is ever shown; the one subprocess run of the script feeds it stdin that is not
a notification event, so it never reaches `osascript` either.

The hook checks that `/usr/bin/osascript` exists and offers no seam for that
check, so the test that expects a call runs only where the file exists.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from harness_fixtures import REAL, clean_env

pytestmark = pytest.mark.harness

SCRIPT = REAL / ".claude" / "hooks" / "notify_prompt.py"
OSASCRIPT = Path("/usr/bin/osascript")
MESSAGE = 'say "hi"'
WORKTREE = "h28-away-prompt"


def hook() -> ModuleType:
    """The hook module, loaded in the test body so its absence fails the test, not collection."""
    if not SCRIPT.exists():
        pytest.fail(".claude/hooks/notify_prompt.py is missing from the checkout")
    spec = importlib.util.spec_from_file_location("h28_notify_prompt", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Run:
    """Stands in for `subprocess.run`: records each call, raises if told to."""

    def __init__(self, raises: BaseException | None = None) -> None:
        self.calls: list[tuple[list[str], dict[str, Any]]] = []
        self._raises = raises

    def __call__(self, argv: list[str], **kwargs: Any) -> subprocess.CompletedProcess[str]:
        self.calls.append((list(argv), kwargs))
        if self._raises is not None:
            raise self._raises
        return subprocess.CompletedProcess(argv, 0, "", "")


def notification(kind: str = "permission_prompt") -> dict[str, Any]:
    return {
        "hook_event_name": "Notification",
        "notification_type": kind,
        "message": MESSAGE,
        "title": "Permission needed",
        "session_id": "s1",
        "transcript_path": "/x/s1.jsonl",
        "cwd": f"/x/{WORKTREE}",
    }


@pytest.mark.skipif(not OSASCRIPT.exists(), reason="/usr/bin/osascript exists only on a Mac")
def test_a_permission_prompt_runs_osascript_once_with_the_text_as_arguments() -> None:
    run = Run()
    argv = hook().notify(notification(), run=run, platform="darwin")
    [(called, kwargs)] = run.calls
    assert argv == called
    assert called[0] == str(OSASCRIPT)
    assert called[-2:] == [MESSAGE, WORKTREE]
    scripts = [called[i + 1] for i, part in enumerate(called[:-1]) if part == "-e"]
    assert scripts, called
    assert all(MESSAGE not in script and "hi" not in script for script in scripts), scripts
    assert kwargs.get("timeout") == 5, kwargs


@pytest.mark.parametrize(
    ("kind", "platform"),
    [("idle_prompt", "darwin"), ("permission_prompt", "linux")],
    ids=["idle-prompt", "not-a-mac"],
)
def test_no_banner_for_another_notification_or_another_platform(kind: str, platform: str) -> None:
    run = Run()
    assert hook().notify(notification(kind), run=run, platform=platform) is None
    assert run.calls == []


@pytest.mark.parametrize(
    "error",
    [OSError("no osascript"), subprocess.TimeoutExpired(str(OSASCRIPT), 5)],
    ids=["os-error", "timeout"],
)
def test_a_failing_osascript_is_swallowed(error: BaseException) -> None:
    run = Run(raises=error)
    hook().notify(notification(), run=run, platform="darwin")
    if OSASCRIPT.exists():  # otherwise the hook stops before running anything
        assert len(run.calls) == 1, run.calls


@pytest.mark.parametrize("stdin", ["{not json", "[1, 2]", ""], ids=["not-json", "list", "empty"])
def test_the_script_exits_0_on_stdin_that_is_not_an_event(stdin: str) -> None:
    if not SCRIPT.exists():
        pytest.fail(".claude/hooks/notify_prompt.py is missing from the checkout")
    result = subprocess.run(
        [sys.executable, str(SCRIPT)],
        input=stdin,
        capture_output=True,
        text=True,
        env=clean_env(),
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
