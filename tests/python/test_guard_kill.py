"""`.claude/hooks/guard_kill.py`: a kill that can hit the wrong process is refused.

On 2026-11-17 a hand-written `kill <pid>` hit the session's own Claude Code
process; on 2026-11-07 a `pkill -f` hit Ola's applications. The hook refuses a
`kill` of its own ancestor, of a `claude` process, of a process group, and every
`pkill` and `killall`; a `kill` of any other pid, or of one it cannot read, passes.
The hook runs from the checkout by path, so its ancestors include this test.
"""

from __future__ import annotations

import os
import subprocess
from collections.abc import Iterator
from pathlib import Path

import pytest

from harness_fixtures import REAL, bash_event, pretool_decision, run_script

pytestmark = pytest.mark.harness

HOOK = ".claude/hooks/guard_kill.py"


def decide(command: str) -> tuple[str, str] | None:
    return pretool_decision(run_script(REAL, HOOK, bash_event(REAL, command)))


@pytest.fixture
def child() -> Iterator[int]:
    """A harmless process this test started."""
    with subprocess.Popen(["sleep", "60"]) as proc:
        yield proc.pid
        proc.kill()


@pytest.fixture
def claude(tmp_path: Path) -> Iterator[int]:
    """A process whose command is `claude` (a link to sleep), not an ancestor of the hook."""
    link = tmp_path / "claude"
    link.symlink_to("/bin/sleep")
    with subprocess.Popen([str(link), "60"]) as proc:
        yield proc.pid
        proc.kill()


@pytest.mark.parametrize(
    "form",
    ["kill {}", "kill -9 {}", "kill -TERM {}", "/bin/kill -s KILL {}", "true && kill {}",
     "echo $(kill -9 {})", "echo {} | xargs kill", "sh -c 'kill {}'"],
)  # fmt: skip
def test_own_ancestor_is_refused(form: str) -> None:
    verdict = decide(form.format(os.getpid()))
    assert verdict is not None and verdict[0] == "deny"
    assert str(os.getpid()) in verdict[1] and "ancestor" in verdict[1]


def test_claude_process_is_refused(claude: int) -> None:
    verdict = decide(f"kill -9 {claude}")
    assert verdict is not None and verdict[0] == "deny"
    assert str(claude) in verdict[1] and "claude" in verdict[1]


@pytest.mark.parametrize(
    "command",
    ["pkill -f sleep", "killall sleep", "/usr/bin/pkill x", "ls; nohup killall -9 node",
     "sh -c 'pkill x'"],
)  # fmt: skip
def test_pkill_and_killall_are_refused(command: str) -> None:
    verdict = decide(command)
    assert verdict is not None and verdict[0] == "deny"
    assert "tools/reap.py --kill" in verdict[1]


@pytest.mark.parametrize("command", ["kill 0", "kill -- -1", "kill -9 -1", "kill -- -{}"])
def test_process_groups_are_refused(command: str, child: int) -> None:
    verdict = decide(command.format(child))
    assert verdict is not None and verdict[0] == "deny"


@pytest.mark.parametrize("form", ["kill {}", "kill -9 {}", "kill -s TERM {}", "kill %1"])
def test_harmless_child_is_allowed(form: str, child: int) -> None:
    assert decide(form.format(child)) is None


@pytest.mark.parametrize(
    "command",
    ["grep kill file", "echo skill", "echo 'pkill -f x'", "git log --grep=killall",
     "kill $PID", "kill $(cat x.pid)", "kill -l"],
)  # fmt: skip
def test_unrelated_or_unreadable_is_allowed(command: str) -> None:
    assert decide(command) is None
