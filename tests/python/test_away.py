"""`tools/away.py`: Ola enters and leaves unattended mode (h3 U1, T12-T17).

The spec is `docs/increments/h3-unattended-u1.md` §3.2, §3.3 and §4. Entering
needs a real terminal and a typed `y`; leaving needs neither. The in-process
tests drive `main(argv, *, root, now, open_tty, spawn, stop)` with fakes, and
`root` is always a temporary repository, so the real repository's harness
state is never touched and no real keep-awake process is started.

Seams, as this suite pins them (§3.3 names them, not their shapes):
`open_tty()` returns a text-mode file object open for reading and writing, as
`open("/dev/tty", "r+")` would; `spawn(argv, **kwargs)` returns an object with
a `pid`, as `subprocess.Popen` does; `stop(pid)` is called with the previous
flag's `keep_awake_pid`. The default stopper is `stop_keep_awake(pid)`.
"""

from __future__ import annotations

import errno
import json
import shutil
import signal
import subprocess
import sys
from collections.abc import Callable
from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from harness_fixtures import (
    Tool,
    clean_env,
    flag_on,
    flag_path,
    iso,
    make_repo,
    state_dir,
)

NOW = datetime(2026, 9, 30, 18, 0, 0, tzinfo=UTC)
FAKE_PID = 4242


away = Tool("away")


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    return make_repo(tmp_path / "repo")


class FakeTty:
    """A terminal that answers one line and records everything written to it."""

    def __init__(self, answer: str) -> None:
        self._answer = answer + "\n"
        self.written = ""
        self.closed = False

    def write(self, text: str) -> int:
        self.written += text
        return len(text)

    def flush(self) -> None:
        pass

    def readline(self, size: int = -1) -> str:
        answer, self._answer = self._answer, ""
        return answer

    def read(self, size: int = -1) -> str:
        return self.readline()

    def close(self) -> None:
        self.closed = True

    def __enter__(self) -> FakeTty:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


class Recorder:
    """Stands in for `spawn` or `stop`; records calls, returns or raises as told."""

    def __init__(self, result: Any = None, raises: BaseException | None = None) -> None:
        self.calls: list[tuple[tuple[Any, ...], dict[str, Any]]] = []
        self._result = result
        self._raises = raises

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        self.calls.append((args, kwargs))
        if self._raises is not None:
            raise self._raises
        return self._result


def no_tty() -> Any:
    raise OSError(errno.ENXIO, "Device not configured")


def tty_answering(answer: str) -> tuple[FakeTty, Callable[[], FakeTty]]:
    tty = FakeTty(answer)
    return tty, lambda: tty


def spawner() -> Recorder:
    return Recorder(result=SimpleNamespace(pid=FAKE_PID))


def windows(repo: Path) -> list[dict[str, Any]]:
    path = state_dir(repo) / "windows.jsonl"
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def read_flag(repo: Path) -> dict[str, Any]:
    loaded: dict[str, Any] = json.loads(flag_path(repo).read_text())
    return loaded


def run_main(argv: list[str], **seams: Any) -> int:
    """`away.main`'s exit status, whether it returns it or raises SystemExit.

    §3.3 pins the status, not the mechanism: an argparse error on `-1h`
    raises SystemExit(2) from inside `main`.
    """
    try:
        return int(away.main(argv, **seams))
    except SystemExit as exit_:
        return exit_.code if isinstance(exit_.code, int) else 1


def enter(
    repo: Path,
    argv: list[str],
    answer: str = "y",
    spawn: Recorder | None = None,
    stop: Recorder | None = None,
) -> tuple[int, FakeTty, Recorder, Recorder]:
    tty, open_tty = tty_answering(answer)
    spawn = spawn or spawner()
    stop = stop or Recorder()
    code = run_main(argv, root=repo, now=NOW, open_tty=open_tty, spawn=spawn, stop=stop)
    return code, tty, spawn, stop


# ---------------------------------------------------------------- T12


def test_away_without_a_terminal_exits_3_and_writes_nothing(repo: Path) -> None:
    script = repo / "tools" / "away.py"
    assert script.exists(), "tools/away.py does not exist yet"
    result = subprocess.run(
        [sys.executable, str(script), "8h"],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        cwd=repo,
        env=clean_env(),
        start_new_session=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 3
    assert "needs your own terminal" in result.stderr
    assert not flag_path(repo).exists()
    assert not (state_dir(repo) / "windows.jsonl").exists()


def test_away_without_a_terminal_starts_nothing(repo: Path) -> None:
    spawn = spawner()
    code = run_main(["8h"], root=repo, now=NOW, open_tty=no_tty, spawn=spawn, stop=Recorder())
    assert code == 3
    assert spawn.calls == []
    assert not flag_path(repo).exists()


# ---------------------------------------------------------------- T13


@pytest.mark.parametrize(
    ("duration", "minutes"), [("8h", 480), ("90m", 90), ("1h30m", 90), ("72h", 4320)]
)
def test_a_valid_duration_sets_a_window_of_it_times_the_buffer(
    repo: Path, duration: str, minutes: int
) -> None:
    code, _, spawn, _ = enter(repo, [duration])
    assert code == 0
    flag = read_flag(repo)
    since = datetime.fromisoformat(flag["since"])
    until = datetime.fromisoformat(flag["until"])
    assert since == NOW
    assert until - since == timedelta(seconds=minutes * 72)  # minutes x 60 s x 1.2
    [(args, _)] = spawn.calls
    assert args[0] == ["caffeinate", "-is", "-t", str(minutes * 72)]


@pytest.mark.parametrize("duration", ["72h1m", "0h", "0m", "", "8", "8x", "-1h", "1.5h"])
def test_an_invalid_duration_exits_2_before_opening_the_terminal(repo: Path, duration: str) -> None:
    opened = Recorder(raises=AssertionError("the terminal was opened"))
    spawn = spawner()
    code = run_main([duration], root=repo, now=NOW, open_tty=opened, spawn=spawn, stop=Recorder())
    assert code == 2
    assert opened.calls == []
    assert spawn.calls == []
    assert not flag_path(repo).exists()


# ---------------------------------------------------------------- T14


def test_entering_for_8h_writes_the_flag_and_starts_keep_awake(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code, tty, spawn, _ = enter(repo, ["8h"])
    assert code == 0
    flag = read_flag(repo)
    since = datetime.fromisoformat(flag["since"])
    until = datetime.fromisoformat(flag["until"])
    assert until - since == timedelta(hours=9, minutes=36)
    assert since.utcoffset() == timedelta(0)
    assert flag["keep_awake_pid"] == FAKE_PID
    assert flag_path(repo).stat().st_mode & 0o777 == 0o600
    assert state_dir(repo).stat().st_mode & 0o777 == 0o700
    [(args, kwargs)] = spawn.calls
    assert args[0] == ["caffeinate", "-is", "-t", "34560"]
    assert kwargs.get("start_new_session") is True
    assert windows(repo) == [{"event": "enter", "since": flag["since"], "until": flag["until"]}]
    assert "Type y to confirm" in tty.written
    assert "UNATTENDED until" in capsys.readouterr().out


@pytest.mark.parametrize("answer", ["n", "", "maybe"])
def test_anything_but_yes_changes_nothing(
    repo: Path, answer: str, capsys: pytest.CaptureFixture[str]
) -> None:
    code, tty, spawn, _ = enter(repo, ["8h"], answer=answer)
    assert code == 1
    assert not flag_path(repo).exists()
    assert spawn.calls == []
    assert windows(repo) == []
    captured = capsys.readouterr()
    assert "Not confirmed; nothing changed." in captured.out + captured.err + tty.written


@pytest.mark.parametrize("answer", ["YES", " y ", "Yes"])
def test_yes_is_read_case_insensitively_and_stripped(repo: Path, answer: str) -> None:
    code, _, _, _ = enter(repo, ["8h"], answer=answer)
    assert code == 0
    assert flag_path(repo).exists()


# ---------------------------------------------------------------- T15


def test_a_missing_keep_awake_still_enters_and_says_so(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    missing = Recorder(raises=FileNotFoundError(2, "No such file", "caffeinate"))
    code, _, _, _ = enter(repo, ["8h"], spawn=missing)
    assert code == 0
    assert read_flag(repo)["keep_awake_pid"] is None
    assert "not started" in capsys.readouterr().out


def test_re_entering_while_on_stops_the_old_keep_awake_and_says_it_replaces(repo: Path) -> None:
    flag_on(repo, now=NOW - timedelta(minutes=10), keep_awake_pid=111)
    code, tty, _, stop = enter(repo, ["2h"])
    assert code == 0
    assert stop.calls == [((111,), {})]
    assert "replaces the window" in tty.written
    assert read_flag(repo)["keep_awake_pid"] == FAKE_PID


# ---------------------------------------------------------------- T16


def test_the_stopper_never_signals_a_process_that_is_not_caffeinate() -> None:
    child = subprocess.Popen(["sleep", "30"])
    try:
        away.stop_keep_awake(child.pid)
        with pytest.raises(subprocess.TimeoutExpired):
            child.wait(timeout=0.5)
    finally:
        child.kill()
        child.wait()


@pytest.mark.parametrize("pid", [None, "123", 12.0, [123]])
def test_the_stopper_ignores_a_pid_that_is_not_an_int(pid: object) -> None:
    away.stop_keep_awake(pid)


@pytest.mark.skipif(shutil.which("caffeinate") is None, reason="caffeinate is macOS only")
def test_the_stopper_terminates_a_caffeinate_it_started() -> None:
    child = subprocess.Popen(["caffeinate", "-t", "30"])
    try:
        away.stop_keep_awake(child.pid)
        assert child.wait(timeout=10) == -signal.SIGTERM
    finally:
        if child.poll() is None:
            child.kill()
            child.wait()


# ---------------------------------------------------------------- T17

BACK_AT = NOW + timedelta(hours=5)

QUEUED = (
    {"at": "2026-09-30T19:00:00+00:00", "branch": "feat-a", "hook": "guard_push",
     "agent_type": None, "agent_id": None, "act": "git push", "why": "w", "cwd": "/r"},
    {"at": "2026-09-30T20:00:00+00:00", "branch": "feat-b", "hook": "guard_governance",
     "agent_type": "tester", "agent_id": "a1", "act": "Edit /r/CLAUDE.md", "why": "w",
     "cwd": "/r"},
    {"at": "2026-09-30T21:00:00+00:00", "branch": "feat-a", "hook": "guard_push",
     "agent_type": "developer", "agent_id": "a2", "act": "gh pr merge 1", "why": "w",
     "cwd": "/r"},
)  # fmt: skip


def write_queue(repo: Path, entries: tuple[dict[str, Any], ...]) -> None:
    state_dir(repo).mkdir(parents=True, exist_ok=True)
    with (state_dir(repo) / "queue.jsonl").open("a") as queue:
        for entry in entries:
            queue.write(json.dumps(entry) + "\n")


def back(repo: Path, stop: Recorder | None = None) -> tuple[int, Recorder]:
    stop = stop or Recorder()
    code = run_main(["--back"], root=repo, now=BACK_AT, open_tty=no_tty, spawn=spawner(), stop=stop)
    return code, stop


def archive(repo: Path) -> Path:
    return state_dir(repo) / "queue-2026-09-30.jsonl"


def test_back_prints_the_queue_by_branch_and_archives_it(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    flag_on(repo, now=NOW, keep_awake_pid=777)
    until = read_flag(repo)["until"]
    write_queue(repo, QUEUED)
    tasks = repo / ".claude" / "current-task"
    tasks.mkdir(parents=True)
    (tasks / "session.md").write_text("round: U1\nASK OLA: may I push?\n")

    code, stop = back(repo)

    assert code == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith(f"Back. Unattended since {iso(NOW)} until {until}")
    at = lines.index("Queued while away (3):")
    assert lines[at + 1 : at + 7] == [
        "  feat-a",
        "    2026-09-30T19:00:00+00:00 guard_push main: git push",
        "    2026-09-30T21:00:00+00:00 guard_push developer: gh pr merge 1",
        "  feat-b",
        "    2026-09-30T20:00:00+00:00 guard_governance tester: Edit /r/CLAUDE.md",
        "ASK OLA lines in .claude/current-task/:",
    ]
    assert "  session.md: ASK OLA: may I push?" in lines
    assert not flag_path(repo).exists()
    assert not (state_dir(repo) / "queue.jsonl").exists()
    assert len(archive(repo).read_text().splitlines()) == 3
    [window] = windows(repo)
    assert window["event"] == "back"
    assert datetime.fromisoformat(window["until"]) == BACK_AT
    assert stop.calls == [((777,), {})]


def test_a_second_back_appends_to_the_same_archive(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    flag_on(repo, now=NOW, keep_awake_pid=777)
    write_queue(repo, QUEUED)
    back(repo)
    write_queue(repo, QUEUED[:1])
    capsys.readouterr()

    code, _ = back(repo)

    assert code == 0
    assert "No unattended flag was set." in capsys.readouterr().out
    assert len(archive(repo).read_text().splitlines()) == 4
    assert not (state_dir(repo) / "queue.jsonl").exists()
    assert [w["event"] for w in windows(repo)] == ["back"]


def test_back_with_no_flag_and_no_queue_says_so(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    code, stop = back(repo)
    assert code == 0
    out = capsys.readouterr().out
    assert out.splitlines()[0] == "No unattended flag was set."
    assert "(none)" in out
    assert stop.calls == []
    assert windows(repo) == []


def test_back_with_a_broken_flag_still_deletes_it(repo: Path) -> None:
    flag_path(repo).parent.mkdir(parents=True)
    flag_path(repo).write_text("{not json")
    code, _ = back(repo)
    assert code == 0
    assert not flag_path(repo).exists()


def test_back_needs_no_terminal_even_as_a_detached_subprocess(repo: Path) -> None:
    flag_on(repo)
    script = repo / "tools" / "away.py"
    assert script.exists(), "tools/away.py does not exist yet"
    result = subprocess.run(
        [sys.executable, str(script), "--back"],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        cwd=repo,
        env=clean_env(),
        start_new_session=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert not flag_path(repo).exists()
