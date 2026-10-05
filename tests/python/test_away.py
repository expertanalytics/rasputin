"""`tools/away.py`: Ola enters and leaves unattended mode (h3 U1, T12-T17, T22-T24).

The spec is `docs/increments/h3-unattended-u1.md` §3.2, §3.3 and §4. Entering
needs a real terminal and a typed `y`; leaving needs neither. The in-process
tests drive `main(argv, *, root, now, open_tty, spawn, stop)` with fakes, and
`root` is always a temporary repository, so the real repository's harness
state is never touched and no real keep-awake process is started.

Seams, as this suite pins them (§3.3 names them, not their shapes):
`open_tty()` returns a text-mode file object open for reading and writing
(the default is the non-seeking opener of §3.3 step 2, and T23 is the one
test that does not inject it); `spawn(argv, **kwargs)` returns an object with
a `pid`, as `subprocess.Popen` does; `stop(pid)` is called with the previous
flag's `keep_awake_pid`. The default stopper is `stop_keep_awake(pid)`.
"""

from __future__ import annotations

import errno
import io
import json
import os
import select
import shlex
import shutil
import signal
import subprocess
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from harness_fixtures import (
    Tool,
    add_worktree,
    clean_env,
    commit_at,
    flag_on,
    flag_path,
    flag_window,
    git,
    iso,
    make_repo,
    state_dir,
    write_flag,
)

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

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
    assert script.exists(), "tools/away.py is missing from the copy"
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


@pytest.mark.skipif(os.geteuid() == 0, reason="root writes into a 0o500 directory")
def test_a_flag_that_cannot_be_written_stops_the_keep_awake_and_changes_nothing(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # §3.3 step 6 (review 1): the flag write is one unit; its failure undoes step 5.
    state = state_dir(repo)
    state.mkdir(mode=0o700)
    state.chmod(0o500)
    try:
        code, _, spawn, stop = enter(repo, ["8h"])
    finally:
        state.chmod(0o700)
    assert code == 1
    assert len(spawn.calls) == 1
    assert stop.calls == [((FAKE_PID,), {})]
    captured = capsys.readouterr()
    assert "could not write the flag" in captured.err
    assert "keep-awake stopped" in captured.err
    assert "UNATTENDED until" not in captured.out
    assert not flag_path(repo).exists()
    assert not (state / "windows.jsonl").exists()
    assert list(state.iterdir()) == []  # no temporary file left behind


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
DECISIONS_HEADING = "ASK OLA lines, every worktree's .claude/current-task/:"

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
        "ASK OLA lines, every worktree's .claude/current-task/:",
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
    assert script.exists(), "tools/away.py is missing from the copy"
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


@pytest.mark.parametrize("branch", ["null", "missing"])
def test_back_lists_a_null_or_missing_branch_under_no_branch(
    repo: Path, capsys: pytest.CaptureFixture[str], branch: str
) -> None:
    # §3.3 (review 1): the same text the recap uses (§3.8).
    entry = {key: value for key, value in QUEUED[0].items() if key != "branch"}
    if branch == "null":
        entry["branch"] = None
    write_queue(repo, (entry,))
    code, _ = back(repo)
    assert code == 0
    lines = capsys.readouterr().out.splitlines()
    at = lines.index("Queued while away (1):")
    assert lines[at + 1 : at + 3] == [
        "  (no branch)",
        "    2026-09-30T19:00:00+00:00 guard_push main: git push",
    ]
    assert not any("None" in line for line in lines)


# ---------------------------------------------------------------- T22


def test_back_survives_a_broken_session_state(repo: Path) -> None:
    # §3.3 (review 1): the recap module is optional to --back; the flag, the
    # keep-awake and the queue are not.
    (repo / "tools" / "session_state.py").write_text('raise RuntimeError("planted")\n')
    flag_on(repo)
    write_queue(repo, QUEUED[:1])
    result = subprocess.run(
        [sys.executable, str(repo / "tools" / "away.py"), "--back"],
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
    lines = result.stdout.splitlines()
    at = lines.index("ASK OLA lines, every worktree's .claude/current-task/:")
    assert lines[at + 1] == "  (unavailable: RuntimeError: planted)"
    assert not flag_path(repo).exists()
    assert not (state_dir(repo) / "queue.jsonl").exists()
    [archived] = state_dir(repo).glob("queue-*.jsonl")
    assert len(archived.read_text().splitlines()) == 1
    assert [w["event"] for w in windows(repo)] == ["back"]


# ---------------------------------------------------------------- T23

PROMPT = b"Type y to confirm: "
PTY_DEADLINE = 10.0  # seconds, measured in the parent (§4, T23)

#: The child's program: the fixture's copy of `away.py`, imported by path and
#: run with no `open_tty`, so the default opener meets the pty as /dev/tty.
PTY_CHILD = """
import importlib.util, sys
from pathlib import Path
from types import SimpleNamespace
root = Path(sys.argv[1])
spec = importlib.util.spec_from_file_location("away_on_pty", root / "tools" / "away.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
spawn = lambda *args, **kwargs: SimpleNamespace(pid=4242)
sys.exit(module.main(["8h"], root=root, spawn=spawn, stop=lambda pid: None))
"""


@dataclass(frozen=True)
class PtyRun:
    status: int
    before_input: bytes  # what the master read before the answer was written
    output: bytes  # everything the master read
    answered: bool


def read_master(master: int, deadline: float) -> bytes | None:
    """One chunk from the pty master; b"" at end of file; None if nothing yet.

    The slave's last close reads as EOF on macOS and as EIO on Linux.
    """
    wait = max(0.0, min(0.1, deadline - time.monotonic()))
    ready, _, _ = select.select([master], [], [], wait)
    if not ready:
        return None
    try:
        return os.read(master, 4096)
    except OSError as error:
        if error.errno == errno.EIO:
            return b""
        raise


def run_on_pty(repo: Path, answer: bytes) -> PtyRun:
    """Run `PTY_CHILD` on a fresh pseudo-terminal, answering the prompt once.

    The child is killed and the test failed if it has not exited within
    `PTY_DEADLINE` seconds of the fork, so a read that never returns cannot
    hang the suite.
    """
    pty = pytest.importorskip("pty")
    script = repo / "tools" / "away.py"
    assert script.exists(), "tools/away.py is missing from the copy"
    pid, master = pty.fork()
    if pid == 0:  # pragma: no cover - the child; it never returns
        for key in [key for key in os.environ if key.startswith("GIT_")]:
            del os.environ[key]
        os.execv(sys.executable, [sys.executable, "-c", PTY_CHILD, str(repo)])
    deadline = time.monotonic() + PTY_DEADLINE
    output, before_input, answered = b"", b"", False
    try:
        while time.monotonic() < deadline:
            chunk = read_master(master, deadline)
            if chunk == b"":
                break
            output += chunk or b""
            if not answered and PROMPT in output:
                before_input, answered = output, True
                os.write(master, answer)
        while time.monotonic() < deadline:
            done, status = os.waitpid(pid, os.WNOHANG)
            if done:
                return PtyRun(
                    status=os.waitstatus_to_exitcode(status),
                    before_input=before_input or output,
                    output=output,
                    answered=answered,
                )
            time.sleep(0.05)
        os.kill(pid, signal.SIGKILL)
        os.waitpid(pid, 0)
        pytest.fail(f"away.py on a pty did not exit within {PTY_DEADLINE} s; read {output!r}")
    finally:
        os.close(master)


@pytest.mark.skipif(sys.platform == "win32", reason="pty is POSIX only")
def test_the_default_terminal_opener_confirms_on_a_real_pty(repo: Path) -> None:
    run = run_on_pty(repo, b"y\n")
    assert run.status == 0, f"exit {run.status}; the pty read {run.output!r}"
    assert run.answered and PROMPT in run.before_input
    assert b"needs your own terminal" not in run.output
    assert flag_path(repo).exists()
    assert [w["event"] for w in windows(repo)] == ["enter"]


@pytest.mark.skipif(sys.platform == "win32", reason="pty is POSIX only")
def test_the_default_terminal_opener_declines_on_a_real_pty(repo: Path) -> None:
    run = run_on_pty(repo, b"n\n")
    assert run.status == 1, f"exit {run.status}; the pty read {run.output!r}"
    assert run.answered and PROMPT in run.before_input
    assert b"needs your own terminal" not in run.output
    assert not flag_path(repo).exists()


# ---------------------------------------------------------------- T24


def failing_tty(error: OSError) -> Callable[[], Any]:
    def open_tty() -> Any:
        raise error

    return open_tty


def test_a_terminal_failure_without_errno_is_named_by_its_type(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # §3.3 step 2 (tty fix): errno None is a defect in away.py, not a missing
    # terminal, and the parenthesis must not pass it off as an OS condition.
    spawn = spawner()
    opener = failing_tty(io.UnsupportedOperation("x"))
    code = run_main(["8h"], root=repo, now=NOW, open_tty=opener, spawn=spawn, stop=Recorder())
    assert code == 3
    stderr = capsys.readouterr().err
    assert "needs your own terminal" in stderr
    assert "(UnsupportedOperation: x)" in stderr
    assert spawn.calls == []
    assert not flag_path(repo).exists()


def test_a_terminal_failure_with_errno_is_named_by_the_errno(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    spawn = spawner()
    code = run_main(["8h"], root=repo, now=NOW, open_tty=no_tty, spawn=spawn, stop=Recorder())
    assert code == 3
    assert "(ENXIO)" in capsys.readouterr().err
    assert spawn.calls == []


# ---------------------------------------------------------------- h8, test 6
#
# docs/increments/h8-window-and-recap.md §3.4: the longest interval between
# consecutive points of `since`, the commit times inside [since, end] sorted,
# and `end`. Pure; commits are (time, short hash).

SINCE = datetime(2026, 10, 2, 20, 0, tzinfo=UTC)
END = SINCE + timedelta(hours=9)


def at(hours: float) -> datetime:
    return SINCE + timedelta(hours=hours)


def quiet(start: datetime, end: datetime, opener: str | None, count: int) -> Any:
    return away.Quiet(start=start, end=end, opener=opener, count=count)


def test_quiet_is_a_frozen_dataclass() -> None:
    found = quiet(SINCE, END, None, 0)
    with pytest.raises(AttributeError):
        found.count = 1


def test_with_no_commits_the_whole_window_is_quiet() -> None:
    assert away.longest_quiet(SINCE, END, []) == quiet(SINCE, END, None, 0)


def test_the_longest_gap_follows_the_last_commit() -> None:
    commits = [(at(1), "aaa1111"), (at(2), "bbb2222")]
    assert away.longest_quiet(SINCE, END, commits) == quiet(at(2), END, "bbb2222", 2)


def test_commits_outside_the_window_are_ignored() -> None:
    commits = [(at(-1), "before0"), (at(1), "aaa1111"), (at(2), "bbb2222"), (at(10), "after00")]
    assert away.longest_quiet(SINCE, END, commits) == quiet(at(2), END, "bbb2222", 2)


def test_the_input_order_does_not_matter() -> None:
    commits = [(at(2), "bbb2222"), (at(8.5), "ccc3333"), (at(1), "aaa1111")]
    assert away.longest_quiet(SINCE, END, commits) == quiet(at(2), at(8.5), "bbb2222", 3)


def test_a_tie_takes_the_earliest_interval() -> None:
    commits = [(at(6), "bbb2222"), (at(3), "aaa1111")]
    assert away.longest_quiet(SINCE, END, commits) == quiet(SINCE, at(3), None, 2)


def test_the_gap_from_the_window_start_wins_when_longest() -> None:
    assert away.longest_quiet(SINCE, END, [(at(8), "aaa1111")]) == quiet(SINCE, at(8), None, 1)


def test_the_window_is_closed_at_both_ends() -> None:
    # Points: since, c1 (= since), c2 (= end), end. The 9 h interval opens at c1.
    commits = [(END, "c2c2c2c"), (SINCE, "c1c1c1c")]
    assert away.longest_quiet(SINCE, END, commits) == quiet(SINCE, END, "c1c1c1c", 2)


# ---------------------------------------------------------------- h8, test 7
#
# §3.4: `--back` prints the longest stretch without a commit on any ref, right
# after the `Back.` line. The window is [since, min(now, until)].

LONGEST = "Longest stretch without a commit (any ref): "


def hhmm(moment: datetime) -> str:
    return f"{moment.astimezone(UTC):%Y-%m-%d %H:%M}"


def back_lines(repo: Path, capsys: pytest.CaptureFixture[str], now: datetime) -> list[str]:
    capsys.readouterr()
    code = run_main(
        ["--back"], root=repo, now=now, open_tty=no_tty, spawn=spawner(), stop=Recorder()
    )
    assert code == 0
    return capsys.readouterr().out.splitlines()


def commit_on_branches(repo: Path, since: datetime) -> tuple[str, str]:
    """One commit on master at since + 10 min, one on `side` at since + 40 min 30 s."""
    first = commit_at(repo, since + timedelta(minutes=10), "on master")
    git(repo, "checkout", "-q", "-b", "side")
    second = commit_at(repo, since + timedelta(minutes=40, seconds=30), "on side")
    git(repo, "checkout", "-q", "master")
    return first, second


def test_back_prints_the_longest_quiet_stretch_on_any_ref(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    since = BACK_AT - timedelta(hours=3)
    flag_window(repo, since, BACK_AT + timedelta(hours=1))
    _, side = commit_on_branches(repo, since)
    lines = back_lines(repo, capsys, BACK_AT)
    assert lines[0].startswith("Back. Unattended since ")
    # side opens the longest interval, 2 h 19 min 30 s up to now; minutes floored.
    gap_start = since + timedelta(minutes=40, seconds=30)
    assert lines[1] == (
        f"{LONGEST}2 h 19 min, {hhmm(gap_start)} to {hhmm(BACK_AT)} UTC, after {side}; "
        "2 commit(s) in the window."
    )


def test_back_with_an_expired_flag_ends_the_window_at_until(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    since = BACK_AT - timedelta(hours=5)
    until = BACK_AT - timedelta(hours=1)
    flag_window(repo, since, until)
    first = commit_at(repo, since + timedelta(hours=1), "inside")
    commit_at(repo, until + timedelta(minutes=30), "after the window")
    lines = back_lines(repo, capsys, BACK_AT)
    assert lines[1] == (
        f"{LONGEST}3 h 00 min, {hhmm(since + timedelta(hours=1))} to {hhmm(until)} UTC, "
        f"after {first}; 1 commit(s) in the window."
    )


def test_back_with_no_commit_in_the_window_counts_from_its_start(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # make_repo's root commit carries today's date, after BACK_AT: outside.
    since = BACK_AT - timedelta(hours=2, minutes=5)
    flag_window(repo, since, BACK_AT + timedelta(hours=1))
    lines = back_lines(repo, capsys, BACK_AT)
    assert lines[1] == (
        f"{LONGEST}2 h 05 min, {hhmm(since)} to {hhmm(BACK_AT)} UTC, from the window's start; "
        "0 commit(s) in the window."
    )


def test_back_with_an_unreadable_flag_says_the_stretch_is_unknown(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    write_flag(repo, "{not json")
    lines = back_lines(repo, capsys, BACK_AT)
    assert lines[1] == (
        "Longest stretch without a commit: unknown (the flag has no readable since)."
    )


def test_back_with_no_flag_prints_no_stretch(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    lines = back_lines(repo, capsys, BACK_AT)
    assert lines[0] == "No unattended flag was set."
    assert not any(line.startswith("Longest stretch") for line in lines)


def test_back_when_git_log_fails_says_so_and_still_clears_the_state(
    repo: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # A ref naming a missing object: `git log --all` fails (fatal: bad object),
    # while `rev-parse --git-common-dir`, which finds the state dir, still works.
    flag_on(repo, now=NOW, keep_awake_pid=777)
    write_queue(repo, QUEUED[:1])
    (repo / ".git" / "refs" / "heads" / "broken").write_text("0" * 39 + "1\n")
    lines = back_lines(repo, capsys, BACK_AT)
    assert lines[1] == "Longest stretch without a commit: unknown (git log failed)."
    assert not flag_path(repo).exists()
    assert not (state_dir(repo) / "queue.jsonl").exists()
    assert len(archive(repo).read_text().splitlines()) == 1


def test_commits_between_reads_every_ref_and_filters_to_the_window(repo: Path) -> None:
    since = BACK_AT - timedelta(hours=3)
    first, second = commit_on_branches(repo, since)
    commit_at(repo, since - timedelta(minutes=1), "before")
    found = away.commits_between(repo, since, BACK_AT)
    assert sorted(found) == [
        (since + timedelta(minutes=10), first),
        (since + timedelta(minutes=40, seconds=30), second),
    ]


def test_commits_between_on_git_failure_is_none(repo: Path) -> None:
    (repo / ".git" / "refs" / "heads" / "broken").write_text("0" * 39 + "1\n")
    assert away.commits_between(repo, BACK_AT - timedelta(hours=3), BACK_AT) is None


#: `git log --format=%h %cI` output that exits 0 but does not parse: a stamp
#: that is not a date, a line with no stamp, and a time with no offset (naive,
#: which cannot be compared with the window's aware times).
MALFORMED_LOG = ("abc1234 not-a-date", "garbage", "abc1234 2026-09-30T20:00:00")


def fake_git_log(tmp_path: Path, output: str) -> Path:
    """A `git` that prints `output` for `log` and runs the real git otherwise."""
    real = shutil.which("git")
    assert real is not None
    bin_dir = tmp_path / "fake-git-bin"
    bin_dir.mkdir()
    script = bin_dir / "git"
    script.write_text(
        "#!/bin/sh\n"
        'for arg in "$@"; do\n'
        f"  if [ \"$arg\" = log ]; then printf '%s\\n' {shlex.quote(output)}; exit 0; fi\n"
        "done\n"
        f'exec {shlex.quote(real)} "$@"\n'
    )
    script.chmod(0o755)
    return bin_dir


@pytest.mark.parametrize("output", MALFORMED_LOG)
def test_back_with_a_malformed_git_log_line_still_clears_the_state(
    repo: Path,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
    output: str,
) -> None:
    # Code review round 1: a non-zero git exit was handled, a line that does
    # not parse ended --back before the queue was archived. §3.3: the flag and
    # the queue are not optional to --back; §3.4: the stretch is then unknown.
    flag_on(repo, now=NOW, keep_awake_pid=777)
    write_queue(repo, QUEUED[:1])
    monkeypatch.setenv("PATH", f"{fake_git_log(tmp_path, output)}{os.pathsep}{os.environ['PATH']}")
    lines = back_lines(repo, capsys, BACK_AT)
    assert lines[0].startswith("Back. Unattended since ")
    assert lines[1].startswith("Longest stretch without a commit: unknown (")
    assert "Queued while away (1):" in lines
    assert DECISIONS_HEADING in lines
    assert not flag_path(repo).exists()
    assert not (state_dir(repo) / "queue.jsonl").exists()
    assert len(archive(repo).read_text().splitlines()) == 1


# ---------------------------------------------------------------- h8, test 8
#
# §3.1: run from a worktree, `--back` lists the main checkout's decisions (the
# fault evidence §1b reproduced), and every worktree's, uncapped.

DECISIONS = DECISIONS_HEADING


def test_back_from_a_worktree_lists_the_main_checkouts_decisions(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    main = make_repo(tmp_path.resolve() / "main")
    tree = add_worktree(main, main / ".claude" / "worktrees" / "w1", "w1")
    main_tasks = main / ".claude" / "current-task"
    main_tasks.mkdir(parents=True)
    (main_tasks / "session.md").write_text(
        "NOW: h8\nQUEUE: x\n" + "".join(f"ASK OLA: main question {n}?\n" for n in range(7))
    )
    tree_tasks = tree / ".claude" / "current-task"
    tree_tasks.mkdir(parents=True)
    (tree_tasks / "tester-101010.md").write_text("ASK OLA: worktree question?\n")

    lines = back_lines(tree, capsys, BACK_AT)

    at_ = lines.index(DECISIONS)
    assert lines[at_ + 1 : at_ + 9] == [
        *(f"  session.md: ASK OLA: main question {n}?" for n in range(7)),
        "  .claude/worktrees/w1: tester-101010.md: ASK OLA: worktree question?",
    ]
