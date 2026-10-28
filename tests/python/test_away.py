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
import hashlib
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
    Granted,
    PrivateGranted,
    Tool,
    add_worktree,
    clean_env,
    commit_at,
    commit_file,
    flag_on,
    flag_path,
    flag_window,
    git,
    grants_used_path,
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
    # h20 §3.1: the enter entry carries the grants, `[]` without --grant.
    assert windows(repo) == [
        {"event": "enter", "since": flag["since"], "until": flag["until"], "grants": []}
    ]
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
DECISIONS_HEADING = "Waiting on Ola:"


def ask(repo: Path, *texts: str) -> None:
    """Open questions in the store `tools/questions.py` reads, written without its checks."""
    with (state_dir(repo).parent / "questions.jsonl").open("a") as out:
        for n, text in enumerate(texts, 1):
            out.write(json.dumps({"id": f"q{n}", "kind": "question", "text": text,
                                  "default": "no", "choices": [], "item": None,
                                  "source": "main", "state": "open", "words": None,
                                  "at": "2026-09-30T19:00:00+00:00"}) + "\n")  # fmt: skip


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
    ask(repo, "May I push?")

    code, stop = back(repo)

    assert code == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].startswith(f"Back. Unattended since {iso(NOW)} until {until}")
    at = lines.index("Queued while away (3):")
    assert lines[at + 1 : at + 9] == [
        "  feat-a",
        "    2026-09-30T19:00:00+00:00 guard_push main: git push",
        "    2026-09-30T21:00:00+00:00 guard_push developer: gh pr merge 1",
        "  feat-b",
        "    2026-09-30T20:00:00+00:00 guard_governance tester: Edit /r/CLAUDE.md",
        "Push grants: (none)",
        DECISIONS_HEADING,
        "  1. May I push? Default: no.",
    ]
    # It is the question tool's `list`: Ola's "1" now means the question he saw.
    shown = json.loads((state_dir(repo).parent / "questions-shown.json").read_text())
    assert shown == {"1": "q1"}
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


def test_back_survives_a_broken_question_tool(repo: Path) -> None:
    # §3.3 (review 1): the recap module is optional to --back; the flag, the
    # keep-awake and the queue are not.
    (repo / "tools" / "questions.py").write_text('raise RuntimeError("planted")\n')
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
    at = lines.index(DECISIONS_HEADING)
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
# Run from a worktree, `--back` lists the questions of the store every
# checkout shares, uncapped.

DECISIONS = DECISIONS_HEADING


def test_back_from_a_worktree_lists_the_shared_questions_uncapped(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    main = make_repo(tmp_path.resolve() / "main")
    tree = add_worktree(main, main / ".claude" / "worktrees" / "w1", "w1")
    state_dir(main).mkdir(parents=True, exist_ok=True)
    ask(main, *(f"Main question {n}?" for n in range(8)))

    lines = back_lines(tree, capsys, BACK_AT)

    at_ = lines.index(DECISIONS)
    assert lines[at_ + 1 : at_ + 9] == [
        f"  {n + 1}. Main question {n}? Default: no." for n in range(8)
    ]


# ---------------------------------------------------------------- h20 §4.1
#
# docs/increments/h20-away-push-grants.md §3.1 and §3.6: `--grant <branch>`
# names the branches an away window may push, open and enqueue; `--back`
# reports each grant. Pinned beyond the design's words: a `grants-used.jsonl`
# line's `act` is one of the act names of §3.2's table ("push", "pull
# request", "enqueue"); `--back` names a head by its first 7 hex digits.

GRANTED_TEXT = "Granted while away: push, pull request and one enqueue per reviewed head of: "


def test_1_each_grant_is_written_to_the_flag_and_the_window_and_shown(repo: Path) -> None:
    git(repo, "branch", "feat")
    git(repo, "branch", "other")
    code, tty, _, _ = enter(repo, ["8h", "--grant", "feat", "--grant", "other"])
    assert code == 0
    assert read_flag(repo)["grants"] == ["feat", "other"]
    [window] = windows(repo)
    assert window["grants"] == ["feat", "other"]
    assert GRANTED_TEXT + "feat, other." in tty.written


def test_2_no_grant_writes_an_empty_list(repo: Path) -> None:
    code, _, _, _ = enter(repo, ["8h"])
    assert code == 0
    assert read_flag(repo)["grants"] == []
    assert windows(repo)[0]["grants"] == []


@pytest.mark.parametrize("name", ["master", "main", "HEAD", "a..b", "~x", "=x", "a*"])
def test_3_a_grant_that_is_not_a_plain_branch_name_exits_2_before_the_terminal(
    repo: Path, name: str
) -> None:
    opened = Recorder(raises=AssertionError("the terminal was opened"))
    spawn = spawner()
    code = run_main(
        ["8h", "--grant", name], root=repo, now=NOW, open_tty=opened, spawn=spawn, stop=Recorder()
    )
    assert code == 2
    assert opened.calls == []
    assert spawn.calls == []
    assert not flag_path(repo).exists()


@pytest.mark.parametrize("argv", [["--back", "--grant", "feat"], ["--grant", "feat"]])
def test_4_a_grant_without_a_duration_exits_2(repo: Path, argv: list[str]) -> None:
    flag_on(repo, now=NOW)
    code = run_main(argv, root=repo, now=NOW, open_tty=no_tty, spawn=spawner(), stop=Recorder())
    assert code == 2
    assert flag_path(repo).exists(), "--back ran despite the --grant"


def test_5_a_grant_for_a_branch_that_does_not_exist_yet_is_accepted(repo: Path) -> None:
    code, tty, _, _ = enter(repo, ["8h", "--grant", "later"])
    assert code == 0
    assert read_flag(repo)["grants"] == ["later"]
    assert "does not exist yet" in tty.written


def test_6_re_entering_replaces_the_grants(repo: Path) -> None:
    enter(repo, ["8h", "--grant", "a"])
    code, _, _, _ = enter(repo, ["8h", "--grant", "b"])
    assert code == 0
    assert read_flag(repo)["grants"] == ["b"]


def used(at: str, window: str, act: str, head: str, review: tuple[str, str]) -> dict[str, Any]:
    return {"at": at, "window": window, "branch": "feat", "act": act, "head": head,
            "review": list(review)}  # fmt: skip


def test_7_back_reports_this_windows_granted_acts_and_each_unused_grant(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    built = Granted(tmp_path / "repo", flag=False)
    repo, head = built.repo, built.head
    git(repo, "branch", "good", head)
    git(repo, "checkout", "-q", "-b", "bad", head)
    late = commit_file(repo, "tools/b.py", "b = 1\n", "a code commit after the review")
    git(repo, "checkout", "-q", "feat")
    write_flag(repo, {"since": iso(NOW), "until": iso(NOW + timedelta(hours=8)),
                      "set_by": "away", "keep_awake_pid": None,
                      "grants": ["feat", "good", "bad"]})  # fmt: skip
    review = (built.x, "2026-09-30T18:30:00+00:00")
    earlier = "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    log = [
        used("2026-09-29T21:00:00+00:00", iso(NOW - timedelta(days=1)), "push", earlier, review),
        used("2026-09-30T19:10:00+00:00", iso(NOW), "push", head, review),
        used("2026-09-30T19:20:00+00:00", iso(NOW), "enqueue", head, review),
    ]
    grants_used_path(repo).write_text("".join(json.dumps(line) + "\n" for line in log))

    lines = back_lines(repo, capsys, BACK_AT)

    at = lines.index("Push grants (3):")
    assert at > next(i for i, line in enumerate(lines) if line.startswith("Queued while away"))
    assert at < next(i for i, line in enumerate(lines) if line == DECISIONS_HEADING)
    feat, good, bad = lines[at + 1 : at + 4]
    assert feat.startswith("  feat: ")
    for part in (f"pushed {head[:7]} at 19:10 UTC", f"enqueued {head[:7]} at 19:20 UTC",
                 f"review of {built.x[:7]}", "recorded 18:30 UTC"):  # fmt: skip
        assert part in feat, feat
    assert earlier[:7] not in feat, "an earlier window's act was reported"
    assert good == f"  good: not used; head {head[:7]} reviewed: yes"
    assert bad.startswith(f"  bad: not used; head {late[:7]} reviewed: no ("), bad
    assert "tools/b.py" in bad, bad
    assert "not judged" not in bad, "a review reason was reported as an unjudged head"


def test_7_back_with_no_grants_says_none(repo: Path, capsys: pytest.CaptureFixture[str]) -> None:
    flag_on(repo, now=NOW)
    assert "Push grants: (none)" in back_lines(repo, capsys, BACK_AT)


# ---------------------------------------------------------------- h22 §7
#
# docs/increments/h22-harness-fixes.md §4.4: `--back` judges a private
# branch's review against private/master, and names a private merge.


def private_back_lines(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> tuple[PrivateGranted, str, str, str]:
    """Fixture A with prv merged with a newer private/master after its review, a
    flag granting prv and `done`, and a logged private merge of `done`; returns
    (fixture, prv's head, prv's report line, done's report line)."""
    built = PrivateGranted(tmp_path / "repo")
    repo = built.repo
    git(repo, "checkout", "-q", "-b", "scratch-private", "refs/remotes/private/master")
    p2 = commit_file(repo, "src/p2.py", "p = 2\n", "P2 on private/master")
    git(repo, "update-ref", "refs/remotes/private/master", p2)
    git(repo, "checkout", "-q", "prv")
    git(repo, "branch", "-q", "-D", "scratch-private")
    git(repo, "merge", "-q", "--no-ff", "--no-edit", "refs/remotes/private/master")
    write_flag(repo, {"since": iso(NOW), "until": iso(NOW + timedelta(hours=8)),
                      "set_by": "away", "keep_awake_pid": None,
                      "grants": ["prv", "done"]})  # fmt: skip
    merged = {"at": "2026-09-30T19:40:00+00:00", "window": iso(NOW), "branch": "done",
              "act": "merge", "head": built.prv, "remote": "private",
              "review": [built.prv, "2026-09-30T18:30:00+00:00"]}  # fmt: skip
    grants_used_path(repo).write_text(json.dumps(merged) + "\n")
    lines = back_lines(repo, capsys, BACK_AT)
    at = lines.index("Push grants (2):")
    return built, built.head("prv"), lines[at + 1], lines[at + 2]


def test_h22_back_reviews_a_private_branch_against_private_master(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A clean merge of private/master after the review passes against private/master;
    against origin/master it would not be a master merge."""
    _, head, prv, _ = private_back_lines(tmp_path, capsys)
    assert prv == f"  prv: not used; head {head[:7]} reviewed: yes", prv


def test_h22_back_names_a_private_merge(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    built, _, _, done = private_back_lines(tmp_path, capsys)
    assert done.startswith("  done: "), done
    assert f"merged {built.prv[:7]} at 19:40 UTC" in done, done


def test_h22_back_says_not_judged_when_the_branchs_repository_cannot_be_told(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Review round 1 follow-up: with refs/remotes/private/master deleted, home()
    gives a reason; `--back` reports the head as not judged, never as unreviewed."""
    built = PrivateGranted(tmp_path / "repo")
    built.delete_private_ref()
    write_flag(built.repo, {"since": iso(NOW), "until": iso(NOW + timedelta(hours=8)),
                            "set_by": "away", "keep_awake_pid": None,
                            "grants": ["prv", "pub"]})  # fmt: skip
    missing = ("the private remote is set up but refs/remotes/private/master is missing, so "
               "it cannot be told whether this branch holds private commits; run git fetch "
               "private")  # fmt: skip
    lines = back_lines(built.repo, capsys, BACK_AT)
    at = lines.index("Push grants (2):")
    prv, pub = lines[at + 1 : at + 3]
    assert prv == f"  prv: not used; head {built.prv[:7]} not judged ({missing})", prv
    assert pub == f"  pub: not used; head {built.pub[:7]} not judged ({missing})", pub


# ---------------------------------------------------------------- h28 §2, §6 tests 5 to 11
#
# `away.py` refuses to enter while a recorded permission prompt has no
# tool_result in its session's transcript. Every transcript here is a fake of
# two to four lines under tmp_path, in the shape Claude Code writes: a
# top-level `type` of `assistant` or `user`, and the blocks under
# `message.content`.

PROMPT_INPUT = {"command": "git push", "description": "Push the branch"}
PROMPT_SESSION = "5e55102a-0000-4000-8000-000000000001"
PROMPT_AT = NOW - timedelta(minutes=5)
NOT_JSON = '{"type": "assistant", "tool_use" cut off'


def input_sha(supplied: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(supplied, sort_keys=True).encode()).hexdigest()


def tool_use(block_id: str, supplied: dict[str, Any] | None = None) -> dict[str, Any]:
    block = {"type": "tool_use", "id": block_id, "name": "Bash",
             "input": PROMPT_INPUT if supplied is None else supplied}  # fmt: skip
    return {"type": "assistant", "message": {"role": "assistant", "content": [block]}}


def tool_result(block_id: str, is_error: bool = False) -> dict[str, Any]:
    text = "The user doesn't want to proceed with this tool use" if is_error else "pushed"
    block = {"type": "tool_result", "tool_use_id": block_id, "content": text,
             "is_error": is_error}  # fmt: skip
    return {"type": "user", "message": {"role": "user", "content": [block]}}


def unrelated() -> dict[str, Any]:
    return {"type": "user", "message": {"role": "user", "content": "Push it, please."}}


def write_transcript(path: Path, entries: list[dict[str, Any] | str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [e if isinstance(e, str) else json.dumps(e) for e in entries]
    path.write_text("".join(line + "\n" for line in lines))
    return path


def main_transcript(tmp_path: Path) -> Path:
    return tmp_path / "projects" / f"{PROMPT_SESSION}.jsonl"


def record_prompt(
    repo: Path, transcript: Path, at: datetime = PROMPT_AT,
    agent_id: str | None = None, agent_type: str | None = None,
) -> None:  # fmt: skip
    """One line of `prompts.jsonl`, as the permission hook writes it (h28 §2)."""
    entry = {"at": iso(at), "session_id": PROMPT_SESSION, "agent_id": agent_id,
             "agent_type": agent_type, "transcript_path": str(transcript), "tool_name": "Bash",
             "input_sha": input_sha(PROMPT_INPUT), "act": "Bash: git push"}  # fmt: skip
    state_dir(repo).mkdir(parents=True, exist_ok=True)
    with (state_dir(repo) / "prompts.jsonl").open("a") as out:
        out.write(json.dumps(entry) + "\n")


class TtyOpener:
    """`open_tty` that counts its calls and answers `y`."""

    def __init__(self) -> None:
        self.calls = 0

    def __call__(self) -> FakeTty:
        self.calls += 1
        return FakeTty("y")


@dataclass
class Entered:
    code: int
    opener: TtyOpener
    spawn: Recorder
    err: str
    out: str


def enter_with_prompts(repo: Path, capsys: pytest.CaptureFixture[str], *extra: str) -> Entered:
    opener, spawn = TtyOpener(), spawner()
    capsys.readouterr()
    code = run_main(["8h", *extra], root=repo, now=NOW, open_tty=opener, spawn=spawn,
                    stop=Recorder())  # fmt: skip
    captured = capsys.readouterr()
    return Entered(code, opener, spawn, captured.err, captured.out)


def assert_entered(repo: Path, entered: Entered) -> None:
    assert entered.code == 0, entered.err
    assert flag_path(repo).exists()
    assert entered.opener.calls == 1


def assert_refused(repo: Path, entered: Entered) -> None:
    assert entered.code == 1, f"exit {entered.code}; stderr: {entered.err}"
    assert "Bash: git push" in entered.err, entered.err
    assert "--despite-prompts" in entered.err, entered.err
    assert PROMPT_SESSION[:8] in entered.err, entered.err
    assert not flag_path(repo).exists()
    assert not (state_dir(repo) / "windows.jsonl").exists()
    assert entered.opener.calls == 0, "the terminal was opened before the prompt check"
    assert entered.spawn.calls == [], "a keep-awake was started"


def test_h28_an_open_prompt_refuses_the_window_before_the_terminal(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    transcript = write_transcript(main_transcript(tmp_path), [unrelated(), tool_use("toolu_1")])
    record_prompt(repo, transcript)
    assert_refused(repo, enter_with_prompts(repo, capsys))


@pytest.mark.parametrize("is_error", [False, True], ids=["answered", "rejected"])
def test_h28_a_prompt_whose_call_has_a_result_does_not_refuse(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str], is_error: bool
) -> None:
    """A rejected prompt writes a tool_result with is_error true; it is not open.
    A line that does not parse is skipped."""
    transcript = write_transcript(
        main_transcript(tmp_path),
        [NOT_JSON, tool_use("toolu_1"), tool_result("toolu_1", is_error), unrelated()],
    )
    record_prompt(repo, transcript)
    assert_entered(repo, enter_with_prompts(repo, capsys))


def test_h28_a_missing_transcript_means_the_session_is_gone(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    record_prompt(repo, main_transcript(tmp_path))
    assert_entered(repo, enter_with_prompts(repo, capsys))


def test_h28_a_transcript_without_the_call_is_taken_as_lagging(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Only a tool_use with the same input counts; another command's answered call does not."""
    other = {"command": "git status"}
    transcript = write_transcript(
        main_transcript(tmp_path),
        [unrelated(), tool_use("toolu_9", other), tool_result("toolu_9")],
    )
    record_prompt(repo, transcript)
    assert_refused(repo, enter_with_prompts(repo, capsys))


@pytest.mark.parametrize(
    ("answered", "refuses"), [("toolu_1", True), ("toolu_2", False)],
    ids=["first-answered-last-open", "first-open-last-answered"],
)  # fmt: skip
def test_h28_the_last_matching_call_decides(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str], answered: str, refuses: bool
) -> None:
    transcript = write_transcript(
        main_transcript(tmp_path), [tool_use("toolu_1"), tool_use("toolu_2"), tool_result(answered)]
    )
    record_prompt(repo, transcript)
    entered = enter_with_prompts(repo, capsys)
    if refuses:
        assert_refused(repo, entered)
    else:
        assert_entered(repo, entered)


LAST_ENTER = NOW - timedelta(hours=1)


@pytest.mark.parametrize(
    ("at", "refuses"),
    [(LAST_ENTER - timedelta(minutes=1), False), (LAST_ENTER, True),
     (LAST_ENTER + timedelta(minutes=1), True)],
    ids=["before-the-last-enter", "at-the-last-enter", "after-the-last-enter"],
)  # fmt: skip
def test_h28_only_prompts_since_the_last_window_start_count(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str], at: datetime, refuses: bool
) -> None:
    windows_path = state_dir(repo) / "windows.jsonl"
    windows_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        {"event": "enter", "since": iso(NOW - timedelta(hours=3)),
         "until": iso(NOW - timedelta(hours=2)), "grants": []},
        {"event": "enter", "since": iso(LAST_ENTER),
         "until": iso(LAST_ENTER + timedelta(minutes=30)), "grants": []},
        {"event": "back", "since": iso(LAST_ENTER),
         "until": iso(LAST_ENTER + timedelta(minutes=20))},
    ]  # fmt: skip
    windows_path.write_text("".join(json.dumps(line) + "\n" for line in lines))
    transcript = write_transcript(main_transcript(tmp_path), [tool_use("toolu_1")])
    record_prompt(repo, transcript, at=at)
    entered = enter_with_prompts(repo, capsys)
    if refuses:
        assert entered.code == 1, entered.err
        assert "Bash: git push" in entered.err, entered.err
    else:
        assert_entered(repo, entered)


SUBAGENT_ID = "a1b2c3d4e5"


def subagent_transcript(tmp_path: Path) -> Path:
    main = main_transcript(tmp_path)
    return main.with_suffix("") / "subagents" / f"agent-{SUBAGENT_ID}.jsonl"


SUBAGENT_CASES = {
    # (main transcript, subagent transcript or None for absent, refuses)
    "result-only-in-subagent-file": ([unrelated()],
                                     [tool_use("toolu_1"), tool_result("toolu_1")], False),
    "result-in-neither": ([unrelated()], [tool_use("toolu_1")], True),
    "subagent-file-wins-over-main": ([tool_use("toolu_1"), tool_result("toolu_1")],
                                     [tool_use("toolu_1")], True),
    "no-subagent-file-falls-back-to-main": ([tool_use("toolu_1"), tool_result("toolu_1")],
                                            None, False),
}  # fmt: skip


@pytest.mark.parametrize("case", SUBAGENT_CASES.values(), ids=SUBAGENT_CASES.keys())
def test_h28_a_subagents_prompt_is_judged_in_its_own_transcript(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str], case: Any
) -> None:
    main_lines, sub_lines, refuses = case
    transcript = write_transcript(main_transcript(tmp_path), main_lines)
    if sub_lines is not None:
        write_transcript(subagent_transcript(tmp_path), sub_lines)
    record_prompt(repo, transcript, agent_id=SUBAGENT_ID, agent_type="tester")
    entered = enter_with_prompts(repo, capsys)
    if refuses:
        assert_refused(repo, entered)
        assert "tester" in entered.err, entered.err
    else:
        assert_entered(repo, entered)


def test_h28_despite_prompts_enters_and_warns(
    repo: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    transcript = write_transcript(main_transcript(tmp_path), [tool_use("toolu_1")])
    record_prompt(repo, transcript)
    entered = enter_with_prompts(repo, capsys, "--despite-prompts")
    assert_entered(repo, entered)
    assert "Bash: git push" in entered.err + entered.out, (entered.err, entered.out)
    assert [w["event"] for w in windows(repo)] == ["enter"]
