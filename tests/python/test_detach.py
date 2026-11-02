"""`tools/detach.py`: a heavy run survives the agent that launched it.

The increment is docs/increments/h40-detached-runs.md. A heavy run is launched
as a *detached job* -- its own session (setsid), writing artifacts to disk --
and re-attached later by reading those artifacts. Liveness never rides on the
launching agent's monitored stream, so a laptop suspend (SIGSTOP/SIGCONT on the
job) no longer loses the run.

These are the red tests (T1-T6 of the increment's "Testability" section),
written before `tools/detach.py` exists. They pin:

* T1  start detaches: the child lands in its own session/process group.
* T2  the child outlives its launcher (the lost-run condition, no suspend).
* T3  suspend stand-in: SIGSTOP keeps `status` reading `running`, SIGCONT
      lets it finish and flip to `done`.
* T4  the four status states -- running / done / failed / unknown -- each
      named by its own test; `unknown` is manifest host != this machine.
* T4b the done-marker contract: `done` only when the marker exists AND parses;
      a truncated or `.tmp` marker is treated as not-done.
* T5  atomic write: a marker left as `.tmp` (replace never happened) is never
      read as done.
* T6  pid-reuse guard: a foreign process that reused the pid is reported
      `failed`, not `running`, by comparing start times normalized to UTC
      within a tolerance this suite fixes (see T6_TOLERANCE_S).

Standard library only; no C++, no `@perf` acceptance (touches neither refine
nor mesh code). Fast: SIGSTOP/SIGCONT and short-lived children, no real sleeps
or suspends in the assertions.
"""

from __future__ import annotations

import contextlib
import importlib.util
import json
import os
import signal
import sys
import time
from pathlib import Path
from types import ModuleType

import pytest

# h17 §4b: a harness test; CI runs it in the `harness` job, not the product legs.
pytestmark = pytest.mark.harness

TOOL = Path(__file__).resolve().parents[2] / "tools" / "detach.py"

# T6 -- the pid-reuse clock tolerance. The increment leaves T6's tolerance to
# @tester (h40 lines 208-213): start times are compared normalized to UTC and
# must agree within this bound. It is >= 1 s because `ps -o lstart=` is 1 s
# granular (macOS) and /proc start time is clock-tick granular (Linux); a
# tighter bound would flap on those sources. We fix it at 2 s: one tick of
# slack over the 1 s granularity floor, still far below any real pid-reuse gap.
# ARCHITECT RULING NEEDED before green (choice beyond the design).
T6_TOLERANCE_S = 2.0


def _load() -> ModuleType:
    """Load tools/detach.py by path, as the repo loads every non-package tool.

    Red reason: until `tools/detach.py` is written this raises, and the message
    names the missing tool -- not an import of something unrelated.
    """
    assert TOOL.exists(), f"tools/detach.py does not exist yet: {TOOL}"
    spec = importlib.util.spec_from_file_location("detach", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------


@pytest.fixture
def detach() -> ModuleType:
    return _load()


@pytest.fixture
def run_dir(tmp_path: Path) -> Path:
    d = tmp_path / "run"
    d.mkdir()
    return d


def _wait_for(pred, timeout: float = 5.0, step: float = 0.01) -> bool:
    """Poll pred() until true or the timeout; a bounded wait on a child, not a sleep.

    Used only to let a short-lived child reach a state (write its marker, exit);
    the assertions themselves do not sleep. Bounded so a hung child fails fast.
    """
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if pred():
            return True
        time.sleep(step)
    return pred()


def _manifest(run_dir: Path, name: str) -> dict:
    return json.loads((run_dir / f"{name}.detached.json").read_text())


def _child_that_writes(marker: Path, extra: str = "") -> list[str]:
    """A /bin/sh command that writes `marker` atomically after a short beat.

    The beat is a 0.05 s shell sleep inside the *child* (not the test), so the
    launcher can return and (T2) even exit before the marker lands.
    """
    tmp = f"{marker}.tmp"
    return [
        "/bin/sh",
        "-c",
        f"sleep 0.05; {extra} printf '{{}}' > '{tmp}'; mv '{tmp}' '{marker}'",
    ]


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _reap(pid: int) -> None:
    """Reap a child we signalled, so the test leaves no zombie or live process."""
    try:
        os.kill(pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        return
    with contextlib.suppress(ChildProcessError, OSError):
        os.waitpid(pid, 0)


# ---------------------------------------------------------------------------
# T1 -- start detaches the child into its own session / process group
# ---------------------------------------------------------------------------


def test_t1_start_puts_child_in_its_own_session(detach: ModuleType, run_dir: Path) -> None:
    marker = run_dir / "t1.done.json"
    pid = detach.start(
        run_dir=run_dir,
        name="t1",
        cmd=_child_that_writes(marker),
        done_marker="t1.done.json",
    )
    try:
        # The whole detachment: the child is session leader of a new session,
        # so its pgid is itself and differs from this test's process group.
        assert os.getpgid(pid) == pid
        assert os.getpgid(pid) != os.getpgid(0)
        m = _manifest(run_dir, "t1")
        assert m["pid"] == pid
        assert m["pgid"] == pid
        assert m["state"] == "running"
        assert m["done_marker"] == "t1.done.json"
    finally:
        _wait_for(lambda: marker.exists())
        _reap(pid)


def test_t1_start_returns_immediately_without_waiting(detach: ModuleType, run_dir: Path) -> None:
    # The launcher must not block on the job (it hands back at once). The child
    # sleeps 0.05 s before writing; start() must return well before that lands.
    marker = run_dir / "t1b.done.json"
    t0 = time.monotonic()
    pid = detach.start(
        run_dir=run_dir, name="t1b", cmd=_child_that_writes(marker), done_marker="t1b.done.json"
    )
    elapsed = time.monotonic() - t0
    try:
        assert elapsed < 0.05, f"start() waited {elapsed:.3f}s on the job"
        assert not marker.exists(), "marker already written: start() waited on the child"
    finally:
        _wait_for(lambda: marker.exists())
        _reap(pid)


def test_t1_start_creates_missing_run_dir(detach: ModuleType, tmp_path: Path) -> None:
    missing = tmp_path / "nope" / "deep"
    marker = missing / "t1c.done.json"
    pid = detach.start(
        run_dir=missing, name="t1c", cmd=_child_that_writes(marker), done_marker="t1c.done.json"
    )
    try:
        assert missing.is_dir()
        assert (missing / "t1c.detached.json").exists()
    finally:
        _wait_for(lambda: marker.exists())
        _reap(pid)


def test_t1_stdin_is_devnull_so_job_never_blocks_on_input(
    detach: ModuleType, run_dir: Path
) -> None:
    # The 2026-10-30 reap case: a job blocked on input. A detached job reads
    # EOF at once, so a child that reads stdin still finishes.
    marker = run_dir / "t1d.done.json"
    tmp = f"{marker}.tmp"
    cmd = ["/bin/sh", "-c", f"cat; printf '{{}}' > '{tmp}'; mv '{tmp}' '{marker}'"]
    pid = detach.start(run_dir=run_dir, name="t1d", cmd=cmd, done_marker="t1d.done.json")
    try:
        assert _wait_for(lambda: marker.exists()), "job blocked on stdin (not /dev/null)"
    finally:
        _reap(pid)


# ---------------------------------------------------------------------------
# T2 -- the child outlives its launcher (the lost-run condition, no suspend)
# ---------------------------------------------------------------------------


def test_t2_child_survives_launcher_process_exit(run_dir: Path) -> None:
    # Reproduce the lost-run condition: the *launcher* dies, the job lives on
    # and writes its marker. The launcher stand-in is a separate Python process
    # that calls detach.start() and exits 0 at once, not this test process.
    marker = run_dir / "t2.done.json"
    launcher = (
        "import importlib.util,sys;"
        f"spec=importlib.util.spec_from_file_location('detach', {str(TOOL)!r});"
        "assert spec and spec.loader, 'tools/detach.py missing';"
        "m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m);"
        f"pid=m.start(run_dir={str(run_dir)!r}, name='t2',"
        f" cmd=['/bin/sh','-c',\"sleep 0.2; printf '{{}}' > '{marker}.tmp';"
        f" mv '{marker}.tmp' '{marker}'\"], done_marker='t2.done.json');"
        "print(pid)"
    )
    import subprocess

    proc = subprocess.run(
        [sys.executable, "-c", launcher], capture_output=True, text=True, timeout=10
    )
    assert proc.returncode == 0, f"launcher stand-in failed: {proc.stderr}"
    child_pid = int(proc.stdout.strip())
    try:
        # Launcher has already exited (subprocess.run returned). The marker is
        # not there yet (child sleeps 0.2 s); it must appear anyway.
        assert not marker.exists(), "marker present before launcher exit: not detached"
        assert _wait_for(lambda: marker.exists()), "child died with its launcher"
    finally:
        _reap(child_pid)


# ---------------------------------------------------------------------------
# T3 -- suspend stand-in: SIGSTOP keeps `running`, SIGCONT lets it finish
# ---------------------------------------------------------------------------


def test_t3_sigstop_reads_running_then_sigcont_reads_done(
    detach: ModuleType, run_dir: Path
) -> None:
    marker = run_dir / "t3.done.json"
    pid = detach.start(
        run_dir=run_dir,
        name="t3",
        cmd=_child_that_writes(marker, extra="sleep 0.15;"),
        done_marker="t3.done.json",
    )
    try:
        # A laptop suspend is exactly SIGSTOP on the process. status() must
        # still read `running` -- pid alive, no done-marker.
        os.kill(pid, signal.SIGSTOP)
        states = detach.status(run_dir=run_dir)
        assert states["t3"] == "running", f"stopped job not running: {states}"
        assert not marker.exists()
        # Resume (SIGCONT, the wake); the job completes and status flips to done.
        os.kill(pid, signal.SIGCONT)
        assert _wait_for(lambda: detach.status(run_dir=run_dir).get("t3") == "done"), (
            "job did not flip to done after SIGCONT"
        )
        assert marker.exists()
    finally:
        _reap(pid)


# ---------------------------------------------------------------------------
# T4 -- the four status states, each named
# ---------------------------------------------------------------------------


def test_t4_running_when_pid_alive_and_no_done_marker(detach: ModuleType, run_dir: Path) -> None:
    marker = run_dir / "run.done.json"
    pid = detach.start(
        run_dir=run_dir,
        name="running",
        cmd=_child_that_writes(marker, extra="sleep 0.3;"),
        done_marker="run.done.json",
    )
    try:
        assert detach.status(run_dir=run_dir)["running"] == "running"
    finally:
        _reap(pid)


def test_t4_done_when_marker_present_and_parses(detach: ModuleType, run_dir: Path) -> None:
    marker = run_dir / "d.done.json"
    pid = detach.start(
        run_dir=run_dir, name="d", cmd=_child_that_writes(marker), done_marker="d.done.json"
    )
    try:
        assert _wait_for(lambda: marker.exists())
        # A present+parsing done-marker reads `done` regardless of the pid:
        # the done-marker check wins over liveness (design §"Partial vs
        # complete"), so the job need not have exited for `done` to hold.
        assert detach.status(run_dir=run_dir)["d"] == "done"
    finally:
        _reap(pid)


def test_t4_failed_when_no_done_marker_and_pid_gone(detach: ModuleType, run_dir: Path) -> None:
    # Kill a job mid-run: no done-marker, pid gone -> failed (crashed/killed).
    marker = run_dir / "f.done.json"
    pid = detach.start(
        run_dir=run_dir,
        name="f",
        cmd=_child_that_writes(marker, extra="sleep 5;"),
        done_marker="f.done.json",
    )
    os.kill(pid, signal.SIGKILL)
    with contextlib.suppress(ChildProcessError, OSError):
        os.waitpid(pid, 0)
    assert _wait_for(lambda: not _alive(pid))
    assert not marker.exists()
    assert detach.status(run_dir=run_dir)["f"] == "failed"


def test_t4_unknown_when_manifest_host_is_not_this_machine(
    detach: ModuleType, run_dir: Path
) -> None:
    # A manifest written on another host: liveness cannot be judged here, so
    # the pid is NOT probed and the state is `unknown`. We write the manifest
    # by hand with a foreign host and a pid that is alive on THIS machine (our
    # own), to prove the host guard wins over any kill(pid,0) result.
    name = "u"
    manifest = {
        "name": name,
        "pid": os.getpid(),
        "pgid": os.getpgid(0),
        "cmd": ["python3", "x.py"],
        "cwd": str(run_dir),
        "started": "2026-11-02T10:46:00Z",
        "host": "some-other-host.example",
        "stdout": f"{name}.stdout",
        "stderr": f"{name}.stderr",
        "done_marker": f"{name}.done.json",
        "state": "running",
    }
    (run_dir / f"{name}.detached.json").write_text(json.dumps(manifest))
    assert detach.status(run_dir=run_dir)[name] == "unknown"


# ---------------------------------------------------------------------------
# T4b -- the done-marker contract: present AND parses
# ---------------------------------------------------------------------------


def test_t4b_truncated_done_marker_is_not_done(detach: ModuleType, run_dir: Path) -> None:
    # The marker file exists but does not parse (a writer exiting mid-write).
    # The design: a truncated marker is treated as not-done; pid-liveness then
    # decides. Here the pid is gone, so the state is `failed`, never `done`.
    name = "trunc"
    manifest = {
        "name": name,
        "pid": 2147480000,  # a pid that is not running
        "pgid": 2147480000,
        "cmd": ["python3", "x.py"],
        "cwd": str(run_dir),
        "started": "2026-11-02T10:46:00Z",
        "host": _load().this_host(),
        "stdout": f"{name}.stdout",
        "stderr": f"{name}.stderr",
        "done_marker": f"{name}.done.json",
        "state": "running",
    }
    (run_dir / f"{name}.detached.json").write_text(json.dumps(manifest))
    (run_dir / f"{name}.done.json").write_text('{"partial": ')  # truncated JSON
    state = detach.status(run_dir=run_dir)[name]
    assert state != "done", "truncated marker read as done"
    assert state == "failed"


# ---------------------------------------------------------------------------
# T5 -- atomic write: a `.tmp` marker is never read as done
# ---------------------------------------------------------------------------


def test_t5_tmp_marker_is_never_read_as_done(detach: ModuleType, run_dir: Path) -> None:
    # Interrupt before os.replace: the `.tmp` file is left and the final name
    # never appears. status() keys only on the final done_marker name, so the
    # `.tmp` is invisible to the done-check.
    name = "atom"
    manifest = {
        "name": name,
        "pid": 2147480000,  # not running
        "pgid": 2147480000,
        "cmd": ["python3", "x.py"],
        "cwd": str(run_dir),
        "started": "2026-11-02T10:46:00Z",
        "host": _load().this_host(),
        "stdout": f"{name}.stdout",
        "stderr": f"{name}.stderr",
        "done_marker": f"{name}.done.json",
        "state": "running",
    }
    (run_dir / f"{name}.detached.json").write_text(json.dumps(manifest))
    # A fully valid JSON, but written to the .tmp name -- replace never ran.
    (run_dir / f"{name}.done.json.tmp").write_text('{"done": true}')
    assert not (run_dir / f"{name}.done.json").exists()
    state = detach.status(run_dir=run_dir)[name]
    assert state != "done", ".tmp marker read as done -- half-written result exposed"
    assert state == "failed"


def test_t5_atomic_write_helper_leaves_no_partial_final_name(
    detach: ModuleType, run_dir: Path
) -> None:
    # The write path detach.py uses for its own manifest: write to <name>.tmp
    # in the SAME directory, then os.replace. Pin that the final name and the
    # tmp name share a directory (cross-filesystem replace is not atomic).
    target = run_dir / "payload.json"
    detach.atomic_write_json(target, {"k": 1})
    assert target.exists()
    assert json.loads(target.read_text()) == {"k": 1}
    assert not (run_dir / "payload.json.tmp").exists()


# ---------------------------------------------------------------------------
# T6 -- pid-reuse guard: a foreign process that reused the pid reads `failed`
# ---------------------------------------------------------------------------


def test_t6_reused_pid_with_mismatched_start_reads_failed(
    detach: ModuleType, run_dir: Path
) -> None:
    # Point the manifest at THIS test's own pid -- a process that is alive but
    # is NOT the original job -- with a `started` far from this process's real
    # start. A naive kill(pid,0) would say `running`; the start-time guard must
    # say `failed`. This is the pid-reuse false positive the guard forbids.
    name = "reuse"
    manifest = {
        "name": name,
        "pid": os.getpid(),  # alive, but a foreign process (us)
        "pgid": os.getpgid(0),
        "cmd": ["python3", "x.py"],
        "cwd": str(run_dir),
        # Far in the past, well beyond T6_TOLERANCE_S from our real start time.
        "started": "2001-01-01T00:00:00Z",
        "host": _load().this_host(),
        "stdout": f"{name}.stdout",
        "stderr": f"{name}.stderr",
        "done_marker": f"{name}.done.json",
        "state": "running",
    }
    (run_dir / f"{name}.detached.json").write_text(json.dumps(manifest))
    assert detach.status(run_dir=run_dir)[name] == "failed"


def test_t6_matching_start_within_tolerance_reads_running(
    detach: ModuleType, run_dir: Path
) -> None:
    # The other side of the guard: when the running pid's start time matches the
    # manifest's `started` within T6_TOLERANCE_S, it IS the original job and
    # reads `running`. Use this process: read its own start via detach's own
    # start-time reader so the test and the tool agree on the clock.
    name = "match"
    pid = os.getpid()
    started = detach.pid_start_utc(pid)  # tool's own reader -> an ISO-8601 UTC string
    manifest = {
        "name": name,
        "pid": pid,
        "pgid": os.getpgid(0),
        "cmd": ["python3", "x.py"],
        "cwd": str(run_dir),
        "started": started,
        "host": _load().this_host(),
        "stdout": f"{name}.stdout",
        "stderr": f"{name}.stderr",
        "done_marker": f"{name}.done.json",
        "state": "running",
    }
    (run_dir / f"{name}.detached.json").write_text(json.dumps(manifest))
    assert detach.status(run_dir=run_dir)[name] == "running"


def test_t6_tolerance_is_at_least_one_second(detach: ModuleType) -> None:
    # The guard tolerance must be >= 1 s because `ps -o lstart=` is 1 s granular
    # and /proc start time is clock-tick granular (h40 lines 208-213). Pin the
    # tool's own tolerance so a future tightening below the source granularity
    # (which would flap) turns this red.
    assert detach.START_TOLERANCE_S >= 1.0
    assert detach.START_TOLERANCE_S == T6_TOLERANCE_S


# ---------------------------------------------------------------------------
# Degeneracy -- the increment's "Size, exclusions, degeneracy" bullet
# ---------------------------------------------------------------------------


def test_missing_run_dir_status_is_empty_not_a_crash(detach: ModuleType, tmp_path: Path) -> None:
    assert detach.status(run_dir=tmp_path / "absent") == {}


def test_second_start_same_name_refuses_without_force(detach: ModuleType, run_dir: Path) -> None:
    marker = run_dir / "dup.done.json"
    pid = detach.start(
        run_dir=run_dir, name="dup", cmd=_child_that_writes(marker), done_marker="dup.done.json"
    )
    try:
        with pytest.raises(FileExistsError):
            detach.start(
                run_dir=run_dir,
                name="dup",
                cmd=_child_that_writes(marker),
                done_marker="dup.done.json",
            )
    finally:
        _wait_for(lambda: marker.exists())
        _reap(pid)
