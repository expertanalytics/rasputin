"""The suite's own reads: no test may depend on a prose file (h11, T5).

`docs/increments/h11-ci-path-filter.md` §2.1 lets a pull request that changes
only prose skip the C++ and Python jobs. That is safe only while no suite reads
a file `tools/ci_changes.py` calls prose. This hook enforces it: an audit hook
records every file the test process opens, and at session end the session
fails if any of them is a prose file of this repository, naming the file and
saying to add it to `NOT_PROSE` in `tools/ci_changes.py` (or to stop reading
it).

The repository root is found from this file's own location, so the hook also
works when loaded as a plugin from elsewhere (`-p conftest` with
`tests/python` on the path, as T5's planted run does).

The tool is imported lazily, at session end, not when this file is imported:
collection and the tests themselves run whether or not the tool exists. If it
is missing, the session fails and says so: a check that silently stops running
is worse than a red session.

Limit: only opens made by the pytest process itself are seen. A tool the suite
runs as a subprocess against the real tree is not; no suite does that today
(they run tools against fixture repositories).
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ci_changes.py"

#: Absolute paths the process opened, as seen at the time of the open. Resolved
#: once at session end, so the audit hook itself stays cheap.
_opened: set[str] = set()


def _record_open(event: str, args: tuple[Any, ...]) -> None:
    if event != "open" or not args:
        return
    path = args[0]
    if isinstance(path, bytes):
        path = os.fsdecode(path)
    if isinstance(path, str):  # an int is a file descriptor: no path to judge
        _opened.add(os.path.abspath(path))


# An audit hook cannot be removed; installed once per process.
if not getattr(sys, "_tin_engine_prose_hook", False):
    sys.addaudithook(_record_open)
    sys._tin_engine_prose_hook = True  # type: ignore[attr-defined]


def _load_tool() -> ModuleType | None:
    if not TOOL.exists():
        return None
    spec = importlib.util.spec_from_file_location("h11_ci_changes", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_verdict: list[str] = []


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    tool = _load_tool()
    if tool is None:
        _verdict.append(
            f"tools/ci_changes.py is missing from {ROOT}, so the check that no "
            "test reads a prose file (NOT_PROSE, h11 T5) could not run."
        )
    else:
        paths = sorted({os.path.realpath(p) for p in _opened})
        for name in tool.prose_reads(paths, ROOT):
            _verdict.append(
                f"The test session read {name}, a prose file: a pull request that "
                "changes only prose skips this suite. Add it to NOT_PROSE in "
                "tools/ci_changes.py, or stop reading it."
            )
    if _verdict and session.exitstatus == pytest.ExitCode.OK:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED


def pytest_terminal_summary(terminalreporter: Any) -> None:
    if _verdict:
        terminalreporter.section("prose files read by the suite (h11 T5)", red=True)
        for line in _verdict:
            terminalreporter.write_line(line, red=True)
