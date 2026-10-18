"""The whole `tin_engine` package carries no project history in its docstrings or comments.

`tools/check_code_history.py` checks only the modules listed under
`[tool.check_code_history]` in `pyproject.toml`. This test runs it over the
whole package, so the modules outside that list are held to the same rule
until the list covers the package.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

# A harness test: it runs a tool and imports nothing from tin_engine.
pytestmark = pytest.mark.harness

ROOT = Path(__file__).resolve().parents[2]


def test_the_history_gate_finds_nothing_in_the_whole_package() -> None:
    run = subprocess.run(
        [sys.executable, "tools/check_code_history.py", "src_python/tin_engine"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert run.stdout == "", run.stdout
    assert run.returncode == 0
