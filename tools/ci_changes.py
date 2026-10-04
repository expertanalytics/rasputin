"""Decide whether a change needs the code jobs in CI (harness increment h11).

`docs/increments/h11-ci-path-filter.md` §2.1 is the specification. Deny by
default: a path is prose only if it is a Markdown file at the repository root
or under `docs/`, and not one the build or the suites read (`NOT_PROSE`).
Everything else is code. A wrong "code" costs one full run; a wrong "prose"
would let a broken change land green.

Usage, in the workflow's `changes` job:

    python3 tools/ci_changes.py <base-sha> <head-sha> >> "$GITHUB_OUTPUT"

It prints exactly one line, `code=true` or `code=false`, and exits 0. Whenever
it cannot tell (no base, a git failure, an empty diff) it prints `code=true`,
the full run, and says why on stderr.

Stdlib only: the governance job runs it with no install step.
"""

from __future__ import annotations

import subprocess
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path, PurePath

#: Markdown files the build or the suites read, so changing one is a code change.
#: `README.md`: pyproject.toml's `readme`, so `pip install -e`. `CLAUDE.md`:
#: test_spawn_rule_text.py and test_session_state.py. The two `docs/` files:
#: test_session_state.py copies them. The T5 hook in tests/python/conftest.py
#: fails the suite on a read of any other prose file.
NOT_PROSE: frozenset[str] = frozenset(
    {"README.md", "CLAUDE.md", "docs/increments/README.md", "docs/PRINCIPLES.md"}
)


def is_prose(path: str) -> bool:
    """True for a root or `docs/` Markdown file not in `NOT_PROSE` (case-sensitive)."""
    if not path.endswith(".md") or path in NOT_PROSE:
        return False
    return "/" not in path or path.startswith("docs/")


def needs_full_ci(paths: Sequence[str]) -> bool:
    """True unless every path is prose; an empty list is a full run."""
    return not paths or not all(is_prose(path) for path in paths)


def changed_paths(base: str, head: str, repo: Path) -> list[str]:
    """Paths changed between two commits; raises on any git failure.

    `--no-renames` lists both sides of a rename, so moving a source file to a
    prose path still counts as a code change. `-z` keeps unusual names unquoted.
    """
    done = subprocess.run(
        ["git", "diff", "--name-only", "-z", "--no-renames", base, head],
        cwd=repo,
        capture_output=True,
        check=True,
    )
    return [name for name in done.stdout.decode().split("\0") if name]


def prose_reads(paths: Iterable[str], root: Path) -> list[str]:
    """The paths inside `root` that are prose, `/`-separated relative to `root`.

    Pure path arithmetic; paths outside `root` are dropped. The T5 hook calls it
    with every file the test session opened.
    """
    found = []
    for path in paths:
        pure = PurePath(path)
        if pure.is_relative_to(root):
            relative = pure.relative_to(root).as_posix()
            if is_prose(relative):
                found.append(relative)
    return found


def _decide(base: str, head: str) -> tuple[bool, str]:
    """Whether to run the code jobs, and why if it is a fallback."""
    if not base:
        return True, "no base commit (a hand-started run): running every job."
    try:
        paths = changed_paths(base, head, Path.cwd())
    except (OSError, subprocess.CalledProcessError) as error:
        return True, f"git could not list the changes ({error}): running every job."
    if not paths:
        return True, "git found no changed files: running every job."
    return needs_full_ci(paths), ""


def main(argv: list[str]) -> int:
    """Print `code=true` or `code=false` for `<base> <head>`; always exit 0."""
    if len(argv) != 2:
        code, why = True, "expected two arguments, <base> <head>: running every job."
    else:
        code, why = _decide(argv[0], argv[1])
    if why:
        print(f"ci_changes: {why}", file=sys.stderr)
    print(f"code={'true' if code else 'false'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
