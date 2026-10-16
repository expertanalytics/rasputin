#!/usr/bin/env python3
"""Rule text in words, with the change since the last retrospective.

Spec: docs/increments/h8-window-and-recap.md §3.6 and §8 (points 13, 14, 19).
"Now" is the working tree: the files on disk that a session loads, tracked or
not, so an uncommitted rule edit shows. The reference is derived from git, never
stored: the newest commit reachable from HEAD that added a dated retrospective
file. Its counts come from `git ls-tree` at that commit, so a file removed since
is still found. When git fails, the current counts print without a reference.

Usage: python3 tools/rule_sizes.py
"""

from __future__ import annotations

import subprocess
from fnmatch import fnmatch
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RULE_FILES = (
    "CLAUDE.md",
    ".claude/REQUIRED-READING.md",
    "docs/increments/README.md",
    "docs/PRINCIPLES.md",
    ".claude/briefs/common.md",  # h9: the brief template
)
RETROSPECTIVES = "docs/retrospectives/2???-??-??-*.md"
TREES = (*RULE_FILES, ".claude/agents", ".claude/skills")


def words(text: str) -> int:
    """What `wc -w` counts for ASCII text."""
    return len(text.split())


def in_set(path: str) -> bool:
    """A rule file, an agent file, or any Markdown file under the skills."""
    if path in RULE_FILES or (path.startswith(".claude/skills/") and path.endswith(".md")):
        return True
    return fnmatch(path, ".claude/agents/*.md") and path.count("/") == 2


def ordered(paths: list[str]) -> list[str]:
    """Rule files in RULE_FILES order, then agents, then skills, each group sorted."""

    def key(p: str) -> tuple[int, int, str]:
        if p in RULE_FILES:
            return 0, RULE_FILES.index(p), ""
        return 1 + p.startswith(".claude/skills/"), 0, p

    return sorted(set(paths), key=key)


def current(root: Path) -> dict[str, int]:
    """Counts of the file set as it is on disk."""
    found = [p.relative_to(root).as_posix() for name in RULE_FILES if (p := root / name).is_file()]
    for pattern in (".claude/agents/*.md", ".claude/skills/**/*.md"):
        found += [p.relative_to(root).as_posix() for p in root.glob(pattern) if p.is_file()]
    return {p: words((root / p).read_text(errors="replace")) for p in ordered(found)}


def _git(root: Path, *args: str) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(root), *args], capture_output=True, text=True, check=False
    )
    return result.stdout if result.returncode == 0 else None


def reference(root: Path) -> tuple[str, dict[str, int]] | str:
    """(label, counts at the last retrospective), or why there is none."""
    log = _git(
        root, "log", "-1", "--diff-filter=A", "--format=%h", "--name-only", "--", RETROSPECTIVES
    )
    if log is None:
        return "no reference: git failed"
    if not log.strip():
        return "no retrospective found"
    rev, *added = log.split()
    listed = _git(root, "ls-tree", "-r", "--name-only", rev, "--", *TREES)
    if listed is None:
        return "no reference: git failed"
    then = {}
    for path in ordered([p for p in listed.splitlines() if in_set(p)]):
        then[path] = words(_git(root, "show", f"{rev}:{path}") or "")
    return f"{Path(max(added)).name} ({rev})", then


def _change(delta: int) -> str:
    return f"{delta:+d}" if delta else "0"


def table(now: dict[str, int], then: dict[str, int] | None, label: str) -> list[str]:
    """The table of §3.6: one line per file, then the total."""
    total = sum(now.values())
    if then is None:
        rows = [f"  {n:5} {path}" for path, n in now.items()]
        return ["Rule text in words (no retrospective found):", *rows, f"  {total:5} total"]
    lines = [f"Rule text in words, change since {label}:"]
    for path, n in now.items():
        lines.append(f"  {n:5} {_change(n - then[path]) if path in then else 'new':>6} {path}")
    lines += [f"  {0:5} removed (was {n}) {path}" for path, n in then.items() if path not in now]
    return [*lines, f"  {total:5} {_change(total - sum(then.values())):>6} total"]


def report(root: Path) -> list[str]:
    """The table for the checkout at `root`."""
    now, found = current(root), reference(root)
    if isinstance(found, tuple):
        return table(now, found[1], found[0])
    return [f"Rule text in words ({found}):", *table(now, None, "")[1:]]


def main() -> int:
    print("\n".join(report(ROOT)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
