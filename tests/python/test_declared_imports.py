"""Every third-party module the package imports is a requirement it declares.

`docs/increments/48-declared-imports.md`, *6. The tests that pin it*, T5.
A module imported by name must be installed by name: one that arrives only
through another package's own requirements (deptry's DEP003) breaks the day
that package drops it, and one declared nowhere (DEP001) breaks a fresh
install at once. Files under `tin_engine/mcp/` may also use the `mcp` extra;
every other file only the core dependencies.

Static: reads `pyproject.toml` and the source text, imports nothing from the
package and needs no extra installed. Import names and requirement names are
compared after PEP 503 normalisation (lower case, runs of `-`, `_` and `.`
to one `-`).
"""

from __future__ import annotations

import ast
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = ROOT / "src_python" / "tin_engine"
REQUIREMENT_NAME = re.compile(r"[A-Za-z0-9._-]+")


def normalised(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def requirement_names(requirements: list[str]) -> set[str]:
    found = (REQUIREMENT_NAME.match(r.strip()) for r in requirements)
    return {normalised(m.group()) for m in found if m}


def third_party_imports(path: Path) -> list[tuple[int, str]]:
    """(line, top-level name) of every absolute import outside the standard library."""
    found: list[tuple[int, str]] = []
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"), filename=str(path))):
        if isinstance(node, ast.Import):
            names = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names = [node.module]
        else:
            continue
        for name in names:
            top = name.partition(".")[0]
            if top not in sys.stdlib_module_names and top != "tin_engine":
                found.append((node.lineno, top))
    return found


def test_every_third_party_import_is_declared() -> None:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    core = requirement_names(project["dependencies"])
    with_mcp = core | requirement_names(project["optional-dependencies"]["mcp"])
    files = sorted(PACKAGE.rglob("*.py"))
    assert len(files) >= 60, f"only {len(files)} files under {PACKAGE}: is the walk rooted right?"
    offenders = [
        f"{path.relative_to(PACKAGE)}:{line} {name}"
        for path in files
        for line, name in third_party_imports(path)
        if normalised(name)
        not in (with_mcp if path.relative_to(PACKAGE).parts[0] == "mcp" else core)
    ]
    assert offenders == [], "imported but not declared:\n" + "\n".join(offenders)
