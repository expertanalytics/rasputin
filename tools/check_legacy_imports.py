"""Assert that the archived tree under legacy/ is internally consistent.

legacy/ is the porting reference: @migration-expert reads it to recover intent
before rebuilding a module on the new backend. It is excluded from ruff
(pyproject.toml), from mypy and from pytest's testpaths, so nothing else in the
project reads it at all -- an archive can rot silently and only be discovered
when someone needs it.

That is not hypothetical. The foundation reset archived most of src/rasputin/
by rename but deleted py2js.py and web_visualize.py outright, which left
legacy/rasputin/geometry.py importing a module that no longer existed and
transitively broke three files under legacy/tests/. No gate caught it.

This resolves every `rasputin.*` import in the tree against the modules
actually present. It parses rather than imports: the archived code needs CGAL
and a compiled extension that deliberately no longer build, so importing it is
not an option, and parsing is enough to catch a missing module.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEGACY = ROOT / "legacy"
PACKAGE = LEGACY / "rasputin"

# triangulate_dem is the compiled CGAL extension module (see legacy/bindings.cpp:
# PYBIND11_MODULE(triangulate_dem, m)). It was never a .py file, so it cannot be
# resolved against the source tree. Its absence is the point of the migration,
# not a defect in the archive.
COMPILED_MODULES = {"triangulate_dem"}


def available_modules() -> set[str]:
    """Names importable as `rasputin.<name>` from the archived package."""
    names = {p.stem for p in PACKAGE.glob("*.py") if p.stem != "__init__"}
    names |= {p.name for p in PACKAGE.iterdir() if p.is_dir() and (p / "__init__.py").exists()}
    # Names bound in __init__.py are importable as `from rasputin import X` even
    # though no module file carries them (rasputin_data_dir is one).
    init = PACKAGE / "__init__.py"
    if init.exists():
        tree = ast.parse(init.read_text(), str(init))
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign):
                names |= {t.id for t in node.targets if isinstance(t, ast.Name)}
            elif isinstance(node, ast.ImportFrom | ast.Import):
                names |= {(a.asname or a.name).split(".")[0] for a in node.names}
    return names | COMPILED_MODULES


def referenced(path: Path) -> set[tuple[str, int]]:
    """`rasputin.*` module names this file imports, with line numbers."""
    tree = ast.parse(path.read_text(), str(path))
    found: set[tuple[str, int]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            parts = node.module.split(".")
            if parts[0] != "rasputin":
                continue
            if len(parts) > 1:
                found.add((parts[1], node.lineno))
            else:
                # `from rasputin import x, y` -- each name must resolve.
                found |= {(a.name, node.lineno) for a in node.names}
        elif isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if parts[0] == "rasputin" and len(parts) > 1:
                    found.add((parts[1], node.lineno))
    return found


def main() -> int:
    if not PACKAGE.is_dir():
        print(f"error: {PACKAGE} does not exist", file=sys.stderr)
        return 1

    known = available_modules()
    broken: list[str] = []
    checked = 0

    for path in sorted(LEGACY.rglob("*.py")):
        checked += 1
        try:
            refs = referenced(path)
        except SyntaxError as exc:
            broken.append(f"{path.relative_to(ROOT)}:{exc.lineno}: does not parse: {exc.msg}")
            continue
        for name, lineno in sorted(refs):
            if name not in known:
                rel = path.relative_to(ROOT)
                broken.append(
                    f"{rel}:{lineno}: imports rasputin.{name}, which is not in the archive"
                )

    if broken:
        print(f"Archive integrity check FAILED ({len(broken)} unresolved):\n", file=sys.stderr)
        for line in broken:
            print(f"  {line}", file=sys.stderr)
        print(
            "\nlegacy/ is the porting reference and must stay self-consistent.\n"
            "Archive the missing module rather than deleting it, or record the drop\n"
            "in project_structure.md and remove the dangling import.",
            file=sys.stderr,
        )
        return 1

    print(f"Archive integrity OK: {checked} files parsed, all rasputin.* imports resolve.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
