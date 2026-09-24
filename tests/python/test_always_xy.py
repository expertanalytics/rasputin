"""Ruling 8 (`11-raster-ingestion.md` §9): `always_xy=True`, enforced rather than asserted.

Every `pyproj.Transformer.from_crs` in `src_python/` passes `always_xy=True`
as a literal, and no CRS there is spelled `+init=`. Without `always_xy`, the
same transform puts a point 6 000 km away and says nothing (prior art §5.1).

Today `src_python/` has no transformer at all, so the check against the real
tree passes vacuously. Per `PRINCIPLES.md` A3 the checker is therefore also run
against planted violations in `tmp_path`, which it must catch; and a
compliant plant, which it must pass, so it cannot pass by flagging everything.
It was also shown failing against a violation planted inside `src_python/`
itself, in the commit that added it (see the commit message).

The design says "grep". This is a grep for `+init=` and an AST walk for the
call sites, because a call's keywords can span lines and a line grep cannot
tell `always_xy=True` from `always_xy=False` or from a comment. Any mention of
`from_crs` that is not a direct call with a literal `always_xy=True` is a
violation, which catches aliasing (`f = Transformer.from_crs`) and `**kwargs`.
`legacy/` is out of scope: frozen, and every CRS there is `+init=`.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SRC_PYTHON = Path(__file__).resolve().parents[2] / "src_python"


def _is_literal_true(node: ast.expr | None) -> bool:
    return isinstance(node, ast.Constant) and node.value is True


def _call_complies(call: ast.Call) -> bool:
    keywords = {k.arg: k.value for k in call.keywords}
    return None not in keywords and _is_literal_true(keywords.get("always_xy"))


def violations(root: Path) -> list[str]:
    """Every breach of ruling 8 under `root`, as `path:line: reason`."""
    found: list[str] = []
    for path in sorted([*root.rglob("*.py"), *root.rglob("*.pyi")]):
        text = path.read_text(encoding="utf-8")
        for number, line in enumerate(text.splitlines(), start=1):
            if "+init=" in line:
                found.append(f"{path}:{number}: '+init=' CRS spelling")
        tree = ast.parse(text, filename=str(path))
        called: set[int] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                callee = node.func
                name = callee.attr if isinstance(callee, ast.Attribute) else None
                name = callee.id if isinstance(callee, ast.Name) else name
                if name == "from_crs":
                    called.add(id(callee))
                    if not _call_complies(node):
                        found.append(f"{path}:{node.lineno}: from_crs without always_xy=True")
        for node in ast.walk(tree):
            mentioned = (isinstance(node, ast.Attribute) and node.attr == "from_crs") or (
                isinstance(node, ast.Name) and node.id == "from_crs"
            )
            if mentioned and id(node) not in called:
                found.append(f"{path}:{node.lineno}: from_crs referenced, not called")
    return found


def test_src_python_obeys_ruling_8() -> None:
    assert SRC_PYTHON.is_dir()
    assert violations(SRC_PYTHON) == []


COMPLIANT = """
from pyproj import Transformer
t = Transformer.from_crs("EPSG:4326", "EPSG:25833", always_xy=True)
u = Transformer.from_crs(
    "EPSG:25833",
    "EPSG:4326",
    always_xy=True,
)
"""


@pytest.mark.parametrize(
    "planted",
    [
        'Transformer.from_crs("EPSG:4326", "EPSG:25833")',
        'Transformer.from_crs("EPSG:4326", "EPSG:25833", always_xy=False)',
        'Transformer.from_crs("EPSG:4326", "EPSG:25833", always_xy=flag)',
        'Transformer.from_crs("EPSG:4326", "EPSG:25833", **options)',
        'from_crs("EPSG:4326", "EPSG:25833")',
        "make = Transformer.from_crs",
        'crs = CRS("+init=epsg:32633")',
    ],
    ids=[
        "missing",
        "false",
        "not_a_literal",
        "kwargs_splat",
        "bare_name",
        "aliased",
        "init_spelling",
    ],
)
def test_checker_catches_a_planted_violation(tmp_path: Path, planted: str) -> None:
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "ok.py").write_text(COMPLIANT, encoding="utf-8")
    (tmp_path / "pkg" / "bad.py").write_text(f"{planted}\n", encoding="utf-8")
    found = violations(tmp_path)
    assert len(found) == 1, found
    assert "bad.py:1:" in found[0]


def test_checker_passes_compliant_calls(tmp_path: Path) -> None:
    (tmp_path / "ok.py").write_text(COMPLIANT, encoding="utf-8")
    assert violations(tmp_path) == []
