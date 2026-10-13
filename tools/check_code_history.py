"""Check that docstrings, comments and CLI help text carry no project history.

The rule (`.claude/skills/python-development/SKILL.md`): docstrings, comments
and `help=` strings describe the code as it is. No increment numbers or labels,
no design-file IDs, no review or PR names, no record of who decided, no commit
hashes, and a pointer to a design file only once, in the module docstring, in
the form `` `<path>.md`, *<heading>* ``.

    python3 tools/check_code_history.py [PATH ...]
    python3 tools/check_code_history.py --same-code BASE

The first form reads every `.py` file under the given paths, or under
`[tool.check_code_history] paths` in `pyproject.toml` when none is given, and
prints `path:line: rule: text` for each match; it exits 1 when it printed any.
It reads comments, module, class and function docstrings, and the string
values of `help=`, `short_help=` and `epilog=` keyword arguments, nothing else.
There is no inline escape: a false positive is reworded.

The second form checks that every `.py` file differing from the commit `BASE`
differs only in prose: the same `ast.dump` once every bare string statement is
removed, and the same tool directives (`noqa`, `type: ignore`, `ruff:` ...) in
its comments, in order. It names each file that differs and exits 1.
"""

from __future__ import annotations

import argparse
import ast
import io
import re
import subprocess
import sys
import tokenize
import tomllib
from collections.abc import Iterator
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

RULES = {
    "increment": re.compile(r"\b[Ii]ncrements?\s+\d{1,2}[a-f]?\b"),
    # Not inside a path or a format spec (`.0f`); scientific notation is
    # dropped below.
    "label": re.compile(r"(?<![\w.$/-])(?:[1-9]|[1-4]\d)[a-f](?:-(?:\d|i{1,3}))?(?![\w.]|-\d)"),
    # D8 is also the hydrology flow-direction scheme, so it is not a design ID.
    "ruling-id": re.compile(r"(?<![\w-])(?!D8(?!\w))(?:[RDQBCSUJMFT]\d{1,2}|[RD]\d\.\d)(?!\w)"),
    "review": re.compile(r"\baudit PR\b|\bPR\s?(?:\d|[A-D])\b|\breview round\b|\bround \d\b"),
    "person": re.compile(r"\bOla\b|\bruled\b|\bruling\b"),
    "hash": re.compile(r"`[0-9a-f]{7,12}`"),
}
SCIENTIFIC = re.compile(r"\d+[eE]-?\d+")
# The two patterns of tools/check_citations.py, copied so neither gate depends
# on the other.
MD_PATH = r"(?<![\w/.-])(?!/)([\w./-]+\.md)"
MD = re.compile(MD_PATH + r"\b")
HEADING_CITATION = re.compile("`" + MD_PATH + r"`, \*([^*\n]+)\*")
HELP_KEYWORDS = {"help", "short_help", "epilog"}
# What a tool reads from a comment, matched anywhere in it and in any case
# (ruff obeys `NOQA`). A prefixed directive is compared to the comment's end;
# the prose after a `noqa` code list is not part of it.
DIRECTIVE = re.compile(
    r"\b(?:ruff|mypy|pyright|isort):.*"
    r"|type:\s*ignore\S*|noqa(?::\s*[A-Z]+\d+(?:,\s*[A-Z]+\d+)*)?|fmt:\s*\w+|pragma\S*",
    re.IGNORECASE,
)


def prose(source: str) -> Iterator[tuple[int, str, str]]:
    """Every comment, docstring line and help string, as (line, text, kind).

    The kind is `comment`, `module` (the module docstring), `docstring` or
    `help`. A help string's constant parts are joined and reported at the line
    where its value starts.
    """
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.COMMENT:
            yield token.start[0], token.string, "comment"
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Module | ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef):
            first = node.body[0] if node.body else None
            if first is not None and _bare_string(first):
                assert isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant)
                kind = "module" if isinstance(node, ast.Module) else "docstring"
                for offset, line in enumerate(str(first.value.value).splitlines()):
                    yield first.lineno + offset, line, kind
        if isinstance(node, ast.keyword) and node.arg in HELP_KEYWORDS:
            value = node.value
            parts = value.values if isinstance(value, ast.JoinedStr) else [value]
            text = "".join(
                p.value for p in parts if isinstance(p, ast.Constant) and isinstance(p.value, str)
            )
            if text:
                yield value.lineno, text, "help"


def matches(source: str) -> Iterator[tuple[int, str, str]]:
    """Every rule match in `source`'s prose, as (line, rule, matched text)."""
    allowed = 1  # the one heading-form pointer a module docstring may hold
    for number, line, kind in prose(source):
        for name, rule in RULES.items():
            for match in rule.finditer(line):
                if not SCIENTIFIC.fullmatch(match.group(0)):
                    yield number, name, match.group(0)
        heading = {m.start() + 1 for m in HEADING_CITATION.finditer(line)}
        for match in MD.finditer(line):
            if kind == "module" and match.start() in heading and allowed:
                allowed -= 1
            else:
                yield number, "pointer", match.group(0)


def scope(paths: list[str]) -> list[Path]:
    """The `.py` files under `paths`, or under the configured scope."""
    if not paths:
        config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
        paths = config.get("tool", {}).get("check_code_history", {}).get("paths", [])
    missing = [p for p in paths if not (ROOT / p).exists()]
    if missing:
        sys.exit(f"check_code_history: no such path: {', '.join(missing)}")
    files: list[Path] = []
    for p in paths:
        path = ROOT / p
        files += sorted(path.rglob("*.py")) if path.is_dir() else [path]
    return files


def code_shape(source: str) -> tuple[str, list[str]]:
    """The AST without bare strings, and the directives in the comments."""
    tree = ast.parse(source)
    for node in ast.walk(tree):
        for field in ("body", "orelse", "finalbody"):
            statements = getattr(node, field, None)
            if isinstance(statements, list):
                setattr(node, field, [s for s in statements if not _bare_string(s)])
    directives = [
        m.group(0).strip()
        for t in tokenize.generate_tokens(io.StringIO(source).readline)
        if t.type == tokenize.COMMENT
        for m in DIRECTIVE.finditer(t.string)
    ]
    return ast.dump(tree), directives


def _bare_string(statement: ast.stmt) -> bool:
    return (
        isinstance(statement, ast.Expr)
        and isinstance(statement.value, ast.Constant)
        and isinstance(statement.value.value, str)
    )


def _git(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)


def same_code(base: str) -> int:
    """Name every `.py` file whose code differs from `base`'s; 1 if any does."""
    changed = _git("diff", "--name-only", base, "--", "*.py").stdout.split()  # untracked: unseen
    differ = []
    for name in sorted(set(changed)):
        old = _git("show", f"{base}:{name}")
        new = ROOT / name
        if (
            old.returncode
            or not new.exists()
            or code_shape(old.stdout) != code_shape(new.read_text(encoding="utf-8"))
        ):
            differ.append(name)
    for name in differ:
        print(f"{name}: code differs from {base}")
    return 1 if differ else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("paths", nargs="*", help="files or directories to check")
    parser.add_argument("--same-code", metavar="BASE", help="check only prose changed")
    args = parser.parse_args()
    if args.same_code:
        return same_code(args.same_code)
    found = 0
    for path in scope(args.paths):
        shown = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
        for number, rule, text in matches(path.read_text(encoding="utf-8")):
            print(f"{shown}:{number}: {rule}: {text}")
            found += 1
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
