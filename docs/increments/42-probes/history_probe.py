"""Increment 42's measurements: history labels in docstrings and comments, and same-code.

Run from the repository root:

    python3 docs/increments/42-probes/history_probe.py scan src_python/tin_engine
    python3 docs/increments/42-probes/history_probe.py same OLD.py NEW.py

`scan` prints one line per match (file:line, rule, matched text, the line) and a
total on stderr. It reads comments, module, class and function docstrings, and
the string values of `help=`, `short_help=` and `epilog=` keyword arguments.
`same` exits 0 when the two files have the same code: the same `ast.dump` once
every bare string statement (docstrings included) is removed from every
statement list (`body`, `orelse`, `finalbody`), and the same tool directives in
their comments, in order. A directive is any comment holding `noqa`, `type:`,
`fmt:`, `pragma`, or a `ruff:`, `mypy:`, `pyright:` or `isort:` prefix.
This is the prototype of the design's checker, not the checker.
"""

from __future__ import annotations

import ast
import io
import re
import sys
import tokenize
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

RULES = {
    "increment": re.compile(r"\b[Ii]ncrements?\s+\d{1,2}[a-f]?\b"),
    "label": re.compile(r"(?<![\w.$/-])(?:[1-9]|[1-4]\d)[a-f](?:-(?:\d|i{1,3}))?(?![\w.]|-\d)"),
    "ruling-id": re.compile(r"(?<![\w-])(?:[RDQBCSUJMFT]\d{1,2}|[RD]\d\.\d)(?!\w)"),
    "review": re.compile(r"\baudit PR\b|\bPR\s?(?:\d|[A-D])\b|\breview round\b|\bround \d\b"),
    "person": re.compile(r"\bOla\b|\bruled\b|\bruling\b"),
    "hash": re.compile(r"`[0-9a-f]{7,12}`"),
}
SCIENTIFIC = re.compile(r"\d+[eE]-?\d+")
# A directive: what a tool reads from a comment, matched anywhere in it, so a
# file-level ruff suppression counts as well as an inline one. A prefixed one
# (`ruff:`, `mypy:`, `pyright:`, `isort:`) is compared whole, to the comment's
# end; after a `noqa` code list the prose is not part of it.
DIRECTIVE = re.compile(
    r"\b(?:ruff|mypy|pyright|isort):.*"
    r"|type:\s*ignore\S*|noqa(?::\s*[A-Z]+\d+(?:,\s*[A-Z]+\d+)*)?|fmt:\s*\w+|pragma\S*"
)
MD = re.compile(r"(?<![\w/.-])(?!/)[\w./-]+\.md\b")
HEADING_POINTER = re.compile(r"`(?<![\w/.-])(?!/)[\w./-]+\.md`, \*[^*\n]+\*")
HELP_KEYWORDS = {"help", "short_help", "epilog"}


def prose(source: str) -> Iterator[tuple[int, str, str]]:
    """Every comment, docstring line and help string, as (line, text, kind).

    The kind is `comment`, `module` (the module docstring), `docstring` (a class
    or function docstring) or `help`.
    """
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.COMMENT:
            yield token.start[0], token.string, "comment"
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Module | ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef):
            first = node.body[0] if node.body else None
            if (
                isinstance(first, ast.Expr)
                and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)
            ):
                kind = "module" if isinstance(node, ast.Module) else "docstring"
                for offset, line in enumerate(first.value.value.splitlines()):
                    yield first.lineno + offset, line, kind
        if isinstance(node, ast.keyword) and node.arg in HELP_KEYWORDS:
            parts = node.value.values if isinstance(node.value, ast.JoinedStr) else [node.value]
            text = "".join(
                p.value for p in parts if isinstance(p, ast.Constant) and isinstance(p.value, str)
            )
            if text:
                yield node.value.lineno, text, "help"


def pointers(lines: list[tuple[int, str, str]]) -> Iterator[tuple[int, str, str]]:
    """The `.md` mentions the rule does not allow.

    Allowed: one per module docstring, in the heading form. Anything else is
    reported: a second one, one not in the heading form, one in a comment, a
    class or function docstring, or a help string.
    """
    allowed = 1
    for number, line, kind in lines:
        heading = {m.start() for m in HEADING_POINTER.finditer(line)}
        for match in MD.finditer(line):
            if kind == "module" and match.start() - 1 in heading and allowed:
                allowed -= 1
                continue
            yield number, match.group(0), line


def one_line(text: str) -> str:
    return " ".join(text.split())[:120]


def scan(root: Path) -> int:
    counts: Counter[str] = Counter()
    for path in sorted(root.rglob("*.py")):
        lines = list(prose(path.read_text(encoding="utf-8")))
        for number, line, _ in lines:
            for name, rule in RULES.items():
                for match in rule.finditer(line):
                    if SCIENTIFIC.fullmatch(match.group(0)):
                        continue
                    counts[name] += 1
                    print(f"{path}:{number}\t{name}\t{match.group(0)}\t{one_line(line)}")
        for number, text, line in pointers(lines):
            counts["pointer"] += 1
            print(f"{path}:{number}\tpointer\t{text}\t{one_line(line)}")
    print(f"TOTAL {sum(counts.values())} {dict(counts)}", file=sys.stderr)
    return 0


def code_shape(source: str) -> tuple[str, list[str]]:
    tree = ast.parse(source)
    for node in ast.walk(tree):
        for field in ("body", "orelse", "finalbody"):
            statements = getattr(node, field, None)
            if isinstance(statements, list):
                setattr(
                    node,
                    field,
                    [
                        s
                        for s in statements
                        if not (
                            isinstance(s, ast.Expr)
                            and isinstance(s.value, ast.Constant)
                            and isinstance(s.value.value, str)
                        )
                    ],
                )
    directives = [
        m.group(0).strip()
        for t in tokenize.generate_tokens(io.StringIO(source).readline)
        if t.type == tokenize.COMMENT
        for m in DIRECTIVE.finditer(t.string)
    ]
    return ast.dump(tree), directives


def same(old: Path, new: Path) -> int:
    equal = code_shape(old.read_text(encoding="utf-8")) == code_shape(
        new.read_text(encoding="utf-8")
    )
    print("same code" if equal else "CODE DIFFERS")
    return 0 if equal else 1


if __name__ == "__main__":
    if sys.argv[1] == "scan":
        sys.exit(scan(Path(sys.argv[2])))
    sys.exit(same(Path(sys.argv[2]), Path(sys.argv[3])))
