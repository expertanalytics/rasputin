"""Increment 42's measurements: history labels in docstrings and comments, and same-code.

Run from the repository root:

    python3 docs/increments/42-probes/history_probe.py scan src_python/tin_engine
    python3 docs/increments/42-probes/history_probe.py same OLD.py NEW.py

`scan` prints one line per match (file:line, rule, matched text, the line) and a
total on stderr. `same` exits 0 when the two files have the same code: the same
`ast.dump` once every bare string statement (docstrings included) is removed, and
the same tool directives (`type:`, `noqa`, `fmt:`, `pragma`) in their comments.
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
DIRECTIVE = re.compile(
    r"#\s*(type:\s*ignore\S*|noqa(?::\s*[A-Z]+\d+(?:,\s*[A-Z]+\d+)*)?|fmt:\s*\w+|pragma\S*)"
)


def prose(source: str) -> Iterator[tuple[int, str]]:
    """Every comment, and every line of every module, class or function docstring."""
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.COMMENT:
            yield token.start[0], token.string
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Module | ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef):
            first = node.body[0] if node.body else None
            if (
                isinstance(first, ast.Expr)
                and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)
            ):
                for offset, line in enumerate(first.value.value.splitlines()):
                    yield first.lineno + offset, line


def scan(root: Path) -> int:
    counts: Counter[str] = Counter()
    for path in sorted(root.rglob("*.py")):
        for number, line in prose(path.read_text(encoding="utf-8")):
            for name, rule in RULES.items():
                for match in rule.finditer(line):
                    if SCIENTIFIC.fullmatch(match.group(0)):
                        continue
                    counts[name] += 1
                    print(f"{path}:{number}\t{name}\t{match.group(0)}\t{line.strip()[:120]}")
    print(f"TOTAL {sum(counts.values())} {dict(counts)}", file=sys.stderr)
    return 0


def code_shape(source: str) -> tuple[str, list[str]]:
    tree = ast.parse(source)
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if isinstance(body, list):
            node.body = [  # type: ignore[attr-defined]
                s
                for s in body
                if not (
                    isinstance(s, ast.Expr)
                    and isinstance(s.value, ast.Constant)
                    and isinstance(s.value.value, str)
                )
            ]
    directives = [
        m.group(1)
        for t in tokenize.generate_tokens(io.StringIO(source).readline)
        if t.type == tokenize.COMMENT
        for m in [DIRECTIVE.search(t.string)]
        if m
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
