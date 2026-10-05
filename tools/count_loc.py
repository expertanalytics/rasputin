#!/usr/bin/env python3
"""The one counter for the production-line ceiling of `CLAUDE.md` §2 (h16 T1).

    python3 tools/count_loc.py <base> [<head>]

Counts the lines `git diff <base>...<head>` adds and removes, under §2's
exclusions: a line counts unless it is blank, a comment, a docstring, or
inside a raw literal's body. Added lines are judged in the file at `<head>`,
removed lines at the merge base. Prints one tab-separated row per counted file
(path, added, removed, net), a `total` row, then one `not counted:` line per
changed file it skipped, with the reason. Committed revisions only.

Spec: docs/increments/h16-harness-fixes.md §2, T1. Exit 2 with one
`count_loc:` line on stderr on a git error.
"""

from __future__ import annotations

import ast
import io
import os
import re
import subprocess
import sys
import tokenize
from typing import Literal

Kind = Literal["python", "c++", "cmake", "shell"]
Row = tuple[str, int, int, int]

#: User and system config, replace refs: none may change what is counted.
GIT_ENV = {"GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1",
           "GIT_NO_REPLACE_OBJECTS": "1"}  # fmt: skip
#: Config given in the environment (`GIT_CONFIG_COUNT` with its `_KEY_n`/`_VALUE_n`,
#: `GIT_CONFIG_PARAMETERS`): dropped from what git sees.
ENV_CONFIG = re.compile(r"GIT_CONFIG_(?:COUNT|PARAMETERS|KEY_\d+|VALUE_\d+)$")
#: Explicit prefixes, hunk context, `--no-relative` and myers (git's default),
#: so the repository's own config (`diff.noprefix`, `diff.srcPrefix`,
#: `diff.interHunkContext`, `diff.relative`, `diff.algorithm`) cannot move the
#: headers, fuse hunks, narrow the diff to the directory the counter runs from,
#: or pair the lines differently.
DIFF = ("diff", "-M", "--no-ext-diff", "--no-color", "--src-prefix=a/", "--dst-prefix=b/",
        "--inter-hunk-context=0", "--no-relative", "--diff-algorithm=myers")  # fmt: skip
SKIPPED_DIRS = ("tests/", "docs/")
SUFFIXES: dict[str, Kind] = {
    ".py": "python", ".pyi": "python", ".h": "c++", ".hpp": "c++", ".cpp": "c++",
    ".cc": "c++", ".cxx": "c++", ".cmake": "cmake", ".sh": "shell",
}  # fmt: skip
#: The prefix of a raw string literal, `R`, `LR`, `uR`, `UR` or `u8R`, not inside a name.
RAW_PREFIX = re.compile(r"[^A-Za-z0-9_](?:u8|[LuU])?R$")
HUNK = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+(\d+)(?:,(\d+))? @@")
NOT_CODE = {tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT,
            tokenize.ENDMARKER, tokenize.ENCODING}  # fmt: skip
#: Tokens opening and closing an f-string (3.12) or t-string (3.14), which span lines.
OPENS = {getattr(tokenize, n, -1) for n in ("FSTRING_START", "TSTRING_START")}
CLOSES = {getattr(tokenize, n, -1) for n in ("FSTRING_END", "TSTRING_END")}


def kind_of(path: str) -> Kind | None:
    """The kind of a counted file, or None for a file of no kind (Markdown, data)."""
    name = path.rsplit("/", 1)[-1]
    if name == "CMakeLists.txt":
        return "cmake"
    _, dot, suffix = name.rpartition(".")
    return SUFFIXES.get(f".{suffix}") if dot else None


def _counted(path: str) -> Kind | None:
    """The kind of `path` when the counter counts it: of a kind, outside `SKIPPED_DIRS`."""
    return None if path.startswith(SKIPPED_DIRS) else kind_of(path)


def _docstring_spans(text: str) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    """(start, end) positions of every module, class and function docstring."""
    spans = []
    for node in ast.walk(ast.parse(text)):
        scoped = (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)
        if isinstance(node, scoped) and node.body:
            first = node.body[0]
            if (isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant)
                    and isinstance(first.value.value, str)):  # fmt: skip
                end = (first.end_lineno or first.lineno, first.end_col_offset or 0)
                spans.append(((first.lineno, first.col_offset), end))
    return spans


def _python_lines(text: str, name: str) -> frozenset[int]:
    """Lines holding a token other than comments and layout; a docstring counts
    on none of its lines, any other string (f-strings too) on its first only.
    Text that does not tokenize or parse counts every non-blank line, with a
    warning naming `name` on stderr."""
    try:
        spans = _docstring_spans(text)
        tokens = list(tokenize.generate_tokens(io.StringIO(text).readline))
    except (SyntaxError, tokenize.TokenError, ValueError) as exc:
        print(f"count_loc: {name} does not tokenize ({exc}); counting its non-blank lines",
              file=sys.stderr)  # fmt: skip
        return frozenset(n for n, s in enumerate(text.splitlines(), 1) if s.strip())
    out: set[int] = set()
    depth = 0  # inside an f-string: its later lines are part of one literal
    for tok in tokens:
        if tok.type in NOT_CODE:
            continue
        inside_doc = tok.type == tokenize.STRING and any(
            start <= tok.start and tok.end <= end for start, end in spans
        )
        if not inside_doc and depth == 0:
            out.add(tok.start[0])
        depth += (tok.type in OPENS) - (tok.type in CLOSES)
    return frozenset(out)


def _raw_open(text: str, i: int) -> int:
    """At a `"` that opens a raw string `R"d(`, the index of its `(`; else -1."""
    window = (" " if i < 4 else "") + text[max(0, i - 4) : i]
    if not RAW_PREFIX.search(window):
        return -1
    paren = text.find("(", i)
    return paren if 0 <= paren - i - 1 <= 16 and "\n" not in text[i:paren] else -1


def _cpp_lines(text: str) -> frozenset[int]:
    """Lines with a character outside comments and outside a raw string's body."""
    out: set[int] = set()
    line, i, n = 1, 0, len(text)
    while i < n:
        c, start = text[i], i
        if c == "\n":
            line, i = line + 1, i + 1
            continue
        if text.startswith("//", i):
            end = text.find("\n", i)
            i = n if end < 0 else end
            continue
        if text.startswith("/*", i):
            end = text.find("*/", i + 2)
            i = n if end < 0 else end + 2
        elif c == '"' and (paren := _raw_open(text, i)) >= 0:
            out.add(line)  # the opening line counts; the body and `)d"` do not
            close = ")" + text[i + 1 : paren] + '"'
            end = text.find(close, paren)
            i = n if end < 0 else end + len(close)
        elif c in "\"'" and not (c == "'" and i and (text[i - 1].isalnum() or text[i - 1] == "_")):
            out.add(line)  # an ordinary literal: comment markers inside it are text
            i += 1
            while i < n and text[i] not in (c, "\n"):
                i += 2 if text[i] == "\\" else 1
            i += 1 if i < n and text[i] == c else 0
        else:
            if not c.isspace():
                out.add(line)
            i += 1
        line += text.count("\n", start, i)
    return frozenset(out)


def counted_lines(text: str, kind: Kind, name: str = "<text>") -> frozenset[int]:
    """The 1-based lines of `text` that count under §2; `name` is for warnings."""
    if kind == "python":
        return _python_lines(text, name)
    if kind == "c++":
        return _cpp_lines(text)
    rows = enumerate(text.splitlines(), 1)
    return frozenset(n for n, s in rows if s.strip() and not s.lstrip().startswith("#"))


def hunks(diff: str) -> dict[str, tuple[list[int], list[int]]]:
    """From `git diff -U0` text: path -> (added lines at the new side, removed
    lines at the old). The path is the new one, or the old for a deletion.
    Only the file headers and the `@@` lines are read, so a changed line that
    looks like a header (`+++i;`) is never taken for one."""
    out: dict[str, tuple[list[int], list[int]]] = {}
    old = new = ""
    in_header = False
    for text in diff.splitlines():
        if text.startswith("diff --git "):
            in_header, old, new = True, "", ""
        elif in_header and text.startswith("--- "):
            old = text[4:].removeprefix("a/")
        elif in_header and text.startswith("+++ "):
            new = text[4:].removeprefix("b/")
        elif (match := HUNK.match(text)) is not None:
            in_header = False
            path = old if new == "/dev/null" else new
            added, removed = out.setdefault(path, ([], []))
            a, b, c, d = (int(g) if g is not None else 1 for g in match.groups())
            removed.extend(range(a, a + b))
            added.extend(range(c, c + d))
    return out


def tally(
    changes: list[tuple[str, str, str]],
    lines: dict[str, tuple[list[int], list[int]]],
    read: dict[tuple[str, str], frozenset[int]],
) -> tuple[list[Row], list[str]]:
    """Join the three: `changes` holds (status, old path, new path) per changed
    file, `lines` is `hunks`' output, and `read[(side, path)]` the counted lines
    of a file at side `old` or `new`. Returns the rows and the skipped lines."""
    rows: list[Row] = []
    skipped: list[str] = []
    for status, old, new in changes:
        if status == "R" and (_counted(old) is None) != (_counted(new) is None):
            # Across the counted boundary: in full, as an added or a deleted file.
            inward = _counted(new) is not None
            n = len(read.get(("new", new) if inward else ("old", old), frozenset()))
            rows.append((new, n, 0, n) if inward else (old, 0, n, -n))
            continue
        path = old if status == "D" else new
        why = next((d for d in SKIPPED_DIRS if path.startswith(d)), None)
        if why is None and kind_of(path) is None:
            why = "not code"
        if why is not None:
            skipped.append(f"not counted: {path} ({why})")
            continue
        added, removed = lines.get(path, ([], []))
        plus = sum(n in read.get(("new", new), frozenset()) for n in added)
        minus = sum(n in read.get(("old", old), frozenset()) for n in removed)
        rows.append((path, plus, minus, plus - minus))
    return sorted(rows), sorted(skipped)


def _env() -> dict[str, str]:
    """This process's environment, less config given in it, plus `GIT_ENV`."""
    return {**{k: v for k, v in os.environ.items() if not ENV_CONFIG.match(k)}, **GIT_ENV}


def _git(*args: str) -> str:
    done = subprocess.run(["git", *args], capture_output=True, text=True, errors="replace",
                          env=_env(), check=False)  # fmt: skip
    if done.returncode != 0:
        first = (done.stderr.strip().splitlines() or ["failed"])[0]
        raise RuntimeError(f"git {args[0]}: {first}")
    return done.stdout


def _read(rev: str, path: str, kind: Kind) -> frozenset[int]:
    return counted_lines(_git("cat-file", "blob", f"{rev}:{path}"), kind, f"{path} at {rev[:12]}")


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) not in (1, 2):
        print("count_loc: usage: count_loc.py <base> [<head>]", file=sys.stderr)
        return 2
    try:
        given = args[1] if len(args) == 2 else "HEAD"
        head = _git("rev-parse", "--verify", "--end-of-options", f"{given}^{{commit}}").strip()
        base = _git("merge-base", args[0], head).strip()
        fields = _git(*DIFF, "--name-status", "-z", base, head).split("\0")
        changes: list[tuple[str, str, str]] = []
        while len(fields) > 1:
            status = fields.pop(0)[0]
            old = fields.pop(0)
            new = fields.pop(0) if status in "RC" else old
            changes.append((status, old, new))
        lines = hunks(_git("-c", "core.quotePath=false", *DIFF, "-U0", base, head))
        read: dict[tuple[str, str], frozenset[int]] = {}
        for status, old, new in changes:
            if status != "A" and (kind := _counted(old)) is not None:
                read[("old", old)] = _read(base, old, kind)
            if status != "D" and (kind := _counted(new)) is not None:
                read[("new", new)] = _read(head, new, kind)
    except RuntimeError as exc:
        print(f"count_loc: {exc}", file=sys.stderr)
        return 2
    rows, skipped = tally(changes, lines, read)
    for path, plus, minus, net in rows:
        print(f"{path}\t{plus}\t{minus}\t{net}")
    totals = (sum(r[1] for r in rows), sum(r[2] for r in rows), sum(r[3] for r in rows))
    print("total\t" + "\t".join(map(str, totals)))
    for line in skipped:
        print(line)
    return 0


if __name__ == "__main__":
    sys.exit(main())
