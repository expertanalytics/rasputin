"""The C++ core's comments carry no project history.

The rule is the one `tools/check_code_history.py` holds the Python package to,
applied to the `//` and `/* */` comments of every `.hpp`, `.cpp` and `.h` file
under `include/`, `bindings/` and `src/`: no increment numbers or labels, no
design-file IDs, no review names, no record of who decided, no commit hashes.
A pointer to a design file is allowed once per file, in the file's opening
comment (the comments before the first line of code, `#pragma once` aside),
and only in the heading form `` `<path>.md`, *<heading>* ``.

The patterns are the gate's own, imported from the tool, so this test and the
gate cannot drift apart. String and character literals, raw strings included,
are skipped, so a `//` inside a literal is not read as a comment.
"""

from __future__ import annotations

import importlib.util
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

import pytest

# A harness test: it reads tools/ and the C++ sources, and imports nothing from tin_engine.
pytestmark = pytest.mark.harness

ROOT = Path(__file__).resolve().parents[2]
SOURCES = ("include", "bindings", "src")
SUFFIXES = {".hpp", ".cpp", ".h"}


def _load_gate() -> ModuleType:
    path = ROOT / "tools" / "check_code_history.py"
    spec = importlib.util.spec_from_file_location("cpp_comment_history_gate", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


GATE = _load_gate()
RULES = GATE.RULES
SCIENTIFIC = GATE.SCIENTIFIC
HEADING_CITATION = GATE.HEADING_CITATION
MD = GATE.MD


@dataclass(frozen=True)
class Comment:
    line: int
    text: str
    opening: bool  # part of the file's opening comment


def comments(source: str) -> Iterator[Comment]:
    """Every comment line of a C++ source, with string and character literals skipped."""
    i, line, n = 0, 1, len(source)
    code_seen = False  # a line of code other than `#pragma once` has been passed
    while i < n:
        c = source[i]
        if c == "\n":
            line += 1
            i += 1
        elif source.startswith("//", i):
            end = source.find("\n", i)
            end = n if end < 0 else end
            yield Comment(line, source[i + 2 : end], not code_seen)
            i = end
        elif source.startswith("/*", i):
            end = source.find("*/", i + 2)
            end = n if end < 0 else end
            for offset, text in enumerate(source[i + 2 : end].split("\n")):
                yield Comment(line + offset, text, not code_seen)
            line += source.count("\n", i, end)
            i = end + 2
        elif c.isspace():
            i += 1
        elif source.startswith("#pragma once", i):
            i += len("#pragma once")
        else:
            code_seen = True
            end = _skip_code(source, i)
            line += source.count("\n", i, end)
            i = end


def _skip_code(source: str, i: int) -> int:
    """Past one token of code, a whole literal if one starts at `i`."""
    c = source[i]
    if c == "R" and source.startswith('R"', i) and not _ident_char(source, i - 1):
        open_paren = source.index("(", i)
        delimiter = source[i + 2 : open_paren]
        return source.index(")" + delimiter + '"', open_paren) + len(delimiter) + 2
    if c == '"' or (c == "'" and not _ident_char(source, i - 1)):
        j = i + 1
        while source[j] != c:
            j += 2 if source[j] == "\\" else 1
        return j + 1
    return i + 1


def _ident_char(source: str, i: int) -> bool:
    """True where `source[i]` continues an identifier or number (so `'` is a digit separator)."""
    return i >= 0 and (source[i].isalnum() or source[i] == "_")


def history(source: str) -> Iterator[tuple[int, str, str]]:
    """Every history match in `source`'s comments, as (line, rule, matched text)."""
    allowed = 1  # the one heading-form pointer the opening comment may hold
    for comment in comments(source):
        for name, rule in RULES.items():
            for match in rule.finditer(comment.text):
                if not SCIENTIFIC.fullmatch(match.group(0)):
                    yield comment.line, name, match.group(0)
        heading = {m.start() + 1 for m in HEADING_CITATION.finditer(comment.text)}
        for match in MD.finditer(comment.text):
            if comment.opening and match.start() in heading and allowed:
                allowed -= 1
            else:
                yield comment.line, "pointer", match.group(0)


def cpp_sources() -> list[Path]:
    return sorted(
        path
        for top in SOURCES
        for path in (ROOT / top).rglob("*")
        if path.suffix in SUFFIXES and path.is_file()
    )


# --- the reader itself, on inline sources ---------------------------------


def test_a_slash_pair_inside_a_string_or_raw_string_is_not_a_comment() -> None:
    source = 'auto a = "see // increment 14";\nauto b = R"x(// R5 )x";\nchar c = \'"\';  // plain\n'
    assert [(c.line, c.text) for c in comments(source)] == [(3, " plain")]


def test_a_digit_separator_is_not_a_character_literal() -> None:
    source = "int n = 1'000;  // increment 14\n"
    assert list(history(source)) == [(1, "increment", "increment 14")]


def test_a_block_comment_reports_each_of_its_lines() -> None:
    source = "int x;\n/* first\n   R5 here */ int y;\n"
    assert list(history(source)) == [(3, "ruling-id", "R5")]


def test_one_heading_pointer_is_allowed_in_the_opening_comment_only() -> None:
    pointer = "// `docs/increments/14-adaptive-refinement.md`, *Rounds*\n"
    assert list(history("#pragma once\n\n" + pointer + "\nint x;\n")) == []
    assert list(history("#pragma once\n" + pointer + pointer)) == [
        (3, "pointer", "docs/increments/14-adaptive-refinement.md")
    ]
    assert list(history("int x;\n" + pointer)) == [
        (2, "pointer", "docs/increments/14-adaptive-refinement.md")
    ]


def test_a_bare_pointer_is_flagged_even_in_the_opening_comment() -> None:
    source = "#pragma once\n// (docs/increments/15f-edge-strip.md)\n"
    assert list(history(source)) == [(2, "pointer", "docs/increments/15f-edge-strip.md")]


# --- the pin -------------------------------------------------------------


def test_the_cpp_core_comments_carry_no_history() -> None:
    files = cpp_sources()
    assert files, f"no C++ sources found under {SOURCES}"
    found = [
        f"{path.relative_to(ROOT)}:{line}: {rule}: {text}"
        for path in files
        for line, rule, text in history(path.read_text(encoding="utf-8"))
    ]
    lines = {entry.split(": ")[0] for entry in found}
    assert found == [], f"{len(found)} matches on {len(lines)} lines:\n" + "\n".join(found)
