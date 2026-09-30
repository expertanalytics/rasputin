#!/usr/bin/env python3
"""Check the citations the prose and the code comments make, and flag those at risk.

Governance and increment docs cite each other and the source by line number,
as `docs/increments/<file>.md:<n>` or `<path>.hpp:<a>-<b>`. Line numbers drift,
and worse, the *text* at a line can be rewritten to say the opposite while the
number still resolves. On 2026-09-17 a sweep resolved all 46 citations in
`kernel-sufficiency-audit.md` and was still wrong three times, because three of
them pointed at text the same branch had rewritten. `@reviewer` put it exactly:
resolving a line number is not the same as resolving the quotation.

So this script does not claim to find wrong citations. The design is
`docs/increments/h2-citation-pinning.md`. It knows five forms:

  unpinned  -- `<path>:<n>`, resolved in the working tree (`legacy/` against
               the `legacy-archive` tag);
  pinned    -- `<path>@<rev>:<n>`, resolved with `git cat-file` at a commit
               sha or a tag. Branches and HEAD move, so they are not pins;
  ID        -- `<path>.md` followed by a section sign and an ID, resolved
               against the ID labels the target's headings declare;
  heading   -- a backticked `<path>.md`, a comma and an italic heading name.

and reports, in this order of precedence:

  broken    -- the file, the rev, the line or the heading is missing.
  unpinned  -- an unpinned line citation into a governed rule file. Rule files
               are rewritten in place, so such a citation is pinned or cited by
               ID or heading instead.
  at-risk   -- an unpinned citation into a file this branch modifies: the
               number may still resolve while the quotation no longer holds.
               Re-read these against the new text; no script can.

Exit status is 1 for `broken` or `unpinned`. `at-risk` is a worklist.

Only prose (`.md`), comments and docstrings are scanned, never string literals.

Usage: python3 tools/check_citations.py [--base master] [--paths docs .claude]
"""

from __future__ import annotations

import argparse
import ast
import io
import re
import subprocess
import sys
import tokenize
from functools import cache
from pathlib import Path, PurePosixPath

REPO = Path(__file__).resolve().parent.parent

# `legacy/` left the working tree with release hygiene, but the increment
# records cite it by line 147 times as their evidence. An unpinned citation
# whose path starts with this prefix is read as pinned to the archive tag,
# never against the working tree -- a stray or regrown `legacy/` in the
# checkout must not answer for the archive. See release-hygiene.md section 3.
LEGACY_PREFIX = "legacy/"
LEGACY_TAG = "legacy-archive"

# A path is a token with an extension, which keeps ratios and timestamps out.
# One starting with `/` is an absolute path or the tail of a URL, never a
# citation.
PATH = r"(?<![\w/.-])(?!/)([\w./-]+\.(?:md|py|hpp|cpp|h|yaml|yml|toml|txt|cmake))"
LINES = r":(\d+)(?:-(\d+))?"
PINNED = re.compile(PATH + r"@([\w./-]+)" + LINES)
CITATION = re.compile(PATH + LINES)
MD_PATH = r"(?<![\w/.-])(?!/)([\w./-]+\.md)"
ID_CITATION = re.compile(MD_PATH + r"`?[ \t]*§([\w.]*\w)")
ID = re.compile(r"\d+(?:\.\d+)*[A-Z]?|[A-Z]{1,2}\d+|[A-Z]")
HEADING_CITATION = re.compile("`" + MD_PATH + r"`, \*([^*\n]+)\*")
SHA = re.compile(r"[0-9a-f]{7,40}")

HEADING = re.compile(r"(#{1,6})\s+(.*?)\s*$")
FENCE = re.compile(r" {0,3}(`{3,}|~{3,})")
NUMERIC_LABEL = re.compile(r"(\d+(?:\.\d+)*[A-Z]?)(?:\.|\s|$)")
ALNUM_LABEL = re.compile(r"([A-Z]{1,2}\d+)(?:[.:]|\s|$)")
LETTER_LABEL = re.compile(r"([A-Z])\.(?:\s|$)")
NAME_ENDS = (":", " [", " —", " (")

# Where a partial path may live. More than one hit is reported as ambiguous
# rather than picked between.
SEARCH_ROOTS = ("docs", ".claude", "include", "src", "src_python", "tests", "tools", "lib")

SCAN_SUFFIXES = (".md", ".py", ".h", ".hpp", ".cpp", ".cmake")
C_SUFFIXES = (".h", ".hpp", ".cpp")
SKIPPED_PREFIXES = ("lib/", "legacy/", ".claude/worktrees/")

# The rule files. An unpinned line citation into one of these is a failure.
# Step 4 of the generic harness replaces this and guard_governance.py's copy.
GOVERNED = (
    "CLAUDE.md",
    "docs/PRINCIPLES.md",
    ".claude/REQUIRED-READING.md",
    "docs/increments/README.md",
    "testing.md",
)
GOVERNED_PREFIXES = (".claude/agents/", ".claude/skills/", ".claude/hooks/", "tools/check_")

# Sessions check out git worktrees here, each a full copy of the repo. Walking
# into them made every bare basename ambiguous and every doc scanned twice
# (229 "broken" on 2026-09-29, the first commit with the gates hook active).
WORKTREES = REPO / ".claude" / "worktrees"


def in_repo_proper(path: Path) -> bool:
    return not path.is_relative_to(WORKTREES)


def git(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)


def changed_files(base: str) -> set[str] | None:
    """Repo-relative paths this branch modifies, or None if base is unknown."""
    for ref in (f"{base}...HEAD", base):
        result = git("diff", "--name-only", ref)
        if result.returncode == 0:
            return {line for line in result.stdout.split("\n") if line}
    return None


@cache
def repo_files() -> list[str]:
    """Tracked and untracked, unignored files that exist, repo-relative."""
    listed = git("ls-files", "-z", "--cached", "--others", "--exclude-standard").stdout
    return sorted({rel for rel in listed.split("\0") if rel and (REPO / rel).is_file()})


def scannable(rel: str) -> bool:
    name = PurePosixPath(rel).name
    wanted = PurePosixPath(rel).suffix in SCAN_SUFFIXES or name == "CMakeLists.txt"
    return wanted and not rel.startswith(SKIPPED_PREFIXES)


def resolve(cited: str) -> Path | list[Path] | None:
    """A cited path, as a file, an ambiguity, or nothing.

    A path that is not a file at the root is matched by its trailing path
    components. An ambiguous path is returned as the candidate list rather
    than silently resolved to the first hit: guessing which README.md was
    meant would make this script commit the error it exists to catch.
    """
    direct = REPO / cited
    if direct.is_file():
        return direct
    parts = PurePosixPath(cited).parts
    hits = [
        REPO / rel
        for rel in repo_files()
        if rel.startswith(tuple(f"{root}/" for root in SEARCH_ROOTS))
        and PurePosixPath(rel).parts[-len(parts) :] == parts
        and in_repo_proper(REPO / rel)
    ]
    if len(hits) == 1:
        return hits[0]
    return hits or None


def governed(rel: str) -> bool:
    return rel in GOVERNED or rel.startswith(GOVERNED_PREFIXES)


@cache
def lines_of(path: Path) -> list[str]:
    return path.read_text(errors="replace").splitlines()


@cache
def refs() -> list[str]:
    return [line.split(" ", 1)[-1] for line in git("show-ref").stdout.splitlines()]


@cache
def pin_error(rev: str) -> str | None:
    """Why `rev` is not a pin, or None when it names a present commit or tag."""
    spec = rev if SHA.fullmatch(rev) else f"refs/tags/{rev}"
    if git("rev-parse", "--verify", "--quiet", f"{spec}^{{commit}}").returncode == 0:
        return None
    if rev == "HEAD" or any(
        ref == f"refs/heads/{rev}" or (ref.startswith("refs/remotes/") and ref.endswith(f"/{rev}"))
        for ref in refs()
    ):
        return f"'{rev}' is a branch or HEAD, which moves: not a pin; pin a sha or a tag"
    if git("rev-parse", "--is-shallow-repository").stdout.strip() == "true":
        return (
            f"no commit or tag '{rev}' in this clone, which is shallow; "
            f"it needs full history (fetch-depth: 0)"
        )
    return f"no commit or tag '{rev}'"


@cache
def lines_at(rev: str, path: str) -> int | None:
    result = git("cat-file", "-p", f"{rev}:{path}")
    return len(result.stdout.splitlines()) if result.returncode == 0 else None


def check_pinned(path: str, rev: str, last: int, legacy: bool = False) -> str | None:
    """The broken message for a pinned line citation, or None when it resolves."""
    # The `legacy/` shorthand keeps its own wording, which names the tag.
    archive = f"not in the working tree; resolved against tag '{rev}', where the"
    if (error := pin_error(rev)) is not None:
        return f"{archive} tag is missing" if legacy else error
    total = lines_at(rev, path)
    if total is None:
        return f"{archive} file is missing" if legacy else f"no '{path}' at '{rev}'"
    if last > total or last == 0:
        return f"file has {total} lines" if legacy else f"'{path}' has {total} lines at '{rev}'"
    return None


def headings(path: Path) -> list[tuple[int, str]]:
    """The ATX headings of `path` outside fenced blocks, as (level, text)."""
    found: list[tuple[int, str]] = []
    fenced = False
    for line in lines_of(path):
        if FENCE.match(line):
            fenced = not fenced
        elif not fenced and (match := HEADING.match(line)):
            found.append((len(match.group(1)), match.group(2)))
    return found


@cache
def declared_ids(path: Path) -> set[str]:
    """The IDs `path`'s headings declare (spec section 2.3)."""
    ids: set[str] = set()
    stack: list[tuple[int, str | None]] = []  # enclosing headings, with numeric labels
    for level, text in headings(path):
        text = text.removeprefix("§")
        while stack and stack[-1][0] >= level:
            stack.pop()
        numeric = NUMERIC_LABEL.match(text)
        other = ALNUM_LABEL.match(text) or LETTER_LABEL.match(text)
        if numeric:
            ids.add(numeric.group(1))
        elif other:
            ids.add(other.group(1))
            if LETTER_LABEL.match(text) and stack and stack[-1][1]:
                ids.add(f"{stack[-1][1]}{other.group(1)}")
        stack.append((level, numeric.group(1) if numeric else None))
    return ids


def heading_named(path: Path, name: str) -> bool:
    return any(
        text == name or any(text.startswith(name + end) for end in NAME_ENDS)
        for _, text in headings(path)
    )


def python_text(text: str) -> dict[int, list[str]]:
    """Comments and docstrings of a Python source, by line. Raises on bad syntax."""
    out: dict[int, list[str]] = {}
    lines = text.splitlines()
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type == tokenize.COMMENT:
            out.setdefault(token.start[0], []).append(token.string)
    for node in ast.walk(ast.parse(text)):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            first = node.body[0] if node.body else None
            if (
                isinstance(first, ast.Expr)
                and isinstance(first.value, ast.Constant)
                and isinstance(first.value.value, str)
            ):
                for number in range(first.lineno, (first.end_lineno or first.lineno) + 1):
                    out.setdefault(number, []).append(lines[number - 1])
    return out


RAW_PREFIX = re.compile(r"(?:^|[^\w])(?:u8|[uUL])?R$")
RAW_OPEN = re.compile(r'"([^()\\\s]{0,16})\(')


def c_text(text: str, cmake: bool = False) -> dict[int, list[str]]:
    """Comment text of a C-family (or, with `cmake`, a CMake) source, by line."""
    out: dict[int, list[str]] = {}

    def emit(start: int, end: int, line: int) -> None:
        for k, part in enumerate(text[start:end].split("\n")):
            out.setdefault(line + k, []).append(part)

    i, n, line = 0, len(text), 1
    while i < n:
        c, j = text[i], i + 1
        eol = text.find("\n", i)
        eol = n if eol < 0 else eol
        if (cmake and c == "#") or (not cmake and text.startswith("//", i)):
            emit(i, eol, line)
            j = eol
        elif not cmake and text.startswith("/*", i):
            close = text.find("*/", i + 2)
            j = n if close < 0 else close + 2
            emit(i, j, line)
        elif (
            c == '"'
            and not cmake
            and (raw := RAW_OPEN.match(text, i))
            and RAW_PREFIX.search(text[max(0, i - 3) : i])
        ):
            close = text.find(f'){raw.group(1)}"', raw.end())
            j = n if close < 0 else close + len(raw.group(1)) + 2
        elif c == '"' or (c == "'" and not cmake):
            between = 0 < i < n - 1 and text[i - 1].isalnum() and text[i + 1].isalnum()
            if not (c == "'" and between):  # a digit separator, as in 0x9513'0000
                stop = n if cmake else eol  # CMake quotes span lines; C's end at the line
                while j < stop and text[j] != c:
                    j += 2 if text[j] == "\\" else 1
                j = min(j + 1, n) if j < stop else stop
        line += text.count("\n", i, j)
        i = j
    return out


def scanned_text(path: Path) -> dict[int, list[str]]:
    """The text a citation is sought in, by line (spec section 4)."""
    text = path.read_text(errors="replace")
    if path.suffix == ".py":
        return python_text(text)
    if path.suffix in C_SUFFIXES:
        return c_text(text)
    if path.suffix == ".cmake" or path.name == "CMakeLists.txt":
        return c_text(text, cmake=True)
    return {number: [line] for number, line in enumerate(text.splitlines(), start=1)}


def default_sources(paths: list[str] | None) -> list[Path]:
    """The files to scan: named files whatever their type, walked ones filtered."""
    if paths is None:
        return [REPO / rel for rel in repo_files() if scannable(rel)]
    sources: list[Path] = []
    for entry in paths:
        target = REPO / entry
        if target.is_file():
            sources.append(target)
        elif target.is_dir():
            sources.extend(
                p
                for p in sorted(target.rglob("*"))
                if p.is_file()
                and in_repo_proper(p)
                and (not p.is_relative_to(REPO) or scannable(p.relative_to(REPO).as_posix()))
            )
    return sources


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", default="master")
    parser.add_argument("--paths", nargs="*", default=None)
    args = parser.parse_args()

    maybe_touched = changed_files(args.base)
    diffable = maybe_touched is not None
    if maybe_touched is None:
        print(f"warning: cannot diff against '{args.base}'; at-risk check disabled")
    touched = maybe_touched or set()

    broken: list[str] = []
    unpinned: list[str] = []
    at_risk: list[str] = []
    warnings: list[str] = []

    def resolved(where: str, cited: str) -> Path | None:
        target = resolve(cited)
        if target is None:
            broken.append(f"{where}: cites '{cited}' -- no such file")
        elif isinstance(target, list):
            names = ", ".join(str(p.relative_to(REPO)) for p in target)
            broken.append(
                f"{where}: cites '{cited}' -- ambiguous ({names}); cite the path, not the basename"
            )
        else:
            return target
        return None

    for source in default_sources(args.paths):
        rel_source = source.relative_to(REPO) if source.is_relative_to(REPO) else source
        try:
            text_by_line = scanned_text(source)
        except (SyntaxError, tokenize.TokenError, ValueError):
            warnings.append(f"{rel_source}: not scanned, does not parse")
            continue
        for number, fragments in sorted(text_by_line.items()):
            where = f"{rel_source}:{number}"
            for fragment in fragments:
                rest = fragment
                for m in PINNED.finditer(fragment):
                    path, rev, start, end = m.groups()
                    last = max(int(start), int(end or 0))
                    if (error := check_pinned(path, rev, last)) is not None:
                        broken.append(f"{where}: cites '{path}@{rev}:{start}' -- {error}")
                    rest = rest[: m.start()] + " " * len(m.group(0)) + rest[m.end() :]
                for m in ID_CITATION.finditer(rest):
                    cited, ident = m.groups()
                    if not ID.fullmatch(ident):
                        continue
                    target = resolved(where, cited)
                    if target is not None and ident not in declared_ids(target):
                        broken.append(
                            f"{where}: cites '{cited}' §{ident} -- no heading declares §{ident}"
                        )
                for m in HEADING_CITATION.finditer(rest):
                    cited, name = m.groups()
                    target = resolved(where, cited)
                    if target is not None and not heading_named(target, name):
                        broken.append(f"{where}: cites '{cited}', *{name}* -- no such heading")
                for m in CITATION.finditer(rest):
                    cited, start, end = m.groups()
                    last = max(int(start), int(end or 0))
                    if cited.startswith(LEGACY_PREFIX):
                        if (
                            error := check_pinned(cited, LEGACY_TAG, last, legacy=True)
                        ) is not None:
                            broken.append(f"{where}: cites '{cited}:{start}' -- {error}")
                        continue
                    target = resolved(where, cited)
                    if target is None:
                        continue
                    total = len(lines_of(target))
                    if last > total or last == 0:
                        broken.append(f"{where}: cites '{cited}:{start}' -- file has {total} lines")
                        continue
                    rel_target = target.relative_to(REPO).as_posix()
                    if governed(rel_target):
                        unpinned.append(
                            f"{where}: cites '{cited}:{start}' in rule file {rel_target} -- "
                            f"pin it (path@<sha>:<n>) or cite a heading"
                        )
                    elif rel_target in touched:
                        at_risk.append(
                            f"{where}: cites '{cited}:{start}', and this branch edits "
                            f"{rel_target} -- re-read the quotation, not the number"
                        )

    for title, items in (
        ("at risk", at_risk),
        ("unpinned", unpinned),
        ("broken", broken),
        ("warnings", warnings),
    ):
        if items:
            print(f"\n== {title} ({len(items)}) ==")
            for item in items:
                print(f"  {item}")
    if broken or unpinned:
        return 1
    if not at_risk:
        print(
            "All citations resolve"
            + (
                ", and none point into a file this branch edits."
                if diffable
                else ". At-risk was not evaluated -- see the warning above."
            )
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
