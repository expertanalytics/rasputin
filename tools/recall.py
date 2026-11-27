#!/usr/bin/env python3
"""What is already known on a topic: ledger lines, saved memories, benchmark verdicts.

    python3 tools/recall.py <topic words>... [--limit N] [--memory <folder>]

An entry (ledger line, memory file, benchmark record) matches when its text
holds every topic word, case-insensitive; each source prints at most
`--limit` entries, newest first, then "+N more".
1. The ledger: the file in `FINDINGS_LEDGER`, else `docs/findings.jsonl` on
   `private/master`, read with `git show`; numbered as `findings.py` numbers it.
2. Memories: `<folder>/*.md` but the index `MEMORY.md`, with their frontmatter
   name, description and date (else the file's), and each body line quoting
   Ola with a date.
3. Benchmark records, `docs/benchmarks/<date>/<name>/README.md` on
   `private/master`, matched on their whole text. The verdict line printed is
   the first line that starts bold or whose label is Verdict or Result; a bare
   heading such as "## Results" gives the next line of text instead.
Each run appends its time and words, as `tool: recall`, to the record
`findings.py search` keeps, so it counts as a ledger search too;
`questions.py add` and `findings.py add` refuse through `unrecalled` until a
run since question mode began (`questions.py mode enter`) named a word of the
question or ledger line. Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent))  # appended: the stdlib wins (h16 G4)
import findings  # beside this file

ROOT = Path(__file__).resolve().parents[1]
REF = "private/master"
#: Claude Code keys a project's memory folder on its path, each other character a "-".
MEMORY = Path.home() / ".claude/projects"
RECORD = re.compile(r"docs/benchmarks/([^/]+)/[^/]+/README\.md")
OLA = re.compile(r"\bOla,? (\d{4}-\d{2}-\d{2})[^\"]*?\"([^\"]+)\"")
VERDICT = re.compile(r"(?:Verdict|Results?)\b")
#: Chosen: words that name no topic, so a recall on them opens no question.
SMALL = {"the", "and", "for", "with", "from", "that", "this", "what", "which", "shall",
         "should", "does", "are", "was", "not", "but", "all", "any", "can", "will", "into",
         "our", "its", "how", "why", "when", "who", "has", "have", "now"}  # fmt: skip
#: Chosen: a few dated Ola lines per memory; the rest are counted.
QUOTES = 3


def unrecalled(searches: Path, text: str) -> str | None:
    """The refusal unless a recall since question mode began named a word of `text`."""
    import questions  # here, not at the top: it imports this module

    since = questions.mode_since(searches.parent)
    if since is None:
        return "enter question mode first: python3 tools/questions.py mode enter"
    try:
        lines = searches.read_text().splitlines()
    except OSError:
        lines = []
    named: set[str] = set()
    for line in lines:
        try:
            run = json.loads(line)
            if run.get("tool") == "recall" and datetime.fromisoformat(run["at"]) >= since:
                named.update(re.findall(r"\w+", " ".join(run["words"]).casefold()))
        except (ValueError, KeyError, TypeError, AttributeError):
            continue
    named = {w for w in named if len(w) > 2 and w not in SMALL}
    if any(t.startswith(w) for t in re.findall(r"\w+", text.casefold()) for w in named):
        return None
    return ("recall what is already known on this topic first: "
            "python3 tools/recall.py <a word of the question or ledger line>")  # fmt: skip


def _git(repo: Path, *argv: str, data: bytes | None = None) -> bytes:
    return subprocess.run(["git", "-C", str(repo), *argv], input=data, capture_output=True,
                          check=True).stdout  # fmt: skip


def _ledger(words: list[str], repo: Path, ref: str, ledger: Path | None) -> list[str]:
    try:
        shown = "" if ledger else _git(repo, "show", f"{ref}:docs/findings.jsonl").decode()
        entries = findings._load(ledger) if ledger else \
            [json.loads(line) for line in shown.splitlines() if line.strip()]  # fmt: skip
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        return [f"unreadable: {error}"]
    return [findings._line(n, e) for n, e in reversed(list(enumerate(entries, 1)))
            if all(w in " ".join([e["text"], *e["topics"]]).casefold() for w in words)]  # fmt: skip


def _memories(words: list[str], folder: Path) -> list[str]:
    found: list[tuple[str, str]] = []  # (date, entry)
    for path in sorted(set(folder.glob("*.md")) - {folder / "MEMORY.md"}):  # not the index
        text = path.read_text(errors="replace")
        if not all(w in text.casefold() for w in words):
            continue
        head, _, body = text.removeprefix("---\n").partition("\n---\n")
        field = dict(re.findall(r"^(name|description): *(.*)$", head, re.M))
        about = field.get("description", "")
        about = json.loads(about) if about.startswith('"') else about
        date = re.search(r"^\s*modified: *(\d{4}-\d{2}-\d{2})", head, re.M)
        stamp = datetime.fromtimestamp(path.stat().st_mtime).date().isoformat()
        day = date.group(1) if date else stamp
        quotes = [f'Ola, {m[1]}: "{m[2]}"' for m in map(OLA.search, body.splitlines()) if m]
        more = [f"+{len(quotes) - QUOTES} more Ola lines in the file"] * (len(quotes) > QUOTES)
        lines = "\n      ".join([f"{path}, {day}: {field.get('name', path.stem)}: {about}",
                              *quotes[:QUOTES], *more])  # fmt: skip
        found.append((day, lines))
    return [lines for _, lines in sorted(found, key=lambda f: f[0], reverse=True)]


def verdict(text: str) -> str:
    """The record's verdict line, by the rule in the module's docstring."""
    lines = [line.strip() for line in text.splitlines()]
    for n, line in enumerate(lines):
        label = line.lstrip("#* ")
        if line.startswith("**") or (VERDICT.match(label) and line[:1] in "#*"):
            if line.startswith("#") and VERDICT.fullmatch(label.rstrip(":. ")):
                return next((r for r in lines[n + 1 :] if r and r[0] != "#"), label)
            return label.replace("**", "")
    return "(no bold or verdict line)"


def _benchmarks(words: list[str], repo: Path, ref: str) -> list[str]:
    try:
        listed = _git(repo, "ls-tree", "-r", "--name-only", ref, "--", "docs/benchmarks")
        paths = [p for p in listed.decode().splitlines() if RECORD.fullmatch(p)]
        out = _git(repo, "cat-file", "--batch", data="".join(f"{ref}:{p}\n" for p in paths)
                   .encode())  # fmt: skip
    except subprocess.CalledProcessError as error:
        return [f"unreadable: {error.stderr.decode().strip()}"]
    found = []
    for path in paths:
        header, _, out = out.partition(b"\n")
        size = int(header.split()[2])
        text, out = out[:size].decode(errors="replace"), out[size + 1 :]
        if all(w in f"{path}\n{text}".casefold() for w in words):
            found.append(f"{path}, {RECORD.fullmatch(path)[1]}: {verdict(text)}")  # type: ignore[index]
    return sorted(found, reverse=True)


def _section(head: str, found: list[str], limit: int, numbered: bool = True) -> list[str]:
    if found[:1] and found[0].startswith("unreadable"):
        return [f"{head}:", f"   {found[0]}"]
    shown = [f"   {f'{n}. ' * numbered}{line}" for n, line in enumerate(found[:limit], 1)]
    more = [f"   +{len(found) - limit} more: rerun with --limit {len(found)}"] * \
        (len(found) > limit)  # fmt: skip
    return [f"{head}, {len(shown)} of {len(found)}:", *(shown or ["   no match"]), *more]


def command(words: list[str], *, repo: Path, ref: str, ledger: Path | None, memory: Path,
            searches: Path, now: datetime, limit: int) -> list[str]:  # fmt: skip
    """Recall `words` from the three sources; record the run; return the lines to print."""
    folded = [w.casefold() for w in words]
    findings._append(searches, {"at": now.isoformat(), "words": words, "tool": "recall"})
    where = str(ledger) if ledger else f"{ref}:docs/findings.jsonl"
    return "\n".join(
        [  # a memory's entry spans lines
            f"Recall: {' '.join(words)} (every word, any case)",
            *_section(f"Ledger ({where})", _ledger(folded, repo, ref, ledger), limit, False),
            *_section(f"Memories ({memory})", _memories(folded, memory), limit),
            *_section(
                f"Benchmark verdicts ({ref}:docs/benchmarks)", _benchmarks(folded, repo, ref), limit
            ),
        ]
    ).splitlines()


def main(argv: list[str] | None = None) -> int:
    import questions

    p = argparse.ArgumentParser(description="What is already known on a topic.")
    p.add_argument("words", nargs="+")
    p.add_argument("--limit", type=int, default=10, help="lines per source; the rest counted")
    p.add_argument("--memory", type=Path, help="the memory folder; default the main session's")
    args = p.parse_args(argv)
    variable = os.environ.get("FINDINGS_LEDGER")
    searches = questions.store_of(ROOT).with_name(findings.SEARCHES)
    main_repo = searches.parent.parent  # the common git directory's parent, for every worktree
    memory = args.memory or MEMORY / re.sub(r"[^A-Za-z0-9]", "-", str(main_repo)) / "memory"
    print("\n".join(command(args.words, repo=ROOT, ref=REF, memory=memory,
                            ledger=Path(variable) if variable else None, searches=searches,
                            now=datetime.now(UTC), limit=args.limit)))  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
