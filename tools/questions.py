#!/usr/bin/env python3
"""The questions for Ola: checked before they are stored, numbered once, answered verbatim.

The store is untracked, `<git-common-dir>/questions.jsonl`, one JSON record per
line, appended only; the newest line for an id is its record, so every worktree
sees one list. Beside it, `questions-shown.json` holds which number the last
`list` gave to which id, so Ola's "2" keeps meaning what he saw. Standard
library only: the recap imports this on the system `python3`.
Spec: `docs/increments/h35-question-tool.md`, *1. What it is*
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Iterable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
#: Chosen: two lines of the recap's 120 characters.
LIMIT = 240
NO_DEFAULT = {"", "none", "tbd", "?", "-"}
ELSEWHERE = re.compile(r"see above|as above|as before|previous question|question \d+", re.I)
COUNTED = re.compile(r"\b(?:increment|row) \d+", re.I)
TREE = ("tools/", "docs/", "src_python/", "tests/", "include/", ".claude/", ".github/")
# fmt: off
KNOWN = {"CI", "PR", "DEM", "CRS", "GIS", "MCP", "CPU", "GPU", "UTC", "CLI", "API", "URL", "CSV",
         "JSON", "EPSG", "QGIS", "NVE", "UTM", "AC", "CORINE", "OK", "PDF", "TIN", "VTK", "HTTP"}
# fmt: on
PUNCT = "?.,;:!()\"'"
HASH = re.compile(r"[0-9a-f]{7,40}")
ABBREVIATION = re.compile(r"[A-Z]{2,6}")


class RefusedError(Exception):
    """A question the checks refuse, with the one plain sentence saying why."""


def load(store: Path) -> tuple[list[dict[str, Any]], int]:
    """The newest record per id, in first-seen order, and the count of unreadable lines."""
    if not store.exists():
        return [], 0
    latest: dict[str, dict[str, Any]] = {}
    bad = 0
    for line in store.read_text().splitlines():
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            record = None
        if not isinstance(record, dict) or not isinstance(record.get("id"), str):
            bad += 1
            continue
        latest[record["id"]] = record
    return list(latest.values()), bad


def _key(text: str) -> str:
    return " ".join(text.split()).casefold()


def _append(store: Path, *records: dict[str, Any]) -> None:
    store.parent.mkdir(parents=True, exist_ok=True)
    with store.open("a") as out:
        out.writelines(json.dumps(r) + "\n" for r in records)


def _shown_path(store: Path) -> Path:
    return store.with_name("questions-shown.json")


def _yes_no(record: dict[str, Any]) -> bool:
    return record["kind"] in ("question", "act") and not record["choices"]


def _clash(new: dict[str, Any], others: Iterable[dict[str, Any]]) -> dict[str, Any] | None:
    """The first open question the same short answer would fit, if any."""
    words = {c.casefold() for c in new["choices"]}
    for other in others:
        if other["state"] != "open" or other["kind"] == "taken" or other["id"] == new["id"]:
            continue
        if _yes_no(new) and _yes_no(other):
            return other
        if words & {c.casefold() for c in other["choices"]}:
            return other
    return None


def _plain(text: str, terms: set[str], act: bool) -> None:
    """Refuse a label, number, hash, handle, branch, repository path or unknown abbreviation."""
    if found := COUNTED.search(text):
        raise RefusedError(f"say what is meant, not the label '{found.group()}'")
    spelled = set(re.findall(r"\(([A-Z]{2,6})\)", text))
    personas = {p.stem for p in (ROOT / ".claude" / "agents").glob("*.md")}
    for raw in text.split():
        tok = raw.strip(PUNCT)
        if not tok or tok.startswith(("/", "~", "$")) or "://" in tok:
            continue
        kind = None
        if tok.startswith(TREE) or ("/" in tok and re.search(r"\.\w+$", tok)):
            kind = "a path inside the repository"
        elif re.fullmatch(r"h\d+|[A-Z]\d{1,2}", tok):
            kind = "a design label"
        elif re.fullmatch(r"#\d+", tok):
            kind = "a pull request number"
        elif HASH.fullmatch(tok) and not tok.isdigit() and not tok.isalpha():
            kind = "a commit hash"
        elif tok.startswith("@") and tok[1:] not in personas:
            kind = "a handle Ola does not know"
        elif tok.startswith("worktree-") and not act:
            kind = "a branch name"
        elif ABBREVIATION.fullmatch(tok) and tok not in KNOWN | terms | spelled:
            kind = "an abbreviation; spell it out, or name a proper name with --term"
        if kind:
            raise RefusedError(f"say it in plain words: '{tok}' is {kind}")


def _check(args: argparse.Namespace, records: list[dict[str, Any]], kind: str) -> None:
    """Refuse a text the checks reject; return quietly otherwise."""
    text, default = args.text, args.default
    if default.strip().casefold() in NO_DEFAULT:
        raise RefusedError("a question needs a default: what happens if Ola is silent")
    if kind != "taken":
        if text.count("?") > 1:
            raise RefusedError("one decision per question: the text has more than one '?'")
        if not text.rstrip().endswith("?"):
            raise RefusedError("a question must end with a question mark")
        key = _key(text)
        for r in records:
            if r["state"] in ("open", "held", "answered") and _key(r["text"]) == key:
                said = f': Ola, {r["at"][:10]}: "{r["words"]}"' if r["state"] == "answered" else ""
                raise RefusedError(f"already asked as {r['id']}{said}")
    if ELSEWHERE.search(text):
        raise RefusedError("each question must read alone: it points elsewhere")
    _plain(text, set(args.term or ()), kind == "act")
    if len(text) + len(default) > LIMIT:
        raise RefusedError(f"at most {LIMIT} characters; say the rest in the increment file")


def _line(n: int, r: dict[str, Any]) -> str:
    if r["kind"] == "taken":
        return f"{n}. Taken: {r['text'].rstrip('.')}. Veto by number."
    choices = f" Choices: {', '.join(r['choices'])}." if r["choices"] else ""
    return f"{n}. {r['text']}{choices} Default: {r['default'].rstrip('.')}."


def _shown(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The records `list` numbers: open questions, open acts, then defaults taken not yet shown."""
    open_ = [r for r in records if r["state"] == "open"]
    return (
        [r for r in open_ if r["kind"] == "question"]
        + [r for r in open_ if r["kind"] == "act"]
        + [r for r in open_ if r["kind"] == "taken" and not r.get("shown")]
    )


def list_lines(records: list[dict[str, Any]]) -> list[str]:
    """The block for Ola, numbered 1..n in one sequence; reads nothing and writes nothing."""
    lines = [_line(n, r) for n, r in enumerate(_shown(records), 1)]
    if held := sum(r["state"] == "held" for r in records):
        noun = "question" if held == 1 else "questions"
        lines.append(f"and {held} more yes/no {noun} held until one above is answered")
    return lines


def _numbered(store: Path, records: list[dict[str, Any]], n: int) -> dict[str, Any]:
    """The record the last `list` printed as number `n`, still open."""
    path = _shown_path(store)
    ids: dict[str, str] = json.loads(path.read_text()) if path.exists() else {}
    record = next((r for r in records if r["id"] == ids.get(str(n))), None)
    if record is None or record["state"] != "open":
        raise RefusedError(f"the last list printed no open question numbered {n}")
    return record


def _close(
    store: Path, records: list[dict[str, Any]], record: dict[str, Any], now: datetime, **fields: Any
) -> list[str]:
    """Close `record` and release the oldest held question nothing open clashes with."""
    closed = {**record, **fields, "at": now.isoformat()}
    rest = [closed if r["id"] == record["id"] else r for r in records]
    released = next((r for r in rest if r["state"] == "held" and _clash(r, rest) is None), None)
    out = [closed]
    lines = []
    if released:
        out.append({**released, "state": "open", "at": now.isoformat()})
        lines.append(f"Now open: {released['text']} Default: {released['default']}")
    _append(store, *out)
    return lines


def command(args: argparse.Namespace, store: Path, now: datetime) -> list[str]:
    """Run one command on `store` at time `now`; return the lines to print."""
    records, _ = load(store)
    if args.cmd in ("add", "taken"):
        kind = "taken" if args.cmd == "taken" else "act" if args.act else "question"
        _check(args, records, kind)
        number = max((int(r["id"][1:]) for r in records if r["id"][1:].isdigit()), default=0)
        choices = [c for c in (getattr(args, "choices", "") or "").split("|") if c]
        record = {
            "id": f"q{number + 1}",
            "kind": kind,
            "text": args.text,
            "default": args.default,
            "choices": choices,
            "item": getattr(args, "item", None),
            "source": args.source,
            "state": "open",
            "words": None,
            "at": now.isoformat(),
        }
        clash = _clash(record, records) if kind != "taken" else None
        if clash:
            record["state"] = "held"
        _append(store, record)
        if clash:
            why = "the same short answer would fit both"
            return [f"held until question {clash['id']} is answered: {why}"]
        return [f"stored as {record['id']}"]
    if args.cmd == "list":
        lines = list_lines(records)
        shown = _shown(records)
        _shown_path(store).parent.mkdir(parents=True, exist_ok=True)
        _shown_path(store).write_text(json.dumps({str(n): r["id"] for n, r in enumerate(shown, 1)}))
        taken = [{**r, "shown": True} for r in shown if r["kind"] == "taken"]
        if taken:
            _append(store, *taken)
        return lines or ["No open questions."]
    record = _numbered(store, records, args.n)
    if args.cmd == "answer":
        released = _close(store, records, record, now, state="answered", words=args.words)
        return [f'Ola, {now.date().isoformat()}: "{args.words}"', *released]
    released = _close(store, records, record, now, state="dropped", reason=args.reason)
    return [f"Dropped: {record['text']}", *released]


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Questions for Ola, checked, numbered and answered.")
    sub = p.add_subparsers(dest="cmd", required=True)
    term_help = "a proper name such as Bane NOR that reads as an abbreviation; not for jargon"
    for name in ("add", "taken"):
        s = sub.add_parser(
            name, help="ask a question" if name == "add" else "record a default taken"
        )
        s.add_argument("text")
        s.add_argument("--default", default="", help="what happens if Ola is silent")
        s.add_argument("--source", default="main", help="who asks: a persona or main")
        s.add_argument("--term", action="append", help=term_help)
        if name == "add":
            s.add_argument("--choices", default="", help="the answers, separated by |")
            s.add_argument("--act", action="store_true", help="a publishing act awaiting a yes")
            s.add_argument("--item", help="the work item the question is about")
    sub.add_parser("list", help="print the numbered block for Ola")
    a = sub.add_parser("answer", help="close question N with Ola's words, verbatim")
    a.add_argument("n", type=int)
    a.add_argument("words")
    d = sub.add_parser("drop", help="close question N without an answer")
    d.add_argument("n", type=int)
    d.add_argument("reason")
    return p


def main(argv: list[str] | None = None) -> int:
    common = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "--path-format=absolute", "--git-common-dir"],
        capture_output=True,
        text=True,
        check=True,
    )
    try:
        lines = command(
            parser().parse_args(argv),
            Path(common.stdout.strip()) / "questions.jsonl",
            datetime.now(UTC),
        )
    except RefusedError as error:
        print(error, file=sys.stderr)
        return 1
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
