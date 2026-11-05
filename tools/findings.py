#!/usr/bin/env python3
"""The ledger of settled results, so a known result is found before it is called new.

    python3 tools/findings.py add "<one-line result>" --commit <sha> --topic <word>... --ola "<..>"
    python3 tools/findings.py search <words>...
    python3 tools/findings.py list

    python3 tools/findings.py replace <n> "<new line>" --ola "<Ola's words>" | drop <n>

The ledger is data, `docs/findings.jsonl` or the path in `FINDINGS_LEDGER`, one
JSON entry per fact; `--commit` names a commit of the repository holding the
ledger, which may be private. A missing file is an empty ledger, and git
history keeps what `replace` and `drop` remove.
`add` and `replace` need Ola's agreement to the wording: `--ola`, words found in
this session's transcript as `brief.py --ola` finds them, or `--question <id>`, a
question `questions.py` stored whose text holds the line and whose answer is
one of `AGREE`, alone or with the box's `RECOMMENDED` suffix. Each `search`
appends its time and words to `<git-common-dir>/findings-searches.jsonl`,
untracked; `questions.py add` and design or probe briefs (`brief.py`) refuse
through `unsearched` until a search ran in the last hour. Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

sys.path.append(str(Path(__file__).resolve().parent))  # appended: the stdlib wins (h16 G4)

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "docs" / "findings.jsonl"
SEARCHES = "findings-searches.jsonl"
#: Chosen: one line of a result, as questions.py caps a question at two.
LIMIT = 160
WINDOW = timedelta(hours=1)
#: Chosen: whole answers that agree; anything else, "Yes, but ..." too, is refused.
AGREE = {"yes", "agree", "agreed", "ok", "approve", "approved"}
#: The question box's suffix on a click of the recommended option, as in "Yes (Recommended)".
RECOMMENDED = " (Recommended)"


class RefusedError(Exception):
    """An entry the checks refuse, with the one plain sentence saying why."""


def unsearched(searches: Path, now: datetime) -> str | None:
    """The refusal when no search of the settled results ran in the hour before `now`."""
    try:
        last = datetime.fromisoformat(json.loads(searches.read_text().splitlines()[-1])["at"])
    except (OSError, IndexError, ValueError, KeyError, TypeError):
        last = None
    if last is not None and now - last <= WINDOW:
        return None
    return ("search the settled results first, in case this is already known: "
            "python3 tools/findings.py search <words>")  # fmt: skip


def _load(ledger: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in ledger.read_text().splitlines() if line.strip()] \
        if ledger.exists() else []  # fmt: skip


def _append(path: Path, record: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as out:
        out.write(json.dumps(record) + "\n")


def _line(n: int, e: dict[str, Any]) -> str:
    return f"{n}. {e['text']} (commit {e['commit'][:8]}, {e['date']}; {', '.join(e['topics'])})"


def _answered(args: argparse.Namespace, store: Path) -> str:
    """Ola's agreeing answer to question `args.question`, whose text holds the line."""
    import questions

    qid = f"q{args.question}" if args.question.isdigit() else args.question
    record = next((r for r in questions.load(store)[0] if r["id"] == qid), None)
    if record is None:
        raise RefusedError(f"there is no question {qid}")
    if record["state"] != "answered":
        raise RefusedError(f"question {qid} is not answered")
    if record["words"].removesuffix(RECOMMENDED).strip(" .!").casefold() not in AGREE:
        raise RefusedError(f'Ola did not agree to question {qid}: "{record["words"]}"')
    if " ".join(args.text.split()) not in " ".join(record["text"].split()):
        raise RefusedError(f"question {qid} does not hold this line")
    return f'"{record["words"]}" (question {qid}, {record["at"]})'


def _worded(args: argparse.Namespace, searches: Path) -> str:
    """Ola's agreement with the wording, as brief.py quotes it, once the wording is plain."""
    import brief  # here, not at the top: these import this module
    import questions

    try:
        questions._plain(args.text, set(args.term or ()), False)
        quoted = _answered(args, searches.with_name("questions.jsonl")) if args.question \
            else brief._ola([args.ola])[1].strip()  # fmt: skip
    except (questions.RefusedError, brief.BriefError) as error:
        raise RefusedError(str(error)) from error
    if len(args.text) > LIMIT:
        raise RefusedError(f"at most {LIMIT} characters: one line of result")
    return quoted


def command(args: argparse.Namespace, ledger: Path, searches: Path, now: datetime) -> list[str]:
    """Run one command on `ledger` at time `now`; return the lines to print."""
    entries = _load(ledger)
    if args.cmd in ("replace", "drop") and not 0 < args.n <= len(entries):
        raise RefusedError(f"the ledger has no entry {args.n}")
    if args.cmd == "add":
        ola = _worded(args, searches)
        home = next(d for d in ledger.parents if d.is_dir())  # the ledger's own repository
        shown = subprocess.run(["git", "-C", str(home), "show", "-s", "--format=%H %cs",
                                f"{args.commit}^{{commit}}"], capture_output=True, text=True,
                               check=False).stdout.split()  # fmt: skip
        if len(shown) != 2:
            raise RefusedError(f"git has no commit {args.commit}")
        entries.append({"text": args.text, "commit": shown[0], "date": shown[1],
                        "topics": args.topic, "ola": ola})  # fmt: skip
    elif args.cmd == "replace":
        entries[args.n - 1].update(text=args.text, ola=_worded(args, searches))
    elif args.cmd == "drop":
        del entries[args.n - 1]
    if args.cmd in ("add", "replace", "drop"):
        ledger.parent.mkdir(parents=True, exist_ok=True)
        ledger.write_text("".join(json.dumps(e) + "\n" for e in entries))
        return [f"{args.cmd}: done; the ledger holds {len(entries)}"]
    if args.cmd == "list":
        return [_line(n, e) for n, e in enumerate(entries, 1)] or ["The ledger is empty."]
    words = [w.casefold() for w in args.words]
    _append(searches, {"at": now.isoformat(), "words": args.words})

    def hits(e: dict[str, Any]) -> int:
        return sum(w in " ".join([e["text"], *e["topics"]]).casefold() for w in words)

    found = sorted((e for e in entries if hits(e)), key=hits, reverse=True)
    return [_line(n, e) for n, e in enumerate(found, 1)] or [
        f"No settled result matches: {' '.join(args.words)}"
    ]


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="The ledger of settled results, one entry per fact.")
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("add", help="record one settled result")
    r = sub.add_parser("replace", help="rewrite entry N in place")
    r.add_argument("n", type=int)
    for s in (a, r):
        s.add_argument("text")
        agreed = s.add_mutually_exclusive_group(required=True)
        agreed.add_argument("--ola", help="Ola's words agreeing to the wording")
        agreed.add_argument("--question", help="an answered question holding the line, as q7")
        s.add_argument("--term", action="append", help="a proper name read as an abbreviation")
    a.add_argument("--commit", required=True, help="the commit that settled it")
    a.add_argument("--topic", action="append", required=True, help="a word to find it by")
    sub.add_parser("drop", help="remove entry N").add_argument("n", type=int)
    sub.add_parser("search", help="print matching results; lets questions through for an hour") \
        .add_argument("words", nargs="+")  # fmt: skip
    sub.add_parser("list", help="print every result")
    return p


def main(argv: list[str] | None = None) -> int:
    import questions

    try:
        searches = questions.store_of(ROOT).with_name(SEARCHES)
        ledger = Path(os.environ.get("FINDINGS_LEDGER", LEDGER))
        lines = command(parser().parse_args(argv), ledger, searches, datetime.now(UTC))
    except RefusedError as error:
        print(error, file=sys.stderr)
        return 1
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
