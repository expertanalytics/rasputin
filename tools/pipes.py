#!/usr/bin/env python3
"""Named work lists ("pipes"), the item a free writer takes next, and stuck items.

A pipe's definition is tracked, `.claude/pipes/<name>.toml` (`name`, `purpose`),
read from the working tree and from `<remote>/master` of every remote but
`origin`, never fetching: the recap runs in a hook, which must not touch the
network. A pipe's live list is untracked, `<git-common-dir>/pipes/<name>.jsonl`,
one JSON record per line, appended only; the newest line for an id is its
record. A pipe's own record (id = the pipe's name, state `pipe`) holds its
priority and any open question about the pipe. Standard library only: the
recap imports this on the system `python3`.
Spec: `docs/increments/h31-pipes.md`, *7. Blueprint: modules, data flow, and what the recap prints*
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import subprocess
import sys
import tomllib
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

sys.path.append(str(Path(__file__).resolve().parent))  # appended: the stdlib wins
import harness_mode

ROOT = Path(__file__).resolve().parents[1]
#: Chosen: about three of Ola's daily check cycles, which allows a weekend.
STUCK_AFTER = timedelta(days=3)
CAP, PIPE_CAP = 3, 2  # items in flight across all pipes, and per pipe
RANK = {"review": 0, "green": 1, "red": 2, "designing": 3}  # nearest to merging first
STATES = ("queued", *RANK, "waiting", "merged", "parked", "dropped")
IDLE = "IDLE: no runnable item in any pipe"
NOTE = re.compile(r"^[^-]+-(.+)-\d{6}\.md$")
Run = Callable[[list[str]], tuple[int, str]]


@dataclass(frozen=True)
class Item:
    id: str
    state: str
    pipe: str
    since: datetime
    title: str = ""
    next: str = ""
    branch: str | None = None
    pr: int | None = None
    question: str | None = None
    default: str | None = None
    plain: bool = False
    fallback: bool = False
    rank: float = 0  # lower is earlier in the pipe; `top` lowers it
    before: str | None = None  # the state `ask` left, which `answer` restores
    priority: int | None = None  # on a pipe's own record only
    priority_words: str | None = None


@dataclass(frozen=True)
class Pipe:
    name: str
    purpose: str
    priority: int | None
    source: str


@dataclass(frozen=True)
class Choice:
    item: Item | None
    reason: str


FIELDS = {f.name for f in dataclasses.fields(Item)}


def fold(records: list[dict[str, Any]]) -> dict[str, Item]:
    """The newest record per id, in the order ids first appear."""
    items: dict[str, Item] = {}
    for r in records:
        known = {k: v for k, v in r.items() if k in FIELDS}
        items[r["id"]] = Item(**{**known, "since": datetime.fromisoformat(r["since"])})
    return items


def _records(path: Path) -> tuple[list[dict[str, Any]], int]:
    """The readable records of a live list, and how many lines are not."""
    good, bad = [], 0
    text = path.read_text() if path.is_file() else ""
    for line in filter(str.strip, text.splitlines()):
        try:
            record = json.loads(line)
            fold([record])
            good.append(record)
        except (ValueError, TypeError, KeyError, AttributeError):
            bad += 1
    return good, bad


def load_items(path: Path) -> dict[str, Item]:
    return fold(_records(path)[0])


def _git(args: list[str]) -> tuple[int, str]:
    done = subprocess.run(["git", *args], capture_output=True, text=True, check=False)
    return done.returncode, done.stdout


def _parse(text: str, source: str, where: str, warnings: list[str]) -> Pipe | None:
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        warnings.append(f"pipes: {where}: not TOML ({error})")
        return None
    extra = sorted(set(data) - {"name", "purpose"})
    warnings += [f"pipes: {where}: unknown key {key}" for key in extra]
    name = str(data.get("name", Path(where).stem))
    return Pipe(name=name, purpose=str(data.get("purpose", "")), priority=None, source=source)


def _definitions(repo: Path, run: Run) -> tuple[list[Pipe], list[str], list[str]]:
    """Pipes, warnings, and one note per remote read for the header."""
    found: list[Pipe | None] = []
    warnings: list[str] = []
    notes: list[str] = []
    for path in sorted((repo / ".claude" / "pipes").glob("*.toml")):
        found.append(_parse(path.read_text(), "tree", str(path.relative_to(repo)), warnings))
    git = ["-C", str(repo)]
    for remote in run([*git, "remote"])[1].split():
        if remote == "origin":
            continue
        ref = f"{remote}/master"
        code, sha = run([*git, "rev-parse", "--verify", "-q", "--short", ref])
        if code != 0:
            warnings.append(f"pipes: no {ref} locally; git fetch {remote} to read its pipes")
            continue
        date = run([*git, "log", "-1", "--format=%cs", ref])[1].strip()
        source = f"{ref} {sha.strip()}"
        notes.append(f"pipes read from {source}, committed {date}")
        names = run([*git, "ls-tree", "--name-only", ref, ".claude/pipes/"])[1].split()
        if not names:
            notes[-1] += ": no pipe definitions"
        for name in names:
            if name.endswith(".toml") and Path(name).stem not in {p.name for p in found if p}:
                text = run([*git, "show", f"{ref}:{name}"])[1]
                found.append(_parse(text, source, name, warnings))
    return [p for p in found if p], warnings, notes


def definitions(repo: Path, run: Run = _git) -> tuple[list[Pipe], list[str]]:
    """The tree's pipe definitions, then those only other remotes' masters hold."""
    found, warnings, _ = _definitions(repo, run)
    return found, warnings


def _priorities(pipes: list[Pipe]) -> dict[str, int]:
    return {p.name: 99 if p.priority is None else p.priority for p in pipes}


def in_flight(items: dict[str, Item]) -> int:
    return sum(i.state in RANK for i in items.values())


def stuck(items: dict[str, Item], now: datetime) -> list[Item]:
    """`queued` and `waiting` items unchanged for longer than STUCK_AFTER, oldest first."""
    waits = [i for i in items.values() if i.state in ("queued", "waiting")]
    return sorted((i for i in waits if now - i.since > STUCK_AFTER), key=lambda i: i.since)


def order(pipes: list[Pipe], items: dict[str, Item], now: datetime, unattended: bool,
          busy: set[str]) -> Choice:  # fmt: skip
    """The item the next free writer takes, and why; `Choice(None, IDLE)` when none."""
    pri, place = _priorities(pipes), {key: n for n, key in enumerate(items)}

    def runnable(i: Item) -> bool:
        return i.state != "waiting" and (i.plain or not unattended)

    flying = [i for i in items.values() if i.state in RANK]
    free = [i for i in flying if runnable(i) and i.branch not in busy]
    if free:
        first = min(free, key=lambda i: (RANK[i.state], pri.get(i.pipe, 99), i.since))
        return Choice(first, "nearest to merging")
    if len(flying) >= CAP:
        return Choice(None, IDLE)
    queued = [i for i in items.values() if i.state == "queued" and runnable(i)]
    plain_queue = [i for i in queued if not i.fallback]
    open_pipes = {i.pipe for i in plain_queue}
    held = {p: sum(i.pipe == p for i in flying) for p in open_pipes}
    ready = [i for i in plain_queue if held[i.pipe] < PIPE_CAP or open_pipes == {i.pipe}]
    old = [i for i in ready if now - i.since > STUCK_AFTER]
    if old:
        first = min(old, key=lambda i: i.since)
        return Choice(first, f"stuck {(now - first.since).days} d")
    if ready:
        first = min(ready, key=lambda i: (pri.get(i.pipe, 99), i.rank, place[i.id]))
        top = pri.get(first.pipe, 99) == min(pri.values(), default=99)
        return Choice(first, f"{first.pipe} is {'first' if top else 'next by priority'}")
    backup = [i for i in queued if i.fallback and i.branch not in busy]
    if backup:
        return Choice(min(backup, key=lambda i: (pri.get(i.pipe, 99), place[i.id])), "fallback")
    return Choice(None, IDLE)


def _notes(tasks: Path) -> dict[str, str]:
    """Branch `worktree-<name>` to the note file `<persona>-<name>-<HHMMSS>.md` naming it."""
    names = [p.name for p in tasks.iterdir()] if tasks.is_dir() else []
    return {f"worktree-{m[1]}": n for n in names if (m := NOTE.match(n))}


def busy_branches(tasks: Path) -> set[str]:
    return set(_notes(tasks))


@dataclass
class State:
    """Everything the commands and the table read, gathered once."""

    pipes: list[Pipe]
    items: dict[str, Item]  # every record, pipes' own included
    warnings: list[str]
    notes: list[str]
    unreadable: int
    lists: Path
    tasks: Path


def gather(repo: Path) -> State:
    common = Path(_git(["-C", str(repo), "rev-parse", "--path-format=absolute",
                        "--git-common-dir"])[1].strip())  # fmt: skip
    pipes, warnings, notes = _definitions(repo, _git)
    lists, bad, records = common / "pipes", 0, list[dict[str, Any]]()
    for path in sorted(lists.glob("*.jsonl")):
        good, n = _records(path)
        records, bad = records + good, bad + n
    items = fold(records)
    pipes = [dataclasses.replace(p, priority=items[p.name].priority) if p.name in items else p
             for p in pipes]  # fmt: skip
    main = common.parent if common.name == ".git" else repo
    return State(pipes, items, warnings, notes, bad, lists, main / ".claude" / "current-task")


def _open(items: dict[str, Item]) -> dict[str, Item]:
    return {k: i for k, i in items.items() if i.state != "pipe"}


def _choose(state: State, now: datetime) -> tuple[Choice, dict[str, str], bool]:
    mode = harness_mode.read_mode(state.lists.parent / "harness", now)
    unattended = mode.state != "off"
    notes = _notes(state.tasks)
    return order(state.pipes, _open(state.items), now, unattended, set(notes)), notes, unattended


def _next_line(choice: Choice) -> str:
    i = choice.item
    return f"next: {i.pipe}/{i.id} {i.state} ({choice.reason})" if i else IDLE


def _cell(text: str, width: int) -> str:
    return (text if len(text) < width else text[: width - 4] + "...").ljust(width)


WIDTHS = (5, 8, 4, 30, 10, 12, 11)
HEADS = ("", "PIPE", "PRI", "ITEM", "STEP", "WAITING ON", "SINCE")


def show_lines(repo: Path, now: datetime) -> list[str]:
    """The recap's table of every open item, its questions, and the next choice."""
    state = gather(repo)
    choice, notes, unattended = _choose(state, now)
    items, pri = _open(state.items), _priorities(state.pipes)
    late = {i.id: n for n, i in enumerate(stuck(items, now))}
    place = {key: n for n, key in enumerate(items)}

    def key(i: Item) -> tuple[Any, ...]:
        group = RANK.get(i.state, 4 if i.state == "waiting" else 5)
        return (late.get(i.id, len(late)), pri.get(i.pipe, 99), group, i.rank,
                i.since if i.state in RANK else now, place[i.id])  # fmt: skip

    rows = sorted((i for i in items.values() if i.state in (*RANK, "waiting", "queued")), key=key)
    asks = [i for i in rows if i.state == "waiting"]
    asks += [i for i in state.items.values() if i.state == "pipe" and i.question]
    unreadable = [f"{state.unreadable} unreadable"] * bool(state.unreadable)
    head = "; ".join([f"in flight {in_flight(items)} of {CAP}", *state.notes, *unreadable])
    lines = [f"Pipes ({head})"[:120], *(w[:120] for w in state.warnings)]
    lines.append("".join(_cell(h, w) for h, w in zip(HEADS, WIDTHS, strict=True)) + "NEXT")
    for i in rows:
        age = now - i.since
        since = f"{age.days} d" if age >= timedelta(days=1) else f"{age // timedelta(hours=1)} h"
        who = ("Ola" if i.state == "waiting" else "Ola present" if unattended and not i.plain
               else "its agent" if i.branch in notes else "free agent")  # fmt: skip
        nxt = (f"question {asks.index(i) + 1} below" if i in asks
               else " ".join(filter(None, [i.next, f"#{i.pr}" if i.pr else ""])))  # fmt: skip
        since += " stuck" if i.id in late else ""
        cells = ("next" if i == choice.item else "", i.pipe, str(pri.get(i.pipe, "-")),
                 f"{i.id} {i.title}", i.state, who, since)  # fmt: skip
        lines.append(("".join(map(_cell, cells, WIDTHS)) + nxt)[:120].rstrip())
    if asks:
        lines.append("Questions for Ola:")
    for n, i in enumerate(asks, 1):
        where = i.pipe if i.state == "pipe" else f"{i.pipe}/{i.id}"
        lines.append(f"  {n}. {where}: {str(i.question)[:60]}  default: {i.default}"[:120])
    lines.append(_next_line(choice))
    return lines[:40] + ([f"... {len(lines) - 40} more"] if len(lines) > 40 else [])


class RefusedError(Exception):
    pass


def _append(st: State, item: Item, now: datetime, **change: Any) -> None:
    """Append `item` with `change` applied; `since` moves only when the state does."""
    new = dataclasses.replace(item, **change)
    if new.state != item.state:
        new = dataclasses.replace(new, since=now.replace(microsecond=0))
    record = {**dataclasses.asdict(new), "since": new.since.isoformat()}
    st.lists.mkdir(parents=True, exist_ok=True)
    with (st.lists / f"{new.pipe}.jsonl").open("a") as out:
        out.write(json.dumps(record) + "\n")


def _item(state: State, ident: str) -> Item:
    found = state.items.get(ident)
    if found is None or found.state == "pipe":
        raise RefusedError(f"pipes: no item {ident}")
    return found


def _pipe_record(state: State, name: str, now: datetime) -> Item:
    if name not in {p.name for p in state.pipes}:
        raise RefusedError(f"pipes: no pipe {name} defined")
    return state.items.get(name) or Item(id=name, state="pipe", pipe=name, since=now)


def command(args: argparse.Namespace, repo: Path, now: datetime) -> list[str]:
    """Apply one command; the lines it prints."""
    state = gather(repo)
    if args.cmd == "add":
        if args.fallback and args.governed:
            raise RefusedError("pipes: --fallback with --governed: a fallback is plain")
        if args.id in state.items:
            raise RefusedError(f"pipes: id {args.id} is taken in pipe {state.items[args.id].pipe}")
        _pipe_record(state, args.pipe, now)
        _append(state, Item(id=args.id, state="queued", pipe=args.pipe, since=now), now,
                title=args.title, next=args.next or "", plain=args.plain or args.fallback,
                fallback=args.fallback)  # fmt: skip
    elif args.cmd == "set":
        item, change = _item(state, args.id), {"state": args.state}
        change |= {k: v for k in ("branch", "pr", "next") if (v := getattr(args, k)) is not None}
        if args.plain or args.governed:
            change["plain"] = args.plain
        _append(state, item, now, **change)
    elif args.cmd == "ask":
        if args.id in {p.name for p in state.pipes}:
            _append(state, _pipe_record(state, args.id, now), now, question=args.question,
                    default=args.default)  # fmt: skip
        else:
            item = _item(state, args.id)
            _append(state, item, now, state="waiting", before=item.before or item.state,
                    question=args.question, default=args.default)  # fmt: skip
    elif args.cmd == "answer":
        item = state.items.get(args.id) or _item(state, args.id)
        back = item.state if item.state == "pipe" else item.before or "queued"
        _append(state, item, now, state=back, before=None, question=None, default=None)
    elif args.cmd == "priority":
        for pair in args.pairs:
            name, _, value = pair.partition("=")
            _append(state, _pipe_record(state, name, now), now, priority=int(value),
                    priority_words=args.words)  # fmt: skip
    elif args.cmd == "top":
        item = _item(state, args.id)
        ranks = [i.rank for i in _open(state.items).values() if i.pipe == item.pipe]
        _append(state, item, now, rank=min(ranks) - 1)
    elif args.cmd == "next":
        choice, notes, _ = _choose(state, now)
        flying = [i for i in _open(state.items).values() if i.state in RANK and i.branch in notes]
        return [f"busy: {i.pipe}/{i.id} ({notes[str(i.branch)]})" for i in flying] + [
            _next_line(choice)
        ]
    else:
        return show_lines(repo, now)
    return []


def parser() -> argparse.ArgumentParser:
    top = argparse.ArgumentParser(description="Pipes: named work lists and the next item.")
    sub = top.add_subparsers(dest="cmd")
    add = sub.add_parser("add", help="queue a new item in a pipe")
    add.add_argument("pipe")
    add.add_argument("id")
    add.add_argument("title")
    add.add_argument("--fallback", action="store_true", help="the pipe's standing fallback")
    put = sub.add_parser("set", help="change an item's state")
    put.add_argument("id")
    put.add_argument("state", choices=STATES)
    put.add_argument("--branch")
    put.add_argument("--pr", type=int)
    for p in (add, put):
        p.add_argument("--next", help="the next step, in plain words")
        p.add_argument("--plain", action="store_true", help="writes no governed path")
        p.add_argument("--governed", action="store_true", help="writes a governed path")
    ask = sub.add_parser("ask", help="put a question for Ola on an item or a pipe")
    ask.add_argument("id")
    ask.add_argument("question")
    ask.add_argument("--default", required=True, help="what is taken if Ola is silent")
    answer = sub.add_parser("answer", help="record Ola's ruling and close the question")
    answer.add_argument("id")
    answer.add_argument("words")
    pri = sub.add_parser("priority", help="set pipe priorities, 1 highest")
    pri.add_argument("pairs", nargs="+", metavar="PIPE=N")
    pri.add_argument("--words", required=True, help="Ola's words")
    sub.add_parser("top", help="move an item to the head of its pipe").add_argument("id")
    sub.add_parser("next", help="the item the next free writer takes")
    sub.add_parser("show", help="the table the recap prints")
    return top


def main(argv: list[str] | None = None) -> int:
    try:
        lines = command(parser().parse_args(argv), ROOT, datetime.now(UTC))
    except RefusedError as error:
        print(error, file=sys.stderr)
        return 1
    if lines:
        print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
