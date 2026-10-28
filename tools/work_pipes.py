#!/usr/bin/env python3
"""Named work lists ("pipes"), the item a free writer takes next, and stuck items.

A pipe's definition is tracked, `.claude/pipes/<name>.toml` (`name`, `purpose`),
read from the working tree and from `<remote>/master` of every remote but
`origin`, never fetching: the recap runs in a hook, which must not touch the
network. A pipe's live list is untracked, `<git-common-dir>/pipes/<name>.jsonl`,
one JSON record per line, appended only; the newest line for an id is its
record. A pipe's own record (id = the pipe's name, state `pipe`) holds its
priority. Questions for Ola live in `tools/questions.py` alone: an item waits
on Ola while a question there names it with `--item`. Standard library only:
the recap imports this on the system `python3`.
Spec: `docs/increments/h31-pipes.md`, *7. Blueprint: modules, data flow, and what the recap prints*
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import re
import subprocess
import sys
import textwrap
import tomllib
from collections.abc import Callable
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any

sys.path.append(str(Path(__file__).resolve().parent))  # appended: the stdlib wins
import harness_mode
import questions

ROOT = Path(__file__).resolve().parents[1]
#: Chosen: about three of Ola's daily check cycles, which allows a weekend.
STUCK_AFTER = timedelta(days=3)
CAP, PIPE_CAP = 3, 2  # items in flight across all pipes, and per pipe
RANK = {"review": 0, "green": 1, "red": 2, "designing": 3}  # nearest to merging first
STATES = ("queued", *RANK, "waiting", "merged", "parked", "dropped")
IDLE = "IDLE: no runnable item in any pipe"
NOTE = re.compile(r"^[^-]+-(.+)-\d{6}\.md$")
NUMBER = re.compile(r"\d+\.\d+")  # `<pipe>.<position>`, as `numbers` gives them
#: Chosen: the pipes Ola names most, in his order; any other follows, by name.
PIPE_NUMBERS = {"mesh": 1, "lsm": 2, "harness": 3}
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
    plain: bool = False
    fallback: bool = False
    rank: float = 0  # lower is earlier in the pipe; `top` lowers it
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
    asked: list[dict[str, Any]]  # the question tool's records


def _asked(state: State) -> dict[str, int]:
    """Item id to the number the question tool gives its first open question."""
    found: dict[str, int] = {}
    for n, record in enumerate(questions.open_records(state.asked), 1):
        if record.get("item") and record["kind"] != "taken":
            found.setdefault(record["item"], n)
    return found


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
    asked = questions.load(common / "questions.jsonl")[0]
    return State(pipes, items, warnings, notes, bad, lists, main / ".claude" / "current-task",
                 asked)  # fmt: skip


def _open(items: dict[str, Item]) -> dict[str, Item]:
    return {k: i for k, i in items.items() if i.state != "pipe"}


def _choose(state: State, now: datetime) -> tuple[Choice, dict[str, str], bool]:
    mode = harness_mode.read_mode(state.lists.parent / "harness", now)
    unattended = mode.state != "off"
    notes, asked = _notes(state.tasks), _asked(state)
    items = {k: dataclasses.replace(i, state="waiting") if k in asked else i
             for k, i in _open(state.items).items()}  # fmt: skip
    return order(state.pipes, items, now, unattended, set(notes)), notes, unattended


def _owed(repo: Path) -> list[str]:
    """The duties line `tools/duties.py` gives the prompt hook, as a list; none on a fault."""
    import duties

    try:
        found = duties.line(repo, "hook")
    except Exception:  # the duties are advice here: a fault never stops `next`
        return []
    return [found] if found else []


def _next_line(choice: Choice, nums: dict[str, str]) -> str:
    i = choice.item
    return f"next: {nums[i.id]} {i.pipe}/{i.id} {i.state} ({choice.reason})" if i else IDLE


def _cell(text: str, width: int) -> str:
    return (text if len(text) < width else text[: width - 4] + "...").ljust(width)


WIDTHS = (6, 5, 8, 4, 30, 10, 12, 11)
HEADS = ("#", "", "PIPE", "PRI", "ITEM", "STEP", "WAITING ON", "SINCE")


def _show_order(state: State, now: datetime) -> list[Item]:
    """The open items in flight, waiting or queued: stuck first, then by pipe and step."""
    items, pri = _open(state.items), _priorities(state.pipes)
    late = {i.id: n for n, i in enumerate(stuck(items, now))}
    place = {key: n for n, key in enumerate(items)}

    def key(i: Item) -> tuple[Any, ...]:
        group = RANK.get(i.state, 4 if i.state == "waiting" else 5)
        return (late.get(i.id, len(late)), pri.get(i.pipe, 99), group, i.rank,
                i.since if i.state in RANK else now, place[i.id])  # fmt: skip

    return sorted((i for i in items.values() if i.state in (*RANK, "waiting", "queued")), key=key)


def numbers(state: State, now: datetime) -> dict[str, str]:
    """Item id to `<pipe>.<position>`: the position counts its pipe's rows in `show` order.

    A number names an item for as long as nothing moves or closes; it is a handle
    for typing, never stored.
    """
    rows = _show_order(state, now)
    others = sorted(({p.name for p in state.pipes} | {i.pipe for i in rows}) - set(PIPE_NUMBERS))
    pipe_no = PIPE_NUMBERS | {name: len(PIPE_NUMBERS) + n for n, name in enumerate(others, 1)}
    found, count = dict[str, str](), dict.fromkeys(pipe_no, 0)
    for i in rows:
        count[i.pipe] += 1
        found[i.id] = f"{pipe_no[i.pipe]}.{count[i.pipe]}"
    return found


def show_lines(repo: Path, now: datetime) -> list[str]:
    """The recap's table of every open item, its questions, and the next choice."""
    state = gather(repo)
    choice, notes, unattended = _choose(state, now)
    items, pri = _open(state.items), _priorities(state.pipes)
    late = {i.id for i in stuck(items, now)}
    rows, nums = _show_order(state, now), numbers(state, now)
    asked = _asked(state)
    unreadable = [f"{state.unreadable} unreadable"] * bool(state.unreadable)
    head = "; ".join([f"in flight {in_flight(items)} of {CAP}", *state.notes, *unreadable])
    lines = [f"Pipes ({head})"[:120], *(w[:120] for w in state.warnings)]
    lines.append("".join(_cell(h, w) for h, w in zip(HEADS, WIDTHS, strict=True)) + "NEXT")
    for i in rows:
        age = now - i.since
        since = f"{age.days} d" if age >= timedelta(days=1) else f"{age // timedelta(hours=1)} h"
        who = ("Ola" if i.state == "waiting" or i.id in asked
               else "Ola present" if unattended and not i.plain
               else "its agent" if i.branch in notes else "free agent")  # fmt: skip
        nxt = (f"question {asked[i.id]} below" if i.id in asked
               else " ".join(filter(None, [i.next, f"#{i.pr}" if i.pr else ""])))  # fmt: skip
        since += " stuck" if i.id in late else ""
        cells = (nums[i.id], "next" if i == choice.item else "", i.pipe, str(pri.get(i.pipe, "-")),
                 f"{i.id} {i.title}", i.state, who, since)  # fmt: skip
        lines.append(("".join(map(_cell, cells, WIDTHS)) + nxt)[:120].rstrip())
    if block := questions.list_lines(state.asked):
        lines.append("Questions for Ola:")
    for line in block:
        lines += textwrap.wrap(line, 120, subsequent_indent="   ")
    lines.append(_next_line(choice, nums))
    return lines[:40] + ([f"... {len(lines) - 40} more"] if len(lines) > 40 else [])


def _md(text: str) -> str:
    """`text` safe inside a Markdown table cell: one line, `|` escaped."""
    return " ".join(text.split()).replace("|", "\\|")


def _waits_on(state: State, item: Item, asked: dict[str, int], running: set[str]) -> str:
    """What `item` really waits for: Ola's answer, a running agent, or a free writer."""
    if item.id in asked:
        return f"Ola (question {asked[item.id]})"
    held = [r for r in state.asked if r["state"] == "held" and r.get("item") == item.id]
    if held or item.state == "waiting":
        return "Ola (a held question)" if held else "Ola"
    return "a running agent" if item.id in running else "a free writer"


def markdown_lines(repo: Path, now: datetime) -> list[str]:
    """The verbose view in GitHub-flavoured Markdown: every open item, nothing truncated."""
    state = gather(repo)
    choice, _, _ = _choose(state, now)
    items, pri = _open(state.items), _priorities(state.pipes)
    asked = _asked(state)
    running = {branch.removeprefix("worktree-") for branch in _notes(state.tasks)}
    nxt, nums = choice.item, numbers(state, now)
    show_place = {i.id: n for n, i in enumerate(_show_order(state, now))}
    taken = IDLE
    if nxt:
        taken = f"{nums[nxt.id]} **{nxt.id}** {_md(nxt.title)} ({nxt.pipe}, {choice.reason})"
    lines = [f"In flight {in_flight(items)} of {CAP}. A free writer takes next: {taken}."]
    lines += [f"> {w}" for w in state.warnings]
    group = {**RANK, "waiting": 4, "queued": 5}
    names = sorted({p.name for p in state.pipes} | {i.pipe for i in items.values()},
                   key=lambda name: (pri.get(name, 99), name))  # fmt: skip
    for name in names:
        rows = sorted((i for i in items.values() if i.pipe == name and i.state in group),
                      key=lambda i: show_place[i.id])  # fmt: skip
        lines += ["", f"### {name}", ""]
        if not rows:
            lines.append("No open items.")
            continue
        lines += ["| # | Item | Step | Waiting on | Next |", "|---|---|---|---|---|"]
        for i in rows:
            todo = " ".join(filter(None, [i.next, f"#{i.pr}" if i.pr else ""]))
            cells = (nums[i.id], f"**{_md(i.id)}** {_md(i.title)}", i.state,
                     _waits_on(state, i, asked, running), _md(todo))  # fmt: skip
            lines.append("| " + " | ".join(cells) + " |")
    block = questions.list_lines(state.asked)
    numbered = [line for line in block if line[:1].isdigit()]
    held = sum(r["state"] == "held" for r in state.asked)
    lines += ["", "### Questions for Ola", "", *(numbered or ["No open questions."])]
    if held:
        lines += ["", block[-1][0].upper() + block[-1][1:] + "."]
    return lines


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


def _item(state: State, ident: str, now: datetime) -> Item:
    """The open item with id `ident`, or with number `ident` in `show` order."""
    if ident not in state.items and NUMBER.fullmatch(ident):
        ident = {n: k for k, n in numbers(state, now).items()}.get(ident, ident)
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
        if NUMBER.fullmatch(args.id):
            raise RefusedError(f"pipes: id {args.id} has the shape of an item number")
        if args.id in state.items:
            raise RefusedError(f"pipes: id {args.id} is taken in pipe {state.items[args.id].pipe}")
        _pipe_record(state, args.pipe, now)
        _append(state, Item(id=args.id, state="queued", pipe=args.pipe, since=now), now,
                title=args.title, next=args.next or "", plain=args.plain or args.fallback,
                fallback=args.fallback)  # fmt: skip
    elif args.cmd == "set":
        item, change = _item(state, args.id, now), {"state": args.state}
        change |= {k: v for k in ("branch", "pr", "next") if (v := getattr(args, k)) is not None}
        if args.plain or args.governed:
            change["plain"] = args.plain
        _append(state, item, now, **change)
    elif args.cmd == "priority":
        for pair in args.pairs:
            name, _, value = pair.partition("=")
            _append(state, _pipe_record(state, name, now), now, priority=int(value),
                    priority_words=args.words)  # fmt: skip
    elif args.cmd == "top":
        item = _item(state, args.id, now)
        ranks = [i.rank for i in _open(state.items).values() if i.pipe == item.pipe]
        _append(state, item, now, rank=min(ranks) - 1)
    elif args.cmd == "next":
        choice, notes, _ = _choose(state, now)
        flying = [i for i in _open(state.items).values() if i.state in RANK and i.branch in notes]
        nums = numbers(state, now)
        return [f"busy: {nums[i.id]} {i.pipe}/{i.id} ({notes[str(i.branch)]})" for i in flying] + [
            _next_line(choice, nums),
            *([] if choice.item else _owed(repo)),
        ]
    elif getattr(args, "markdown", False):
        return markdown_lines(repo, now)
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
    pri = sub.add_parser("priority", help="set pipe priorities, 1 highest")
    pri.add_argument("pairs", nargs="+", metavar="PIPE=N")
    pri.add_argument("--words", required=True, help="Ola's words")
    sub.add_parser("top", help="move an item to the head of its pipe").add_argument("id")
    sub.add_parser("next", help="the item the next free writer takes")
    show = sub.add_parser("show", help="the table the recap prints")
    show.add_argument("--markdown", action="store_true", help="every open item, as Markdown")
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
