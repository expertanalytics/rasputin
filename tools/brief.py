#!/usr/bin/env python3
"""The fixed part of a persona's brief, as one hashed block (h9).

Spec: docs/increments/h9-spawn-briefs.md §3.1 to §3.4. The main session runs

    python3 tools/brief.py <persona> --worktree <path> --beside <none|persona:path>...
                           [--increment <file>] [--no-build] [--untimed] [--ola <text>]...
                           [--lane fast] [--item <id>] [--area <path>...]

and pastes the block, unchanged, before the task. `.claude/hooks/guard_spawn.py`
reads it back with `find_blocks` and `check`, so the format has one owner. A
refusal is exit 2 with the reason on stderr and no block.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import string
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent))  # appended: the stdlib wins (h16 G4)
import session_state  # beside this file, as the hooks import it

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / ".claude" / "briefs" / "common.md"

#: The write limit per persona, h6's table (docs/increments/h6-role-limits.md §3).
#: Read-only is derived from an empty entry (§3.1a). When h6 lands, its ROLES replace this.
WRITES: dict[str, tuple[str, ...]] = {
    "architect": ("docs/ except docs/retrospectives/", "ROADMAP.md", "CLAUDE.md",
                  ".claude/ files ending in .md", "root *.md files"),
    "developer": ("src_python/", "include/", "src/", "bindings/", "tools/", ".claude/hooks/",
                  ".github/", "CMakeLists.txt", "pyproject.toml"),
    "orchestrator": ("docs/retrospectives/", "the lessons log ../rasputin_scratch/lessons.md"),
    "perf": ("docs/benchmarks/",),
    "reviewer": (),
    "tester": ("tests/",),
}  # fmt: skip
PERSONAS = tuple(WRITES)
#: The fast lane's developer writes the test, and rule text is fast-lane work.
FAST_WRITES = ("tests/", "CLAUDE.md", ".claude/ files ending in .md", "docs/increments/README.md")
NEEDS_INCREMENT = ("tester", "developer", "reviewer")
#: The line .claude/hooks/record_review.py records (h29 §3.1); the reviewer's brief demands it.
REVIEW_LINE = (
    "End the handback with exactly one line `Review: <APPROVED | APPROVED WITH WORDINGS | "
    "CHANGES REQUESTED> <40-hex sha of the head you reviewed>`, and start no other line with "
    "`Review:`; without it no review is recorded (docs/increments/README.md step 4)."
)
QUOTED = re.compile(r"(?i)invariant-critical|mutation|@perf|acceptance run|tools/bench\.py")
HEADER = re.compile(
    r"^<<<BRIEF persona=([a-z]+) worktree=(/\S+) head=([0-9a-f]{40}) hash=([0-9a-f]{12})>>>[ \t]*$"
)
END = re.compile(r"^<<<END BRIEF ([0-9a-f]{12})>>>[ \t]*$")
MAX_QUOTED = 8
MAX_WRITERS = 3  # pipe 3.6; the rule text below says "three"


class BriefError(Exception):
    """A brief that cannot be made: exit 2, the reason on stderr."""


class MalformedError(ValueError):
    """Text with a marker that is not a well-formed block."""


@dataclass(frozen=True)
class Block:
    persona: str
    worktree: str
    head: str
    hash: str
    body: str


def _normalise(body: str) -> str:
    lines = [line.rstrip() for line in body.replace("\r\n", "\n").split("\n")]
    while lines and not lines[0]:
        lines.pop(0)
    while lines and not lines[-1]:
        lines.pop()
    return "\n".join(lines)


def block_hash(persona: str, worktree: str, head: str, body: str) -> str:
    """The first 12 hex digits of SHA-256 over persona, worktree, head and the body (§3.2)."""
    text = f"{persona}\n{worktree}\n{head}\n{_normalise(body)}"
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def find_blocks(text: str) -> list[Block]:
    """Every block in `text`; a `MalformedError` for a marker that does not make one."""
    found: list[Block] = []
    opened: tuple[re.Match[str], list[str]] | None = None
    for line in text.replace("\r\n", "\n").split("\n"):
        header, end = HEADER.match(line), END.match(line)
        if header:
            if opened:
                raise MalformedError("a header inside a block")
            opened = header, []
        elif end:
            if opened is None or end.group(1) != opened[0].group(4):
                raise MalformedError("an END line without its header, or with another hash")
            persona, worktree, head, digest = opened[0].groups()
            found.append(Block(persona, worktree, head, digest, _normalise("\n".join(opened[1]))))
            opened = None
        elif "<<<BRIEF" in line or "<<<END BRIEF" in line:
            raise MalformedError(
                f"a marker line that is not a header or an END line: {line[:120]!r}"
            )
        elif opened:
            opened[1].append(line)
    if opened:
        raise MalformedError("a header without its END line")
    return found


def _git(path: Path | str, *args: str) -> str | None:
    done = subprocess.run(["git", "-C", str(path), *args], capture_output=True, text=True,
                          check=False)  # fmt: skip
    return done.stdout.strip() if done.returncode == 0 else None


def check(block: Block, root: Path) -> str | None:
    """Why `block` is refused (edited, not a checkout, stale), or None. `root` is
    the checkout asking, from which git runs."""
    if block_hash(block.persona, block.worktree, block.head, block.body) != block.hash:
        return "brief: edited; paste the output of tools/brief.py unchanged"
    head = _git(block.worktree, "rev-parse", "HEAD") if Path(block.worktree).is_dir() else None
    if head is None:
        return f"brief: worktree {block.worktree} is not a checkout; run brief.py again"
    if head != block.head:
        return f"brief: stale: {block.worktree} is at {head[:12]}, not {block.head[:12]}; " \
               "run brief.py again"  # fmt: skip
    return None


def _checkout(given: str) -> Path:
    """`given`, resolved, when it is the top of a checkout of this repository."""
    if any(c.isspace() for c in given):
        raise BriefError(f"the path {given!r} contains whitespace")
    path = Path(given).resolve()
    top = _git(path, "rev-parse", "--show-toplevel") if path.is_dir() else None
    if top is None or Path(top).resolve() != path:
        raise BriefError(f"{given} is not the top of a checkout of this repository")
    return path


def _inside(path: str, prefix: str) -> bool:
    return path == prefix or path.startswith(prefix.rstrip("/") + "/")


def _area(worktree: str, declared: list[str]) -> set[str]:
    """Pipe 3.6: what the branch changes against its base, plus `declared` before its
    first commit. The base is where it leaves origin/master (master without a remote)."""
    base = _git(worktree, "merge-base", "HEAD", "origin/master") or \
        _git(worktree, "merge-base", "HEAD", "master")  # fmt: skip
    changed = set((_git(worktree, "diff", "--name-only", base, "HEAD") or "").splitlines()
                  if base else ())  # fmt: skip
    return changed | set(declared) if base in (None, _git(worktree, "rev-parse", "HEAD")) \
        else changed  # fmt: skip


def _running(worktree: Path) -> tuple[Path, list[dict[str, str | list[str]]]]:
    """`areas.json` in the git common dir, and its entries for live writers elsewhere:
    the note file and the worktree still there."""
    store = Path(_git(worktree, "rev-parse", "--path-format=absolute", "--git-common-dir")
                 or worktree) / "areas.json"  # fmt: skip
    try:
        entries = json.loads(store.read_text())
    except (OSError, ValueError):
        entries = []
    return store, [e for e in entries if e["worktree"] != str(worktree)
                   and Path(e["note"]).is_file() and Path(e["worktree"]).is_dir()]  # fmt: skip


def _overlap(mine: dict[str, str | list[str]], running: list[dict[str, str | list[str]]]) -> None:
    """A refusal naming both writers when `mine` would write a running writer's paths."""

    def name(e: dict[str, str | list[str]]) -> str:
        return f"@{e['persona']} ({e['item'] or 'no item'}, in {e['worktree']})"

    for other in running:
        theirs = _area(str(other["worktree"]), list(other["declared"]))
        shared = sorted({b if _inside(b, a) else a for a in mine["area"] for b in theirs
                         if _inside(a, b) or _inside(b, a)})  # fmt: skip
        if shared:
            raise BriefError(f"{name(mine)} and {name(other)} would write the same files: "
                             f"{', '.join(shared)}")  # fmt: skip


def _concurrency(persona: str, worktree: Path, beside: list[str], no_build: bool,
                 untimed: bool = False, running: Sequence[str] = ()) -> str:  # fmt: skip
    """§3.1a: the generated line, or a refusal naming the rule."""
    if untimed and persona != "perf" and not any(e.partition(":")[0] == "perf" for e in beside):
        raise BriefError("--untimed marks a @perf run")
    too_many = BriefError("at most three writers run at once")
    if WRITES[persona] and len({str(worktree), *running}) > MAX_WRITERS:
        raise too_many
    if beside == ["none"]:
        alone = "Concurrency: you run alone"
        return f"{alone}; nothing runs beside a timing run." if persona == "perf" else f"{alone}."
    if "none" in beside:
        raise BriefError("--beside none goes alone, not with other --beside entries")
    others = []
    for entry in beside:
        name, _, place = entry.partition(":")
        if name not in PERSONAS or not place:
            raise BriefError(
                f"--beside {entry}: give <persona>:<worktree>, a persona of {PERSONAS}"
            )
        others.append((name, _checkout(place)))
    everyone = [(persona, worktree), *others]
    if any(name == "perf" for name, _ in everyone) and not untimed:
        raise BriefError(
            "nothing runs beside a timing run (@perf), read-only agents included; "
            "brief its untimed stages with --untimed"
        )
    writers = [place for name, place in everyone if WRITES[name]]
    if len({*map(str, writers), *running}) > MAX_WRITERS:
        raise too_many
    if len(set(writers)) < len(writers):
        raise BriefError("two writers share one worktree")
    for who, where in everyone:
        if not WRITES[who] and where in writers:
            raise BriefError(f"@{who} is read-only and shares {where} with a writer")
    kinds = {name: "writer" if WRITES[name] else "read-only" for name, _ in others}
    if untimed:
        kinds["perf"] = "writer, untimed run"
    listed = ", ".join(f"@{name} ({kinds[name]}, in {place})" for name, place in others)
    build = "You may not build C++ in this run." if no_build else "You may build C++ in this run."
    tail = " This run is untimed." if untimed and persona == "perf" else ""
    return (f"Concurrency: also running: {listed}. At most three writers run at once; read-only "
            "agents do not count toward the three, as long as no writer changes the files they "
            f"read.{tail}\n{build}")  # fmt: skip


def _from_increment(path: Path, shown: str) -> list[str]:
    """§3.1: the Status line and the quoted lines outside ## Review."""
    lines = path.read_text(errors="replace").splitlines()
    starts = [n for n, line in enumerate(lines) if line.startswith("## Review")]
    review = range(0)
    if starts:
        after = [n for n in range(starts[0] + 1, len(lines)) if lines[n].startswith("## ")]
        review = range(starts[0], after[0] if after else len(lines))
    out = [f"From {shown}, its own lines:"]
    out += [f"  {line}" for line in lines if line.startswith("Status:")][:1]
    hits = [(n + 1, line.strip()) for n, line in enumerate(lines)
            if n not in review and QUOTED.search(line)]  # fmt: skip
    out += [f"  {n}: {line}" for n, line in hits[:MAX_QUOTED]]
    if len(hits) > MAX_QUOTED:
        out.append(f"  ... {len(hits) - MAX_QUOTED} more; read them in the file")
    return out


def _ola(quotations: Sequence[str]) -> list[str]:
    """Each quotation, found in a human turn of this session's transcript, or a refusal."""
    session = os.environ.get("CLAUDE_CODE_SESSION_ID")
    path = session_state.TRANSCRIPTS / f"{session}.jsonl"
    turns = session_state.human_turns(path) if session and path.is_file() else []
    out = ["Ola, verbatim, checked against this session's transcript:"]
    for given in quotations:
        text = " ".join(given.split())
        stamp = next((when for when, _, said in turns if text in said), None)
        if stamp is None:
            raise BriefError(
                f"--ola {text!r} is not in any human turn of this session's transcript"
            )
        out.append(f'  "{text}" ({stamp})')
    return out


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:  # type: ignore[override]
        raise BriefError(message)


def _brief(argv: list[str] | None) -> str:
    parser = _Parser(prog="brief.py", description=__doc__)
    parser.add_argument("persona")
    parser.add_argument("--worktree", required=True)
    parser.add_argument("--beside", action="append", required=True)
    parser.add_argument("--increment")
    parser.add_argument("--no-build", action="store_true")
    parser.add_argument("--untimed", action="store_true")
    parser.add_argument("--ola", action="append", default=[])
    parser.add_argument("--lane", choices=("fast",))  # h34 §3: test and fix by one developer
    parser.add_argument("--item", default="")  # pipe 3.6: the pipe item a writer works on
    parser.add_argument("--area", nargs="+", default=[])  # narrows the write limit
    args = parser.parse_args(argv)
    persona, fast = args.persona, args.lane == "fast"
    if persona not in PERSONAS:
        raise BriefError(f"{persona} is not a persona; one of {', '.join(PERSONAS)}")
    if fast and persona != "developer":
        raise BriefError("--lane fast is for @developer only")
    worktree = _checkout(args.worktree)
    increment, parts = "none named", []
    if args.increment is not None:
        path = worktree / args.increment  # relative: the worktree's copy; absolute: as given
        if path.is_file():
            increment, parts = args.increment, _from_increment(path, args.increment)
        elif persona == "architect":
            increment = f"{args.increment} (new: you create it)"
        else:
            raise BriefError(f"--increment {args.increment} does not exist")
    elif persona in NEEDS_INCREMENT and not fast:
        raise BriefError(f"@{persona} needs --increment")
    # <persona>-<worktree>-<HHMMSS>.md, numbered -2, -3 while the name is taken (h16 T3)
    tasks = session_state.main_checkout(ROOT) / ".claude" / "current-task"
    stem = f"{persona}-{worktree.name}-{datetime.now().strftime('%H%M%S')}"
    note, n = tasks / f"{stem}.md", 1
    while note.exists():
        n += 1
        note = tasks / f"{stem}-{n}.md"
    try:
        fixed = string.Template(TEMPLATE.read_text()).substitute(
            persona=persona, worktree=worktree, increment=increment, note=note
        )
    except (KeyError, ValueError) as exc:
        raise BriefError(f"{TEMPLATE}: the template names an unknown variable ({exc})") from exc
    limits = WRITES[persona] + FAST_WRITES * fast
    prefixes = [e.split()[0] for e in limits if e.split()[0] not in ("root", "the")]
    for path in args.area:
        if not any(_inside(path, p) for p in prefixes):
            raise BriefError(f"--area {path} is outside the write limit")
    store, running = _running(worktree)
    declared = args.area or prefixes
    mine = {"persona": persona, "worktree": str(worktree), "item": args.item,
            "branch": _git(worktree, "rev-parse", "--abbrev-ref", "HEAD") or "", "note": str(note),
            "declared": declared, "area": sorted(_area(str(worktree), declared))}  # fmt: skip
    if WRITES[persona]:
        _overlap(mine, running)
    limit = ", ".join(limits) or "nothing"
    body = [fixed, _concurrency(persona, worktree, args.beside, args.no_build, args.untimed,
                                [str(e["worktree"]) for e in running]),
            f"Write limit: {limit}, and your note file.", *parts]  # fmt: skip
    if persona == "reviewer":
        body.insert(3, REVIEW_LINE)
    if args.ola:
        body += _ola(args.ola)
    text, head = _normalise("\n".join(body)), _git(worktree, "rev-parse", "HEAD") or ""
    digest = block_hash(persona, str(worktree), head, text)
    if WRITES[persona]:  # pipe 3.6: the note file marks the writer live until it is deleted
        note.parent.mkdir(parents=True, exist_ok=True)
        note.write_text(f"persona: @{persona}; made by tools/brief.py, yours to overwrite\n")
        temporary = store.with_name(f"areas.json.{os.getpid()}")
        temporary.write_text(json.dumps([*running, mine], indent=1))
        temporary.replace(store)
    return f"<<<BRIEF persona={persona} worktree={worktree} head={head} hash={digest}>>>\n" \
           f"{text}\n<<<END BRIEF {digest}>>>\n"  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    try:
        sys.stdout.write(_brief(argv))
    except BriefError as exc:
        print(f"brief.py: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
