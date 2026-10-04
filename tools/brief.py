#!/usr/bin/env python3
"""The fixed part of a persona's brief, as one hashed block (h9).

Spec: docs/increments/h9-spawn-briefs.md §3.1 to §3.4. The main session runs

    python3 tools/brief.py <persona> --worktree <path> --beside <none|persona:path>...
                           [--increment <file>] [--no-build] [--ola <text>]...

and pastes the block, unchanged, before the task. `.claude/hooks/guard_spawn.py`
reads it back with `find_blocks` and `check`, so the format has one owner. A
refusal is exit 2 with the reason on stderr and no block.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import string
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import session_state  # beside this file, as the hooks import it

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE = ROOT / ".claude" / "briefs" / "common.md"

#: The write limit per persona, h6's table (docs/increments/h6-role-limits.md §3).
#: Read-only is derived from an empty entry (§3.1a). When h6 lands, its ROLES replace this.
WRITES: dict[str, tuple[str, ...]] = {
    "architect": ("docs/ except docs/retrospectives/", "ROADMAP.md", "CLAUDE.md",
                  ".claude/ files ending in .md"),
    "developer": ("src_python/", "include/", "src/", "bindings/", "tools/", ".claude/hooks/",
                  ".github/", "CMakeLists.txt", "pyproject.toml"),
    "orchestrator": ("docs/retrospectives/",),
    "perf": ("docs/benchmarks/",),
    "reviewer": (),
    "tester": ("tests/",),
}  # fmt: skip
PERSONAS = tuple(WRITES)
NEEDS_INCREMENT = ("tester", "developer", "reviewer")
QUOTED = re.compile(r"(?i)invariant-critical|mutation|@perf|acceptance run|tools/bench\.py")
HEADER = re.compile(
    r"^<<<BRIEF persona=([a-z]+) worktree=(/\S+) head=([0-9a-f]{40}) hash=([0-9a-f]{12})>>>[ \t]*$"
)
END = re.compile(r"^<<<END BRIEF ([0-9a-f]{12})>>>[ \t]*$")
MAX_QUOTED, MAX_VERDICT = 8, 200


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


def _concurrency(persona: str, worktree: Path, beside: list[str], no_build: bool) -> str:
    """§3.1a: the generated line, or a refusal naming the rule."""
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
    if any(name == "perf" for name, _ in everyone):
        raise BriefError("nothing runs beside a timing run (@perf), read-only agents included")
    writers = [place for name, place in everyone if WRITES[name]]
    if len(writers) > 2:
        raise BriefError("at most two writers run at once")
    if len(set(writers)) < len(writers):
        raise BriefError("two writers share one worktree")
    for name, place in everyone:
        if not WRITES[name] and place in writers:
            raise BriefError(f"@{name} is read-only and shares {place} with a writer")
    kinds = {name: "writer" if WRITES[name] else "read-only" for name, _ in others}
    listed = ", ".join(f"@{name} ({kinds[name]}, in {place})" for name, place in others)
    build = "You may not build C++ in this run." if no_build else "You may build C++ in this run."
    return (f"Concurrency: also running: {listed}. At most two writers run at once; read-only "
            "agents do not count toward the two, as long as no writer changes the files they "
            f"read.\n{build}")  # fmt: skip


def _from_increment(path: Path, shown: str) -> list[str]:
    """§3.1a: the Status line, the quoted lines outside ## Review, the review count."""
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
    rounds = [lines[n] for n in review if "APPROVED" in lines[n] or "CHANGES REQUESTED" in lines[n]]
    if not starts:
        out.append("Review rounds recorded: none")
    else:
        last = f"; last: {rounds[-1][:MAX_VERDICT]}" if rounds else ""
        out.append(f"Review rounds recorded: {len(rounds)}{last}")
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
    parser.add_argument("--ola", action="append", default=[])
    args = parser.parse_args(argv)
    persona = args.persona
    if persona not in PERSONAS:
        raise BriefError(f"{persona} is not a persona; one of {', '.join(PERSONAS)}")
    worktree = _checkout(args.worktree)
    increment, parts = "none named", []
    if args.increment is not None:
        path = ROOT / args.increment
        if path.is_file():
            increment, parts = args.increment, _from_increment(path, args.increment)
        elif persona == "architect":
            increment = f"{args.increment} (new: you create it)"
        else:
            raise BriefError(f"--increment {args.increment} does not exist")
    elif persona in NEEDS_INCREMENT:
        raise BriefError(f"@{persona} needs --increment")
    stamp = datetime.now().strftime("%H%M%S")
    note = session_state.main_checkout(ROOT) / ".claude" / "current-task" / f"{persona}-{stamp}.md"
    try:
        fixed = string.Template(TEMPLATE.read_text()).substitute(
            persona=persona, worktree=worktree, increment=increment, note=note
        )
    except (KeyError, ValueError) as exc:
        raise BriefError(f"{TEMPLATE}: the template names an unknown variable ({exc})") from exc
    limit = ", ".join(WRITES[persona]) or "nothing"
    body = [fixed, _concurrency(persona, worktree, args.beside, args.no_build),
            f"Write limit: {limit}, and your note file.", *parts]  # fmt: skip
    if args.ola:
        body += _ola(args.ola)
    text, head = _normalise("\n".join(body)), _git(worktree, "rev-parse", "HEAD") or ""
    digest = block_hash(persona, str(worktree), head, text)
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
