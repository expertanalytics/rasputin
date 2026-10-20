"""The `@orchestrator` runs owed, computed from git and the away log, as one line.

Run with no arguments as the `UserPromptSubmit` hook (the push item left out),
with `--status` as the status line (the session's JSON on stdin), and imported
by `tools/session_state.py` for the recap. It writes nothing and always exits
0: a hook exiting 2 would block the prompt, so any fault becomes a
`Duties: unknown` line instead, never silence, which would read as "none owed".
Spec: `docs/increments/h23-duties.md`, *2. The duties, and how each is computed*.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path
from typing import Any, Literal, TextIO

sys.path.append(str(Path(__file__).resolve().parent))  # appended: the stdlib wins
import harness_mode as hm

ROOT = Path(__file__).resolve().parents[1]
#: Chosen: an away window this long owes a check whether it ran by night or by day.
NIGHT = timedelta(hours=4)
WEEK = timedelta(days=7)
RETROS = "docs/retrospectives/"
#: A check commit is a commit by `@orchestrator`, known by this subject suffix.
CHECK = "(@orchestrator)"
Place = Literal["hook", "recap", "status"]
# The git argument vectors, one per query; the long two are packed so each reads as one command.
GIT = ("git", "-C")
ORIGIN = ("rev-parse", "--verify", "--quiet", "refs/remotes/origin/master")
CHECKS = ("log", "--diff-filter=A", "--name-only", "--format=C %cI %s", "--branches")
MERGES = ("log", "--first-parent", "--merges", "--diff-merges=first-parent",
          "--name-only", "--format=M %h %cI %s")  # fmt: skip
UNMERGED = ("log", "--no-merges", "--diff-filter=A", "--source",
            "--name-only", "--format=S %S", "--branches", "--not")  # fmt: skip


@dataclass(frozen=True)
class Window:
    start: datetime
    end: datetime


@dataclass(frozen=True)
class Merge:
    when: datetime
    label: str  # `#N` from a pull request merge's subject, else the short sha
    paths: tuple[str, ...]


@dataclass(frozen=True)
class Inputs:
    retro: datetime | None  # the latest first add of a retrospective by a check commit
    research: datetime | None  # the same, over files whose name contains `-research`
    merges: list[Merge]  # first-parent merges on the base, newest first
    windows: list[Window]
    unmerged: list[str]  # local branches adding a retrospective the base lacks


def window_ends(events: list[dict[str, Any]], now: datetime) -> list[Window]:
    """The away windows ended by `now`, in log order; a malformed event is skipped.

    A `back` ends the window it names; an `enter` with no `back` ends at its own
    `until`, or at the next `enter`'s `since` when that is earlier.
    """
    ended: list[Window] = []
    pending: Window | None = None
    for event in events:
        try:
            kind = event["event"]
            since, until = (datetime.fromisoformat(event[k]) for k in ("since", "until"))
        except (KeyError, TypeError, ValueError):
            continue
        if kind not in ("enter", "back"):
            continue
        if pending and kind == "enter":
            ended.append(Window(pending.start, min(pending.end, since)))
        elif kind == "back":
            ended.append(Window(since, until))
        pending = Window(since, until) if kind == "enter" else None
    ended += [pending] if pending else []
    return [window for window in ended if window.end <= now]


def _minute(moment: datetime) -> str:
    return moment.astimezone(UTC).strftime("%Y-%m-%d %H:%M UTC")


def _recording(paths: tuple[str, ...]) -> bool:
    """A merge that records a check: a retrospective, beside ROADMAP and increment files only."""
    other = ("ROADMAP.md", "docs/increments/", RETROS)
    return any(p.startswith(RETROS) for p in paths) and all(p.startswith(other) for p in paths)


def _listed(names: list[str], shown: int) -> str:
    extra = len(names) - shown
    return ", ".join(names[:shown]) + (f", +{extra}" if extra > 0 else "")


def owed(inputs: Inputs, now: datetime) -> list[str]:
    """The owed items, in output order: post-merge, morning, research, push."""
    retro, items = inputs.retro, []
    work = [
        m for m in inputs.merges if (retro is None or m.when > retro) and not _recording(m.paths)
    ]
    if work and retro is None:
        items.append("post-merge check (no retrospective on record)")
    elif work and retro is not None:
        labels = _listed([m.label for m in work], 3)
        plural = "s" * (len(work) > 1)
        since = f"since the last retrospective, {_minute(retro)}"
        items.append(f"post-merge check ({len(work)} merge{plural} {since}: {labels})")
    nights = [
        w
        for w in inputs.windows
        if w.end <= now and w.end - w.start >= NIGHT and (retro is None or w.end > retro)
    ]
    if nights:
        latest = _minute(max(w.end for w in nights))
        plural = "s" * (len(nights) > 1)
        items.append(
            f"morning check ({len(nights)} away window{plural} ended since, latest {latest})"
        )
    if inputs.research is None:
        items.append("research round (none on record)")
    elif now - inputs.research >= WEEK:
        items.append(f"research round (last {inputs.research.astimezone(UTC).date().isoformat()})")
    names = sorted(set(inputs.unmerged))
    if names:
        items.append(
            f"push {len(names)} retrospective branch{'es' * (len(names) > 1)} ({_listed(names, 2)})"
        )
    return items


def render(items: list[str], place: Place) -> str | None:
    """The one line for `place`: the hook drops the push item, which only Ola can act on."""
    if place == "hook":
        items = [item for item in items if not item.startswith("push ")]
    if items:
        return "Duties owed: " + "; ".join(items)
    return "Duties: none owed" if place == "status" else None


def _git(repo: Path, *args: str) -> str:
    result = subprocess.run([*GIT, str(repo), *args], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"git {args[0]}: {result.stderr.strip()}")
    return result.stdout


def _blocks(text: str, mark: str) -> list[tuple[str, list[str]]]:
    """`git log --name-only` output, as (header after `mark`, its paths) pairs."""
    blocks: list[tuple[str, list[str]]] = []
    for line in text.splitlines():
        if line.startswith(mark):
            blocks.append((line[len(mark) :], []))
        elif line and blocks:
            blocks[-1][1].append(line)
    return blocks


def gather(repo: Path, unmerged: bool) -> Inputs:
    """Read git and the away log; the unmerged-branch query only when `unmerged`."""
    origin = subprocess.run([*GIT, str(repo), *ORIGIN], capture_output=True, check=False)
    base, spec = ("origin/master" if origin.returncode == 0 else "master"), RETROS + "*.md"
    first: dict[str, datetime] = {}
    for header, paths in _blocks(_git(repo, *CHECKS, base, "--", spec), "C "):
        stamp, _, subject = header.partition(" ")
        for path in paths if subject.endswith(CHECK) else []:
            added = datetime.fromisoformat(stamp)
            first[path] = min(first.get(path, added), added)
    retro = max(first.values(), default=None)
    research = max((t for p, t in first.items() if "-research" in Path(p).name), default=None)
    since = [f"--since={retro.isoformat()}"] if retro else []
    merges = []
    for header, paths in _blocks(_git(repo, *MERGES, *since, base), "M "):
        sha, stamp, subject = [*header.split(" ", 2), ""][:3]
        number = re.match(r"Merge pull request (#\d+)", subject)
        label = number[1] if number else sha
        merges.append(Merge(datetime.fromisoformat(stamp), label, tuple(paths)))
    names: set[str] = set()
    if unmerged:
        held = set(_git(repo, "ls-tree", "--name-only", base, RETROS).splitlines())
        for header, paths in _blocks(_git(repo, *UNMERGED, base, "--", spec), "S "):
            if set(paths) - held:
                names.add(header.removeprefix("refs/heads/"))
    state = hm.state_dir(repo)
    log = state / "windows.jsonl" if state else None
    lines = log.read_text().splitlines() if log and log.exists() else []
    events = [e for e in map(_json, lines) if isinstance(e, dict)]
    return Inputs(retro, research, merges, window_ends(events, datetime.now(UTC)), sorted(names))


def _json(line: str) -> Any:
    try:
        return json.loads(line)
    except json.JSONDecodeError:
        return None


def line(repo: Path, place: Place) -> str | None:
    """The duties line for `place`, read from the repository at `repo` now."""
    return render(owed(gather(repo, place != "hook"), datetime.now(UTC)), place)


def main(argv: list[str], stdin: TextIO) -> int:
    """Print the line for the hook, or with `--status` for the status line; always 0."""
    try:
        if argv[:1] == ["--status"]:
            text = line(Path(json.loads(stdin.read())["workspace"]["project_dir"]), "status")
        else:
            text = line(ROOT, "hook")
    except Exception as error:
        message = str(error).replace("\n", " ")[:120]
        text = f"Duties: unknown ({type(error).__name__}: {message})"
    if text:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:], sys.stdin))
