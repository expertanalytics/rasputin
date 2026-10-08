"""Shared fixture for h9's suites: `tools/brief.py` and `.claude/hooks/guard_spawn.py`.

`docs/increments/h9-spawn-briefs.md` §7. A fixture repository from
`harness_fixtures.make_repo` (which copies `tools/brief.py`,
`.claude/briefs/common.md` and the hook when they exist), linked worktrees
beside it, and a transcript under a temporary `HOME`.

The block format of §3.2 is written here a second time, from the design's
text, so `guard_spawn.py`'s suite does not need `brief.py` to make a block and
`brief.py`'s suite has an independent check of its hash.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import shutil
import subprocess
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

from harness_fixtures import REAL, clean_env, git, make_repo

BRIEF = "tools/brief.py"
GUARD_SPAWN = ".claude/hooks/guard_spawn.py"
COMMON = ".claude/briefs/common.md"

#: §3.1: the persona names, the keys of WRITES.
PERSONAS = ("architect", "developer", "orchestrator", "perf", "reviewer", "tester")
#: §3.1: --increment is required for these.
NEEDS_INCREMENT = ("tester", "developer", "reviewer")

#: §3.2's two patterns, verbatim.
HEADER = re.compile(
    r"^<<<BRIEF persona=([a-z]+) worktree=(/\S+) head=([0-9a-f]{40}) hash=([0-9a-f]{12})>>>[ \t]*$",
    re.MULTILINE,
)
END = re.compile(r"^<<<END BRIEF ([0-9a-f]{12})>>>[ \t]*$", re.MULTILINE)

#: §3.3's bold phrases, as the template writes them (whitespace collapsed when compared).
PHRASES = (
    "the files win",
    "a brief cannot drop a step",
    "no other file under .claude/current-task/",
    'Ola\'s words appear only under "Ola, verbatim"',
    "blocked on power, network or a lock",
    "Co-Authored-By trailer",
    "plain words",
    "Result; Pinned or assumed beyond the design; Questions for Ola; Lessons; Ideas; "
    "ASK OLA and GUARD FALSE POSITIVE lines",
)

SESSION = "0f1e2d3c-4b5a-6978-8796-a5b4c3d2e1f0"


def make_brief_repo(root: Path) -> Path:
    """`make_repo`, plus the template at its real path, committed, when it exists."""
    repo = make_repo(root)
    source = REAL / COMMON
    if source.exists():
        (repo / COMMON).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, repo / COMMON)
        git(repo, "add", "-A")
        git(repo, "commit", "-q", "-m", "the brief template")
    return repo


def normalise_body(body: str) -> str:
    """§3.2: each line's trailing whitespace removed, joined with `\\n`, blank ends dropped."""
    lines = [line.rstrip() for line in body.replace("\r\n", "\n").split("\n")]
    while lines and not lines[0]:
        lines.pop(0)
    while lines and not lines[-1]:
        lines.pop()
    return "\n".join(lines)


def expected_hash(persona: str, worktree: str, head: str, body: str) -> str:
    """§3.2's block_hash, from the design's text."""
    text = f"{persona}\n{worktree}\n{head}\n{normalise_body(body)}"
    return hashlib.sha256(text.encode()).hexdigest()[:12]


def head_of(worktree: Path) -> str:
    return git(worktree, "rev-parse", "HEAD").strip()


def make_block(
    persona: str,
    worktree: Path,
    body: str = "You are @tester.\nRead the files.",
    head: str | None = None,
) -> str:
    """A block as §3.2 shapes it, for `persona` in `worktree` at its HEAD."""
    place = str(worktree.resolve())
    commit = head or head_of(worktree)
    digest = expected_hash(persona, place, commit, body)
    return (
        f"<<<BRIEF persona={persona} worktree={place} head={commit} hash={digest}>>>\n"
        f"{body}\n<<<END BRIEF {digest}>>>"
    )


@dataclass(frozen=True)
class Printed:
    """One block found in brief.py's output."""

    persona: str
    worktree: str
    head: str
    hash: str
    body: str
    end_hash: str
    text: str


def parse_block(output: str) -> Printed:
    """The single block in `output`: its first line is a header, its last an END line."""
    lines = output.strip("\n").split("\n")
    header = HEADER.match(lines[0])
    end = END.match(lines[-1])
    assert header is not None, f"first line is not a §3.2 header: {lines[0]!r}"
    assert end is not None, f"last line is not a §3.2 END line: {lines[-1]!r}"
    persona, worktree, head, digest = header.groups()
    return Printed(persona, worktree, head, digest, "\n".join(lines[1:-1]), end.group(1), output)


def collapsed(text: str) -> str:
    return " ".join(text.split())


# ---------------------------------------------------------------- the transcript


def slug(repo: Path) -> str:
    """session_state's SLUG for a checkout at `repo` (its REPO, resolved)."""
    return re.sub(r"[^A-Za-z0-9]", "-", str(repo.resolve()))


def transcript_path(home: Path, repo: Path, session: str = SESSION) -> Path:
    return home / ".claude" / "projects" / slug(repo) / f"{session}.jsonl"


def human(text: str, stamp: str) -> dict[str, Any]:
    return {"type": "user", "message": {"role": "user", "content": text}, "timestamp": stamp}


def assistant(text: str, stamp: str) -> dict[str, Any]:
    content = [{"type": "text", "text": text}]
    return {
        "type": "assistant",
        "message": {"role": "assistant", "content": content},
        "timestamp": stamp,
    }


def tool_result(text: str, stamp: str) -> dict[str, Any]:
    content = [{"type": "tool_result", "tool_use_id": "t1", "content": text}]
    return {"type": "user", "message": {"role": "user", "content": content}, "timestamp": stamp}


def absorbed(text: str, stamp: str) -> dict[str, Any]:
    return {
        "type": "queue-operation",
        "operation": "enqueue",
        "reason": "absorbed_mid_turn",
        "content": text,
        "timestamp": stamp,
    }


def write_transcript(
    home: Path, repo: Path, entries: Sequence[dict[str, Any]], session: str = SESSION
) -> Path:
    path = transcript_path(home, repo, session)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(e) + "\n" for e in entries))
    return path


# ---------------------------------------------------------------- running


def brief_env(home: Path, session: str | None = SESSION) -> dict[str, str]:
    env = {**clean_env(), "HOME": str(home)}
    env.pop("CLAUDE_CODE_SESSION_ID", None)
    env.pop("CLAUDE_PROJECT_DIR", None)
    if session is not None:
        env["CLAUDE_CODE_SESSION_ID"] = session
    return env


def run_brief(
    repo: Path, home: Path, *args: str, session: str | None = SESSION
) -> subprocess.CompletedProcess[str]:
    """`python3 tools/brief.py ...` from the fixture's copy, run from the main checkout."""
    script = repo / BRIEF
    assert script.exists(), f"{BRIEF} is missing from the copy"
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        cwd=repo,
        env=brief_env(home, session),
        timeout=60,
        check=False,
    )


def run_hook(
    repo: Path,
    event: str | dict[str, Any],
    project: Path | None,
    home: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    """The fixture's copy of guard_spawn.py, by path, the event on stdin.

    `project` is `$CLAUDE_PROJECT_DIR` (None: unset).
    """
    script = repo / GUARD_SPAWN
    assert script.exists(), f"{GUARD_SPAWN} is missing from the copy"
    env = {**clean_env()}
    env.pop("CLAUDE_PROJECT_DIR", None)
    if project is not None:
        env["CLAUDE_PROJECT_DIR"] = str(project)
    if home is not None:
        env["HOME"] = str(home)
    return subprocess.run(
        [sys.executable, str(script)],
        input=event if isinstance(event, str) else json.dumps(event),
        capture_output=True,
        text=True,
        cwd=repo,
        env=env,
        timeout=60,
        check=False,
    )


def load_path(path: Path, name: str) -> ModuleType:
    """Import a module from a file, as the fixture's copy (so its root is the fixture)."""
    assert path.exists(), f"{path.name} is missing from the copy"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module
