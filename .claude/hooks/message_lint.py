#!/usr/bin/env python3
"""Stop: send the main session's last message back when it is not self-contained.

Three mechanical rules, kept conservative so a false block stays rare:

1. A repository path is absolute. `tools/x.py`, `docs/...`, `.claude/...` and
   the other top-level directories, or a bare file name with a line number, are flagged unless
   an absolute prefix leads to them. Checked in prose and inline code spans;
   fenced code blocks and URLs are skipped.
2. Work is not named by a bare number: `increment 150`, `design 180`,
   `item 2.58`, `q421`, `the 150 tool`.
3. A short run label is not a name: `P`, `L`, `PL`, `A14`, `HC15`, `F91`.

Rules 2 and 3 skip code spans, code blocks and URLs, and allow a number or
label inside parentheses or square brackets, where it follows a description.

Subagents (an `agent_id` in the event) pass untouched. A message is sent back
at most once (its hash is kept in `<git-common-dir>/harness/message_lint.jsonl`),
and never more than `STREAK` times in a row, so a session cannot loop. A crash,
or a hash that cannot be recorded, lets the stop through.
"""

import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

STREAK = 3
LEDGER = "message_lint.jsonl"

TOP_DIRS = r"(?:docs|src|src_python|tests|tools|include|bindings|\.claude|\.github)"
REPO_PATH = re.compile(rf"(?<![\w./~\-]){TOP_DIRS}/[\w./\-]*\.\w+(?::\d+)?")
BARE_LINE = re.compile(
    r"(?<![\w./~\-])[\w\-]+\.(?:py|cpp|hpp|h|md|toml|json|jsonl|yaml|yml|txt|sh|cmake):\d+"
)
WORK_NUMBER = re.compile(
    r"\b(?:[Ii]ncrement \d+|[Dd]esign \d+|[Ii]tem \d+\.\d+|[qQ]\d{3,4}"
    r"|the \d{2,3} (?:tool|record|design|run|probe|increment))\b"
)
RUN_LABEL = re.compile(r"(?<![\w\-./+':])(?:PL|P|L|[A-Z]{1,2}\d{1,3})(?![\w\-/'+]|\s*[=<>])")
#: Units, standards, map products (AR5, AR50, N50) and common terms, never run names.
ALLOWED_LABELS = re.compile(
    r"P\d0|P9[59]|L[12]|CO2|Q[1-4]|A[3-5]|M[1-5]|B\d\d|H[1-6]|AR\d{1,3}|N\d{2,4}"
)

FENCE = re.compile(r"^\s*(```|~~~).*?^\s*\1[^\n]*$", re.MULTILINE | re.DOTALL)
CODE_SPAN = re.compile(r"`[^`\n]+`")
URL = re.compile(r"\b(?:https?|file)://\S+")

RULES = {
    "path": "a relative path; write it absolute (starting /Users/ or /private/)",
    "number": "work named by a bare number; say what it is, the number in parentheses after",
    "label": "a short run label used as a name; say what the run is, the label in parentheses",
}


def blank(pattern: re.Pattern[str], text: str) -> str:
    """`text` with each match of `pattern` replaced by as many spaces (offsets kept)."""
    return pattern.sub(lambda m: re.sub(r"[^\n]", " ", m.group(0)), text)


def bracketed(line: str, at: int) -> bool:
    """Whether offset `at` of `line` lies inside parentheses or square brackets."""
    depth = 0
    for char in line[:at]:
        if char in "([":
            depth += 1
        elif char in ")]":
            depth = max(0, depth - 1)
    return depth > 0


def findings(text: str) -> list[tuple[str, str]]:
    """(fragment, rule) for each violation in `text`, in order, each fragment once."""
    prose = blank(URL, blank(FENCE, text))
    found: list[tuple[str, str]] = []
    for line in prose.splitlines():
        for pattern in (REPO_PATH, BARE_LINE):
            found += [(m.group(0), "path") for m in pattern.finditer(line)]
    for line in blank(CODE_SPAN, prose).splitlines():
        for pattern, rule in ((WORK_NUMBER, "number"), (RUN_LABEL, "label")):
            for m in pattern.finditer(line):
                if bracketed(line, m.start()):
                    continue
                if rule == "label" and ALLOWED_LABELS.fullmatch(m.group(0)):
                    continue
                found.append((m.group(0), rule))
    return list(dict.fromkeys(found))


def last_message(event: dict[str, Any]) -> str:
    """The text of the main session's last assistant message."""
    given = event.get("last_assistant_message")
    if isinstance(given, str):
        return given
    path = event.get("transcript_path")
    if not isinstance(path, str):
        return ""
    entries = []
    for raw in Path(path).read_text().splitlines():
        try:
            entries.append(json.loads(raw))
        except ValueError:
            continue
    assistant = [e for e in entries if isinstance(e, dict) and e.get("type") == "assistant"]
    if not assistant:
        return ""
    last_id = (assistant[-1].get("message") or {}).get("id")
    texts = []
    for entry in assistant:
        message = entry.get("message") or {}
        if message.get("id") != last_id:
            continue
        content = message.get("content")
        if isinstance(content, str):
            texts.append(content)
        elif isinstance(content, list):
            texts += [b.get("text", "") for b in content if b.get("type") == "text"]
    return "\n".join(texts)


def reason(found: list[tuple[str, str]]) -> str:
    lines = [f"- `{fragment}`: {RULES[rule]}" for fragment, rule in found]
    return (
        "message_lint: this message is not self-contained, as Ola has ruled messages must "
        "be. Rewrite it and send it again, fixing each fragment:\n" + "\n".join(lines)
    )


def settle(event: dict[str, Any]) -> dict[str, str] | None:
    """The block for a main-session stop whose last message breaks a rule, or None."""
    if event.get("hook_event_name") != "Stop" or event.get("agent_id"):
        return None
    text = last_message(event)
    found = findings(text)
    if not found:
        return None
    sys.path.append(str(Path(__file__).resolve().parents[2] / "tools"))  # the stdlib wins
    import harness_mode

    state = harness_mode.state_dir(harness_mode.ROOT)
    if state is None:
        return None
    sha = hashlib.sha256(text.encode()).hexdigest()
    session = event.get("session_id")
    seen = []
    if (state / LEDGER).exists():
        for raw in (state / LEDGER).read_text().splitlines():
            try:
                seen.append(json.loads(raw))
            except ValueError:
                continue
    if any(entry.get("sha") == sha for entry in seen):
        return None
    tail = seen[-STREAK:]
    if event.get("stop_hook_active") and len(tail) == STREAK and all(
        entry.get("session_id") == session for entry in tail
    ):
        return None
    entry = {"session_id": session, "sha": sha}
    if harness_mode.append_queue(state, entry, LEDGER) is not None:
        return None
    return {"decision": "block", "reason": reason(found)}


def main() -> int:
    try:
        event = json.load(sys.stdin)
        output = settle(event) if isinstance(event, dict) else None
    except Exception as error:  # a crash lets the stop through
        print(f"message_lint failed: {type(error).__name__}: {error}", file=sys.stderr)
        return 0
    if output is not None:
        print(json.dumps(output))
    return 0


if __name__ == "__main__":
    sys.exit(main())
