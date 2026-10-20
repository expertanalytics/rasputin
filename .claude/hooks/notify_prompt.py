#!/usr/bin/env python3
"""Notification: a Mac banner when a permission prompt waits unanswered.

Claude Code sends the `permission_prompt` notification about six seconds after
a prompt opened, so Ola sees it with the terminal hidden. The message and the
worktree name reach osascript as arguments, never spliced into the script.
Always exits 0: Claude Code ignores this hook's exit code and stderr.
"""

import contextlib
import json
import subprocess
import sys
from pathlib import Path

OSASCRIPT = "/usr/bin/osascript"
SCRIPT = ("on run argv", "display notification (item 1 of argv) with title "
          '"rasputin: a prompt is waiting" subtitle (item 2 of argv)', "end run")  # fmt: skip


def notify(event: dict, run=subprocess.run, platform=sys.platform) -> list[str] | None:
    """The osascript argv for a permission_prompt notification, run once; else None."""
    if event.get("notification_type") != "permission_prompt" or platform != "darwin":
        return None
    if not Path(OSASCRIPT).exists():
        return None
    argv = [OSASCRIPT, *(part for line in SCRIPT for part in ("-e", line)),
            str(event.get("message", "")), Path(str(event.get("cwd", ""))).name]  # fmt: skip
    with contextlib.suppress(OSError, subprocess.TimeoutExpired):
        run(argv, timeout=5, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    return argv


def main() -> int:
    try:
        event = json.load(sys.stdin)
    except ValueError:
        return 0
    if isinstance(event, dict):
        notify(event)
    return 0


if __name__ == "__main__":
    sys.exit(main())
