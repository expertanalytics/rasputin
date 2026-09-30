#!/usr/bin/env python3
"""U0 probe (b): does a command run with `!` at the prompt get a terminal?

docs/increments/h3-unattended-u1.md, section 2. Not production code.
Prints the result and appends it to the probe log (hook_event_name "tty-probe").
Exit 0: /dev/tty opened and "y" was read. Exit 1: opened, no "y" within 20 s.
Exit 3: /dev/tty could not be opened (no terminal).
"""

import errno
import json
import os
import select
import sys
from datetime import UTC, datetime
from pathlib import Path

LOG = Path("/Users/skavhaug/projects/rasputin_scratch/u0-probe/log.jsonl")


def main() -> int:
    result: dict[str, object] = {
        "t": datetime.now(UTC).isoformat(timespec="seconds"),
        "hook_event_name": "tty-probe",
        "stdin_isatty": sys.stdin.isatty(),
    }
    code = 3
    try:
        fd = os.open("/dev/tty", os.O_RDWR)
    except OSError as err:
        result["dev_tty"] = errno.errorcode.get(err.errno or 0, str(err))
    else:
        result["dev_tty"] = "opened"
        os.write(fd, b"U0: type y and Enter within 20 s\n")
        ready, _, _ = select.select([fd], [], [], 20)
        answer = os.read(fd, 64).decode(errors="replace").strip() if ready else None
        os.close(fd)
        result["read"] = answer
        code = 0 if answer and answer.lower() == "y" else 1
    result["exit"] = code
    print(json.dumps(result))
    try:
        LOG.parent.mkdir(parents=True, exist_ok=True)
        with LOG.open("a") as out:
            out.write(json.dumps(result) + "\n")
    except OSError as err:
        print(f"could not write {LOG}: {err}", file=sys.stderr)
    return code


if __name__ == "__main__":
    sys.exit(main())
