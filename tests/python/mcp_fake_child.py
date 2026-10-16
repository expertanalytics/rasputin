"""A stand-in for ``python -m tin_engine`` under the MCP server's child runner.

Increment 37 (`docs/increments/37-mcp-server.md`, sections 6.4 and 9). The
server runs ``[*child, <command>, <options>...]``; a test gives it
``[sys.executable, this file, SPEC.json]`` as ``child``, so ``sys.argv[2:]`` is
exactly the argv the server built from the tool arguments. The spec, written
by `mcp_driver.fake_child`, says what to do, in this order:

``argv_file``
    write ``sys.argv[2:]`` here as JSON (T4: what the server passed).
``pid_file``
    write this process's pid (T8: is it dead after a cancel?).
``create_record``
    write ``{}`` to the ``--record`` path found in argv (T5: a sidecar the
    child did create, which the server must still remove).
``stdout``
    lines to print to stdout, at once.
``stderr``
    ``[seconds, line]`` pairs: sleep, then print the line to stderr, flushed
    (T9: each line is one progress notification).
``sleep``
    seconds to sleep before exiting (T8: a run to cancel; T10: a run to queue).
``interval_file``
    append ``{"start": t, "end": t}`` wall-clock seconds, one JSON line (T10:
    did two runs overlap?).
``signal``
    end by sending this signal number to itself (T5: ``killed by signal``).
``exit``
    the exit code, default 0.

Test support only; nothing in `src_python` runs it.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path


def main() -> None:
    spec = json.loads(Path(sys.argv[1]).read_text())
    argv = sys.argv[2:]
    started = time.time()
    if "argv_file" in spec:
        Path(spec["argv_file"]).write_text(json.dumps(argv))
    if "pid_file" in spec:
        Path(spec["pid_file"]).write_text(str(os.getpid()))
    if spec.get("create_record") and "--record" in argv:
        Path(argv[argv.index("--record") + 1]).write_text("{}")
    for line in spec.get("stdout", []):
        print(line, flush=True)
    for seconds, line in spec.get("stderr", []):
        time.sleep(seconds)
        print(line, file=sys.stderr, flush=True)
    time.sleep(spec.get("sleep", 0.0))
    if "interval_file" in spec:
        with open(spec["interval_file"], "a") as f:
            f.write(json.dumps({"start": started, "end": time.time()}) + "\n")
    if "signal" in spec:
        os.kill(os.getpid(), int(spec["signal"]))
        time.sleep(5.0)  # never reached: the signal ends the process
    sys.exit(int(spec.get("exit", 0)))


if __name__ == "__main__":
    main()
