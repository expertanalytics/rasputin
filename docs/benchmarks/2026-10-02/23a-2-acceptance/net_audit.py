"""Run a script under an audit hook that counts network events.

Usage: python net_audit.py [--control] SCRIPT ARGS...
Runs SCRIPT as __main__ and prints, on stderr at exit,
NET_EVENTS urllib.Request=<n> socket.connect=<n> socket.__new__=<n> and the
range of every urllib request (from its Range header), so the requests a
fetch makes are counted outside the code under test.
--control opens and connects nothing but creates one socket first, so a hook
that never fires is caught.
"""

import atexit
import collections
import runpy
import sys

counts: collections.Counter[str] = collections.Counter()
ranges: list[str] = []


def _hook(name: str, args: tuple[object, ...]) -> None:
    if name.startswith("socket.") or name == "urllib.Request":
        counts[name] += 1
    if name == "urllib.Request":
        headers = args[2] if len(args) > 2 else {}
        if isinstance(headers, dict):
            ranges.append(str(headers.get("Range", headers.get("range", "-"))))


def _report() -> None:
    keys = ("urllib.Request", "socket.connect", "socket.__new__")
    line = " ".join(f"{k}={counts.get(k, 0)}" for k in keys)
    print(f"NET_EVENTS {line}", file=sys.stderr)
    for r in ranges:
        print(f"RANGE {r}", file=sys.stderr)


sys.addaudithook(_hook)
atexit.register(_report)
argv = sys.argv[1:]
if argv and argv[0] == "--control":
    import socket

    socket.socket().close()
    argv = argv[1:]
sys.argv = argv
runpy.run_path(argv[0], run_name="__main__")
