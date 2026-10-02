"""Run a command line under an audit hook that records every socket event.

Usage: python socket_audit.py [--control] SCRIPT ARGS...
Runs SCRIPT as __main__ and prints SOCKET_EVENTS <n> on stderr at exit.
--control opens one socket first, so a hook that never fires is caught.
"""

import atexit
import runpy
import sys

events: list[str] = []


def _hook(name: str, args: tuple[object, ...]) -> None:
    if name.startswith("socket."):
        events.append(name)


sys.addaudithook(_hook)
atexit.register(lambda: print(f"SOCKET_EVENTS {len(events)} {sorted(set(events))}", file=sys.stderr))
argv = sys.argv[1:]
if argv and argv[0] == "--control":
    import socket

    socket.socket().close()
    argv = argv[1:]
sys.argv = argv
runpy.run_path(argv[0], run_name="__main__")
