"""Compare each cache's NOTICE.txt with sources.notice(SOURCES[id]).

Usage: python check_notice.py CACHE_ROOT... (each holding one source directory)
"""

import sys
from pathlib import Path

from tin_engine.sources import SOURCES, notice

bad = 0
for root in map(Path, sys.argv[1:]):
    for d in sorted(p for p in root.iterdir() if p.is_dir()):
        on_disk = (d / "NOTICE.txt").read_text(encoding="utf-8")
        same = on_disk == notice(SOURCES[d.name])
        bad += not same
        print(("PASS " if same else "FAIL ") + f"{d.name}: NOTICE.txt equals notice(SOURCES[{d.name!r}])")
sys.exit(1 if bad else 0)
