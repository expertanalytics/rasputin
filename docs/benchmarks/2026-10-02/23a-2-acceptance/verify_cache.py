"""Check a fetched cache against a second cache, the resume log and the COG.

Usage: python verify_cache.py CACHE_A CACHE_B PRESENT_AFTER_KILL RESUME_STDERR [N_SAMPLE]

Independent of tin_engine: the header is parsed with tifffile, the blocks are
named blocks/<i // across>/<i % across>.bin as the design says, and the COG is
read with curl, not with tin_engine's RangeClient. Prints one line per check
and exits 1 if any fails.
"""

import hashlib
import io
import json
import math
import random
import re
import subprocess
import sys
from pathlib import Path

import tifffile

OBJ = "anadem-v1/anadem_v1_compressed_COG"
cache_a, cache_b, present_file, resume_err = (Path(p) for p in sys.argv[1:5])
n_sample = int(sys.argv[5]) if len(sys.argv) > 5 else 12
failed = False


def check(ok: bool, what: str) -> None:
    global failed
    failed |= not ok
    print(("PASS " if ok else "FAIL ") + what)


page = tifffile.TiffFile(io.BytesIO((cache_a / OBJ / "header.bin").read_bytes())).pages[0]
across = math.ceil(page.imagewidth / page.tilewidth)
offsets, counts = page.dataoffsets, page.databytecounts
manifest = json.loads((cache_a / "anadem-v1" / "manifest.json").read_text())
url = manifest["objects"]["anadem_v1_compressed_COG"]["url"]


def files(root: Path) -> dict[int, Path]:
    out = {}
    for p in (root / OBJ / "blocks").glob("*/*.bin"):
        out[int(p.parent.name) * across + int(p.stem)] = p
    return out


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


fa, fb = files(cache_a), files(cache_b)
print(f"page {page.imagewidth} x {page.imagelength}, tile {page.tilewidth}, across {across}")
check(set(fa) == set(fb), f"same block set: A {len(fa)}, B {len(fb)}")
same = [i for i in fa if i in fb and sha(fa[i]) == sha(fb[i])]
check(len(same) == len(fa), f"sha256 equal block by block: {len(same)} of {len(fa)}")
sizes = [fa[i].stat().st_size == counts[i] for i in fa]
check(all(sizes), f"every block's size is the header's byte count: {sum(sizes)} of {len(fa)}")
total = sum(p.stat().st_size for p in fa.values())
print(f"block bytes A {total:,}, B {sum(p.stat().st_size for p in fb.values()):,}")
check(sha(cache_a / OBJ / "header.bin") == sha(cache_b / OBJ / "header.bin"), "header.bin equal")

# Resume: what the resumed run asked for is exactly what the killed run lacked.
present = set()
for line in present_file.read_text().split():
    row, col = line.split("/")[-2:]
    present.add(int(row) * across + int(col.removesuffix(".bin")))
asked = set()
for m in re.finditer(r"RANGE bytes=(\d+)-(\d+)", resume_err.read_text()):
    lo, hi = int(m[1]), int(m[2]) + 1
    asked |= {i for i in fa if offsets[i] >= lo and offsets[i] + counts[i] <= hi}
check(not (asked & present), f"resume asked no present block: overlap {len(asked & present)}")
check(asked | present == set(fa), f"present {len(present)} + asked {len(asked)} = {len(fa)}")

# The COG itself, read by curl.
rng = random.Random(20261002)
for i in sorted(rng.sample(sorted(fa), min(n_sample, len(fa)))):
    lo, hi = offsets[i], offsets[i] + counts[i] - 1
    body = subprocess.run(
        ["curl", "-sf", "-r", f"{lo}-{hi}", url], check=True, capture_output=True
    ).stdout
    ok = hashlib.sha256(body).hexdigest() == sha(fa[i])
    check(ok, f"block {i} (row {i // across}, col {i % across}) bytes {lo}-{hi} equal to the COG")
sys.exit(1 if failed else 0)
