"""Medians of the land cover row and the total from raw/stats/*_stats.md, plus the vtk identity check."""
import re
import statistics
from pathlib import Path

R = Path(__file__).resolve().parent.parent / "raw" / "stats"
ROW = re.compile(r"^\| \**(land cover|total)\** \| \**([0-9.]+)\**")


def rows(path: Path) -> dict[str, float]:
    out = {}
    for line in path.read_text().splitlines():
        m = ROW.match(line)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


sha = dict(line.split() for line in (R / "vtk_sha256.txt").read_text().splitlines())
for c in ("numedalslagen", "skiensvassdraget"):
    med = {}
    for side in ("base", "branch"):
        vals = [rows(R / f"{side}_{c}_r{r}_stats.md") for r in (1, 2, 3)]
        for k in ("land cover", "total"):
            xs = [v[k] for v in vals]
            med[side, k] = statistics.median(xs)
            print(f"{c} {side} {k}: {xs} median {med[side, k]:.3f}")
    lc = med["branch", "land cover"] / med["base", "land cover"]
    print(f"{c} land cover branch/base {lc:.3f} (speed-up {1 / lc:.2f}x) gate<=1/3: {lc <= 1 / 3}")
    print(f"{c} total branch/base {med['branch', 'total'] / med['base', 'total']:.3f}")
    hashes = {sha[f"{s}_{c}_r{r}"] for s in ("base", "branch") for r in (1, 2, 3)}
    print(f"{c} vtk sha256 distinct over 6 runs: {len(hashes)} {sorted(hashes)}")
