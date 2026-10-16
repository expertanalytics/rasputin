"""Increment 36 acceptance (b): the GeoPackage writer's time and memory on the
built command, `docs/increments/36-geopackage-output.md` section 9.

For each case, `rasputin mesh ... --out X.gpkg --stats -` and the same command
with `--out X.vtk --binary` (the subtraction baseline), each under
`/usr/bin/time -l`, one warm-up then `RUNS` timed runs. While a command runs,
the rasputin child's RSS is sampled with `ps` every 0.1 s, so the writer's own
share is the RSS growth inside the `write: geopackage` window (the last
`write: geopackage` seconds of the process; it is the last phase of the run),
not the whole command's peak, which the meshing sets. Medians are reported.

Usage, from the worktree root:
    RASPUTIN_DATA=... .venv/bin/python docs/benchmarks/2026-10-11/gpkg-out/scripts/writer_time.py \
        <out dir> [case ...]
Writes <out dir>/writer-<case>-<fmt>-<n>.log (the command's stderr+stdout and the
/usr/bin/time block), <out dir>/writer-<case>-<fmt>-<n>.rss (t_s rss_kB samples),
<out dir>/writer-results.json and prints the table.
"""

from __future__ import annotations

import json
import os
import re
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path

RUNS = 3
TREE = Path(__file__).resolve().parents[5]  # scripts/gpkg-out/2026-10-11/benchmarks/docs/<tree>
DATA = os.environ.get("RASPUTIN_DATA", str(TREE.parent / "rasputin_data"))
RASPUTIN = TREE / ".venv" / "bin" / "rasputin"
CASES = {
    "numedalslagen": [
        "mesh", "--dem", f"{DATA}/DTM10_UTM33_20260925",
        "--domain", f"{DATA}/numedalslagen_outline_nve.geojson",
        "--features", f"{DATA}/corine2018_dtm10_utm33.gpkg", "--features-layer", "corine2018",
        "--features-map", "corine", "--tolerance", "10",
    ],
    "romsdal-slope": [
        "mesh", "--dem", f"{DATA}/DTM10_UTM33_20260925/6901_3_10m_z33.tif",
        "--tolerance", "10", "--tolerance-slope", "2", "25", "35",
    ],
}  # fmt: skip
FORMATS = {"gpkg": ["--out", "mesh.gpkg"], "vtk": ["--out", "mesh.vtk", "--binary"]}


def sample(proc: subprocess.Popen[bytes], out: list[tuple[float, int]], t0: float) -> None:
    child: str | None = None
    while proc.poll() is None:
        if child is None:
            kids = subprocess.run(["pgrep", "-P", str(proc.pid)], capture_output=True, text=True)
            child = (kids.stdout.split() or [None])[0]
        if child is not None:
            ps = subprocess.run(["ps", "-o", "rss=", "-p", child], capture_output=True, text=True)
            if ps.stdout.strip():
                out.append((time.perf_counter() - t0, int(ps.stdout.split()[0])))
        time.sleep(0.1)


def one(case: str, fmt: str, n: int, out_dir: Path, work: Path) -> dict[str, float]:
    target = work / FORMATS[fmt][1]
    argv = ["/usr/bin/time", "-l", str(RASPUTIN), *CASES[case], "--out", str(target),
            *FORMATS[fmt][2:], "--stats", "-"]  # fmt: skip
    for stale in work.glob("mesh.*"):
        stale.unlink()
    samples: list[tuple[float, int]] = []
    t0 = time.perf_counter()
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    thread = threading.Thread(target=sample, args=(proc, samples, t0))
    thread.start()
    stdout, stderr = proc.communicate()
    t_end = time.perf_counter() - t0
    thread.join()
    text = stdout.decode() + "\n--- stderr ---\n" + stderr.decode()
    (out_dir / f"writer-{case}-{fmt}-{n}.log").write_text(text)
    (out_dir / f"writer-{case}-{fmt}-{n}.rss").write_text(
        "".join(f"{t:.2f} {kb}\n" for t, kb in samples)
    )
    if proc.returncode != 0:
        raise SystemExit(f"{case} {fmt} run {n}: exited {proc.returncode}\n{text[-3000:]}")
    real = float(re.search(r"([\d.]+) real", text).group(1))
    peak = int(re.search(r"(\d+)\s+maximum resident set size", text).group(1))
    # the --stats Markdown table: `| write: geopackage | 3.980 | 40.0 % |`
    phases = {m.group(1): float(m.group(2))
              for m in re.finditer(r"\|\s*(write: [a-z]+(?: [a-z]+)*)\s*\|\s*([\d.]+)\s*\|", text)}  # fmt: skip
    write_s = phases.get("write: geopackage", phases.get("write: encode", 0.0) + phases.get("write: disk", 0.0))
    window = [kb for t, kb in samples if t >= t_end - write_s - 0.05]
    result = {
        "real_s": real, "peak_rss_bytes": peak, "write_s": write_s,
        "size_bytes": target.stat().st_size,
        "writer_rss_start_bytes": (window[0] if window else 0) * 1024,
        "writer_rss_peak_bytes": (max(window) if window else 0) * 1024,
        "samples": len(samples),
    }  # fmt: skip
    result["writer_share_bytes"] = result["writer_rss_peak_bytes"] - result["writer_rss_start_bytes"]
    for phase, seconds in phases.items():
        result[phase] = seconds
    return result


def main() -> None:
    out_dir = Path(sys.argv[1]).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cases = sys.argv[2:] or list(CASES)
    work = Path(os.environ.get("TMPDIR", "/tmp")) / "rasputin-writer-time"
    work.mkdir(exist_ok=True)
    results: dict[str, dict[str, object]] = {}
    for case in cases:
        for fmt in FORMATS:
            one(case, fmt, 0, out_dir, work)  # warm-up
            runs = [one(case, fmt, n, out_dir, work) for n in range(1, RUNS + 1)]
            med = {k: statistics.median(r[k] for r in runs) for k in runs[0]}
            results[f"{case}/{fmt}"] = {"runs": runs, "median": med}
            print(f"{case} {fmt}: real {med['real_s']:.2f} s, write {med['write_s']:.3f} s, "
                  f"peak RSS {med['peak_rss_bytes'] / 1e9:.2f} GB, writer window "
                  f"{med['writer_rss_start_bytes'] / 1e9:.2f} -> {med['writer_rss_peak_bytes'] / 1e9:.2f} GB "
                  f"(share {med['writer_share_bytes'] / 1e6:.0f} MB), file {med['size_bytes'] / 1e6:.1f} MB",
                  flush=True)  # fmt: skip
    (out_dir / "writer-results.json").write_text(json.dumps(results, indent=1) + "\n")


if __name__ == "__main__":
    main()
