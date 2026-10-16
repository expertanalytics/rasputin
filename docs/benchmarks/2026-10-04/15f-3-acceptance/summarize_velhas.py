"""Tables for 15f-3-acceptance.md part (3), Velhas (@perf). python summarize_velhas.py > tables-velhas.md
Reads velhas/<label>/results.json (velhas.py run / check) and velhas-check/*.json (strip_check.py).
Timings: 4 timed --binary runs a side (2 per run, runs base-r1 15f3-r1 15f3-r2 base-r2)."""
import json
import re
import statistics as st
from pathlib import Path

HERE = Path(__file__).resolve().parent
V = HERE / "velhas"
SIDES = {"base": ["base-r1", "base-r2"], "15f3": ["15f3-r1", "15f3-r2"]}
ROWS = ["max_error_m", "resampled_grid_max_error_m", "dem_check_points_inserted", "dem_check_rounds",
        "line_points_checked", "line_points_on_nodata", "line_points_inserted", "line_points_refused",
        "line_points_refused_max_error_m", "line_max_error_m", "line_points_duplicate", "refinement_rounds",
        "points_inserted"]
PHASES = ["decode", "resample", "refine", "edge strip: generate", "final check: scan (parallel)",
          "final check: split + flip (serial)", "write: encode", "other"]


def phases(p: Path) -> dict[str, float]:
    t = p.read_text()
    out = {m[1]: float(m[2]) for m in re.finditer(r"^\| ([^|*]+?) \| ([\d.]+) \| [\d.]+ % \|$", t, re.M)}
    tot = re.search(r"^\| \*\*total\*\* \| \*\*([\d.]+)\*\*", t, re.M)
    out["total"] = float(tot[1]) if tot else float("nan")
    return out


def named(p: Path) -> dict[str, str]:
    t = p.read_text()
    return {m[2]: m[1].strip() for m in re.finditer(r"^\|[^|]*\| ([^|]*) \| `([a-z_0-9]+)` \|$", t, re.M)}


res = {k: json.loads((V / k / "results.json").read_text()) for v in SIDES.values() for k in v}
for i, tol in enumerate((20.0, 10.0, 5.0)):
    print(f"### Velhas, EPSG:31983, tolerance {tol:g} m\n")
    print("| item | base c193cb1 | 15f-3 83c7fd2 |\n|---|---:|---:|")
    col = {}
    for side, ks in SIDES.items():
        recs = [next(r for r in res[k]["runs"] if r["tolerance"] == tol) for k in ks]
        timed = [t for r in recs for t in r["timed"]]
        logs = [V / k / "logs" for k in ks]
        ph = [phases(lg / f"t{tol:g}_run{j}.stats.md") for lg in logs for j in (1, 2)]
        col[side] = dict(
            real=st.median(t["real_s"] for t in timed), rr=f"{min(t['real_s'] for t in timed):.2f}-{max(t['real_s'] for t in timed):.2f}",
            rss=max(t["max_rss_bytes"] for t in timed) / 2**30, named=named(logs[0] / f"t{tol:g}_ascii.stats.md"),
            ph={p: st.median(x.get(p, 0.0) for x in ph) for p in PHASES + ["total"]},
            q=recs[0]["quality"], sha=recs[0]["mesh_sha256"], check=recs[0].get("check"), control=recs[0].get("control"),
            power={r["pmset_before"].splitlines()[0] for r in recs} | {r["pmset_after"].splitlines()[0] for r in recs})
    b, n = col["base"], col["15f3"]
    print(f"| process wall, median of 4 (s) | {b['real']:.2f} ({b['rr']}) | {n['real']:.2f} ({n['rr']}) |")
    print(f"| peak RSS, max of 4 (GiB) | {b['rss']:.2f} | {n['rss']:.2f} |")
    for k in ROWS:
        print(f"| `{k}` | {b['named'].get(k, '-')} | {n['named'].get(k, '-')} |")
    for k in PHASES + ["total"]:
        print(f"| phase: {k} (s, median) | {b['ph'][k]:.3f} | {n['ph'][k]:.3f} |")
    print(f"| quality: max_error_m <= tolerance | {float(b['named']['max_error_m']) <= tol} | "
          f"{float(n['named']['max_error_m']) <= tol} |")  # velhas.py passed NaN to quality(): its own flag is not used
    for k in ("worst_angle", "max_degree", "share_under_1", "delaunay_violations", "delaunay_checked"):
        print(f"| quality: {k} | {b['q'][k]} | {n['q'][k]} |")
    for part in ("interior", "strip"):
        f = lambda c: f"{c[part]['over_tolerance']} of {c[part]['nodes']}" if c else "-"
        print(f"| 15c source-node check, {part}: over | {f(b['check'])} | {f(n['check'])} |")
    if n["control"]:
        print(f"| control (mesh +30 m east), interior over | - | {n['control']['interior']['over_tolerance']} |")
    cb = json.loads((HERE / "velhas-check" / f"base-r1_t{tol:g}.json").read_text())
    cn = json.loads((HERE / "velhas-check" / f"15f3-r1_t{tol:g}.json").read_text())
    for part in ("crossings", "midpoints"):
        print(f"| strip check, {part}: over of checked (max m) | {cb[part]['over']} of {cb[part]['checked']} "
              f"({cb[part]['max_error_m']:.4f}) | {cn[part]['over']} of {cn[part]['checked']} ({cn[part]['max_error_m']:.4f}) |")
    print(f"| power (pmset first line) | {'; '.join(sorted(b['power']))} | {'; '.join(sorted(n['power']))} |")
    print(f"| mesh sha256 | `{b['sha'][:16]}` | `{n['sha'][:16]}` |")
    print()
