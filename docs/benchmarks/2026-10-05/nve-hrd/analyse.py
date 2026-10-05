"""The acceptance run's tables (@perf), from the committed files only.

Usage: analyse.py HERE  (HERE holds results.csv, summary.json, batch.log,
rss.tsv.gz and, once checks.py has run, checks.csv); writes Markdown on stdout.

Per-station peak memory: rss.tsv.gz samples the batch process's resident
set every 0.2 s; batch.log time-stamps each station's line as it is written,
after the station ends. A station's peak is the largest sample after the
previous station's line and up to its own: a sampled figure (a peak shorter
than 0.2 s can be missed), and memory the allocator keeps from an earlier
station is counted again.
"""

from __future__ import annotations

import csv
import gzip
import json
import statistics
import sys
from pathlib import Path

BANDS = ("under 10", "10-100", "100-1000", "over 1000")
CLASSES = ("match", "close", "miss", "uncertain", "refused")


def band(area: float | None) -> str:
    if area is None:
        return "no polygon"
    return BANDS[sum(area >= b for b in (10.0, 100.0, 1000.0))]


def fnum(v: str, fmt: str = "{:.1f}", scale: float = 1.0) -> str:
    return "" if v in ("", None) else fmt.format(float(v) * scale)


def main(here: Path) -> None:
    rows = list(csv.DictReader(open(here / "results.csv", encoding="utf-8")))
    summary = json.loads((here / "summary.json").read_text())
    log = []
    for line in (here / "batch.log").read_text(encoding="utf-8").splitlines():
        t, _, text = line.partition("\t")
        log.append((float(t), text))
    rss = [
        (float(a), int(b))
        for a, b in (ln.split("\t") for ln in gzip.open(here / "rss.tsv.gz", "rt") if ln.strip())
    ]
    ids = {r["station"] for r in rows}
    ends = [(t, text.split(" ", 1)[0]) for t, text in log if text.split(" ", 1)[0] in ids]
    peak: dict[str, float] = {}
    start = log[0][0]
    for t, sid in ends:
        inside = [kb for ts, kb in rss if start < ts <= t + 0.25]
        peak[sid] = max(inside) / 1024 / 1024 if inside else float("nan")
        start = t
    p = print
    p("## Classes\n")
    p("| class | all | river-seeded | lake-seeded | no seed (refused before) |")
    p("|---|---:|---:|---:|---:|")
    for c in CLASSES:
        sel = [r for r in rows if r["class"] == c]
        n = lambda k: sum(r["seeded_by"] == k for r in sel)  # noqa: E731
        p(f"| {c} | {len(sel)} | {n('river')} | {n('lake')} | {n('')} |")
    p(f"| total | {len(rows)} | | | |\n")
    p(f"`match_by`: {summary['match_by']}\n")
    p("## By size band (NVE's polygon area, km²)\n")
    p("| band | stations | match | close | miss | uncertain | refused | share uncertain |")
    p("|---|---:|---:|---:|---:|---:|---:|---:|")
    for b in (*BANDS, "no polygon"):
        sel = [r for r in rows if band(float(r["reference_area_km2"]) if r["reference_area_km2"] else None) == b]
        if not sel:
            continue
        cnt = [sum(r["class"] == c for r in sel) for c in CLASSES]
        assessed = len(sel) - cnt[4]
        share = f"{100 * cnt[3] / assessed:.0f} %" if assessed else ""
        p(f"| {b} | {len(sel)} | " + " | ".join(map(str, cnt)) + f" | {share} |")
    p("\n## By tile count (tiles our fine outline meets)\n")
    p("| tiles | stations | match | close | miss | uncertain | share uncertain |")
    p("|---|---:|---:|---:|---:|---:|---:|")
    for g, grp in summary["by_tiles"].items():
        c = grp["classes"]
        share = "" if grp["uncertain_share"] is None else f"{100 * grp['uncertain_share']:.0f} %"
        p(f"| {g} | {grp['stations']} | {c['match']} | {c['close']} | {c['miss']} | {c['uncertain']} | {share} |")
    p("\n## Scored stations (match, close, miss): percentiles\n")
    p("| measure | min | p10 | p25 | p50 | p75 | p90 | max |")
    p("|---|---:|---:|---:|---:|---:|---:|---:|")
    for k, m in summary["scored"].items():
        if m is None:
            continue
        p(f"| {k} | " + " | ".join(f"{v:.3f}" for v in m.values()) + " |")
    p(f"\nUncertain causes (a station can have several): {summary['uncertain_causes']}\n")
    p(f"Refusal causes: {summary['refusal_causes']}\n")
    p("## Refusals\n")
    p("| station | name | cause | NVE km² | message |")
    p("|---|---|---|---:|---|")
    for r in rows:
        if r["class"] == "refused":
            msg = r["refusal_message"].replace("|", "/")
            if r["grid_tiles"]:
                msg += f" (tiles: {r['grid_tiles'].replace(';', ', ')})"
            p(f"| {r['station']} | {r['name']} | {r['refusal_cause']} | {fnum(r['reference_area_km2'])} | {msg} |")
    p("\n## Every station\n")
    p("| station | name | class | by | seed | NVE km² | ours km² | NVE's in ours % | ours in NVE's % | ratio | offset m | causes | tiles | s | peak GB |")
    p("|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|")
    for r in rows:
        p(
            f"| {r['station']} | {r['name']} | {r['class']} | {r['match_by']} | {r['seeded_by']} | "
            f"{fnum(r['reference_area_km2'])} | {fnum(r['fine_area_km2'], '{:.2f}')} | "
            f"{fnum(r['nve_in_ours'], scale=100)} | {fnum(r['ours_in_nve'], scale=100)} | "
            f"{fnum(r['area_ratio'], '{:.3f}')} | {fnum(r['divide_offset_m'])} | "
            f"{r['causes'].replace(';', ', ')} | {r['tiles']} | {fnum(r['seconds'])} | "
            f"{peak.get(r['station'], float('nan')):.2f} |"
        )
    secs = [float(r["seconds"]) for r in rows]
    p("\n## Time and memory\n")
    p(f"Stations' own seconds: total {sum(secs):.0f} s, median {statistics.median(secs):.1f} s, "
      f"max {max(secs):.0f} s ({rows[secs.index(max(secs))]['station']}).")
    big = max(peak, key=lambda k: peak[k] if peak[k] == peak[k] else -1)
    p(f"Sampled peak RSS per station: median {statistics.median(v for v in peak.values() if v == v):.2f} GB, "
      f"max {peak[big]:.2f} GB ({big}).")
    tail = [text for _, text in log if text.strip().startswith(("real", "maximum resident"))
            or "maximum resident set size" in text or text.strip().endswith("real")]
    p("Whole process (`/usr/bin/time -l`): " + "; ".join(t.strip() for t in tail))
    chk = here / "checks.csv"
    if chk.exists():
        c = list(csv.DictReader(open(chk)))
        river = [r for r in c if r["seeded_by"] == "river"]
        unc = {r["station"] for r in rows if r["class"] == "uncertain"}
        p("\n## Step 6 checks and increment 22's outline guarantees\n")
        p(f"Stations checked: {len(c)} not refused ({len(river)} river-seeded, {len(c) - len(river)} lake-seeded).\n")
        bad = lambda key, ok: [r["station"] for r in c if r[key] != "" and r[key] != ok]  # noqa: E731
        orc = bad("oracle", "True")
        mono = bad("monotone", "True")
        drn = bad("drains", "True")
        off = [r["station"] for r in river if int(r["lowered_off_chain"]) != 0]
        simple = bad("simple", "True")
        inside = bad("seed_inside", "True")
        area = [r["station"] for r in c if abs(float(r["area_diff_m2"])) > 1e-9 * float(r["fine_area_m2"])]
        held = sorted(set(mono) | set(drn) | set(off))
        p("| check | failures | stations |")
        p("|---|---:|---|")
        p(f"| placed node's count = catchment nodes (oracle; a failure is a defect) | {len(orc)} | {', '.join(orc)} |")
        p(f"| counts strictly rising downstream along the chain, end closed (`monotone`) | {len(mono)} | {', '.join(mono)} |")
        p(f"| each chain node drains into the next (`drains`) | {len(drn)} | {', '.join(drn)} |")
        p(f"| burn lowered no node off the chain (first window) | {len(off)} | {', '.join(off)} |")
        p(f"| reduced outline simple | {len(simple)} | {', '.join(simple)} |")
        p(f"| seed strictly inside the reduced outline | {len(inside)} | {', '.join(inside)} |")
        worst = max(abs(float(r["area_diff_m2"])) for r in c)
        p(f"| area kept (reduced − fine) | {len(area)} over 1e-9 of the area | largest difference {worst:.3g} m² |")
        p(f"\nStations whose burn did not hold (monotone or drains failed, or a node off the chain lowered): "
          f"{len(held)}; of the {len(unc)} `uncertain`, they explain {len(set(held) & unc)}.\n")
        lk = [r["station"] for r in river if r["lake_above_p"] == "True"]
        p(f"River rows with a lake on the reach above P: {len(lk)}\n")
        p("| station | name | class | NVE's in ours % | ours in NVE's % | causes |")
        p("|---|---|---|---:|---:|---|")
        by = {r["station"]: r for r in rows}
        for s in lk:
            r = by[s]
            p(f"| {s} | {r['name']} | {r['class']} | {fnum(r['nve_in_ours'], scale=100)} | "
              f"{fnum(r['ours_in_nve'], scale=100)} | {r['causes'].replace(';', ', ')} |")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
