#!/usr/bin/env python3
"""Mesh every SMHI SVAR 2022 sub-catchment of four Swedish river systems with `rasputin mesh` (increment 39).

    se2_batch.py prepare
    se2_batch.py run --level baro|fine|huvud [--jobs N] [--only SYSTEM] [--retry] [--timeout S]
    se2_batch.py report [--copy-to DIR]

Inputs (never written): `$RASPUTIN_DATA/sweden_smhi_svar/se2/*_delaro_svar2022_3006.geojson`
(one file per main catchment, EPSG:3006), the GLO-30 cache and the CORINE
GeoPackage the Lagan quick-check case uses. Everything else lives under
`$RASPUTIN_SCRATCH/sweden/se2/` (default `<main checkout>/../rasputin_scratch`):

    units/<system>/<ARO_UUID>.geojson    baro/<system>/<BARONR>.geojson    huvud/<system>/<HARO>.geojson
    index_<level>.csv                    one row per polygon: system, id, name, area_km2, path, note
    out/<system>/<level>/<id>.gpkg (+ .vtk for baro and huvud; _stats.md, _record.json beside)
    summary_<level>.csv                  the ledger, one row appended as each subprocess exits
    report.md, overview_<system>.svg, summary_<system>.geojson

`run` meshes the polygons of a level not yet in its ledger, largest SVAR area
first, `--jobs` subprocesses at once; a killed batch restarts with the same
command and skips what is done. A non-zero exit of `mesh` is a `refused` row
(its message in `mesh`'s own words) and the batch goes on; a crash of this
script stops it. `--retry` drops the `refused`/`error` rows first and runs them
again. Standard library plus shapely; no plotting library (the maps are SVG
from shapely's own `svg()`).
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from shapely.geometry import mapping, shape
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parent
TOLERANCE_M = 10.0
CORINE = "corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg"
CORINE_LAYER = "U2018_CLC2018_V2020_20u1"
SYSTEMS = {"umealven": 28000, "angermanalven": 38000, "indalsalven": 40000, "ljungan": 42000}
LEVELS = {"fine": ("gpkg",), "baro": ("gpkg", "vtk"), "huvud": ("gpkg", "vtk")}
INDEX_FIELDS = ["system", "level", "id", "name", "area_km2", "units", "path", "note"]
LEDGER_FIELDS = ["system", "level", "id", "name", "area_km2", "triangles", "vertices", "max_error_m",
                 "seconds", "peak_gb", "status", "message", "vtk_seconds"]  # fmt: skip
CRS_MEMBER = {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::3006"}}
BOX_CHARS = re.compile(r"[─-╿]")
NOISE = re.compile(r"^(<frozen|Usage:|Try ')")
STATUS_FILL = {"ok": "#4c9a6a", "refused": "#d1495b", "error": "#d1495b", "not run": "#d9d9d9"}


def main_checkout() -> Path:
    done = subprocess.run(["git", "-C", str(HERE), "rev-parse", "--path-format=absolute", "--git-common-dir"],
                          capture_output=True, text=True, check=False)  # fmt: skip
    return Path(done.stdout.strip()).parent if done.returncode == 0 else HERE.parents[3]


def env_path(var: str, default: Path) -> Path:
    return Path(os.environ[var]).expanduser() if os.environ.get(var) else default


class Paths:
    def __init__(self, args: argparse.Namespace) -> None:
        main = main_checkout()
        self.data = env_path("RASPUTIN_DATA", main.parent / "rasputin_data")
        self.root = env_path("RASPUTIN_SCRATCH", main.parent / "rasputin_scratch") / "sweden" / "se2"
        self.svar = self.data / "sweden_smhi_svar" / "se2"
        self.cache = self.data / "sweden_glo30_cache"
        self.corine = self.data / CORINE
        self.rasputin = Path(args.rasputin) if args.rasputin else HERE.parents[3] / ".venv" / "bin" / "rasputin"
        for what in (self.svar, self.cache, self.corine, self.rasputin):
            if not what.exists():
                sys.exit(f"refused: {what} is missing")
        self.root.mkdir(parents=True, exist_ok=True)

    def index(self, level: str) -> Path:
        return self.root / f"index_{level}.csv"

    def ledger(self, level: str) -> Path:
        return self.root / f"summary_{level}.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    tmp = path.with_suffix(".part")
    with tmp.open("w", newline="") as f:
        w = csv.DictWriter(f, fields)
        w.writeheader()
        w.writerows(rows)
    tmp.replace(path)


def write_polygon(path: Path, geom: object, props: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    feature = {"type": "Feature", "properties": props, "geometry": mapping(geom)}
    path.write_text(json.dumps({"type": "FeatureCollection", "crs": CRS_MEMBER, "features": [feature]}))


def prepare(p: Paths, _args: argparse.Namespace) -> None:
    """One GeoJSON per polygon at each level, and one index per level; nothing is fetched."""
    index: dict[str, list[dict[str, object]]] = {level: [] for level in LEVELS}
    for system, haro in SYSTEMS.items():
        src = p.svar / f"{system}_{haro}_delaro_svar2022_3006.geojson"
        feats = json.loads(src.read_text())["features"]
        by_baro: dict[int, list[dict[str, object]]] = {}
        for f in feats:
            props, uid = f["properties"], f["properties"]["ARO_UUID"]
            path = p.root / "units" / system / f"{uid}.geojson"
            write_polygon(path, shape(f["geometry"]), props)
            index["fine"].append({"system": system, "level": "fine", "id": uid, "name": props.get("NAME") or "",
                                  "area_km2": f"{props['AREA'] / 1e6:.4f}", "units": 1, "path": path, "note": ""})  # fmt: skip
            by_baro.setdefault(int(props["BARONR"]), []).append(f)
        for baro, members in sorted(by_baro.items()):
            geom, note = unary_union([shape(f["geometry"]) for f in members]).buffer(0), ""
            if geom.geom_type != "Polygon":
                parts = sorted(geom.geoms, key=lambda g: g.area, reverse=True)
                note = f"dissolve gave {len(parts)} parts; kept the largest ({parts[0].area / 1e6:.2f} km2)"
                geom = parts[0]
            area = sum(f["properties"]["AREA"] for f in members) / 1e6
            path = p.root / "baro" / system / f"{baro}.geojson"
            write_polygon(path, geom, {"BARONR": baro, "HARO": haro, "system": system, "units": len(members),
                                       "AREA_km2": round(area, 4), "note": note})  # fmt: skip
            index["baro"].append({"system": system, "level": "baro", "id": str(baro), "name": f"BARO {baro}",
                                  "area_km2": f"{area:.4f}", "units": len(members), "path": path, "note": note})  # fmt: skip
        whole = json.loads((p.svar / f"huvud_{system}_{haro}_svar2022_3006.geojson").read_text())["features"][0]
        path = p.root / "huvud" / system / f"{haro}.geojson"
        write_polygon(path, shape(whole["geometry"]), whole["properties"])
        index["huvud"].append({"system": system, "level": "huvud", "id": str(haro), "name": whole["properties"]["NAME"],
                               "area_km2": f"{whole['properties']['AREA']:.4f}", "units": len(feats), "path": path, "note": ""})  # fmt: skip
        print(f"{system}: {len(feats)} units, {len(by_baro)} BARO areas", flush=True)
    for level, rows in index.items():
        write_csv(p.index(level), INDEX_FIELDS, rows)


def message_of(stderr: str) -> str:
    """`mesh`'s refusal in its own words: the Rich box and the interpreter's warnings stripped."""
    lines = [BOX_CHARS.sub(" ", ln).strip() for ln in stderr.splitlines()]
    text = " ".join(ln for ln in lines if ln and not NOISE.match(ln))
    text = re.sub(r"\s+", " ", text).removeprefix("Error ").strip()
    return text[:300]


def mesh_once(p: Paths, domain: Path, out: Path, timeout: float) -> tuple[subprocess.CompletedProcess[str], float, float]:
    """One `rasputin mesh`, as the Lagan quick-check case runs it; returns (process, wall s, peak GB)."""
    stem = out.with_suffix("") if out.suffix == ".gpkg" else out.with_name(f"{out.stem}_{out.suffix[1:]}")
    cmd = [str(p.rasputin), "mesh", "--dem", "glo30", "--cache", str(p.cache), "--out-crs", "EPSG:3006",
           "--domain", str(domain), "--features", str(p.corine), "--features-layer", CORINE_LAYER,
           "--features-crs", "EPSG:3035", "--features-map", "corine", "--tolerance", str(TOLERANCE_M),
           "--out", str(out), "--stats", f"{stem}_stats.md", "--record", f"{stem}_record.json"]  # fmt: skip
    t0 = time.perf_counter()
    with (stem.parent / f"{stem.name}.stderr").open("w+") as err:
        proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=err, text=True)
        watchdog = threading.Timer(timeout, proc.kill)  # a hung unit is one refused row, not a hung night
        watchdog.start()
        try:
            _, status, usage = os.wait4(proc.pid, 0)  # wait4: this child's own peak memory
            proc.returncode = os.waitstatus_to_exitcode(status)
        except BaseException:
            proc.kill()
            raise
        finally:
            watchdog.cancel()
        err.seek(0)
        stderr = err.read()
    if proc.returncode < 0:
        stderr += f"\nkilled by signal {-proc.returncode} after {time.perf_counter() - t0:.0f} s (timeout {timeout:g} s)"
    peak = usage.ru_maxrss / (1e9 if sys.platform == "darwin" else 1e6)
    return subprocess.CompletedProcess(cmd, proc.returncode, "", stderr), time.perf_counter() - t0, peak


def run_one(p: Paths, row: dict[str, str], timeout: float) -> dict[str, object]:
    out_dir = p.root / "out" / row["system"] / row["level"]
    out_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, object] = {k: row[k] for k in ("system", "level", "id", "name", "area_km2")}
    result.update(triangles="", vertices="", max_error_m="", seconds="", peak_gb="", vtk_seconds="",
                  status="error", message=row["note"])  # fmt: skip
    try:
        out = out_dir / f"{row['id']}.gpkg"
        proc, seconds, peak = mesh_once(p, Path(row["path"]), out, timeout)
        result.update(seconds=f"{seconds:.2f}", peak_gb=f"{peak:.2f}")
        if proc.returncode != 0 or not out.exists():
            result.update(status="refused", message=message_of(proc.stderr) or f"exit {proc.returncode}")
            return result
        record = json.loads(out.with_name(f"{row['id']}_record.json").read_text())
        stats = out.with_name(f"{row['id']}_stats.md").read_text()
        sizes = dict(re.findall(r"^\| (output triangles|output vertices) \| (\d+) \|$", stats, re.M))
        result.update(triangles=sizes.get("output triangles", ""), vertices=sizes.get("output vertices", ""),
                      max_error_m=record.get("max_error_m", ""), status="ok")  # fmt: skip
        if float(record.get("max_error_m", 0)) > TOLERANCE_M:
            result["message"] = f"{result['message']} over tolerance".strip()
        if "vtk" in LEVELS[row["level"]]:
            vtk, vtk_s, _ = mesh_once(p, Path(row["path"]), out_dir / f"{row['id']}.vtk", timeout)
            result["vtk_seconds"] = f"{vtk_s:.2f}"
            if vtk.returncode != 0:
                result["message"] = f"{result['message']} vtk refused: {message_of(vtk.stderr)}".strip()
    except Exception as exc:  # noqa: BLE001  (an error row; the batch goes on)
        result.update(status="error", message=f"{type(exc).__name__}: {exc}"[:300])
    return result


def run(p: Paths, args: argparse.Namespace) -> None:
    index = read_csv(p.index(args.level))
    ledger = p.ledger(args.level)
    done_rows = read_csv(ledger) if ledger.exists() else []
    if args.retry:
        done_rows = [r for r in done_rows if r["status"] == "ok"]
        write_csv(ledger, LEDGER_FIELDS, done_rows)
    done = {r["id"] for r in done_rows}
    todo = [r for r in index if r["id"] not in done and (not args.only or r["system"] == args.only)]
    todo.sort(key=lambda r: float(r["area_km2"]), reverse=True)
    print(f"{args.level}: {len(index)} polygons, {len(done)} done, {len(todo)} to run, {args.jobs} at once", flush=True)
    if not ledger.exists():
        write_csv(ledger, LEDGER_FIELDS, [])
    t0, n = time.perf_counter(), 0
    with ledger.open("a", newline="") as f, ThreadPoolExecutor(max_workers=args.jobs) as pool:
        writer = csv.DictWriter(f, LEDGER_FIELDS)
        futures = [pool.submit(run_one, p, row, args.timeout) for row in todo]
        for future in as_completed(futures):
            row = future.result()  # an exception here is the runner's own bug: the batch stops
            writer.writerow(row)
            f.flush()
            n += 1
            print(f"{n}/{len(todo)} {row['system']} {row['id']} {row['area_km2']} km2 {row['status']} "
                  f"{row['triangles'] or '-'} tri {row['seconds']} s {row['message']}"[:200], flush=True)  # fmt: skip
    print(f"{args.level}: {n} run in {time.perf_counter() - t0:.0f} s", flush=True)


def latest_rows(p: Paths, level: str) -> dict[str, dict[str, str]]:
    ledger = p.ledger(level)
    return {r["id"]: r for r in read_csv(ledger)} if ledger.exists() else {}


def svg_map(path: Path, polygons: list[tuple[object, str]], title: str) -> None:
    """Outlines simplified to 100 m, filled by status, drawn with shapely's own svg()."""
    geoms = [g for g, _ in polygons]
    minx, miny, maxx, maxy = unary_union(geoms).bounds
    w, h = maxx - minx, maxy - miny
    width = 800
    scale = width / w
    height = int(h * scale) + 40
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
             f'viewBox="0 0 {width} {height}" font-family="sans-serif" font-size="12">',
             f'<text x="8" y="14">{title}</text>',
             f'<g transform="translate(0,{height}) scale({scale},-{scale}) translate({-minx},{-miny})">']  # fmt: skip
    for geom, status in polygons:
        body = geom.simplify(100).svg(scale_factor=0, fill_color=STATUS_FILL[status], opacity=1.0)
        parts.append(re.sub(r'stroke-width="[^"]*"', f'stroke-width="{1 / scale:.1f}"', body))  # one pixel
    parts.append("</g>")
    x = 8
    for status, fill in STATUS_FILL.items():
        parts.append(f'<rect x="{x}" y="22" width="12" height="12" fill="{fill}"/>'
                     f'<text x="{x + 16}" y="32">{status}</text>')  # fmt: skip
        x += 90
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def report(p: Paths, args: argparse.Namespace) -> None:
    lines = ["# SE2 catchments: the batch's record", "",
             f"Tolerance {TOLERANCE_M:g} m, GLO-30, CORINE 2018. Källa: SMHI (SVAR 2022 polygons as outlines; "
             "BARO areas dissolved from the units, a modification).", ""]  # fmt: skip
    for level in LEVELS:
        if not p.index(level).exists():
            continue
        index, rows = read_csv(p.index(level)), latest_rows(p, level)
        lines += [f"## Level `{level}`", "",
                  "| system | polygons | ok | refused | error | not run | triangles | mesh s | largest km² | slowest s | over 10 m |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]  # fmt: skip
        refusals: list[str] = []
        for system in SYSTEMS:
            mine = [r for r in index if r["system"] == system]
            if not mine:
                continue
            got = [rows.get(r["id"]) for r in mine]
            count = {s: sum(1 for g in got if g and g["status"] == s) for s in ("ok", "refused", "error")}
            not_run = sum(1 for g in got if g is None)
            oks = [g for g in got if g and g["status"] == "ok"]
            tri = sum(int(g["triangles"] or 0) for g in oks)
            secs = sum(float(g["seconds"] or 0) + float(g["vtk_seconds"] or 0) for g in got if g)
            largest = max(mine, key=lambda r: float(r["area_km2"]))
            slowest = max(oks, key=lambda g: float(g["seconds"]), default=None)
            over = sum(1 for g in oks if float(g["max_error_m"] or 0) > TOLERANCE_M)
            lines.append(f"| {system} | {len(mine)} | {count['ok']} | {count['refused']} | {count['error']} | {not_run} | "
                         f"{tri:,} | {secs:,.0f} | {float(largest['area_km2']):,.1f} | "
                         f"{(slowest['seconds'] + ' (' + slowest['id'][:8] + ')') if slowest else '-'} | {over} |")  # fmt: skip
            refusals += [f"- {system} `{g['id']}` {g['name']} ({float(g['area_km2']):.2f} km²): {g['status']}, {g['message']}"
                         for g in got if g and g["status"] != "ok"]  # fmt: skip
            polygons = [(shape(json.loads(Path(r["path"]).read_text())["features"][0]["geometry"]),
                         rows[r["id"]]["status"] if r["id"] in rows else "not run") for r in mine]  # fmt: skip
            svg_map(p.root / f"overview_{system}_{level}.svg", polygons, f"{system}, level {level}: {len(mine)} polygons")
            features = [{"type": "Feature", "geometry": mapping(g), "properties": {**r, **(rows.get(r["id"]) or {})}}
                        for (g, _), r in zip(polygons, mine, strict=True)]  # fmt: skip
            (p.root / f"summary_{system}_{level}.geojson").write_text(
                json.dumps({"type": "FeatureCollection", "crs": CRS_MEMBER, "features": features}))
        lines += ["", f"Refusals and errors ({len(refusals)}):", *(refusals or ["- none"]), ""]
    (p.root / "report.md").write_text("\n".join(lines))
    print("\n".join(lines))
    if args.copy_to:
        dest = Path(args.copy_to)
        dest.mkdir(parents=True, exist_ok=True)
        for f in [*p.root.glob("overview_*.svg"), *p.root.glob("summary_*.csv"), p.root / "report.md"]:
            (dest / f.name).write_bytes(f.read_bytes())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rasputin", help="the rasputin executable (default: this worktree's .venv/bin/rasputin)")
    sub = ap.add_subparsers(dest="step", required=True)
    sub.add_parser("prepare")
    r = sub.add_parser("run")
    r.add_argument("--level", choices=LEVELS, required=True)
    r.add_argument("--jobs", type=int, default=1)
    r.add_argument("--only", choices=SYSTEMS)
    r.add_argument("--retry", action="store_true")
    r.add_argument("--timeout", type=float, default=3600.0, help="seconds per mesh run; longer is killed and a refused row")
    rep = sub.add_parser("report")
    rep.add_argument("--copy-to")
    args = ap.parse_args()
    {"prepare": prepare, "run": run, "report": report}[args.step](Paths(args), args)


if __name__ == "__main__":
    main()
