"""Placement figures, before and after (@perf, 2026-10-05; increment 29, PR 2).

`docs/increments/29-nve-reference-catchments.md`, "Placement figures, before
and after (PR 2)". Composes the package's public functions and adds no logic
to the placement: `read_stations`, `read_segments`, `read_references`,
`gauge.place`, `catchment.delineate`, and for the flow paths one window read
round the reach, `burn.burn_reach` on it and `accumulate` on the raw and the
burnt arrays. Matplotlib is used here only (Ola's ruling, 2026-10-04).

Usage, from the repository root, the venv's python with matplotlib on its path:

    python render.py survey  DATA DEM WORK   # place and burn every station, first window only
    python render.py choose  DATA DEM WORK   # pick the four cases, confirming by full runs
    python render.py draw    DATA DEM WORK OUT [STATION ...]   # the cases, then any extras

Placement is `catchment --rivers`'s (the nearest line, no watercourse number),
as the design's confirming run is that command; `--tiers` anywhere uses the
batch's tiered placement instead (and suffixes the PNGs with `_tiers`);
`--extras-only` draws only the stations named after OUT.

DATA holds stations.geojson, rivers.geojson and reference.geojson (as
`rasputin fetch-stations nve-hrd` writes them); DEM is a tile directory;
WORK takes survey.csv, cases.json and one runs/<station>.json per full run.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import shapely  # noqa: E402
from matplotlib.colors import LightSource, LogNorm  # noqa: E402
from shapely.geometry import LineString, Polygon  # noqa: E402

from tin_engine._core import accumulate  # noqa: E402
from tin_engine.burn import BurnRefusal, burn_reach  # noqa: E402
from tin_engine.catchment import WINDOW_MARGIN_M, CatchmentError, CatchmentRequest, delineate  # noqa: E402
from tin_engine.dem_input import repository_for  # noqa: E402
from tin_engine.gauge import Gauge, place  # noqa: E402
from tin_engine.io.rivers import read_segments  # noqa: E402
from tin_engine.io.station_set import read_references, read_stations  # noqa: E402
from tin_engine.mosaic import Bounds, MosaicError, assemble, plan_mosaic  # noqa: E402
from tin_engine.raster import to_core  # noqa: E402
from tin_engine.sensitivity import SWING_MAX, assess  # noqa: E402

STREAM_KM2 = 0.1  # flow paths drawn: nodes draining at least this
MAP_PAD_M = 300.0  # the map: the reach's bounds plus this


class Inputs:
    def __init__(self, data: Path, dem: Path, tiers: bool) -> None:
        self.tiers = tiers
        self.stations, crs = read_stations(data / "stations.geojson")
        self.segments, rcrs, _ = read_segments(data / "rivers.geojson")
        self.refs, fcrs = read_references(data / "reference.geojson")
        assert crs == rcrs == fcrs, (crs, rcrs, fcrs)
        self.crs = crs
        self.repo, _ = repository_for((dem,))
        self.foot = self.repo.footprints()
        self.by_id = {s.station: s for s in self.stations}

    def placement(self, sid: str):
        """`catchment --rivers`'s placement, the nearest line; with --tiers the
        batch's (own watercourse number, prefix, river name, any)."""
        s = self.by_id[sid]
        more = {"watercourse": s.watercourse, "river": s.river} if self.tiers else {}
        return place(Gauge(x=s.x, y=s.y, **more), self.segments)

    def window(self, reach):  # stage A's first window: the reach's bounds plus corridor and margin
        line, pad = np.asarray(reach.line), reach.corridor + WINDOW_MARGIN_M
        (x0, y0), (x1, y1) = line.min(axis=0) - pad, line.max(axis=0) + pad
        plan = plan_mosaic(self.foot, Bounds(x_min=x0, y_min=y0, x_max=x1, y_max=y1))
        return assemble(plan, self.repo.load).tile


def first_window(inp: Inputs, p):
    """The survey's view of one station: burn and assess on the first window."""
    tile = inp.window(p.reach)
    burnt, path = burn_reach(tile, p.reach)
    acc = accumulate(to_core(burnt))
    m, line = tile.meta, np.asarray(p.reach.line)
    down = float(np.hypot(*np.diff(line, axis=0).T).sum()) - p.reach.at
    s = assess(acc.count, acc.reach, acc.flow_to, path, m.delta_x * m.delta_y / 1e6,
               p.reach.uncertainty, down)  # fmt: skip
    return tile, burnt, path, s


def survey(inp: Inputs, work: Path) -> None:
    cols = ["station", "name", "nve_poly_km2", "placed_on", "distance_m", "objectid", "lake",
            "confluence_near", "lowered_nodes", "lowered_max_m", "node_offset_m", "swing",
            "direction_ok", "causes", "well_posed_first", "refusal"]  # fmt: skip
    with (work / "survey.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, cols)
        w.writeheader()
        for st in inp.stations:
            ref = inp.refs.get(st.station)
            row: dict = {"station": st.station, "name": st.name,
                         "nve_poly_km2": "" if ref is None else round(ref.area / 1e6, 4)}  # fmt: skip
            p = inp.placement(st.station)
            if p is None:
                w.writerow({**row, "refusal": "no river line within 500 m"})
                continue
            row |= {"placed_on": p.placed_on, "distance_m": round(p.distance_m, 1),
                    "objectid": p.objectid, "lake": p.lake, "confluence_near": p.confluence_near}  # fmt: skip
            try:
                _, _, path, s = first_window(inp, p)
            except (MosaicError, BurnRefusal, ValueError) as exc:
                w.writerow({**row, "refusal": str(exc).replace("\n", " ")[:200]})
                continue
            causes = [*s.causes, *([] if path.direction_ok else ["direction"])]
            w.writerow({**row, "lowered_nodes": path.lowered_nodes,
                        "lowered_max_m": round(path.lowered_max_m, 3),
                        "node_offset_m": round(path.node_offset_m, 1), "swing": round(s.swing, 4),
                        "direction_ok": path.direction_ok, "causes": ";".join(causes),
                        "well_posed_first": not causes})  # fmt: skip
            print(st.station, "survey", ";".join(causes) or "well posed", file=sys.stderr)


def overlap(fine: Polygon, ref, meta) -> dict:
    """The design's agreement numbers on the final window's node lattice."""
    x0, y0, x1, y1 = shapely.union(fine, ref).bounds
    dx, dy = meta.delta_x, meta.delta_y
    c = np.arange(math.floor((x0 - meta.x_min) / dx), math.ceil((x1 - meta.x_min) / dx) + 1)
    xs, ours, refs, both = meta.x_min + c * dx, 0, 0, 0
    shapely.prepare(fine), shapely.prepare(ref)
    for r in range(math.floor((meta.y_max - y1) / dy), math.ceil((meta.y_max - y0) / dy) + 1):
        y = np.full_like(xs, meta.y_max - r * dy)
        a = shapely.contains_xy(fine, xs, y)
        b = shapely.contains_xy(ref, xs, y)
        ours, refs, both = ours + int(a.sum()), refs + int(b.sum()), both + int((a & b).sum())
    offset = (refs + ours - 2 * both) * dx * dy / ref.length
    return {"nve_in_ours": both / refs, "ours_in_nve": both / ours,
            "area_ratio": fine.area / ref.area, "divide_offset_m": offset}  # fmt: skip


def full(inp: Inputs, sid: str, work: Path) -> dict:
    """`delineate` with the reach (what `catchment --rivers` runs), cached."""
    out = work / "runs" / f"{sid}.json"
    if out.exists():
        return json.loads(out.read_text())
    st, p = inp.by_id[sid], inp.placement(sid)
    rec: dict = {"station": sid, "name": st.name}
    t0 = time.perf_counter()
    try:
        res = delineate(CatchmentRequest(seed=(st.x, st.y), seed_crs=inp.crs, reach=p.reach), inp.repo)
    except CatchmentError as exc:
        rec["refusal"] = str(exc)
    else:
        g, s = res.gauge, res.gauge.sensitivity
        causes = [*s.causes, *([] if g.direction_ok else ["direction"])]
        rec |= {"node": g.node, "chain": g.chain, "lowered_nodes": g.lowered_nodes,
                "lowered_max_m": g.lowered_max_m, "node_offset_m": g.node_offset_m,
                "end_extended_m": g.end_extended_m, "sens": {k: getattr(s, k) for k in s.__slots__},
                "causes": causes, "well_posed": not causes, "area_km2": res.fine_area / 1e6,
                "fine": list(res.fine.exterior.coords), "windows": len(res.windows),
                **overlap(res.fine, inp.refs[sid], res.meta)}  # fmt: skip
    rec |= {"seconds": time.perf_counter() - t0, "placed_on": p.placed_on,
            "distance_m": p.distance_m, "objectid": p.objectid, "lake": p.lake,
            "confluence_near": p.confluence_near, "U": p.reach.uncertainty}  # fmt: skip
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(rec))
    return rec


def choose(inp: Inputs, work: Path) -> None:
    rows = list(csv.DictReader((work / "survey.csv").open()))
    good = [r for r in rows if r["well_posed_first"] == "True"]
    order = {s.station: k for k, s in enumerate(inp.stations)}
    # Pinned (@perf): a confluence step is `swing` as the only cause (the chain
    # drains, its end is closed, the downstream side was read, so the gain is
    # real), with one step between neighbouring chain nodes over 5 % by itself.
    confluence_step = lambda f: f.get("causes") == ["swing"] and f["sens"]["largest_step"] > SWING_MAX * f["sens"]["a0"]  # noqa: E731
    rules = {
        "1 placed straight onto the flow path": (
            sorted((r for r in good if r["lowered_nodes"] == "0"), key=lambda r: float(r["nve_poly_km2"])),
            lambda f: f.get("well_posed", False)),
        "2 the line burnt in": (
            sorted(good, key=lambda r: -float(r["lowered_max_m"])), lambda f: f.get("well_posed", False)),
        "3 marked uncertain (swing, confluence step)": (
            [r for r in rows if r["confluence_near"] == "True"], confluence_step),
        "4 a lake gauge": ([r for r in rows if r["lake"] == "True"], lambda f: f.get("well_posed", False)),
    }  # fmt: skip
    cases: dict = {}
    for name, (cands, ok) in rules.items():
        cands = sorted(cands, key=lambda r: order[r["station"]]) if name[0] in "34" else cands
        tried = []
        for r in cands:
            f = full(inp, r["station"], work)
            why = f.get("refusal") or ";".join(f.get("causes", [])) or "well posed"
            tried.append({"station": r["station"], "full_run": why, "chosen": ok(f)})
            print(name, r["station"], why, file=sys.stderr)
            if ok(f):
                break
        cases[name] = {"station": tried[-1]["station"] if tried and tried[-1]["chosen"] else None,
                       "candidates_run": tried}  # fmt: skip
    (work / "cases.json").write_text(json.dumps(cases, indent=1))


def _base(ax, tile, extent, acc_count, title):
    m = tile.meta
    z = np.asarray(tile.array, dtype=float)
    z = np.where((z == m.nodata) | np.isnan(z), np.nan, z) if m.nodata is not None else z
    r0, r1 = int((m.y_max - extent[3]) / m.delta_y), int((m.y_max - extent[2]) / m.delta_y) + 1
    c0, c1 = int((extent[0] - m.x_min) / m.delta_x), int((extent[1] - m.x_min) / m.delta_x) + 1
    img = (m.x_min + c0 * m.delta_x, m.x_min + (c1 - 1) * m.delta_x,
           m.y_max - (r1 - 1) * m.delta_y, m.y_max - r0 * m.delta_y)  # fmt: skip
    zz = z[r0:r1, c0:c1]
    hs = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(zz, nan=np.nanmin(zz)),
                                                     vert_exag=2, dx=m.delta_x, dy=m.delta_y)  # fmt: skip
    ax.imshow(hs, cmap="gray", extent=img, origin="upper", vmin=0, vmax=1, alpha=0.8)
    cnt = np.asarray(acc_count, dtype=float)[r0:r1, c0:c1]
    least = math.ceil(STREAM_KM2 * 1e6 / (m.delta_x * m.delta_y))
    flow = np.ma.masked_less(cnt, least)
    ax.imshow(flow, cmap="Blues", norm=LogNorm(least, max(least * 10, flow.max())), extent=img,
              origin="upper", interpolation="nearest")  # fmt: skip
    ax.set_xlim(extent[0], extent[1]), ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal"), ax.set_title(title, fontsize=10)
    ax.ticklabel_format(useOffset=False, style="plain"), ax.tick_params(labelsize=7)
    ax.set_xlabel(f"x, m ({tile.meta.crs})", fontsize=8), ax.set_ylabel("y, m", fontsize=8)


def _panel(ax, lines: list[str]) -> None:
    ax.axis("off")
    ax.text(0, 1, "\n".join(lines), va="top", fontsize=8, family="monospace", transform=ax.transAxes)


def draw(inp: Inputs, work: Path, out: Path, extra: list[str]) -> None:
    """Two PNGs per case, and per extra station named on the command line."""
    cases = json.loads((work / "cases.json").read_text())
    todo = [(c["station"], name) for name, c in cases.items() if c["station"] is not None]
    todo = [] if "--extras-only" in sys.argv else todo
    for sid, name in [*todo, *((s, "extra, not one of the four cases") for s in extra)]:
        f, st, p = full(inp, sid, work), inp.by_id[sid], inp.placement(sid)
        tile, burnt, path, _ = first_window(inp, p)
        m, raw = tile.meta, np.asarray(tile.array)
        xy = [(m.x_min + c * m.delta_x, m.y_max - r * m.delta_y) for r, c in path.chain.tolist()]
        print(sid, "first-window chain equals the full run's:", xy == [tuple(v) for v in f["chain"]],
              file=sys.stderr)  # fmt: skip
        line = np.asarray(p.reach.line)
        (x0, y0), (x1, y1) = line.min(axis=0) - MAP_PAD_M, line.max(axis=0) + MAP_PAD_M
        ext, box = (x0, x1, y0, y1), shapely.box(x0, y0, x1, y1)
        others = [s for s in inp.segments if s.objectid != p.objectid and LineString(s.line).intersects(box)]
        seg = next(s for s in inp.segments if s.objectid == p.objectid)
        verdict = "well defined" if f["well_posed"] else "uncertain: " + ", ".join(f["causes"])
        for when, arr in (("before", tile), ("after", burnt)):
            fig, (ax, side) = plt.subplots(1, 2, figsize=(13, 7.5), dpi=120, width_ratios=[2.4, 1])
            how = "batch's tiered placement" if inp.tiers else "nearest line, as catchment --rivers"
            _base(ax, tile, ext, accumulate(to_core(arr)).count, f"{sid} {st.name}, case {name}: {when}\n({how})")
            for k, s in enumerate(others):
                ax.plot(*np.asarray(s.line).T, color="tab:olive", lw=0.8, ls="--",
                        label="other NVE river lines" if k == 0 else None)  # fmt: skip
            if when == "before":
                ax.plot(*line.T, color="tab:green", lw=1.4, label="NVE river line: the reach")
                ax.plot(*np.asarray(seg.line).T, color="tab:green", lw=3.5, alpha=0.45,
                        label=f"chosen segment {p.objectid}")  # fmt: skip
                ax.plot(*p.position, "o", mfc="none", mec="k", ms=10, mew=1.5, label="mapped position P")
                _panel(side, [
                    "Before: the raw DEM.", "Blue: DEM flow paths, the nodes",
                    f"draining >= {STREAM_KM2} km2, darker as more drain.", "",
                    f"station  {sid} {st.name}", f"number   {st.watercourse}", f"river    {st.river}", "",
                    f"chosen   segment {p.objectid}", f"         {seg.name or '(no name)'}, {seg.vassdragsnr}",
                    f"         {seg.kind} line, tier: {p.placed_on}", "",
                    f"P is {p.distance_m:.1f} m from the station", f"U = {p.reach.uncertainty:.0f} m",
                    f"another river line within 100 m of P: {p.confluence_near}"])  # fmt: skip
            else:
                ch = np.asarray(f["chain"])
                arc = np.r_[0.0, np.cumsum(np.hypot(*np.diff(ch, axis=0).T))]
                arc -= arc[int(np.argmin(np.hypot(*(ch - f["node"]).T)))]
                low = [k for k, (r, c) in enumerate(path.chain.tolist()) if burnt.array[r, c] < raw[r, c]]
                ax.plot(ch[low, 0], ch[low, 1], ".", color="darkred", ms=4,
                        label=f"lowered chain nodes ({f['lowered_nodes']})")  # fmt: skip
                ax.plot(*ch.T, color="tab:orange", lw=1.2, label="burnt chain")
                inw = np.abs(arc) <= f["U"] + 1e-6
                ax.plot(*ch[inw].T, color="magenta", lw=5, alpha=0.45,
                        label=f"sensitivity window, U = {f['U']:.0f} m each way")  # fmt: skip
                ax.plot(*f["node"], "*", color="yellow", mec="k", ms=16, label="placed node", zorder=5)
                ins = side.inset_axes([0.0, 0.0, 1.0, 0.5])
                for part in shapely.get_parts(inp.refs[sid]):
                    ins.fill(*part.exterior.xy, color="tab:green", alpha=0.35, lw=0)
                ins.plot(*np.asarray(f["fine"]).T, color="k", lw=0.8)
                ins.plot(*f["node"], "*", color="yellow", mec="k", ms=10)
                ins.set_aspect("equal"), ins.set_xticks([]), ins.set_yticks([])
                ins.set_title("ours (black line) over NVE's polygon (green)", fontsize=8)
                s = f["sens"]
                _panel(side, [
                    "After: the reach burnt in.", "Blue: flow paths of the burnt DEM.", "",
                    f"placed node {f['node_offset_m']:.1f} m from P",
                    f"lowered {f['lowered_nodes']} nodes, at most {f['lowered_max_m']:.2f} m",
                    f"chain end extended {f['end_extended_m']:.0f} m", "",
                    f"area at the node {s['a0']:.4g} km2", f"swing within U   {100 * s['swing']:.1f} %",
                    f"largest step {s['largest_step']:.4g} km2, at {s['largest_step_at_m']:.0f} m",
                    f"verdict: {verdict}", "",
                    f"NVE's polygon {inp.refs[sid].area / 1e6:.4g} km2, ours {f['area_km2']:.4g} km2",
                    f"NVE's in ours {100 * f['nve_in_ours']:.1f} %, ours in NVE's {100 * f['ours_in_nve']:.1f} %",
                    f"area ratio {f['area_ratio']:.3g}"])  # fmt: skip
            ax.plot(st.x, st.y, "r^", ms=10, label="NVE station point", zorder=5)
            ax.legend(loc="lower left", fontsize=7, framealpha=0.85)
            fig.savefig(out / f"{sid}{'_tiers' if inp.tiers else ''}_{when}.png", bbox_inches="tight")
            plt.close(fig)
        print("drew", sid, file=sys.stderr)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a not in ("--tiers", "--extras-only")]
    cmd, data, dem, work = args[0], Path(args[1]), Path(args[2]), Path(args[3])
    work.mkdir(parents=True, exist_ok=True)
    inp = Inputs(data, dem, "--tiers" in sys.argv)
    {"survey": lambda: survey(inp, work), "choose": lambda: choose(inp, work),
     "draw": lambda: draw(inp, work, Path(args[4]), args[5:])}[cmd]()  # fmt: skip
