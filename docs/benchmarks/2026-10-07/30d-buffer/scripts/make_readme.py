"""Write ../README.md from README.template.md and the files under ../raw/.

Every table in the README is generated here; the prose is the template's.
Usage: python make_readme.py   (no arguments; reads and writes beside itself)
"""

import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RAW = ROOT / "raw"
RUNS = RAW / "runs"
CATCH = {"lagan": "Lagan", "ljungan_flasjo": "Ljungan above Flåsjö", "numedalslagen": "Numedalslågen"}
INPUT = {"geojson": "GeoJSON window", "gpkg": "European GeoPackage", "gpkg33": "UTM33 GeoPackage"}
PAIRS = [("lagan", "geojson"), ("lagan", "gpkg"), ("ljungan_flasjo", "geojson"),
         ("ljungan_flasjo", "gpkg"), ("numedalslagen", "gpkg33"), ("numedalslagen", "gpkg")]  # fmt: skip
PHASES = ["decode", "features read", "features clip", "total"]
VARIANTS = [("master", "master"), ("p1", "1: buffer once"), ("p12", "1+2: + region from the hull"),
            ("p123k1000", "1+2+3, pieces of 1000"), ("p123k500", "1+2+3, pieces of 500"),
            ("p123k250", "1+2+3, pieces of 250")]  # fmt: skip


def stats(path: Path) -> dict[str, object]:
    """Phase seconds and quality figures from one --stats file."""
    text = path.read_text()
    out: dict[str, object] = {}
    for m in re.finditer(r"^\| (?:\*\*)?([a-z: +()]+?)(?:\*\*)? \| (?:\*\*)?([0-9.]+)(?:\*\*)? \| ", text, re.M):
        out.setdefault(m.group(1), float(m.group(2)))
    m = re.search(r"^\| minimum angle \| ([0-9.]+)° \| ([0-9.]+) % \| ([0-9.]+) % \| ([0-9.]+)° \|", text, re.M)
    out["worst angle"] = m.group(4) if m else "?"
    m = re.search(r"^\| vertex degree \(triangles\) \| \d+ \| \d+ \| (\d+) \|", text, re.M)
    out["max degree"] = m.group(1) if m else "?"
    m = re.search(r"^\| Largest height error \| ([0-9.]+) \|", text, re.M)
    out["max error"] = m.group(1) if m else "?"
    m = re.search(r"^\| Self-check: DEM nodes with data left outside the mesh \| (\d+) \|", text, re.M)
    out["outside"] = m.group(1) if m else "?"
    m = re.search(r"^\| output triangles \| (\d+) \|", text, re.M)
    out["triangles"] = m.group(1) if m else "?"
    return out


def power_ok(stem: str) -> bool:
    lines = (RUNS / f"{stem}_power.txt").read_text().splitlines()
    return len(lines) == 2 and all("AC Power" in x for x in lines)


def exit_ok(stem: str) -> bool:
    return (RUNS / f"{stem}_stderr.txt").read_text().rstrip().endswith("exit 0")


SHA = dict(line.split() for line in (RUNS / "vtk_sha256.txt").read_text().splitlines())


def runs(tag: str, c: str, f: str) -> list[dict[str, object]]:
    out = []
    for r in ("r1", "r2", "r3"):
        stem = f"{tag}_{r}_{c}_{f}"
        assert power_ok(stem), f"{stem}: not on AC power before and after"
        assert exit_ok(stem), f"{stem}: did not exit 0"
        out.append(stats(RUNS / f"{stem}_stats.md"))
    return out


def med(rows: list[dict[str, object]], key: str) -> float:
    return statistics.median(float(r[key]) for r in rows)  # type: ignore[arg-type]


def cell(rows: list[dict[str, object]], key: str) -> str:
    vals = ", ".join(f"{float(r[key]):.2f}" for r in rows)  # type: ignore[arg-type]
    return f"**{med(rows, key):.2f}** ({vals})"


def table_master() -> str:
    head = "| catchment | features | " + " | ".join(PHASES) + " |\n|---|---|" + "---|" * len(PHASES)
    body = []
    for c, f in PAIRS:
        rows = runs("master", c, f)
        body.append(f"| {CATCH[c]} | {INPUT[f]} | " + " | ".join(cell(rows, p) for p in PHASES) + " |")
    return "\n".join([head, *body])


def trace(c: str, f: str) -> list[dict[str, object]]:
    return [json.loads(x) for x in (RAW / "trace" / f"{c}_{f}.jsonl").read_text().splitlines()]


def table_trace() -> str:
    head = ("| catchment | features | call site | calls | distance | join | vertices in | seconds |\n"
            "|---|---|---|---|---|---|---|---|")  # fmt: skip
    body = []
    for c, f in PAIRS:
        agg: dict[tuple[str, str, object, str], list[float]] = defaultdict(list)
        verts: dict[tuple[str, str, object, str], set[int]] = defaultdict(set)
        for r in trace(c, f):
            if r["call"] != "buffer" and r["call"] != "query_features" and r["call"] != "read_json":
                continue
            join = (r.get("kwargs") or {}).get("join_style", "round") if r["call"] == "buffer" else ""
            k = (str(r["call"]), str(r["where"]), r.get("distance"), str(join))
            agg[k].append(float(r["s"]))  # type: ignore[arg-type]
            verts[k].add(int(r.get("vertices") or r.get("rows") or 0))  # type: ignore[arg-type]
        for k, ss in sorted(agg.items(), key=lambda kv: -sum(kv[1])):
            call, where, d, join = k
            name = f"`{where}`" if call == "buffer" else f"`{call}` (`{where}`)"
            dist = f"{d:.4g}" if isinstance(d, float) else ""
            v = "/".join(str(x) for x in sorted(verts[k])[:3]) + ("…" if len(verts[k]) > 3 else "")
            if call == "query_features":
                v += " rows"
            if call == "read_json":
                v = ""
            body.append(f"| {CATCH[c]} | {INPUT[f]} | {name} | {len(ss)} | {dist} | {join} | {v} | {sum(ss):.3f} |")
        st = stats(RUNS / f"trace_{c}_{f}_stats.md")
        bufs = defaultdict(float)
        for r in trace(c, f):
            if r["call"] == "buffer":
                site = str(r["where"])
                phase = "features read" if site.startswith("feature_input") else "decode"
                bufs[phase] += float(r["s"])  # type: ignore[arg-type]
        for phase in ("decode", "features read"):
            body.append(f"| {CATCH[c]} | {INPUT[f]} | *{phase}* in this run: {float(st[phase]):.3f} s, "  # type: ignore[arg-type]
                        f"of which buffers {bufs[phase]:.3f} s | | | | | "
                        f"rest {float(st[phase]) - bufs[phase]:.3f} |")  # type: ignore[arg-type]  # fmt: skip
    return "\n".join([head, *body])


def probes(name: str) -> list[dict[str, object]]:
    p = RAW / "probes" / f"{name}.jsonl"
    power = (RAW / "probes" / f"{name}_power.txt").read_text().splitlines()
    assert len(power) == 2 and all("AC Power" in x for x in power), f"{name}: power"
    assert (RAW / "probes" / f"{name}_stderr.txt").read_text().rstrip().endswith("exit 0"), name
    return [json.loads(x) for x in p.read_text().splitlines()]


def secs(r: dict[str, object]) -> str:
    vals = ", ".join(f"{x:.3f}" for x in r["runs"])  # type: ignore[union-attr]
    return f"**{float(r['median']):.3f}** ({vals})"  # type: ignore[arg-type]


def table_isolated() -> str:
    head = "| catchment | step | seconds, median (3 runs) |\n|---|---|---|"
    body, sizes = [], []
    for c in CATCH:
        for r in probes(f"isolated_{c}"):
            if "step" in r:
                extra = f" (d = {float(r['distance']):.4g}, {r['vertices']} vertices)" if "distance" in r else ""  # type: ignore[arg-type]
                body.append(f"| {CATCH[c]} | {r['step']}{extra} | {secs(r)} |")
            elif "source" in r:
                sizes.append(f"| {CATCH[c]} | {INPUT[str(r['source'])]} | {r['rows']} | {int(r['vertices']):,} |")  # type: ignore[arg-type]
    rows = "| catchment | features | rows read | vertices read |\n|---|---|---|---|"
    return "\n".join([head, *body]) + "\n\nWhat each input hands to the clip:\n\n" + "\n".join([rows, *sizes])


def table_candidates() -> str:
    head = ("| catchment | features | variant | decode | features read | features clip | total | "
            "`.vtk` vs master | worst angle | max degree | largest error m | DEM nodes outside |\n"
            "|---|---|---|---|---|---|---|---|---|---|---|---|")  # fmt: skip
    body = []
    for c, f in PAIRS[:4]:  # the patched runs were made on Lagan and Ljungan only
        base = {SHA[f"master_{r}_{c}_{f}"] for r in ("w0", "r1", "r2", "r3")}
        assert len(base) == 1, f"master {c} {f}: meshes differ between repeats"
        for tag, label in VARIANTS:
            rows = runs(tag, c, f)
            shas = {SHA[f"{tag}_{r}_{c}_{f}"] for r in ("w0", "r1", "r2", "r3")}
            same = "identical" if shas == base else f"**differs** ({len(shas)} hash(es))"
            q = rows[0]
            body.append(f"| {CATCH[c]} | {INPUT[f]} | {label} | " + " | ".join(cell(rows, p) for p in PHASES)
                        + f" | {same} | {q['worst angle']}° | {q['max degree']} | {q['max error']} | {q['outside']} |")  # fmt: skip
    return "\n".join([head, *body])


def table_pieces() -> str:
    head = ("| catchment | method | seconds, median (3 runs) | sym. difference m² | largest boundary distance m "
            "| pokes out of GEOS's m | GEOS's pokes out m | bounds, largest difference m | every vertex within 1e-6 m |\n"
            "|---|---|---|---|---|---|---|---|---|")  # fmt: skip
    body = []
    for c in CATCH:
        for r in probes(f"candidates_{c}_pieces"):
            if r["method"] == "GEOS buffer":
                body.append(f"| {CATCH[c]} ({r['edges']} edges, d = {float(r['distance']):.4g} m) | GEOS buffer, mitre | {secs(r)} | | | | | | |")  # type: ignore[arg-type]
                continue
            far = max(float(r["far_got_to_exact_m"]), float(r["far_exact_to_got_m"]))  # type: ignore[arg-type]
            body.append(f"| {CATCH[c]} | {r['method']} | {secs(r)} | {float(r['sym_diff_m2']):.3g} | {far:.3g} | "  # type: ignore[arg-type]
                        f"{float(r['got_outside_exact_m']):.3g} | {float(r['exact_outside_got_m']):.3g} | "  # type: ignore[arg-type]
                        f"{float(r['bounds_max_diff_m']):.3g} | {'yes' if r['vertices_within_1e6'] else 'no'} |")  # type: ignore[arg-type]  # fmt: skip
    return "\n".join([head, *body])


def table_region() -> str:
    head = ("| catchment | master: hull of the grown domain, s | candidate: the hull grown, s | sym. difference m² "
            "| largest boundary distance m | bounds, largest difference m | candidate covers the domain grown by 99.5 m |\n"
            "|---|---|---|---|---|---|---|")  # fmt: skip
    body = []
    for c in CATCH:
        m, k = probes(f"candidates_{c}_region")
        far = max(float(k["far_cand_to_master_m"]), float(k["far_master_to_cand_m"]))  # type: ignore[arg-type]
        body.append(f"| {CATCH[c]} | {secs(m)} | {secs(k)} | {float(k['sym_diff_m2']):.4g} | {far:.3g} | "  # type: ignore[arg-type]
                    f"{float(k['bounds_max_diff_m']):.3g} | {'yes' if k['cand_covers_domain_grown_99_5'] else 'no'} |")  # type: ignore[arg-type]  # fmt: skip
    return "\n".join([head, *body])


def table_scaling() -> str:
    head = "| catchment | shape | join | distance m | vertices | seconds, median (3 runs) |\n|---|---|---|---|---|---|"
    body = []
    for c in CATCH:
        for r in probes(f"candidates_{c}_scaling"):
            body.append(f"| {CATCH[c]} | {r['shape']} | {r['join']} | {float(r['distance']):.4g} | {r['vertices']} | {secs(r)} |")  # type: ignore[arg-type]
    return "\n".join([head, *body])


def main() -> None:
    text = (HERE / "README.template.md").read_text()
    for name, f in {"MASTER": table_master, "TRACE": table_trace, "ISOLATED": table_isolated,
                    "CANDIDATES": table_candidates, "PIECES": table_pieces, "REGION": table_region,
                    "SCALING": table_scaling}.items():  # fmt: skip
        key = "{{" + name + "}}"
        if key in text:
            text = text.replace(key, f())
    (ROOT / "README.md").write_text(text)


if __name__ == "__main__":
    main()
