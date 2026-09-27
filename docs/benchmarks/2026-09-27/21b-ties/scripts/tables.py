"""tables.py data/ties_quarter_t1.jsonl data/ties_tile_t1.jsonl > data/tables.md

Summarises instrument.patch's RASPUTIN_TIES_OUT records (one JSON line per
refine call) into the README's tables. Phase p2 is the refine loop, p0
legalise_all, p1 the quality pass.
"""
import json
import re
import sys
from pathlib import Path


def load(path):
    return json.loads(Path(path).read_text().splitlines()[-1])


def shape_group(name):
    if name.startswith("axis_square"):
        return "axis-aligned square"
    if name.startswith("axis_rectangle"):
        return "axis-aligned rectangle, not square"
    return {
        "rotated_square": "rotated square",
        "rotated_rectangle": "rotated rectangle, not square",
        "isosceles_trapezoid_axis": "isosceles trapezoid, parallel sides on rows or columns",
        "isosceles_trapezoid_other": "isosceles trapezoid, other direction",
        "other_cyclic": "other cyclic quad (no parallel sides)",
        "collinear_triple": "three corners collinear",
    }[name]


def pct(n, d):
    return f"{100.0 * n / d:.3f} %" if d else "-"


def main(paths):
    recs = {Path(p).stem.split("_")[1]: load(p) for p in paths}
    names = list(recs)
    print("| measure | " + " | ".join(names) + " |")
    print("|---|" + "---:|" * len(names))
    rows = []
    for ph, label in (("p2", "refine loop"), ("p0", "legalise_all"), ("p1", "quality pass")):
        def g(r, k, ph=ph):
            return r["counts"].get(f"{ph}.{k}", 0)

        rows += [
            (f"{label}: incircle calls", lambda r, g=g: f"{g(r, 'calls'):,}"),
            (f"{label}: exact path", lambda r, g=g: f"{g(r, 'exact.calls'):,} ({pct(g(r, 'exact.calls'), g(r, 'calls'))})"),
            (f"{label}: exact path, Cocircular", lambda r, g=g: f"{g(r, 'exact.result.cocircular'):,}"),
            (f"{label}: exact path, four node corners", lambda r, g=g: f"{g(r, 'exact.nodes4'):,}"),
            (f"{label}: exact path, lattice det = 0", lambda r, g=g: f"{g(r, 'exact.n4.lattice_det_zero'):,}"),
            (f"{label}: exact path, lattice det != 0", lambda r, g=g: f"{g(r, 'exact.n4.exact_path_lattice_nonzero'):,}"),
            (f"{label}: exact path, QW2 conditions hold", lambda r, g=g: f"{g(r, 'exact.n4.qw2_qualifies'):,}"),
            (f"{label}: filtered path, four node corners", lambda r, g=g: f"{g(r, 'filtered.nodes4'):,}"),
            (f"{label}: filtered path, QW2 conditions hold", lambda r, g=g: f"{g(r, 'filtered.n4.qw2_qualifies'):,}"),
            (f"{label}: all calls QW2 would answer", lambda r, g=g: (
                lambda q: f"{q:,} ({pct(q, g(r, 'calls'))})")(
                g(r, 'exact.n4.qw2_qualifies') + g(r, 'filtered.n4.qw2_qualifies'))),
            (f"{label}: calls with fewer than four node corners", lambda r, g=g: f"{sum(g(r, f'{p}.nodes{k}') for p in ('exact', 'filtered') for k in range(4)):,}"),
            (f"{label}: lattice sign differs from kernel (four nodes)", lambda r, g=g: f"{g(r, 'exact.n4.lattice_sign_differs') + g(r, 'filtered.n4.lattice_sign_differs'):,}"),
        ]
    rows += [
        ("max spread from d, four-node calls, refine loop (nodes)", lambda r: f"{r['max_spread_all4_refine']}"),
        ("max spread from d, exact-path calls, refine loop (nodes)", lambda r: f"{r['max_spread_exact4_refine']}"),
        ("dx, dy; rows x cols", lambda r: f"{r['dx']:g}, {r['dy']:g}; {r['rows']} x {r['cols']}"),
        ("frame-exact sufficient condition (bits of dx + bit_width) <= 53", lambda r: f"{r['dx_significant_bits']} + {r['index_bit_width']} = {r['dx_significant_bits'] + r['index_bit_width']}: {r['global_frame_exact_condition']}"),
        ("rounds, inserted, flips", lambda r: f"{r['rounds']}, {r['inserted']:,}, {r['flips']:,}"),
    ]
    for label, f in rows:
        print(f"| {label} | " + " | ".join(f(recs[n]) for n in names) + " |")

    print()
    print("Shapes of the exact-path (tie) quads, refine loop:")
    print()
    print("| shape | " + " | ".join(names) + " |")
    print("|---|" + "---:|" * len(names))
    groups = {}
    for n in names:
        tot = recs[n]["counts"].get("p2.exact.calls", 0)
        for k, v in recs[n]["counts"].items():
            m = re.fullmatch(r"p2\.exact\.n4\.shape\.(.+)", k)
            if m:
                groups.setdefault(shape_group(m.group(1)), {}).setdefault(n, 0)
                groups[shape_group(m.group(1))][n] += v
    order = sorted(groups, key=lambda s: -groups[s].get(names[0], 0))
    for s in order:
        cells = []
        for n in names:
            v = groups[s].get(n, 0)
            cells.append(f"{v:,} ({pct(v, recs[n]['counts'].get('p2.exact.calls', 0))})")
        print(f"| {s} | " + " | ".join(cells) + " |")
    print()
    print("Most frequent axis-aligned sizes (short x long side, nodes), refine loop exact path:")
    print()
    for n in names:
        sizes = sorted(((v, k.split(".")[-1]) for k, v in recs[n]["counts"].items()
                        if k.startswith("p2.exact.n4.shape.axis_")), reverse=True)[:6]
        print(f"- {n}: " + ", ".join(f"{s.removeprefix('axis_')} {v:,}" for v, s in sizes))


if __name__ == "__main__":
    main(sys.argv[1:])
