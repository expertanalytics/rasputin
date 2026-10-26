"""Fix 3 probe, round 3: the picture docs/increments/106-probes/kultsjon-budget-borders.png.
Reads budget_band.py's OUT (src-NAME, sub-N, dom-N .geojson) and OUT/mN.gpkg, the mesh of
each budget N (`rasputin mesh --out OUT/mN.gpkg --dem glo30 --out-crs EPSG:3006
--domain OUT/dom-N.geojson --features OUT/sub-N.geojson --triangles N`); prints the
triangles, the border vertices and the mesh vertices on the outline, and draws SMHI's
outline, the simplified border and the mesh. Needs matplotlib (the main checkout's venv).
Run:  python -I docs/increments/106-probes/budget_borders_plot.py OUT NAME PNG N [N ...]
"""

import json
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np
import shapely
from shapely.geometry import shape

from tin_engine.feature_input import read_source

mpl.use("Agg")
import matplotlib.pyplot as plt

d, name, png = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
budgets = [int(x) for x in sys.argv[4:]]


def polygons(path: Path) -> list:
    return [shape(f["geometry"]) for f in json.loads(path.read_text())["features"]]


def outline_vertices(tris: list) -> int:
    edges: dict[tuple, int] = {}
    for t in tris:
        c = [tuple(np.round(p, 3)) for p in shapely.get_coordinates(t)[:3]]
        for a, b in ((c[0], c[1]), (c[1], c[2]), (c[2], c[0])):
            k = (min(a, b), max(a, b))
            edges[k] = edges.get(k, 0) + 1
    return len({p for k, v in edges.items() if v == 1 for p in k})


def draw(ax, geoms: list, **style) -> None:
    for g in geoms:
        for p in shapely.get_parts(g):
            x, y = p.exterior.xy
            ax.plot(np.array(x) / 1e3, np.array(y) / 1e3, **style)


src = polygons(d / f"src-{name}.geojson")
fig, axs = plt.subplots(1, len(budgets), figsize=(7 * len(budgets), 5.6))
for ax, n in zip(axs, budgets, strict=True):
    rows = read_source(d / f"m{n}.gpkg", "triangles", None, lambda _c: (-1e12, -1e12, 1e12, 1e12))
    tris = [shapely.force_2d(g) for _, g, _ in rows.rows]
    sub = polygons(d / f"sub-{n}.geojson")
    sv = np.unique(np.concatenate([shapely.get_coordinates(s) for s in sub]), axis=0)
    print(
        f"{name} N {n}: triangles {len(tris)}, simplified border vertices {len(sv)}, "
        f"mesh vertices on the outline {outline_vertices(tris)}"
    )
    draw(ax, tris, color="#c8c8c8", lw=0.4, zorder=1)
    draw(ax, src, color="#1f5fa8", lw=1.6, zorder=2, label="SMHI's outline (SVAR 2022)")
    draw(ax, sub, color="#d4561a", lw=1.2, zorder=3, label="border used as constraint")
    ax.plot(sv[:, 0] / 1e3, sv[:, 1] / 1e3, "o", color="#d4561a", ms=3, zorder=4)
    ax.set_title(f"--triangles {n}: {len(tris)} triangles, {len(sv)} border vertices")
    ax.set_aspect("equal")
    ax.set_xlabel("km east (EPSG:3006)")
    ax.set_ylabel("km north")
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles, strict=True))
    ax.legend(unique.values(), unique.keys(), loc="lower left", fontsize=9)
fig.suptitle("Kultsjön, SMHI's main sub-basin: border band = triangle edge / 8; grey = mesh")
fig.tight_layout()
fig.savefig(png, dpi=110)
