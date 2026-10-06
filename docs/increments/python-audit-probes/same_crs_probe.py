"""The probe set for `same_crs` (python-audit.md, section 9, "the rule" and its
probe table), run against the rule `5197f9a` implements and the rule ruled
after code review round 1. Both rules are copied here, so the script does not
depend on which of them `src_python/` holds.

Run from the repository root with the project venv (about 5 minutes):

    .venv/bin/python docs/increments/python-audit-probes/same_crs_probe.py

Prints, per probe pair: both rules' answers in both argument orders, and the
largest move over a 15 by 15 grid on the EPSG code's area of use. Then, over
every CRS in the set: ordered same pairs, asymmetric pairs, transitivity
breaks, and the CRSs not the same as themselves. Then the tolerance: one
parameter of EPSG:25833 or 31287 perturbed, and whether the rule still calls
it the same.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable
from typing import Any

import numpy as np
from pyproj import CRS, Transformer
from pyproj.exceptions import ProjError

warnings.simplefilter("ignore")


def no_id(code: int) -> dict[str, Any]:
    j = CRS.from_epsg(code).to_json_dict()
    j.pop("id", None)
    return j


def wkt_no_id(code: int) -> str:
    return CRS.from_json_dict(no_id(code)).to_wkt()


def axes_swapped(code: int) -> str:
    j = no_id(code)
    axes = j["coordinate_system"]["axis"]
    j["coordinate_system"]["axis"] = [axes[1], axes[0]]
    return CRS.from_json_dict(j).to_wkt()


def in_grads(code: int) -> str:
    j = no_id(code)
    for axis in j["coordinate_system"]["axis"]:
        axis["unit"] = {"type": "AngularUnit", "name": "grad", "conversion_factor": math.pi / 200}
    return CRS.from_json_dict(j).to_wkt()


def proj4(code: int) -> str:
    return CRS.from_epsg(code).to_proj4()


U = "+proj=utm +zone=33 +ellps=GRS80 +units=m +no_defs"
L = (
    "+proj=lcc +lat_1=46 +lat_2=49 +lat_0=47.5 +lon_0=13.33333333333333 "
    "+x_0=400000 +y_0=400000 +ellps=bessel +units=m +no_defs"
)
TM = "+proj=tmerc +lat_0=0 +lon_0=15 +k=0.9996 +x_0=500000 +y_0=0 +ellps=GRS80 +units=m +no_defs"
G = "+proj=longlat +ellps=GRS80 +no_defs"
LL = "+proj=longlat +datum=WGS84 +no_defs"
TOWGS84 = "+towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232"
LONG_ISLAND = (
    "+proj=lcc +lat_0=40.1666666666667 +lon_0=-74 +lat_1=41.0333333333333 "
    "+lat_2=40.6666666666667 +x_0=300000 +y_0=0 +datum=NAD83 +units=us-ft +no_defs"
)
DYNAMIC_2020 = CRS.from_epsg(9057).to_wkt().replace("FRAMEEPOCH[2005]", "FRAMEEPOCH[2020]")

Pair = tuple[str | int, str | int]
ROWS: list[tuple[str, list[Pair]]] = [
    (
        "axis order",
        [
            (CRS.from_epsg(4326).to_wkt("WKT1_GDAL"), 4326),
            (axes_swapped(25833), 25833),
            (axes_swapped(31287), 31287),
            (axes_swapped(3035), 3035),
        ],
    ),
    ("axis order, no datum", [(U + " +axis=neu", 25833)]),
    (
        "axis direction",
        [
            (U + " +axis=wnu", 25833),
            (U + " +axis=esu", 25833),
            (U + " +axis=wsu", 25833),
            (L + " +axis=wnu", 31287),
            (G + " +axis=wnu", 4258),
        ],
    ),
    (
        "linear unit",
        [
            (U.replace("+units=m", "+units=us-ft"), 25833),
            (U.replace("+units=m", "+units=km"), 25833),
            (U.replace("+units=m", "+to_meter=1.0000001"), 25833),
            (U.replace("+units=m", "+to_meter=1.000000001"), 25833),
            (L.replace("+units=m", "+units=us-ft"), 31287),
            ("EPSG:2263", 32118),
        ],
    ),
    ("angular unit", [(in_grads(4326), 4326)]),
    (
        "prime meridian",
        [(U + " +pm=paris", 25833), (G + " +pm=ferro", 4258), ("EPSG:31251", 31254)],
    ),
    (
        "conversion parameter held in the remark",
        [
            (LL.replace("+no_defs", "+lon_0=10 +no_defs"), 4326),
            (LL.replace("+no_defs", "+lon_0=-3 +no_defs"), 4326),
            ("+proj=longlat +datum=NAD83 +lon_0=10 +no_defs", 4269),
            (
                LL.replace("+no_defs", "+lon_0=10 +no_defs"),
                LL.replace("+no_defs", "+lon_0=20 +no_defs"),
            ),
            (LL.replace("+no_defs", "+lon_0=10 +no_defs"), "OGC:CRS84"),
        ],
    ),
    (
        "datum",
        [
            ("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs", 25833),
            (U.replace("GRS80", "WGS84"), 25833),
            (U + " +towgs84=0,0,0,0,0,0,0", 25833),
            (U + " +towgs84=0,0,0", 25833),
            (L + " +nadgrids=@null", 31287),
            (L + " " + TOWGS84, 31287),
            (G, 4258),
            (G, 4269),
            (G, 4283),
            ("EPSG:4269", 4258),
            ("EPSG:4258", 4326),
            ("EPSG:26917", 6346),
        ],
    ),
    (
        "no datum named",
        [
            (L, 31287),
            (L.replace("+lat_1=46 +lat_2=49", "+lat_1=49 +lat_2=46"), 31287),
            (proj4(31287), 31287),
            (proj4(25833), 25833),
            (proj4(3035), 3035),
        ],
    ),
    (
        "method, no datum",
        [(TM, 25833), (TM.replace("tmerc", "etmerc"), 25833), (TM + " +approx", 25833)],
    ),
    (
        "parameters, no datum",
        [
            (L.replace("13.33333333333333", "13.33333333533333"), 31287),
            (TM.replace("0.9996", repr(0.9996 * (1 + 2e-10))), 25833),
            (TM.replace("+x_0=500000", "+x_0=500000.00001"), 25833),
            (U + " +k=0.5", 25833),
            (U + " +x_0=1000", 25833),
        ],
    ),
    ("parameters", [(L.replace("13.33333333333333", "13.5"), 31287), ("EPSG:25832", 25833)]),
    ("hemisphere", [(U + " +south", 25833)]),
    ("dimension", [("EPSG:25833+5941", 25833), ("EPSG:4937", 4258), ("EPSG:4936", 4258)]),
    ("vertical", [(U + " +vunits=ft", 25833), (U + " +geoidgrids=@egm96_15.gtx", 25833)]),
    ("longitude range", [(G + " +lon_wrap=180", 4258), (G + " +over", 4258)]),
    (
        "derived",
        [
            (
                "+proj=ob_tran +o_proj=longlat +o_lat_p=40 +o_lon_p=0 +lon_0=10 "
                "+ellps=GRS80 +no_defs",
                4258,
            )
        ],
    ),
    ("epoch", [(DYNAMIC_2020, 9057)]),
    (
        "same by definition",
        [
            (wkt_no_id(31287), 31287),
            (wkt_no_id(25833), 25833),
            (CRS.from_epsg(3035).to_wkt("WKT1_GDAL"), 3035),
            (CRS.from_epsg(3035).to_wkt("WKT1_ESRI"), 3035),
            (CRS.from_epsg(25833).to_wkt("WKT1_GDAL"), 25833),
            (CRS.from_epsg(25833).to_wkt("WKT1_ESRI"), 25833),
            (CRS.from_epsg(31287).to_wkt("WKT1_GDAL"), 31287),
            ("OGC:CRS84", 4326),
            ("EPSG:3045", 25833),
            ("epsg:31287", 31287),
            ("ESRI:102100", 3857),
            ("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs", 32633),
            (LL, 4326),
            ("+proj=longlat +datum=NAD83 +no_defs", 4269),
            (LONG_ISLAND, 2263),
        ],
    ),
    ("other body", [("+proj=longlat +a=3396190 +b=3376200 +no_defs", 4326)]),
]


def name(x: str | int) -> str:
    return f"EPSG:{x}" if isinstance(x, int) else x


def transformer(a: str, b: str) -> Transformer | None:
    try:
        t = Transformer.from_crs(CRS(a), CRS(b), always_xy=True)
    except ProjError:
        return None
    return None if t.source_crs is None or t.target_crs is None else t


def rule_now(a: str, b: str) -> bool:
    """Equivalent once x-then-y, and the operation is PROJ's noop."""
    t = transformer(a, b)
    if t is None or t.source_crs is None or t.target_crs is None:
        return False
    return bool(t.source_crs.equals(t.target_crs)) and t.definition.split(" ", 1)[0] == "proj=noop"


def _same_frame(a: CRS, b: CRS) -> bool:
    if len(a.axis_info) != len(b.axis_info):
        return False
    for u, v in zip(a.axis_info, b.axis_info, strict=True):
        if u.direction.lower() != v.direction.lower() or not math.isclose(
            u.unit_conversion_factor, v.unit_conversion_factor, rel_tol=1e-12
        ):
            return False
    pm = [
        0.0
        if c.prime_meridian is None
        else c.prime_meridian.longitude * float(c.prime_meridian.unit_conversion_factor)
        for c in (a, b)
    ]
    return math.isclose(pm[0], pm[1], rel_tol=0.0, abs_tol=1e-12)


def rule_5197f9a(a: str, b: str) -> bool:
    """Equivalent once x-then-y, or an EPSG code in common at 70 in the same frame."""
    t = transformer(a, b)
    if t is None or t.source_crs is None or t.target_crs is None:
        return False
    if t.source_crs.equals(t.target_crs):
        return True
    codes = [{m.code for m in CRS(x).list_authority("EPSG", 70)} for x in (a, b)]
    return bool(codes[0] & codes[1]) and _same_frame(t.source_crs, t.target_crs)


def moves(a: str, b: str) -> str:
    base = CRS(b if b.startswith(("EPSG", "OGC", "ESRI")) else a)
    if base.area_of_use is None:
        return "n/a"
    w, s, e, n = base.area_of_use.bounds
    lon, lat = np.meshgrid(np.linspace(w, e, 15), np.linspace(s, n, 15))
    try:
        geo = base.geodetic_crs
        xa, ya = Transformer.from_crs(geo, CRS(a), always_xy=True).transform(
            lon.ravel(), lat.ravel()
        )
        xb, yb = Transformer.from_crs(CRS(a), CRS(b), always_xy=True).transform(xa, ya)
    except ProjError:
        return "no operation"
    scale = [
        ax.unit_conversion_factor
        * (6_371_000 if ax.unit_name in ("degree", "grad", "radian") else 1)
        for ax in CRS(b).axis_info[:2]
    ]
    d = np.hypot((np.asarray(xb) - xa) * scale[0], (np.asarray(yb) - ya) * scale[1])
    d = d[np.isfinite(d)]
    return "no image" if d.size == 0 else f"{float(d.max()):.3g} m"


def short(x: str) -> str:
    return x if len(x) < 60 else x[:57] + "..."


def relation(rule: Callable[[str, str], bool], crss: list[str]) -> None:
    same = {(i, j) for i, x in enumerate(crss) for j, y in enumerate(crss) if i != j and rule(x, y)}
    asym = [p for p in same if (p[1], p[0]) not in same]
    breaks = [
        (i, j, k) for i, j in same for j2, k in same if j2 == j and k != i and (i, k) not in same
    ]
    print(f"{rule.__name__}: same pairs {len(same)}, asymmetric {len(asym)}, breaks {len(breaks)}")
    for i, j, k in breaks[:3]:
        print(f"    {short(crss[i])} ~ {short(crss[j])} ~ {short(crss[k])}")
    print(f"    not the same as itself: {[short(c) for c in crss if not rule(c, c)]}")


def tolerance() -> None:
    for code, key in (
        (25833, "Scale factor"),
        (25833, "False easting"),
        (31287, "Longitude of false origin"),
        (31287, "Latitude of 1st"),
    ):
        for rel in (1.1e-10, 2e-10, 1e-9):
            j = no_id(code)
            for p in j["conversion"]["parameters"]:
                if p["name"].startswith(key):
                    p["value"] *= 1 + rel
            perturbed, own = CRS.from_json_dict(j), CRS.from_epsg(code)
            assert own.area_of_use is not None
            w, s, e, n = own.area_of_use.bounds
            lon, lat = np.meshgrid(np.linspace(w, e, 15), np.linspace(s, n, 15))
            geo = own.geodetic_crs
            xa, ya = Transformer.from_crs(geo, perturbed, always_xy=True).transform(
                lon.ravel(), lat.ravel()
            )
            xb, yb = Transformer.from_crs(geo, own, always_xy=True).transform(
                lon.ravel(), lat.ravel()
            )
            mm = float(np.nanmax(np.hypot(np.subtract(xa, xb), np.subtract(ya, yb)))) * 1000
            same = rule_now(perturbed.to_wkt(), f"EPSG:{code}")
            print(f"EPSG:{code} {key} x(1+{rel:g}): same={same} moves {mm:.3f} mm")


def main() -> None:
    crss: set[str] = set()
    for component, pairs in ROWS:
        print(f"## {component}")
        for x, y in pairs:
            a, b = name(x), name(y)
            crss |= {a, b}
            old = (rule_5197f9a(a, b), rule_5197f9a(b, a))
            now = (rule_now(a, b), rule_now(b, a))
            print(f"  {short(a):60} vs {short(b):22} 5197f9a={old} now={now} moves={moves(a, b)}")
    ordered = sorted(crss)
    print(f"{len(ordered)} CRSs")
    relation(rule_now, ordered)
    relation(rule_5197f9a, ordered)
    tolerance()


if __name__ == "__main__":
    main()
