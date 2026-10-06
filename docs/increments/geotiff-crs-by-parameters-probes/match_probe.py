"""Probe for `docs/increments/geotiff-crs-by-parameters.md`: the matching rule on
hand-built cases, and the sweep over EPSG's own codes.

Run from the repository root with the project venv:

    .venv/bin/python docs/increments/geotiff-crs-by-parameters-probes/match_probe.py [--sweep]

It needs pyproj only, and loads `src_python/tin_engine/crs.py` by path for
`same_crs` (B's rule). The CRS is built the way the design says the reader
builds it: PROJJSON, the base geographic CRS from its EPSG code, the
parameters in the EPSG method's own order, axes east then north in metres.
"""

from __future__ import annotations

import collections
import importlib.util
import sys
import time
from pathlib import Path
from typing import Any

from pyproj import CRS, Transformer
from pyproj.database import query_crs_info
from pyproj.enums import PJType
from pyproj.exceptions import CRSError, ProjError

ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location("crs", ROOT / "src_python/tin_engine/crs.py")
assert _spec is not None and _spec.loader is not None
crs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(crs)

LCC2 = {"name": "Lambert Conic Conformal (2SP)", "id": {"authority": "EPSG", "code": 9802}}
TM = {"name": "Transverse Mercator", "id": {"authority": "EPSG", "code": 9807}}
AXES = [
    {"name": "Easting", "abbreviation": "E", "direction": "east", "unit": "metre"},
    {"name": "Northing", "abbreviation": "N", "direction": "north", "unit": "metre"},
]


def p(name: str, code: int, value: float, unit: str) -> dict[str, Any]:
    return {"name": name, "value": value, "unit": unit, "id": {"authority": "EPSG", "code": code}}


def build(geog: int | dict[str, Any], method: dict[str, Any], params: list[dict[str, Any]]) -> CRS:
    base = CRS.from_epsg(geog).to_json_dict() if isinstance(geog, int) else geog
    return CRS.from_json_dict(
        {
            "type": "ProjectedCRS",
            "name": "unknown",
            "base_crs": base,
            "conversion": {"name": "unknown", "method": method, "parameters": params},
            "coordinate_system": {"subtype": "Cartesian", "axis": AXES},
        }
    )


def equivalent_xy(a: CRS, b: CRS) -> bool:
    """Leg 1 of `same_crs`: PROJ equivalence on the transformer's x-then-y copies."""
    try:
        t = Transformer.from_crs(a, b, always_xy=True)
    except ProjError:
        return False
    s, d = t.source_crs, t.target_crs
    return s is not None and d is not None and s.equals(d)


def matches(built: CRS) -> list[int]:
    """The design's rule: PROJ's EPSG candidates at every confidence, kept when on
    the file's own base geographic CRS and equivalent once x-then-y; lowest first."""
    base = built.geodetic_crs.to_epsg(min_confidence=100) if built.geodetic_crs else None
    found = []
    for m in built.list_authority(auth_name="EPSG", min_confidence=0):
        other = CRS.from_epsg(int(m.code))
        on_base = (
            other.geodetic_crs is not None
            and other.geodetic_crs.to_epsg(min_confidence=100) == base
        )
        if on_base and equivalent_xy(built, other):
            found.append(int(m.code))
    return sorted(found)


def lcc(geog: int, lat0: float, lon0: float, p1: float, p2: float, fe: float, fn: float) -> CRS:
    return build(
        geog,
        LCC2,
        [
            p("Latitude of false origin", 8821, lat0, "degree"),
            p("Longitude of false origin", 8822, lon0, "degree"),
            p("Latitude of 1st standard parallel", 8823, p1, "degree"),
            p("Latitude of 2nd standard parallel", 8824, p2, "degree"),
            p("Easting at false origin", 8826, fe, "metre"),
            p("Northing at false origin", 8827, fn, "metre"),
        ],
    )


def tm(
    geog: int,
    lon0: float,
    k: float = 0.9996,
    fe: float = 500000.0,
    fn: float = 0.0,
    lat0: float = 0.0,
) -> CRS:
    return build(
        geog,
        TM,
        [
            p("Latitude of natural origin", 8801, lat0, "degree"),
            p("Longitude of natural origin", 8802, lon0, "degree"),
            p("Scale factor at natural origin", 8805, k, "unity"),
            p("False easting", 8806, fe, "metre"),
            p("False northing", 8807, fn, "metre"),
        ],
    )


AUSTRIA_LON = 13.33333333300013  # the file's ProjFalseOriginLongGeoKey (3084)


def cases() -> None:
    table = {
        "Austrian file (MGI 4312)": lcc(4312, 47.5, AUSTRIA_LON, 46.0, 49.0, 400000.0, 400000.0),
        "Austrian file, parallels 49 then 46": lcc(
            4312, 47.5, AUSTRIA_LON, 49.0, 46.0, 400000.0, 400000.0
        ),
        "Austrian, false easting +1e-5 m": lcc(
            4312, 47.5, AUSTRIA_LON, 46.0, 49.0, 400000.00001, 400000.0
        ),
        "Austrian, false easting +1 mm": lcc(
            4312, 47.5, AUSTRIA_LON, 46.0, 49.0, 400000.001, 400000.0
        ),
        "Austrian on ETRS89 4258": lcc(4258, 47.5, AUSTRIA_LON, 46.0, 49.0, 400000.0, 400000.0),
        "Austrian on WGS 84 4326": lcc(4326, 47.5, AUSTRIA_LON, 46.0, 49.0, 400000.0, 400000.0),
        "TM 15E on ETRS89": tm(4258, 15.0),
        "TM 15E on WGS 84": tm(4326, 15.0),
        "TM 15E k=1 on ETRS89": tm(4258, 15.0, k=1.0),
        "MGI GK M31": tm(4312, 13.333333333333334, k=1.0, fe=450000.0, fn=-5000000.0),
        "TM 27E on ETRS89 (UTM 35N)": tm(4258, 27.0),
        "Xian 1980 GK CM 75E": tm(4610, 75.0, k=1.0),
    }
    for name, built in table.items():
        t0 = time.perf_counter()
        found = matches(built)
        ms = (time.perf_counter() - t0) * 1000
        same = (
            all(crs.same_crs(f"EPSG:{found[0]}", f"EPSG:{n}") for n in found[1:]) if found else None
        )
        offered = len(built.list_authority(auth_name="EPSG", min_confidence=0))
        print(f"{name:40s} PROJ offers {offered}, matches={found} mutually same={same} {ms:.0f} ms")
    built = table["Austrian file (MGI 4312)"]
    print("B's same_crs(built Austrian, EPSG:31287):", crs.same_crs(built, "EPSG:31287"))
    print(
        "pyproj to_epsg() at its default 70:",
        built.to_epsg(),
        "; PROJ's offer:",
        [(m.code, m.confidence) for m in built.list_authority(auth_name="EPSG", min_confidence=0)],
    )


def sweep() -> None:
    """Every non-deprecated EPSG projected CRS by 9807 or 9802 on a geographic 2D
    base with an EPSG code, Greenwich, parameters in degrees, metres or unity."""
    units = {"degree", "metre", "unity"}
    counts: collections.Counter[str] = collections.Counter()
    groups: list[list[int]] = []
    t0 = time.time()
    for info in query_crs_info(
        auth_name="EPSG", pj_types=PJType.PROJECTED_CRS, allow_deprecated=False
    ):
        try:
            c = CRS.from_epsg(int(info.code))
        except CRSError:
            continue
        d = c.to_json_dict()
        conv = d.get("conversion", {})
        method = conv.get("method", {})
        if method.get("id", {}).get("code") not in (9807, 9802):
            continue
        base = d["base_crs"]
        if base.get("id") is None or c.geodetic_crs is None or len(c.geodetic_crs.axis_info) != 2:
            continue
        unit_names = {
            q["unit"] if isinstance(q["unit"], str) else q["unit"].get("name")
            for q in conv["parameters"]
        }
        if not unit_names <= units or any(a.unit_name != "metre" for a in c.axis_info):
            counts["skipped: units"] += 1
            continue
        if c.prime_meridian is None or c.prime_meridian.longitude != 0:
            counts["skipped: prime meridian"] += 1
            continue
        params = [{k: q[k] for k in ("name", "value", "unit", "id")} for q in conv["parameters"]]
        found = matches(build(base, method, params))
        code = int(info.code)
        if not found:
            counts["matches nothing"] += 1
        elif found == [code]:
            counts["itself only"] += 1
        elif len(found) == 1:
            counts[
                "one other code, same_crs to it"
                if crs.same_crs(f"EPSG:{code}", f"EPSG:{found[0]}")
                else "one other code, NOT same"
            ] += 1
        else:
            groups.append(found)
            counts["several codes"] += 1
    pairs = collections.Counter(
        crs.same_crs(f"EPSG:{a}", f"EPSG:{b}")
        for g in groups
        for i, a in enumerate(g)
        for b in g[i + 1 :]
    )
    print(dict(counts), f"{time.time() - t0:.0f} s")
    print(
        "pairs within several-code groups, same_crs:",
        dict(pairs),
        "; group sizes:",
        dict(collections.Counter(len(g) for g in groups)),
        "; e.g.",
        groups[:3],
    )


if __name__ == "__main__":
    cases()
    if "--sweep" in sys.argv:
        sweep()
