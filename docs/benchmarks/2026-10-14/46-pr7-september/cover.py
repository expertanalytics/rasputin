"""Increment 46 PR 7: snow cover of a `rasputin snow` EnSight output.

For each output directory: the frame count, the snow-covered fraction on the
frames asked for (by map area and by triangle count; one bin, so sca is 0 or
1), the domain mean SWE, and the melt-out day (`melt_out_day`, the first
frame's day with sca < 0.5, -1 if never): the last day a triangle melts out,
and how many triangles never do.

    python cover.py <out.json> <frame> [<frame> ...] -- <label>=<dir> [...]
"""

import json
import sys
from pathlib import Path

import numpy as np


def text(b: bytes, at: int) -> tuple[str, int]:
    return b[at:at + 80].rstrip(b"\0 ").decode(), at + 80


def geometry(d: Path) -> tuple[np.ndarray, np.ndarray]:
    b = (d / "snow.geo").read_bytes()
    at = 6 * 80 + 4 + 80 + 80
    n = int(np.frombuffer(b, "<i4", 1, at)[0]); at += 4
    xyz = np.frombuffer(b, "<f4", 3 * n, at).reshape(3, n).astype(np.float64); at += 12 * n
    assert text(b, at)[0] == "tria3"; at += 80
    m = int(np.frombuffer(b, "<i4", 1, at)[0]); at += 4
    tri = np.frombuffer(b, "<i4", 3 * m, at).reshape(m, 3) - 1
    return xyz, tri


def var(path: Path, m: int) -> np.ndarray:
    v = np.frombuffer(path.read_bytes(), "<f4", offset=3 * 80 + 4)
    assert len(v) == m
    return v.astype(np.float64)


def main() -> None:
    out, rest = sys.argv[1], sys.argv[2:]
    cut = rest.index("--")
    frames, dirs = [int(f) for f in rest[:cut]], dict(a.split("=", 1) for a in rest[cut + 1:])
    result = {}
    for label, d in dirs.items():
        d = Path(d)
        xyz, tri = geometry(d)
        x, y = xyz[0][tri], xyz[1][tri]
        area = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
        m = len(tri)
        n_frames = len(list(d.glob("swe[0-9]*.var")))
        r = {"triangles": m, "frames": n_frames, "map_area_km2": area.sum() / 1e6, "frames_asked": {}}
        for f in frames:
            sca = var(d / f"sca{f:04d}.var", m)
            swe = var(d / f"swe{f:04d}.var", m)
            r["frames_asked"][f] = {
                "day": f, "sca_values": sorted({float(v) for v in np.unique(sca)})[:5],
                "covered_by_area": float((sca * area).sum() / area.sum()),
                "covered_by_count": float(sca.mean()),
                "swe_mean_mm_by_area": float((swe * area).sum() / area.sum()),
                "swe_max_mm": float(swe.max()),
            }
        mo = var(d / "melt_out_day.var", m)
        melted = mo >= 0
        r["melt_out"] = {
            "never": int((~melted).sum()),
            "never_area_share": float(area[~melted].sum() / area.sum()),
            "last_day": float(mo[melted].max()) if melted.any() else None,
            "median_day_by_count": float(np.median(mo[melted])) if melted.any() else None,
        }
        result[label] = r
    if "on" in result and "off" in result:
        don, doff = Path(dirs["on"]), Path(dirs["off"])
        xyz, tri = geometry(don)
        m = len(tri)
        mo_on, mo_off = var(don / "melt_out_day.var", m), var(doff / "melt_out_day.var", m)
        result["on_minus_off"] = {f: {"swe_diff_range_mm": [float(v) for v in (
            (var(don / f"swe{f:04d}.var", m) - var(doff / f"swe{f:04d}.var", m)).min(),
            (var(don / f"swe{f:04d}.var", m) - var(doff / f"swe{f:04d}.var", m)).max())]}
            for f in frames}
        both = (mo_on >= 0) & (mo_off >= 0)
        result["on_minus_off"]["melt_out_later_with_shadows"] = int((mo_on[both] > mo_off[both]).sum())
        result["on_minus_off"]["melt_out_max_delay_days"] = float((mo_on[both] - mo_off[both]).max())
    Path(out).write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
