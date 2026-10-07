"""Increment 30b, section 4.3: when does GEOS's intersection of a line with a
polygon that covers it return the line unchanged?

`docs/increments/30b-clip-speed.md`. Prints, per hand-made case, whether the
domain covers the line, whether `shapely.intersection(line, domain)` returns
the line itself (same WKB, one piece), and the pieces it returns otherwise.
With a `tin_engine` that has `feature_input._untouched` (30b's green tree), it
also prints that helper's answer, which must be True exactly where the line
comes back unchanged. Run with any interpreter that has shapely:

    python docs/increments/30b-probes/geos_cases.py
"""

from __future__ import annotations

import shapely
from shapely.geometry import LineString, Polygon

try:
    from tin_engine.feature_input import _untouched  # type: ignore[attr-defined]
except ImportError:
    _untouched = None

X0, Y0 = 500_000.0, 6_600_000.0


def at(x: float, y: float) -> tuple[float, float]:
    return (X0 + x, Y0 + y)


# A square with an extra vertex at (100, 50) on its east side, and a holed square.
DOMAIN = Polygon([at(0, 0), at(100, 0), at(100, 50), at(100, 100), at(0, 100)])
HOLED = Polygon(
    [at(0, 0), at(100, 0), at(100, 100), at(0, 100)],
    [[at(40, 40), at(60, 40), at(60, 60), at(40, 60)]],
)
CASES = {
    "plain inside": (DOMAIN, [at(10, 10), at(20, 30), at(50, 20)]),
    "repeated vertex": (DOMAIN, [at(10, 10), at(20, 30), at(20, 30), at(50, 20)]),
    "closed ring inside": (DOMAIN, [at(10, 10), at(30, 10), at(30, 30), at(10, 30), at(10, 10)]),
    "closed ring touching at its start": (
        DOMAIN, [at(0, 10), at(30, 10), at(30, 30), at(0, 30), at(0, 10)]),
    "closed ring touching mid": (
        DOMAIN, [at(10, 10), at(30, 10), at(30, 30), at(0, 20), at(10, 10)]),
    "domain vertex on the line": (DOMAIN, [at(50, 50), at(100, 40), at(100, 60), at(60, 70)]),
    "along the boundary": (DOMAIN, [at(100, 10), at(100, 90)]),
    "along, then inside": (DOMAIN, [at(100, 10), at(100, 90), at(50, 90)]),
    "self-crossing": (DOMAIN, [at(10, 10), at(50, 50), at(50, 10), at(10, 50)]),
    "self-touching at a vertex": (
        DOMAIN, [at(10, 10), at(50, 10), at(30, 30), at(10, 10), at(10, 50)]),
    "doubling back": (DOMAIN, [at(10, 10), at(50, 10), at(30, 10)]),
    "collinear middle vertex": (DOMAIN, [at(10, 10), at(20, 10), at(30, 10)]),
    "touching a hole": (HOLED, [at(10, 10), at(40, 50), at(10, 90)]),
    "zero length": (DOMAIN, [at(10, 10), at(10, 10)]),
}  # fmt: skip


def main() -> None:
    print(f"shapely {shapely.__version__}, GEOS {shapely.geos_version_string}")
    for name, (domain, coords) in CASES.items():
        line = LineString(coords)
        parts = [p for p in shapely.get_parts(shapely.intersection(line, domain))
                 if isinstance(p, LineString) and p.length > 0]  # fmt: skip
        same = len(parts) == 1 and shapely.to_wkb(parts[0]) == shapely.to_wkb(line)
        shapely.prepare(domain)
        helper = "" if _untouched is None else f" untouched={_untouched(line, domain)!s:5}"
        print(f"{name:34s} covers={domain.covers(line)!s:5} unchanged={same!s:5} "
              f"pieces={len(parts)}{helper}")  # fmt: skip


if __name__ == "__main__":
    main()
