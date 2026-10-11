"""Increment 34 design probe: section 3's slope at the grid border and next to
NoData (docs/increments/34-probes/README.md, item 1b; design review round 1, B1).

On planes, Horn's slope by round 1's rule (a missing neighbour takes the node's
own z) against section 3's rule (its reflection through the node, 2 z(node) -
z(opposite), when the opposite is there). Then the two rules on the Romsdalen
window: where they differ.
"""

import numpy as np
from slope_stats import CASES, horn, horn_edge, window


def plane(rows, cols, dx, dy, deg, azimuth_deg):
    """z rising at `deg` towards `azimuth_deg` (0: +x, 90: -row, i.e. north)."""
    r, c = np.indices((rows, cols), dtype=np.float64)
    x, y = c * dx, -r * dy
    a = np.radians(azimuth_deg)
    return np.tan(np.radians(deg)) * (x * np.cos(a) + y * np.sin(a))


def round1(z, dx, dy, valid):
    """Round 1's rule: every missing neighbour (outside, or NoData) takes z(node)."""
    p = np.pad(np.where(valid, z, np.nan), 1, mode="constant", constant_values=np.nan)
    rows, cols = z.shape

    def nb(di, dj):
        v = p[1 + di : 1 + di + rows, 1 + dj : 1 + dj + cols]
        return np.where(np.isnan(v), z, v)

    a, b, c, d, f = nb(-1, -1), nb(-1, 0), nb(-1, 1), nb(0, -1), nb(0, 1)
    g, h, i = nb(1, -1), nb(1, 0), nb(1, 1)
    gx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8 * dx)
    gy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8 * dy)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def classes(s):
    return np.ceil(2.0 * s)


def main():
    for deg, az, dx, dy in (
        (40.0, 0.0, 10.0, 10.0),
        (40.0, 30.0, 10.0, 10.0),
        (37.3, 30.0, 10.0, 5.0),
    ):
        z = plane(21, 21, dx, dy, deg, az)
        valid = np.ones(z.shape, bool)
        valid[10, 10] = False  # one NoData node in the middle
        e, h = round1(z, dx, dy, valid), horn(z, dx, dy, valid)
        print(f"plane {deg} deg towards {az} deg, {dx:g} x {dy:g} m:")
        for name, (i, j) in {
            "interior": (5, 5),
            "border": (0, 5),
            "corner": (0, 0),
            "next to NoData": (10, 11),
            "corner next to NoData": (9, 9),
        }.items():
            print(
                f"  {name:22s} round 1 {e[i, j]:7.3f}  section 3 {h[i, j]:7.3f}  "
                f"class {classes(h[i, j]):.0f}"
            )
        ok = valid.copy()
        print(
            f"  section 3, every valid node: {h[ok].min():.9f} .. {h[ok].max():.9f} deg, "
            f"classes {sorted({int(v) for v in classes(h[ok])})}"
        )
    # Where reflection cannot help: both neighbours along an axis missing.
    strip = plane(1, 21, 10.0, 10.0, 40.0, 30.0)
    print(
        f"one-row grid, plane 40 deg towards 30 deg: {horn(strip, 10.0)[0, 5]:.3f} deg "
        "(the north-south part is lost)"
    )
    # NoData above and below a node, its four corner neighbours valid: the edge
    # neighbours fall back to z(node), so only the corners carry that axis.
    z = plane(21, 21, 10.0, 10.0, 40.0, 90.0)
    valid = np.ones(z.shape, bool)
    valid[9, 10] = valid[11, 10] = False
    print(
        f"NoData above and below, corners valid, plane 40 deg facing north: "
        f"{horn(z, 10.0, 10.0, valid)[10, 10]:.3f} deg (the north-south part at half weight)"
    )
    z = window(*CASES["romsdalen"])
    a, b = classes(horn_edge(z)), classes(horn(z))
    ring = np.zeros(z.shape, bool)
    ring[[0, -1], :] = ring[:, [0, -1]] = True
    diff = a != b
    print(
        f"romsdalen: classes differ at {diff.sum()} nodes, all on the outer ring: "
        f"{bool((diff & ~ring).sum() == 0)}; ring {ring.sum()} nodes; "
        f"at 30 deg or more on the ring: round 1 probe (numpy edge padding) "
        f"{np.mean(a[ring] >= 60) * 100:.1f}%, section 3 {np.mean(b[ring] >= 60) * 100:.1f}%, "
        f"interior {np.mean(b[~ring] >= 60) * 100:.1f}%"
    )


if __name__ == "__main__":
    main()
