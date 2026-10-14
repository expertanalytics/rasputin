# ruff: noqa: E501 -- a one-off design probe; its formulas and report lines are long
"""Increment 46 addendum, design review round 3: the figures section 10 gained.

Run from the worktree root with the main checkout's .venv:
    python docs/increments/46-probes/probe_review3.py > docs/increments/46-probes/probe_review3.out

1. W: the direct beam at the sun elevations where the far field flipped a lit
   state (probe_shadow.out, B3), with probe_shadow.py's beam formula.
2. F: the bound on what terrain beyond 60 km can do to Bygdin's horizon.
3. M: PR 5's mutation targets, each with the fixture that kills it, on a
   one-dimensional model of the radiation design's pruned horizon search (its
   R4): a ray due east along a grid row, so bilinear sampling is linear
   interpolation along the row; a max-pyramid over 2^L-node blocks; the bound
   a/d_lo - c*d_lo for a >= 0 and a/d_hi - c*d_lo for a < 0; the ray's early
   stop with the larger of both rasters' maxima over the rest of the ray;
   reach found first, from coverage alone, and the march ended there. The
   exhaustive loop is the reference. A mutant is killed when its int16 table
   entry or its reach differs from the reference's. The model fixes no more of
   the kernel than R4 does; @tester's mutation run on the real kernel is the
   check that counts.
"""

import datetime as dt
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from noaa_sun import sun as noaa

C_EARTH = (1 - 0.13) / (2 * 6371000.0)

# ---- 1. the beam at the flip elevations ----
for day, e_flip in ((dt.datetime(2027, 4, 15, 12), 1.02), (dt.datetime(2027, 5, 15, 12), 1.90), (dt.datetime(2027, 6, 21, 12), 2.97)):
    dist = noaa(day, 61.3694, 8.4852)[2]
    b = 1368 / dist**2 * 0.75 ** (1 / math.sin(math.radians(e_flip)))
    print(f"W {day:%Y-%m-%d}: R {dist:.5f} AU; beam at e {e_flip:.2f} deg {b:.2f} W/m2 (at 1 AU {1368 * 0.75 ** (1 / math.sin(math.radians(e_flip))):.2f})")
for e in (1.0, 1.116, 2.0, 2.97, 3.0):
    print(f"W beam at 1 AU, e {e:.2f} deg: {1368 * 0.75 ** (1 / math.sin(math.radians(e))):.3f} W/m2")

# ---- 2. the bound beyond 60 km ----
# tan h <= (z_max - z_obs)/d - c d for any terrain at ground distance d; it falls with d, so d = 60 km bounds all d >= 60 km.
for z_max, what in ((2469.0, "Galdhopiggen, Norway's highest summit"), (2468.0, "B0's near-raster maximum")):
    for z_obs in (1054.0, 1500.0, 2240.0):
        t = (z_max - z_obs) / 60000.0 - C_EARTH * 60000.0
        print(f"F z_max {z_max:.0f} m ({what}), receiver at {z_obs:.0f} m: horizon from beyond 60 km at most {math.degrees(math.atan(t)):.3f} deg")


# ---- 3. the pruned search, modelled in one dimension, and its mutants ----
def at(prof, h, x):
    """Linear interpolation along the row; None off the profile or on a NaN node (a missing sample)."""
    i = math.floor(x / h)
    f = x / h - i
    if i < 0 or i >= len(prof) or (f > 0 and i + 1 >= len(prof)):
        return None
    a = prof[i]
    if f == 0:
        return None if math.isnan(a) else a
    b = prof[i + 1]
    if math.isnan(a) or math.isnan(b):
        return None
    return a * (1 - f) + b * f


def ray(fx, mut):
    """The ray's samples as (d, raster) pairs: near from r_i at half a near cell to D_near, far after it to D_far."""
    r_i, dn, df = fx["r_i"], fx["D_near"], fx["D_far"]
    if fx["far"] is None and mut != "none_is_zeros":
        df = dn  # no far raster: D_far = D_near
    start = 0.0 if mut == "start_at_0" else r_i + (fx["hn"] / 2 if mut == "start_late" else 0.0)
    out, j = [], 0
    near_end = df if mut == "near_past_D_near" else dn
    while start + j * fx["hn"] / 2 <= near_end + 1e-9:
        out.append((start + j * fx["hn"] / 2, "n"))
        j += 1
    if mut == "far_from_r_i":
        out = [(d, "n") for d, _ in out if d <= 0] + [(start + k * fx["hf"] / 2, "f") for k in range(int((df - start) / (fx["hf"] / 2)) + 1)]
        return out
    j = 1
    while dn + j * fx["hf"] / 2 <= df + 1e-9:
        out.append((dn + j * fx["hf"] / 2, "f"))
        j += 1
    return out


def prof_of(fx, rk, mut):
    if rk == "f":
        if fx["far"] is None:
            return np.zeros(4096), 1000.0  # only reached under the none_is_zeros mutant
        return fx["far"], fx["hf"]
    return fx["near"], fx["hn"]


def reach_of(fx, samples, mut):
    """The distance searched before the first missing sample, uint8 km, floor; 255 when none is missing."""
    for d, rk in samples:
        p, h = prof_of(fx, rk, mut)
        if at(p, h, d) is None:
            if mut == "reach_255":
                return 255, d
            return (round(d / 1000) if mut == "reach_round" else math.floor(d / 1000)), d
    return 255, math.inf


def value(fx, d, z):
    return (z - fx["z_obs"]) / d - fx["c"] * d


def exhaustive(fx, mut=None):
    s = ray(fx, mut)
    reach, d_miss = reach_of(fx, s, mut)
    best = -math.inf
    for d, rk in s:
        if d >= d_miss and mut != "past_nan":
            break
        if d <= 0:
            continue
        p, h = prof_of(fx, rk, mut)
        z = at(p, h, d)
        if z is not None:
            best = max(best, value(fx, d, z))
    return best, reach


def pyramid(p):
    lv = [np.where(np.isnan(p), -np.inf, p)]
    while len(lv[-1]) > 1:
        q = lv[-1]
        lv.append(np.array([q[2 * i : 2 * i + 2].max() for i in range((len(q) + 1) // 2)]))
    return lv


def blockmax(p, lv, lev, b):
    """Max over nodes b*2^L .. (b+1)*2^L inclusive (a sample in the block reads at most those)."""
    hi = min(len(p) - 1, (b + 1) * 2**lev)
    m = lv[lev][b] if b < len(lv[lev]) else -np.inf
    return max(m, -np.inf if np.isnan(p[hi]) else p[hi])


def bound(fx, a, d_lo, d_hi, mut):
    c = fx["c"]
    if a >= 0:
        return a / (d_hi if mut == "pos_branch_d_hi" else d_lo) - c * d_lo
    return a / (d_lo if mut == "neg_branch_d_lo" else d_hi) - c * (d_hi if mut == "curv_d_hi" else d_lo)


def bound_pos(fx, a, d_lo, d_hi, mut):
    if mut == "curv_d_hi" and a >= 0:
        return a / d_lo - fx["c"] * d_hi
    return bound(fx, a, d_lo, d_hi, mut)


def pruned(fx, mut=None):
    s = ray(fx, mut)
    reach, d_miss = reach_of(fx, s, mut)
    if mut != "past_nan":
        s = [(d, rk) for d, rk in s if d < d_miss]
    pyr = {rk: pyramid(prof_of(fx, rk, mut)[0]) for rk in ("n", "f")}
    gmax = {rk: pyr[rk][-1][0] for rk in pyr}
    best, j = -math.inf, 0
    while j < len(s):
        d, rk = s[j]
        if d <= 0:
            j += 1
            continue
        rest = {r for _, r in s[j:]}
        zr = gmax[rk] if mut == "stop_current_raster" else max(gmax[r] for r in rest)
        if bound_pos(fx, zr - fx["z_obs"], d, s[-1][0], mut) <= best:
            break  # the ray's early stop
        p, h = prof_of(fx, rk, mut)
        lv = pyr[rk]
        i0 = math.floor(d / h)
        took = False
        for lev in range(len(lv) - 1, 0, -1):
            b = i0 // 2**lev
            k = j
            while k + 1 < len(s) and s[k + 1][1] == rk and s[k + 1][0] / h <= (b + 1) * 2**lev:
                k += 1
            zmax = blockmax(p, lv, lev, b)
            if bound_pos(fx, zmax - fx["z_obs"], d, s[k][0], mut) <= best:
                j = k + 1
                took = True
                break
        if took:
            continue
        z = at(p, h, d)
        if z is not None:
            best = max(best, value(fx, d, z))
        j += 1
    return best, reach


def q(t):
    """int16 in 0.01 deg, rounded half away from zero; -90 deg with no sample."""
    if t == -math.inf:
        return -9000
    hdeg = 100 * math.degrees(math.atan(t))
    return int(math.copysign(math.floor(abs(hdeg) + 0.5), hdeg))


def flat(n, z=0.0):
    return np.full(n, z)


def tent(p, i, z):
    p = p.copy()
    p[i] = z
    return p


FX = {
    # receiver above all terrain: a < 0 everywhere; the near ridge sets a best that the far plain beats only far out
    "F1 above all terrain": dict(near=tent(flat(2101), 1000, 2000.0), hn=10.0, far=flat(2001, 1000.0), hf=100.0, z_obs=3000.0, r_i=10.0, D_near=20000.0, D_far=200000.0, c=C_EARTH),
    # a lone far peak just past D_near that beats the near best by less than the curvature over the far raster
    "F2 far peak at 25 km": dict(near=flat(2101), hn=10.0, far=tent(flat(2001), 250, 100.0), hf=100.0, z_obs=0.0, r_i=10.0, D_near=20000.0, D_far=200000.0, c=C_EARTH),
    # a spike on the node at 10 m; r_i = 7.5 m is not a multiple of half a cell
    "F4 spike at r_i": dict(near=tent(flat(401), 1, 100.0), hn=10.0, far=None, hf=100.0, z_obs=20.0, r_i=7.5, D_near=2000.0, D_far=2000.0, c=0.0),
    # near and far rasters disagree on purpose: a near-only ridge past D_near, a far-only ridge inside it
    "F5 near-far switch": dict(near=tent(flat(401), 300, 1000.0), hn=10.0, far=tent(tent(flat(61), 15, 1000.0), 40, 400.0), hf=100.0, z_obs=0.0, r_i=10.0, D_near=2000.0, D_far=5000.0, c=0.0),
    # a NaN node at 3.7 km and a ridge behind it
    "F6 coverage gap at 3.7 km": dict(near=tent(tent(flat(2101), 370, math.nan), 500, 1000.0), hn=10.0, far=None, hf=100.0, z_obs=0.0, r_i=10.0, D_near=20000.0, D_far=20000.0, c=0.0),
    # far raster None, D_far requested 60 km; receiver 100 m above a flat near raster
    "F7 far raster None": dict(near=flat(2101), hn=10.0, far=None, hf=100.0, z_obs=100.0, r_i=10.0, D_near=20000.0, D_far=60000.0, c=C_EARTH),
}
MUT = {
    "neg_branch_d_lo": "the a < 0 bound with d_lo for d_hi",
    "pos_branch_d_hi": "the a >= 0 bound with d_hi for d_lo",
    "curv_d_hi": "the curvature term at d_hi for d_lo",
    "stop_current_raster": "the early stop with the current raster's maximum only",
    "start_late": "the search starting at r_i + half a cell",
    "start_at_0": "the search starting at 0 (d_j = j * half a cell)",
    "near_past_D_near": "the near raster read past D_near",
    "far_from_r_i": "the far raster read from r_i",
    "past_nan": "the march not ended at the first missing sample",
    "reach_round": "reach rounded, not floored",
    "reach_255": "reach not recorded (255)",
    "none_is_zeros": "a far raster of None read as zeros to D_far",
}
for name, fx in FX.items():
    ref = exhaustive(fx)
    pr = pruned(fx)
    print(f"{name}: exhaustive h {math.degrees(math.atan(ref[0])) if ref[0] > -math.inf else -90:.4f} deg, int16 {q(ref[0])}, reach {ref[1]}; pruned int16 {q(pr[0])}, reach {pr[1]}; equal {(q(ref[0]), ref[1]) == (q(pr[0]), pr[1])}")
for m, what in MUT.items():
    kills = []
    for name, fx in FX.items():
        ref = exhaustive(fx)
        got = pruned(fx, m)
        if (q(got[0]), got[1]) != (q(ref[0]), ref[1]):
            kills.append(f"{name.split()[0]} ({q(got[0])}, {got[1]})")
    print(f"M {m} ({what}): killed by {', '.join(kills) if kills else 'NONE'}")
print(f"M quantisation: V-valley h(45) 22.2077 deg -> round {q(math.tan(math.radians(22.2077)))}, truncated {int(2220.77)}; F7 -> round {q(exhaustive(FX['F7 far raster None'])[0])}, floor {math.floor(100 * math.degrees(math.atan(exhaustive(FX['F7 far raster None'])[0])))}")

# ---- 4. the flat plane's sunrise at Bygdin, NOAA standing in for TR1 (PR 6's test S) ----
for day in (dt.datetime(2027, 4, 15), dt.datetime(2027, 6, 21)):
    h_flat = -math.degrees(math.atan(C_EARTH * 10.0))
    target = -0.8333 + h_flat  # true centre elevation at which the apparent upper limb meets the flat horizon (standard refraction 0.5667)
    lo, hi = day, day + dt.timedelta(hours=12)
    for _ in range(60):
        mid = lo + (hi - lo) / 2
        if noaa(mid, 61.3694, 8.4852)[0] < target:
            lo = mid
        else:
            hi = mid
    rate = (noaa(hi + dt.timedelta(seconds=30), 61.3694, 8.4852)[0] - noaa(hi - dt.timedelta(seconds=30), 61.3694, 8.4852)[0]) / 60.0
    print(f"S {day:%Y-%m-%d}: flat-plane horizon at r_i 10 m {h_flat:.7f} deg; sunrise (upper limb, standard refraction) {hi:%H:%M:%S} UTC, azimuth {noaa(hi, 61.3694, 8.4852)[1]:.2f} deg; de/dt {rate * 60:.4f} deg/min, so 0.003 deg is {0.003 / rate:.1f} s")
