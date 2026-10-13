# ruff: noqa: E501 -- a one-off design probe; its formulas and report lines are long
"""Increment 46 design probe: the synthetic spring in NumPy, standing in for PR 1-3.

Run from the worktree root with its venv:
    .venv/bin/python docs/increments/46-probes/probe_spring.py WORKDIR > docs/increments/46-probes/probe_spring.out
WORKDIR receives synth.tif and mesh.vtk (not committed). The sun is NOAA's
solar calculator (noaa_sun.py), standing in for TR1. Everything printed is quoted by
docs/increments/46-hbv-snow.md, sections 4.5, 4.7, 6 and 7.
"""

import datetime as dt
import math
import pathlib
import statistics
import subprocess
import sys
import time

import numpy as np
import tifffile
from pyproj import CRS, Proj

work = pathlib.Path(sys.argv[1])
work.mkdir(parents=True, exist_ok=True)
rasputin = pathlib.Path(sys.executable).parent / "rasputin"

# ---- terrain (section 4.5) ----
n_side, step, x0, ytop = 1201, 10.0, 494000.0, 6806000.0
xs = (np.arange(n_side) * step) - 6000.0
ys = 6000.0 - (np.arange(n_side) * step)
X, Y = np.meshgrid(xs, ys)
r, phi = np.hypot(X, Y), np.arctan2(Y, X)
z = 300 + 1600 * np.exp(-((r / 2500) ** 2)) * (
    1 + 0.15 * (1 - np.exp(-((r / 1500) ** 2))) * np.cos(5 * phi)
)
z += 250 * np.exp(-((X - 3500) ** 2 + (Y + 3000) ** 2) / 900**2)
rng = np.random.default_rng(46)
k = np.hypot(*np.meshgrid(np.fft.fftfreq(n_side, step), np.fft.fftfreq(n_side, step)))
amp = np.zeros_like(k)
keep = (k > 0) & (k <= 1 / 40.0)
amp[keep] = k[keep] ** -1.5
noise = np.real(np.fft.ifft2(amp * np.exp(2j * np.pi * rng.random(k.shape))))
z += noise * (25.0 / noise.std())
gy, gx = np.gradient(z, step)
print(
    f"terrain z {z.min():.1f} .. {z.max():.1f} m; grid slope mean "
    f"{np.degrees(np.arctan(np.hypot(gx, gy))).mean():.1f} max {np.degrees(np.arctan(np.hypot(gx, gy))).max():.1f} deg"
)
geokeys = (1, 1, 0, 3, 1024, 0, 1, 1, 1025, 0, 1, 2, 3072, 0, 1, 25833)
tifffile.imwrite(
    work / "synth.tif",
    z.astype(np.float32),
    extratags=[
        (33550, "d", 3, (step, step, 0.0), False),
        (33922, "d", 6, (0, 0, 0, x0, ytop, 0), False),
        (34735, "H", len(geokeys), geokeys, False),
    ],
)

# ---- mesh (2 m tolerance) ----
t = time.perf_counter()
subprocess.run(
    [
        str(rasputin),
        "mesh",
        "--dem",
        str(work / "synth.tif"),
        "--tolerance",
        "2",
        "--out",
        str(work / "mesh.vtk"),
    ],
    check=True,
    capture_output=True,
)
data = (work / "mesh.vtk").read_bytes()
i = data.index(b"POINTS ")
line = data[i : data.index(b"\n", i)].split()
nv = int(line[1])
p0 = data.index(b"\n", i) + 1
pts = np.frombuffer(data, ">f8", nv * 3, p0).reshape(-1, 3).astype(float)
i = data.index(b"POLYGONS ", p0 + nv * 24)
line = data[i : data.index(b"\n", i)].split()
tri = (
    np.frombuffer(data, ">i4", int(line[2]), data.index(b"\n", i) + 1)
    .reshape(-1, 4)[:, 1:]
    .astype(np.int64)
)
print(
    f"mesh --tolerance 2: {len(tri)} triangles, {nv} vertices, {time.perf_counter() - t:.2f} s incl. write"
)

# ---- geometry ----
a, b, c = pts[tri[:, 0]], pts[tri[:, 1]], pts[tri[:, 2]]
nrm = np.cross(b - a, c - a)
nrm[nrm[:, 2] < 0] *= -1
area3 = 0.5 * np.linalg.norm(nrm, axis=1)
nrm /= np.linalg.norm(nrm, axis=1)[:, None]
gam = math.radians(Proj(CRS.from_epsg(25833)).get_factors(15.0, 61.334).meridian_convergence)
ce, se = math.cos(gam), math.sin(gam)  # grid -> true ENU (gamma = 0 on the central meridian)
nrm = np.column_stack([ce * nrm[:, 0] - se * nrm[:, 1], se * nrm[:, 0] + ce * nrm[:, 1], nrm[:, 2]])
zc = (a[:, 2] + b[:, 2] + c[:, 2]) / 3
slope = np.degrees(np.arccos(np.clip(nrm[:, 2], -1, 1)))
aspect = np.degrees(np.arctan2(nrm[:, 0], nrm[:, 1])) % 360
band = (zc >= 900) & (zc <= 1100) & (slope >= 25) & (slope <= 35)
south = band & (np.abs(aspect - 180) <= 22.5)
north = band & ((aspect <= 22.5) | (aspect >= 337.5))
print(
    f"gamma at centre {math.degrees(gam):.4f} deg; 900-1100 m, 25-35 deg: south {south.sum()}, north {north.sum()} triangles"
)

# ---- sun: NOAA's solar calculator (noaa_sun.py), standing in for TR1 ----
from noaa_sun import sun as noaa  # noqa: E402


def sun(tt):
    e_deg, a_deg, dist = noaa(tt, 61.334, 15.0)
    el, az = math.radians(e_deg), math.radians(a_deg)
    u = np.array([math.cos(el) * math.sin(az), math.cos(el) * math.cos(az), math.sin(el)])
    return u, (1368 / dist**2 * 0.75 ** (1 / math.sin(el)) if el > 0 else 0.0)


def plane(sl, asp):
    s, a_ = math.radians(sl), math.radians(asp)
    return np.array([math.sin(s) * math.sin(a_), math.sin(s) * math.cos(a_), math.cos(s)])


# R1: one sample, sun at 30 deg due south, beam 800
u30 = np.array([0, -math.cos(math.radians(30)), math.sin(math.radians(30))])
print(
    "R1 anomaly:",
    {
        nm: round(800 * max(0, plane(30, a_) @ u30) - 800 * u30[2], 1)
        for nm, a_ in (("S", 180), ("N", 0), ("E", 90))
    },
)
acc = dict.fromkeys("HSNEW", 0.0)
for q in range(96):
    u, B = sun(dt.datetime(2027, 4, 15) + dt.timedelta(minutes=15 * q + 7.5))
    acc["H"] += B * max(0, u[2])
    for nm, a_ in zip("SNEW", (180, 0, 90, 270), strict=False):
        acc[nm] += B * max(0, plane(30, a_) @ u)
print(
    "R2 day-mean beam 2027-04-15:",
    {k_: round(v / 96, 1) for k_, v in acc.items()},
    "anomaly:",
    {k_: round((acc[k_] - acc["H"]) / 96, 1) for k_ in "SNEW"},
)


# ---- D1 ----
def lognormal(cv, n):
    if cv == 0:
        return np.ones(1)
    s = math.sqrt(math.log(1 + cv * cv))
    nd = statistics.NormalDist()
    f = np.array([math.exp(s * nd.inv_cdf((q + 0.5) / n)) for q in range(n)])
    return f / f.mean()


print("D1 lognormal(0.4, 10):", np.round(lognormal(0.4, 10), 4).tolist())

# ---- the spring ----
TX, TS, CX, LW, CFR, RR, DD = 0.0, 0.0, 3.0, 0.1, 0.05, 0.015, 0.02
start, steps, dts = dt.datetime(2027, 4, 1), 91 * 24, 3600.0
runs = {
    "none": ("none", 0.4, 10),
    "melt_factor": ("melt_factor", 0.4, 10),
    "threshold": ("threshold", 0.4, 10),
    "none_one_bin": ("none", 0.0, 1),
    "melt_factor_one_bin": ("melt_factor", 0.0, 1),
}
state = {}
for name, (mode, cv, nb) in runs.items():
    f = lognormal(cv, nb)
    state[name] = dict(
        mode=mode,
        f=f,
        sp=np.tile(600.0 * f, (len(tri), 1)),
        sw=np.zeros((len(tri), nb)),
        out=np.zeros(len(tri)),
        sca_hist=[],
        swe_hist=[],
    )
for s_ in state.values():
    s_["sca_hist"].append(np.ones(len(tri)))
    s_["swe_hist"].append(np.full(len(tri), 600.0))
anomaly_sum = np.zeros(len(tri))
t_step = []
dtd = dts / 86400
for st in range(steps):
    t0 = time.perf_counter()
    tm = start + dt.timedelta(seconds=dts * (st + 0.5))
    hours = st + 0.5
    hsol = (tm.hour + tm.minute / 60 + 1.0) % 24  # mean solar time at 15 deg E
    T = (-3 + 17 * hours / (91 * 24) + 4 * math.cos(2 * math.pi * (hsol - 14) / 24)) - 0.006 * (
        zc - 300
    )
    anomaly = np.zeros(len(tri))
    for q in range(4):
        u, B = sun(start + dt.timedelta(seconds=dts * st + dts / 4 * (q + 0.5)))
        if B > 0:
            anomaly += B * (np.maximum(0, nrm @ u) - max(0, u[2]))
    anomaly /= 4
    anomaly_sum += anomaly
    for s_ in state.values():
        sp, sw, mode = s_["sp"], s_["sw"], s_["mode"]
        m = np.full(len(tri), CX)
        th = np.full(len(tri), TS)
        if mode == "melt_factor":
            m = np.maximum(0, CX + RR * anomaly)
        if mode == "threshold":
            th = TS - DD * anomaly
        P = 0.0  # dry spring: steps 1, 2 and 4 of section 4.2 are no-ops
        Tm, mm, thm = T[:, None], m[:, None], th[:, None]
        melt = np.where(Tm > thm, np.minimum(mm * dtd * (Tm - thm), sp), 0.0)
        refr = np.where((Tm <= thm) & (Tm < TS), np.minimum(CFR * CX * dtd * (TS - Tm), sw), 0.0)
        sp += refr - melt
        sw += melt - refr
        o = np.maximum(0, sw - LW * sp)
        sw -= o
        s_["out"] += o.mean(axis=1)
    if (st + 1) % 24 == 0:
        for s_ in state.values():
            s_["sca_hist"].append((s_["sp"] > 0).mean(axis=1))
            s_["swe_hist"].append((s_["sp"] + s_["sw"]).mean(axis=1))
    t_step.append(time.perf_counter() - t0)

print(
    f"cost: {np.median(t_step) * 1e3:.1f} ms median per step, all five runs plus dI, {sum(t_step):.1f} s total"
)
w = area3 / area3.sum()
for name, s_ in state.items():
    sca, swe = np.array(s_["sca_hist"]), np.array(s_["swe_hist"])
    below = sca < 0.5
    mo = np.where(below.any(axis=0), below.argmax(axis=0), -1)
    first = np.where((sca < 1).any(axis=0), (sca < 1).argmax(axis=0), -1)
    bare = np.where((sca == 0).any(axis=0), (sca == 0).argmax(axis=0), -1)
    both = (first >= 0) & (bare >= 0)
    mass = abs(600.0 - swe[-1] - s_["out"]).max()
    print(
        f"{name}: domain SCA d30 {w @ sca[30]:.4f} d60 {w @ sca[60]:.4f} d90 {w @ sca[90]:.4f}; "
        f"SWE d30 {w @ swe[30]:.2f} d60 {w @ swe[60]:.2f} d90 {w @ swe[90]:.2f} mm; "
        f"never below 0.5: {(mo < 0).mean():.4f} of triangles"
    )
    for cls, sel in (("south", south), ("north", north)):
        ws = area3[sel] / area3[sel].sum()
        csca = sca[:, sel] @ ws
        half = int(np.argmax(csca < 0.5)) if (csca < 0.5).any() else -1
        print(
            f"   {cls} 30deg 900-1100 m: median melt-out day {np.median(mo[sel]):.1f}, "
            f"class SCA d60 {csca[60]:.3f}, class SCA < 0.5 first at day {half}, "
            f"median days first-bare-bin..bare {np.median(bare[sel & both] - first[sel & both]):.1f}"
        )
    print(f"   I1 over the run: max |600 - swe_end - outflow| = {mass:.3e} mm")
print(
    f"radiation_anomaly_mean: S class {(anomaly_sum / steps)[south].mean():.1f}, N class {(anomaly_sum / steps)[north].mean():.1f} W/m2"
)
