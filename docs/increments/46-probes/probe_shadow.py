# ruff: noqa: E501 -- a one-off design probe; its formulas and report lines are long
"""Increment 46 addendum probe (section 10): cast shadows by the radiation design's horizon.

Run from the worktree root with a venv that has numpy, tifffile, pyproj and
tin_engine (the main checkout's .venv was used):
    python docs/increments/46-probes/probe_shadow.py DTM10_DIR > docs/increments/46-probes/probe_shadow.out
The horizon is the radiation design's reference definition (its R4: true
azimuth A_k = k*360/N, samples from the exclusion radius r_i at half a cell,
tan h = max_j (z(d_j) - z_obs - c d_j^2) / d_j, bilinear, a sample off the
raster or on NoData is no obstacle), written exhaustively in NumPy. The lit
test is its R5's without refraction (e > h(A) interpolated linearly between
directions, and n.s > 0). The sun is NOAA's (noaa_sun.py), standing in for TR1.
Everything printed is quoted by docs/increments/46-hbv-snow.md, section 10.
"""

import datetime as dt
import math
import pathlib
import resource
import sys
import time

import numpy as np
from pyproj import CRS, Proj, Transformer

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from noaa_sun import sun as noaa

N = 64
C_EARTH = (1 - 0.13) / (2 * 6371000.0)


def bilinear(z, x0, ytop, step, x, y):
    """z at (x, y); NaN off the raster. Node (0, 0) is at (x0, ytop)."""
    col, row = (x - x0) / step, (ytop - y) / step
    c0, r0 = np.floor(col).astype(np.int64), np.floor(row).astype(np.int64)
    ok = (c0 >= 0) & (r0 >= 0) & (c0 < z.shape[1] - 1) & (r0 < z.shape[0] - 1)
    c0c, r0c = np.where(ok, c0, 0), np.where(ok, r0, 0)
    fc, fr = col - c0c, row - r0c
    v = (
        z[r0c, c0c] * (1 - fc) * (1 - fr)
        + z[r0c, c0c + 1] * fc * (1 - fr)
        + z[r0c + 1, c0c] * (1 - fc) * fr
        + z[r0c + 1, c0c + 1] * fc * fr
    )
    return np.where(ok, v, np.nan)


def horizon(z, x0, ytop, step, rx, ry, rz, r_i, d_from, d_to, delta, c, gamma=0.0, chunk=400):
    """tan h, shape (n, N), searched from max(r_i, d_from) to d_to at `delta`; -inf where no sample."""
    out = np.full((len(rx), N), -np.inf)
    for k in range(N):
        a = math.radians(k * 360.0 / N) - gamma  # true azimuth -> grid azimuth (A_grid = A_true - gamma)
        ux, uy = math.sin(a), math.cos(a)
        for s in range(0, len(rx), chunk):
            sl = slice(s, s + chunk)
            start = np.maximum(r_i[sl], d_from)[:, None]
            d = start + delta * np.arange(int((d_to - d_from) / delta) + 1)[None, :]
            d = np.where(d <= d_to, d, np.nan)
            zz = bilinear(z, x0, ytop, step, rx[sl, None] + d * ux, ry[sl, None] + d * uy)
            t = (zz - rz[sl, None] - c * d * d) / d
            out[sl, k] = np.max(np.where(np.isnan(t), -np.inf, t), axis=1)
    return out


def h_at(h_deg, az_deg):
    """Horizon (deg) at true azimuth az, linear between the two neighbouring directions; h_deg (n, N)."""
    f = (az_deg % 360.0) / (360.0 / N)
    k0 = math.floor(f) % N
    w = f - math.floor(f)
    return h_deg[:, k0] * (1 - w) + h_deg[:, (k0 + 1) % N] * w


def deg(t):
    return np.degrees(np.arctan(t))


print(f"probe_shadow.py, run {dt.datetime.now(dt.UTC):%Y-%m-%dT%H:%MZ}")

# ---- 1. V-valley, the radiation design's TR3a closed form ----
step = 10.0
xs = np.arange(-3000, 3000 + step, step)
X, Y = np.meshgrid(xs, -xs)
tb = math.tan(math.radians(30))
zv = tb * np.minimum(np.abs(X), 1000.0)
one = np.array([0.0])
for c in (0.0, C_EARTH):
    t = horizon(zv, -3000.0, 3000.0, step, one, one, one * 0, one + 5.0, 0.0, 4300.0, 5.0, c)
    closed = np.abs(np.sin(np.radians(np.arange(N) * 360.0 / N))) * tb - c * 5.0
    print(
        f"V1 V-valley beta 30, W 1000 m, 10 m grid, r_i 5 m, c {c:.3e}: max |tan h - closed form| {np.max(np.abs(t[0] - closed)):.2e};"
        f" h(0) {deg(t[0, 0]):.6f}, h(45) {deg(t[0, 8]):.4f}, h(90) {deg(t[0, 16]):.4f}, h(180) {deg(t[0, 32]):.6f} deg"
    )
hv = deg(horizon(zv, -3000.0, 3000.0, step, one, one, one * 0, one + 5.0, 0.0, 4300.0, 5.0, 0.0))
print(
    f"V2 lit on the axis floor (c = 0): h at A 90 = {h_at(hv, 90)[0]:.4f}; h at A 120, interpolated between k 21 and 22 = {h_at(hv, 120)[0]:.4f}"
    f" (closed form atan(sin 120 tan 30) = {math.degrees(math.atan(math.sin(math.radians(120)) * tb)):.4f}); k 21 {hv[0, 21]:.4f}, k 22 {hv[0, 22]:.4f}"
)
for e in (29.99, 30.01):
    print(f"V2 sun A 90, e {e}: lit {bool(e > h_at(hv, 90)[0])}")
for e in (20.0, 40.0):
    s = math.sin(math.radians(e))
    lit = e > h_at(hv, 90)[0]
    print(f"V3 snow anomaly, floor, sun A 90 e {e}, beam 800: b {lit * s:.4f}, reference {s:.4f}, dI {800 * (lit * s - s):.1f} W/m2")

# ---- 2. the synthetic massif (section 4.5's terrain, rebuilt) ----
n_side, x0, ytop = 1201, 494000.0, 6806000.0
gx_ = (np.arange(n_side) * step) - 6000.0
X, Y = np.meshgrid(gx_, 6000.0 - np.arange(n_side) * step)
r, phi = np.hypot(X, Y), np.arctan2(Y, X)
zs = 300 + 1600 * np.exp(-((r / 2500) ** 2)) * (1 + 0.15 * (1 - np.exp(-((r / 1500) ** 2))) * np.cos(5 * phi))
zs += 250 * np.exp(-((X - 3500) ** 2 + (Y + 3000) ** 2) / 900**2)
rng = np.random.default_rng(46)
k_ = np.hypot(*np.meshgrid(np.fft.fftfreq(n_side, step), np.fft.fftfreq(n_side, step)))
amp = np.zeros_like(k_)
keep = (k_ > 0) & (k_ <= 1 / 40.0)
amp[keep] = k_[keep] ** -1.5
noise = np.real(np.fft.ifft2(amp * np.exp(2j * np.pi * rng.random(k_.shape))))
zs += noise * (25.0 / noise.std())
print(f"synthetic terrain z {zs.min():.1f} .. {zs.max():.1f} m (section 4.5: 234.0 .. 1924.4)")


def receivers(z, x0, ytop, rows, cols):
    gy, gx = np.gradient(z, step)
    n_e, n_n = -gx[rows, cols], gy[rows, cols]  # rows run south, so dz/dnorth = -gy
    n = np.column_stack([n_e, n_n, np.ones_like(n_e)])
    n /= np.linalg.norm(n, axis=1)[:, None]
    return x0 + cols * step, ytop - rows * step, z[rows, cols].astype(float), n


def days(lat, lon, n, h_deg, label, dates):
    slope = np.degrees(np.arccos(n[:, 2]))
    aspect = np.degrees(np.arctan2(n[:, 0], n[:, 1])) % 360
    northish = (slope >= 25) & ((aspect <= 45) | (aspect >= 315))
    for day in dates:
        b_un, b_sh, ref = np.zeros(len(n)), np.zeros(len(n)), 0.0
        for q in range(96):
            e, a, dist = noaa(day + dt.timedelta(minutes=15 * q + 7.5), lat, lon)
            if e <= 0:
                continue
            el, az = math.radians(e), math.radians(a)
            u = np.array([math.cos(el) * math.sin(az), math.cos(el) * math.cos(az), math.sin(el)])
            beam = 1368 / dist**2 * 0.75 ** (1 / math.sin(el))
            cosi = np.maximum(0, n @ u)
            b_un += beam * cosi
            b_sh += beam * cosi * (e > h_at(h_deg, a))
            ref += beam * u[2]
        lost = 1 - b_sh.sum() / b_un.sum()
        d_di = (b_sh - b_un) / 96
        print(
            f"{label} {day:%Y-%m-%d}: beam on the slopes removed by cast shadow {100 * lost:.1f} %;"
            f" daily-mean dI change: mean {d_di.mean():.1f}, 5th pct {np.percentile(d_di, 5):.1f} W/m2;"
            f" receivers losing > 20 W/m2 {100 * np.mean(d_di < -20):.1f} %;"
            f" north-facing >= 25 deg ({northish.sum()}): mean {d_di[northish].mean():.1f} W/m2;"
            f" horizontal reference {ref / 96:.1f} W/m2"
        )


rows = rng.integers(50, n_side - 50, 3000)
cols = rng.integers(50, n_side - 50, 3000)
rx, ry, rz, nrm = receivers(zs, x0, ytop, rows, cols)
t0 = time.perf_counter()
hs = deg(horizon(zs, x0, ytop, step, rx, ry, rz, np.full(len(rx), 10.0), 0.0, 17000.0, 5.0, C_EARTH))
print(f"S0 synthetic: 3000 random nodes, N 64, to 17 km at 5 m, exhaustive NumPy: {time.perf_counter() - t0:.1f} s; median of the 64-direction mean horizon {np.median(hs.mean(axis=1)):.2f} deg, max horizon {hs.max():.2f} deg")
days(61.334, 15.0, nrm, hs, "S1 synthetic", [dt.datetime(2027, 4, 15), dt.datetime(2027, 5, 15), dt.datetime(2027, 6, 21)])

# ---- 3. Bygdin on DTM10 through the existing mesh path's open_dem ----
if len(sys.argv) > 1:
    from shapely import contains_xy

    from tin_engine.dem_input import DemRequest, open_dem
    from tin_engine.io.domain_file import read_domain
    from tin_engine.io.models import Bounds

    dtm = pathlib.Path(sys.argv[1])
    dom = read_domain(pathlib.Path("docs/benchmarks/2026-09-29/bygdin/bygdin_reduced_t20.geojson"))
    bx = dom.polygon.bounds

    def grown(g):
        return Bounds(x_min=bx[0] - g, y_min=bx[1] - g, x_max=bx[2] + g, y_max=bx[3] + g)

    t0 = time.perf_counter()
    near = open_dem(DemRequest(sources=(dtm,), bounds=grown(20000.0)))
    m = near.tile.meta
    zn = np.asarray(near.tile.array, dtype=np.float32)
    print(
        f"B0 Bygdin near raster: open_dem(DemRequest(sources=(DTM10,), bounds=catchment box + 20 km)): {m.cols} x {m.rows} nodes,"
        f" {zn.nbytes / 1e6:.0f} MB float32, {time.perf_counter() - t0:.1f} s; NaN {100 * np.isnan(zn).mean():.2f} %; z {np.nanmin(zn):.0f} .. {np.nanmax(zn):.0f} m;"
        f" step {m.delta_x}, pixel_is_area {m.pixel_is_area}"
    )
    lon, lat = Transformer.from_crs(CRS.from_epsg(25833), CRS.from_epsg(4326), always_xy=True).transform((bx[0] + bx[2]) / 2, (bx[1] + bx[3]) / 2)
    gam = Proj(CRS.from_epsg(25833)).get_factors(lon, lat).meridian_convergence
    gw = Proj(CRS.from_epsg(25833)).get_factors(*Transformer.from_crs(25833, 4326, always_xy=True).transform(bx[0], (bx[1] + bx[3]) / 2)).meridian_convergence
    ge = Proj(CRS.from_epsg(25833)).get_factors(*Transformer.from_crs(25833, 4326, always_xy=True).transform(bx[2], (bx[1] + bx[3]) / 2)).meridian_convergence
    print(f"B0 catchment centre {lat:.4f} N {lon:.4f} E; convergence gamma {gam:.3f} deg at the centre, {gw:.3f} (west edge) .. {ge:.3f} (east edge)")
    zf = np.where(np.isnan(zn), np.nan, zn)
    rows_all, cols_all = np.meshgrid(np.arange(1, m.rows - 1), np.arange(1, m.cols - 1), indexing="ij")
    xa, ya = m.x_min + cols_all * m.delta_x, m.y_max - rows_all * m.delta_y
    inside = contains_xy(dom.polygon, xa, ya) & ~np.isnan(zf[rows_all, cols_all])
    pick = rng.choice(np.flatnonzero(inside), 2000, replace=False)
    rows, cols = rows_all.ravel()[pick], cols_all.ravel()[pick]
    rx, ry, rz, nrm = receivers(zf, m.x_min, m.y_max, rows, cols)
    g = math.radians(gam)
    # grid -> true ENU: azimuths turn clockwise by gamma (A_true = A_grid + gamma, the radiation design's P4 sign)
    nrm = np.column_stack([math.cos(g) * nrm[:, 0] + math.sin(g) * nrm[:, 1], -math.sin(g) * nrm[:, 0] + math.cos(g) * nrm[:, 1], nrm[:, 2]])
    t0 = time.perf_counter()
    tn = horizon(zf, m.x_min, m.y_max, m.delta_x, rx, ry, rz, np.full(len(rx), 10.0), 0.0, 20000.0, 5.0, C_EARTH, gamma=g)
    hn = deg(tn)
    print(f"B1 Bygdin: 2000 random nodes in the catchment, N 64, to 20 km at 5 m: {time.perf_counter() - t0:.1f} s; z {rz.min():.0f} .. {rz.max():.0f} m; median mean horizon {np.median(hn.mean(axis=1)):.2f} deg, max {hn.max():.2f} deg")
    days(lat, lon, nrm, hn, "B2 Bygdin", [dt.datetime(2027, 4, 15), dt.datetime(2027, 5, 15), dt.datetime(2027, 6, 21)])
    del near, zn, zf
    # far field: 20-60 km, from a 10 m canvas block-averaged to 100 m (TR3b's far raster, built the slow way)
    t0 = time.perf_counter()
    try:
        far = open_dem(DemRequest(sources=(dtm,), bounds=grown(60000.0)))
        fm = far.tile.meta
        za = np.asarray(far.tile.array, dtype=np.float32)
        rr, cc = (fm.rows // 10) * 10, (fm.cols // 10) * 10
        with np.errstate(invalid="ignore"):
            zb = np.nanmean(za[:rr, :cc].reshape(rr // 10, 10, cc // 10, 10), axis=(1, 3))
        print(f"B3 far canvas: {fm.cols} x {fm.rows} nodes at 10 m ({za.nbytes / 1e6:.0f} MB), block means to 100 m, {time.perf_counter() - t0:.1f} s; NaN at 100 m {100 * np.isnan(zb).mean():.1f} %")
        x0b, ytb = fm.x_min + 4.5 * fm.delta_x, fm.y_max - 4.5 * fm.delta_y
        del far, za
        tf = horizon(zb, x0b, ytb, 100.0, rx, ry, rz, np.zeros(len(rx)), 20000.0, 60000.0, 50.0, C_EARTH, gamma=g)
        hf = deg(tf)
        above = hf > hn
        print(f"B3 far field (20-60 km) above the near horizon in {100 * above.mean():.2f} % of (receiver, direction) pairs; where above, far horizon median {np.median(hf[above]) if above.any() else float('nan'):.2f}, max {hf[above].max() if above.any() else float('nan'):.2f} deg; near horizon there median {np.median(hn[above]) if above.any() else float('nan'):.2f} deg")
        hb = np.maximum(hn, hf)
        for day in (dt.datetime(2027, 4, 15), dt.datetime(2027, 5, 15), dt.datetime(2027, 6, 21)):
            flips, lost, tot, emax = 0, 0.0, 0.0, -90.0
            for q in range(96):
                e, a, dist = noaa(day + dt.timedelta(minutes=15 * q + 7.5), lat, lon)
                if e <= 0:
                    continue
                el, az = math.radians(e), math.radians(a)
                u = np.array([math.cos(el) * math.sin(az), math.cos(el) * math.cos(az), math.sin(el)])
                beam = 1368 / dist**2 * 0.75 ** (1 / math.sin(el))
                cosi = np.maximum(0, nrm @ u)
                ln, lb = e > h_at(hn, a), e > h_at(hb, a)
                f = (ln != lb) & (cosi > 0)
                flips += int(f.sum())
                if f.any():
                    emax = max(emax, e)
                lost += float((beam * cosi * (ln & ~lb)).sum())
                tot += float((beam * cosi * ln).sum())
            print(f"B3 {day:%Y-%m-%d}: far field flips {flips} (receiver, 15-min) lit states; highest sun elevation at a flip {emax:.2f} deg; beam energy it removes {100 * lost / tot:.4f} % of the near-only lit beam")
    except Exception as exc:  # a probe reports and goes on
        print(f"B3 far canvas refused: {type(exc).__name__}: {exc}"[:400])
print(f"peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30:.2f} GiB")
