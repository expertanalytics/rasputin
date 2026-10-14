# ruff: noqa: E501, N803, N806, F821 -- a one-off design probe: long report lines, the kernel's S, D, W, T, P, M names, and names defined by the executed probe_spring.py
"""Increment 46 PR 8 design probe: linear bins (a continuous SCA), in NumPy.

Run from the worktree root with a venv whose `rasputin` is built:
    .venv/bin/python docs/increments/46-probes/probe_sca.py WORKDIR > docs/increments/46-probes/probe_sca.out
It runs probe_spring.py's terrain, mesh, geometry and sun (the same code, executed
up to its D1 section), then section 12's kernel as a NumPy reference: the hand
fixtures L1-L10, the distribution figures D3, a flat ten-bin run that must
reproduce probe_spring.out's first row (the reference is today's kernel on flat
bins), and the ten-bin rows at cv 0.4. Everything printed is quoted by section 12 of
docs/increments/46-hbv-snow.md.
"""

import math
import pathlib
import statistics
import sys

import numpy as np

here = pathlib.Path(__file__).parent
sys.path.insert(0, str(here))
spring_src = (here / "probe_spring.py").read_text("utf-8")
exec(compile(spring_src.split("# ---- D1 ----")[0], "probe_spring.py", "exec"))


# ---- the distribution (section 12.4) ----
def lognormal_linear(cv, n):
    """Bin means and deep-edge factors: each band of a mean-1 log-normal matched in mean and variance by a uniform."""
    if cv == 0:
        return np.ones(1), np.ones(1)
    s = math.sqrt(math.log1p(cv * cv))
    mu = -s * s / 2
    nd = statistics.NormalDist()

    def cdf(z, shift):
        return 0.0 if z == -math.inf else 1.0 if z == math.inf else nd.cdf(z - shift)

    f, h = [], []
    for k in range(n):
        za = nd.inv_cdf(k / n) if k > 0 else -math.inf
        zb = nd.inv_cdf((k + 1) / n) if k + 1 < n else math.inf
        m1 = math.exp(mu + s * s / 2) * (cdf(zb, s) - cdf(za, s)) * n
        m2 = math.exp(2 * mu + 2 * s * s) * (cdf(zb, 2 * s) - cdf(za, 2 * s)) * n
        sd = math.sqrt(max(0.0, m2 - m1 * m1))
        f.append(m1)
        h.append(m1 + math.sqrt(3) * sd)
    f, h = np.array(f), np.array(h)
    mean = math.fsum(f) / n
    return f / mean, h / mean


f10, h10 = lognormal_linear(0.4, 10)
print("D3 lognormal(0.4, 10) factors:", np.round(f10, 4).tolist())
print("D3 lognormal(0.4, 10) deep:   ", np.round(h10, 4).tolist())
print("D3 lognormal(0.4, 10) shallow:", np.round(2 * f10 - h10, 4).tolist())
f1, h1 = lognormal_linear(0.4, 1)
print(f"D3 lognormal(0.4, 1): factors {f1.tolist()}, deep {h1[0]:.9f}, shallow {2 - h1[0]:.9f}")
for n_ in (1, 2, 4, 10, 64):
    lo_, hi_ = 0.0, 3.0
    for _ in range(60):
        mid = (lo_ + hi_) / 2
        ff, hh = lognormal_linear(mid, n_)
        lo_, hi_ = (mid, hi_) if (2 * ff - hh).min() >= 0 else (lo_, mid)
    print(f"D3 largest cv with no negative shallow edge, {n_} bins: {lo_:.4f}")


# ---- the kernel (section 12.3), vectorised over triangles x bins ----
def cover(S, D):
    full = 2 * S >= D
    return np.where(S == 0, 0.0, np.where(full, 1.0, np.divide(2 * S, D, out=np.zeros_like(S), where=D > 0)))


def step(st, f, h, T, P, M, th, ts, cx, cfr, lw, days, tx=0.0):
    """One step. M is the melt factor m (n,), th the melt threshold (n,); T, P (n,)."""
    S, D, W = st["S"], st["D"], st["W"]
    Tm, Pm = T[:, None], P[:, None]
    snow = np.where(Tm < tx, Pm, 0.0)
    rain = Pm - snow
    respread = (snow > 0) & ((S == 0) | (2 * S < D))
    D[:] = np.where(respread, 2 * S, D)
    S += f * snow
    D += h * snow
    c = cover(S, D)
    W += c * rain
    out = ((1 - c) * rain).sum(axis=1)
    melting = Tm > th[:, None]
    pot = np.where(melting, M[:, None] * days * (Tm - th[:, None]), 0.0)
    has = melting & (S > 0)
    bare = has & (pot >= D)
    full = 2 * S >= D
    stay = has & ~bare & full & (pot <= 2 * S - D)
    part = has & ~bare & ~stay
    Dn = D - pot
    with np.errstate(divide="ignore", invalid="ignore"):
        Sp = np.where(full, Dn * Dn / (4 * (D - S)), S * (Dn / D) ** 2)
    Snew = np.where(bare, 0.0, np.where(stay, S - pot, np.where(part, Sp, S)))
    Dnew = np.where(bare, 0.0, np.where(stay | part, Dn, D))
    Dnew = np.where(Snew == 0, 0.0, Dnew)
    delta = np.where(bare, S, np.where(stay, pot, np.where(part, S - Snew, 0.0)))
    W += delta
    freezing = ~melting & (Tm < ts)
    F = np.where(freezing, np.minimum(c * cfr * cx * days * (ts - Tm), W), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        Dref = np.where(full | (S == 0), Dnew + F, Dnew * (Snew + F) / Snew)
    Dnew = np.where(F > 0, Dref, Dnew)
    Snew = Snew + F
    W -= F
    o = np.maximum(0.0, W - lw * Snew)
    W -= o
    S[:], D[:] = Snew, Dnew
    return out + o.sum(axis=1)


def fixture(name, st, f, h, T, P, M=3.0, th=0.0, ts=0.0, cfr=0.05, lw=0.1, tx=0.0):
    out = step(st, np.array(f), np.array(h), np.array([T]), np.array([P]), np.array([M]), np.array([th]), ts, M if M else 3.0, cfr, lw, 1.0, tx)[0] / len(f)
    S, D, W = st["S"][0], st["D"][0], st["W"][0]
    print(f"{name}: sp {np.round(S, 6).tolist()} deep {np.round(D, 6).tolist()} sw {np.round(W, 6).tolist()} outflow {out:.6f} swe {np.mean(S + W):.6f} sca {cover(st['S'], st['D'])[0].mean():.6f}")


def new(swe, f=(1.0,), h=(1.5,)):
    return {"S": np.array([[swe * x for x in f]]), "D": np.array([[swe * x for x in h]]), "W": np.zeros((1, len(f)))}


one, deep = (1.0,), (1.5,)
s = new(100.0)
print(f"L0 initial: sp 100 deep 150 sca {cover(s['S'], s['D'])[0, 0]:.6f}")
fixture("L1 melt 15, stays full", s, one, deep, T=5.0, P=0.0)
s = new(100.0)
fixture("L2 melt 60, full to partial", s, one, deep, T=20.0, P=0.0)
fixture("L3 melt 30 from partial", s, one, deep, T=10.0, P=0.0)
fixture("L4 melt 60, out", s, one, deep, T=20.0, P=0.0)
s = new(100.0)
fixture("L5a lw 2: melt 60", s, one, deep, T=20.0, P=0.0, lw=2.0)
s5 = {k: v.copy() for k, v in s.items()}
fixture("L5 rain 10 on the partial bin, T = tx = ts = 0", s5, one, deep, T=0.0, P=10.0, lw=2.0)
fixture("L6 refreeze on the partial bin, T = -10", s, one, deep, T=-10.0, P=0.0, lw=2.0)
s = new(100.0)
fixture("L7a melt 60", s, one, deep, T=20.0, P=0.0, cfr=0.0)
fixture("L7b melt 30", s, one, deep, T=10.0, P=0.0, cfr=0.0)
fixture("L7 snow 10 on the partial bin, cfr 0", s, one, deep, T=-1.0, P=10.0, cfr=0.0)
fixture("L7c melt 10 after it", s, one, deep, T=10.0 / 3.0, P=0.0, cfr=0.0)
s = new(0.0)
fixture("L8 snow 10 on a bare bin", s, one, deep, T=-1.0, P=10.0, cfr=0.0)
fixture("L8b melt 7.5", s, one, deep, T=2.5, P=0.0, cfr=0.0)
s = new(100.0)
fixture("L9a melt 60", s, one, deep, T=20.0, P=0.0)
fixture("L9 T = ts on the partial bin", s, one, deep, T=0.0, P=0.0)
# L10: melt-only depletion curve, one bin, deep 1.5, cfr 0
s, worst = new(100.0), 0.0
cum = 0.0
for k_ in range(400):
    Tk = 0.5 + 1.0 * (k_ % 7)
    step(s, np.array(one), np.array(deep), np.array([Tk]), np.zeros(1), np.array([3.0]), np.zeros(1), 0.0, 3.0, 0.0, 0.1, 1 / 24)
    cum += 3.0 / 24 * Tk
    want = min(1.0, max(0.0, (150.0 - cum) / 100.0))
    worst = max(worst, abs(cover(s["S"], s["D"])[0, 0] - want))
print(f"L10 melt-only depletion: 400 hourly steps, cumulative melt {cum:.4f} mm; max |sca - clamp((150 - melt)/100, 0, 1)| = {worst:.3e}")

# ---- the spring, as probe_spring.py, with the kernel above ----
TX, TS, CX, LW, CFR, RR, DD = 0.0, 0.0, 3.0, 0.1, 0.05, 0.015, 0.02
start, steps, dts = dt.datetime(2027, 4, 1), 91 * 24, 3600.0
sigma_ = math.sqrt(math.log(1 + 0.16))
fl10 = np.array([math.exp(sigma_ * statistics.NormalDist().inv_cdf((q + 0.5) / 10)) for q in range(10)])
fl10 /= fl10.mean()
runs = {
    "flat_none_N10 (must match probe_spring.out row none)": ("none", fl10, fl10),
    "linear_none_N10": ("none", f10, h10),
    "linear_melt_factor_N10": ("melt_factor", f10, h10),
    "linear_threshold_N10": ("threshold", f10, h10),
    "linear_none_N1": ("none", f1, h1),
    "linear_melt_factor_N1": ("melt_factor", f1, h1),
}
state = {}
for name, (mode, f, h) in runs.items():
    state[name] = dict(mode=mode, f=f, h=h, S=np.tile(600.0 * f, (len(tri), 1)), D=np.tile(600.0 * h, (len(tri), 1)),
                       W=np.zeros((len(tri), len(f))), out=np.zeros(len(tri)), sca=[np.ones(len(tri))], swe=[np.full(len(tri), 600.0)], mono=0.0)
dtd = dts / 86400
for st_ in range(steps):
    tm = start + dt.timedelta(seconds=dts * (st_ + 0.5))
    hours = st_ + 0.5
    hsol = (tm.hour + tm.minute / 60 + 1.0) % 24
    T = (-3 + 17 * hours / (91 * 24) + 4 * math.cos(2 * math.pi * (hsol - 14) / 24)) - 0.006 * (zc - 300)
    anomaly = np.zeros(len(tri))
    for q in range(4):
        u, B = sun(start + dt.timedelta(seconds=dts * st_ + dts / 4 * (q + 0.5)))
        if B > 0:
            anomaly += B * (np.maximum(0, nrm @ u) - max(0, u[2]))
    anomaly /= 4
    for s_ in state.values():
        m = np.full(len(tri), CX)
        th = np.full(len(tri), TS)
        if s_["mode"] == "melt_factor":
            m = np.maximum(0, CX + RR * anomaly)
        if s_["mode"] == "threshold":
            th = TS - DD * anomaly
        before = cover(s_["S"], s_["D"]).mean(axis=1)
        s_["out"] += step(s_, s_["f"], s_["h"], T, np.zeros(len(tri)), m, th, TS, CX, CFR, LW, dtd) / len(s_["f"])
        after = cover(s_["S"], s_["D"]).mean(axis=1)
        s_["mono"] = max(s_["mono"], float((after - before).max()))
    if (st_ + 1) % 24 == 0:
        for s_ in state.values():
            s_["sca"].append(cover(s_["S"], s_["D"]).mean(axis=1))
            s_["swe"].append((s_["S"] + s_["W"]).mean(axis=1))
w = area3 / area3.sum()
for name, s_ in state.items():
    sca, swe = np.array(s_["sca"]), np.array(s_["swe"])
    below = sca < 0.5
    mo = np.where(below.any(axis=0), below.argmax(axis=0), -1)
    first = np.where((sca < 1).any(axis=0), (sca < 1).argmax(axis=0), -1)
    bare = np.where((sca == 0).any(axis=0), (sca == 0).argmax(axis=0), -1)
    both = (first >= 0) & (bare >= 0)
    frac = ((sca > 0) & (sca < 1)).mean()
    print(f"{name}: domain SCA d30 {w @ sca[30]:.4f} d60 {w @ sca[60]:.4f} d90 {w @ sca[90]:.4f}; SWE d30 {w @ swe[30]:.2f} d60 {w @ swe[60]:.2f} d90 {w @ swe[90]:.2f} mm; "
          f"daily frames with 0 < sca < 1: {frac:.4f}; largest SCA rise in a step {s_['mono']:.3e}")
    for cls, sel in (("south", south), ("north", north)):
        print(f"   {cls}: median melt-out day {np.median(mo[sel]):.1f}, median days first-bare..bare {np.median(bare[sel & both] - first[sel & both]) if (sel & both).any() else float('nan'):.1f}")
    print(f"   I1 over the run: max |600 - swe_end - outflow| = {abs(600.0 - swe[-1] - s_['out']).max():.3e} mm")
