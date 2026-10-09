"""Design probe for increment 33: triangles by distance band on the Geilo-Ål section.

Reads sec_t{T}.ply (ASCII, `rasputin mesh ... --ascii`) for T in 20, 10, 5, 2, 1
and section.pkl from section.py; estimates the ramp's triangle count.
"""

import pickle

import numpy as np
import shapely

EDGES = [0, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 1e9]
MIDS = np.array([50, 175, 375, 750, 1250, 1750, 2250, 2750, 4000.0])
TOLS = (1, 2, 5, 10, 20)


def load(tol: int) -> tuple[np.ndarray, np.ndarray]:
    with open(f"sec_t{tol}.ply") as f:
        nv = nf = 0
        while True:
            line = f.readline()
            if line.startswith("element vertex"):
                nv = int(line.split()[2])
            if line.startswith("element face"):
                nf = int(line.split()[2])
            if line.startswith("end_header"):
                break
        v = np.loadtxt(f, max_rows=nv)
        faces = np.loadtxt(f, max_rows=nf, dtype=np.int64)[:, 1:4]
    return v, faces


def areas(p: np.ndarray) -> np.ndarray:
    ux, uy = p[:, 1, 0] - p[:, 0, 0], p[:, 1, 1] - p[:, 0, 1]
    wx, wy = p[:, 2, 0] - p[:, 0, 0], p[:, 2, 1] - p[:, 0, 1]
    return 0.5 * np.abs(ux * wy - wx * uy)


def ramp(d: float, r0: float, r1: float, near: float = 1.0, far: float = 20.0) -> float:
    if d <= r0:
        return near
    if d >= r1:
        return far
    return near + (far - near) * (d - r0) / (r1 - r0)


def main() -> None:
    with open("section.pkl", "rb") as f:
        sec, _dom = pickle.load(f)
    res = {}
    for tol in (20, 10, 5, 2, 1):
        v, faces = load(tol)
        p = v[faces][:, :, :2]
        d = shapely.distance(shapely.points(p.mean(axis=1)), sec)
        k = np.digitize(d, EDGES) - 1
        n = len(EDGES) - 1
        res[tol] = (np.bincount(k, minlength=n), np.bincount(k, weights=areas(p), minlength=n))
        print(tol, "triangles", len(faces), "per band", res[tol][0].tolist())
    ar1 = res[1][1]
    print("band areas km2 (t=1 mesh)", np.round(ar1 / 1e6, 2).tolist())
    dens = {t: res[t][0] / res[t][1] for t in TOLS}
    logs = np.log(TOLS)
    cases = {"linear 0-3000": (0, 3000), "flat 0-100 then linear to 3000": (100, 3000)}
    for name, (r0, r1) in cases.items():
        tot = 0.0
        for b in range(len(MIDS)):
            t = ramp(float(MIDS[b]), r0, r1)
            tot += np.exp(np.interp(np.log(t), logs, [np.log(dens[s][b]) for s in TOLS])) * ar1[b]
        print(name, "estimated triangles", int(tot))
    step = sum(dens[1][b] * ar1[b] for b in range(8)) + dens[20][8] * ar1[8]
    print("step 1 m to 3 km, then 20", int(step))
    print("density per km2 at 1 m, band 0-100:", round(dens[1][0] * 1e6), " band >3000:",
          round(dens[1][8] * 1e6), "; at 20 m band >3000:", round(dens[20][8] * 1e6))  # fmt: skip


if __name__ == "__main__":
    main()
