"""Design probe for increment 33: triangles by distance band on the corridor.

Reads corr_t{T}.ply (binary, `rasputin mesh ... --out corr_tT.ply`) for T in
20, 10, 5, 2, 1 and line_3.pkl from probe_line.py; samples 400 000 triangles
per mesh (seed 1) and estimates the ramp's triangle count.
"""

import pickle

import numpy as np
import shapely

EDGES = [0, 100, 250, 500, 1000, 1500, 2000, 2500, 3000, 1e9]
MIDS = np.array([50, 175, 375, 750, 1250, 1750, 2250, 2750, 4000.0])
TOLS = [1, 2, 5, 10, 20]


def load(tol: int) -> tuple[np.ndarray, np.ndarray]:
    with open(f"corr_t{tol}.ply", "rb") as f:
        b = f.read()
    h = b.index(b"end_header\n") + len(b"end_header\n")
    hdr = b[:h].decode()
    nv = int(hdr.split("element vertex ")[1].split()[0])
    nf = int(hdr.split("element face ")[1].split()[0])
    v = np.frombuffer(b, dtype="<f8", count=3 * nv, offset=h).reshape(nv, 3)
    face = np.dtype([("n", "u1"), ("i", "<u4", 3)])
    return v, np.frombuffer(b, dtype=face, count=nf, offset=h + 24 * nv)["i"]


def areas(p: np.ndarray) -> np.ndarray:
    ux, uy = p[:, 1, 0] - p[:, 0, 0], p[:, 1, 1] - p[:, 0, 1]
    wx, wy = p[:, 2, 0] - p[:, 0, 0], p[:, 2, 1] - p[:, 0, 1]
    return 0.5 * np.abs(ux * wy - wx * uy)


def main() -> None:
    with open("line_3.pkl", "rb") as f:
        ml, _corr = pickle.load(f)
    segs = []
    for g in ml.simplify(1.0, preserve_topology=False).geoms:
        c = np.asarray(g.coords)
        segs.append(shapely.linestrings(np.stack([c[:-1], c[1:]], 1)))
    tree = shapely.STRtree(np.concatenate(segs))
    print("segments after 1 m simplification", sum(len(s) for s in segs))
    rng = np.random.default_rng(1)
    res = {}
    for tol in (20, 10, 5, 2, 1):
        v, f = load(tol)
        n = len(f)
        idx = rng.choice(n, min(n, 400000), replace=False)
        scale = n / len(idx)
        p = v[f[idx].astype(np.int64)][:, :, :2]
        centres = shapely.points(p.mean(1))
        _, d = tree.query_nearest(centres, return_distance=True, all_matches=False)
        k = np.digitize(d, EDGES) - 1
        count = np.bincount(k, minlength=9) * scale
        res[tol] = (count, np.bincount(k, weights=areas(p), minlength=9) * scale)
        print(tol, "triangles", n, "per band", np.round(count).astype(int).tolist(), flush=True)
    ar1 = res[1][1]
    print("band areas km2", np.round(ar1 / 1e6, 1).tolist())
    dens = {t: res[t][0] / res[t][1] for t in TOLS}

    def estimate(r0: float, r1: float, near: float = 1.0, far: float = 20.0) -> int:
        tot = 0.0
        for b in range(9):
            d = float(MIDS[b])
            t = near if d <= r0 else far if d >= r1 else near + (far - near) * (d - r0) / (r1 - r0)
            log_d = np.interp(np.log(t), np.log(TOLS), [np.log(dens[s][b]) for s in TOLS])
            tot += np.exp(log_d) * ar1[b]
        return int(tot)

    print("estimate linear 0-3000:", estimate(0, 3000), " flat to 100 then linear:",
          estimate(100, 3000), " flat to 500:", estimate(500, 3000))  # fmt: skip


if __name__ == "__main__":
    main()
