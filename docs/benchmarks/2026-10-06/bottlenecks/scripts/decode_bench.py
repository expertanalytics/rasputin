"""decode, decomposed: open_dem's steps timed with perf_counter (wrappers add
one call each), node counts, and two candidates checked for an identical canvas
and identical seams."""
import sys, time, json, pathlib, functools, hashlib
import numpy as np, shapely
import tin_engine.dem_input as di, tin_engine.mosaic as mo, tin_engine.io.repository as rp, tin_engine.io.cog as cog
from tin_engine.domain import read_domain
c = sys.argv[1]
dom = read_domain(pathlib.Path(f"/Users/skavhaug/projects/rasputin_scratch/norway/{c}/{c}_outline_nve.geojson"))
DEM = pathlib.Path("/Users/skavhaug/projects/rasputin_data/DTM10_UTM33_20260925")
T, N = {}, {}
def timed(mod, name, key, count=None):
    f = getattr(mod, name)
    @functools.wraps(f)
    def w(*a, **k):
        t = time.perf_counter(); r = f(*a, **k); T[key] = T.get(key, 0) + time.perf_counter() - t
        if count: N[count] = N.get(count, 0) + count_of(a, r)
        return r
    setattr(mod, name, w); return f
count_of = lambda a, r: int(np.size(r)) if r is not ... else 0
orig_cov = timed(mo, "_covered", "assemble: _covered (seam report mask)", "overlap nodes tested")
timed(mo, "_decide", "assemble: _decide")
timed(mo, "_loaded", "assemble: load windows (decode)")
timed(mo, "_seam", "assemble: _seam")
orig_dw = rp.decode_window
timed(di, "_domain_plan", "domain plan")
_fp = rp.TiffDemRepository.footprints
def fp(self):
    t = time.perf_counter(); r = _fp(self); T["footprints (all 254 headers)"] = T.get("footprints (all 254 headers)", 0) + time.perf_counter() - t; return r
rp.TiffDemRepository.footprints = fp
def run():
    T.clear(); N.clear()
    req = di.DemRequest(sources=[DEM], domain=dom)
    t = time.perf_counter(); out = di.open_dem(req); T["open_dem total"] = time.perf_counter() - t
    return out
# footprints alone, on a fresh repository
t = time.perf_counter(); repo = rp.TiffDemRepository.from_directory(DEM); fps = repo.footprints(); tf = time.perf_counter() - t
res = {"catchment": c, "footprints (all headers)_s": round(tf, 3), "tiles in directory": len(fps)}
base = [run() for _ in range(3)]; res["baseline"] = {k: round(v, 3) for k, v in T.items()}; res["baseline"].update(N)
ref = base[-1]; arr = np.asarray(ref.tile.array)
res["canvas"] = dict(rows=arr.shape[0], cols=arr.shape[1], dtype=str(arr.dtype), MB=round(arr.nbytes/1e6),
                     tiles=len(ref.plan.tiles), sha256=hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()[:16])
# Candidate 1: decode_window with 10 threads instead of 4
rp.decode_window = functools.partial(orig_dw, threads=10)
r1 = [run() for _ in range(3)]; res["cand: decode threads 10"] = {k: round(v, 3) for k, v in T.items()}
res["cand: decode threads 10"]["identical canvas"] = bool(np.array_equal(np.asarray(r1[-1].tile.array), arr, equal_nan=True))
rp.decode_window = orig_dw
# Candidate 2: the seam report's domain mask tested only at nodes whose gap can count
import inspect
src = inspect.getsource(mo.assemble)
old = """            kept = _covered(box, plan.meta, needed)
            first, second = strips[a.name, b.name][kept], strips[b.name, a.name][kept]
"""
new = """            sa, sb = strips[a.name, b.name], strips[b.name, a.name]
            nd = plan.meta.nodata
            cand = _valid(sa, nd) & _valid(sb, nd)
            cand[cand] = np.abs(sa[cand].astype(np.float64) - sb[cand].astype(np.float64)) >= SEAM_THRESHOLD
            rr, cc = np.nonzero(cand)
            if needed is not None:
                xs = plan.meta.x_min + (box.col0 + cc) * plan.meta.delta_x
                ys = plan.meta.y_max - (box.row0 + rr) * plan.meta.delta_y
                cov = shapely.intersects_xy(needed, xs, ys)
                rr, cc = rr[cov], cc[cov]
            first, second = sa[rr, cc], sb[rr, cc]
"""
assert old in src
ns = dict(vars(mo)); exec(src.replace(old, new), ns)
real_assemble = mo.assemble
di.assemble = ns["assemble"]
r2 = [run() for _ in range(3)]
di.assemble = real_assemble
res["cand: seam mask on qualifying nodes only"] = {k: round(v, 3) for k, v in T.items()}
res["cand: seam mask on qualifying nodes only"]["identical canvas"] = bool(np.array_equal(np.asarray(r2[-1].tile.array), arr, equal_nan=True))
res["cand: seam mask on qualifying nodes only"]["identical seams"] = r2[-1].seams == ref.seams
res["seams"] = [str(x) for x in ref.seams]
# Interleaved medians of open_dem total, 5 rounds
import statistics
lean = ns["assemble"]
variants = {"baseline": (real_assemble, orig_dw), "decode threads 10": (real_assemble, functools.partial(orig_dw, threads=10)),
            "seam mask on qualifying nodes": (lean, orig_dw), "both": (lean, functools.partial(orig_dw, threads=10))}
tot = {k: [] for k in variants}; parts = {k: [] for k in variants}
for _ in range(5):
    for k, (asm, dw) in variants.items():
        di.assemble, rp.decode_window = asm, dw
        out = run(); tot[k].append(T["open_dem total"]); parts[k].append(dict(T))
        assert np.array_equal(np.asarray(out.tile.array), arr, equal_nan=True) and out.seams == ref.seams
di.assemble, rp.decode_window = real_assemble, orig_dw
res["interleaved medians (5 rounds, canvas and seams identical every run)"] = {
    k: {"open_dem total": round(statistics.median(v), 3), "runs": [round(x, 3) for x in v],
        **{kk: round(statistics.median(p.get(kk, 0) for p in parts[k]), 3) for kk in parts[k][0] if kk != "open_dem total"}}
    for k, v in tot.items()}
print(json.dumps(res, indent=1, default=str))
