"""Prints the DEM tile array refine is handed on the bench's tile run (no
--domain), and times refine there: dtype, flags, base chain, and whether its
pages were touched before refine (@perf).
python tile_probe.py PKG DEM OUT"""
import sys
import time

pkg, dem, out = sys.argv[1:4]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)
import tin_engine.cli as cli  # noqa: E402

real_core, real_refine = cli.to_core, cli.refine
seen = []


def to_core(tile):
    if not seen:
        a = tile.array
        chain, b = [], a.base
        while b is not None and len(chain) < 5:
            chain.append(f"{type(b).__name__}{getattr(b, 'shape', '')}")
            b = getattr(b, "base", None)
        print("array", a.dtype, a.shape, a.strides, "C" if a.flags.c_contiguous else "-",
              "own" if a.flags.owndata else "view", "w" if a.flags.writeable else "ro",
              "aligned" if a.flags.aligned else "unaligned", hex(a.ctypes.data % 4096),
              "base:", chain, file=sys.stderr)  # fmt: skip
    seen.append(1)
    return real_core(tile)


def refine(*a, **k):
    t0 = time.perf_counter()
    r = real_refine(*a, **k)
    print("refine_s", round(time.perf_counter() - t0, 4), file=sys.stderr)
    return r


cli.to_core, cli.refine = to_core, refine
sys.argv = ["rasputin", "mesh", "--dem", dem, "--tolerance", "1", "--out", out, "--binary"]
cli.app()
