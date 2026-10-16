"""Increment 40's probe: how the flood's time grows with the window.

Run as `native_flood.py` is. The ANADEM fixture is mirrored into a seamless
640 x 640 block and tiled to 2560 and 6400 nodes a side; `upstream` (from the
centre node) and `accumulate` are timed once each on a unit lattice. The core
is today's: increment 40 changes no C++.
"""

# ruff: noqa: E402 -- the imports follow the sys.path set-up on purpose

import sys
import time
from pathlib import Path

sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, sys.argv[1])

import numpy as np

from tin_engine import _core
from tin_engine.dem_input import repository_for

repo, _ = repository_for((Path("tests/fixtures/velhas/anadem_velhas.tif"),))
a = np.asarray(repo.load(repo.footprints()[0].name).array)
block = np.block([[a, a[:, ::-1]], [a[::-1], a[::-1, ::-1]]])
for k in (1, 4, 10):
    z = np.ascontiguousarray(np.tile(block, (k, k)) if k > 1 else a)
    z.flags.writeable = False
    v = _core.raster_view(z, x_min=0.0, y_max=0.0, delta_x=1.0, delta_y=1.0, nodata=-9999.0)
    s = np.zeros(z.shape, np.uint8)
    s[z.shape[0] // 2, z.shape[1] // 2] = 1
    t = time.perf_counter()
    _core.upstream(v, s)
    t_up = time.perf_counter() - t
    t = time.perf_counter()
    _core.accumulate(v)
    t_acc = time.perf_counter() - t
    print(f"{z.shape[0]} x {z.shape[1]} = {z.size / 1e6:.2f} M nodes: upstream {t_up:.3f} s "
          f"({z.size / t_up / 1e6:.1f} M nodes/s), accumulate {t_acc:.3f} s")
