"""Dump the polygons and region ring land_cover_map._fixed sees, then stop."""
import sys, numpy as np, shapely
from tin_engine import land_cover_map as lcm
def dump(pairs, xy):
    np.savez(sys.argv[1], wkb=np.array(shapely.to_wkb([g for _, g in pairs]), dtype=object),
             codes=np.array([c for c, _ in pairs]), xy=xy, allow_pickle=True)
    raise SystemExit(0)
lcm._fixed = dump
from tin_engine.cli import app
app(sys.argv[2:], standalone_mode=False)
