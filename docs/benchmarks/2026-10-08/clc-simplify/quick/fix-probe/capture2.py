import sys, numpy as np
import tin_engine.feature_input as fi
from tin_engine.border_simplify import BorderResult
fi.simplify_borders = lambda polys, band: BorderResult(tuple(polys), None)
sys.argv = [sys.argv[0], "clip_nosimp.pkl"]
exec(open("capture.py").read())
