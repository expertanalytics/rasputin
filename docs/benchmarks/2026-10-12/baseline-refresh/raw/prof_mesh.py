"""Run `rasputin mesh ...` from a bench.py pkg tree under cProfile.

usage: python prof_mesh.py <pkg> <out.prof> mesh <args...>
"""

import cProfile
import sys

pkg, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, pkg)
import tin_engine  # noqa: E402
import tin_engine._core  # noqa: E402

for module in (tin_engine, tin_engine._core):
    if not str(module.__file__).startswith(pkg):
        raise SystemExit(f"{module.__file__} is not under {pkg}")
from tin_engine.cli import app  # noqa: E402

sys.argv = ["rasputin", *sys.argv[3:]]
profile = cProfile.Profile()
profile.enable()
try:
    app(standalone_mode=False)
except SystemExit:
    pass
profile.disable()
profile.dump_stats(out)
