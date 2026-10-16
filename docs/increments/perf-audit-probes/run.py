"""Run a copy of tin_engine (here: master 91262c73's src_python/tin_engine, with
the main venv's _core copied in beside it) instead of the editable install.

Layout: this file beside pkg/tin_engine/. Usage:
    python run.py [--profile OUT] -- <rasputin args...>
"""

import contextlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# The editable install's finder comes first in sys.meta_path; drop it so the copy loads.
sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in type(f).__module__]
sys.path.insert(0, str(HERE / "pkg"))

import tin_engine  # noqa: E402

print("tin_engine from", tin_engine.__file__, file=sys.stderr)

args = sys.argv[1:]
prof = None
if args and args[0] == "--profile":
    prof, args = args[1], args[2:]
if args and args[0] == "--":
    args = args[1:]

from tin_engine.cli import app  # noqa: E402

sys.argv = ["rasputin", *args]
if prof:
    import cProfile
    import pstats

    p = cProfile.Profile()
    with contextlib.suppress(SystemExit):
        p.runcall(app)
    p.dump_stats(prof)
    st = pstats.Stats(prof)
    st.sort_stats("cumulative").print_stats(60)
    st.sort_stats("tottime").print_stats(40)
else:
    app()
