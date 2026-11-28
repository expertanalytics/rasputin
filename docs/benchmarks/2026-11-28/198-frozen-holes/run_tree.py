"""Run `rasputin` from a given src_python folder, past the venv's editable finder.

Usage: python run_tree.py SRC_PYTHON mesh ...  The folder must hold tin_engine
with a built _core beside it (increment 198 changes no C++, so the head's
_core serves master's Python too).
"""

import sys

src = sys.argv.pop(1)
sys.meta_path[:] = [f for f in sys.meta_path if "editable" not in type(f).__module__.lower()]
sys.path.insert(0, src)
for name in [m for m in sys.modules if m == "tin_engine" or m.startswith("tin_engine.")]:
    del sys.modules[name]
from tin_engine.plugins import main  # noqa: E402

assert __import__("tin_engine").__file__.startswith(src), __import__("tin_engine").__file__
sys.argv[0] = "rasputin"
sys.exit(main())
