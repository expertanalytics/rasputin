"""Run `rasputin <argv>` under cProfile, writing the profile to argv[1]."""
import cProfile, sys
out, sys.argv = sys.argv[1], ["rasputin", *sys.argv[2:]]
from tin_engine.cli import app
prof = cProfile.Profile()
try:
    prof.runcall(app)
except SystemExit as exc:
    print("exit", exc.code)
prof.dump_stats(out)
