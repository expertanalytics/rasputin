import sys
sys.meta_path[:] = [f for f in sys.meta_path if "skbc" not in type(f).__module__ and "editable" not in type(f).__module__.lower()]
sys.path.insert(0, "/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/master/src_python")
for m in [k for k in sys.modules if k.startswith("tin_engine")]:
    del sys.modules[m]
import tin_engine, tin_engine._core
assert tin_engine.__file__.startswith("/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/master"), tin_engine.__file__
assert tin_engine._core.__file__.startswith("/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/66e9d031-e53c-4911-b4a9-03c4d27bffa6/scratchpad/master"), tin_engine._core.__file__
from tin_engine.plugins import main
sys.argv[0] = "rasputin"
sys.exit(main())
