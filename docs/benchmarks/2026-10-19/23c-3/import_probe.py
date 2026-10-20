"""D8's two facts, per package: the time to import tin_engine.cli in a fresh
process, and which cut-path modules are loaded after the bench's uncut quarter
run (--domain, --tolerance 1). One line to stdout (@perf).
python import_probe.py PKG DEM QUARTER OUT"""
import sys
import time

pkg, dem, quarter, out = sys.argv[1:5]
sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
sys.path.insert(0, pkg)
t0 = time.perf_counter()
import tin_engine.cli as cli  # noqa: E402

imp = time.perf_counter() - t0
assert cli.__file__.startswith(pkg)
cut = ("tin_engine.basin_run", "tin_engine.pieces", "tin_engine.io.mesh_index")
before = [m for m in cut if m in sys.modules]
cli.app(args=["mesh", "--dem", dem, "--domain", quarter, "--tolerance", "1", "--out", out, "--binary"],
        prog_name="rasputin", standalone_mode=False)  # fmt: skip
after = [m for m in cut if m in sys.modules]
print(f"import_s {imp:.4f} cut_modules_after_import {before} after_run {after}")
