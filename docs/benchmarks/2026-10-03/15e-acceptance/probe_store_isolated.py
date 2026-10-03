"""probe_store.py's store (same seed, rows, columns and points), with every
input batch generated *before* the baseline reading, so the footprint rise
over add + freeze is the store's alone: no freed input array is left in the
process's charge inside the measured interval. A measurement script, not
production code. Run with the venv that has build-bench/pkg on its path."""
import ctypes, os, numpy as np
from tin_engine._core import CheckPoints
_libc, _buf = ctypes.CDLL("/usr/lib/libSystem.B.dylib"), ctypes.create_string_buffer(512)
def fp() -> float:
    _libc.proc_pid_rusage(os.getpid(), 4, _buf)
    return int.from_bytes(_buf.raw[72:80], "little")
R = C = 10000; N = 100_000_000; h = 30.0
rng = np.random.default_rng(1)
batches = []
for k in range(20):
    xy = rng.uniform(0, (C - 1) * h, size=(N // 20, 2)); z = rng.uniform(0, 1000, N // 20).astype(np.float32)
    batches.append((xy, z))
s = CheckPoints(x_min=0.0, y_max=R * h, spacing=h, rows=R, cols=C)
f0 = fp()
for xy, z in batches:
    s.add(xy, z)
f1 = fp()
s.freeze()
f2 = fp()
n = s.size
print(f"points {n}  baseline fp {f0/1e9:.3f} GB")
print(f"after add: store {(f1-f0)/1e9:.3f} GB = {(f1-f0)/n:.2f} B/point")
print(f"after freeze: store {(f2-f0)/1e9:.3f} GB = {(f2-f0)/n:.2f} B/point (live 16 B/point = {16*n/1e9:.3f} GB)")
