"""Probe for README.md (basin-phases); a measurement script, not production code.
Run from this directory with the venv that has build-bench/pkg on its path."""
import ctypes, os, numpy as np
from tin_engine._core import CheckPoints
_libc, _buf = ctypes.CDLL("/usr/lib/libSystem.B.dylib"), ctypes.create_string_buffer(512)
def u():
    _libc.proc_pid_rusage(os.getpid(), 4, _buf); q=lambda o:int.from_bytes(_buf.raw[o:o+8],"little")/1e9
    return f"res {q(64):.2f} fp {q(72):.2f} lifemax {q(240):.2f}"
R=C=10000; N=100_000_000; h=30.0
s=CheckPoints(x_min=0.0, y_max=R*h, spacing=h, rows=R, cols=C)
print("empty", u())
rng=np.random.default_rng(1)
for k in range(20):
    xy=rng.uniform(0, (C-1)*h, size=(N//20,2)); z=rng.uniform(0,1000,N//20).astype(np.float32)
    s.add(xy, z)
del xy, z
print("added", s.size, u())
s.freeze(); print("frozen", u())
print("size", s.size)
import subprocess
out = subprocess.run(["vmmap", "-summary", str(os.getpid())], capture_output=True, text=True).stdout
print("\n".join(l for l in out.splitlines() if l.startswith(("Malloc", "Physical footprint", "REGION TYPE", "TOTAL", "DefaultMallocZone", "MallocHelperZone"))))
