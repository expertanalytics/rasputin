"""Probe for README.md (basin-phases); a measurement script, not production code.
Run from this directory with the venv that has build-bench/pkg on its path."""
import ctypes, os, numpy as np
_libc, _buf = ctypes.CDLL("/usr/lib/libSystem.B.dylib"), ctypes.create_string_buffer(512)
def fp():
    _libc.proc_pid_rusage(os.getpid(), 4, _buf); return round(int.from_bytes(_buf.raw[72:80], "little")/1e9, 2)
a = np.ones((20000, 25000), np.float32); b = np.empty_like(a); b[:] = a
print("alive", fp(), end="; ")
del a, b
print("freed", fp(), end="; ")
_libc.malloc_zone_pressure_relief(None, 0); print("relief", fp())
