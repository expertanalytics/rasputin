import re, sys, numpy as np
b = open(sys.argv[1], "rb").read()
m = re.search(rb"\nPOINTS (\d+) (\w+)\n", b)
n, t = int(m.group(1)), m.group(2).decode()
dt = {"double": ">f8", "float": ">f4"}[t]
p = np.frombuffer(b, dtype=dt, count=3 * n, offset=m.end()).reshape(n, 3)
for name, (x, y) in {"open1": (413602.500, 6331638.254), "open2": (413602.501, 6331638.244), "far": (413677.295, 6331664.081)}.items():
    d = np.hypot(p[:, 0] - x, p[:, 1] - y); i = d.argmin()
    print(f"{name}: nearest vertex {d[i]:.4f} m away ({p[i,0]:.3f} {p[i,1]:.3f}); vertex within 1 mm: {d[i] < 1e-3}")
