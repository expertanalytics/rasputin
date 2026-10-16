import re, sys
def body(f):
    b = open(f, "rb").read(); i = b.index(b"\nPOINTS "); return b[i:]
for x, y in zip(sys.argv[1::2], sys.argv[2::2]): print(x.split('/')[-1], y.split('/')[-1], "geometry and cells bit for bit:", body(x) == body(y))
