"""Cut the first N bytes of a GML FeatureCollection at a featureMember end."""

import sys

src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(src, "rb") as f:
    data = f.read(n)
start = data.index(b"<gml:featureMember>")
end = data.rindex(b"</gml:featureMember>") + len(b"</gml:featureMember>")
with open(dst, "wb") as out:
    out.write(data[:end] + b"\n</gml:FeatureCollection>\n")
print(start, end, data[:end].count(b"<gml:featureMember>"))
