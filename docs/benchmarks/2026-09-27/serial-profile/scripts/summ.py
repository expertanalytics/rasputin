import json, sys, statistics as st
from collections import defaultdict
d = defaultdict(list)
for line in open(sys.argv[1]):
    p, js = line.split(" ", 1); r = json.loads(js); d[r["threads"]].append(r)
b = None
print("| threads | n | refine s | speed-up | scan s | scan speed-up | split s | rest s | serial share |")
print("|---|---|---|---|---|---|---|---|---|")
for t in sorted(d):
    rs = d[t]; m = {k: st.median(x[k] for x in rs) for k in ("refine_s","scan_s","split_s","rest_s","legalise_s","quality_s")}
    if b is None: b = m
    ser = m["split_s"]+m["rest_s"]+m["legalise_s"]+m["quality_s"]
    print(f"| {t} | {len(rs)} | {m['refine_s']:.3f} | {b['refine_s']/m['refine_s']:.2f}x | {m['scan_s']:.3f} | {b['scan_s']/m['scan_s']:.2f}x | {m['split_s']:.3f} | {m['rest_s']:.3f} | {ser/m['refine_s']*100:.0f} % |")
assert len({(r['inserted'], r['flips'], r['triangles']) for rs in d.values() for r in rs}) == 1, "output differs"
print("identical counts across all samples:", {(r['inserted'], r['flips'], r['triangles'], r['rounds']) for rs in d.values() for r in rs})
