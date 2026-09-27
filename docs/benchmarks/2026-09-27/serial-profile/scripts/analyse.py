"""analyse.py INSTR_TXT [--chunks]: per-round tables from the instrumented run.

Record lines (instrument.patch): R round active scan_s nchunks; C begin end
start stop nodes; I col row on_edge flips | written | read; S marked
deferred_touched deferred_edge triangles split_s incircle exact exact_zero
orient orient_exact.

Footprint of an insertion = written slots (t, appended slots, u) plus every
slot legalise_around read (each popped triangle and its neighbour across the
tested edge; a flip writes exactly such a pair). Two insertions of the same
round conflict when their footprints share a slot.
"""
import sys
from collections import defaultdict

import numpy as np

path = sys.argv[1]
rounds = []
cur = None
for line in open(path):
    tag = line[0]
    if tag == "R":
        _, r, active, scan_s, nch = line.split()
        cur = dict(round=int(r), active=int(active), scan_s=float(scan_s), chunks=[], ins=[])
        rounds.append(cur)
    elif tag == "C":
        _, b, e, s0, s1, nodes = line.split()
        cur["chunks"].append((int(b), int(e), float(s0), float(s1), int(nodes)))
    elif tag == "I":
        head, written, read, fw = line[2:].split("|")
        col, row, on_edge, flips = head.split()
        cur["ins"].append((float(col), float(row), int(on_edge), int(flips),
                           [int(x) for x in written.split()], [int(x) for x in read.split()],
                           [int(x) for x in fw.split()]))
    elif tag == "S":
        f = line.split()[1:]
        cur.update(marked=int(f[0]), def_touched=int(f[1]), def_edge=int(f[2]), tris=int(f[3]),
                   split_s=float(f[4]), inc=int(f[5]), inc_exact=int(f[6]), inc_zero=int(f[7]),
                   ori=int(f[8]), ori_exact=int(f[9]))

if "--chunks" in sys.argv:
    print("| round | active | nodes scanned | scan ms | chunks | max chunk ms | mean chunk ms | imbalance max/mean | last start ms | join tail ms |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    tot = defaultdict(float)
    for r in rounds:
        ch = r["chunks"]
        durs = [c[3] - c[2] for c in ch]
        nodes = sum(c[4] for c in ch)
        mx, mean = max(durs), sum(durs) / len(durs)
        last_start = max(c[2] for c in ch)
        tail = r["scan_s"] - max(c[3] for c in ch)
        tot["scan"] += r["scan_s"]; tot["max"] += mx; tot["mean"] += mean
        tot["start"] += last_start; tot["tail"] += tail; tot["nodes"] += nodes
        print(f"| {r['round']} | {r['active']} | {nodes} | {1e3*r['scan_s']:.2f} | {len(ch)} | {1e3*mx:.2f} | {1e3*mean:.2f} | {mx/mean:.2f} | {1e3*last_start:.3f} | {1e3*tail:.3f} |")
    print(f"| sum | | {tot['nodes']:.0f} | {1e3*tot['scan']:.1f} | | {1e3*tot['max']:.1f} | {1e3*tot['mean']:.1f} | {tot['max']/tot['mean']:.2f} | {1e3*tot['start']:.2f} | {1e3*tot['tail']:.2f} |")
    sys.exit(0)

BLOCK = 64  # nodes; a coarse partition only to describe spatial spread
print("| round | active | marked | inserted | deferred (touched) | deferred (edge) | flips/ins mean | footprint mean | p50 | p99 | max | write set mean | ins. conflicting (fp) | conflict pairs (fp) | mean/max conflict degree (fp) | ins. conflicting (write-write) | earlier-overlap (fp) | median NN dist (nodes) | 64-blocks hit | max ins/block |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
T = defaultdict(float)
allfp = []
allw = []
for r in rounds:
    ins = r["ins"]
    if not ins:
        print(f"| {r['round']} | {r['active']} | {r['marked']} | 0 |" + " |" * 16)
        continue
    fps = [set(w) | set(rd) for (_, _, _, _, w, rd, _) in ins]
    wps = [set(w) | set(fw) for (_, _, _, _, w, _, fw) in ins]
    sizes = np.array([len(f) for f in fps])
    allfp += list(sizes)
    flips = np.array([x[3] for x in ins])
    def conflicts(sets):
        owners = defaultdict(list)
        for k, f in enumerate(sets):
            for x in f:
                owners[x].append(k)
        pairs = set()
        for lst in owners.values():
            for i1 in range(len(lst)):
                for i2 in range(i1 + 1, len(lst)):
                    pairs.add((lst[i1], lst[i2]))
        deg = np.zeros(len(sets), int)
        for i1, i2 in pairs:
            deg[i1] += 1
            deg[i2] += 1
        return pairs, deg
    pairs, deg = conflicts(fps)
    wpairs, wdeg = conflicts(wps)
    involved = deg > 0
    wsizes = np.array([len(w) for w in wps])
    # overlap with an earlier insertion of the round, in the serial order
    seen, earlier = set(), 0
    for f in fps:
        earlier += 1 if f & seen else 0
        seen |= f
    pts = np.array([(x[0], x[1]) for x in ins])
    # nearest-neighbour distance by a grid hash (no scipy in the venv)
    cell = 16.0
    grid = defaultdict(list)
    for k, (c, rw) in enumerate(pts):
        grid[(int(c // cell), int(rw // cell))].append(k)
    nn = np.full(len(pts), np.inf)
    for k, (c, rw) in enumerate(pts):
        gx, gy = int(c // cell), int(rw // cell)
        best = np.inf
        for ring in range(0, 400):
            # every point in rings > ring is at least ring * cell away
            if best <= ring * cell:
                break
            cand = [j for dx in range(-ring, ring + 1) for dy in range(-ring, ring + 1)
                    if max(abs(dx), abs(dy)) == ring for j in grid.get((gx + dx, gy + dy), []) if j != k]
            if cand:
                best = min(best, float(np.min(np.hypot(pts[cand, 0] - c, pts[cand, 1] - rw))))
        nn[k] = best
    blocks = defaultdict(int)
    for c, rw in pts:
        blocks[(int(c // BLOCK), int(rw // BLOCK))] += 1
    n = len(ins)
    T["marked"] += r["marked"]; T["ins"] += n; T["dt"] += r["def_touched"]; T["de"] += r["def_edge"]
    T["pairs"] += len(pairs); T["winv"] += (wdeg > 0).sum(); T["wpairs"] += len(wpairs); allw.extend(wsizes); T["inv"] += involved.sum(); T["earlier"] += earlier
    print(f"| {r['round']} | {r['active']} | {r['marked']} | {n} | {r['def_touched']} | {r['def_edge']} | {flips.mean():.2f} | {sizes.mean():.1f} | {np.percentile(sizes,50):.0f} | {np.percentile(sizes,99):.0f} | {sizes.max()} | {wsizes.mean():.1f} | {100*involved.mean():.0f} % | {len(pairs)} | {deg.mean():.1f}/{deg.max()} | {100*(wdeg>0).mean():.0f} % | {100*earlier/n:.0f} % | {np.median(nn[np.isfinite(nn)]) if np.isfinite(nn).any() else float('nan'):.1f} | {len(blocks)} | {max(blocks.values())} |")
allfp = np.array(allfp)
print(f"\nTotals: marked {T['marked']:.0f}, inserted {T['ins']:.0f} ({100*T['ins']/T['marked']:.0f} %), "
      f"deferred touched {T['dt']:.0f} ({100*T['dt']/T['marked']:.0f} %), deferred edge {T['de']:.0f} ({100*T['de']/T['marked']:.0f} %); "
      f"insertions in >=1 conflict {T['inv']:.0f} ({100*T['inv']/T['ins']:.0f} %), conflict pairs {T['pairs']:.0f}, "
      f"write-write: insertions in >=1 conflict {T['winv']:.0f} ({100*T['winv']/T['ins']:.0f} %), pairs {T['wpairs']:.0f}; "
      f"overlap an earlier insertion {T['earlier']:.0f} ({100*T['earlier']/T['ins']:.0f} %)")
allw = np.array(allw)
print(f"write set slots over all insertions: mean {allw.mean():.2f}, p50 {np.percentile(allw,50):.0f}, p90 {np.percentile(allw,90):.0f}, p99 {np.percentile(allw,99):.0f}, max {allw.max()}")
print(f"footprint slots over all insertions: mean {allfp.mean():.2f}, p50 {np.percentile(allfp,50):.0f}, p90 {np.percentile(allfp,90):.0f}, p99 {np.percentile(allfp,99):.0f}, max {allfp.max()}")
inc = sum(r.get("inc", 0) for r in rounds); ie = sum(r.get("inc_exact", 0) for r in rounds)
iz = sum(r.get("inc_zero", 0) for r in rounds); oc = sum(r.get("ori", 0) for r in rounds); oe = sum(r.get("ori_exact", 0) for r in rounds)
fl = sum(x[3] for r in rounds for x in r["ins"])
print(f"split-phase predicates: incircle {inc} (exact fallback {ie} = {100*ie/inc:.1f} %, of which cocircular {iz}), orient2d {oc} (exact {oe}); flips {fl}")
