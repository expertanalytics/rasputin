"""Summarise a macOS `sample` call graph under one function (@perf).
python sample_tree.py SAMPLE.txt PATTERN [DEPTH] [MIN_SHARE]
Finds the first node whose name contains PATTERN and prints its subtree to
DEPTH levels: each node's samples and its share of PATTERN's samples. Names are
cut to the function and shortened (templates elided)."""
import re
import sys

path, pat = sys.argv[1], sys.argv[2]
depth = int(sys.argv[3]) if len(sys.argv) > 3 else 3
min_share = float(sys.argv[4]) if len(sys.argv) > 4 else 0.01
line_re = re.compile(r"^([ +!:|]*?)(\d+) (.*)$")
nodes = []  # (column, count, name)
started = False
for raw in open(path, encoding="utf-8", errors="replace"):
    if raw.startswith("Call graph:"):
        started = True
        continue
    if started and raw.strip() == "":
        if nodes:
            break
        continue
    if not started:
        continue
    m = line_re.match(raw.rstrip("\n"))
    if m:
        nodes.append((len(m.group(1)), int(m.group(2)), m.group(3)))


def short(name: str) -> str:
    name = re.sub(r"\s+\(in [^)]*\).*$", "", name)
    out, level = [], 0
    for ch in name:  # drop template arguments
        if ch == "<":
            level += 1
            if level == 1:
                out.append("<…>")
            continue
        if ch == ">" and level:
            level -= 1
            continue
        if level == 0:
            out.append(ch)
    s = "".join(out)
    s = re.sub(r"\((?:[^()]|\([^()]*\))*\)", "()", s)  # drop parameter lists
    return s[:150]


# Every node matching PATTERN at the shallowest column (all call sites).
roots = [i for i, (c, _, n) in enumerate(nodes) if pat in n]
col0 = min(nodes[i][0] for i in roots)
roots = [i for i in roots if nodes[i][0] == col0]
total = sum(nodes[i][1] for i in roots)
print(f"{pat}: {total} samples in {len(roots)} call site(s)")
agg: dict[tuple, int] = {}
for r in roots:
    stack = [(nodes[r][0], ())]
    for c, n, name in nodes[r + 1 :]:
        if c <= nodes[r][0]:
            break
        while stack and c <= stack[-1][0]:
            stack.pop()
        key = stack[-1][1] + (short(name),)
        if len(key) <= depth:
            agg[key] = agg.get(key, 0) + n
        stack.append((c, key))
for key in sorted(agg, key=lambda k: (k[:-1], -agg[k])):
    pass
def show(prefix: tuple, level: int) -> None:
    kids = sorted((k for k in agg if len(k) == level + 1 and k[:level] == prefix), key=lambda k: -agg[k])
    for k in kids:
        if agg[k] / total >= min_share:
            print(f"{'  ' * level}{agg[k]:6d} {100 * agg[k] / total:5.1f} %  {k[-1]}")
            show(k, level + 1)
show((), 0)
