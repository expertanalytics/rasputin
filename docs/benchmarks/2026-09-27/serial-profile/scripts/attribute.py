"""attribute.py SAMPLE_TXT DSYM : self-sample attribution of a sample(1) report.

Parses sample(1)'s call graph, computes each node's self count (its count
minus its children's), and resolves every _core address with `atos -i`
against the dSYM to its inline chain (innermost first, with file:line). Each
self sample inside the refine call is put in one category: frames are walked
innermost to outermost (inline chain, then the sampled stack) and the first
frame that matches a rule decides. Allocation (malloc/free/memmove) is its
own category, and is also charged to the next matching frame as its owner.
Samples not under refine (Python, idle threads) are dropped.
"""
import re
import subprocess
import sys
from collections import Counter, defaultdict

path, dsym = sys.argv[1], sys.argv[2]
text = open(path).read()
lines = text.split("\n")
start = next(i for i, l in enumerate(lines) if l.startswith("Call graph:")) + 1
end = next(i for i, l in enumerate(lines[start:], start) if not l.strip())
load = int(re.search(r"(0x[0-9a-f]+) -\s+0x[0-9a-f]+ \+_core\.cpython", text).group(1), 16)

node_re = re.compile(r"^([ +!:|]*)(\d+) (.*)$")
nodes = []  # depth, count, name, addrs, parent, in_core, self
stack = []
for l in lines[start:end]:
    mm = node_re.match(l)
    if not mm:
        continue
    depth, count, rest = len(mm.group(1)), int(mm.group(2)), mm.group(3)
    am = re.search(r"\[([0-9a-fx,]+)\]", rest)
    addrs = [int(a, 16) for a in am.group(1).split(",")] if am else []
    while stack and nodes[stack[-1]][0] >= depth:
        stack.pop()
    parent = stack[-1] if stack else None
    nodes.append([depth, count, rest.split("  (in ")[0], addrs, parent, "(in _core." in rest, 0])
    stack.append(len(nodes) - 1)
child_sum = defaultdict(int)
for n in nodes:
    if n[4] is not None:
        child_sum[n[4]] += n[1]
for i, n in enumerate(nodes):
    n[6] = n[1] - child_sum[i]

chains = {}
for a in sorted({a for n in nodes if n[5] for a in n[3]}):
    out = subprocess.run(["atos", "-i", "-o", dsym, "-l", hex(load), hex(a)],
                         capture_output=True, text=True).stdout
    chains[a] = [x.strip() for x in out.split("\n") if x.strip()]


def frames_of(i):
    """Every frame from node i to the root, innermost first: a _core node
    contributes the full inline chain of its first address."""
    out = []
    while i is not None:
        n = nodes[i]
        out += chains.get(n[3][0], [n[2]]) if n[5] and n[3] else [n[2]]
        i = n[4]
    return out


ALLOC = r"malloc|_free\b|free_tc|memmove|memset|bzero|operator new|operator delete|_xzm|madvise|vm_"
RULES = [
    ("exact incircle (adaptive fallback)", r"incircleadapt|expansion_|DetriaExact::incircle"),
    ("filtered predicates (incircle/orient2d)", r"FilteredKernel|orient2d|incircle|orient_sign"),
    ("must_flip (edge lookup, frame points)", r"must_flip|lawson\.hpp:(5\d|6\d|7\d|8\d)\)"),
    ("topology writes: flip", r"LatticeMesh::flip"),
    ("topology writes: split", r"split_inside|split_edge|add_vertex"),
    ("topology writes: put/repoint (callers above)", r"LatticeMesh::put|repoint"),
    ("legalise_around loop (stack, touched marks)", r"legalise_around|lawson\.hpp:1[0-2]\d\)"),
    ("legalise_all (start mesh)", r"legalise_all"),
    ("start quality (improve)", r"improve<|quality\.hpp"),
    ("scan (parallel phase)", r"refinement::scan<|for_each_row|row_spans|row_segments|scan\.hpp|refine\.hpp:29[0-4]\)|chunks\.hpp"),
    ("rebuild active: collect+sort+unique", r"__introsort|__bitset_partition|__partial_sort|__insertion_sort|__sort|refine\.hpp:36[0-6]\)"),
    ("split loop bookkeeping (results read, touched, skipped)", r"refine\.hpp:(30\d|31\d|32\d|33\d|34[0-6]|35\d)\)"),
    ("results.resize / round setup", r"refine\.hpp:28[6-9]\)"),
    ("refine body, no line info (refine.hpp:0)", r"refine\.hpp:0\)"),
    ("output: vertices, z, triangles, constraint_edges", r"constraint_edges|refine\.hpp:3[7-9]\d\)|bilinear"),
    ("setup: to_lattice, LatticeMesh::build", r"to_lattice|LatticeMesh::build|__hash_table|__tree|refine\.hpp:(9\d|1[0-3]\d|2[4-6]\d|28[0-4])\)"),
    ("pybind: outcome conversion", r"pybind11::|type_caster|cast_op")
]


def strip(frame):
    """A frame reduced to its function name and file:line: balanced <...>
    and (...) groups removed, so template arguments cannot match a rule."""
    loc = re.search(r"\(([\w.+-]+:\d+)\)$", frame)
    name = frame.split(" (in ")[0]
    out, depth = [], 0
    for ch in name:
        if ch in "<(":
            depth += 1
        elif ch in ">)" and depth:
            depth -= 1
        elif depth == 0:
            out.append(ch)
    return "".join(out) + (f" ({loc.group(1)})" if loc else "")


def classify(frames):
    for f in frames:
        s = strip(f)
        for name, rx in RULES:
            if re.search(rx, s):
                return name
    return "other: " + (strip(frames[0])[:80] if frames else "?")


cats, owners, leaves = Counter(), Counter(), Counter()
total = 0.0
for i, n in enumerate(nodes):
    if n[6] <= 0:
        continue
    rest = frames_of(n[4]) if n[4] is not None else []
    if not any("refinement::refine<" in f for f in rest + [n[2]]):
        continue
    addrs = n[3] if (n[5] and n[3]) else [None]
    per = n[6] / len(addrs)
    for a in addrs:
        frames = (chains.get(a, [n[2]]) if a is not None else [n[2]]) + rest
        total += per
        if re.search(ALLOC, strip(frames[0])):
            cats["allocation (malloc/free/memmove)"] += per
            owners[classify(frames[1:])] += per
        else:
            cats[classify(frames)] += per
        leaves[frames[0][:160]] += per

print(f"self samples under refine: {total:.0f} (1 sample = 1 ms)")
for c, v in cats.most_common():
    print(f"{v:8.0f}  {100 * v / total:5.1f} %  {c}")
print("\nallocation charged to its owner:")
for c, v in owners.most_common():
    print(f"{v:8.0f}  {100 * v / total:5.1f} %  {c}")
print("\ntop innermost frames:")
for c, v in leaves.most_common(40):
    print(f"{v:8.0f}  {100 * v / total:5.1f} %  {c}")
