"""Medians over the repeats of each case, from the --stats reports in stats/ and
the peak RSS in logs/*.log. Prints one Markdown table per group."""
import re, statistics as S, sys
from collections import defaultdict
from pathlib import Path

D = Path(sys.argv[1])
rss = defaultdict(list)
for log in sorted((D / "logs").glob("*.log")):
    case = None
    for line in log.read_text().splitlines():
        if m := re.match(r"### (\S+) r(\d+)", line):
            case = m.group(1)
        elif "maximum resident set size" in line and case:
            rss[case].append(int(line.split()[0]))

def cells(md, header):
    lines = md.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("| " + header))
    return dict(zip([c.strip() for c in lines[i].strip("|").split("|")], [c.strip() for c in lines[i + 2].strip("|").split("|")]))

runs = defaultdict(list)
for f in sorted((D / "stats").glob("*.md")):
    case = re.sub(r"-r\d+\.md$", "", f.name)
    md = f.read_text()
    phases = {m.group(1): float(m.group(2)) for m in re.finditer(r"^\| ([^|*]+?) \| ([0-9.]+) \|", md, re.M)}
    total = float(re.search(r"\*\*total\*\* \| \*\*([0-9.]+)\*\*", md).group(1))
    size = {m.group(1): m.group(2) for m in re.finditer(r"^\| ([a-z ]+) \| (\d+) \|$", md, re.M)}
    q = cells(md, "metric | median | < 1°")
    deg = cells(md, "metric | median | p99")
    ref = cells(md, "tolerance | achieved")
    runs[case].append(dict(phases=phases, total=total, tri=int(size["output triangles"]),
        start=int(size["start vertices"]), q=q, deg=deg, ref=ref))

def med(case, key):
    vals = [r["phases"].get(key, 0.0) for r in runs[case]]
    return S.median(vals)

print("| case | n | triangles | start vertices | min angle median | < 1° | worst | max degree | achieved max error | quality inserted | feet | features read | features clip | node | start quality | refine (scan / split) | encode | total | peak RSS |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|")
for case in sorted(runs, key=lambda c: (c.split("-")[0], c)):
    rs = runs[case]; r0 = rs[0]
    same = all(r["tri"] == r0["tri"] for r in rs)
    def t(k): return f"{med(case, k):.3f}"
    print(f"| {case} | {len(rs)} | {r0['tri']}{'' if same else ' (differs)'} | {r0['start']} | {r0['q']['median']} | {r0['q']['< 1°']} | {r0['q']['worst']} | {r0['deg']['max']} | {r0['ref']['achieved max error']} | "
          f"{r0['ref'].get('quality inserted', '-')} | {r0['ref'].get('feet', '-')} | {t('features read')} | {t('features clip')} | {t('start mesh: node')} | {t('refine: start quality')} | "
          f"{t('refine')} ({t('refine: scan (parallel)')} / {t('refine: split + flip (serial)')}) | {t('write: encode')} | {S.median(r['total'] for r in rs):.2f} | "
          f"{S.median(rss[case]) / 1e9:.2f} GB |")
