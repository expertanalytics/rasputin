"""summarize.py DOMAIN label... -> Markdown table rows from logs/."""
import sys, json, re, statistics as st
B = "/private/tmp/claude-501/-Users-skavhaug-projects-rasputin/d378b7a0-c4f3-40c9-be13-20c779636779/scratchpad/bench1m/logs"
dom = sys.argv[1]
print("| label | refine s (med of 3) | app s | process s | triangles | achieved max err m | rounds | flips | min angle median | % < 1° | n < 1° | n < 0.1° | worst ° | max degree | deg ≥ 12 | deg ≥ 20 | CDT violations (edges) |")
print("|" + "---|" * 17)
for l in sys.argv[2:]:
    v = {"refine_s": [], "app_s": [], "proc_s": []}
    for i in (1, 2, 3):
        for line in open(f"{B}/{dom}_{l}.t{i}.out"):
            m = re.match(r"BENCH (\w+) ([\d.]+)", line)
            if m and m.group(1) in v: v[m.group(1)].append(float(m.group(2)))
    out = open(f"{B}/{dom}_{l}.ascii.out").read()
    g = lambda pat: (re.search(pat, out) or [None, "?"])[1]
    q = json.loads(open(f"{B}/{dom}_{l}.q.json").read())
    med = lambda k: f"{st.median(v[k]):.3f}" if v[k] else "?"
    spread = lambda k: f" ({min(v[k]):.3f}–{max(v[k]):.3f})" if v[k] else ""
    err = g(r"achieved max error ([\d.e-]+)")
    err = f"{float(err):.6f}" if err != "?" else err
    print(f"| {l} | {med('refine_s')}{spread('refine_s')} | {med('app_s')} | {med('proc_s')} | {q['triangles']} | {err} | {g(r'(\d+) rounds')} | {g(r'(\d+) flips')} | "
          f"{q['minang_median']:.2f} | {q['pct_lt1']:.3f} | {q['n_lt1']} | {q['n_lt01']} | {q['worst']:.4f} | {q['deg_max']} | {q['deg_ge12']} | {q['deg_ge20']} | {q['cdt_violations']} / {q['cdt_edges_checked']} |")
