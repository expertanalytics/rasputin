"""a0eo.py DATA SYNC_TXT PHASES -- tables and model for README "A0 with evaluate-once".

DATA is data/a0eo (run_a0eo.sh's output): today and A0 once per domain,
rep{1,2,3}/ and verify/ for the three dirty rules (A0eot, A0eo, A0eofp), and
plant_*/. SYNC_TXT is data/a0/sync_bench.txt; PHASES is
data/a0eo/split_phases_4_8_16.txt (split_sweep.sh, THREADS="4 8 16").

Everything under "Model" is a MODEL, not a measurement. The split phase at T threads is

    sum over steps of  ceil(items / T) * per-item cost  +  syncs * barrier cost

with per-item costs measured serially by the simulation's pass timers (median
of three runs), today's production split per insertion for a commit (21c
section 6), and the std::barrier cost of a persistent team (sync_bench.txt).
Refine = modelled split + the scan and rest measured at T threads (PHASES).
Writes Markdown to stdout.
"""
import json
import math
import re
import statistics
import sys
from pathlib import Path

data, sync_txt, phases_txt = (Path(a) for a in sys.argv[1:4])
RULES = ("A0eot", "A0eo", "A0eofp")
DOMS = ("quarter", "tile")
THREADS = (1, 4, 8, 16)
SPLIT1 = {"quarter": 92.6, "tile": 92.6}  # ms, today's serial split, 21c section 6 (7f688aa, battery)
SYNCS_EVAL, SYNCS_COMMIT, SYNCS_RENUM = 1, 2, 1  # 3 per sub-round, as in model21c.py

barrier: dict[int, list[float]] = {}
cur = 0
for line in sync_txt.read_text().splitlines():
    if m := re.match(r"threads (\d+)", line):
        cur = int(m.group(1))
    elif m := re.match(r"barrier\s+median\s+([\d.]+) us", line):
        barrier.setdefault(cur, []).append(float(m.group(1)))
BARRIER_US = {t: statistics.median(v) for t, v in barrier.items()}

phase: dict[tuple[str, int], list[dict]] = {}
for line in phases_txt.read_text().splitlines():
    dom, js = line.split(" ", 1)
    j = json.loads(js)
    phase.setdefault((dom, j["threads"]), []).append(j)
PH = {k: {x: statistics.median(r[x + "_s"] for r in v) * 1e3 for x in ("refine", "scan", "split", "rest")}
      for k, v in phase.items()}


def load(rel: str) -> dict:
    return json.loads((data / rel).read_text())


def rounds_equal(a: dict, b: dict, keys: tuple[str, ...]) -> int:
    """Rounds whose listed fields differ, plus any difference in round count."""
    ra, rb = a["sim"]["rounds"], b["sim"]["rounds"]
    return sum(any(x.get(k) != y.get(k) for k in keys) for x, y in zip(ra, rb)) + abs(len(ra) - len(rb))


def checks(d: dict, today: dict, a0: dict) -> str:
    s, q = d["sim"], d["quality"]
    cons = sum(r["cavity_violations"] for r in s["rounds"])
    return (f"`{d['vtk_sha256'][:8]}` | {d['vtk_triangles']:,} | {d['child']['rounds']} | {d['child']['inserted']:,} | "
            f"{d['child']['flips']:,} | {rounds_equal(d, today, ('hash',))} | "
            f"{rounds_equal(d, a0, ('subround_winners', 'sr_pending')) if d['mode'] != 'today' else '-'} | "
            f"{s['full_rescan_max_error']:.5f} / {s['full_rescan_needs_split']} | {q['delaunay_violations']} | "
            f"{cons} | {s['order_violations']}")


print("### Identity and checks (measured)\n")
print("Round hash: rounds whose lattice-mesh hash differs from today's. Sub-rounds: rounds whose "
      "per-sub-round winner and pending lists differ from A0's (same build).\n")
print("| domain | run | sha256 | triangles | rounds | inserted | flips | round hash differs | sub-rounds differ "
      "from A0 | full rescan max / over | Delaunay violations | consistency | order violations |")
print("|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|")
for dom in DOMS:
    today, a0 = load(f"{dom}_today.json"), load(f"{dom}_A0.json")
    print(f"| {dom} | today | {checks(today, today, a0)} |")
    print(f"| {dom} | A0 | {checks(a0, today, a0)} |")
    for rule in RULES:
        for sub in ("rep1", "rep2", "rep3", "verify"):
            print(f"| {dom} | {rule} {sub} | {checks(load(f'{sub}/{dom}_{rule}.json'), today, a0)} |")

print("\n### Evaluations (measured; counts are equal in all three repetitions and the verify run)\n")
print("| domain | rule | marks | pending visits / mark | full evaluations / mark | first | re-evaluated: slot "
      "dirty | re-evaluated: `footed` | clean re-bids | stale re-uses (verify) | re-evaluations that gave the "
      "same footprint (verify) | slots loaded by the dirty check / visit | version stores / commit |")
print("|---|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|")
for dom in DOMS:
    a0 = load(f"{dom}_A0.json")["sim"]
    marks = sum(r["marked"] for r in a0["rounds"])
    visits = sum(sum(r.get("sr_pending", [])) for r in a0["rounds"])
    full = visits - sum(r["skipped_touched"] for r in a0["rounds"])
    print(f"| {dom} | A0 (re-evaluate every sub-round) | {marks:,} | {visits / marks:.2f} | {full / marks:.2f} | "
          f"| | | | | | | |")
    for rule in RULES:
        e = load(f"rep1/{dom}_{rule}.json")["sim"]["eo"]
        v = load(f"verify/{dom}_{rule}.json")["sim"]["eo"]
        ins = load(f"rep1/{dom}_{rule}.json")["child"]["inserted"]
        print(f"| {dom} | {rule} | {marks:,} | {e['n_check'] / marks:.2f} | **{e['n_eval'] / marks:.2f}** | "
              f"{e['first']:,} | {e['reeval_slot']:,} | {e['reeval_footed']} | {e['rebid_clean']:,} | "
              f"{v['stale']} of {v['verified']:,} | {v['reeval_same']:,} of {v['reeval_slot']:,} | "
              f"{e['n_loads'] / e['n_check']:.1f} | {e['n_marks'] / ins:.1f} |")

print("\n### Per-item costs, serial (measured; median of 3 runs, min-max)\n")
print("| domain | rule | check: skip rule + dirty check, per visit | evaluation: plan + cavity + ring, "
      "per evaluation | bid: reserve + winner check + reset, per bid | clean re-bid = check + bid | "
      "version stores, per commit |")
print("|---|---|---|---|---|---|---|")
COST: dict[tuple[str, str], dict[str, float]] = {}


def med(xs: list[float]) -> str:
    return f"{statistics.median(xs):.0f} ({min(xs):.0f}-{max(xs):.0f})"


for dom in DOMS:
    a0t = load(f"{dom}_A0.json")["sim"]["timers"]
    for rule in RULES:
        es = [load(f"{r}/{dom}_{rule}.json")["sim"] for r in ("rep1", "rep2", "rep3")]
        ck = [s["eo"]["t_check"] / s["eo"]["n_check"] * 1e9 for s in es]
        ev = [s["eo"]["t_eval"] / s["eo"]["n_eval"] * 1e9 for s in es]
        bd = [s["eo"]["t_bid"] / s["eo"]["n_bid"] * 1e9 for s in es]
        mk = [s["eo"]["t_mark"] / s["timers"]["n_commit"] * 1e9 for s in es]
        rb = [c * s["eo"]["n_check"] / s["eo"]["n_bid"] + b for c, b, s in zip(ck, bd, es)]
        COST[(dom, rule)] = {k: statistics.median(v) * 1e-9 for k, v in
                             (("check", ck), ("eval", ev), ("bid", bd), ("mark", mk))}
        print(f"| {dom} | {rule} | {med(ck)} ns | {med(ev)} ns | {med(bd)} ns | {med(rb)} ns | {med(mk)} ns |")
    COST[(dom, "A0")] = {"eval": a0t["t_eval"] / a0t["n_eval"]}
    print(f"| {dom} | A0 (1 run, this build) | | {a0t['t_eval'] / a0t['n_eval'] * 1e9:.0f} ns per pending visit, "
          f"everything included | | | |")


def steps(dom: str, opt: str) -> list[tuple[list[tuple[int, float]], int]]:
    """Per step: ([(items, per-item cost s)], syncs)."""
    commit = SPLIT1[dom] * 1e-3 / load(f"{dom}_today.json")["child"]["inserted"]
    src = load(f"{dom}_A0.json" if opt == "A0" else f"rep1/{dom}_{opt}.json")["sim"]
    c = COST[(dom, opt)]
    out: list[tuple[list[tuple[int, float]], int]] = []
    for r in src["rounds"]:
        pend, win = r.get("sr_pending", []), r.get("subround_winners", [])
        ev, bids = r.get("sr_evals", [0] * len(pend)), r.get("sr_bids", [0] * len(pend))
        for p, w, e, b in zip(pend, win, ev, bids):
            if opt == "A0":
                out.append(([(p, c["eval"])], SYNCS_EVAL))
                out.append(([(w, commit)], SYNCS_COMMIT))
            else:
                out.append(([(p, c["check"]), (e, c["eval"]), (b, c["bid"])], SYNCS_EVAL))
                out.append(([(w, commit + c["mark"])], SYNCS_COMMIT))
    per_round = src["timers"]["t_renum"] / len(src["rounds"])
    out += [([(1_000_000, per_round / 1_000_000)], SYNCS_RENUM)] * len(src["rounds"])
    return out


def predict(st: list, t: int) -> float:
    b = BARRIER_US.get(t, 0.0) * 1e-6 if t > 1 else 0.0
    return sum(sum(math.ceil(n / t) * c for n, c in items) + k * b for items, k in st) * 1e3


print("\n### Model: split and refine at T threads, team + barrier (MODEL, not a measurement)\n")
print("Split in ms, modelled. Refine = modelled split + scan + rest measured at T threads "
      "(PHASES, battery). 'today' rows are measured, not modelled.\n")
print("| domain | T | barrier (us) | today split (measured) | today refine (measured) | A0 split | A0eot split "
      "| A0eo split | A0eofp split | A0 refine | A0eot refine | A0eo refine | A0eofp refine | syncs |")
print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for dom in DOMS:
    for t in THREADS:
        ph = PH.get((dom, t))
        sp = {o: predict(steps(dom, o), t) for o in ("A0",) + RULES}
        nsync = sum(k for _, k in steps(dom, "A0eot"))
        ref = {o: (f"{v + ph['scan'] + ph['rest']:.0f}" if ph else "") for o, v in sp.items()}
        today_split = f"{ph['split']:.1f}" if ph else f"{SPLIT1[dom]} (1 thread, section 6)"
        today_ref = f"{ph['refine']:.1f}" if ph else ""
        print(f"| {dom} | {t} | {BARRIER_US.get(t, 0.0) if t > 1 else 0:.2f} | {today_split} | {today_ref} | "
              f"{sp['A0']:.0f} | **{sp['A0eot']:.0f}** | {sp['A0eo']:.0f} | {sp['A0eofp']:.0f} | {ref['A0']} | "
              f"{("**" + ref["A0eot"] + "**") if ref["A0eot"] else ""} | {ref['A0eo']} | {ref['A0eofp']} | {nsync} |")
print("\nBarrier costs used (us, median of runs' medians): "
      + ", ".join(f"T={t}: {v:.2f}" for t, v in sorted(BARRIER_US.items())))
print("Scan and rest used (ms, medians of 15): " + "; ".join(
    f"{d} T={t}: scan {v['scan']:.1f}, rest {v['rest']:.1f}" for (d, t), v in sorted(PH.items())))

print("\n### Plants (quarter circle, A0eot, verify on)\n")
print("| plant | sha256 | round hash differs | sub-rounds differ from A0 | consistency | order violations | "
      "stale re-uses | full rescan over | Delaunay violations |")
print("|---|---|---:|---:|---:|---:|---|---:|---:|")
today, a0 = load("quarter_today.json"), load("quarter_A0.json")
for pl in ("eo_nodirty", "eo_ringless", "eo_nofooted"):
    f = data / f"plant_{pl}" / "quarter_A0eot.json"
    d = json.loads(f.read_text())
    s, e = d["sim"], d["sim"]["eo"]
    print(f"| `{pl}` | `{d['vtk_sha256'][:8]}` | {rounds_equal(d, today, ('hash',))} | "
          f"{rounds_equal(d, a0, ('subround_winners', 'sr_pending'))} | "
          f"{sum(r['cavity_violations'] for r in s['rounds']):,} | {s['order_violations']} | "
          f"{e['stale']:,} of {e['verified']:,} | {s['full_rescan_needs_split']} | {d['quality']['delaunay_violations']} |")
print("| `eo_stalecommit` | no mesh: the process crashes (see plant_eo_stalecommit/outcome.txt) | | | | | | | |")

print("\n### Where the modelled split goes, A0eot at T = 8 (MODEL)\n")
print("| domain | check | evaluation | bid | commit + version stores | renumbering | barriers | total |")
print("|---|---:|---:|---:|---:|---:|---:|---:|")
for dom in DOMS:
    src = load(f"rep1/{dom}_A0eot.json")["sim"]
    c = COST[(dom, "A0eot")]
    commit = SPLIT1[dom] * 1e-3 / load(f"{dom}_today.json")["child"]["inserted"]
    t, parts = 8, {"check": 0.0, "eval": 0.0, "bid": 0.0, "commit": 0.0, "renum": 0.0, "sync": 0.0}
    for r in src["rounds"]:
        for p, w, e, b in zip(r.get("sr_pending", []), r.get("subround_winners", []), r.get("sr_evals", []),
                              r.get("sr_bids", [])):
            parts["check"] += math.ceil(p / t) * c["check"]
            parts["eval"] += math.ceil(e / t) * c["eval"]
            parts["bid"] += math.ceil(b / t) * c["bid"]
            parts["commit"] += math.ceil(w / t) * (commit + c["mark"])
            parts["sync"] += (SYNCS_EVAL + SYNCS_COMMIT) * BARRIER_US[t] * 1e-6
        parts["renum"] += src["timers"]["t_renum"] / len(src["rounds"]) / t + SYNCS_RENUM * BARRIER_US[t] * 1e-6
    print(f"| {dom} | " + " | ".join(f"{v * 1e3:.1f}" for v in parts.values()) + f" | {sum(parts.values()) * 1e3:.1f} |")
