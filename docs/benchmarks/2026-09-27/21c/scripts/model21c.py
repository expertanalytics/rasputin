"""model21c.py SIMDIR SYNC_TXT -- the scaling model of README "A0 and the scaling model".

A MODEL, not a measurement. From the per-step work the simulations recorded
(A0/A1: marks evaluated and winners per sub-round; C: marks per batch, candidate
edges and flips per flip round), it predicts the split phase at T threads as

    sum over steps of  ceil(items / T) * per-item cost  +  syncs * sync cost

Per-item costs are serial measurements: the simulation's own phase timers
(`timers` in each JSON, sim code, 1 thread) and, for a commit (split plus
Lawson), today's production split phase per insertion (21c section 6).
Sync costs come from scripts/sync_bench.cpp (SYNC_TXT, median of the runs'
medians): "spawn" = start and join T fresh std::jthreads per step, "pool" = one
std::barrier wait per step in a persistent team.
Writes the tables as Markdown to stdout.
"""
import json
import math
import re
import statistics
import sys
from pathlib import Path

sim_dir, sync_txt = Path(sys.argv[1]), Path(sys.argv[2])

# 21c section 6, battery, 7f688aa: today's serial split (1 thread) and the
# 8-thread scan and rest, ms.
SPLIT1 = {"quarter": 92.6, "tile": 92.6}
SCAN8 = {"quarter": 47.5, "tile": 60.6}
REST8 = {"quarter": 16.7, "tile": 27.9}
THREADS = (1, 4, 8, 16)
SYNCS_PER_SUBROUND = 3   # A0/A1: eval+reserve | winner check + count | prefix sum + commit
SYNCS_PER_BATCH = 2      # C: plan + ownership reduce | prefix sum + split
SYNCS_PER_FLIPROUND = 3  # C: test candidates | select by min key | flip

# sync_bench.txt: blocks starting "threads N", lines "name median X us ..."
sync: dict[int, dict[str, list[float]]] = {}
cur = 0
for line in sync_txt.read_text().splitlines():
    if m := re.match(r"threads (\d+)", line):
        cur = int(m.group(1))
    elif m := re.match(r"(\w+)\s+median\s+([\d.]+) us", line):
        sync.setdefault(cur, {}).setdefault(m.group(1), []).append(float(m.group(2)))
SYNC_US = {t: {k: statistics.median(v) for k, v in d.items()} for t, d in sync.items()}


def load(name: str) -> dict:
    return json.loads((sim_dir / f"{name}.json").read_text())


def steps(dom: str, opt: str, costs: dict) -> list[tuple[int, float, int]]:
    """(items, per-item cost in s, syncs) per step."""
    once = opt.endswith("-once")  # lower bound: each mark evaluated once, losers not re-evaluated
    opt = opt.removesuffix("-once")
    d = load(f"{dom}_{opt}")["sim"]
    out: list[tuple[int, float, int]] = []
    if opt in ("A0", "A1"):
        for r in d["rounds"]:
            for k, (pend, win) in enumerate(zip(r.get("sr_pending", []), r.get("subround_winners", []))):
                out.append((pend if not once or k == 0 else 0, costs["eval"], 1))
                out.append((win, costs["commit"], SYNCS_PER_SUBROUND - 1))
        if opt == "A0":  # the end-of-round renumbering, whole mesh, parallel over T
            per_round = d["timers"]["t_renum"] / len(d["rounds"])
            for _ in d["rounds"]:
                out.append((1_000_000, per_round / 1_000_000, 1))
    else:
        t = d["timers"]
        recorded = sum(sum(r.get("fr_cand", [])) for r in d["rounds"])
        unrec = t["n_ctest"] - recorded  # the last, flip-free test round of each refine round
        rounds_with_flips = [r for r in d["rounds"] if r.get("fr_cand")]
        for r in d["rounds"]:
            if r["marked"]:
                out.append((r["marked"], costs["cbatch"], SYNCS_PER_BATCH))
            for c, f in zip(r.get("fr_cand", []), r.get("fr_flips", [])):
                out.append((c, costs["ctest"], 1))
                out.append((f, costs["cflip"], SYNCS_PER_FLIPROUND - 1))
            if r.get("fr_cand"):
                out.append((math.ceil(unrec / len(rounds_with_flips)), costs["ctest"], 1))
    return out


def predict(st: list[tuple[int, float, int]], threads: int, sync_us: float) -> float:
    """ms"""
    s = sum(math.ceil(n / threads) * c + k * sync_us * 1e-6 for n, c, k in st)  # a step with no items still syncs
    return s * 1e3


print("| domain | option | per-item costs (ns, serial) | steps | syncs | T | split, no sync | split, spawn | split, pool | refine at 8, spawn | refine at 8, pool |")
print("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
for dom in ("quarter", "tile"):
    today = load(f"{dom}_today")
    ins_today = today["child"]["inserted"]
    for opt in ("A0", "A1", "C", "A0-once", "A1-once"):
        d = load(f"{dom}_{opt.removesuffix('-once')}")["sim"]
        t = d["timers"]
        # Evaluation cost measured on the A0 run of the same domain (same routine as A1's).
        a0t = load(f"{dom}_A0")["sim"]["timers"]
        costs = {
            "eval": a0t["t_eval"] / a0t["n_eval"],
            "commit": SPLIT1[dom] * 1e-3 / ins_today,
        }
        if opt == "C":
            costs = {"cbatch": t["t_cbatch"] / t["n_cbatch"], "ctest": t["t_ctest"] / t["n_ctest"],
                     "cflip": t["t_cflip"] / t["n_cflip"]}
        st = steps(dom, opt, costs)
        nsync = sum(k for n, c, k in st)
        cs = ", ".join(f"{k} {v * 1e9:.0f}" for k, v in costs.items())
        for T in THREADS:
            s0 = predict(st, T, 0.0)
            sp = predict(st, T, SYNC_US.get(T, {}).get("spawn", 0.0)) if T > 1 else s0
            sb = predict(st, T, SYNC_US.get(T, {}).get("barrier", 0.0)) if T > 1 else s0
            r8 = (f"{SCAN8[dom] + REST8[dom] + sp:.0f}", f"{SCAN8[dom] + REST8[dom] + sb:.0f}") if T == 8 else ("", "")
            print(f"| {dom} | {opt} | {cs if T == 1 else ''} | {len([s for s in st if s[0]]) if T == 1 else ''} | "
                  f"{nsync if T == 1 else ''} | {T} | {s0:.1f} | {sp:.1f} | {sb:.1f} | {r8[0]} | {r8[1]} |")
print()
print("Sync costs used (us, median of runs' medians): " + "; ".join(
    f"T={t}: spawn {v.get('spawn', 0):.1f}, for_each_block {v.get('feb', 0):.1f}, "
    f"std::barrier {v.get('barrier', 0):.2f}, spin {v.get('spin', 0):.2f}" for t, v in sorted(SYNC_US.items())))
print(f"Today, measured: split {SPLIT1['quarter']} ms (quarter) / {SPLIT1['tile']} ms (tile) at 1 thread; "
      f"refine at 8 threads 161.7 / 186.3 ms (21c section 6).")
