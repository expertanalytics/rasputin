"""Increment 45, design review round 2 probe: the quick check's input rule
compared per case. Sketches of three rules against today's bench_quick and the
stored baseline; not product code. Run: .venv/bin/python -I probe45_r2.py <tree>.

- today: the union of both records' input names (bench_quick._mismatch);
- round 1: the new record's input names only;
- round 2: per case that ran and is in the baseline, the union of the two
  records' input lists for that case, each name present in both lists with
  the same fingerprint (inputs stored per case in the record, CaseResult.inputs);
- option A (not taken): per case that ran and is in the baseline, today's
  cases.toml list for that case, each name present in the baseline with the
  same fingerprint."""

import json
import sys
from pathlib import Path

sys.path.insert(0, sys.argv[1] + "/tools")
import bench_quick as q
from pydantic import Field


class CaseResult2(q.CaseResult):
    inputs: list[str] = Field(default_factory=list)


class QuickRecord2(q.QuickRecord):
    cases: dict[str, CaseResult2]  # type: ignore[assignment]


def machine_power(new, base):
    if new.power not in ("ac", "battery") or new.power != base.power:
        return f"power: {new.power} vs {base.power}"
    old_machine = base.machine.model_dump()
    for key, value in new.machine.model_dump().items():
        if value != old_machine[key]:
            return f"machine {key}: {value} vs {old_machine[key]}"
    return None


def mismatch_r1(new, base):
    if (reason := machine_power(new, base)) is not None:
        return reason
    for key in sorted(new.inputs):
        if key in base.inputs and new.inputs[key] != base.inputs[key]:
            return f"input {key} changed"
    return None


def mismatch_r2(new, base):
    if (reason := machine_power(new, base)) is not None:
        return reason
    for name in sorted(new.cases.keys() & base.cases.keys()):
        mine, theirs = set(new.cases[name].inputs), set(base.cases[name].inputs)
        for key in sorted(mine | theirs):
            if key not in mine & theirs or new.inputs.get(key) != base.inputs.get(key):
                return f"input {key} of {name} changed"
    return None


def mismatch_a(new, base, by_case):
    if (reason := machine_power(new, base)) is not None:
        return reason
    for name in sorted(new.cases.keys() & base.cases.keys()):
        for key in sorted(by_case[name.split(" ")[0]]):
            if new.inputs.get(key) != base.inputs.get(key):
                return f"input {key} of {name} changed"
    return None


def judged(rule, new, base):
    today = q._mismatch
    q._mismatch = rule
    try:
        v = q.verdict(new, base)
    finally:
        q._mismatch = today
    return f"{v.status}, exit {v.exit_code}"


tree = Path(sys.argv[1])
stored = json.loads((tree / "docs/benchmarks/quick/baseline-ac.json").read_text())
old = QuickRecord2.model_validate(stored)
print("== the stored baseline under the round 2 model (CaseResult.inputs defaults to [])")
print("  loads:", type(old).__name__, "| cases:", list(old.cases))
print("  inputs per case as loaded:", {k: v.inputs for k, v in old.cases.items()})

by_case = {c.name: c.inputs for c in q.load_cases(tree / "docs/benchmarks/quick/cases.toml")}
# The baseline as a round 2 save would write it: each case's cases.toml list.
for name, result in stored["cases"].items():
    result["inputs"] = by_case[name.split(" ")[0]]
base = QuickRecord2.model_validate(stored)
print("  a round 2 save of it holds per case, e.g. quarter:", base.cases["quarter"].inputs)

QUARTER = "$TREE/docs/benchmarks/2026-09-26/quarter.geojson"
OTHER = "$TREE/docs/benchmarks/2026-09-26/quarter_v2.geojson"
TIF = "$TREE/tests/fixtures/dem_archive/7908_3_10m_z33.tif"
in_tree = [n for n, i in by_case.items() if all("$TREE" in x for x in i)]


def run_of(keep, lists=None, changed=None):
    """A new record: the cases in ``keep`` ran, with ``lists`` overriding a
    case's cases.toml inputs, and ``changed`` inputs refingerprinted."""
    lists = {**by_case, **(lists or {})}
    inputs = {i: base.inputs.get(i, "file:1:new") for n in keep for i in lists[n]}
    for key in changed or ():
        inputs[key] = "file:0:changed"
    cases = {
        k: v.model_copy(update={"inputs": lists[k.split(" ")[0]]})
        for k, v in base.cases.items()
        if k.split(" ")[0] in keep
    }
    return base.model_copy(update={"inputs": inputs, "cases": cases}), lists


print("== the input rules, against the baseline as a round 2 save writes it")
for label, keep, lists, changed in [
    ("nothing skipped", list(by_case), None, None),
    ("only the two $TREE cases ran (the skip)", in_tree, None, None),
    ("nothing skipped, the DEM tile's content changed", list(by_case), None, [TIF]),
    (
        "quarter's .geojson replaced by another file in cases.toml (the reviewer's case)",
        list(by_case),
        {"quarter": [TIF, OTHER]},
        None,
    ),
    (
        "the same, with only the two $TREE cases run",
        in_tree,
        {"quarter": [TIF, OTHER]},
        None,
    ),
    (
        "quarter's .geojson dropped from its inputs in cases.toml",
        list(by_case),
        {"quarter": [TIF]},
        None,
    ),
]:
    new, lists = run_of(keep, lists, changed)
    print(f"  {label}: ran {keep}")
    print(f"    today:    {q._mismatch(new, base)}  ->  {judged(q._mismatch, new, base)}")
    print(f"    round 1:  {mismatch_r1(new, base)}  ->  {judged(mismatch_r1, new, base)}")
    print(f"    round 2:  {mismatch_r2(new, base)}  ->  {judged(mismatch_r2, new, base)}")
    print(f"    option A: {mismatch_a(new, base, lists)}")

print("== a new run against the stored baseline as it is (no per-case lists)")
new, _ = run_of(list(by_case))
print(f"  round 2: {mismatch_r2(new, old)}  ->  {judged(mismatch_r2, new, old)}")
