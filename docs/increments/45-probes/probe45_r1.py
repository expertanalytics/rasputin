"""Increment 45, design review round 1 probe: the quick check's input rule
when a case is skipped, and the Linux power reader restricted to five named
values. Sketches of the rules, run against today's bench_quick and fixture
trees; not product code. Run: .venv/bin/python -I probe45_r1.py <tree>."""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, sys.argv[1] + "/tools")
import bench_quick as q

# ------------------------------------------------ the skip and _mismatch
stored = json.loads((Path(sys.argv[1]) / "docs/benchmarks/quick/baseline-ac.json").read_text())
base = q.QuickRecord.model_validate(stored)
print("== _mismatch when a case is skipped (the stored baseline, its own inputs)")
print("  baseline cases:", list(base.cases), "| inputs:", len(base.inputs))


def mismatch_r1(new, base):
    """Round 1: the machine and power as today; inputs only of the cases that ran."""
    if new.power not in ("ac", "battery") or new.power != base.power:
        return f"power: {new.power} vs {base.power}"
    old_machine = base.machine.model_dump()
    for key, value in new.machine.model_dump().items():
        if value != old_machine[key]:
            return f"machine {key}: {value} vs {old_machine[key]}"
    for key in sorted(new.inputs):
        if key in base.inputs and new.inputs[key] != base.inputs[key]:
            return f"input {key} changed"
    return None


cases = q.load_cases(Path(sys.argv[1]) / "docs/benchmarks/quick/cases.toml")
by_case = {c.name: c.inputs for c in cases}
first = next(iter(base.inputs))
for label, keep, changed in [
    ("nothing skipped", list(by_case), None),
    ("one case skipped (the first)", list(by_case)[1:], None),
    (
        "only the two $TREE cases ran",
        [n for n, i in by_case.items() if all("$TREE" in x for x in i)],
        None,
    ),
    ("nothing skipped, one input changed", list(by_case), first),
]:
    inputs = {i: base.inputs[i] for n in keep for i in by_case[n]}
    if changed:
        inputs[changed] = "file:0:changed"
    new = base.model_copy(
        update={
            "inputs": inputs,
            "cases": {k: v for k, v in base.cases.items() if k.split(" ")[0] in keep},
        }
    )
    print(f"  {label}: ran {keep}")
    print(f"    today's _mismatch: {q._mismatch(new, base)}")
    print(f"    round 1 rule:      {mismatch_r1(new, base)}")
    today = q._mismatch
    q._mismatch = mismatch_r1
    v = q.verdict(new, base)
    q._mismatch = today
    print(f"    verdict under the round 1 rule: {v.status}, exit {v.exit_code}")

# ---------------------------------------------- the power reader, round 1
VALUES = ("type", "scope", "online", "status", "capacity")
MAINS = ("Mains", "USB", "Wireless")
NONE = "no power supply under /sys/class/power_supply"


def read_value(path):
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return None
    line = text.splitlines()[0].strip() if text.strip() else ""
    return line or None


def ac_by_absence(rows, mains, batts):
    """Rule 3: no mains-type supply, no battery, and no supply of unknown type."""
    return not mains and not batts and all("type" in a for _, a in rows)


def ac_by_batteries(batts):
    """Rule 4: at least one battery, and every one Charging, Full or Not charging."""
    return bool(batts) and all(
        a.get("status") in ("Charging", "Full", "Not charging") for a in batts
    )


def power(root):
    top = root / "sys/class/power_supply"
    try:
        entries = sorted(p for p in top.iterdir() if p.is_dir())
    except FileNotFoundError:
        entries = []
    except OSError as exc:
        return "unknown", None, f"cannot list /sys/class/power_supply: {exc.strerror}"
    rows = []
    for d in entries:
        a = {k: v for k in VALUES if (v := read_value(d / k)) is not None}
        if a.get("scope") == "Device":
            continue
        rows.append((d.name, a))
    mains = [a for _, a in rows if a.get("type") in MAINS]
    batts = [a for _, a in rows if a.get("type") == "Battery"]
    if any(a.get("online") in ("1", "2") for a in mains):
        state = "ac"
    elif any(a.get("status") == "Discharging" for a in batts):
        state = "battery"
    elif ac_by_absence(rows, mains, batts) or ac_by_batteries(batts):
        state = "ac"
    else:
        state = "unknown"
    cap = batts[0].get("capacity", "") if batts else ""
    pct = int(cap) if cap.isdigit() else None
    raw = "\n".join(
        f"{n}: " + (", ".join(f"{k}={v}" for k, v in sorted(a.items())) or "nothing readable")
        for n, a in rows
    )
    return state, pct, raw or NONE


def tree(root, supplies, make_dir=True):
    ps = root / "sys/class/power_supply"
    if make_dir:
        ps.mkdir(parents=True)
    for name, attrs in supplies.items():
        (ps / name).mkdir()
        for k, v in attrs.items():
            if v is IS_DIR:
                (ps / name / k).mkdir()
            elif isinstance(v, bytes):
                (ps / name / k).write_bytes(v)
            else:
                (ps / name / k).write_text(v + "\n")


IS_DIR = object()
UEVENT = (
    "POWER_SUPPLY_NAME=BAT0\nPOWER_SUPPLY_TYPE=Battery\n"
    "POWER_SUPPLY_STATUS=Discharging\nPOWER_SUPPLY_CAPACITY=64"
)
ROWS = [
    ("no /sys/class/power_supply at all", {}, False),
    ("an empty /sys/class/power_supply", {}, True),
    (
        "AC online 1; BAT0 Charging 81; mouse scope Device Discharging 40",
        {
            "AC": {"type": "Mains", "online": "1"},
            "BAT0": {"type": "Battery", "status": "Charging", "capacity": "81"},
            "hidpp_battery_0": {
                "type": "Battery",
                "scope": "Device",
                "status": "Discharging",
                "capacity": "40",
            },
        },
        True,
    ),
    (
        "AC online 0; BAT0 Discharging 64",
        {
            "AC": {"type": "Mains", "online": "0"},
            "BAT0": {"type": "Battery", "status": "Discharging", "capacity": "64"},
        },
        True,
    ),
    (
        "BAT0 Full 100, no mains entry",
        {"BAT0": {"type": "Battery", "status": "Full", "capacity": "100"}},
        True,
    ),
    (
        "AC online 0; BAT0 Unknown 50",
        {
            "AC": {"type": "Mains", "online": "0"},
            "BAT0": {"type": "Battery", "status": "Unknown", "capacity": "50"},
        },
        True,
    ),
    ("ups0 UPS Discharging only", {"ups0": {"type": "UPS", "status": "Discharging"}}, True),
    (
        "BAT0 files Full 100, its multi-line uevent saying Discharging 64",
        {
            "BAT0": {
                "type": "Battery",
                "status": "Full",
                "capacity": "100",
                "uevent": UEVENT,
                "manufacturer": "SMP",
                "voltage_now": "12541000",
            }
        },
        True,
    ),
    (
        "AC online unreadable (a directory); BAT0 status not UTF-8, capacity unreadable",
        {
            "AC": {"type": "Mains", "online": IS_DIR},
            "BAT0": {"type": "Battery", "status": b"\xff\xfe\n", "capacity": IS_DIR},
        },
        True,
    ),
    ("ACAD type unreadable, online 1; no battery", {"ACAD": {"type": IS_DIR, "online": "1"}}, True),
]
print("== Linux power, round 1 reader (five values; unreadable = missing)")
for label, supplies, make_dir in ROWS:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        tree(root, supplies, make_dir)
        state, pct, raw = power(root)
        print(f"  {label}: {state} {pct}")
        print("    raw: " + raw.replace("\n", " | "))

with tempfile.TemporaryDirectory() as d:
    root = Path(d)
    (root / "sys/class").mkdir(parents=True)
    (root / "sys/class/power_supply").write_text("not a directory\n")
    print("  /sys/class/power_supply a file, not a directory:", power(root))
