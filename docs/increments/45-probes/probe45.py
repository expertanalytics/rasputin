"""Increment 45 design probe: sketches of the cache rule and the Linux readers,
run against fixture trees, to fix the figures the tests pin. Not product code."""

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, sys.argv[1] + "/tools")
import bench


def cache_root(option, environ):
    if option is not None:
        return option
    if environ.get("RASPUTIN_DATA"):
        return Path(environ["RASPUTIN_DATA"]) / "cache"
    xdg = environ.get("XDG_CACHE_HOME", "")
    if xdg and Path(xdg).is_absolute():
        return Path(xdg) / "rasputin"
    home = environ.get("HOME", "")
    return Path(home) / ".cache" / "rasputin" if home and Path(home).is_absolute() else None


print("== cache_root")
for opt, env in [
    (Path("/given"), {"RASPUTIN_DATA": "/data", "XDG_CACHE_HOME": "/xdg", "HOME": "/home/u"}),
    (None, {"RASPUTIN_DATA": "/data", "XDG_CACHE_HOME": "/xdg", "HOME": "/home/u"}),
    (None, {"RASPUTIN_DATA": "", "XDG_CACHE_HOME": "/xdg", "HOME": "/home/u"}),
    (None, {"XDG_CACHE_HOME": "/xdg", "HOME": "/home/u"}),
    (None, {"XDG_CACHE_HOME": "rel/cache", "HOME": "/home/u"}),
    (None, {"XDG_CACHE_HOME": "", "HOME": "/home/u"}),
    (None, {"HOME": "/home/u"}),
    (None, {"HOME": ""}),
    (None, {}),
]:
    print(f"  {opt!s:8s} {json.dumps(env):70s} -> {cache_root(opt, env)}")

CPUINFO_X86 = "".join(
    f"processor\t: {p}\nvendor_id\t: GenuineIntel\n"
    "model name\t: Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz\n"
    f"physical id\t: {p // 4}\ncore id\t\t: {(p // 2) % 2}\n\n"
    for p in range(8)
)
CPUINFO_ARM = "".join(
    f"processor\t: {p}\nBogoMIPS\t: 243.75\nCPU implementer\t: 0x41\nCPU part\t: 0xd0c\n\n"
    for p in range(4)
)
MEMINFO = "MemTotal:       16374624 kB\nMemFree:         1000000 kB\n"
OS_RELEASE = 'NAME="Ubuntu"\nPRETTY_NAME="Ubuntu 24.04.3 LTS"\nID=ubuntu\n'


def tree(root, cpuinfo, topo, supplies):
    (root / "proc").mkdir(parents=True)
    (root / "proc/cpuinfo").write_text(cpuinfo)
    (root / "proc/meminfo").write_text(MEMINFO)
    (root / "etc").mkdir()
    (root / "etc/os-release").write_text(OS_RELEASE)
    for cpu, (pkg, core) in enumerate(topo):
        t = root / f"sys/devices/system/cpu/cpu{cpu}/topology"
        t.mkdir(parents=True)
        (t / "physical_package_id").write_text(f"{pkg}\n")
        (t / "core_id").write_text(f"{core}\n")
    ps = root / "sys/class/power_supply"
    ps.mkdir(parents=True)
    for name, attrs in supplies.items():
        (ps / name).mkdir()
        for k, v in attrs.items():
            (ps / name / k).write_text(v + "\n")


def field(text, key):
    for line in text.splitlines():
        k, _, v = line.partition(":")
        if k.strip() == key:
            return v.strip()
    return None


def machine(root):
    info = (root / "proc/cpuinfo").read_text()
    brand = field(info, "model name") or field(info, "Hardware") or "unknown"
    pairs = {
        (
            (d / "topology/physical_package_id").read_text().strip(),
            (d / "topology/core_id").read_text().strip(),
        )
        for d in (root / "sys/devices/system/cpu").glob("cpu[0-9]*")
        if (d / "topology/core_id").exists()
    }
    cores = len(pairs) or info.count("processor\t:")
    mem = int(field((root / "proc/meminfo").read_text(), "MemTotal").split()[0]) * 1024
    pretty = (
        field((root / "etc/os-release").read_text().replace("=", ":"), "PRETTY_NAME") or "Linux"
    ).strip('"')
    return brand, cores, mem, pretty


def power(root):
    rows = []
    for d in sorted((root / "sys/class/power_supply").iterdir()):
        a = {f.name: f.read_text().strip() for f in d.iterdir() if f.is_file()}
        if a.get("scope", "System") == "Device":
            continue
        rows.append((d.name, a))
    mains = [a for _, a in rows if a.get("type") in ("Mains", "USB", "Wireless")]
    batts = [a for _, a in rows if a.get("type") == "Battery"]
    if any(a.get("online", "0") != "0" for a in mains):
        state = "ac"
    elif any(a.get("status") == "Discharging" for a in batts):
        state = "battery"
    elif (not mains and not batts) or (
        batts and all(a.get("status") in ("Charging", "Full", "Not charging") for a in batts)
    ):
        state = "ac"
    else:
        state = "unknown"
    pct = next((int(a["capacity"]) for a in batts if a.get("capacity", "").isdigit()), None)
    raw = "".join(
        f"{n}: " + ", ".join(f"{k}={v}" for k, v in sorted(a.items())) + "\n" for n, a in rows
    )
    return state, pct, raw or "no power supply under /sys/class/power_supply\n"


X86_TOPO = [(p // 4, (p // 2) % 2) for p in range(8)]
LAPTOP_AC = {
    "AC": {"type": "Mains", "online": "1"},
    "BAT0": {"type": "Battery", "status": "Charging", "capacity": "81"},
    "hidpp_battery_0": {
        "type": "Battery",
        "scope": "Device",
        "status": "Discharging",
        "capacity": "40",
    },
}
LAPTOP_BATT = {
    "AC": {"type": "Mains", "online": "0"},
    "BAT0": {"type": "Battery", "status": "Discharging", "capacity": "64"},
}
NO_ADP_FULL = {"BAT0": {"type": "Battery", "status": "Full", "capacity": "100"}}
BATT_UNKNOWN = {
    "AC": {"type": "Mains", "online": "0"},
    "BAT0": {"type": "Battery", "status": "Unknown", "capacity": "50"},
}
UPS_ONLY = {"ups0": {"type": "UPS", "status": "Discharging"}}

print("== Linux machine and power")
for label, cpuinfo, topo, supplies in [
    ("x86 server, 2 sockets x 2 cores x 2 threads, no supply", CPUINFO_X86, X86_TOPO, {}),
    ("ARM, no topology dir, laptop on AC + mouse", CPUINFO_ARM, [], LAPTOP_AC),
    ("laptop on battery", CPUINFO_X86, X86_TOPO, LAPTOP_BATT),
    ("battery Full, no mains entry", CPUINFO_X86, X86_TOPO, NO_ADP_FULL),
    ("mains offline, battery status Unknown", CPUINFO_X86, X86_TOPO, BATT_UNKNOWN),
    ("UPS only", CPUINFO_X86, X86_TOPO, UPS_ONLY),
]:
    with tempfile.TemporaryDirectory() as d:
        root = Path(d)
        tree(root, cpuinfo, topo, supplies)
        print(f"  {label}: machine={machine(root)} power={power(root)[:2]}")
        if label.startswith("ARM"):
            print("    raw:", power(root)[2].replace("\n", " | "))

print("== the stored quick baseline under a Machine with system/os_version (alias macos)")
from pydantic import AliasChoices, BaseModel, Field  # noqa: E402


class Machine(BaseModel):
    cpu_brand: str
    p_cores: int
    e_cores: int
    memory_bytes: int
    system: str = "macOS"
    os_version: str = Field(validation_alias=AliasChoices("os_version", "macos"))
    python: str
    numpy: str
    shapely: str | None = None
    geos: str | None = None
    pyproj: str | None = None
    proj: str | None = None


stored = json.loads((Path(sys.argv[1]) / "docs/benchmarks/quick/baseline-ac.json").read_text())
old = Machine.model_validate(stored["machine"])
print("  loads:", old.system, old.os_version, "| dump keys:", list(old.model_dump()))
runs = sorted((Path(sys.argv[1]) / "docs/benchmarks").rglob("run.json"))
ok = sum(1 for p in runs if Machine.model_validate(json.loads(p.read_text())["machine"]))
print(f"  stored bench run.json files: {len(runs)}, machine block loads in {ok}")
print("  today's Machine fields:", list(bench.Machine.model_fields))
