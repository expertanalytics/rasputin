"""`tools/bench.py`'s hosts: the machine and the power state on macOS and Linux.

`docs/increments/45-linux-defaults.md`, sections 4.3 and 6 (T4, T5, T6, T8).
`LinuxHost(root)` reads text files under `root` only, so every Linux case here
is a fixture tree of `/proc`, `/sys` and `/etc` text under `tmp_path`, written
from the kernel's documented formats (`Documentation/ABI/testing/sysfs-class-power`,
`Documentation/admin-guide/cputopology.rst`, `proc(5)`, `os-release(5)`), not
captured from a machine. The one read of a real kernel is the last test, which
runs on Linux only.

The fixture conventions: a sysfs value is its text and a newline, as the
kernel writes it; an *unreadable* value is a directory in the file's place,
which raises `IsADirectoryError` whoever reads it (root included, so a root
container and CI's runner user agree); a value that is *not UTF-8* is the
bytes `ff fe 0a`.

`tools/` is not a package, so `bench.py` is loaded from its path in a fixture,
as `test_bench.py` loads it: while `LinuxHost` is missing each test fails on
its own.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from collections.abc import Iterator, Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

TOOL = Path(__file__).resolve().parents[2] / "tools" / "bench.py"


@pytest.fixture(scope="module")
def bench() -> Iterator[ModuleType]:
    spec = importlib.util.spec_from_file_location("bench", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["bench"] = module  # pydantic resolves annotations through it
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop("bench", None)


class NoProcess:
    """A runner that fails the test on any call: `LinuxHost` starts no process."""

    def run(self, argv: Any, cwd: Any = None, timeout: Any = None) -> Any:
        raise AssertionError(f"a process was started: {argv}")


@pytest.fixture
def no_process(bench: ModuleType, monkeypatch: pytest.MonkeyPatch) -> None:
    """Every way `bench.py` starts a process fails the test instead."""

    def refuse(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError(f"a process was started: {args}")

    monkeypatch.setattr(bench, "make_runner", NoProcess)
    monkeypatch.setattr(subprocess, "Popen", refuse)
    monkeypatch.setattr(subprocess, "run", refuse)


# ------------------------------------------------------------- machine (T4)

XEON = "Intel(R) Xeon(R) Platinum 8370C CPU @ 2.80GHz"
MEMINFO = "MemTotal:       16374624 kB\nMemFree:         9081232 kB\nMemAvailable:   12763920 kB\n"
OS_RELEASE = (
    'PRETTY_NAME="Ubuntu 24.04.3 LTS"\nNAME="Ubuntu"\nVERSION_ID="24.04"\n'
    'VERSION="24.04.3 LTS (Noble Numbat)"\nID=ubuntu\n'
)
#: Two sockets of two cores of two threads: eight `cpuN`, four cores.
XEON_TOPOLOGY = [(0, 0), (0, 0), (0, 1), (0, 1), (1, 0), (1, 0), (1, 1), (1, 1)]


def x86_cpuinfo(count: int) -> str:
    return "".join(
        f"processor\t: {n}\nvendor_id\t: GenuineIntel\ncpu family\t: 6\nmodel\t\t: 106\n"
        f"model name\t: {XEON}\nstepping\t: 6\ncpu MHz\t\t: 2793.437\n"
        f"physical id\t: {n // 4}\ncore id\t\t: {n // 2 % 2}\nflags\t\t: fpu vme de pse\n\n"
        for n in range(count)
    )


def arm_cpuinfo(count: int, hardware: str | None = None) -> str:
    """aarch64's `/proc/cpuinfo`: no `model name`; `Hardware` only on some boards."""
    body = "".join(
        f"processor\t: {n}\nBogoMIPS\t: 50.00\nFeatures\t: fp asimd evtstrm aes\n"
        "CPU implementer\t: 0x41\nCPU architecture: 8\nCPU variant\t: 0x3\n"
        "CPU part\t: 0xd0c\nCPU revision\t: 1\n\n"
        for n in range(count)
    )
    return body + (f"Hardware\t: {hardware}\n" if hardware else "")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def linux_root(
    root: Path,
    cpuinfo: str,
    topology: list[tuple[int, int]] | None,
    os_release: str | None = OS_RELEASE,
) -> Path:
    """`/proc/cpuinfo`, `/proc/meminfo`, the CPU topology (None: no `cpuN`
    directories at all) and `/etc/os-release` (None: absent)."""
    write_text(root / "proc" / "cpuinfo", cpuinfo)
    write_text(root / "proc" / "meminfo", MEMINFO)
    if topology is not None:
        cpus = root / "sys" / "devices" / "system" / "cpu"
        for n, (package, core) in enumerate(topology):
            write_text(cpus / f"cpu{n}" / "topology" / "physical_package_id", f"{package}\n")
            write_text(cpus / f"cpu{n}" / "topology" / "core_id", f"{core}\n")
        # Beside the cpuN directories, as on a real kernel: not CPUs.
        (cpus / "cpufreq").mkdir()
        (cpus / "cpuidle").mkdir()
        write_text(cpus / "online", f"0-{len(topology) - 1}\n")
    if os_release is not None:
        write_text(root / "etc" / "os-release", os_release)
    return root


def test_an_x86_server_counts_cores_not_threads(
    bench: ModuleType, tmp_path: Path, no_process: None
) -> None:
    root = linux_root(tmp_path, x86_cpuinfo(8), XEON_TOPOLOGY)
    machine = bench.LinuxHost(root).machine()
    assert machine.cpu_brand == XEON
    assert (machine.p_cores, machine.e_cores) == (4, 0)
    assert machine.memory_bytes == 16_767_614_976  # 16374624 kB times 1024
    assert machine.system == "Linux"
    assert machine.os_version == "Ubuntu 24.04.3 LTS"


def test_an_arm_board_without_topology_counts_processor_lines(
    bench: ModuleType, tmp_path: Path, no_process: None
) -> None:
    root = linux_root(tmp_path, arm_cpuinfo(4), None)
    machine = bench.LinuxHost(root).machine()
    assert machine.cpu_brand == "unknown"
    assert (machine.p_cores, machine.e_cores) == (4, 0)
    assert machine.memory_bytes == 16_767_614_976


def test_an_arm_board_naming_its_hardware_is_named_by_it(
    bench: ModuleType, tmp_path: Path, no_process: None
) -> None:
    root = linux_root(tmp_path, arm_cpuinfo(4, hardware="BCM2835"), None)
    assert bench.LinuxHost(root).machine().cpu_brand == "BCM2835"


def test_no_os_release_is_plain_linux(bench: ModuleType, tmp_path: Path, no_process: None) -> None:
    root = linux_root(tmp_path, x86_cpuinfo(8), XEON_TOPOLOGY, os_release=None)
    machine = bench.LinuxHost(root).machine()
    assert machine.os_version == "Linux"
    assert machine.system == "Linux"


# --------------------------------------------------------------- power (T5)

IS_DIR = object()  # a directory in the value's place: unreadable
NOT_UTF8 = b"\xff\xfe\n"
NO_SUPPLY = "no power supply under /sys/class/power_supply"
UEVENT = (
    "POWER_SUPPLY_NAME=BAT0\nPOWER_SUPPLY_TYPE=Battery\n"
    "POWER_SUPPLY_STATUS=Discharging\nPOWER_SUPPLY_CAPACITY=64"
)


def supplies(root: Path, entries: Mapping[str, Mapping[str, object]]) -> Path:
    """`<root>/sys/class/power_supply/<name>/<value>` for each entry."""
    top = root / "sys" / "class" / "power_supply"
    top.mkdir(parents=True)
    for name, values in entries.items():
        (top / name).mkdir()
        for key, value in values.items():
            path = top / name / key
            if value is IS_DIR:
                path.mkdir()
            elif isinstance(value, bytes):
                path.write_bytes(value)
            else:
                path.write_text(f"{value}\n")
    return root


AC_ON = {"type": "Mains", "online": "1"}
AC_OFF = {"type": "Mains", "online": "0"}


def battery(status: str, capacity: str) -> dict[str, object]:
    return {"type": "Battery", "status": status, "capacity": capacity}


#: The ten rows with supplies listed (the eleventh case, an unlistable
#: directory, is its own test): name, supplies (None: no directory at all),
#: state, percent, raw lines.
ROWS = [
    ("no-power-supply-directory", None, "ac", None, [NO_SUPPLY]),
    ("an-empty-power-supply-directory", {}, "ac", None, [NO_SUPPLY]),
    (
        "mains-online-with-a-device-scope-mouse-battery-dropped",
        {
            "AC": AC_ON,
            "BAT0": battery("Charging", "81"),
            "hidpp_battery_0": {**battery("Discharging", "40"), "scope": "Device"},
        },
        "ac", 81,
        ["AC: online=1, type=Mains", "BAT0: capacity=81, status=Charging, type=Battery"],
    ),
    (
        "mains-offline-battery-discharging",
        {"AC": AC_OFF, "BAT0": battery("Discharging", "64")},
        "battery", 64,
        ["AC: online=0, type=Mains", "BAT0: capacity=64, status=Discharging, type=Battery"],
    ),
    (
        "battery-full-and-no-mains-entry",
        {"BAT0": battery("Full", "100")},
        "ac", 100,
        ["BAT0: capacity=100, status=Full, type=Battery"],
    ),
    (
        "mains-offline-battery-status-unknown",
        {"AC": AC_OFF, "BAT0": battery("Unknown", "50")},
        "unknown", 50,
        ["AC: online=0, type=Mains", "BAT0: capacity=50, status=Unknown, type=Battery"],
    ),
    (
        "a-ups-alone-decides-nothing",
        {"ups0": {"type": "UPS", "status": "Discharging"}},
        "ac", None,
        ["ups0: status=Discharging, type=UPS"],
    ),
    (
        "uevent-and-vendor-files-are-never-read",
        {"BAT0": {**battery("Full", "100"), "uevent": UEVENT, "manufacturer": "SMP",
                  "voltage_now": "12541000"}},
        "ac", 100,
        ["BAT0: capacity=100, status=Full, type=Battery"],
    ),
    (
        "unreadable-and-non-utf8-values-are-missing",
        {"AC": {"type": "Mains", "online": IS_DIR},
         "BAT0": {"type": "Battery", "status": NOT_UTF8, "capacity": IS_DIR}},
        "unknown", None,
        ["AC: type=Mains", "BAT0: type=Battery"],
    ),
    (
        "a-supply-of-unreadable-type-may-be-a-battery",
        {"ACAD": {"type": IS_DIR, "online": "1"}},
        "unknown", None,
        ["ACAD: online=1"],
    ),
]  # fmt: skip


@pytest.mark.parametrize(
    ("entries", "state", "percent", "raw"),
    [pytest.param(*row[1:], id=row[0]) for row in ROWS],
)
def test_linux_power(
    bench: ModuleType,
    tmp_path: Path,
    no_process: None,
    entries: Mapping[str, Mapping[str, object]] | None,
    state: str,
    percent: int | None,
    raw: list[str],
) -> None:
    """Rules 1 to 5 in order: mains online is AC; else a discharging battery
    is battery; else no supply of any (readable) type is AC; else batteries
    all charging, full or not charging is AC; else unknown. The `uevent` row
    fails a reader of `uevent`, which says Discharging 64."""
    root = tmp_path if entries is None else supplies(tmp_path, entries)
    power = bench.LinuxHost(root).power()
    assert (power.state, power.percent) == (state, percent)
    assert power.raw.splitlines() == raw


def test_an_unlistable_power_supply_directory_is_unknown_with_the_reason(
    bench: ModuleType, tmp_path: Path, no_process: None
) -> None:
    """A file in the directory's place: listing raises `NotADirectoryError`,
    root or not."""
    write_text(tmp_path / "sys" / "class" / "power_supply", "not a directory\n")
    power = bench.LinuxHost(tmp_path).power()
    assert (power.state, power.percent) == ("unknown", None)
    assert power.raw == "cannot list /sys/class/power_supply: Not a directory"


def test_a_supply_is_found_through_its_link_as_on_a_real_kernel(
    bench: ModuleType, tmp_path: Path, no_process: None
) -> None:
    """On a kernel each entry of `/sys/class/power_supply` is a symbolic link
    into `/sys/devices`; the link is followed."""
    devices = tmp_path / "sys" / "devices" / "LNXSYSTM:00" / "power_supply"
    for name, values in (("AC", AC_OFF), ("BAT0", battery("Discharging", "64"))):
        for key, value in values.items():
            write_text(devices / name / key, f"{value}\n")
    top = tmp_path / "sys" / "class" / "power_supply"
    top.mkdir(parents=True)
    for name in ("AC", "BAT0"):
        (top / name).symlink_to(devices / name, target_is_directory=True)
    power = bench.LinuxHost(tmp_path).power()
    assert (power.state, power.percent) == ("battery", 64)


# ------------------------------------------------------------ factory (T6)


def test_darwin_gets_the_mac_host(bench: ModuleType) -> None:
    assert isinstance(bench.make_host(NoProcess(), "darwin"), bench.MacHost)


def test_linux_gets_the_linux_host_rooted_at_slash(bench: ModuleType) -> None:
    host = bench.make_host(NoProcess(), "linux")
    assert isinstance(host, bench.LinuxHost)
    assert Path("/") in vars(host).values(), vars(host)


def test_another_system_is_refused_naming_it(bench: ModuleType) -> None:
    with pytest.raises(bench.BenchError, match="win32"):
        bench.make_host(NoProcess(), "win32")


def test_the_mac_host_runs_what_the_mac_reader_ran(bench: ModuleType) -> None:
    """`sysctl`, `sw_vers` and `pmset -g batt`, read as before."""
    canned = {
        "machdep.cpu.brand_string": "Apple M1 Max", "hw.perflevel0.physicalcpu": "8",
        "hw.perflevel1.physicalcpu": "2", "hw.memsize": str(32 * 2**30),
    }  # fmt: skip
    pmset = "Now drawing from 'AC Power'\n -InternalBattery-0 (id=1)\t97%; charged; present: true\n"

    class Canned:
        def run(self, argv: Any, cwd: Any = None, timeout: Any = None) -> Any:
            name = Path(argv[0]).name
            out = {"sysctl": canned.get(argv[-1], "0") + "\n", "sw_vers": "27.0\n",
                   "pmset": pmset}[name]  # fmt: skip
            return bench.Completed(returncode=0, stdout=out, stderr="", wall_s=0.01)

    host = bench.MacHost(Canned())
    machine = host.machine()
    assert (machine.cpu_brand, machine.p_cores, machine.e_cores) == ("Apple M1 Max", 8, 2)
    assert (machine.system, machine.os_version) == ("macOS", "27.0")
    assert host.power() == bench.parse_pmset(pmset)


# ------------------------------------------------------- the real kernel (T8)


@pytest.mark.skipif(sys.platform != "linux", reason="reads this machine's /proc and /sys")
def test_the_real_kernel(bench: ModuleType) -> None:
    host = bench.LinuxHost()
    machine = host.machine()
    assert machine.p_cores >= 1
    assert machine.memory_bytes > 0
    assert machine.system == "Linux"
    assert host.power().state in ("ac", "battery", "unknown")
