"""`tools/bench_quick.py`, the 15-minute quick check (`docs/increments/perf-quick-check.md`).

The design's section 8 names the seams tested here: `verdict(new, base)` and
`hotspots(record, ruled)`, pure over records built in the test; the directory
fingerprint over a `tmp_path` tree; the deadline, with a fake runner whose
child times out; and one real child (skipped without `_core`) proving that the
`PhaseClock` replacement in `bench.py _child` sees `mesh`'s clock.

Names, fields and strings the design leaves open are pinned here and listed in
the red step's handback ("Pinned or assumed beyond the design"):

* a record is `QuickRecord`, built with `model_validate` from the dict shape
  `record()` below returns; a stored baseline has the same shape;
* the measure `total` is the clock's total (`total_s`); every other measure is
  a `--stats` phase row by its name;
* a case run at more than one thread count is recorded once per count, the
  default (0) as `<name>`, another as `<name> t=<n>`;
* `fingerprint(path)` returns a `str`, compared only for equality;
* `cases.toml` is `[[case]]` tables of `name`, `threads`, `runs`, `warmup`,
  `args` (the `mesh` argv, without `--out` and `--binary`, which the tool adds)
  and `inputs`; `hotspots.toml` is `[[hotspot]]` tables of `case`, `phase`,
  `share` (percent), `date` and `ruling`;
* `QUICK` is the directory holding `cases.toml`, `baseline-<power>.json` and
  `hotspots.toml`, and `make_runner` the factory the tests replace, both
  module attributes of `bench_quick`.

`tools/` is not a package, so the module is loaded from its path in a fixture,
as `test_bench.py` loads `bench.py`: while `bench_quick.py` is missing each
test fails on its own instead of one collection error stopping the suite.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import os
import subprocess
import sys
from collections.abc import Callable, Iterator, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any, NamedTuple

import numpy as np
import pytest
from typer.testing import CliRunner

from geotiff_fixtures import micro_tiff

REPO = Path(__file__).resolve().parents[2]
TOOLS = REPO / "tools"
TOOL = TOOLS / "bench_quick.py"
BENCH = TOOLS / "bench.py"
SHIPPED = REPO / "docs" / "benchmarks" / "quick"

cli = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})


@pytest.fixture(scope="module")
def quick() -> Iterator[ModuleType]:
    """``bench_quick`` loaded from its path, with ``tools/`` importable for the
    ``bench`` module it builds on; both are removed from ``sys.modules`` after."""
    sys.path.insert(0, str(TOOLS))
    spec = importlib.util.spec_from_file_location("bench_quick", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["bench_quick"] = module  # pydantic resolves annotations through it
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop("bench_quick", None)
        sys.modules.pop("bench", None)
        sys.path.remove(str(TOOLS))


# ------------------------------------------------------------------ records

MACHINE = {
    "cpu_brand": "Apple M4 Pro",
    "p_cores": 10,
    "e_cores": 4,
    "memory_bytes": 48 * 2**30,
    "macos": "26.0",
    "python": "3.14.0",
    "numpy": "2.3.0",
}


def spread(median: float, low: float | None = None, high: float | None = None) -> dict[str, float]:
    return {"median": median, "min": median if low is None else low,
            "max": median if high is None else high}  # fmt: skip


def lagan(**measures: dict[str, float]) -> dict[str, Any]:
    """Lagan after 30d (section 3): 16.97 s, `features clip` 40.7 % of it."""
    base = {"total": spread(16.97, 16.90, 17.10), "features clip": spread(6.89, 6.80, 6.95),
            "decode": spread(2.00, 1.98, 2.02)}  # fmt: skip
    return {
        "tolerance": 10.0,
        "max_error": 9.97,
        "mesh_sha256": "c" * 64,
        "measures": {**base, **measures},
    }


def record(**cases: dict[str, Any]) -> dict[str, Any]:
    return {
        "commit": "a" * 40,
        "machine": copy.deepcopy(MACHINE),
        "power": "ac",
        "inputs": {"lagan dem": "f" * 64, "lagan features": "e" * 64},
        "cases": cases or {"lagan": lagan()},
    }


def make(quick: ModuleType, data: dict[str, Any]) -> Any:
    return quick.QuickRecord.model_validate(data)


def judge(quick: ModuleType, new: dict[str, Any], base: dict[str, Any] | None) -> Any:
    return quick.verdict(make(quick, new), None if base is None else make(quick, base))


def findings(v: Any, kind: str) -> list[str]:
    return [line for line in v.lines if line.startswith(f"{kind}: ")]


# ------------------------------------------------------------------- verdict


def test_identical_runs_are_no_change(quick: ModuleType) -> None:
    v = judge(quick, record(), record())
    assert (v.status, v.exit_code) == ("NO CHANGE", 0)
    assert findings(v, "SLOWER") == findings(v, "FASTER") == []


def test_a_branch_commit_against_the_master_baseline_is_compared(quick: ModuleType) -> None:
    """The baseline is a master run and the check runs on a branch head."""
    new = record()
    new["commit"] = "b" * 40
    assert judge(quick, new, record()).status == "NO CHANGE"


def test_a_slower_phase_names_case_phase_times_size_and_band(quick: ModuleType) -> None:
    """Section 3's example line, to the character."""
    v = judge(quick, record(lagan=lagan(**{"features clip": spread(8.10)})), record())
    assert (v.status, v.exit_code) == ("SLOWER", 1)
    assert findings(v, "SLOWER") == [
        "SLOWER: lagan features clip 6.89 -> 8.10 s (+17.6 %, band 5.0 %)"
    ]


def test_a_slower_total_is_named_total(quick: ModuleType) -> None:
    v = judge(quick, record(lagan=lagan(total=spread(18.50))), record())
    assert v.status == "SLOWER"
    assert findings(v, "SLOWER") == ["SLOWER: lagan total 16.97 -> 18.50 s (+9.0 %, band 5.0 %)"]


# *Scale* for the band tests: a 10 s total, as in section 3's runs of 1 to
# 20 s; band edges are approached no closer than 1 % (0.1 s), checked at 10 s.


@pytest.mark.parametrize(
    ("new_median", "status"),
    [pytest.param(10.40, "NO CHANGE", id="+4%"), pytest.param(10.60, "SLOWER", id="+6%")],
)
def test_the_band_is_five_percent_when_the_baseline_is_tight(
    quick: ModuleType, new_median: float, status: str
) -> None:
    base = record(lagan=lagan(total=spread(10.0, 9.9, 10.1)))  # spread 2 %
    v = judge(quick, record(lagan=lagan(total=spread(new_median))), base)
    assert v.status == status


@pytest.mark.parametrize(
    ("new_median", "status"),
    [pytest.param(11.50, "NO CHANGE", id="+15%"), pytest.param(12.50, "SLOWER", id="+25%")],
)
def test_a_wide_baseline_spread_widens_the_band(
    quick: ModuleType, new_median: float, status: str
) -> None:
    base = record(lagan=lagan(total=spread(10.0, 9.0, 11.0)))  # (11 - 9) / 10 = 20 %
    v = judge(quick, record(lagan=lagan(total=spread(new_median))), base)
    assert v.status == status
    if status == "SLOWER":
        assert findings(v, "SLOWER") == [
            "SLOWER: lagan total 10.00 -> 12.50 s (+25.0 %, band 20.0 %)"
        ]


def test_the_new_runs_spread_does_not_widen_the_band(quick: ModuleType) -> None:
    """Section 3: one slow outlier in the check cannot hide a slowdown."""
    base = record(lagan=lagan(total=spread(10.0, 9.9, 10.1)))
    v = judge(quick, record(lagan=lagan(total=spread(10.60, 8.0, 14.0))), base)
    assert v.status == "SLOWER"
    assert "band 5.0 %" in findings(v, "SLOWER")[0]


@pytest.mark.parametrize(
    ("new_median", "status"),
    [
        pytest.param(0.54, "NO CHANGE", id="+8%-but-0.04s"),
        pytest.param(0.56, "SLOWER", id="+12%-and-0.06s"),
    ],
)
def test_a_slowdown_must_also_exceed_five_hundredths_of_a_second(
    quick: ModuleType, new_median: float, status: str
) -> None:
    """*Scale:* a 0.5 s phase of a 5 s total, the smallest phase the design
    judges (0.1 s and up); 0.01 s from the 0.05 s floor either side."""
    base = record(lagan=lagan(total=spread(5.0), **{"features clip": spread(0.50)}))
    new = record(lagan=lagan(total=spread(5.0), **{"features clip": spread(new_median)}))
    assert judge(quick, new, base).status == status


def test_a_phase_under_five_percent_of_the_baseline_total_is_not_judged(
    quick: ModuleType,
) -> None:
    base = record(lagan=lagan(total=spread(10.0), trim=spread(0.45)))  # 4.5 %
    v = judge(quick, record(lagan=lagan(total=spread(10.0), trim=spread(2.0))), base)
    assert v.status == "NO CHANGE"
    assert not any("trim" in line for line in v.lines)


def test_a_phase_at_exactly_five_percent_of_the_baseline_total_is_judged(
    quick: ModuleType,
) -> None:
    """Section 3's "at 5 % or more": 0.5 of 10.0 is exactly 5 %."""
    base = record(lagan=lagan(total=spread(10.0), trim=spread(0.50)))
    v = judge(quick, record(lagan=lagan(total=spread(10.0), trim=spread(2.0))), base)
    assert v.status == "SLOWER"
    assert [line for line in findings(v, "SLOWER") if " trim " in line]


def test_faster_is_symmetric_and_exits_0(quick: ModuleType) -> None:
    v = judge(quick, record(lagan=lagan(total=spread(15.00))), record())
    assert (v.status, v.exit_code) == ("FASTER", 0)
    assert findings(v, "FASTER") == ["FASTER: lagan total 16.97 -> 15.00 s (-11.6 %, band 5.0 %)"]


def test_slower_outranks_faster(quick: ModuleType) -> None:
    new = record(lagan=lagan(total=spread(16.97), decode=spread(1.0),
                             **{"features clip": spread(8.10)}))  # fmt: skip
    v = judge(quick, new, record())
    assert (v.status, v.exit_code) == ("SLOWER", 1)
    assert findings(v, "FASTER") and findings(v, "SLOWER")


def test_max_error_over_the_tolerance_is_broken(quick: ModuleType) -> None:
    case = lagan()
    case["max_error"] = 10.2
    v = judge(quick, record(lagan=case), record())
    assert (v.status, v.exit_code) == ("BROKEN", 1)
    assert [line for line in findings(v, "BROKEN") if "lagan" in line and "max_error" in line]


def test_max_error_at_the_tolerance_is_not_broken(quick: ModuleType) -> None:
    case = lagan()
    case["max_error"] = 10.0
    assert judge(quick, record(lagan=case), record()).status == "NO CHANGE"


def test_a_changed_mesh_hash_is_reported_not_judged(quick: ModuleType) -> None:
    case = lagan()
    case["mesh_sha256"] = "d" * 64
    v = judge(quick, record(lagan=case), record())
    assert (v.status, v.exit_code) == ("NO CHANGE", 0)
    assert [line for line in v.lines if "lagan" in line and "mesh" in line]


def test_no_stored_baseline_is_no_baseline_exit_2(quick: ModuleType) -> None:
    v = judge(quick, record(), None)
    assert v.status.startswith("NO BASELINE")
    assert v.exit_code == 2


@pytest.mark.parametrize(
    ("change", "field"),
    [
        pytest.param({"power": "battery"}, "power", id="power"),
        pytest.param({"machine": {**MACHINE, "p_cores": 8}}, "machine", id="machine"),
        pytest.param({"inputs": {"lagan dem": "0" * 64, "lagan features": "e" * 64}},
                     "lagan dem", id="input"),
    ],
)  # fmt: skip
def test_a_mismatch_is_no_baseline_naming_the_field_never_a_comparison(
    quick: ModuleType, change: dict[str, Any], field: str
) -> None:
    new = {**record(lagan=lagan(total=spread(30.0))), **change}  # 77 % slower: never judged
    v = judge(quick, new, record())
    assert v.status.startswith("NO BASELINE: ")
    assert field in v.status
    assert v.exit_code == 2
    assert findings(v, "SLOWER") == []


@pytest.mark.parametrize("power", ["mixed", "unknown"])
def test_an_indeterminate_power_state_is_never_a_baseline(quick: ModuleType, power: str) -> None:
    new, base = record(), record()
    new["power"] = base["power"] = power
    v = judge(quick, new, base)
    assert v.status.startswith("NO BASELINE: ")
    assert v.exit_code == 2


# ------------------------------------------------------------------ hotspots


def ruled(quick: ModuleType, share: float, case: str = "lagan") -> Any:
    return quick.RuledHotspot(case=case, phase="features clip", share=share,
                              date="2026-10-08", ruling="Ola: known, leave it")  # fmt: skip


def hot(quick: ModuleType, clip: float, ruled_list: Sequence[Any] = ()) -> list[str]:
    """`features clip` at ``clip`` s of a 10 s total; `decode` at 2 s (20 %)."""
    rec = make(quick, record(lagan=lagan(total=spread(10.0), decode=spread(2.0),
                                         **{"features clip": spread(clip)})))  # fmt: skip
    lines: list[str] = quick.hotspots(rec, list(ruled_list))
    return lines


def test_a_phase_at_forty_percent_or_more_is_a_hotspot(quick: ModuleType) -> None:
    """Section 3's example line: 40.7 % of the run, unruled."""
    assert hot(quick, 4.07) == [
        "HOTSPOT: lagan features clip 41 % of the run (not in hotspots.toml)"
    ]


def test_exactly_forty_percent_is_a_hotspot(quick: ModuleType) -> None:
    assert len(hot(quick, 4.0)) == 1


def test_under_forty_percent_is_not_a_hotspot(quick: ModuleType) -> None:
    assert hot(quick, 3.99) == []


def test_the_total_is_never_a_hotspot(quick: ModuleType) -> None:
    assert not any(" total " in line for line in hot(quick, 4.07))


def test_a_ruled_hotspot_is_silent(quick: ModuleType) -> None:
    assert hot(quick, 4.07, [ruled(quick, 40.7)]) == []


def test_a_ruled_hotspot_grown_under_ten_points_is_silent(quick: ModuleType) -> None:
    assert hot(quick, 4.97, [ruled(quick, 40.7)]) == []  # +9 points


def test_a_ruled_hotspot_grown_more_than_ten_points_is_raised_again(quick: ModuleType) -> None:
    lines = hot(quick, 5.17, [ruled(quick, 40.7)])  # +11 points
    assert len(lines) == 1
    assert lines[0].startswith("HOTSPOT: lagan features clip 52 % of the run")


def test_a_ruling_on_another_case_does_not_silence(quick: ModuleType) -> None:
    assert len(hot(quick, 4.07, [ruled(quick, 40.7, case="numedalslagen")])) == 1


def test_a_hotspot_does_not_change_the_verdict(quick: ModuleType) -> None:
    v = judge(quick, record(), record())  # lagan's features clip is 40.6 % of 16.97 s
    assert (v.status, v.exit_code) == ("NO CHANGE", 0)


def test_load_hotspots_reads_the_toml_entries(quick: ModuleType, tmp_path: Path) -> None:
    path = tmp_path / "hotspots.toml"
    path.write_text(
        '[[hotspot]]\ncase = "lagan"\nphase = "features clip"\nshare = 40.7\n'
        'date = "2026-10-08"\nruling = "Ola: known, leave it"\n'
    )
    assert quick.load_hotspots(path) == [ruled(quick, 40.7)]


# -------------------------------------------------------- shipped data files


def test_the_shipped_cases_are_sections_3_table(quick: ModuleType) -> None:
    cases = quick.load_cases(SHIPPED / "cases.toml")
    table = [(c.name, list(c.threads), c.warmup, c.runs) for c in cases]
    assert table == [
        ("tile", [0, 1], 0, 5),
        ("quarter", [0], 0, 5),
        ("numedalslagen", [0], 1, 3),
        ("lagan", [0], 1, 3),
    ]
    tolerance = {c.name: c.args[c.args.index("--tolerance") + 1] for c in cases}
    assert {k: float(v) for k, v in tolerance.items()} == {
        "tile": 1.0, "quarter": 1.0, "numedalslagen": 10.0, "lagan": 10.0
    }  # fmt: skip
    lagan_args = next(c.args for c in cases if c.name == "lagan")
    assert lagan_args[lagan_args.index("--out-crs") + 1] == "EPSG:3006"
    for c in cases:
        assert "--out" not in c.args and "--binary" not in c.args, c.name


def test_the_shipped_hotspots_file_loads(quick: ModuleType) -> None:
    ruled_list = quick.load_hotspots(SHIPPED / "hotspots.toml")
    assert all(r.ruling for r in ruled_list)


# ------------------------------------------------------------- fingerprints


def tree(root: Path) -> Path:
    for name, body in (("a.tif", b"aaaa"), ("sub/b.tif", b"bbbbbb"), ("sub/deep/c.tif", b"c")):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(body)
        os.utime(path, ns=(1_700_000_000_000_000_000, 1_700_000_000_000_000_000))
    return root


def test_a_directory_fingerprint_is_stable(quick: ModuleType, tmp_path: Path) -> None:
    root = tree(tmp_path / "dtm10")
    assert quick.fingerprint(root) == quick.fingerprint(root)


def test_a_directory_fingerprint_uses_relative_names(quick: ModuleType, tmp_path: Path) -> None:
    """The same files under another root: the same fingerprint."""
    assert quick.fingerprint(tree(tmp_path / "one")) == quick.fingerprint(tree(tmp_path / "two"))


def rename(path: Path) -> None:
    path.rename(path.with_name("renamed.tif"))


def resize(path: Path) -> None:
    stat = path.stat()
    path.write_bytes(path.read_bytes() + b"x")
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))


def touch(path: Path) -> None:
    os.utime(path, ns=(1_800_000_000_000_000_000, 1_800_000_000_000_000_000))


def remove(path: Path) -> None:
    path.unlink()


@pytest.mark.parametrize("change", [rename, resize, touch, remove])
def test_a_directory_fingerprint_sees_a_change_deep_in_the_tree(
    quick: ModuleType, tmp_path: Path, change: Callable[[Path], None]
) -> None:
    root = tree(tmp_path / "dtm10")
    before = quick.fingerprint(root)
    change(root / "sub" / "deep" / "c.tif")
    assert quick.fingerprint(root) != before


def test_a_directory_fingerprint_reads_no_content(quick: ModuleType, tmp_path: Path) -> None:
    """Same names, sizes and mtimes, other bytes: the same fingerprint."""
    root = tree(tmp_path / "dtm10")
    before = quick.fingerprint(root)
    target = root / "sub" / "b.tif"
    stat = target.stat()
    target.write_bytes(b"BBBBBB")
    os.utime(target, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert quick.fingerprint(root) == before


def test_a_small_file_fingerprint_reads_its_content(quick: ModuleType, tmp_path: Path) -> None:
    path = tmp_path / "outline.geojson"
    path.write_bytes(b"{}  ")
    stat = path.stat()
    before = quick.fingerprint(path)
    path.write_bytes(b"[]  ")  # same size
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert quick.fingerprint(path) != before


def test_a_file_of_100_mb_or_more_is_not_read(quick: ModuleType, tmp_path: Path) -> None:
    """101 MiB, sparse: over the 100 MB line however MB is read. A byte changed
    with size and mtime kept: the fingerprint must not see it."""
    path = tmp_path / "corine.gpkg"
    with path.open("wb") as f:
        f.truncate(101 * 2**20)
    stat = path.stat()
    before = quick.fingerprint(path)
    with path.open("r+b") as f:
        f.seek(2**20)
        f.write(b"\x01")
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    assert quick.fingerprint(path) == before


def test_a_file_under_100_mb_that_is_touched_keeps_its_fingerprint(
    quick: ModuleType, tmp_path: Path
) -> None:
    """Git sets a file's mtime at checkout, so the same fixture has another
    mtime in each checkout: under 100 MB, size and content alone decide."""
    path = tmp_path / "outline.geojson"
    path.write_bytes(b"{}")
    before = quick.fingerprint(path)
    touch(path)
    assert quick.fingerprint(path) == before


def test_a_file_of_100_mb_or_more_that_is_touched_gets_a_new_fingerprint(
    quick: ModuleType, tmp_path: Path
) -> None:
    """101 MiB, sparse, as above: its content is not read, so mtime stays in."""
    path = tmp_path / "corine.gpkg"
    with path.open("wb") as f:
        f.truncate(101 * 2**20)
    before = quick.fingerprint(path)
    touch(path)
    assert quick.fingerprint(path) != before


# ------------------------------------------------------ the run, with a fake

AC = "Now drawing from 'AC Power'\n -InternalBattery-0 (id=1)\t100%; charged; present: true\n"
SYSCTL = {
    "machdep.cpu.brand_string": "Apple M4 Pro",
    "hw.perflevel0.physicalcpu": "10",
    "hw.perflevel1.physicalcpu": "4",
    "hw.memsize": str(48 * 2**30),
}
VTK = b"# vtk DataFile Version 3.0\nrasputin\nBINARY\nDATASET POLYDATA\nPOINTS 0 double\n"


class Done(NamedTuple):
    """What a ``Runner`` returns: ``bench.Completed`` with ``timed_out``."""

    returncode: int
    stdout: str
    stderr: str
    wall_s: float
    timed_out: bool = False


class FakeRunner:
    """Canned output keyed on argv[0]; the child is recognised by ``_child``.

    Each child reports the next of ``totals`` (cycled) as its clock's total,
    with `decode` a fifth of it; a child whose case name is in ``time_out``
    reports a timeout instead. Every call's ``timeout`` is kept.
    """

    def __init__(self, totals: Sequence[float] = (1.0,), time_out: Sequence[str] = ()) -> None:
        self.totals = list(totals)
        self.time_out = set(time_out)
        self.calls: list[tuple[list[str], float | None]] = []
        self.child_count = 0

    def run(self, argv: Sequence[str], cwd: Path | None = None,
            timeout: float | None = None) -> Done:  # fmt: skip
        args = [str(a) for a in argv]
        self.calls.append((args, timeout))
        program = Path(args[0]).name
        if "_child" in args:
            return self._child(args)
        if program == "pmset":
            return Done(0, AC, "", 0.01)
        if program == "sysctl":
            return Done(0, SYSCTL.get(args[-1], "0") + "\n", "", 0.01)
        if program == "sw_vers":
            return Done(0, "26.0\n", "", 0.01)
        if program == "git":
            out = "a" * 40 + "\n" if "rev-parse" in args else ""
            return Done(0, out, "", 0.01)
        if program == "cmake":
            return Done(0, "", "", 0.01)
        return Done(127, "", f"{program}: not faked", 0.01)

    def _child(self, args: list[str]) -> Done:
        mesh = args[args.index("--") + 1 :]
        if "--out" in mesh:
            Path(mesh[mesh.index("--out") + 1]).write_bytes(VTK)
        if any(marker in mesh for marker in self.time_out):
            return Done(-9, "", "", 290.0, timed_out=True)
        total = self.totals[self.child_count % len(self.totals)]
        self.child_count += 1
        phases = {"decode": 0.2 * total, "refine": 0.1 * total}
        values = {"refine_s": 0.1 * total, "app_s": total, "max_error": 0.5, "rounds": 3,
                  "inserted": 10, "flips": 20, "hardening": "none",
                  "phases": phases, "total_s": total}  # fmt: skip
        return Done(0, "", "BENCH " + json.dumps(values) + "\n", total)

    def children(self) -> list[list[str]]:
        return [args for args, _ in self.calls if "_child" in args]

    def started(self, case_input: str) -> bool:
        return any(case_input in args for args in self.children())


@pytest.fixture
def place(tmp_path: Path, quick: ModuleType, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    """A tree whose faked build leaves one ``_core``, a ``QUICK`` directory and
    two cases, `alpha` and `beta`, each with its own input file."""
    tree_dir = tmp_path / "tree"
    (tree_dir / "src_python" / "tin_engine").mkdir(parents=True)
    (tree_dir / "build-bench").mkdir()
    (tree_dir / "build-bench" / "_core.cpython-314-darwin.so").write_bytes(b"never loaded")
    data = tmp_path / "data"
    data.mkdir()
    places = {"tree": tree_dir, "quick": tmp_path / "quick", "data": data}
    places["quick"].mkdir()
    for name in ("alpha", "beta"):
        places[name] = data / f"{name}.tif"
        places[name].write_bytes(name.encode())
    monkeypatch.setattr(quick, "QUICK", places["quick"])
    monkeypatch.setenv("RASPUTIN_DATA", str(data))
    return places


def write_cases(places: dict[str, Path], alpha: str = "", beta: str = "") -> None:
    """Two cases on absolute paths, so no `$RASPUTIN_DATA` convention is pinned."""

    def case(name: str, extra: str) -> str:
        dem = places[name]
        return (f'[[case]]\nname = "{name}"\nthreads = [0]\nruns = 3\nwarmup = 0\n'
                f'args = ["mesh", "--dem", "{dem}", "--tolerance", "1"]\n'
                f'inputs = ["{dem}"]\n{extra}\n')  # fmt: skip

    (places["quick"] / "cases.toml").write_text(case("alpha", alpha) + case("beta", beta))


@pytest.fixture
def use_runner(quick: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Callable[[FakeRunner], None]:
    def install(fake: FakeRunner) -> None:
        monkeypatch.setattr(quick, "make_runner", lambda: fake)

    return install


def invoke(quick: ModuleType, places: dict[str, Path], *extra: str) -> Any:
    return cli.invoke(quick.app, ["run", "--tree", str(places["tree"]), *extra])


def last_line(output: str) -> str:
    return [line for line in output.splitlines() if line.strip()][-1]


def save_baseline(quick: ModuleType, places: dict[str, Path], use: Callable[..., None]) -> Path:
    use(FakeRunner())
    invoke(quick, places, "--save-baseline")
    path = places["quick"] / "baseline-ac.json"
    assert path.is_file(), "--save-baseline wrote no baseline-ac.json"
    return path


def test_a_budget_over_900_seconds_is_refused_before_anything_runs(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    write_cases(place)
    fake = FakeRunner()
    use_runner(fake)
    result = invoke(quick, place, "--budget", "901")
    assert result.exit_code == 3, result.output
    assert "--budget" in result.output
    assert fake.calls == []


def test_save_baseline_writes_median_min_max_without_the_warm_up(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    """One warm-up then three timed runs: the warm-up's 100 s is never seen."""
    write_cases(place, alpha="", beta="")
    text = (place["quick"] / "cases.toml").read_text().replace("warmup = 0", "warmup = 1", 1)
    (place["quick"] / "cases.toml").write_text(text)
    fake = FakeRunner(totals=(100.0, 1.0, 2.0, 4.0))
    use_runner(fake)
    invoke(quick, place, "--save-baseline")
    data = json.loads((place["quick"] / "baseline-ac.json").read_text())
    assert data["power"] == "ac"
    assert data["commit"] == "a" * 40
    assert data["cases"]["alpha"]["measures"]["total"] == {"median": 2.0, "min": 1.0, "max": 4.0}
    assert data["cases"]["alpha"]["measures"]["decode"]["median"] == pytest.approx(0.4)
    assert len(fake.children()) == 4 + 3
    assert all("--binary" in child for child in fake.children())


def test_a_case_at_two_thread_counts_is_recorded_once_per_count(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    write_cases(place)
    path = place["quick"] / "cases.toml"
    path.write_text(path.read_text().replace("threads = [0]", "threads = [0, 1]", 1))
    fake = FakeRunner()
    use_runner(fake)
    invoke(quick, place, "--save-baseline")
    data = json.loads((place["quick"] / "baseline-ac.json").read_text())
    assert sorted(data["cases"]) == ["alpha", "alpha t=1", "beta"]
    alpha = [c for c in fake.children() if str(place["alpha"]) in c]
    threads = sorted(c[c.index("--threads") + 1] for c in alpha)
    assert threads == ["0", "0", "0", "1", "1", "1"]


def test_every_cmake_call_and_child_gets_the_time_left_less_ten_seconds(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    write_cases(place)
    save_baseline(quick, place, use_runner)
    fake = FakeRunner()
    use_runner(fake)
    result = invoke(quick, place, "--budget", "300")
    assert result.exit_code == 0, result.output
    timed = [(args, t) for args, t in fake.calls
             if "_child" in args or Path(args[0]).name == "cmake"]  # fmt: skip
    assert timed and any(Path(a[0]).name == "cmake" for a, _ in timed)
    for args, timeout in timed:
        assert timeout is not None, args
        assert 0.0 < timeout <= 290.0, (args, timeout)


def test_a_timed_out_child_stops_the_run_out_of_time_naming_what_was_not_measured(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    write_cases(place)
    save_baseline(quick, place, use_runner)
    fake = FakeRunner(time_out=[str(place["alpha"])])
    use_runner(fake)
    result = invoke(quick, place, "--budget", "300")
    assert result.exit_code == 4, result.output
    verdict_line = last_line(result.output)
    assert verdict_line.startswith("OUT OF TIME: ")
    assert "alpha" in verdict_line and "beta" in verdict_line
    assert not fake.started(str(place["beta"]))


def test_a_case_whose_baseline_does_not_fit_is_skipped_and_the_rest_judged(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    """`beta`'s baseline median, 1000 s times its 3 runs, exceeds a 300 s budget."""
    write_cases(place)
    path = save_baseline(quick, place, use_runner)
    data = json.loads(path.read_text())
    data["cases"]["beta"]["measures"]["total"] = {"median": 1000.0, "min": 1000.0, "max": 1000.0}
    path.write_text(json.dumps(data))
    fake = FakeRunner()
    use_runner(fake)
    result = invoke(quick, place, "--budget", "300")
    assert result.exit_code == 4, result.output
    verdict_line = last_line(result.output)
    assert verdict_line.startswith("OUT OF TIME: ")
    assert "beta" in verdict_line and "alpha" not in verdict_line
    assert fake.started(str(place["alpha"]))
    assert not fake.started(str(place["beta"]))
    assert not [line for line in result.output.splitlines() if line.startswith("SLOWER: alpha")]


def test_a_slower_case_against_the_saved_baseline_exits_1(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    write_cases(place)
    save_baseline(quick, place, use_runner)
    use_runner(FakeRunner(totals=(1.5,)))
    result = invoke(quick, place)
    assert result.exit_code == 1, result.output
    assert last_line(result.output) == "SLOWER"
    assert [line for line in result.output.splitlines() if line.startswith("SLOWER: alpha total")]


def test_without_a_stored_baseline_the_run_exits_2(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    write_cases(place)
    use_runner(FakeRunner())
    result = invoke(quick, place)
    assert result.exit_code == 2, result.output
    assert last_line(result.output).startswith("NO BASELINE")
    assert not (place["quick"] / "baseline-ac.json").exists()


def test_a_changed_input_file_is_no_baseline(
    quick: ModuleType, place: dict[str, Path], use_runner: Callable[[FakeRunner], None]
) -> None:
    write_cases(place)
    save_baseline(quick, place, use_runner)
    place["beta"].write_bytes(b"another beta")
    use_runner(FakeRunner())
    result = invoke(quick, place)
    assert result.exit_code == 2, result.output
    assert last_line(result.output).startswith("NO BASELINE: ")


# ----------------------------------------------------------- one real child


def test_the_real_child_reports_mesh_clock_phases_and_total(tmp_path: Path) -> None:
    """Section 3: `bench.py _child` replaces `cli.PhaseClock` and keeps its
    instance, so its BENCH line carries `mesh`'s own phases and total
    (`src_python/tin_engine/cli.py@85a3e6dd:787`). A synthetic DEM, not the
    1 m tile: this proves the wiring, not a time."""
    pytest.importorskip("tin_engine._core")
    array = np.random.default_rng(14).uniform(0.0, 50.0, (17, 21)).astype(np.float32)
    dem = tmp_path / "bumpy.tif"
    dem.write_bytes(micro_tiff(array).getvalue())
    argv = [
        sys.executable, str(BENCH), "_child", "--threads", "0", "--",
        "mesh", "--dem", str(dem), "--tolerance", "1", "--out", str(tmp_path / "m.vtk"),
        "--binary",
    ]  # fmt: skip
    proc = subprocess.run(argv, capture_output=True, text=True, timeout=300, check=False)
    assert proc.returncode == 0, proc.stderr
    (line,) = [line for line in proc.stderr.splitlines() if line.startswith("BENCH ")]
    values = json.loads(line[len("BENCH ") :])
    phases = dict(values["phases"])
    assert "refine" in phases and "decode" in phases, phases
    # The clock starts inside `cli.app` and is read just after it returns: at
    # most the app's time, give or take 10 ms of reading order (a 17 x 21 DEM).
    assert 0.0 < values["total_s"] <= values["app_s"] + 0.01
    assert all(0.0 <= seconds <= values["total_s"] for seconds in phases.values()), phases
