"""`tin_engine.resources`: the one module that reads the machine (increment 38).

`docs/increments/38-threads-and-memory.md`, sections 4.1 to 4.3, 6 and 8
(T1 to T5, and T11's `decompose`). The module is fetched inside a fixture, so
its absence fails these tests and leaves the rest of the session collecting.

Figures that depend on the machine (its core count, its memory) are computed
from the module's own readers with no cgroup, never written as literals: the
design's "10 cores, 32 GiB" are this machine's and the suite runs on others.
The cgroup readers take the cgroup root, so a `tmp_path` with a synthetic
`cpu.max` or `memory.max` stands for a container; the scale assumed is one
container with its own cgroup root (the design's stated limit).

T11's `project_structure.md` assertion is not here: that file is prose under
`tools/ci_changes.py`'s rule, and `conftest.py` fails a session in which a
test opens one (h11, T5). `@reviewer` reads the two rows instead.
"""

from __future__ import annotations

import importlib
import importlib.util
import os
from pathlib import Path
from types import ModuleType

import pytest
from pydantic import ValidationError

import tin_engine

PACKAGE = Path(tin_engine.__file__).resolve().parent

#: The tokens by which Python reads the machine; `resources.py` alone may use them.
MACHINE_TOKENS = ("cpu_count", "sched_getaffinity", "sysconf", "/sys/fs/cgroup")

GIB = 1 << 30
#: cgroup v1's "no limit": PAGE_COUNTER_MAX, about 2^63, the figure a host prints.
V1_UNLIMITED = "9223372036854771712"


@pytest.fixture(scope="module")
def res() -> ModuleType:
    return importlib.import_module("tin_engine.resources")


@pytest.fixture
def no_cgroup(tmp_path: Path) -> Path:
    """A cgroup root with no files: the bare host."""
    root = tmp_path / "empty"
    root.mkdir()
    return root


def physical_memory() -> int:
    return os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")


# ---------------------------------------------------------------- T1


def test_t1_only_resources_reads_the_machine() -> None:
    """grep-equivalent over the package: each token occurs in `resources.py` only."""
    readers = {
        str(path.relative_to(PACKAGE))
        for path in PACKAGE.rglob("*.py")
        if any(token in path.read_text(encoding="utf-8") for token in MACHINE_TOKENS)
    }
    assert readers == {"resources.py"}


# ---------------------------------------------------------------- the model


class TestResourcesModel:
    def test_frozen_and_never_zero(self, res: ModuleType) -> None:
        r = res.Resources(threads=3, memory_budget=4096)
        assert (r.threads, r.memory_budget) == (3, 4096)
        with pytest.raises(ValidationError):
            r.threads = 4
        for bad in ({"threads": 0, "memory_budget": 1}, {"threads": 1, "memory_budget": 0}):
            with pytest.raises(ValidationError):
                res.Resources(**bad)

    def test_fallback_is_one_thread_and_half_the_memory(self, res: ModuleType) -> None:
        fallback = res.Resources.fallback()
        assert fallback.threads == 1
        assert fallback.memory_budget == res.machine_memory() // 2


# ---------------------------------------------------------------- T2


class TestT2TheThreadsRule:
    """Section 4.2's table, one case each; `machine_cores()` is the machine's."""

    @pytest.mark.parametrize(
        ("threads", "environ"),
        [
            pytest.param(None, {}, id="nothing given"),
            pytest.param(0, {}, id="--threads 0"),
            pytest.param(None, {"RASPUTIN_THREADS": "0"}, id="RASPUTIN_THREADS=0"),
            pytest.param(None, {"RASPUTIN_THREADS": ""}, id="RASPUTIN_THREADS empty"),
            pytest.param(0, {"RASPUTIN_THREADS": "7"}, id="--threads 0 beats the variable"),
        ],
    )
    def test_every_core(
        self, res: ModuleType, threads: int | None, environ: dict[str, str]
    ) -> None:
        assert res.resources(threads, None, environ).threads == res.machine_cores()

    @pytest.mark.parametrize(
        ("threads", "environ", "expected"),
        [
            pytest.param(3, {"RASPUTIN_THREADS": "7"}, 3, id="the option beats the variable"),
            pytest.param(None, {"RASPUTIN_THREADS": "7"}, 7, id="the variable alone"),
            pytest.param(None, {"RASPUTIN_THREADS": " 4 "}, 4, id="the variable is stripped"),
            pytest.param(64, {}, 64, id="above the core count, honoured"),
            pytest.param(None, {"RASPUTIN_THREADS": "64"}, 64, id="above the cores, by variable"),
        ],
    )
    def test_a_count(
        self, res: ModuleType, threads: int | None, environ: dict[str, str], expected: int
    ) -> None:
        assert res.resources(threads, None, environ).threads == expected

    @pytest.mark.parametrize("text", ["abc", "-1", "1.5", "0x10"])
    def test_a_variable_that_is_not_a_whole_number_is_refused(
        self, res: ModuleType, text: str
    ) -> None:
        with pytest.raises(ValueError, match="RASPUTIN_THREADS") as info:
            res.resources(None, None, {"RASPUTIN_THREADS": text})
        assert text in str(info.value)
        assert "0 means every core" in str(info.value)


# ---------------------------------------------------------------- T3


class TestT3MachineCores:
    """`machine_cores(cgroup)`: the process's count, capped by a cgroup v2 CPU
    quota `ceil(MAX / PERIOD)`; anything that is not two positive whole numbers
    is no quota."""

    def test_the_plain_count_is_the_processs(self, res: ModuleType, no_cgroup: Path) -> None:
        expected = getattr(os, "process_cpu_count", os.cpu_count)() or 1
        assert res.machine_cores(no_cgroup) == expected
        assert res.machine_cores() >= 1

    def test_a_quota_of_two_cpus_caps_at_two(self, res: ModuleType, tmp_path: Path) -> None:
        (tmp_path / "cpu.max").write_text("200000 100000\n")
        assert res.machine_cores(tmp_path) == min(2, res.machine_cores(tmp_path / "none"))

    def test_half_a_cpu_rounds_up_to_one(self, res: ModuleType, tmp_path: Path) -> None:
        (tmp_path / "cpu.max").write_text("50000 100000\n")
        assert res.machine_cores(tmp_path) == 1

    def test_one_and_a_half_cpus_round_up_to_two(self, res: ModuleType, tmp_path: Path) -> None:
        (tmp_path / "cpu.max").write_text("150000 100000\n")
        assert res.machine_cores(tmp_path) == min(2, res.machine_cores(tmp_path / "none"))

    @pytest.mark.parametrize(
        "text",
        [
            pytest.param("max 100000\n", id="no quota"),
            pytest.param("abc 100000\n", id="first word not a number"),
            pytest.param("200000 0\n", id="zero period"),
            pytest.param("200000\n", id="one word"),
            pytest.param("200000 abc\n", id="period not a number"),
            pytest.param("", id="empty file"),
            pytest.param("-200000 100000\n", id="negative quota"),
        ],
    )
    def test_no_quota_gives_the_plain_count(
        self, res: ModuleType, tmp_path: Path, no_cgroup: Path, text: str
    ) -> None:
        (tmp_path / "cpu.max").write_text(text)
        assert res.machine_cores(tmp_path) == res.machine_cores(no_cgroup)

    @pytest.mark.skipif(os.geteuid() == 0, reason="root reads every file")
    def test_an_unreadable_file_is_no_quota(
        self, res: ModuleType, tmp_path: Path, no_cgroup: Path
    ) -> None:
        path = tmp_path / "cpu.max"
        path.write_text("200000 100000\n")
        path.chmod(0)
        try:
            assert res.machine_cores(tmp_path) == res.machine_cores(no_cgroup)
        finally:
            path.chmod(0o600)


# ---------------------------------------------------------------- T4


class TestT4MachineMemory:
    """`machine_memory(cgroup)`: physical memory, capped by `memory.max` (v2)
    or `memory/memory.limit_in_bytes` (v1). Scale: the limits below are 1 GiB
    and 512 MiB, under any machine this suite runs on (asserted)."""

    def test_the_plain_figure_is_physical_memory(self, res: ModuleType, no_cgroup: Path) -> None:
        assert res.machine_memory(no_cgroup) == physical_memory()
        assert res.machine_memory() >= 1

    def test_a_v2_limit_under_the_machine_wins(self, res: ModuleType, tmp_path: Path) -> None:
        assert physical_memory() > GIB, "the fixture assumes a machine over 1 GiB"
        (tmp_path / "memory.max").write_text(f"{GIB}\n")
        assert res.machine_memory(tmp_path) == GIB

    def test_a_v1_limit_under_the_machine_wins(self, res: ModuleType, tmp_path: Path) -> None:
        (tmp_path / "memory").mkdir()
        (tmp_path / "memory" / "memory.limit_in_bytes").write_text(f"{GIB // 2}\n")
        assert res.machine_memory(tmp_path) == GIB // 2

    def test_v1s_no_limit_value_loses_to_the_physical_figure(
        self, res: ModuleType, tmp_path: Path
    ) -> None:
        (tmp_path / "memory").mkdir()
        (tmp_path / "memory" / "memory.limit_in_bytes").write_text(f"{V1_UNLIMITED}\n")
        assert res.machine_memory(tmp_path) == physical_memory()

    def test_v2_wins_when_both_exist(self, res: ModuleType, tmp_path: Path) -> None:
        (tmp_path / "memory.max").write_text(f"{GIB}\n")
        (tmp_path / "memory").mkdir()
        (tmp_path / "memory" / "memory.limit_in_bytes").write_text(f"{GIB // 2}\n")
        assert res.machine_memory(tmp_path) == GIB

    @pytest.mark.parametrize(
        "text",
        [
            pytest.param("max\n", id="max"),
            pytest.param("0\n", id="zero"),
            pytest.param("-1\n", id="negative"),
            pytest.param("abc\n", id="text"),
            pytest.param("", id="empty file"),
        ],
    )
    def test_no_limit_gives_the_physical_figure(
        self, res: ModuleType, tmp_path: Path, text: str
    ) -> None:
        (tmp_path / "memory.max").write_text(text)
        assert res.machine_memory(tmp_path) == physical_memory()


# ---------------------------------------------------------------- T5


class TestT5ParseBytes:
    @pytest.mark.parametrize(
        ("text", "expected"),
        [
            ("8G", 8 * GIB),
            ("512m", 512 << 20),
            ("1T", 1 << 40),
            ("4096", 4096),
            ("2k", 2048),
        ],
    )
    def test_a_whole_number_with_an_optional_1024_suffix(
        self, res: ModuleType, text: str, expected: int
    ) -> None:
        assert res.parse_bytes(text) == expected

    @pytest.mark.parametrize("text", ["0", "8GB", "-1", "abc", "", "1.5G", "G"])
    def test_anything_else_is_refused_naming_the_text(self, res: ModuleType, text: str) -> None:
        with pytest.raises(ValueError) as info:
            res.parse_bytes(text)
        assert (text or repr(text)) in str(info.value)

    def test_the_default_budget_is_half_the_machines_memory(self, res: ModuleType) -> None:
        assert res.resources(None, None, {}).memory_budget == res.machine_memory() // 2

    def test_a_given_budget_is_parsed(self, res: ModuleType) -> None:
        assert res.resources(None, "8G", {}).memory_budget == 8 * GIB
        assert res.resources(None, "4096", {}).memory_budget == 4096

    @pytest.mark.parametrize("text", ["0", "8GB", "-1", "abc"])
    def test_a_bad_budget_is_refused_naming_the_text(self, res: ModuleType, text: str) -> None:
        with pytest.raises(ValueError) as info:
            res.resources(None, text, {})
        assert text in str(info.value)


# ---------------------------------------------------------------- T11


def test_t11_decompose_is_gone() -> None:
    assert importlib.util.find_spec("tin_engine.decompose") is None
