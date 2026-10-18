"""What a run may use of the machine: the one module that reads it.

`docs/increments/38-threads-and-memory.md`, *4. The design*.

The core count, the physical memory and a container's cgroup limits are read
here and nowhere else; `Resources` carries the resolved figures, never 0, never None,
and every pool and `_core` call takes the one integer it uses.

The cgroup root is `/sys/fs/cgroup`, as a container runtime mounts the
container's own cgroup. A limit set on an ancestor cgroup of a process on a
bare host (a systemd slice) is not read. Pydantic and `os` only: no `_core`,
no typer.
"""

from __future__ import annotations

import os
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Self

from pydantic import BaseModel, ConfigDict, Field

CGROUP = Path("/sys/fs/cgroup")
VARIABLE = "RASPUTIN_THREADS"
_SUFFIX = {"": 1, "K": 1 << 10, "M": 1 << 20, "G": 1 << 30, "T": 1 << 40}


class Resources(BaseModel):
    """What a run may use of the machine, resolved: never 0, never None."""

    model_config = ConfigDict(frozen=True)

    #: The count every pool and `_core` call receives.
    threads: int = Field(ge=1)
    #: Bytes one DEM window may take, the flood's or the mosaic's.
    memory_budget: int = Field(ge=1)

    @classmethod
    def fallback(cls) -> Self:
        """A caller that gave nothing: one thread, half the machine's memory."""
        return cls(threads=1, memory_budget=machine_memory() // 2)


def _read(path: Path) -> str:
    """`path`'s text, stripped; empty when it is absent or unreadable."""
    try:
        return path.read_text(encoding="ascii").strip()
    except (OSError, UnicodeDecodeError):
        return ""


def _whole(word: str) -> int | None:
    """`word` as a whole number, decimal digits only, else None."""
    return int(word) if word.isascii() and word.isdigit() else None


def machine_cores(cgroup: Path = CGROUP) -> int:
    """The cores this process may use (affinity counted where Python reports
    it), capped by a cgroup v2 CPU quota: `<cgroup>/cpu.max` as `MAX PERIOD`,
    `ceil(MAX / PERIOD)`, counted only when both are positive whole numbers.
    Never under 1."""
    count = getattr(os, "process_cpu_count", os.cpu_count)() or 1
    words = _read(cgroup / "cpu.max").split()
    if len(words) == 2 and (quota := _whole(words[0])) and (period := _whole(words[1])):
        count = min(count, -(-quota // period))
    return max(1, int(count))


def machine_memory(cgroup: Path = CGROUP) -> int:
    """Physical memory in bytes, capped by the cgroup's limit: `memory.max`
    (v2), else `memory/memory.limit_in_bytes` (v1). A file that is absent,
    unreadable, reads `max` or holds no positive whole number is no limit;
    v1's "no limit" value (about 2^63) loses to the physical figure."""
    physical = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    for name in ("memory.max", "memory/memory.limit_in_bytes"):
        if limit := _whole(_read(cgroup / name)):
            return min(physical, limit)
    return physical


def parse_bytes(text: str) -> int:
    """`text` as bytes: a whole number with an optional suffix K, M, G or T
    (either case), 1024-based; anything else, or 0, is a `ValueError`."""
    found = re.fullmatch(r"([0-9]+)([KMGTkmgt]?)", text)
    if found is None or int(found[1]) == 0:
        raise ValueError(
            f"{text!r} is not a number of bytes (a whole number above 0: 8G, 512M, 4096)"
        )
    return int(found[1]) * _SUFFIX[found[2].upper()]


def resources(
    threads: int | None, memory_budget: str | None, environ: Mapping[str, str]
) -> Resources:
    """Sections 4.2 and 4.3: `--threads` (0 = every core), else `RASPUTIN_THREADS`
    by the same rule (stripped; empty is unset), else every core; and
    `--memory-budget`, else half the machine's memory. A `ValueError` names
    the text it refuses."""
    if threads is None:
        given = environ.get(VARIABLE, "").strip()
        if given and _whole(given) is None:
            raise ValueError(
                f"{VARIABLE}={given} is not a whole number of threads (0 means every core)"
            )
        threads = int(given) if given else 0
    budget = machine_memory() // 2 if memory_budget is None else parse_bytes(memory_budget)
    return Resources(threads=threads or machine_cores(), memory_budget=budget)
