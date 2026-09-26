"""Mesh statistics for ``rasputin mesh --stats``: sizes, quality and timings.

``docs/increments/17-mesh-stats.md`` R1, R3 to R5. Pure: numpy only, no
``_core`` and no typer, so it is testable without the extension. ``cli.py``
builds a :class:`Report` from what it ran and writes :func:`render`'s text.

Results, not configuration, so frozen dataclasses rather than Pydantic, as
``elevation.Trimmed`` does.
"""

# The report is Markdown with a multiplication sign, ≥ and °.
# ruff: noqa: RUF001

from __future__ import annotations

import time
from collections.abc import Callable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

Phases = Sequence[tuple[str, float]]


class PhaseClock:
    """Append-only phase times. A repeated name accumulates; ``phases()`` is in
    first-seen order. ``now`` returns integer nanoseconds (R5)."""

    def __init__(self, now: Callable[[], int] = time.perf_counter_ns) -> None:
        self._now = now
        self._start = now()
        self._seconds: dict[str, float] = {}

    @contextmanager
    def phase(self, name: str) -> Iterator[None]:
        """Time the block as ``name``, recorded even when it raises."""
        t0 = self._now()
        try:
            yield
        finally:
            self.add(name, (self._now() - t0) / 1e9)

    def add(self, name: str, seconds: float) -> None:
        """Record a time measured elsewhere, e.g. inside ``refine``."""
        self._seconds[name] = self._seconds.get(name, 0.0) + seconds

    def elapsed(self) -> float:
        """Seconds since the clock was made: the report's total."""
        return (self._now() - self._start) / 1e9

    def phases(self) -> list[tuple[str, float]]:
        return list(self._seconds.items())


@dataclass(frozen=True, slots=True)
class Quality:
    """R3. Angles in degrees, shares as fractions in [0, 1]; degrees as ints."""

    angle_median: float
    angle_under_1: float
    angle_under_10: float
    angle_worst: float
    degree_median: int
    degree_p99: int
    degree_max: int
    degree_at_least_12: int
    degree_at_least_20: int


def quality(vertices: npt.ArrayLike, triangles: npt.ArrayLike) -> Quality:
    """Minimum angle per triangle in plan view (columns 0 and 1 only), and the
    vertex degree as incident triangles (R3).

    Each angle is ``atan2(|e x f|, e . f)``, exact where ``acos`` of a
    normalised dot loses everything below about 1e-4 rad; a zero-area triangle
    gets 0.
    """
    xy = np.asarray(vertices, dtype=np.float64)[:, :2]
    tris = np.asarray(triangles, dtype=np.int64)
    p = xy[tris]  # (T, 3, 2)
    corners = []
    for k in range(3):
        e = p[:, (k + 1) % 3] - p[:, k]
        f = p[:, (k + 2) % 3] - p[:, k]
        cross = np.abs(e[:, 0] * f[:, 1] - e[:, 1] * f[:, 0])
        corners.append(np.arctan2(cross, np.einsum("ij,ij->i", e, f)))
    angles = np.degrees(np.min(corners, axis=0))
    degree = np.bincount(tris.ravel(), minlength=len(xy))

    def observed(q: float) -> int:
        return int(np.percentile(degree, q, method="inverted_cdf"))

    return Quality(
        angle_median=float(np.median(angles)),
        angle_under_1=float(np.mean(angles < 1.0)),
        angle_under_10=float(np.mean(angles < 10.0)),
        angle_worst=float(angles.min()),
        degree_median=observed(50),
        degree_p99=observed(99),
        degree_max=int(degree.max()),
        degree_at_least_12=int(np.count_nonzero(degree >= 12)),
        degree_at_least_20=int(np.count_nonzero(degree >= 20)),
    )


@dataclass(frozen=True, slots=True)
class Refinement:
    """``refine``'s counters, copied out by ``cli.py`` so no ``_core`` type
    reaches this module. ``carved`` None omits its column, and so does a
    ``quality_*`` None (increment 20's start-quality pass) or a ``feet`` None
    (increment 20b's constraint feet)."""

    tolerance: float
    max_error: float
    rounds: int
    inserted: int
    flips: int
    uncovered: int
    carved: int | None = None
    quality_inserted: int | None = None
    quality_skipped: int | None = None
    feet: int | None = None


@dataclass(frozen=True, slots=True)
class Sizes:
    """R4's Sizes rows; a None row does not apply and is omitted.
    ``dem_nodes`` is (rows, cols), ``dem_spacing`` (dx, dy), ``files``
    (name, bytes)."""

    output_vertices: int
    output_triangles: int
    constraint_edges: int
    files: Sequence[tuple[str, int]]
    dem_nodes: tuple[int, int] | None = None
    dem_spacing: tuple[float, float] | None = None
    domain_vertices: int | None = None
    domain_holes: int | None = None
    start_vertices: int | None = None
    start_triangles: int | None = None
    dropped: int | None = None


@dataclass(frozen=True, slots=True)
class Report:
    command: str
    sizes: Sizes
    quality: Quality
    refinement: Refinement | None
    phases: Phases
    total: float
    stats_seconds: float
    threads: int | None = None


def _exact(value: float) -> str:
    """``value`` short where that loses nothing, else every digit, so a printed
    achieved error can never read as above the tolerance it met."""
    short = f"{value:g}"
    return short if float(short) == value else repr(value)


def _bytes(n: int) -> str:
    """Decimal units, one decimal above bytes."""
    if n < 1000:
        return f"{n} B"
    size = float(n)
    for unit in ("kB", "MB", "GB"):
        size /= 1000
        if size < 1000 or unit == "GB":
            break
    return f"{size:.1f} {unit}"


def _table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    return lines + ["| " + " | ".join(row) + " |" for row in rows]


def _sizes(s: Sizes) -> list[str]:
    rows: list[tuple[str, str]] = []
    if s.dem_nodes is not None and s.dem_spacing is not None:
        dx, dy = s.dem_spacing
        spacing = f"{dx:g}" if dx == dy else f"{dx:g} × {dy:g}"
        rows.append(("DEM nodes", f"{s.dem_nodes[0]} × {s.dem_nodes[1]} ({spacing} m)"))
    if s.domain_vertices is not None:
        holes = s.domain_holes or 0
        plural = "" if holes == 1 else "s"
        rows.append(("domain vertices", f"{s.domain_vertices} (1 ring, {holes} hole{plural})"))
    counts = (
        ("start vertices", s.start_vertices),
        ("start triangles", s.start_triangles),
        ("output vertices", s.output_vertices),
        ("output triangles", s.output_triangles),
        ("constraint edges", s.constraint_edges),
        ("vertices without data dropped", s.dropped),
    )
    rows += [(item, str(n)) for item, n in counts if n is not None]
    rows += [(name, _bytes(n)) for name, n in s.files]
    return _table(("item", "count"), rows)


def _quality(q: Quality) -> list[str]:
    angle = (
        "minimum angle",
        f"{q.angle_median:.2f}°",
        f"{100 * q.angle_under_1:.2f} %",
        f"{100 * q.angle_under_10:.2f} %",
        f"{q.angle_worst:.3g}°",
    )
    counts = [q.degree_median, q.degree_p99, q.degree_max]
    counts += [q.degree_at_least_12, q.degree_at_least_20]
    degree = ("vertex degree (triangles)", *map(str, counts))
    return [
        *_table(("metric", "median", "< 1°", "< 10°", "worst"), [angle]),
        "",
        *_table(("metric", "median", "p99", "max", "≥ 12", "≥ 20"), [degree]),
    ]


def _refinement(r: Refinement) -> list[str]:
    cells = [
        ("tolerance", f"{_exact(r.tolerance)} m"),
        ("achieved max error", f"{_exact(r.max_error)} m"),
        ("rounds", str(r.rounds)),
        ("inserted", str(r.inserted)),
        *([("carved", str(r.carved))] if r.carved is not None else []),
        ("flips", str(r.flips)),
        ("uncovered", str(r.uncovered)),
    ]
    for header, count in (
        ("quality inserted", r.quality_inserted),
        ("quality skipped", r.quality_skipped),
        ("feet", r.feet),
    ):
        if count is not None:
            cells.append((header, str(count)))
    return _table([h for h, _ in cells], [[c for _, c in cells]])


_TIMINGS = """\
Wall clock, `time.perf_counter_ns` (Python) and `std::chrono::steady_clock`
(inside `refine`), one run, no warm-up. Total is the `mesh` command body, from
argument checks to the last file written; interpreter start-up and imports are
not in it."""


def _timings(report: Report) -> list[str]:
    """A phase ``"<parent>: ..."`` with an earlier ``<parent>`` row is a sub-row:
    listed, not added to the total. "other" is what the top-level rows miss."""
    total = report.total
    seen: set[str] = set()
    top = 0.0
    subs = False
    for name, seconds in report.phases:
        if name.split(": ", 1)[0] in seen:
            subs = True
        else:
            top += seconds
        seen.add(name)

    def share(seconds: float) -> str:
        return f"{100 * seconds / total:.1f} %" if total > 0 else "0.0 %"

    rows = [(n, f"{s:.3f}", share(s)) for n, s in (*report.phases, ("other", total - top))]
    rows.append(("**total**", f"**{total:.3f}**", "**100 %**"))
    intro = _TIMINGS
    if report.threads is not None:
        intro += f" Threads: {report.threads} (hardware concurrency)."
    lines = ["## Timings", "", intro, "", *_table(("phase", "seconds", "share"), rows), ""]
    if subs:
        lines += ["Sub-rows sum to their parent and are not added to the total.", ""]
    return [*lines, f"Statistics computed in {report.stats_seconds:.3f} s, not included above."]


def render(report: Report) -> str:
    """The Markdown report (R4). Pure: no clock, no I/O."""
    lines = ["# rasputin mesh — statistics", "", f"`{report.command}`", ""]
    lines += ["## Sizes", "", *_sizes(report.sizes), ""]
    lines += ["## Quality (plan view, x/y)", "", *_quality(report.quality), ""]
    if report.refinement is not None:
        lines += ["## Refinement", "", *_refinement(report.refinement), ""]
    lines += _timings(report)
    return "\n".join(lines) + "\n"
