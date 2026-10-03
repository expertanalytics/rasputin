"""The acceptance benchmark: the 1 m run and the thread-scaling sweep.

Design: ``docs/benchmarks/bench-py.md``. The rule it serves: "Acceptance: an
increment that touches refine or mesh code" in ``docs/increments/README.md``.

``run`` builds the tree's ``_core`` in Release (``<tree>/build-bench``), times
one child process per sample, measures the mesh quality, reads the power state
before and after, writes ``<out-root>/<date>/<label>/`` and judges the run
against the newest comparable stored baseline. ``compare`` re-judges stored
evidence. ``_child`` is the process each sample runs in; it is dispatched before
Typer so the rasputin argv after ``--`` reaches the CLI untouched.

Every subprocess goes through a :class:`Runner`, so the tests replace
:func:`make_runner` and start none.
"""

from __future__ import annotations

import hashlib
import json
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import warnings
from collections.abc import Callable, Mapping, Sequence
from datetime import datetime
from fractions import Fraction
from pathlib import Path
from typing import Annotated, Any, Literal, NamedTuple, Protocol

import numpy as np
import numpy.typing as npt
import typer
from pydantic import BaseModel

REPO = Path(__file__).resolve().parents[1]
BENCH = Path(__file__).resolve()
DEFAULT_DEM = REPO / "tests/fixtures/dem_archive/7908_3_10m_z33.tif"
DEFAULT_QUARTER = REPO / "docs/benchmarks/2026-09-26/quarter.geojson"
#: The 2026-09-26 scaling ceiling each run's is reported against.
CEILING_2026_09_26 = "2.2x"
#: The README's generated part ends at this line; what follows it is @perf's.
MARKER = "<!-- bench.py: generated above this line; hand-written prose below is kept -->"
PMSET = ["pmset", "-g", "batt"]
CHILD_KEYS = ("refine_s", "app_s", "max_error", "rounds", "inserted", "flips")

State = Literal["ac", "battery", "mixed", "unknown"]


# ------------------------------------------------------------ subprocesses
class Completed(NamedTuple):
    returncode: int
    stdout: str
    stderr: str
    wall_s: float


class Runner(Protocol):
    def run(self, argv: Sequence[str], cwd: Path | None = None) -> Completed: ...


class SubprocessRunner:
    def run(self, argv: Sequence[str], cwd: Path | None = None) -> Completed:
        t0 = time.perf_counter()
        proc = subprocess.run(
            [str(a) for a in argv], cwd=cwd, capture_output=True, text=True, check=False
        )
        return Completed(proc.returncode, proc.stdout, proc.stderr, time.perf_counter() - t0)


def make_runner() -> Runner:
    """The factory the tests replace with a fake."""
    return SubprocessRunner()


class BenchError(Exception):
    """A run that cannot produce evidence: exit 3, nothing written."""


# ------------------------------------------------------------------ models
class Power(BaseModel):
    state: State
    percent: int | None
    raw: str


class Sample(BaseModel):
    domain: str
    threads: int
    repeat: int
    refine_s: float
    app_s: float
    proc_s: float
    max_error: float
    rounds: int
    inserted: int
    flips: int


class TimeStats(BaseModel):
    domain: str
    threads: int
    n: int
    median: float
    min: float
    max: float


class MeshQuality(BaseModel):
    worst_angle: float
    angle_median: float
    share_under_1: float
    max_degree: int
    within_tolerance: bool
    delaunay_checked: int
    delaunay_ambiguous: int
    delaunay_violations: int
    mesh_sha256: str = ""


class Tree(BaseModel):
    commit: str
    dirty: bool


class Build(BaseModel):
    no_build: bool
    type: str | None = None
    cxx_flags_release: str | None = None
    compiler: str | None = None
    so_sha256: str | None = None


class Machine(BaseModel):
    cpu_brand: str
    p_cores: int
    e_cores: int
    memory_bytes: int
    macos: str
    python: str
    numpy: str


class Domain(BaseModel):
    name: str
    path: str | None
    sha256: str | None


class Inputs(BaseModel):
    dem: str
    dem_sha256: str
    domains: list[Domain]
    tolerance: float
    extra_args: list[str]


class RunRecord(BaseModel):
    label: str
    started: datetime
    tree: Tree
    bench_blob: str
    child_argv: list[str]
    build: Build
    machine: Machine
    power: Power
    inputs: Inputs
    samples: list[Sample]
    stats: list[TimeStats]
    quality: dict[str, MeshQuality]
    accept_quality: bool = False
    threshold_pct: float = 5.0
    verdict: list[str] = []
    #: The bounds checks the children's ``_core`` reported (increment 24); a
    #: run.json from before it loads as "none", which every such run was.
    hardening: str = "none"


class ChildError(ValueError):
    """A child's stderr without exactly one well-formed ``BENCH`` line."""


class Child(NamedTuple):
    refine_s: float
    app_s: float
    max_error: float
    rounds: int
    inserted: int
    flips: int
    hardening: str = "none"


class Ceiling(NamedTuple):
    top_threads: int
    at_top: float
    best: float
    best_threads: int


class VtkMesh(NamedTuple):
    points: npt.NDArray[np.float64]
    triangles: npt.NDArray[np.int64]
    edges: npt.NDArray[np.int64]
    sha256: str


class Verdict(NamedTuple):
    status: str
    lines: list[str]
    exit_code: int


# -------------------------------------------------------------- pure parts
def parse_pmset(text: str) -> Power:
    """``pmset -g batt``: AC and battery are the two states the rule compares;
    UPS, empty and anything else is ``unknown``."""
    source = re.search(r"Now drawing from '([^']*)'", text)
    state: State = {"AC Power": "ac", "Battery Power": "battery"}.get(  # type: ignore[assignment]
        source.group(1) if source else "", "unknown"
    )
    percent = re.search(r"\t(\d+)%;", text)
    return Power(state=state, percent=int(percent.group(1)) if percent else None, raw=text)


def combine_power(before: Power, after: Power) -> Power:
    state: State = before.state if before.state == after.state else "mixed"
    raw = f"{before.raw}--- after ---\n{after.raw}"
    return Power(state=state, percent=after.percent, raw=raw)


def parse_child(stderr: str, run: str) -> Child:
    lines = [line[6:] for line in stderr.splitlines() if line.startswith("BENCH ")]
    if len(lines) != 1:
        raise ChildError(f"{run}: expected one BENCH line on stderr, found {len(lines)}")
    try:
        values = json.loads(lines[0])
        return Child(
            float(values["refine_s"]), float(values["app_s"]), float(values["max_error"]),
            int(values["rounds"]), int(values["inserted"]), int(values["flips"]),
            str(values.get("hardening", "none")),
        )  # fmt: skip
    except (ValueError, KeyError, TypeError) as exc:
        raise ChildError(f"{run}: malformed BENCH line: {exc!r}") from exc


def median_stats(samples: Sequence[Sample]) -> list[TimeStats]:
    groups: dict[tuple[str, int], list[float]] = {}
    for s in samples:
        groups.setdefault((s.domain, s.threads), []).append(s.refine_s)
    return [
        TimeStats(
            domain=d, threads=t, n=len(v), median=statistics.median(v), min=min(v), max=max(v)
        )
        for (d, t), v in sorted(groups.items())
    ]


def ceiling(medians: Mapping[int, float]) -> Ceiling:
    """Speed-up over 1 thread; key 0, the CLI's default, is not a count."""
    one = medians[1]
    counts = sorted(k for k in medians if k != 0)
    speedups = {k: one / medians[k] for k in counts}
    best_threads = max(counts, key=lambda k: speedups[k])
    top = counts[-1]
    return Ceiling(top, speedups[top], speedups[best_threads], best_threads)


def _incircle_exact(a: Any, b: Any, c: Any, d: Any) -> Fraction:
    (ax, ay), (bx, by), (cx, cy) = (
        (Fraction(p[0]) - Fraction(d[0]), Fraction(p[1]) - Fraction(d[1])) for p in (a, b, c)
    )
    al, bl, cl = ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy
    return ax * (by * cl - bl * cy) - ay * (bx * cl - bl * cx) + al * (bx * cy - by * cx)


def _delaunay(
    xy: npt.NDArray[np.float64], tris: npt.NDArray[np.int64], cons: npt.NDArray[np.int64]
) -> tuple[int, int, int]:
    """Per interior non-constraint edge: is one side's apex strictly inside the
    other side's circumcircle? ``quality.py``'s test, 2026-09-26: a float
    determinant unless within 4x Shewchuk's incircle bound, then exact."""
    a, b, c = xy[tris[:, 0]], xy[tris[:, 1]], xy[tris[:, 2]]
    orient = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    ccw = tris.copy()
    ccw[orient < 0] = ccw[orient < 0][:, [0, 2, 1]]
    u, v, w = ccw.ravel(), ccw[:, [1, 2, 0]].ravel(), ccw[:, [2, 0, 1]].ravel()
    owner = np.repeat(np.arange(len(ccw)), 3)
    key = np.minimum(u, v) * (len(xy) + 1) + np.maximum(u, v)
    order = np.argsort(key, kind="stable")
    same = np.nonzero(key[order][1:] == key[order][:-1])[0]
    h1, h2 = order[same], order[same + 1]
    ckey = np.minimum(cons[:, 0], cons[:, 1]) * (len(xy) + 1) + np.maximum(cons[:, 0], cons[:, 1])
    keep = ~np.isin(key[h1], ckey)
    h1, h2 = h1[keep], h2[keep]
    tri, apex = ccw[owner[h1]], w[h2]
    ad, bd, cd = (xy[tri[:, k]] - xy[apex] for k in range(3))
    al, bl, cl = ((p**2).sum(1) for p in (ad, bd, cd))
    x1 = bd[:, 0] * cd[:, 1] - bd[:, 1] * cd[:, 0]
    x2 = cd[:, 0] * ad[:, 1] - cd[:, 1] * ad[:, 0]
    x3 = ad[:, 0] * bd[:, 1] - ad[:, 1] * bd[:, 0]
    det = al * x1 + bl * x2 + cl * x3
    eps = np.finfo(float).eps / 2
    bound = 4 * (10 + 96 * eps) * eps * (al * np.abs(x1) + bl * np.abs(x2) + cl * np.abs(x3))
    inside = det > bound
    ambiguous = np.nonzero(np.abs(det) <= bound)[0]
    for k in ambiguous.tolist():
        inside[k] = _incircle_exact(xy[tri[k, 0]], xy[tri[k, 1]], xy[tri[k, 2]], xy[apex[k]]) > 0
    return len(h1), len(ambiguous), int(inside.sum())


def quality(
    points: npt.ArrayLike, triangles: npt.ArrayLike, constraint_edges: npt.ArrayLike,
    tolerance: float, max_error: float,
) -> MeshQuality:  # fmt: skip
    """Angles and degree as ``rasputin mesh --stats`` computes them, the
    tolerance check, and the constrained Delaunay check."""
    from tin_engine.stats import quality as stats_quality  # lazy: the child sets sys.path first

    xy = np.asarray(points, dtype=np.float64)[:, :2]
    tris = np.asarray(triangles, dtype=np.int64).reshape(-1, 3)
    cons = np.asarray(constraint_edges).astype(np.int64).reshape(-1, 2)
    q = stats_quality(xy, tris)
    checked, ambiguous, violations = _delaunay(xy, tris, cons)
    return MeshQuality(
        worst_angle=q.angle_worst, angle_median=q.angle_median, share_under_1=q.angle_under_1,
        max_degree=q.degree_max, within_tolerance=max_error <= tolerance,
        delaunay_checked=checked, delaunay_ambiguous=ambiguous, delaunay_violations=violations,
    )  # fmt: skip


def read_vtk_ascii(path: Path) -> VtkMesh:
    """rasputin's legacy ASCII POLYDATA; ``sha256`` covers the bytes from
    ``POINTS`` on, so header fields (CRS, title) do not change it."""
    data = path.read_bytes()
    lines = data.decode("latin-1").split("\n")
    if len(lines) > 2 and lines[2].strip() == "BINARY":
        raise ValueError(f"{path}: a BINARY VTK file; the quality run writes --ascii")
    sections: dict[str, list[str]] = {}
    for i, line in enumerate(lines):
        head = line.split()
        if head and head[0] in ("POINTS", "LINES", "POLYGONS") and head[0] not in sections:
            sections[head[0]] = lines[i + 1 : i + 1 + int(head[1])]
    if "POINTS" not in sections:
        raise ValueError(f"{path}: no POINTS section")
    table = {k: " ".join(sections.get(k, [])).split() for k in ("POINTS", "LINES", "POLYGONS")}
    return VtkMesh(
        points=np.array(table["POINTS"], dtype=np.float64).reshape(-1, 3),
        triangles=np.array(table["POLYGONS"], dtype=np.int64).reshape(-1, 4)[:, 1:4],
        edges=np.array(table["LINES"], dtype=np.int64).reshape(-1, 3)[:, 1:3],
        sha256=hashlib.sha256(data[data.index(b"\nPOINTS ") + 1 :]).hexdigest(),
    )


# --------------------------------------------------------- comparison
def comparable(a: RunRecord, b: RunRecord) -> str | None:
    """None, or the first field that forbids comparing ``a`` with ``b``."""
    if a.power.state not in ("ac", "battery") or a.power.state != b.power.state:
        return f"power: {a.power.state} vs {b.power.state}"
    pairs: list[tuple[str, object, object]] = [
        ("cpu_brand", a.machine.cpu_brand, b.machine.cpu_brand),
        ("p_cores", a.machine.p_cores, b.machine.p_cores),
        ("e_cores", a.machine.e_cores, b.machine.e_cores),
        ("dem_sha256", a.inputs.dem_sha256, b.inputs.dem_sha256),
        ("domains", sorted((d.name, d.sha256) for d in a.inputs.domains),
         sorted((d.name, d.sha256) for d in b.inputs.domains)),
        ("tolerance", a.inputs.tolerance, b.inputs.tolerance),
        ("extra_args", a.inputs.extra_args, b.inputs.extra_args),
        ("hardening", a.hardening, b.hardening),
    ]  # fmt: skip
    for name, x, y in pairs:
        if x != y:
            return f"{name}: {x} vs {y}"
    return None


def find_baseline(
    root: Path, record: RunRecord, is_ancestor: Callable[[str, str], bool]
) -> tuple[Path, RunRecord] | None:
    """The newest ``run.json`` under ``root`` (any depth) older than ``record``,
    comparable with it, and whose commit is an ancestor of ``record``'s."""
    found: tuple[Path, RunRecord] | None = None
    for path in sorted(root.rglob("run.json")):
        try:
            candidate = _load(path.parent)
        except BenchError as exc:  # skipped, not fatal: bench-py.md "Bad input exits 3"
            warnings.warn(f"skipped {exc}", UserWarning, stacklevel=2)
            continue
        if candidate.started >= record.started or comparable(record, candidate) is not None:
            continue
        if found is not None and candidate.started <= found[1].started:
            continue
        if is_ancestor(candidate.tree.commit, record.tree.commit):
            found = (path.parent, candidate)
    return found


def verdict(new: RunRecord, base: RunRecord | None) -> Verdict:
    reason = "no comparable stored run" if base is None else comparable(new, base)
    if reason is not None or base is None:
        return Verdict("NO BASELINE", [f"NO BASELINE: {reason}"], 2)
    lines: list[str] = []
    before = {(s.domain, s.threads): s.median for s in base.stats}
    for s in new.stats:
        old = before.get((s.domain, s.threads))
        if old is None or old <= 0:
            continue
        pct = (s.median / old - 1.0) * 100.0
        # Exact products, not the rounded ratio: exactly the threshold is accepted.
        if s.median * 100.0 > old * (100.0 + new.threshold_pct):
            where = f"{s.domain} refine_s[t={s.threads}]"
            lines.append(f"REGRESSION: {where} {old:.4f} -> {s.median:.4f} (+{pct:.1f} %)")
    for name, q in new.quality.items():
        b = base.quality.get(name)
        if b is None:
            continue
        if not q.within_tolerance:
            lines.append(f"REGRESSION: {name} tolerance {b.within_tolerance} -> False")
        if q.delaunay_violations > 0:
            n = f"{b.delaunay_violations} -> {q.delaunay_violations}"
            lines.append(f"REGRESSION: {name} delaunay_violations {n}")
        waivable = [
            ("worst_angle", b.worst_angle, q.worst_angle, q.worst_angle < b.worst_angle),
            ("max_degree", b.max_degree, q.max_degree, q.max_degree > b.max_degree),
        ]
        for measure, old_q, new_q, worse in waivable:
            if worse:
                prefix = "WAIVED (--accept-quality): " if new.accept_quality else "REGRESSION: "
                lines.append(f"{prefix}{name} {measure} {old_q} -> {new_q}")
    if any(line.startswith("REGRESSION: ") for line in lines):
        return Verdict("REGRESSION", lines, 1)
    return Verdict("ACCEPTED", ["ACCEPTED", *lines], 0)


# ---------------------------------------------------------------- evidence
def _readme(r: RunRecord) -> list[str]:
    m, b, i, p = r.machine, r.build, r.inputs, r.power
    build = (
        "none (--no-build: the installed tin_engine)"
        if b.no_build
        else (f"{b.type}, {b.compiler}, `{b.cxx_flags_release}`, _core sha256 `{b.so_sha256}`")
    )
    out = [
        f"# Benchmark run `{r.label}`", "", "Generated by `tools/bench.py`.", "", "## Verdict", "",
        *[f"- {line}" for line in r.verdict], "", f"Time threshold: {r.threshold_pct} %.",
        *(["`--accept-quality` waived worst angle and max degree loss (named above)."]
          if r.accept_quality else []),
        "", "## Method", "",
        f"- Started {r.started.isoformat()}; tree `{r.tree.commit}`{' (dirty)' * r.tree.dirty}; "
        f"bench.py blob `{r.bench_blob}`. Build: {build}; bounds checks: {r.hardening}.",
        f"- {m.cpu_brand}, {m.p_cores} P + {m.e_cores} E cores, {m.memory_bytes / 2**30:.0f} GiB, "
        f"macOS {m.macos}, Python {m.python}, numpy {m.numpy}.",
        f"- Power **{p.state}** ({p.percent}%), `pmset -g batt` before and after (in run.json).",
        f"- DEM `{i.dem}` (sha256 `{i.dem_sha256}`), tolerance {i.tolerance}, extra mesh args "
        f"`{' '.join(i.extra_args)}`; domains: "
        + ", ".join(f"`{d.name}` ({d.path or 'the whole tile'})" for d in i.domains) + ".",
        f"- One child per sample, `{' '.join(r.child_argv)}`, repeats interleaved over thread "
        "counts; t=0 is the CLI's default, other counts are forced into `refine`. Quality: one "
        "`--ascii` run per domain at t=0, kept out of the repository; rerun the child with "
        "`--ascii --out PATH` to regenerate it.",
    ]  # fmt: skip
    for d in i.domains:
        rows = [s for s in r.stats if s.domain == d.name]
        medians = {s.threads: s.median for s in rows}
        c = ceiling(medians) if 1 in medians else None
        q = r.quality.get(d.name)
        out += [
            "", f"## `{d.name}`", "", "| threads | n | median s | min s | max s |",
            "|---:|---:|---:|---:|---:|",
            *(f"| {s.threads} | {s.n} | {s.median:.4f} | {s.min:.4f} | {s.max:.4f} |"
              for s in rows),
            "", (f"Ceiling: {c.at_top:.2f}x at {c.top_threads} threads over 1, best {c.best:.2f}x "
                 f"at {c.best_threads}" if c else "Ceiling: no 1-thread run")
            + f"; 2026-09-26: about {CEILING_2026_09_26}, flat from about 7.",
        ]  # fmt: skip
        if q is not None:
            out += ["", f"Quality: worst angle {q.worst_angle:.4f} deg, median {q.angle_median:.2f}"
                    f", share under 1 deg {q.share_under_1:.5f}, max degree {q.max_degree}, "
                    f"within tolerance {q.within_tolerance}, Delaunay {q.delaunay_violations} "
                    f"violations of {q.delaunay_checked} edges ({q.delaunay_ambiguous} decided "
                    f"exactly), mesh sha256 `{q.mesh_sha256}`."]  # fmt: skip
    return out


def write_evidence(record: RunRecord, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "run.json").write_text(record.model_dump_json(indent=1) + "\n")
    rows = [
        f"{s.domain}\t{s.threads}\t{s.repeat}\t{s.refine_s:.6f}\t{record.power.state}\n"
        for s in record.samples
    ]
    (directory / "raw.tsv").write_text("".join(rows))
    readme = directory / "README.md"
    kept = readme.read_text().partition(MARKER + "\n")[2] if readme.exists() else ""
    readme.write_text("\n".join([*_readme(record), "", MARKER]) + "\n" + kept)


# -------------------------------------------------------------------- run
def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _cache_value(cache: Path, key: str) -> str | None:
    found = re.search(rf"^{key}:\w+=(.*)$", cache.read_text(), re.M) if cache.exists() else None
    return found.group(1) if found else None


def _checked(runner: Runner, argv: Sequence[str], cwd: Path | None = None) -> Completed:
    done = runner.run(argv, cwd)
    if done.returncode != 0:
        raise BenchError(f"{' '.join(map(str, argv))} failed: {done.stderr.strip()[-2000:]}")
    return done


def build(runner: Runner, tree: Path, hardening: str = "on") -> tuple[Path, Build]:
    """Release ``_core`` in ``<tree>/build-bench``, and ``build-bench/pkg``: the
    tree's ``tin_engine`` as symlinks plus a copy of the fresh ``.so``. The
    hardening option is passed on every configure, so a cached one never decides."""
    out = tree / "build-bench"
    cache = out / "CMakeCache.txt"
    kind = _cache_value(cache, "CMAKE_BUILD_TYPE")
    if cache.exists() and kind != "Release":
        raise BenchError(f"{cache}: CMAKE_BUILD_TYPE is {kind or 'unset'}, not Release; refused")
    _checked(runner, [
        "cmake", "-S", str(tree), "-B", str(out), "-DCMAKE_BUILD_TYPE=Release",
        "-DRASPUTIN_BUILD_PYTHON=ON", "-DRASPUTIN_BUILD_TESTS=OFF",
        f"-DPython_EXECUTABLE={sys.executable}", "-DPYBIND11_FINDPYTHON=ON",
        f"-DRASPUTIN_HARDENING={hardening.upper()}",
    ])  # fmt: skip
    _checked(runner, ["cmake", "--build", str(out), "-j", "--target", "_core"])
    built = sorted(out.glob("_core*.so"))
    if len(built) != 1:
        raise BenchError(f"{out}: expected one _core*.so after the build, found {len(built)}")
    pkg = out / "pkg"
    shutil.rmtree(pkg, ignore_errors=True)
    (pkg / "tin_engine").mkdir(parents=True)
    for entry in sorted((tree / "src_python" / "tin_engine").iterdir()):
        if entry.name != "__pycache__" and not entry.name.startswith("_core."):
            (pkg / "tin_engine" / entry.name).symlink_to(entry)
    shutil.copy2(built[0], pkg / "tin_engine" / built[0].name)
    specs = [p.read_text() for p in sorted(out.glob("CMakeFiles/*/CMakeCXXCompiler.cmake"))]
    compiler = re.findall(r'CMAKE_CXX_COMPILER_(?:ID|VERSION) "([^"]*)"', "".join(specs[-1:]))
    return pkg, Build(
        no_build=False, type=_cache_value(cache, "CMAKE_BUILD_TYPE"),
        cxx_flags_release=_cache_value(cache, "CMAKE_CXX_FLAGS_RELEASE"),
        compiler=" ".join(compiler) or None,
        so_sha256=_sha256(built[0]),
    )  # fmt: skip


def _machine(runner: Runner) -> Machine:
    keys = ["machdep.cpu.brand_string", "hw.perflevel0.physicalcpu"]
    keys += ["hw.perflevel1.physicalcpu", "hw.memsize"]
    brand, *counts = (runner.run(["sysctl", "-n", k]).stdout.strip() for k in keys)
    p, e, mem = (int(c) if c.isdigit() else 0 for c in counts)
    return Machine(
        cpu_brand=brand, p_cores=p, e_cores=e, memory_bytes=mem,
        macos=runner.run(["sw_vers", "-productVersion"]).stdout.strip(),
        python=platform.python_version(), numpy=np.__version__,
    )  # fmt: skip


def _load(directory: Path) -> RunRecord:
    path = directory / "run.json"
    try:
        return RunRecord.model_validate_json(path.read_text())
    except (OSError, ValueError) as exc:  # pydantic's ValidationError is a ValueError
        raise BenchError(f"{path}: not a stored run: {str(exc).splitlines()[0]}") from exc


def _judge(
    runner: Runner, record: RunRecord, baseline: RunRecord | None, root: Path, tree: Path
) -> Verdict:
    if baseline is not None:
        return verdict(record, baseline)
    git = ["git", "-C", str(tree), "merge-base", "--is-ancestor"]
    found = find_baseline(root, record, lambda a, d: runner.run([*git, a, d]).returncode == 0)
    return verdict(record, found[1] if found else None)


def _child(runner: Runner, argv: list[str], run: str) -> tuple[Child, float]:
    """One child process: its parsed ``BENCH`` line and its wall time."""
    done = runner.run(argv)
    if done.returncode != 0:
        raise BenchError(f"child {run} exited {done.returncode}: {done.stderr[-2000:]}")
    try:
        return parse_child(done.stderr, run), done.wall_s
    except ChildError as exc:
        raise BenchError(str(exc)) from exc


def _measure(
    runner: Runner, head: list[str], domains: list[Domain], mesh: list[str],
    threads: list[int], repeats: int, mesh_dir: Path, tolerance: float, label: str,
) -> tuple[list[Sample], dict[str, MeshQuality], set[str]]:  # fmt: skip
    """Timing runs per domain, repeat, thread count (0 first), then that
    domain's ``--ascii`` quality run at t=0; and the hardening modes reported."""
    modes: set[str] = set()
    samples: list[Sample] = []
    qualities: dict[str, MeshQuality] = {}
    for d in domains:
        base = [*mesh, *(["--domain", str(d.path)] if d.path else [])]
        timing = ["--out", str(mesh_dir / "timing.vtk"), "--binary"]
        for r in range(repeats):
            for t in [0, *threads]:
                argv = [*head, "--threads", str(t), "--", *base, *timing]
                c, proc_s = _child(runner, argv, f"{d.name} t={t} r={r}")
                modes.add(c.hardening)
                samples.append(Sample(domain=d.name, threads=t, repeat=r, proc_s=proc_s,
                                      **{k: getattr(c, k) for k in CHILD_KEYS}))  # fmt: skip
        out = mesh_dir / f"{label}_{d.name}.vtk"
        argv = [*head, "--threads", "0", "--", *base, "--out", str(out), "--ascii"]
        c, _ = _child(runner, argv, f"{d.name} quality")
        modes.add(c.hardening)
        vtk = read_vtk_ascii(out)
        q = quality(vtk.points, vtk.triangles, vtk.edges, tolerance, c.max_error)
        qualities[d.name] = q.model_copy(update={"mesh_sha256": vtk.sha256})
    return samples, qualities, modes


app = typer.Typer(help="The acceptance benchmark (docs/benchmarks/bench-py.md).")

OutRoot = Annotated[Path, typer.Option("--out-root", help="Where evidence is stored and searched.")]
BaselineOpt = Annotated[
    Path | None, typer.Option("--baseline", help="A run directory to judge against.")
]


@app.command()
def run(
    label: Annotated[str, typer.Option("--label", help="Names the evidence directory.")],
    extra: Annotated[list[str] | None, typer.Argument(help="After --: extra mesh args.")] = None,
    tree: Annotated[Path, typer.Option("--tree", help="The source tree measured.")] = REPO,
    dem: Annotated[Path, typer.Option("--dem")] = DEFAULT_DEM,
    domain: Annotated[list[str] | None, typer.Option(help="Repeatable; `tile`: whole DEM.")] = None,
    tolerance: Annotated[float, typer.Option("--tolerance")] = 1.0,
    threads: Annotated[str, typer.Option(help="1,2,...")] = ",".join(map(str, range(1, 21))),
    repeats: Annotated[int, typer.Option("--repeats", min=1)] = 5,
    mesh_dir: Annotated[Path | None, typer.Option("--mesh-dir", help="Quality meshes.")] = None,
    no_build: Annotated[bool, typer.Option("--no-build", help="Installed package.")] = False,
    baseline: BaselineOpt = None,
    out_root: OutRoot = REPO / "docs/benchmarks",
    threshold: Annotated[float, typer.Option(help="Percent on the median.")] = 5.0,
    accept_quality: Annotated[bool, typer.Option(help="Waive angle and degree loss.")] = False,
    hardening: Annotated[Literal["on", "off"] | None, typer.Option(help="Default on.")] = None,
) -> None:
    """Build, measure, store the evidence, and judge it against a baseline."""
    runner = make_runner()
    tree, dem = tree.resolve(), dem.resolve()
    counts = [int(t) for t in threads.split(",") if t.strip()]
    names = domain or ["tile", str(DEFAULT_QUARTER)]
    try:  # refuse bad input before the build and any child
        for path in [dem, *(Path(n) for n in names if n != "tile")]:
            if not path.exists():
                raise BenchError(f"{path}: no such file")
        base = _load(baseline) if baseline is not None else None
        if no_build and hardening is not None:
            raise BenchError("--hardening with --no-build: nothing is built")
    except BenchError as exc:
        typer.echo(f"bench: {exc}")
        raise typer.Exit(3) from exc
    domains = [
        Domain(name="tile", path=None, sha256=None) if n == "tile"
        else Domain(name=Path(n).stem, path=str(Path(n).resolve()), sha256=_sha256(Path(n)))
        for n in names
    ]  # fmt: skip
    started = datetime.now().astimezone()
    meshes = mesh_dir or Path(tempfile.gettempdir()) / "rasputin-bench" / label
    meshes.mkdir(parents=True, exist_ok=True)
    try:
        if no_build:
            pkg, built = None, Build(no_build=True)
        else:
            pkg, built = build(runner, tree, hardening or "on")
        head = [sys.executable, str(BENCH), "_child", *(["--pkg", str(pkg)] if pkg else [])]
        mesh = ["mesh", "--dem", str(dem), "--tolerance", f"{tolerance:g}", *(extra or [])]
        commit = _checked(runner, ["git", "-C", str(tree), "rev-parse", "HEAD"]).stdout.strip()
        dirty = bool(runner.run(["git", "-C", str(tree), "status", "--porcelain"]).stdout.strip())
        blob = runner.run(["git", "-C", str(REPO), "hash-object", str(BENCH)]).stdout.strip()
        machine = _machine(runner)
        before = parse_pmset(runner.run(PMSET).stdout)
        samples, qualities, modes = _measure(
            runner, head, domains, mesh, counts, repeats, meshes, tolerance, label
        )
        if len(modes) != 1:
            raise BenchError(f"children reported different hardening modes: {sorted(modes)}")
        (mode,) = modes
        if hardening is not None and (hardening == "on") != (mode != "none"):
            raise BenchError(f"--hardening {hardening}, but the children's _core reports {mode}")
        power = combine_power(before, parse_pmset(runner.run(PMSET).stdout))
    except BenchError as exc:
        typer.echo(f"bench: {exc}")
        raise typer.Exit(3) from exc
    record = RunRecord(
        label=label, started=started, tree=Tree(commit=commit, dirty=dirty), bench_blob=blob,
        child_argv=[*head, "--threads", "<N>", "--", *mesh, "--out", "<path>", "--binary"],
        build=built, machine=machine, power=power,
        inputs=Inputs(dem=str(dem), dem_sha256=_sha256(dem), domains=domains,
                      tolerance=tolerance, extra_args=list(extra or [])),
        samples=samples, stats=median_stats(samples), quality=qualities,
        accept_quality=accept_quality, threshold_pct=threshold, hardening=mode,
    )  # fmt: skip
    v = _judge(runner, record, base, out_root, tree)
    record.verdict = v.lines
    directory = out_root / started.date().isoformat() / label
    write_evidence(record, directory)
    typer.echo(f"power: {power.state}; evidence: {directory}; meshes: {meshes}")
    typer.echo("\n".join(v.lines))
    raise typer.Exit(v.exit_code)


@app.command()
def compare(
    new_dir: Annotated[Path, typer.Argument(help="A stored run directory.")],
    baseline: BaselineOpt = None,
    out_root: OutRoot = REPO / "docs/benchmarks",
) -> None:
    """Re-judge stored evidence with the threshold and waiver it carries."""
    try:
        new, base = _load(new_dir), (_load(baseline) if baseline is not None else None)
    except BenchError as exc:
        typer.echo(f"bench: {exc}")
        raise typer.Exit(3) from exc
    v = _judge(make_runner(), new, base, out_root, REPO)
    typer.echo("\n".join(v.lines))
    raise typer.Exit(v.exit_code)


# ------------------------------------------------------------------ child
def child_main(argv: list[str]) -> int:
    """``_child [--pkg DIR] --threads N -- <rasputin argv>``: one sample."""
    split = argv.index("--")
    own, rasputin = argv[:split], argv[split + 1 :]
    threads = int(own[own.index("--threads") + 1])
    pkg = own[own.index("--pkg") + 1] if "--pkg" in own else None
    if pkg is not None:  # run.py's technique: the pkg tree, not the editable finder
        sys.meta_path[:] = [f for f in sys.meta_path if "ScikitBuild" not in type(f).__name__]
        sys.path.insert(0, pkg)
    import tin_engine
    import tin_engine._core
    import tin_engine.cli as cli

    for module in (tin_engine, tin_engine._core) if pkg is not None else ():
        if not str(module.__file__).startswith(pkg or ""):
            raise SystemExit(f"bench child: {module.__file__} is not under --pkg {pkg}")
    real = vars(cli)["refine"]  # cli imports refine from _core without re-exporting it
    seen: dict[str, Any] = {}

    def refine(*args: Any, **kwargs: Any) -> Any:  # scale.py's technique
        if threads != 0:
            kwargs["threads"] = threads
        t0 = time.perf_counter()
        out = real(*args, **kwargs)
        seen.update(refine_s=time.perf_counter() - t0, out=out)
        return out

    vars(cli)["refine"] = refine
    t0 = time.perf_counter()
    cli.app(args=rasputin, prog_name="rasputin", standalone_mode=False)
    app_s = time.perf_counter() - t0
    if "out" not in seen:
        print("bench child: refine was never called (no --tolerance?)", file=sys.stderr)
        return 1
    out = seen["out"]
    values = dict(refine_s=seen["refine_s"], app_s=app_s, max_error=out.max_error,
                  rounds=out.rounds, inserted=out.inserted, flips=out.flips,
                  hardening=getattr(tin_engine._core, "hardening", "none"))  # fmt: skip
    print("BENCH " + json.dumps(values), file=sys.stderr)
    return 0


if __name__ == "__main__":
    if sys.argv[1:2] == ["_child"]:
        sys.exit(child_main(sys.argv[2:]))
    app()
