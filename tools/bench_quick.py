"""The 15-minute quick check: is ``rasputin mesh`` slower than on master?

Design: ``docs/increments/perf-quick-check.md``. ``run`` builds the tree's
``_core`` in Release (``bench.build``), times the fixed cases in
``QUICK/cases.toml`` through ``bench.py _child``, and judges each case's total
and its ``--stats`` phases against ``QUICK/baseline-<power>.json``. Every
subprocess gets the time left before ``--budget`` less 10 s; a case whose
baseline does not fit in what is left is skipped. Exit 0 no change or faster,
1 slower or broken, 2 no baseline, 3 bad input, 4 out of time.
"""

from __future__ import annotations

import hashlib
import os
import statistics
import sys
import tempfile
import time
import tomllib
from pathlib import Path
from string import Template
from typing import Annotated

import bench
import typer
from bench import BenchError, Machine, Verdict
from pydantic import BaseModel

REPO = Path(__file__).resolve().parents[1]
QUICK = REPO / "docs/benchmarks/quick"
MAX_BUDGET = 900.0
#: Files under this many bytes are fingerprinted by content as well.
READ_UNDER = 100_000_000
make_runner = bench.make_runner


# ------------------------------------------------------------------ models
class Spread(BaseModel):
    median: float
    min: float
    max: float


class CaseResult(BaseModel):
    tolerance: float
    max_error: float
    mesh_sha256: str
    #: ``total`` is the clock's total; every other key a ``--stats`` phase.
    measures: dict[str, Spread]


class QuickRecord(BaseModel):
    commit: str
    machine: Machine
    power: str
    #: Each input as written in cases.toml, to its fingerprint.
    inputs: dict[str, str]
    cases: dict[str, CaseResult]


class Case(BaseModel):
    name: str
    threads: list[int]
    runs: int
    warmup: int
    args: list[str]
    inputs: list[str]


class RuledHotspot(BaseModel):
    case: str
    phase: str
    share: float
    date: str
    ruling: str


def load_cases(path: Path) -> list[Case]:
    return [Case.model_validate(c) for c in tomllib.loads(path.read_text())["case"]]


def load_hotspots(path: Path) -> list[RuledHotspot]:
    entries = tomllib.loads(path.read_text()).get("hotspot", []) if path.exists() else []
    return [RuledHotspot.model_validate(h) for h in entries]


def fingerprint(path: Path) -> str:
    """A file: size, mtime and, under 100 MB, content. A directory: its sorted
    relative file names, sizes and mtimes, recursive; no content read."""
    if path.is_dir():
        rows = sorted(
            f"{p.relative_to(path).as_posix()}\t{p.stat().st_size}\t{p.stat().st_mtime_ns}\n"
            for p in path.rglob("*")
            if p.is_file()
        )
        return "dir:" + hashlib.sha256("".join(rows).encode()).hexdigest()
    stat = path.stat()
    content = hashlib.sha256(path.read_bytes()).hexdigest() if stat.st_size < READ_UNDER else ""
    return f"file:{stat.st_size}:{stat.st_mtime_ns}:{content}"


# ----------------------------------------------------------------- judging
def _mismatch(new: QuickRecord, base: QuickRecord) -> str | None:
    """None, or the first field that forbids comparing ``new`` with ``base``."""
    if new.power not in ("ac", "battery") or new.power != base.power:
        return f"power: {new.power} vs {base.power}"
    old_machine = base.machine.model_dump()
    for key, value in new.machine.model_dump().items():
        if value != old_machine[key]:
            return f"machine {key}: {value} vs {old_machine[key]}"
    for key in sorted(new.inputs.keys() | base.inputs.keys()):
        if new.inputs.get(key) != base.inputs.get(key):
            return f"input {key} changed"
    return None


def verdict(new: QuickRecord, base: QuickRecord | None) -> Verdict:
    reason = "no stored baseline" if base is None else _mismatch(new, base)
    if reason is not None or base is None:
        return Verdict(f"NO BASELINE: {reason}", [], 2)
    lines: list[str] = []
    for name, case in new.cases.items():
        old = base.cases.get(name)
        if old is None:
            lines.append(f"NOT JUDGED: {name} is not in the baseline")
            continue
        if case.max_error > case.tolerance:
            lines.append(f"BROKEN: {name} max_error {case.max_error} over {case.tolerance}")
        if case.mesh_sha256 != old.mesh_sha256:
            lines.append(f"MESH CHANGED: {name} mesh sha256 {old.mesh_sha256[:12]} -> "
                         f"{case.mesh_sha256[:12]} (reported, not judged)")  # fmt: skip
        total = old.measures["total"].median
        for measure, b in old.measures.items():
            m = case.measures.get(measure)
            if m is None or b.median <= 0 or (measure != "total" and b.median * 100 < total * 5):
                continue
            band = max(5.0, (b.max - b.min) / b.median * 100)
            pct = (m.median / b.median - 1) * 100
            found = f"{name} {measure} {b.median:.2f} -> {m.median:.2f} s ({pct:+.1f} %, band "
            if m.median * 100 > b.median * (100 + band) and m.median - b.median > 0.05:
                lines.append(f"SLOWER: {found}{band:.1f} %)")
            elif m.median * 100 < b.median * (100 - band) and b.median - m.median > 0.05:
                lines.append(f"FASTER: {found}{band:.1f} %)")
    for status, code in (("BROKEN", 1), ("SLOWER", 1), ("FASTER", 0)):
        if any(line.startswith(f"{status}: ") for line in lines):
            return Verdict(status, lines, code)
    return Verdict("NO CHANGE", lines, 0)


def hotspots(record: QuickRecord, ruled: list[RuledHotspot]) -> list[str]:
    """A phase at 40 % or more of its case's total, unless ruled at a share
    it has not outgrown by 10 points."""
    lines: list[str] = []
    for name, case in record.cases.items():
        total = case.measures["total"].median
        for phase, spread in case.measures.items():
            if phase == "total" or total <= 0 or spread.median * 100 < total * 40:
                continue
            share = spread.median / total * 100
            rule = next((r for r in ruled if (r.case, r.phase) == (name, phase)), None)
            if rule is None:
                why = "not in hotspots.toml"
            elif share >= rule.share + 10:
                why = f"ruled at {rule.share:.0f} %, grown 10 points or more"
            else:
                continue
            lines.append(f"HOTSPOT: {name} {phase} {share:.0f} % of the run ({why})")
    return lines


# --------------------------------------------------------------------- run
def _spread(values: list[float]) -> Spread:
    return Spread(median=statistics.median(values), min=min(values), max=max(values))


def _mesh_sha256(path: Path) -> str:
    """From ``POINTS`` on, as bench.py hashes, so header fields do not count."""
    data = path.read_bytes() if path.exists() else b""
    return hashlib.sha256(data[data.find(b"\nPOINTS ") + 1 :]).hexdigest()


app = typer.Typer(help="The 15-minute quick check (docs/increments/perf-quick-check.md).")


@app.callback()
def main() -> None:
    """Subcommands: run."""


@app.command()
def run(
    tree: Annotated[Path, typer.Option("--tree", help="The source tree measured.")] = REPO,
    budget: Annotated[float, typer.Option("--budget", help="Seconds, at most 900.")] = MAX_BUDGET,
    save: Annotated[bool, typer.Option("--save-baseline", help="Store this run.")] = False,
) -> None:
    """Build, time the cases, and judge them against the stored baseline."""
    if not 0 < budget <= MAX_BUDGET:
        typer.echo(f"bench_quick: --budget {budget:g}: must be over 0 and at most 900 s")
        raise typer.Exit(3)
    deadline = time.monotonic() + budget
    tree = tree.resolve()
    data = os.environ.get("RASPUTIN_DATA") or str(REPO.parent / "rasputin_data")

    def expand(text: str) -> str:
        return Template(text).safe_substitute(RASPUTIN_DATA=data, TREE=str(tree))

    runner = make_runner()
    try:
        cases = load_cases(QUICK / "cases.toml")
        for path in (Path(expand(i)) for c in cases for i in c.inputs):
            if not path.exists():
                raise BenchError(f"{path}: no such input")
        inputs = {i: fingerprint(Path(expand(i))) for c in cases for i in c.inputs}
        before = bench.parse_pmset(runner.run(bench.PMSET).stdout)
        stored = QUICK / f"baseline-{before.state}.json"
        base = QuickRecord.model_validate_json(stored.read_text()) if stored.exists() else None
        pkg, _ = bench.build(runner, tree, deadline=deadline)
        commit = bench._checked(runner, ["git", "-C", str(tree), "rev-parse", "HEAD"])
        machine = bench._machine(runner)
    except (BenchError, OSError, ValueError) as exc:  # ValidationError is a ValueError
        late = time.monotonic() > deadline - 10
        typer.echo(f"{'OUT OF TIME: the build' if late else 'bench_quick:'} {exc}")
        raise typer.Exit(4 if late else 3) from exc
    head = [sys.executable, str(bench.BENCH), "_child", "--pkg", str(pkg)]
    out = Path(tempfile.mkdtemp(prefix="rasputin-quick-")) / "mesh.vtk"
    plan = [(c, t, c.name if t == 0 else f"{c.name} t={t}") for c in cases for t in c.threads]
    results: dict[str, CaseResult] = {}
    missed: list[str] = []
    for k, (case, threads, name) in enumerate(plan):
        old = base.cases.get(name) if base else None
        if old and old.measures["total"].median * case.runs > deadline - time.monotonic():
            missed.append(name)
            continue
        argv = [*head, "--threads", str(threads), "--", *map(expand, case.args)]
        children: list[bench.Child] = []
        for _ in range(case.warmup + case.runs):
            done = runner.run(
                [*argv, "--out", str(out), "--binary"], None, bench.time_left(deadline)
            )
            if done.timed_out:
                break
            try:
                if done.returncode != 0:
                    raise BenchError(f"exited {done.returncode}: {done.stderr[-2000:]}")
                children.append(bench.parse_child(done.stderr, name))
            except (BenchError, ValueError) as exc:
                typer.echo(f"bench_quick: child {name}: {exc}")
                raise typer.Exit(3) from exc
        else:
            timed = children[case.warmup :]
            measures = {"total": _spread([c.total_s or 0.0 for c in timed])}
            for phase in dict.fromkeys(p for c in timed for p, _ in c.phases):
                measures[phase] = _spread([s for c in timed for p, s in c.phases if p == phase])
            tolerance = float(case.args[case.args.index("--tolerance") + 1])
            results[name] = CaseResult(
                tolerance=tolerance,
                max_error=children[-1].max_error,
                mesh_sha256=_mesh_sha256(out),
                measures=measures,
            )
            continue
        missed += [n for _, _, n in plan[k:]]
        break
    power = bench.combine_power(before, bench.parse_pmset(runner.run(bench.PMSET).stdout))
    record = QuickRecord(commit=commit.stdout.strip(), machine=machine, power=power.state,
                         inputs=inputs, cases=results)  # fmt: skip
    v = verdict(record, base)
    lines = [*v.lines, *hotspots(record, load_hotspots(QUICK / "hotspots.toml"))]
    status, code = (f"OUT OF TIME: {', '.join(missed)}", 4) if missed else (v.status, v.exit_code)
    if save and not missed and power.state in ("ac", "battery"):
        path = QUICK / f"baseline-{power.state}.json"
        path.write_text(record.model_dump_json(indent=1) + "\n")
        lines.append(f"baseline saved: {path}")
    typer.echo("\n".join([*lines, status]))
    raise typer.Exit(code)


if __name__ == "__main__":
    app()
