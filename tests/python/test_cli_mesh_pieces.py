"""`rasputin mesh ... --pieces P --memory-budget B`: increment 23c end to end.

`docs/increments/23-basin-scale.md`, "Choosing the cuts", "The seam
protocol", "Output", K1-K6 and K9, and "Tests @tester can write red", 23c:
DC1 (K1), DC2 (conformity from the files), DC3 (a valid union), DC4
(constrained Delaunay with the seams) and DC5 (the tolerance over the union),
on a synthetic 41 x 41-node DEM at 10 m and an off-node square domain
(`pieces_fixtures`). DC2-DC5 are checked together on every path a run takes
(tester.md section 3D: both oracles on every path), with feet and the start
quality each on and off, and at two tolerances; each oracle has a control
that fails it on a wrong mesh.

PINNED HERE, where the design is silent (handback, choices 2-5 and 7-9):

- `--pieces` and `--memory-budget` are `mesh` options; the budget is a whole
  number of bytes (default 2**34, the design's "16 GB, read as 16 GiB");
  both need `--dem` and `--tolerance`; `--pieces 0` and a budget under one
  node's `b(T)` are usage errors naming the option.
- A cut run writes `<out name>.pieces/` beside `--out`, holding
  `index.json` (`MeshIndex`, `tin_engine.io.mesh_index`) and one file per
  piece in `--out`'s format. Whether `--out` itself is written is 23d's
  (the stitcher) and is not asserted.
- An uncut run (one piece) writes no `.pieces` directory and the same bytes
  as today, with the options given or not: so no field records the budget
  on that path (K1 over "recorded in the file", K5).
- K9 through spies: `node` and `triangulate` run once in a cut run.

HOW THIS FILE GOES RED: `--pieces` and `--memory-budget` are not options of
`mesh`, so every cut run exits 2 with "No such option"; `cut()` asserts the
exit code and the directory. The uncut controls (today's bytes) pass.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from shapely.geometry import Polygon

import pieces_fixtures as pf
from pieces_fixtures import Cut, Union
from test_cli_mesh_dem import USAGE
from tin_engine import _core

DOMAIN = Polygon(pf.SQUARE)
GIB = 1 << 30


@pytest.fixture(scope="module")
def scene(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, np.ndarray]:
    where = tmp_path_factory.mktemp("scene")
    array = pf.hills()
    return (
        pf.write_dem(where / "dem.tif", array),
        pf.write_domain(where / "square.geojson", pf.SQUARE),
        array,
    )


def base(scene: tuple[Path, Path, np.ndarray], tolerance: str = "1") -> tuple[str, ...]:
    dem, domain, _ = scene
    return ("--dem", str(dem), "--domain", str(domain), "--tolerance", tolerance)


def all_problems(done: Cut, array: np.ndarray, tolerance: float) -> list[str]:
    """DC2-DC5 over one cut run, as one list of what is wrong."""
    problems = pf.conformity_problems(done.pieces)
    problems += [f"seam edge off the lattice: {e}" for e in pf.off_lattice_seams(done.pieces)]
    u = pf.union(done.pieces)
    problems += [f"two z at {k}" for k in u.z_conflicts]
    problems += pf.validity_problems(u, DOMAIN)
    problems += [f"triangles {p} overlap" for p in pf.overlapping_pairs(u)]
    violations = pf.delaunay_violations(u.points[:, :2], u.triangles, u.constrained)
    problems += [f"not constrained Delaunay: {v}" for v in violations]
    found = pf.tolerance_check(u, array, DOMAIN, tolerance)
    if found.over or found.outside_mesh:
        problems.append(f"tolerance: {found}")
    return problems


# ------------------------------------------------------------------ options


class TestOptions:
    @pytest.mark.parametrize("value", ["0", "-2"])
    def test_pieces_under_one_is_refused(self, tmp_path: Path, scene: Any, value: str) -> None:
        done = pf.run(tmp_path, *base(scene), "--pieces", value)
        assert done.code == USAGE, done.output
        assert "No such option" not in done.output
        assert "--pieces" in done.output
        assert done.dir is None

    @pytest.mark.parametrize("value", ["0", "-1", "266"])
    def test_a_budget_under_one_node_is_refused(
        self, tmp_path: Path, scene: Any, value: str
    ) -> None:
        """b(1 m) is 267 bytes: a budget below it cannot hold a cell of one node."""
        done = pf.run(tmp_path, *base(scene), "--memory-budget", value)
        assert done.code == USAGE, done.output
        assert "No such option" not in done.output
        assert "--memory-budget" in done.output
        assert done.dir is None

    @pytest.mark.parametrize("option", [("--pieces", "4"), ("--memory-budget", "4096")])
    def test_needs_a_tolerance(self, tmp_path: Path, scene: Any, option: tuple[str, str]) -> None:
        dem, _, _ = scene
        done = pf.run(tmp_path, "--dem", str(dem), *option)
        assert done.code == USAGE, done.output
        assert "No such option" not in done.output
        assert "--tolerance" in done.output
        assert done.dir is None

    def test_needs_a_dem(self, tmp_path: Path) -> None:
        done = pf.run(tmp_path, "square", "--flat", "--pieces", "4")
        assert done.code == USAGE, done.output
        assert "No such option" not in done.output
        assert "--dem" in done.output

    def test_the_help_names_both(self) -> None:
        from typer.testing import CliRunner

        from tin_engine.cli import app

        text = CliRunner().invoke(app, ["mesh", "--help"], env={"COLUMNS": "200"}).output
        assert "--pieces" in text
        assert "--memory-budget" in text


# ------------------------------------------------------------------ DC1, K1


class TestUncutIsTodaysMesh:
    """DC1: a domain under the budget, and with a budget at its estimate,
    writes the same bytes as today's path, and no pieces."""

    def test_the_options_at_one_piece_change_no_byte(self, tmp_path: Path, scene: Any) -> None:
        today = pf.run(tmp_path / "today", *base(scene), "--binary")
        assert today.code == 0, today.output
        same = pf.run(
            tmp_path / "same",
            *base(scene),
            "--binary",
            *("--pieces", "1", "--memory-budget", str(1 << 40)),
        )
        assert same.code == 0, same.output
        assert same.out.read_bytes() == today.out.read_bytes()
        assert (today.dir, same.dir) == (None, None)

    def test_the_budget_boundary_through_the_cli(self, tmp_path: Path, scene: Any) -> None:
        """One piece at `N b = B`, cut at `N b = B + 1`: N and b from a cut run's
        index, so the window is the run's own, not re-derived here."""
        probe = pf.cut(tmp_path / "probe", *base(scene), "--pieces", "4")
        part = probe.index["partition"]
        budget = int(part["cols"]) * int(part["rows"]) * int(part["bytes_per_node"])
        assert int(part["bytes_per_node"]) == 267  # b(1 m)

        today = pf.run(tmp_path / "today", *base(scene))
        whole = pf.run(tmp_path / "whole", *base(scene), "--memory-budget", str(budget))
        assert whole.code == 0, whole.output
        assert whole.dir is None
        assert whole.out.read_bytes() == today.out.read_bytes()

        split = pf.cut(tmp_path / "split", *base(scene), "--memory-budget", str(budget - 1))
        assert int(split.index["partition"]["nx"]) * int(split.index["partition"]["ny"]) > 1


# ------------------------------------------------------------------ the output


class TestTheOutput:
    @pytest.fixture(scope="class")
    def done(self, tmp_path_factory: pytest.TempPathFactory, scene: Any) -> Cut:
        return pf.cut(tmp_path_factory.mktemp("out"), *base(scene), "--pieces", "4")

    def test_the_index_records_the_run(self, done: Cut, tmp_path: Path, scene: Any) -> None:
        today = pf.run(tmp_path, *base(scene))
        assert today.code == 0, today.output
        crs = pf.read_vtk(today.out.read_bytes()).field_data["crs"].values[0]
        index = done.index
        assert index["crs"] == crs
        assert index["tolerance"] == 1.0
        part = index["partition"]
        assert part["pieces"] == 4
        assert part["memory_budget"] == 16 * GIB
        assert part["bytes_per_node"] == 267
        assert (part["nx"], part["ny"]) == (2, 2)
        for key in ("dx", "dy", "cols", "rows"):
            assert isinstance(part[key], int) and part[key] > 0, key
        assert part["dx"] * (part["nx"] - 1) < part["cols"] <= part["dx"] * part["nx"]
        assert part["dy"] * (part["ny"] - 1) < part["rows"] <= part["dy"] * part["ny"]

    def test_one_piece_per_cell_of_a_convex_domain(self, done: Cut) -> None:
        ids = sorted(p.id for p in done.pieces)
        assert ids == [(0, 0, 0), (0, 1, 0), (1, 0, 0), (1, 1, 0)]

    def test_every_listed_file_exists_with_its_sha256(self, done: Cut) -> None:
        assert done.dir is not None
        for entry in done.index["pieces"]:
            path = done.dir / entry["file"]
            assert path.parent == done.dir, entry
            assert path.suffix == ".vtk", entry
            assert entry["sha256"] == pf.sha256(path), entry

    def test_the_index_is_a_mesh_index(self, done: Cut) -> None:
        from tin_engine.io.mesh_index import MeshIndex

        assert done.dir is not None
        text = (done.dir / "index.json").read_text()
        model = MeshIndex.model_validate_json(text)
        assert json.loads(model.model_dump_json()) == json.loads(text)

    def test_each_piece_file_is_a_mesh_in_the_runs_crs(self, done: Cut) -> None:
        for piece in done.pieces:
            assert piece.vtk.field_data["crs"].values[0] == done.index["crs"]
            assert len(piece.triangles) > 0, piece.id
            assert piece.seam > 0

    def test_the_seam_bit_is_not_a_feature_bit(self, done: Cut) -> None:
        names = list(done.pieces[0].vtk.field_data["feature_names"].values)
        bits = [int(b) for b in np.asarray(done.pieces[0].vtk.field_data["feature_bits"].values)]
        assert names.count("seam") == 1
        assert len(set(bits)) == len(bits)

    def test_a_ply_out_writes_ply_pieces(self, tmp_path: Path, scene: Any) -> None:
        done = pf.run(tmp_path, *base(scene), "--pieces", "4", out_name="x.ply")
        assert done.code == 0, done.output
        assert done.dir is not None
        index = json.loads((done.dir / "index.json").read_text())
        files = [done.dir / e["file"] for e in index["pieces"]]
        assert len(files) == 4
        for path in files:
            assert path.suffix == ".ply" and path.read_bytes().startswith(b"ply\n"), path


# ------------------------------------------------------------------ DC2-DC5, section 3D


PATHS = [
    pytest.param((), 1.0, id="defaults"),
    pytest.param(("--no-constraint-feet",), 1.0, id="feet-off"),
    pytest.param(("--start-min-angle", "0"), 1.0, id="quality-off"),
    pytest.param(("--no-constraint-feet", "--start-min-angle", "0"), 1.0, id="both-off"),
    pytest.param((), 0.25, id="tolerance-0.25"),
    pytest.param((), 4.0, id="tolerance-4"),
]


class TestTheUnion:
    @pytest.mark.parametrize(("extra", "tolerance"), PATHS)
    @pytest.mark.parametrize("pieces", ["4", "9"])
    def test_both_oracles_and_conformity_on_every_path(
        self, tmp_path: Path, scene: Any, extra: tuple[str, ...], tolerance: float, pieces: str
    ) -> None:
        _, _, array = scene
        done = pf.cut(tmp_path, *base(scene, repr(tolerance)), "--pieces", pieces, *extra)
        assert len(done.pieces) == int(pieces)
        assert all_problems(done, array, tolerance) == []

    def test_the_union_has_the_uncut_area(self, tmp_path: Path, scene: Any) -> None:
        done = pf.cut(tmp_path / "cut", *base(scene), "--pieces", "4")
        today = pf.run(tmp_path / "today", *base(scene))
        whole = pf.read_vtk(today.out.read_bytes())
        u = pf.union(done.pieces)

        def area(points: np.ndarray, tris: np.ndarray) -> float:
            a, b, c = points[tris[:, 0], :2], points[tris[:, 1], :2], points[tris[:, 2], :2]
            twice = (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (
                c[:, 0] - a[:, 0]
            )
            return float(twice.sum() / 2)

        expected = area(np.asarray(whole.points), pf.polygons_as_array(whole))
        assert area(u.points, u.triangles) == pytest.approx(expected, rel=1e-12)


class TestTheControls:
    """Each oracle fails on a mesh that is wrong in its way."""

    @pytest.fixture(scope="class")
    def done(self, tmp_path_factory: pytest.TempPathFactory, scene: Any) -> Cut:
        return pf.cut(tmp_path_factory.mktemp("controls"), *base(scene), "--pieces", "4")

    def test_dc4_without_the_seams_as_constraints_finds_violations(self, done: Cut) -> None:
        """What the union is: constrained Delaunay with the seams, not without."""
        u = pf.union(done.pieces)
        with_seams = pf.delaunay_violations(u.points[:, :2], u.triangles, u.constrained)
        assert with_seams == []
        assert u.seam_only, "no seam edge in the union"
        without = pf.delaunay_violations(u.points[:, :2], u.triangles, u.constrained - u.seam_only)
        assert without != []
        assert {edge for edge, _ in without} <= u.seam_only

    def test_dc5_with_z_shifted_by_twice_the_tolerance_fails(self, done: Cut, scene: Any) -> None:
        _, _, array = scene
        u = pf.union(done.pieces)
        assert pf.tolerance_check(u, array, DOMAIN, 1.0).over == 0
        shifted = Union(
            points=u.points + np.array([0.0, 0.0, 2.0]),
            triangles=u.triangles,
            constrained=u.constrained,
            seam_only=u.seam_only,
            z_conflicts=u.z_conflicts,
        )
        assert pf.tolerance_check(shifted, array, DOMAIN, 1.0).over > 0

    def test_dc2_a_seam_vertex_moved_in_one_file_is_caught(self, done: Cut) -> None:
        assert pf.conformity_problems(done.pieces) == []
        first = done.pieces[0]
        tampered = pf.Piece(**{**first.__dict__, "points": first.points.copy()})
        (a, _), _m = next(
            (e, m) for e, m in zip(first.edges, first.masks, strict=True) if m & first.seam
        )
        tampered.points[a, 2] = np.nextafter(tampered.points[a, 2], np.inf)
        assert pf.conformity_problems([tampered, *done.pieces[1:]]) != []

    def test_dc3_a_dropped_triangle_is_caught(self, done: Cut) -> None:
        u = pf.union(done.pieces)
        assert pf.validity_problems(u, DOMAIN) == []
        short = Union(
            points=u.points,
            triangles=u.triangles[1:],
            constrained=u.constrained,
            seam_only=u.seam_only,
            z_conflicts=[],
        )
        assert pf.validity_problems(short, DOMAIN) != []


# ------------------------------------------------------------------ K9


def test_k9_noding_and_the_start_triangulation_run_once(
    tmp_path: Path, scene: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One global vector step: `node` and `triangulate` once in a cut run,
    wherever the pieces' code imported them from."""
    calls = {"node": 0, "triangulate": 0}
    for name in calls:
        original = getattr(_core, name)

        def spy(*args: Any, _name: str = name, _f: Any = original, **kwargs: Any) -> Any:
            calls[_name] += 1
            return _f(*args, **kwargs)

        monkeypatch.setattr(_core, name, spy)
        for module in [m for k, m in sys.modules.items() if k.startswith("tin_engine")]:
            if getattr(module, name, None) is original:
                monkeypatch.setattr(module, name, spy)
    pf.cut(tmp_path, *base(scene), "--pieces", "4")
    assert calls == {"node": 1, "triangulate": 1}
