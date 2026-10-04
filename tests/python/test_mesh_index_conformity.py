"""K4's check when the index is written: two pieces' records of one seam edge
must be equal bit for bit, or the run fails naming the seam and the pieces.

`docs/increments/23-basin-scale.md`, "Output" ("Conformity is checked when the
index is written") and K4. The end-to-end suites (`test_cli_mesh_pieces*.py`)
check conformity from the piece files (DC2); a correct run never reaches the
refusal, so this file reaches it directly.

PINNED HERE, where the design names `io/mesh_index.py` and "seam records" but
no surface (handback, choice 6):

    SeamRecord(piece=(j, i, k), edge=<index in the start triangulation's edge
               list>, xyz=<(K, 3) vertex sequence from the edge's lower end>)
    check_conformity(records) -> None
    ConformityError(ValueError) with `.edge` and `.pieces` (the two ids)

Bit for bit means bytes: a NaN z (a seam end with no height) equals a NaN z,
and 0.0 differs from -0.0. One record alone is a seam along the outline
(degeneracy policy) and passes.

HOW THIS FILE GOES RED: `tin_engine.io.mesh_index` does not exist; the `mi`
fixture fails every test on `ModuleNotFoundError`.
"""

from __future__ import annotations

import importlib
import json
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from pydantic import ValidationError

XYZ = (
    (500_000.0, 6_600_000.0, 12.5),
    (500_010.0, 6_600_000.0, 13.25),
    (500_030.0, 6_600_000.0, 9.0),
)


@pytest.fixture
def mi() -> ModuleType:
    return importlib.import_module("tin_engine.io.mesh_index")


def record(mi: ModuleType, piece: tuple[int, int, int], edge: int, xyz: Any = XYZ) -> Any:
    return mi.SeamRecord(piece=piece, edge=edge, xyz=tuple(tuple(map(float, p)) for p in xyz))


def changed(xyz: Any, row: int, col: int, value: float) -> tuple[tuple[float, ...], ...]:
    out = [list(p) for p in xyz]
    out[row][col] = value
    return tuple(tuple(p) for p in out)


class TestConforming:
    def test_two_equal_records_pass(self, mi: ModuleType) -> None:
        mi.check_conformity([record(mi, (0, 0, 0), 3), record(mi, (0, 1, 0), 3)])

    def test_many_edges_in_any_order_pass(self, mi: ModuleType) -> None:
        other = changed(XYZ, 1, 2, 99.0)
        records = [
            record(mi, (0, 1, 0), 3),
            record(mi, (1, 0, 0), 7, other),
            record(mi, (0, 0, 0), 3),
            record(mi, (0, 0, 0), 7, other),
        ]
        mi.check_conformity(records)
        mi.check_conformity(list(reversed(records)))

    def test_a_seam_along_the_outline_has_one_record(self, mi: ModuleType) -> None:
        mi.check_conformity([record(mi, (0, 0, 0), 5)])

    def test_nan_heights_compare_bit_for_bit(self, mi: ModuleType) -> None:
        xyz = changed(XYZ, 0, 2, float("nan"))
        mi.check_conformity([record(mi, (0, 0, 0), 2, xyz), record(mi, (1, 0, 0), 2, xyz)])

    def test_no_records(self, mi: ModuleType) -> None:
        mi.check_conformity([])


class TestRefused:
    @pytest.mark.parametrize(
        ("row", "col", "how"),
        [
            (1, 2, "last bit of z"),
            (1, 0, "last bit of x"),
            (2, 1, "last bit of y"),
            (0, 2, "minus zero"),
        ],
    )
    def test_one_value_differs(self, mi: ModuleType, row: int, col: int, how: str) -> None:
        base = changed(XYZ, 0, 2, 0.0) if how == "minus zero" else XYZ
        value = -0.0 if how == "minus zero" else float(np.nextafter(base[row][col], np.inf))
        other = changed(base, row, col, value)
        with pytest.raises(mi.ConformityError) as caught:
            mi.check_conformity([record(mi, (0, 0, 0), 4, base), record(mi, (0, 1, 0), 4, other)])
        assert isinstance(caught.value, ValueError)
        assert caught.value.edge == 4
        assert set(map(tuple, caught.value.pieces)) == {(0, 0, 0), (0, 1, 0)}
        assert "seam" in str(caught.value)

    def test_a_vertex_more_on_one_side(self, mi: ModuleType) -> None:
        longer = (*XYZ, (500_040.0, 6_600_000.0, 1.0))
        with pytest.raises(mi.ConformityError) as caught:
            mi.check_conformity([record(mi, (1, 1, 0), 9), record(mi, (1, 0, 0), 9, longer)])
        assert caught.value.edge == 9
        assert set(map(tuple, caught.value.pieces)) == {(1, 1, 0), (1, 0, 0)}

    def test_the_failing_edge_is_named_among_good_ones(self, mi: ModuleType) -> None:
        bad = changed(XYZ, 2, 2, 9.5)
        records = [
            record(mi, (0, 0, 0), 1),
            record(mi, (0, 1, 0), 1),
            record(mi, (0, 0, 0), 6),
            record(mi, (1, 0, 0), 6, bad),
        ]
        with pytest.raises(mi.ConformityError) as caught:
            mi.check_conformity(records)
        assert caught.value.edge == 6
        assert set(map(tuple, caught.value.pieces)) == {(0, 0, 0), (1, 0, 0)}


# ------------------------------------------------------------------ MeshIndex drift


def an_index() -> dict[str, Any]:
    """A valid index by hand, with one written piece and one all over NoData."""
    counts = {"triangles": 12, "vertices": 10, "on_frozen": 0}
    window = {"row0": -2, "col0": 0, "rows": 21, "cols": 20}
    return {
        "crs": "EPSG:25833",
        "tolerance": 1.0,
        "source": "refined from DEM nodes",
        "vocabulary": "0" * 64,
        "partition": {
            "pieces": 4,
            "memory_budget": 2**34,
            "bytes_per_node": 267,
            "nx": 2,
            "ny": 1,
            "dx": 20,
            "dy": 39,
            "cols": 39,
            "rows": 39,
        },
        "pieces": [
            {"id": [0, 0, 0], "file": "0-0-0.vtk", "sha256": "a" * 64, "counts": counts,
             "window": window},
            {"id": [0, 1, 0], "file": None, "sha256": None,
             "counts": {"triangles": 0, "vertices": 0, "on_frozen": 0}, "window": window},
        ],
    }  # fmt: skip


class TestMeshIndexDrift:
    """ "Three points from 23c-1's green (3ffae13)", point 2: piece ids are
    strict non-negative integers, and `file` and `sha256` are null together
    or set together."""

    def test_the_hand_built_index_validates(self, mi: ModuleType) -> None:
        index = mi.MeshIndex.model_validate(an_index())
        assert index.pieces[1].file is None and index.pieces[1].sha256 is None

    @pytest.mark.parametrize(
        "piece_id",
        [["1", "2", "3"], [0, "1", 0], [0, -1, 0], [-1, 0, 0], [0, 0, -1], [0.0, 1.0, 0.0]],
        ids=["text", "one-text", "negative-i", "negative-j", "negative-k", "floats"],
    )
    def test_a_piece_id_not_three_non_negative_integers(
        self, mi: ModuleType, piece_id: list[Any]
    ) -> None:
        doc = an_index()
        doc["pieces"][0]["id"] = piece_id
        with pytest.raises(ValidationError):
            mi.MeshIndex.model_validate(doc)
        with pytest.raises(ValidationError):
            mi.MeshIndex.model_validate_json(json.dumps(doc))

    @pytest.mark.parametrize(
        ("file", "sha256"), [("0-0-0.vtk", None), (None, "a" * 64)], ids=["no-sha", "no-file"]
    )
    def test_file_and_sha256_null_only_together(
        self, mi: ModuleType, file: str | None, sha256: str | None
    ) -> None:
        doc = an_index()
        doc["pieces"][0]["file"], doc["pieces"][0]["sha256"] = file, sha256
        with pytest.raises(ValidationError):
            mi.MeshIndex.model_validate(doc)
        with pytest.raises(ValidationError):
            mi.MeshIndex.model_validate_json(json.dumps(doc))
