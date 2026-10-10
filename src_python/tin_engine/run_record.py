"""The record of one ``rasputin mesh`` run, and every word it is printed in.

``docs/increments/25-plain-output.md`` D1 to D7. Pure: no numpy, no ``_core``,
no typer, so the wording is tested without the extension. ``cli.py`` copies
numbers and strings out of what ran into a builder, and prints the one
:class:`RunRecord` four ways: the mesh file's fields (:func:`file_fields`),
the ``--stats`` rows (:func:`stats_rows`), one stderr summary
(:func:`summary`) and the ``--record`` JSON (:func:`as_json`).
"""

from __future__ import annotations

import json
from dataclasses import dataclass

#: A builder argument: a value, or None for "not produced on this path".
Value = str | int | float | None

#: Every entry, in record order, with its ``--stats`` wording ("The fields,
#: for Ola"). The Inputs entries come first (``INPUTS``), then the mesh-file
#: fields in D2's order, then the counts and checks.
WORDING = {
    "dem_grid": "DEM grid size and spacing",
    "dem_tiles": "The tiles or downloaded blocks used",
    "dem_seams": "Overlapping tiles that disagree, and by how much",
    "dem_vertical_unit": "Unit of z",
    "dem_crs": "The DEM's own coordinate system",
    "dem_transform": "Conversion from the DEM's coordinate system",
    "resampled_grid": "The grid the DEM was interpolated onto",
    "domain": "The domain file and shape",
    "domain_crs": "The domain's coordinate system",
    "domain_transform": "Conversion from the domain's coordinate system",
    "features": "Each features file: layer, class map, features, lines, vertices",
    "features_crs": "Each features file's coordinate system",
    "features_transform": "Conversion from each features file's coordinate system",
    "features_repair_m": "Land-cover borders closer than this made one, narrower gaps filled",
    "features_merge_same_class": "Borders between land-cover polygons of one class dropped",
    "features_tolerance_m": "Land-cover borders simplified by this much (0 = off)",
    "features_outline_snap_m": "Land-cover borders this close to the outline moved onto it",
    "start_mesh": "What refinement started from",
    "start_min_angle_deg": "Starting mesh improved to this smallest angle (0 = off)",
    "start_quality_gain_deg": (
        "Starting mesh: a point is added only if it raises the smallest angle around it "
        "by at least this many degrees (negative = always added)"
    ),
    "snap_to_lines": "Points very close to a line were moved onto it",
    "tolerance_near_m": "Tolerance on the tolerance lines",
    "tolerance_ramp_m": "Distance from the lines over which the tolerance rises to its far value",
    "tolerance_lines": "The tolerance lines file, and its segments after simplifying",
    "crs": "Coordinate system",
    "tolerance_m": "Tolerance",
    "max_error_m": "Largest height error",
    "tolerance_slope_m": "Tolerance on steep ground",
    "tolerance_slope_deg": (
        "Slope where the tolerance starts to tighten, and where it reaches the steep value"
    ),
    "slope_nodes_tightened": (
        "DEM nodes held to less than the general tolerance because of their slope"
    ),
    "max_error_slope_share_of_tolerance": (
        "Largest error as a share of the tolerance its slope allows"
    ),
    "dem_source": "DEM",
    "dem_credit": "Credit",
    "licence_note": "Licence",
    "cite": "Please cite",
    "nodata_vertices_removed": "Vertices removed on NoData",
    "heights": "Heights",
    "features_notice": "Credit for the features data",
    "land_cover_codes": "What land_cover_code holds",
    "land_cover_vertices": "Land-cover vertices before and after clean-up",
    "land_cover_area_moved_m2": (
        "Land-cover area inside the outline that the outline rule gave to another polygon, m2"
    ),
    "resampled_grid_max_error_m": "Largest error against the resampled grid",
    "max_error_near_lines_m": "Largest height error where the tolerance is the lines' own",
    "dem_nodes_checked": "Nodes of the original DEM compared with the mesh",
    "dem_check_points_inserted": "Points that comparison added",
    "dem_check_rounds": "Passes of that comparison",
    "dem_nodes_at_vertices": "DEM nodes on a vertex (to rounding), not compared",
    "dem_nodes_at_vertices_max_error_m": "Their largest difference",
    "line_points_checked": "Points along lines (grid-line crossings and halfway between) compared",
    "line_max_error_m": "Their largest difference from the mesh",
    "line_points_on_nodata": "Points along lines on NoData cells, not compared",
    "line_points_refused": "Points along lines that could not be added",
    "line_points_refused_max_error_m": "Their largest difference from the mesh",
    "line_points_inserted": "Points along lines added",
    "line_check_dem_nodes_inserted": "DEM nodes the line check added",
    "line_points_duplicate": "Points along lines dropped as duplicates",
    "dem_nodes_outside_mesh": "Self-check: DEM nodes with data left outside the mesh",
    "refinement_rounds": "Refinement passes",
    "points_inserted": "Points added",
    "points_inserted_on_nodata": "Of them, into triangles with a NoData corner",
    "edge_flips": "Edge swaps",
    "start_quality_points_inserted": "Points added to improve the starting mesh",
    "start_quality_points_skipped": "Tries skipped while improving it",
    "start_quality_points_without_gain": "Tries that would not have improved the angles",
    "start_quality_points_snapped_to_lines": (
        "Points moved onto lines while improving the starting mesh"
    ),
    "start_quality_lines_split": "Land-cover and outline lines split to improve angles",
    "points_snapped_to_lines": "Points refinement moved onto lines",
    "snaps_refused": "Moves onto lines refused",
    "final_check_points_snapped_to_lines": (
        "Points the final check against the DEM moved onto lines"
    ),
    "final_check_snapped_points_added_anyway": "Of them, also added where they were",
    "start_vertices_between_dem_nodes": "Starting-mesh vertices not on a DEM node",
}
_NAMES = list(WORDING)
INPUTS = frozenset(_NAMES[: _NAMES.index("crs")])
#: D2's mesh-file fields, and the CORINE notice kept by Ola's ruling.
#: ``land_cover_codes`` is in the file too, but the writers put it there.
IN_FILE = frozenset(_NAMES[_NAMES.index("crs") : _NAMES.index("features_notice") + 1])
FLAT_HEIGHTS = "none: every z is 0 (--flat)"
SEAMS_AGREE = "none: the tiles agree where they overlap"


@dataclass(frozen=True, slots=True)
class Entry:
    name: str  # plain snake_case words, the unit as a suffix (_m, _deg)
    wording: str  # the label --stats prints
    value: str  # ASCII, already formatted
    number: float | int | None  # the bare number, for --record's JSON; None for text
    in_file: bool  # True only for the mesh-file fields of D2
    in_inputs: bool = False  # --stats "Inputs", else "Result"


@dataclass(frozen=True, slots=True)
class RunRecord:
    entries: tuple[Entry, ...]  # in record order; an omitted entry is absent
    triangles: int = 0  # the output's triangle count, which the summary names


def _exact(value: float) -> str:
    """``value`` short where that loses nothing, else every digit, so a printed
    maximum can never read as above the tolerance it met."""
    short = f"{value:g}"
    return short if float(short) == value else repr(value)


def escaped_ascii(text: str) -> str:
    """`text` with non-ASCII escaped (`\\xe9`), as a file field holds it."""
    return text.encode("ascii", "backslashreplace").decode("ascii")


def plural(n: int, one: str, many: str) -> str:
    return f"{n} {one if n == 1 else many}"


def ordinal(n: int) -> str:
    """``every 2nd DEM node``: the English ordinal; 1 reads ``every DEM node``."""
    if n == 1:
        return "every DEM node"
    suffix = "th" if n % 100 in (11, 12, 13) else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"every {n}{suffix} DEM node"


def _entry(name: str, value: str | int | float) -> Entry:
    """One entry: a float, or a ``_m`` or ``_deg`` number, is a measured float,
    an int a count, anything else text. D3 rule 3: a zero count is not in the
    file."""
    number: float | int | None = None
    if isinstance(value, float) or (name.endswith(("_m", "_deg")) and not isinstance(value, str)):
        number = float(value)
        text = _exact(number)
    elif isinstance(value, int):
        number, text = value, str(value)
    else:
        text = str(value)
    in_file = name in IN_FILE and not (name == "nodata_vertices_removed" and number == 0)
    return Entry(name, WORDING[name], text, number, in_file, name in INPUTS)


def _record(triangles: int, values: dict[str, Value]) -> RunRecord:
    unknown = set(values) - set(WORDING)
    if unknown:
        raise TypeError(f"no record entry named {sorted(unknown)}")
    entries = (_entry(n, v) for n in WORDING if (v := values.get(n)) is not None)
    return RunRecord(entries=tuple(entries), triangles=triangles)


def refined_record(*, triangles: int, **values: Value) -> RunRecord:
    """A run with ``--tolerance``. ``max_error_m`` is the measured figure; D2
    raises it to the difference at DEM nodes on a vertex where that is larger,
    so no DEM node inside the mesh is further from it than ``max_error_m``."""
    at = values.get("dem_nodes_at_vertices_max_error_m")
    measured = values["max_error_m"]
    if isinstance(at, float | int) and isinstance(measured, float | int):
        values["max_error_m"] = max(float(measured), float(at))
    return _record(triangles, values)


def stride_record(*, triangles: int, **values: Value) -> RunRecord:
    """A run without ``--tolerance``: the stride grid, sampled."""
    return _record(triangles, values)


def flat_record(*, triangles: int, crs: str | None = None) -> RunRecord:
    """A gallery fixture written with ``--flat``: z is not real, and says so."""
    return _record(triangles, {"crs": crs or None, "heights": FLAT_HEIGHTS})


def file_fields(record: RunRecord) -> list[tuple[str, str]]:
    """The ``.vtk`` FieldData and the ``.ply`` comments, ``(name, value)``."""
    return [(e.name, e.value) for e in record.entries if e.in_file]


def stats_rows(record: RunRecord) -> list[tuple[str, str, str]]:
    """Every entry as a ``--stats`` row, ``(wording, value, name)``."""
    return [(e.wording, e.value, e.name) for e in record.entries]


def summary(record: RunRecord) -> str:
    """D7: one line for a person, then one ``Warning:`` line per broken promise."""
    by = {e.name: e for e in record.entries}
    said = [f"{record.triangles} triangles."]
    tolerance, top = by.get("tolerance_m"), by.get("max_error_m")
    if tolerance is not None and top is not None:
        who = "node of the original DEM" if "dem_nodes_checked" in by else "DEM node"
        largest = f"{top.number:.5g}"
        if float(top.value) <= float(tolerance.value):
            said.append(
                f"Every {who} inside the mesh is within {tolerance.value} m of it "
                f"(largest difference {largest} m)."
            )
        else:
            said.append(
                f"The largest difference at a {who} inside the mesh is {largest} m, "
                f"above the tolerance of {tolerance.value} m."
            )
    removed = by.get("nodata_vertices_removed")
    if removed is not None and removed.number:
        n = int(removed.value)
        said.append(
            f"{plural(n, 'vertex', 'vertices')} on NoData cells "
            f"{'was removed with its' if n == 1 else 'were removed with their'} triangles."
        )
    if "heights" in by:
        said.append("The heights are not real: every z is 0 (--flat).")
    share = by.get("max_error_slope_share_of_tolerance")
    steep = f"{100 * float(share.value):.4g} % of" if share is not None else ""
    # 34, 4.3: the weight 1 / t is rounded, so up to 1 + 1e-12 is within.
    within = share is not None and float(share.value) <= 1.0 + 1e-12
    if within:
        said.append(
            "Nodes on steep ground are within the tolerance their slope allows "
            f"(largest error {steep} it)."
        )
    lines = [" ".join(said)]
    if share is not None and not within:
        lines.append(
            f"Warning: the largest error on steep ground is {steep} what its slope allows."
        )
    at, count = by.get("dem_nodes_at_vertices_max_error_m"), by.get("dem_nodes_at_vertices")
    if (
        at is not None
        and count is not None
        and tolerance is not None
        and (float(at.value) > float(tolerance.value))
    ):
        n = int(count.value)
        lines.append(
            f"Warning: {plural(n, 'DEM node', 'DEM nodes')} on a vertex "
            f"{'differs' if n == 1 else 'differ'} from it by up "
            f"to {at.value} m, more than the tolerance of {tolerance.value} m."
        )
    refused, worst = by.get("line_points_refused"), by.get("line_points_refused_max_error_m")
    if refused is not None and worst is not None and refused.number:
        n = int(refused.value)
        lines.append(
            f"Warning: {plural(n, 'point', 'points')} along the lines could not be added; "
            f"{'its' if n == 1 else 'their largest'} difference is {worst.value} m."
        )
    outside = by.get("dem_nodes_outside_mesh")
    if outside is not None and outside.number:
        n = int(outside.value)
        lines.append(
            f"Warning: {plural(n, 'DEM node', 'DEM nodes')} with data "
            f"{'lies' if n == 1 else 'lie'} outside the mesh; there should be none."
        )
    return "\n".join(lines)


def as_json(record: RunRecord, version: str, command: str) -> str:
    """D5: ``--record``'s text. The record's order, ASCII, one trailing newline."""
    obj: dict[str, str | float | int] = {"rasputin_version": version, "command": command}
    for e in record.entries:
        obj[e.name] = e.value if e.number is None else e.number
    return json.dumps(obj, indent=1, ensure_ascii=True) + "\n"
