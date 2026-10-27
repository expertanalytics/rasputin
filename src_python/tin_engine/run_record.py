"""The record of one ``rasputin mesh`` run, and every word it is printed in.

`docs/increments/25-plain-output.md`, *The design*. Pure: no numpy, no
``_core``, no typer, so the wording is tested without the extension. ``cli.py`` copies
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

#: Every entry, in record order, with its ``--stats`` wording. The Inputs
#: entries come first (``INPUTS``), then the mesh-file fields in the order the
#: file carries them, then the counts and checks.
WORDING = {
    "dem_grid": "DEM grid size and spacing",
    "dem_tiles": "The tiles or downloaded blocks used",
    "dem_seams": "Overlapping tiles that disagree, and by how much",
    "dem_vertical_unit": "Unit of z",
    "dem_file_crs": "The DEM file's own coordinate system, when it is in feet",
    "dem_crs": "The DEM's own coordinate system",
    "dem_transform": "Conversion from the DEM's coordinate system",
    "resampled_grid": "The grid the DEM was interpolated onto",
    "domain": "The domain file and shape",
    "domain_crs": "The domain's coordinate system",
    "domain_band_m": "Band the domain outline was simplified within, m",
    "domain_transform": "Conversion from the domain's coordinate system",
    "features": "Each features file: layer, meshing rules, features, lines, vertices",
    "features_crs": "Each features file's coordinate system",
    "features_transform": "Conversion from each features file's coordinate system",
    "features_repair_m": "Land-cover borders closer than this made one, narrower gaps filled",
    "features_merge_same_class": "Borders between land-cover polygons of one class dropped",
    "features_tolerance_m": (
        "Land-cover borders simplified, each at most this far from its border after the "
        "repair and the outline rule, each class's area kept (0 = off)"
    ),
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
    "threads": "Threads used",
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
    "land_cover_source": "Land-cover classes: the source and the year",
    "land_cover_credit": "Credit for the land-cover classes",
    "land_cover_licence_note": "Licence of the land-cover classes",
    "land_cover_cite": "Please cite for the land-cover classes",
    "land_cover_prefer_credit": "Credit for the classes taken from a preferred source",
    "land_cover_prefer_licence_note": "Licence of the classes taken from a preferred source",
    "land_cover_prefer_cite": "Please cite for the classes taken from a preferred source",
    "features_notice": "Credit for the features data",
    "features_licence": "Licence of the features data",
    "land_cover_codes": "What land_cover_code holds",
    "cells": "The cells layer",
    "reservoirs": "Reservoir polygons splitting the cells' lake share",
    "land_cover_prefer": (
        "Classes taken from a preferred source: inside the mesh, the land cover's area of the "
        "class, the preferred source's, and what went to neighbouring classes"
    ),
    "land_cover_repaired": "Land-cover polygons that crossed themselves, repaired",
    "reservoirs_repaired": "Reservoir polygons that crossed themselves, repaired",
    "catchments": "Sub-catchment units: their borders constraints, their ids the cells'",
    "catchments_units": "Sub-catchment units whose inner point lies in the domain",
    "catchments_band_m": "Band the sub-catchment borders were simplified within, m",
    "land_cover_vertices": "Land-cover vertices before and after clean-up",
    "land_cover_area_moved_m2": (
        "Land-cover area inside the outline that the outline rule gave to another polygon, m2"
    ),
    "triangles_asked": "Triangles asked for (about this many)",
    "triangles_cap": "Triangles at most",
    "budget_trials": "Tolerances tried for the triangle budget, and the triangles each gave",
    "budget_window_missed": (
        "The triangle budget's window was missed; the nearest count under the cap was written"
    ),
    "rms_error_m": "Root-mean-square height error",
    "error_nodes_measured": "DEM nodes inside the mesh the height errors are over",
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
    "triangles_flat_as_written": (
        "Triangles flat or turned over at the written coordinates; their edges were left "
        "out of the last edge-swap check"
    ),
    "edges_not_delaunay_as_written": (
        "Edges that fail the empty-circle test at the written coordinates and could not be swapped"
    ),
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
#: The mesh-file fields, and the CORINE credit, which the Copernicus data
#: policy asks for on data derived from CORINE.
#: ``land_cover_codes`` is in the file too, but the writers put it there.
IN_FILE = frozenset(_NAMES[_NAMES.index("crs") : _NAMES.index("features_notice") + 1])
SEAMS_AGREE = "none: the tiles agree where they overlap"


@dataclass(frozen=True, slots=True)
class Entry:
    name: str  # plain snake_case words, the unit as a suffix (_m, _deg)
    wording: str  # the label --stats prints
    value: str  # ASCII, already formatted
    number: float | int | None  # the bare number, for --record's JSON; None for text
    in_file: bool  # True only for the fields the mesh file carries
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
    an int a count, anything else text. A zero count of removed NoData
    vertices is left out of the file."""
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
    """A run with ``--tolerance``. ``max_error_m`` is the measured figure,
    raised to the difference at DEM nodes on a vertex where that is larger,
    so no DEM node inside the mesh is further from it than ``max_error_m``."""
    at = values.get("dem_nodes_at_vertices_max_error_m")
    measured = values["max_error_m"]
    if isinstance(at, float | int) and isinstance(measured, float | int):
        values["max_error_m"] = max(float(measured), float(at))
    return _record(triangles, values)


def stride_record(*, triangles: int, **values: Value) -> RunRecord:
    """A run without ``--tolerance``: the stride grid, sampled."""
    return _record(triangles, values)


def file_fields(record: RunRecord) -> list[tuple[str, str]]:
    """The ``.vtk`` FieldData and the ``.ply`` comments, ``(name, value)``."""
    return [(e.name, e.value) for e in record.entries if e.in_file]


def stats_rows(record: RunRecord) -> list[tuple[str, str, str]]:
    """Every entry as a ``--stats`` row, ``(wording, value, name)``."""
    return [(e.wording, e.value, e.name) for e in record.entries]


def summary(record: RunRecord) -> str:
    """One line for a person, then one ``Warning:`` line per broken promise."""
    by = {e.name: e for e in record.entries}
    said = [f"{record.triangles} triangles."]
    tolerance, top = by.get("tolerance_m"), by.get("max_error_m")
    asked, cap, tried = by.get("triangles_asked"), by.get("triangles_cap"), by.get("budget_trials")
    shown = tolerance.value if tolerance is not None else ""
    figures = ".5g"
    if asked is not None and cap is not None and tried is not None and tolerance is not None:
        shown, figures = f"{float(tolerance.value):.3g}", ".3g"
        n, found = record.triangles, f"tolerance {shown} m"
        trials = plural(tried.value.count(";") + 1, "trial", "trials")
        if top is not None and top.number == 0:
            said = [
                f"{n} triangles: every DEM node is a vertex, the most this DEM gives "
                f"(about {asked.value} asked for)."
            ]
        elif "budget_window_missed" in by:
            said = [
                f"{n} triangles, the most under the cap of {cap.value} (about {asked.value} "
                f"asked for; no tolerance gave a count in between): {found}, {trials}."
            ]
        else:
            said = [
                f"{n} triangles, about the {asked.value} asked for (at most {cap.value}): "
                f"{found} found in {trials}."
            ]
    if tolerance is not None and top is not None:
        who = "node of the original DEM" if "dem_nodes_checked" in by else "DEM node"
        largest = format(top.number, figures)
        rms = by.get("rms_error_m")
        mean = f"root-mean-square {float(rms.value):.2f} m" if rms is not None else ""
        if float(top.value) <= float(tolerance.value):
            said.append(
                f"Every {who} inside the mesh is within {shown} m of it "
                f"(largest difference {largest} m{', ' if mean else ''}{mean})."
            )
        else:
            said.append(
                f"The largest difference at a {who} inside the mesh is {largest} m, "
                f"above the tolerance of {shown} m{f' ({mean})' if mean else ''}."
            )
    removed = by.get("nodata_vertices_removed")
    if removed is not None and removed.number:
        n = int(removed.value)
        said.append(
            f"{plural(n, 'vertex', 'vertices')} on NoData cells "
            f"{'was removed with its' if n == 1 else 'were removed with their'} triangles."
        )
    share = by.get("max_error_slope_share_of_tolerance")
    steep = f"{100 * float(share.value):.4g} % of" if share is not None else ""
    # The weight 1 / t is rounded, so up to 1 + 1e-12 is within.
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
    """``--record``'s text. The record's order, ASCII, one trailing newline."""
    obj: dict[str, str | float | int] = {"rasputin_version": version, "command": command}
    for e in record.entries:
        obj[e.name] = e.value if e.number is None else e.number
    return json.dumps(obj, indent=1, ensure_ascii=True) + "\n"
