"""Where the land-cover clean-up spends its time, and whether each candidate
change gives the same bytes. Run as

    .venv/bin/python docs/increments/35-probes/clean_probe.py numedalslagen|lagan

from the tree. It reads the case's CORINE source and domain as `rasputin mesh`
does (the quick check's arguments, docs/benchmarks/quick/cases.toml), builds
the `cover` list `_Tally._take` hands to `_Tally._clean`, then runs the step's
sub-steps one by one with a timer, master's way, and beside each the candidate:

  C1  the reprojection of every polygon in one pyproj call (Lagan only);
  C2  no clip to the read region for a polygon the region contains properly;
  C3  `_add`'s `covers` skipped for a line `_untouched` already passed;
  C4  `_add`'s `intersection(line, dom)` skipped for a line whose interior lies
      in the domain's interior and whose ends lie in its closure (DE-9IM
      "T**F**FF*" with the domain first);
  C5  `_outline_lines`'s STRtree of the outline built once per step;
  C6  the clip to the domain done part by part.

Each candidate's output is compared with master's as WKB (or as float64
bits). Timings are one run each, on whatever power state `pmset -g batt`
reports; they rank the sub-steps, they are not an acceptance figure."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import shapely
from shapely.geometry import MultiPolygon, Polygon

from tin_engine import feature_input as fi
from tin_engine.border_simplify import simplify_borders
from tin_engine.crs import parse_crs, reprojector, same_crs
from tin_engine.feature_input import (
    CLASS_MAPS,
    FeatureRequest,
    FeatureSource,
    _outline_lines,
    _polygonal,
    _untouched,
    read_source,
    snap_to_outline,
    source_region,
)
from tin_engine.io.domain_file import read_domain

DATA = Path("/Users/skavhaug/projects/rasputin_data")
CASES = {
    "numedalslagen": (
        DATA / "numedalslagen_outline_nve.geojson",
        DATA / "corine2018_dtm10_utm33.gpkg",
        "corine2018",
        None,
        "EPSG:25833",
    ),
    "lagan": (
        DATA / "sweden_smhi_svar/lagan_mouth_svar2022_3006.geojson",
        DATA / "corine_sql/u2018_clc2018_v2020_20u1_geoPackage/DATA/U2018_CLC2018_V2020_20u1.gpkg",
        "U2018_CLC2018_V2020_20u1",
        "EPSG:3035",
        "EPSG:3006",
    ),
}
ASK = {"repair_m": 1.0, "merge_same_class": True, "tolerance_m": 50.0, "outline_snap_m": 5.0}


def wkb(geoms):
    return [shapely.to_wkb(g) for g in geoms]


def read_cover(case: str):
    """The domain in the DEM's CRS, that CRS, the read region, and the cover
    list as `_Tally._take` builds it (C1 measured on the way)."""
    domain_path, source_path, layer, crs, dem_crs = CASES[case]
    domain = read_domain(domain_path).to_crs(dem_crs)
    source = FeatureSource(path=source_path, class_map=CLASS_MAPS["corine"], layer=layer, crs=crs)
    FeatureRequest(sources=(source,), **ASK)
    t = time.perf_counter()
    read = read_source(
        source.path, source.layer, source.class_map.attribute,
        lambda own: source_region(domain, dem_crs, own).bounds, source.crs,
    )
    print(f"read: {len(read.rows)} rows in {time.perf_counter() - t:.2f} s")
    src = parse_crs(read.crs)
    move = None if same_crs(src, dem_crs) else reprojector(src, dem_crs)
    region = source_region(domain, dem_crs, dem_crs)
    shapely.prepare(region)
    dom = domain.polygon
    shapely.prepare(dom)

    # --- the cover list, as _take builds it; C1 beside it ---
    rows = [
        (fid, shapely.force_2d(g), v)
        for fid, g, v in read.rows
        if g is not None and not g.is_empty
    ]
    rows = [r for r in rows if r[1].geom_type in ("Polygon", "MultiPolygon")]
    geoms = [g for _, g, _ in rows]
    n_vertices = int(shapely.get_num_coordinates(np.array(geoms, dtype=object)).sum())
    print(f"cover: {len(geoms)} polygons, {n_vertices} vertices before the clip")
    if move is not None:
        t = time.perf_counter()
        per = [shapely.transform(g, move) for g in geoms]
        t_per = time.perf_counter() - t
        t = time.perf_counter()
        batched = list(shapely.transform(np.array(geoms, dtype=object), move))
        t_batch = time.perf_counter() - t
        same_bits = all(
            np.array_equal(shapely.get_coordinates(a), shapely.get_coordinates(b))
            and shapely.to_wkb(a) == shapely.to_wkb(b)
            for a, b in zip(per, batched, strict=True)
        )
        print(f"C1 reprojection: per polygon {t_per:.2f} s, one call {t_batch:.2f} s, "
              f"same bits: {same_bits}")
        moved = per
    else:
        moved = geoms
    cover = [(fid, 0, int(v), g) for (fid, _, v), g in zip(rows, moved, strict=True)]
    return domain, dem_crs, region, cover


def main(case: str) -> None:
    domain, _dem_crs, region, cover = read_cover(case)
    dom = domain.polygon
    ask = FeatureRequest(sources=(), **ASK)

    # --- _clean, sub-step by sub-step ---
    t_step = time.perf_counter()
    t = time.perf_counter()
    clipped = [_polygonal(shapely.intersection(g, region)) for *_, g in cover]
    t_region = time.perf_counter() - t
    # C2: which polygons the region contains properly, and whether their clip is themselves
    t = time.perf_counter()
    inside = shapely.contains_properly(region, np.array([g for *_, g in cover], dtype=object))
    t_inside = time.perf_counter() - t
    t = time.perf_counter()
    c2 = [
        _polygonal(g) if ok else _polygonal(shapely.intersection(g, region))
        for (*_, g), ok in zip(cover, inside, strict=True)
    ]
    t_c2 = time.perf_counter() - t
    same_c2 = wkb(x for x in clipped if x is not None) == wkb(x for x in c2 if x is not None)
    eq_c2 = [
        shapely.equals_exact(a, b, 0) for a, b, ok in zip(clipped, c2, inside, strict=True) if ok
    ]
    print(f"clip to region: {t_region:.2f} s ({len(cover)} calls); "
          f"C2: {int(inside.sum())} inside ({t_inside:.3f} s test), {t_c2:.2f} s, "
          f"same WKB: {same_c2}, same point sets: {all(eq_c2)}")
    items = [c[:3] for c, g in zip(cover, clipped, strict=True) if g is not None]
    polys = np.array([g for g in clipped if g is not None], dtype=object)
    t = time.perf_counter()
    polys = shapely.coverage_clean(
        polys, snapping_distance=ask.repair_m, gap_width=ask.repair_m, merge_strategy="min_area"
    )
    print(f"coverage_clean: {time.perf_counter() - t:.2f} s")
    t = time.perf_counter()
    groups: dict[int, list[int]] = {}
    for i, (_, _, code) in enumerate(items):
        groups.setdefault(code, []).append(i)
    first = sorted(g[0] for g in groups.values())
    merged = [shapely.coverage_union_all(polys[groups[items[i][2]]]) for i in first]
    polys = np.array(merged, dtype=object)
    items = [items[i] for i in first]
    print(f"merge by class: {time.perf_counter() - t:.2f} s ({len(first)} classes)")

    # the outline rule, with its GEOS calls timed through wrappers of coarse calls
    tally: dict[str, float] = {}
    def timed(name, fn, by_line=False):
        def wrap(*a, **k):
            t = time.perf_counter()
            out = fn(*a, **k)
            key = f"{name}@{sys._getframe(1).f_lineno}" if by_line else name  # coarse calls only
            tally[key] = tally.get(key, 0.0) + time.perf_counter() - t
            return out
        return wrap
    saved = {
        "union_all": shapely.union_all, "make_valid": shapely.make_valid,
        "_loops": fi._loops, "_chains": fi._chains,
    }
    shapely.union_all = timed("union_all", shapely.union_all, True)
    shapely.make_valid = timed("make_valid", shapely.make_valid, True)
    fi._loops = timed("_loops", fi._loops)
    fi._chains = timed("_chains", fi._chains)
    ring = fi._Outline.ring
    fi._Outline.ring = timed("_Outline.ring", ring)
    t = time.perf_counter()
    snapped = snap_to_outline(list(polys), dom, ask.outline_snap_m)
    t_snap = time.perf_counter() - t
    shapely.union_all, shapely.make_valid = saved["union_all"], saved["make_valid"]
    fi._loops, fi._chains, fi._Outline.ring = saved["_loops"], saved["_chains"], ring
    inner = {k: round(v, 2) for k, v in tally.items()}
    print(f"snap_to_outline: {t_snap:.2f} s; inside it {inner}")
    # the loops' union and clip alone (line 687)
    polygons = snapped.polygons

    t = time.perf_counter()
    inside_dom = [
        Polygon() if (c := _polygonal(shapely.intersection(g, dom))) is None else c
        for g in polygons
    ]
    t_dom = time.perf_counter() - t
    # C6: part by part
    t = time.perf_counter()
    c6 = []
    n_kept = n_cut = 0
    for g in polygons:
        parts = shapely.get_parts(g)
        keep = shapely.contains_properly(dom, parts)
        cut = ~keep & shapely.intersects(dom, parts)
        n_kept += int(keep.sum())
        n_cut += int(cut.sum())
        rest = None
        if cut.any():
            rest = _polygonal(shapely.intersection(MultiPolygon(list(parts[cut])), dom))
        whole = [*parts[keep], *([] if rest is None else shapely.get_parts(rest))]
        c6.append(Polygon() if not whole else whole[0] if len(whole) == 1 else MultiPolygon(whole))
    t_c6 = time.perf_counter() - t
    same_c6 = wkb(inside_dom) == wkb(c6)
    sets_c6 = all(
        sorted(wkb(shapely.get_parts(a))) == sorted(wkb(shapely.get_parts(b)))
        for a, b in zip(inside_dom, c6, strict=True)
    )
    print(f"clip to domain: {t_dom:.2f} s ({len(polygons)} calls); C6: {t_c6:.2f} s, "
          f"{n_kept} parts kept whole, {n_cut} cut; same WKB: {same_c6}, same part sets: {sets_c6}")

    t = time.perf_counter()
    simplified = simplify_borders(inside_dom, ask.tolerance_m, ask.repair_m).polygons
    print(f"simplify_borders: {time.perf_counter() - t:.2f} s")
    t = time.perf_counter()
    lines_of = tuple(_outline_lines(g, dom) for g in simplified)
    t_ol = time.perf_counter() - t
    # C5: the tree once
    t = time.perf_counter()
    rims = [shapely.get_coordinates(r) for r in shapely.get_rings(dom)]
    segments = shapely.linestrings(
        np.concatenate([np.stack([r[:-1], r[1:]], axis=1) for r in rims])
    )
    shapely.STRtree(segments)
    t_tree = time.perf_counter() - t
    print(f"_outline_lines: {t_ol:.2f} s ({len(simplified)} calls); "
          f"C5: one outline tree costs {t_tree:.3f} s, "
          f"so about {t_tree * (len(simplified) - 1):.2f} s saved")

    # _add's per-line work, master's way and C3/C4
    all_lines = [ln for ls in lines_of for ln in ls]
    t = time.perf_counter()
    untouched = [_untouched(ln, dom) for ln in all_lines]
    t_unt = time.perf_counter() - t
    t = time.perf_counter()
    cut_master = [
        shapely.get_parts(shapely.intersection(ln, dom))
        for ln, u in zip(all_lines, untouched, strict=True)
        if not u
    ]
    t_cut = time.perf_counter() - t
    t = time.perf_counter()
    covers = [dom.covers(ln) for ln in all_lines]
    t_cov = time.perf_counter() - t
    t = time.perf_counter()
    covers_c3 = [True if u else dom.covers(ln) for ln, u in zip(all_lines, untouched, strict=False)]
    t_c3 = time.perf_counter() - t
    print(f"_add lines: {len(all_lines)}; _untouched {t_unt:.2f} s, {sum(untouched)} pass; "
          f"intersection of the rest {t_cut:.2f} s; covers {t_cov:.2f} s; "
          f"C3 covers only where not untouched {t_c3:.2f} s, same: {covers == covers_c3}")
    rest = [ln for ln, u in zip(all_lines, untouched, strict=False) if not u]
    t = time.perf_counter()
    pat = [shapely.relate_pattern(dom, ln, "T**F**FF*") for ln in rest]
    t_pat = time.perf_counter() - t
    t = time.perf_counter()
    pat2 = [shapely.relate_pattern(ln, dom, "TFF**F***") for ln in rest]
    t_pat2 = time.perf_counter() - t
    same_c4 = all(
        (len(parts) == 1 and shapely.to_wkb(parts[0]) == shapely.to_wkb(ln))
        for ln, parts, ok in zip(rest, cut_master, pat, strict=False) if ok
    )
    print(f"C4: of {len(rest)} not untouched, {sum(pat)} match the pattern "
          f"(domain first {t_pat:.2f} s, line first {t_pat2:.2f} s, agree: {pat == pat2}); "
          f"their intersection is the line itself: {same_c4}; "
          f"{sum(1 for p in cut_master if len(p) != 1)} of the rest split into parts")
    print(f"whole step re-run here: {time.perf_counter() - t_step:.2f} s "
          f"(includes the candidates' extra work); area_changed {snapped.area_changed:.3f}")


if __name__ == "__main__":
    main(sys.argv[1])
