"""`rasputin draw`: increment 6b-ii's command, and the only composition root.

Committed RED, before `cli.py` gains the command. `tin_engine.cli` already
exists, so it is imported at module scope and the intended failure is Click's
`No such command 'draw'` -- exit code 2 with no picture. Every test that
asserts a *refusal* also asserts what the refusal says, because an assertion of
the form `exit_code != 0` alone would pass against the missing command and
measure nothing.

This is where the two halves of increment 6 meet. `viz/` never imports `_core`
and `_core` never sees a path, a file or a CRS; `cli.py` is the single module
that knows both, and its job is:

1. look a fixture up in `viz.fixtures.GALLERY`;
2. map the fixture's string roles onto `_core.ChainRole` and call `build_pslg`,
   then `node`, then `triangulate`;
3. hand the **most-processed graph that exists** to `build_scene` as the
   `PslgLike`, together with the closed-role vocabulary THAT graph's roles
   speak -- strings for a fixture, `ChainRole` values for a `NodedPslg`;
4. pass the status name, `describe(status)` and the engine's message to
   `render_svg` as text;
5. resolve and validate the output path, and write the bytes.

Step 3 is the one the design does not state, and `TestFailurePresentation`
below is why it has to be that way: `06-cdt-viewer.md` says the `degenerate`
fixture shows "the `DegenerateGeometry` presentation", but an all-collinear ring
never reaches `triangulate` -- the PSLG validator rejects it first, so there is
no `Pslg` object at all and nothing for `build_scene` to be given unless the
fixture is itself a `PslgLike`. Measured, on this tree -- the third element of
a chain spec is a mask of property bits, so the empty set is `0`, and the
pre-migration `False` raises `ValueError` at the binding instead of returning
diagnostics at all:

    [(d.error, d.message) for d in build_pslg(
        collinear, [(range(4), ChainRole.Outer, 0)]).diagnostics]
    -> [(<PslgError.DegenerateRing: 6>,
         'chain 0 is declared Outer but every vertex is collinear')]

The presentation the design wants is still seen -- the input drawn alone with
the engine's own words in the header band -- but the words come from a
`PslgDiagnostic` rather than from a `CdtStatus`, and the assertions below take
them from the engine at run time rather than hard-coding either.

Path handling is tested as a hostile boundary, per the python skill: traversal
out of an explicitly provided parent, and a symlink that resolves back inside it.

AMENDED AT INCREMENT 5c, which is the increment this file exists to observe:
`rasputin draw not-noded` stops drawing zero triangles and starts drawing a
noded crossing. Step 3 is the part that a status assertion cannot see. A
`NodedPslg`'s roles are `ChainRole` values and `CLOSED_ROLES` was a tuple of
STRINGS; `ChainRole.Outer == "outer"` is False, so passing the noded graph
while leaving the vocabulary alone closes no ring, and every ring's closing
edge becomes a `MASKED_EDGE_WITHOUT_CHAIN` -- the whole gallery in the alarm
colour over a correct mesh, with `CdtStatus.Ok` in the band. `TestNodedCrossing`
below holds that line twice: once at the picture, where findings must be zero,
and once structurally, against the source and the vocabulary the composition
root actually paired.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from typer.testing import CliRunner

import tin_engine.cli as cli
from tin_engine.cli import ROLES, app

SVG_NS = "http://www.w3.org/2000/svg"

GALLERY_NAMES = (
    "catchment",
    "sliver-fan",
    "corner-hole",
    "hole-in-hole",
    "breakline-chain",
    "river",
    "not-noded",
    "degenerate",
)

DEFAULT_LABEL_LIMIT = 500

#: `cli.py`'s default, restated here rather than imported: a test that took the
#: number from the module under test would agree with any value it was given.
#: A millimetre on the ground, the gallery's coordinates being UTM-shaped.
DEFAULT_SNAP_SPACING = 1e-3

#: Fine enough that a UTM easting overflows the snap lattice -- `kMaxGridIndex`
#: is 2**51, so the representable coordinate at spacing s is s * 2.25e15. This
#: is how a NODER refusal is reached from the CLI with no fixture of its own,
#: and it is the lower bound that makes `--snap-spacing` an option rather than
#: a constant hidden in a header.
OVERFLOWING_SNAP_SPACING = 1e-12

runner = CliRunner(env={"NO_COLOR": "1", "TERM": "dumb"})

ANSI = re.compile(r"\x1b\[[0-9;]*m")
BOX = re.compile(r"[─-╿]")


def plain(text: str) -> str:
    """Output as a reader sees it: no colour, no box rule, no line wrapping.

    Typer renders errors inside a Rich panel whose width is the terminal's, so
    a message asserted verbatim would fail on a narrow one and pass on a wide
    one -- a flaky test by construction. Collapsing the box and the whitespace
    makes every assertion below independent of where Rich chose to wrap, while
    leaving words and numbers intact, which is all any of them look at.
    """
    return " ".join(BOX.sub(" ", ANSI.sub("", text)).split())


def invoke(*args: str) -> Any:
    return runner.invoke(app, ["draw", *args])


def written(path: Path) -> ET.Element:
    assert path.is_file(), f"{path} was not written"
    return ET.fromstring(path.read_text(encoding="utf-8"))


def group_text(document: ET.Element, gid: str) -> str:
    groups = [g for g in document.iter(f"{{{SVG_NS}}}g") if g.get("id") == gid]
    assert len(groups) == 1, f"expected exactly one <g id={gid!r}>, found {len(groups)}"
    return " ".join("".join(groups[0].itertext()).split())


def elements(document: ET.Element, gid: str, name: str) -> list[ET.Element]:
    groups = [g for g in document.iter(f"{{{SVG_NS}}}g") if g.get("id") == gid]
    return [e for g in groups for e in g.iter(f"{{{SVG_NS}}}{name}")]


def gallery() -> Any:
    import importlib

    return importlib.import_module("tin_engine.viz.fixtures").GALLERY


def core_verdict(name: str, spacing: float = DEFAULT_SNAP_SPACING) -> tuple[bool, list[str]]:
    """What the engine actually says about a fixture, computed here.

    The oracle for a passthrough claim is the engine's own words, so this
    re-runs the composition the command performs rather than hard-coding a
    status. It is deliberately *not* an independent derivation -- there is
    nothing independent to derive; the claim under test is only that whatever
    the engine said reaches the picture.

    Returns `(drawable, words)`: whether a mesh was produced at all, and the
    strings the header band must carry. Three refusal layers can answer now --
    the validator, the noder and the backend -- and the caller does not get to
    know which, which is the point: the claim under test is only that whatever
    the engine said reaches the picture.
    """
    import tin_engine._core as core

    fixture = gallery()[name]
    roles = {
        "outer": core.ChainRole.Outer,
        "hole": core.ChainRole.Hole,
        "breakline": core.ChainRole.Breakline,
    }
    chains = [
        ([int(i) for i in fixture.indices_of(c)], roles[chain.role], int(chain.properties))
        for c, chain in enumerate(fixture.chains)
    ]
    result = core.build_pslg(np.asarray(fixture.vertices), chains)
    if not result.ok:
        return False, [d.error.name for d in result.diagnostics]
    assert result.pslg is not None
    noded = core.node(result.pslg, spacing)
    if not noded.ok():
        return False, [noded.status.name]
    outcome = core.triangulate(noded.pslg)
    return outcome.ok(), [outcome.status.name]


class TestCommandSurface:
    def test_the_command_is_registered(self) -> None:
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "draw" in plain(result.output)

    def test_list_names_every_fixture(self) -> None:
        result = invoke("--list")
        assert result.exit_code == 0
        listing = plain(result.output)
        for name in GALLERY_NAMES:
            assert name in listing, f"--list does not name {name}"

    def test_list_says_what_each_fixture_is_for(self) -> None:
        # Name alone is not a gallery. The description is the reason the fixture
        # exists, and a listing without it is a list of eight opaque nouns.
        result = invoke("--list")
        listing = plain(result.output)
        for name, fixture in gallery().items():
            first = fixture.description.split()[0]
            assert first in listing, f"--list drops {name}'s description"

    def test_list_needs_no_fixture_name(self) -> None:
        # `--list` is asked before the user knows a name to give, so the
        # argument has to be optional rather than dummy-valued.
        assert invoke("--list").exit_code == 0

    def test_an_unknown_fixture_is_refused_by_name(self) -> None:
        result = invoke("no-such-fixture", "--out", "/dev/null")
        assert result.exit_code != 0
        message = plain(result.output)
        assert "no-such-fixture" in message
        assert any(name in message for name in GALLERY_NAMES), "no alternatives offered"

    def test_a_missing_fixture_name_points_at_the_listing(self) -> None:
        # The name is optional so that `--list` can be asked before the user
        # knows one, which means "no name" is a usage error the command has to
        # write itself. It says how to find a name; asserting that rather than
        # just a non-zero exit is what keeps this test from passing against a
        # `draw` command that does not exist at all.
        result = invoke()
        assert result.exit_code != 0
        message = plain(result.output)
        assert "Usage" in message
        assert "--list" in message

    def test_the_cli_maps_every_role_the_enum_has(self) -> None:
        # The claim `test_viz_svg.py::TestModuleIsolation` cannot make, because
        # that suite may not import `_core`: the mapping this module owns is
        # total both ways. A role the enum gains and `ROLES` does not is a
        # fixture vocabulary that silently cannot express it.
        import tin_engine._core as core

        # `__members__` rather than iteration: a pybind11 enum type is not
        # iterable, so `set(ChainRole)` is a TypeError and not a set.
        assert set(ROLES.values()) == set(core.ChainRole.__members__.values())


class TestOutput:
    def test_out_writes_a_parseable_svg(self, tmp_path: Path) -> None:
        target = tmp_path / "catchment.svg"
        result = invoke("catchment", "--out", str(target))
        assert result.exit_code == 0, plain(result.output)
        assert written(target).tag == f"{{{SVG_NS}}}svg"

    def test_the_catchment_actually_triangulates(self, tmp_path: Path) -> None:
        # Able to fail on its own: if the catchment stopped producing a mesh,
        # the happy path above would still pass on the failure presentation.
        target = tmp_path / "catchment.svg"
        assert invoke("catchment", "--out", str(target)).exit_code == 0
        assert elements(written(target), "triangles", "polygon") != []

    def test_a_well_formed_fixture_draws_no_disagreement(self, tmp_path: Path) -> None:
        # `closed_roles` is the composition root's to supply -- a `Pslg` never
        # stores a ring's closing edge, and `viz/` cannot name `ChainRole` to
        # work out which chains are rings. Forget it and every ring closure
        # becomes a `MASKED_EDGE_WITHOUT_CHAIN`: a picture covered in alarms on
        # input the engine accepted without complaint. Measured against a
        # prototype: with `closed_roles=()` every other test in both new suites
        # still passed, so this is the only assertion that holds that line.
        target = tmp_path / "catchment.svg"
        assert invoke("catchment", "--out", str(target)).exit_code == 0
        document = written(target)
        alarms = [
            line
            for line in elements(document, "edges", "line")
            if "finding" in (line.get("class") or "")
        ]
        assert alarms == []
        assert re.search(r"(?i)\bfindings\b\D{0,3}0\b", group_text(document, "header"))

    def test_the_river_fixture_draws_its_property_stroke(self, tmp_path: Path) -> None:
        # `SvgStyle.property_strokes` defaults to `()` and must stay `()` --
        # a default naming `river` would be policy in the module that declares
        # it holds none -- so the precedence list is the composition root's to
        # supply, exactly as `closed_roles` is. Forget it and the river fixture
        # renders as a plain breakline: a picture that looks entirely correct
        # and has lost the one thing the fixture exists to show.
        target = tmp_path / "river.svg"
        result = invoke("river", "--out", str(target))
        assert result.exit_code == 0, plain(result.output)
        tokens = {
            token
            for line in elements(written(target), "edges", "line")
            for token in (line.get("class") or "").split()
        }
        assert "river" in tokens

    def test_a_fixture_with_no_properties_draws_no_property_stroke(
        self, tmp_path: Path
    ) -> None:
        # Able to fail on its own: a renderer that put the token on every
        # constrained edge would pass the test above. `breakline-chain` is the
        # same picture as `river` with the mask cleared, which is what makes
        # the pair a controlled comparison rather than two unrelated fixtures.
        target = tmp_path / "breakline-chain.svg"
        assert invoke("breakline-chain", "--out", str(target)).exit_code == 0
        tokens = {
            token
            for line in elements(written(target), "edges", "line")
            for token in (line.get("class") or "").split()
        }
        assert "river" not in tokens

    def test_with_no_out_it_prints_the_path_it_wrote(self, tmp_path: Path) -> None:
        # "The common case is show me this now": a temp file and a path on
        # stdout, because a picture the user cannot find is not a picture.
        result = invoke("catchment")
        assert result.exit_code == 0, plain(result.output)
        match = re.search(r"\S+\.svg", plain(result.output))
        assert match is not None, "no path printed"
        assert written(Path(match.group(0))).tag == f"{{{SVG_NS}}}svg"

    def test_the_title_reaches_the_picture(self, tmp_path: Path) -> None:
        # A CRS name goes here, as text, having never been near C++.
        target = tmp_path / "titled.svg"
        assert invoke("catchment", "--out", str(target), "--title", "Gaula, UTM 33N").exit_code == 0
        assert "Gaula, UTM 33N" in group_text(written(target), "header")

    def test_no_delaunay_still_draws(self, tmp_path: Path) -> None:
        # One bound bool, and the cheapest intuition in the increment: the user
        # sees both triangulations of the same input. This asserts only that the
        # picture survives the flag; that the two pictures DIFFER is the test
        # below, which is what actually holds the flag to the backend.
        target = tmp_path / "raw.svg"
        result = invoke("catchment", "--out", str(target), "--no-delaunay")
        assert result.exit_code == 0, plain(result.output)
        assert elements(written(target), "triangles", "polygon") != []

    def test_no_delaunay_draws_a_different_triangulation(self, tmp_path: Path) -> None:
        # The flag has to reach `triangulate`, and nothing above holds it there:
        # hard-wiring `delaunay=True` in `cli.py` left the whole suite passing.
        # The fixtures are fixed data, so this is not a gamble on some quad
        # happening to be flippable -- measured on this tree, `catchment`'s two
        # triangle sets are 26 triangles each with a symmetric difference of 38,
        # and it is the only fixture that both triangulates and differs by
        # enough to be worth asserting on (`corner-hole`'s difference is 0).
        drawn = {}
        for label, extra in (("delaunay", ()), ("raw", ("--no-delaunay",))):
            target = tmp_path / f"{label}.svg"
            result = invoke("catchment", "--out", str(target), *extra)
            assert result.exit_code == 0, plain(result.output)
            polygons = elements(written(target), "triangles", "polygon")
            assert polygons != [], f"{label} drew no triangle"
            drawn[label] = {polygon.get("points") for polygon in polygons}
        assert len(drawn["raw"]) == len(drawn["delaunay"]), "the two fill the same domain"
        assert drawn["raw"] != drawn["delaunay"], "--no-delaunay drew the Delaunay mesh"

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_fixture_in_the_gallery_draws(self, name: str, tmp_path: Path) -> None:
        # Including the three deliberate failures: a fixture that cannot be
        # drawn is a fixture nobody will look at.
        target = tmp_path / f"{name}.svg"
        result = invoke(name, "--out", str(target))
        assert result.exit_code == 0, plain(result.output)
        document = written(target)
        assert document.tag == f"{{{SVG_NS}}}svg"
        assert elements(document, "edges", "line") != [], f"{name} rendered blank"


class TestLabels:
    """Risk 5: `--labels` on a large mesh produces an enormous unreadable file.

    The guard is a hard error naming the count, not a warning, and not a
    silently unlabelled picture.
    """

    def test_labels_are_off_by_default(self, tmp_path: Path) -> None:
        target = tmp_path / "plain.svg"
        assert invoke("catchment", "--out", str(target)).exit_code == 0
        assert elements(written(target), "labels", "text") == []

    def test_labels_draw_one_index_per_vertex(self, tmp_path: Path) -> None:
        target = tmp_path / "labelled.svg"
        result = invoke("catchment", "--out", str(target), "--labels")
        assert result.exit_code == 0, plain(result.output)
        assert elements(written(target), "labels", "text") != []

    def test_above_the_limit_it_refuses_and_names_the_count(self, tmp_path: Path) -> None:
        # The gallery is deliberately small, so the limit is lowered rather than
        # the fixture grown -- the flag exists precisely so the refusal can be
        # exercised without a 40 000-triangle fixture in the repository.
        count = len(np.asarray(gallery()["catchment"].vertices))
        assert count > 3, "the fixture cannot exceed a limit of 3"
        target = tmp_path / "unreadable.svg"
        result = invoke(
            "catchment", "--out", str(target), "--labels", "--label-limit", "3"
        )
        assert result.exit_code != 0
        message = plain(result.output)
        assert str(count) in message, "the refusal does not name the vertex count"
        assert "3" in message, "the refusal does not name the limit"

    def test_the_refusal_writes_nothing(self, tmp_path: Path) -> None:
        # A hard error, so no half-written file and no unreadable one. The
        # message is asserted alongside the absent file because "nothing was
        # written" is also true of a command that was never reached, and a test
        # that cannot tell those apart measures nothing.
        count = len(np.asarray(gallery()["catchment"].vertices))
        target = tmp_path / "unreadable.svg"
        result = invoke("catchment", "--out", str(target), "--labels", "--label-limit", "3")
        assert str(count) in plain(result.output)
        assert not target.exists()

    def test_the_limit_does_not_fire_without_labels(self, tmp_path: Path) -> None:
        # The limit bounds the *labels*, not the picture: an unlabelled mesh of
        # any size is perfectly readable and must still draw.
        target = tmp_path / "big.svg"
        result = invoke("catchment", "--out", str(target), "--label-limit", "1")
        assert result.exit_code == 0, plain(result.output)

    def test_the_default_limit_is_five_hundred(self) -> None:
        result = runner.invoke(app, ["draw", "--help"])
        assert result.exit_code == 0
        assert str(DEFAULT_LABEL_LIMIT) in plain(result.output)


class TestPathHandling:
    """The CLI is the only place a path exists at all, so it is the only place
    path validation can happen. Treated as a hostile boundary per the python
    skill: `resolve()`, an explicitly provided parent, and no symlink following.
    """

    def test_a_path_inside_the_permitted_parent_is_written(self, tmp_path: Path) -> None:
        target = tmp_path / "inside.svg"
        result = invoke(
            "catchment", "--out", str(target), "--out-parent", str(tmp_path)
        )
        assert result.exit_code == 0, plain(result.output)
        assert target.is_file()

    def test_traversal_out_of_the_permitted_parent_is_refused(self, tmp_path: Path) -> None:
        parent = tmp_path / "permitted"
        parent.mkdir()
        escape = parent / ".." / "escaped.svg"
        result = invoke("catchment", "--out", str(escape), "--out-parent", str(parent))
        assert result.exit_code != 0
        assert "outside" in plain(result.output).lower()
        assert not (tmp_path / "escaped.svg").exists()

    def test_an_absolute_path_elsewhere_is_refused(self, tmp_path: Path) -> None:
        parent = tmp_path / "permitted"
        parent.mkdir()
        result = invoke(
            "catchment", "--out", "/tmp/rasputin-escape.svg", "--out-parent", str(parent)
        )
        assert result.exit_code != 0
        assert "outside" in plain(result.output).lower()

    def test_a_symlinked_target_is_refused(self, tmp_path: Path) -> None:
        # The link points back INSIDE the permitted parent, so containment alone
        # accepts it: only a `is_symlink()` check refuses. That is the case the
        # design means by "no symlink following", and the one a containment-only
        # implementation passes by accident.
        parent = tmp_path / "permitted"
        parent.mkdir()
        real = parent / "real.svg"
        real.write_text("not an svg", encoding="utf-8")
        link = parent / "link.svg"
        link.symlink_to(real)
        result = invoke("catchment", "--out", str(link), "--out-parent", str(parent))
        assert result.exit_code != 0
        assert "symlink" in plain(result.output).lower()
        assert real.read_text(encoding="utf-8") == "not an svg"

    def test_a_missing_directory_is_a_message_not_a_traceback(self, tmp_path: Path) -> None:
        result = invoke("catchment", "--out", str(tmp_path / "nope" / "x.svg"))
        assert result.exit_code != 0
        assert result.exception is None or isinstance(result.exception, SystemExit)
        assert "nope" in plain(result.output)


class TestFailurePresentation:
    """The three deliberate failure fixtures, seen rather than assumed.

    `06-cdt-viewer.md`: "a failure presentation nobody has looked at is a failure
    presentation that is wrong." All three exit 0 and produce a picture -- the
    exit code reports whether a picture was drawn, and a drawn failure is still
    a drawn picture. What the engine refused is in the header band.

    TWO of them as of increment 5c, not three: `not-noded` is now drawn as a
    mesh and has its own class below. The two that remain fail at two different
    depths -- `degenerate` never reaches `triangulate` at all (the PSLG
    validator rejects it), `hole-in-hole` is a backend refusal
    (`InvalidTopology`, since `4482649`). A third depth is now reachable and
    has no fixture of its own: the NODER's refusal, exercised through
    `--snap-spacing` in `TestSnapSpacing`.
    """

    def test_the_degenerate_fixture_is_refused_before_triangulation(self) -> None:
        # The design says this fixture shows "the DegenerateGeometry
        # presentation". It does not: the PSLG validator rejects an all-collinear
        # outer ring, so `triangulate` is never called and the words come from a
        # `PslgDiagnostic`. Asserted here so the discrepancy is recorded in the
        # suite rather than only in a report.
        drawable, words = core_verdict("degenerate")
        assert drawable is False
        assert words != []

    def test_the_degenerate_picture_carries_the_validators_own_words(
        self, tmp_path: Path
    ) -> None:
        _, words = core_verdict("degenerate")
        target = tmp_path / "degenerate.svg"
        assert invoke("degenerate", "--out", str(target)).exit_code == 0
        header = group_text(written(target), "header")
        assert any(word in header for word in words), f"{words} not in {header!r}"

    def test_the_hole_in_hole_fixture_really_fails_in_the_backend(self) -> None:
        # The third refusal, and the newest: `4482649` turned this fixture from
        # a drawable mesh into an `InvalidTopology` refusal. The same probe the
        # other two carry, for the same reason -- if it ever triangulated again
        # the blank-page case below would silently stop testing a failure.
        drawable, words = core_verdict("hole-in-hole")
        assert drawable is False
        assert words != ["Ok"]

    @pytest.mark.parametrize("name", ["degenerate", "hole-in-hole"])
    def test_a_failed_fixture_is_never_a_blank_page(
        self, name: str, tmp_path: Path
    ) -> None:
        # The input drawn alone, in its role colours. This is the guard that
        # increment 4's risk 5 asks the viewer to be, and it is worth nothing if
        # the page comes out empty.
        target = tmp_path / f"{name}.svg"
        assert invoke(name, "--out", str(target)).exit_code == 0
        document = written(target)
        assert elements(document, "edges", "line") != []
        assert any(
            "role-" in (line.get("class") or "")
            for line in elements(document, "edges", "line")
        ), "the input is drawn without its role colours"


class TestNodedCrossing:
    """`rasputin draw not-noded`, which is the increment's entire product.

    The fixture keeps its name because the name describes the INPUT, which is
    still not noded; its value is being a regression fixture with two states,
    and the picture changing between two commits under one name is the thing
    worth having.

    Every assertion here would be satisfied by a wrong picture if it only
    checked the status, which is why none of them only checks the status.
    """

    def test_the_crossing_now_triangulates(self) -> None:
        # The probe the rest of the class rests on, and the exact inversion of
        # the case this suite carried until 5c: the same fixture, the same
        # helper, the opposite answer.
        drawable, words = core_verdict("not-noded")
        assert drawable is True
        assert words == ["Ok"]

    def test_the_picture_is_a_mesh(self, tmp_path: Path) -> None:
        target = tmp_path / "not-noded.svg"
        result = invoke("not-noded", "--out", str(target))
        assert result.exit_code == 0, plain(result.output)
        assert elements(written(target), "triangles", "polygon") != []

    def test_the_picture_carries_no_finding(self, tmp_path: Path) -> None:
        # THE assertion of the increment. A `NodedPslg`'s roles are `ChainRole`
        # values; `CLOSED_ROLES` was a tuple of strings, and
        # `ChainRole.Outer == "outer"` is False. Pass the noded graph without
        # changing the vocabulary and no ring closes: one
        # `MASKED_EDGE_WITHOUT_CHAIN` per ring, drawn in the alarm colour over
        # a mesh the engine is perfectly happy with. A test that checked only
        # `CdtStatus.Ok` passes against exactly that picture.
        target = tmp_path / "not-noded.svg"
        assert invoke("not-noded", "--out", str(target)).exit_code == 0
        document = written(target)
        alarms = [
            line
            for line in elements(document, "edges", "line")
            if "finding" in (line.get("class") or "")
        ]
        assert alarms == []
        assert re.search(r"(?i)\bfindings\b\D{0,3}0\b", group_text(document, "header"))

    def test_the_road_and_the_river_are_drawn_apart(self, tmp_path: Path) -> None:
        # The user's sentence, in two colours: "a road crossing a river". Both
        # property strokes must reach the picture, which needs `_PRECEDENCE` to
        # name both and `svg.py` to have a rule for each -- without the second,
        # a road edge draws identically to the row above it and the legend
        # gains a row a reader cannot tell apart.
        target = tmp_path / "not-noded.svg"
        assert invoke("not-noded", "--out", str(target)).exit_code == 0
        tokens = {
            token
            for line in elements(written(target), "edges", "line")
            for token in (line.get("class") or "").split()
        }
        assert {"river", "road"} <= tokens

    def test_water_is_drawn_over_infrastructure(self) -> None:
        # The ordering rule `cli._PRECEDENCE` states (`grep -n _PRECEDENCE
        # src_python/tin_engine/cli.py`), named by FEATURE rather than
        # derived from the vocabulary's numbering: an edge carrying both bits
        # is drawn with exactly one token, and it is the river's.
        assert [stroke.token for stroke in cli.PROPERTY_STROKES] == ["river", "road"]

    def test_the_band_says_the_backend_was_happy(self, tmp_path: Path) -> None:
        target = tmp_path / "not-noded.svg"
        assert invoke("not-noded", "--out", str(target)).exit_code == 0
        assert "Ok" in group_text(written(target), "header")


class TestAttempt:
    """The source and its closed-role vocabulary are ONE value.

    `05c-noder-wiring.md`, "The scene's source": a source without its roles
    must be unrepresentable, because the failure it causes is invisible to
    every status assertion in this file. The picture test above is the
    consequence; this is the mechanism, asserted directly so that a failure
    says "the vocabulary does not match the source" rather than "the picture
    has alarms on it".
    """

    def attempt(self, name: str, spacing: float = DEFAULT_SNAP_SPACING) -> Any:
        # Keyword arguments, so this pins the two knobs by name and not the
        # positional order of a private helper.
        return cli._triangulated(gallery()[name], delaunay=True, spacing=spacing)

    def test_it_is_frozen(self) -> None:
        import dataclasses

        assert dataclasses.is_dataclass(cli.Attempt)
        assert cli.Attempt.__dataclass_params__.frozen is True

    def test_it_carries_the_source_and_the_vocabulary_together(self) -> None:
        import dataclasses

        names = {f.name for f in dataclasses.fields(cli.Attempt)}
        assert {"source", "closed_roles", "mesh", "ok", "status", "message"} <= names

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_ring_of_the_source_speaks_the_vocabulary(self, name: str) -> None:
        # The oracle is the FIXTURE's own declaration -- which chains it called
        # `outer` and `hole` -- and never the source's roles, which is the
        # value under test. Chain order is preserved through `build_pslg` and
        # `node`, so position is the join.
        fixture = gallery()[name]
        attempt = self.attempt(name)
        for c, declared in enumerate(fixture.chains):
            if declared.role not in ("outer", "hole"):
                continue
            role = attempt.source.chains[c].role
            assert role in attempt.closed_roles, (
                f"{name} chain {c} is a ring with role {role!r}, "
                f"which is not in {attempt.closed_roles!r}"
            )

    def test_a_noded_fixture_is_drawn_from_the_noded_graph(self) -> None:
        # "The most-processed graph that exists", with no branch on the CDT's
        # status: when the noder succeeded, its output is the truer picture,
        # because it shows the splits the backend was actually handed.
        import tin_engine._core as core

        attempt = self.attempt("not-noded")
        assert isinstance(attempt.source, core.NodedPslg)
        # And it really did gain a node -- otherwise the claim is vacuous.
        assert len(attempt.source.vertices) > len(np.asarray(gallery()["not-noded"].vertices))

    def test_a_rejected_fixture_is_drawn_from_the_fixture(self) -> None:
        # `degenerate` has no `Pslg` at all, so the fixture is the only
        # `PslgLike` that exists -- and its roles are strings, which is exactly
        # why the vocabulary cannot be a module constant.
        attempt = self.attempt("degenerate")
        assert attempt.source is gallery()["degenerate"]
        assert attempt.closed_roles == ("outer", "hole")


class TestSnapSpacing:
    """Two refusal channels, neither removable in favour of the other.

    `typer` refuses a non-finite or non-positive `--snap-spacing` as a USAGE
    error and exits 2 without drawing anything -- the same line `_destination`
    already draws between a bad `--out` and a refused fixture. The engine
    refuses it as a library precondition, `NodeStatus::InvalidSnapSpacing`, for
    callers that are not this CLI; that half is pinned in `test_core_noding.py`.
    Delete the CLI guard and `--snap-spacing -1` silently produces a picture of
    a failure that is not about the terrain.
    """

    @pytest.mark.parametrize("value", ["0", "-1", "nan", "-inf"])
    def test_a_non_positive_or_non_finite_spacing_is_a_usage_error(
        self, value: str, tmp_path: Path
    ) -> None:
        target = tmp_path / "refused.svg"
        result = invoke("catchment", f"--snap-spacing={value}", "--out", str(target))
        assert result.exit_code == 2, plain(result.output)
        message = plain(result.output)
        assert "snap-spacing" in message, f"the refusal does not name the option: {message!r}"
        assert not target.exists(), "a usage error drew a picture"

    def test_the_option_is_documented_with_its_default(self) -> None:
        result = runner.invoke(app, ["draw", "--help"])
        assert result.exit_code == 0
        assert "--snap-spacing" in plain(result.output)

    def test_the_default_is_a_named_constant_with_its_reasoning(self) -> None:
        # Not a literal inside a `typer.Option`: three of the nine noder
        # statuses point at this number and two point in opposite directions,
        # so it is a policy choice the composition root makes on the record.
        assert cli.DEFAULT_SNAP_SPACING == DEFAULT_SNAP_SPACING
        source = (
            Path(cli.__file__).read_text(encoding="utf-8")
            if cli.__file__
            else ""
        )
        assert "DEFAULT_SNAP_SPACING" in source

    def test_a_coarser_spacing_still_draws(self, tmp_path: Path) -> None:
        target = tmp_path / "coarse.svg"
        result = invoke("not-noded", "--snap-spacing=0.5", "--out", str(target))
        assert result.exit_code == 0, plain(result.output)
        assert elements(written(target), "triangles", "polygon") != []

    def test_the_spacing_reaches_the_engine(self, tmp_path: Path) -> None:
        # Nothing above holds the option to `node()`: a `cli.py` that accepted
        # `--snap-spacing` and passed the constant would satisfy every one of
        # them. At 1e-12 a UTM easting overflows the lattice, so the engine
        # must answer `CoordinateOutOfRange` -- a value only the real argument
        # can produce.
        drawable, words = core_verdict("catchment", OVERFLOWING_SNAP_SPACING)
        assert drawable is False
        assert words == ["CoordinateOutOfRange"]

        target = tmp_path / "overflow.svg"
        result = invoke(
            "catchment", f"--snap-spacing={OVERFLOWING_SNAP_SPACING}", "--out", str(target)
        )
        assert result.exit_code == 0, plain(result.output)
        assert "CoordinateOutOfRange" in group_text(written(target), "header")


class TestNoderRefusalPresentation:
    """What a user sees when the NODER refuses: a diagnosis about the data,
    with the spacing as the lever.

    The input drawn alone in role colours, findings suppressed, exit 0, and a
    band carrying `describe(status)` followed by the engine's own message --
    the same two-part join `_triangulated` already builds for a `CdtStatus`.
    A drawn failure is still a drawn picture, so the exit code stays 0 and a
    caller who wants the verdict reads the band.
    """

    def test_a_refused_noding_draws_the_input_in_role_colours(
        self, tmp_path: Path
    ) -> None:
        target = tmp_path / "refused.svg"
        result = invoke(
            "catchment", f"--snap-spacing={OVERFLOWING_SNAP_SPACING}", "--out", str(target)
        )
        assert result.exit_code == 0, plain(result.output)
        document = written(target)
        lines = elements(document, "edges", "line")
        assert lines != [], "the refusal drew a blank page"
        assert any("role-" in (line.get("class") or "") for line in lines)

    def test_a_refused_noding_draws_no_triangle_and_no_finding(
        self, tmp_path: Path
    ) -> None:
        target = tmp_path / "refused.svg"
        assert invoke(
            "catchment", f"--snap-spacing={OVERFLOWING_SNAP_SPACING}", "--out", str(target)
        ).exit_code == 0
        document = written(target)
        assert elements(document, "triangles", "polygon") == []
        alarms = [
            line
            for line in elements(document, "edges", "line")
            if "finding" in (line.get("class") or "")
        ]
        assert alarms == [], "findings are not suppressed on a refusal"

    def test_the_band_carries_the_engines_own_words(self, tmp_path: Path) -> None:
        import tin_engine._core as core

        target = tmp_path / "refused.svg"
        assert invoke(
            "catchment", f"--snap-spacing={OVERFLOWING_SNAP_SPACING}", "--out", str(target)
        ).exit_code == 0
        header = group_text(written(target), "header")
        sentence = core.describe(core.NodeStatus.CoordinateOutOfRange)
        # The prose is the engine's and is not pinned here; that it CROSSED is.
        assert sentence.split()[0] in header
        assert "CoordinateOutOfRange" in header


class TestCornerGrazePresentation:
    """`NotConverged`, which is the first place a human meets the corner graze.

    No gallery fixture reaches it -- the period-two orbit needs input already
    rounded to the grid's own resolution -- so the outcome is substituted at
    the seam `cli.py` calls. That is the only thing faked: the status and the
    message are the engine's real ones, and everything downstream of
    `_triangulated` is the shipped code.

    It is a diagnosis about the data with the spacing as the lever. It is NOT
    a crash, and it is not a promise that a bigger cap helps -- for the graze
    it provably cannot, the orbit having period two.
    """

    class _Refusal:
        """A `NodeOutcome` the binding cannot construct from Python."""

        def __init__(self, status: Any, message: str) -> None:
            self.pslg = None
            self.status = status
            self.message = message

        def ok(self) -> bool:
            return False

    @pytest.fixture
    def not_converged(self, monkeypatch: pytest.MonkeyPatch) -> str:
        import tin_engine._core as core

        message = (
            "the constraint set did not settle in 4 rounds at spacing 0.001: "
            "a constraint edge still meets a node's cell"
        )
        monkeypatch.setattr(
            cli,
            "node",
            lambda *_args, **_kw: self._Refusal(core.NodeStatus.NotConverged, message),
        )
        return message

    def test_it_is_drawn_rather_than_raised(self, not_converged: str, tmp_path: Path) -> None:
        target = tmp_path / "graze.svg"
        result = invoke("catchment", "--out", str(target))
        assert result.exit_code == 0, plain(result.output)
        assert elements(written(target), "edges", "line") != []
        assert elements(written(target), "triangles", "polygon") == []

    def test_the_band_carries_the_status_the_prose_and_the_message(
        self, not_converged: str, tmp_path: Path
    ) -> None:
        import tin_engine._core as core

        target = tmp_path / "graze.svg"
        assert invoke("catchment", "--out", str(target)).exit_code == 0
        header = group_text(written(target), "header")
        assert "NotConverged" in header
        # The whole sentence, not its first word -- which is "the", and passes
        # on a band printing one article. Derived from `describe` at run time
        # rather than written out, so a reworded row keeps this green: what is
        # pinned is that the band carries the engine's words verbatim, which is
        # the passthrough claim. Normalised on both sides because `group_text`
        # collapses the SVG's whitespace and the row may be laid out with more.
        sentence = " ".join(core.describe(core.NodeStatus.NotConverged).split())
        assert sentence in header
        assert "did not settle in 4 rounds" in header

    def test_the_cli_adds_no_lever_of_its_own(self, not_converged: str, tmp_path: Path) -> None:
        # The band is the engine's two sentences and nothing this module wrote.
        # A CLI-authored sentence here would be a second authority for a C++
        # fact, drifting the first time a row is reworded -- and the one lever
        # it would be tempted to name is the cap, which the CLI does not turn.
        target = tmp_path / "graze.svg"
        assert invoke("catchment", "--out", str(target)).exit_code == 0
        header = group_text(written(target), "header")
        assert "--max-rounds" not in header
        assert "max_rounds" not in header


class TestMaxRoundsIsNotExposed:
    """The cap is a bound on a loop, not a parameter of the answer: at any
    value the outcome is the same mesh or `NotConverged`, never a different
    mesh. For the corner graze `Ok` is unreachable at ANY cap, so a knob whose
    visible effect on the failure a user is most likely to meet is "the same
    refusal, slower" is worse than not having it."""

    def test_the_option_does_not_exist(self, tmp_path: Path) -> None:
        target = tmp_path / "x.svg"
        result = invoke("not-noded", "--max-rounds", "16", "--out", str(target))
        assert result.exit_code == 2
        assert "max-rounds" in plain(result.output)
        assert not target.exists()

    def test_the_help_offers_the_spacing_and_not_the_cap(self) -> None:
        text = plain(runner.invoke(app, ["draw", "--help"]).output)
        assert "--snap-spacing" in text
        assert "max-rounds" not in text


class TestTheDefaultSpacingSuitsTheGallery:
    """A default that refuses a shipped fixture is a broken default, and the
    answer to that is a different default rather than a per-fixture knob --
    `Fixture` has no spacing field, and adding one would put a policy in the
    module that declares it holds none. Measured over the whole gallery rather
    than argued."""

    @pytest.mark.parametrize("name", GALLERY_NAMES)
    def test_every_valid_fixture_nodes_at_the_default(self, name: str) -> None:
        import tin_engine._core as core

        fixture = gallery()[name]
        chains = [
            ([int(i) for i in fixture.indices_of(c)], ROLES[chain.role], int(chain.properties))
            for c, chain in enumerate(fixture.chains)
        ]
        result = core.build_pslg(np.asarray(fixture.vertices), chains)
        if not result.ok:
            # `degenerate` never reaches the noder; that is its own test.
            pytest.skip(f"{name} is refused by the PSLG validator")
        outcome = core.node(result.pslg, DEFAULT_SNAP_SPACING)
        assert outcome.status == core.NodeStatus.Ok, (
            f"{name} is refused at the default spacing: {outcome.message}"
        )
