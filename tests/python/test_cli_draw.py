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
   then `triangulate`;
3. hand the **fixture itself** to `build_scene` as the `PslgLike`, with the
   mesh and the `ok` bit from the backend, and `closed_roles=("outer", "hole")`;
4. pass the status name, `describe(status)` and the backend message to
   `render_svg` as text;
5. resolve and validate the output path, and write the bytes.

Step 3 is the one the design does not state, and `TestFailurePresentation`
below is why it has to be that way: `06-cdt-viewer.md` says the `degenerate`
fixture shows "the `DegenerateGeometry` presentation", but an all-collinear ring
never reaches `triangulate` -- the PSLG validator rejects it first, so there is
no `Pslg` object at all and nothing for `build_scene` to be given unless the
fixture is itself a `PslgLike`. Measured, on this tree:

    build_pslg(collinear, [(range(4), ChainRole.Outer, False)]).diagnostics
    -> [PslgError.DegenerateRing, 'chain 0 is declared Outer but every vertex
        is collinear']

The presentation the design wants is still seen -- the input drawn alone with
the engine's own words in the header band -- but the words come from a
`PslgDiagnostic` rather than from a `CdtStatus`, and the assertions below take
them from the engine at run time rather than hard-coding either.

Path handling is tested as a hostile boundary, per the python skill: traversal
out of an explicitly provided parent, and a symlink that resolves back inside it.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from typer.testing import CliRunner

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


def core_verdict(name: str) -> tuple[bool, list[str]]:
    """What the engine actually says about a fixture, computed here.

    The oracle for a passthrough claim is the engine's own words, so this
    re-runs the composition the command performs rather than hard-coding a
    status. It is deliberately *not* an independent derivation -- there is
    nothing independent to derive; the claim under test is only that whatever
    the engine said reaches the picture.

    Returns `(drawable, words)`: whether a mesh was produced at all, and the
    strings the header band must carry.
    """
    import tin_engine._core as core

    fixture = gallery()[name]
    roles = {
        "outer": core.ChainRole.Outer,
        "hole": core.ChainRole.Hole,
        "breakline": core.ChainRole.Breakline,
    }
    chains = [
        ([int(i) for i in fixture.indices_of(c)], roles[chain.role], chain.is_river)
        for c, chain in enumerate(fixture.chains)
    ]
    result = core.build_pslg(np.asarray(fixture.vertices), chains)
    if not result.ok:
        return False, [d.error.name for d in result.diagnostics]
    assert result.pslg is not None
    outcome = core.triangulate(result.pslg)
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

    They fail at three different depths, which is why three: `degenerate` never
    reaches `triangulate` at all (the PSLG validator rejects it), `not-noded`
    and `hole-in-hole` are backend refusals -- `NotNoded` and, since `4482649`,
    `InvalidTopology`.
    """

    def test_the_not_noded_fixture_really_fails_in_the_backend(self) -> None:
        # The probe that keeps the assertion below from going vacuous: if the
        # fixture ever triangulated cleanly, the non-`Ok` presentation would not
        # be exercised by anything.
        drawable, words = core_verdict("not-noded")
        assert drawable is False
        assert words != ["Ok"]

    def test_the_not_noded_picture_carries_the_backends_own_words(
        self, tmp_path: Path
    ) -> None:
        _, words = core_verdict("not-noded")
        target = tmp_path / "not-noded.svg"
        assert invoke("not-noded", "--out", str(target)).exit_code == 0
        header = group_text(written(target), "header")
        assert any(word in header for word in words), f"{words} not in {header!r}"

    def test_the_not_noded_picture_draws_no_triangle(self, tmp_path: Path) -> None:
        target = tmp_path / "not-noded.svg"
        assert invoke("not-noded", "--out", str(target)).exit_code == 0
        assert elements(written(target), "triangles", "polygon") == []

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

    @pytest.mark.parametrize("name", ["not-noded", "degenerate", "hole-in-hole"])
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
