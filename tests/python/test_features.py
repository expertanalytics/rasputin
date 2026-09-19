"""`tin_engine.features`: the edge-property vocabulary at the boundary.

Committed red at `dc8994b`: `src_python/tin_engine/features.py` did not exist
yet, so the import below raised `ModuleNotFoundError: No module named
'tin_engine.features'` while this file was collected, which was the intended
red for increment 7's commit 1. The module landed in `e5dee4e` and the suite
has been green since.

What outlives that step is the choice the red step made: the import is
deliberately module-scope rather than the per-test import
`test_viz_protocols.py` uses. There the point was that a missing `viz/`
package should name itself in one test instead of aborting collection of a
suite that also tests the bindings; here the whole file is about one module
that either exists or does not, so aborting collection is the honest failure
and remains the right shape now that the module exists.

**This is the only invariant-critical suite in increment 7 and it carries the
increment's one mutation round** (`docs/increments/07-edge-properties.md`, "What
is worth testing"). The reason is not that the code is hard: it is that this is
the only place in the system where a wrong answer is *silent*. The C++ core
merges opaque bits and can never see that two producers disagree about which bit
means "river" -- risk 20. Everything else in this increment fails at a compiler.
The seven mutants the round must kill are named in the design; each is named
again at the test that kills it, in the form ``MUTANT n``, so that deleting a
test makes a named claim disappear rather than a line count drop.

Two rules here are architecture rather than hygiene:

* **`features.py` imports nothing first-party and never imports `_core`.** The
  reason is the one the module's own docstring gives: it is constructible and
  testable with no extension in the process. It is emphatically *not* that
  `viz/` is allowed to depend on it -- `viz/` may not import this module at
  all, and `test_viz_svg.py::TestModuleIsolation` exists to deny exactly that
  permission, pinning `style.py` and `fixtures.py` to zero first-party imports.
  `cli.py`, the composition root, is the only module that imports both this
  one and `viz.style` (`cli.py:39,42`). The rule is checked here by parsing the
  source, not by inspecting `sys.modules`, because `tin_engine/__init__.py`
  imports `_core` itself -- so an import-time check
  would be asserting something about the package rather than about this module.
  Same reasoning, same mechanism, as `test_viz_protocols.py`.
* **The name pattern keeps a vocabulary name usable as a CSS class token, and
  is defence in depth. It is not what closes the injection hole.**
  `viz/svg.py`'s `_edge_classes` (`svg.py:166`) joins
  `style.PropertyStroke.token` values, which `_edges` then interpolates into a
  ``class="..."`` attribute **unescaped** (`svg.py:225`) -- only `_text`
  escapes -- and never an `EdgeProperty` name, because `viz/` cannot import
  this module. The only bridge is `cli.py:84`, which builds a `PropertyStroke`
  from a name and so re-validates through `style.py:51`'s own, deliberately
  wider `^[a-z][a-z0-9_-]*$` -- a CSS class may carry a hyphen and a feature
  name may not. A quote-carrying name dies there whatever
  `^[a-z][a-z0-9_]*$` says. The attribute is closed by
  `PropertyStroke.token`, and the suite over that boundary is
  `test_viz_svg.py::TestPropertyStrokes`.

  What tokens exist today is not an enum, either: `cli.py:79`'s
  `_PRECEDENCE = ("river",)` is a literal tuple of `str`, and there is no enum
  anywhere on this path. A vocabulary read from a configuration file is still
  untrusted input, and a second, narrower gate on it is worth its lines -- so
  every case below stays.
"""

from __future__ import annotations

import ast
import os
import subprocess
import sys
from pathlib import Path

import pydantic
import pytest

from tin_engine.features import DEFAULT_VOCABULARY, EdgeProperty, EdgeVocabulary

REPO_ROOT = Path(__file__).resolve().parents[2]
FEATURES_SOURCE = REPO_ROOT / "src_python" / "tin_engine" / "features.py"

# The seven linear features `DEFAULT_VOCABULARY` ships with, in bit order. The
# tuple is written out rather than derived from the constant, because a test
# that derives its expectation from the thing under test asserts nothing.
DEFAULT_ROWS = (
    ("river", 0),
    ("road", 1),
    ("railway", 2),
    ("coastline", 3),
    ("contour", 4),
    ("wall", 5),
    ("ditch", 6),
)


def vocabulary(*rows: tuple[str, int]) -> EdgeVocabulary:
    """An `EdgeVocabulary` over ``(name, bit)`` pairs, for readability."""
    return EdgeVocabulary(properties=tuple(EdgeProperty(name=n, bit=b) for n, b in rows))


# ---------------------------------------------------------------------------
# EdgeProperty: the bit range and the name pattern.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bit", [0, 1, 15, 30, 31])
def test_edge_property_accepts_every_bit_the_cpp_word_has(bit: int) -> None:
    assert EdgeProperty(name="river", bit=bit).bit == bit


@pytest.mark.parametrize("bit", [-1, -32, 32, 33, 64])
def test_edge_property_rejects_a_bit_outside_the_32_bit_ceiling(bit: int) -> None:
    # `EdgeProperties::kMaxProperties == 32` is a C++ constant no Python caller
    # can see. This is that ceiling restated at the only layer that meets
    # untrusted input, and 32 itself is the boundary that an `le=32` typo moves.
    with pytest.raises(pydantic.ValidationError):
        EdgeProperty(name="river", bit=bit)


@pytest.mark.parametrize("name", ["river", "road", "r", "a0", "land_cover_2", "x_"])
def test_edge_property_accepts_a_lowercase_css_safe_token(name: str) -> None:
    assert EdgeProperty(name=name, bit=0).name == name


@pytest.mark.parametrize(
    "name",
    [
        pytest.param('ri"ver', id="double-quote-ends-the-class-attribute"),
        pytest.param("ri'ver", id="single-quote"),
        pytest.param("ri ver", id="space-splits-into-two-class-tokens"),
        pytest.param('river" onload="alert(1)', id="attribute-injection"),
        pytest.param("<script>", id="angle-brackets"),
        pytest.param("river\n", id="trailing-newline"),
        pytest.param("riv\ner", id="embedded-newline"),
        pytest.param("River", id="uppercase"),
        pytest.param("1river", id="leading-digit"),
        pytest.param("_river", id="leading-underscore"),
        pytest.param("riv-er", id="hyphen"),
        pytest.param("", id="empty"),
        pytest.param("rivér", id="non-ascii"),
    ],
)
def test_edge_property_refuses_a_name_that_is_not_a_css_class_token(name: str) -> None:
    """MUTANT 5: the name pattern relaxed to any ``str``.

    The design names two of these as the minimum -- a quote and a space -- and
    both stay, under the corrected reason: a quote or a space makes a name
    unusable as a CSS class token, which is what this pattern is for. It is not
    what protects the ``class="..."`` attribute. `svg.py:_edge_classes`
    (`svg.py:166`) joins `style.PropertyStroke.token` values and `_edges`
    (`svg.py:225`) interpolates the result unescaped, never a name from this
    module, and `cli.py:84` re-validates any name through that model before it
    can reach the attribute -- see `test_viz_svg.py::TestPropertyStrokes`, which
    refuses the same shapes at the boundary that is load-bearing. Relaxing this
    pattern to any ``str`` therefore does not open the injection hole; it lets a
    name into the vocabulary that the drawing layer would then refuse, which is
    a failure moved to a worse place rather than prevented.

    ``river\\n`` is the one that is not about SVG. Python's ``re`` matches ``$``
    *before* a trailing newline, so a validator hand-written as
    ``re.match(r"^[a-z][a-z0-9_]*$", name)`` accepts it while
    ``Field(pattern=...)`` -- pydantic's Rust engine, where ``$`` is
    end-of-haystack -- refuses it. Measured against pydantic 2.13.5 before this
    case was written down.
    """
    with pytest.raises(pydantic.ValidationError):
        EdgeProperty(name=name, bit=0)


def test_edge_property_is_frozen() -> None:
    prop = EdgeProperty(name="river", bit=0)

    with pytest.raises(pydantic.ValidationError):
        prop.bit = 1  # type: ignore[misc]


def test_edge_property_forbids_an_unknown_field() -> None:
    # extra="forbid", so a typo in a config file is a refusal and not a silent
    # default -- the failure mode being a vocabulary that looks configured and
    # is not.
    with pytest.raises(pydantic.ValidationError):
        EdgeProperty(name="river", bit=0, colour="blue")  # type: ignore[call-arg]


# ---------------------------------------------------------------------------
# EdgeVocabulary: construction-time validation.
# ---------------------------------------------------------------------------


def test_an_empty_vocabulary_is_legal() -> None:
    # A vocabulary with no properties is a pipeline that classifies nothing,
    # which is the default state of the system and not an error.
    empty = EdgeVocabulary(properties=())

    assert empty.properties == ()
    assert empty.mask() == 0
    assert empty.names(0) == ()


def test_a_single_property_vocabulary_is_legal() -> None:
    one = vocabulary(("river", 0))

    assert len(one.properties) == 1
    assert one.mask("river") == 1


def test_vocabulary_is_frozen() -> None:
    vocab = vocabulary(("river", 0))

    with pytest.raises(pydantic.ValidationError):
        vocab.properties = ()  # type: ignore[misc]


def test_vocabulary_forbids_an_unknown_field() -> None:
    with pytest.raises(pydantic.ValidationError):
        EdgeVocabulary(properties=(), crs="EPSG:25833")  # type: ignore[call-arg]


def test_two_properties_on_one_bit_are_refused_at_construction() -> None:
    """MUTANT 2: the bit-uniqueness validator dropped.

    Two names on one bit is risk 20 in its cheapest and most detectable form,
    and construction is the only moment it is detectable for free. The refusal
    names the bit, per this project's "refuses with the value that did not fit"
    rule.
    """
    with pytest.raises(pydantic.ValidationError) as excinfo:
        vocabulary(("river", 3), ("road", 3))

    assert "3" in str(excinfo.value)


def test_two_properties_with_one_name_are_refused_at_construction() -> None:
    """MUTANT 3: the name-uniqueness validator dropped.

    Without it, `mask("river")` depends on which duplicate the lookup reaches
    first -- an order-dependent answer from a frozen model, which is the worst
    shape a defect can take because it reproduces.
    """
    with pytest.raises(pydantic.ValidationError) as excinfo:
        vocabulary(("river", 0), ("river", 4))

    assert "river" in str(excinfo.value)


def test_duplicate_bits_are_refused_even_when_not_adjacent() -> None:
    # A validator that compares each property only with its neighbour passes
    # the adjacent case above and fails here.
    with pytest.raises(pydantic.ValidationError):
        vocabulary(("river", 0), ("road", 1), ("railway", 2), ("ditch", 0))


def test_duplicate_names_are_refused_even_when_not_adjacent() -> None:
    with pytest.raises(pydantic.ValidationError):
        vocabulary(("river", 0), ("road", 1), ("railway", 2), ("river", 9))


def test_distinct_names_on_distinct_bits_are_accepted_in_any_order() -> None:
    # The uniqueness validators must not have become "sorted by bit" validators:
    # a vocabulary is free to list its properties in whatever order its author
    # wrote them.
    descending = vocabulary(("ditch", 6), ("road", 1), ("river", 0))

    assert descending.mask("river", "ditch") == 0b100_0001


# ---------------------------------------------------------------------------
# mask(): names to a word.
# ---------------------------------------------------------------------------


def test_mask_shifts_by_the_declared_bit_and_not_by_the_tuple_index() -> None:
    """MUTANT 1: `mask()` shifting by the tuple index instead of by `.bit`.

    The cheapest killer the design names: a one-property vocabulary whose single
    property is at bit 3, where index-shifting yields 1 and bit-shifting yields
    8.
    """
    assert vocabulary(("road", 3)).mask("road") == 8


def test_mask_shifts_by_the_declared_bit_for_every_member() -> None:
    # The same mutant, over a vocabulary whose bits are a permutation of a
    # contiguous range: index-shifting agrees with bit-shifting on the set but
    # not on the individual answers.
    vocab = vocabulary(("river", 2), ("road", 0), ("railway", 1))

    assert vocab.mask("river") == 0b100
    assert vocab.mask("road") == 0b001
    assert vocab.mask("railway") == 0b010


def test_mask_of_no_names_is_the_empty_set() -> None:
    assert DEFAULT_VOCABULARY.mask() == 0


def test_mask_unions_its_arguments() -> None:
    assert DEFAULT_VOCABULARY.mask("river", "railway") == 0b101


def test_mask_is_commutative_and_idempotent_over_its_arguments() -> None:
    # The same three laws `EdgeProperties::operator|` is held to in C++, at the
    # layer that produces the word the C++ reduces. A `mask()` that summed or
    # xor-ed instead of or-ing passes the disjoint case above and fails here.
    vocab = DEFAULT_VOCABULARY

    assert vocab.mask("river", "road") == vocab.mask("road", "river")
    assert vocab.mask("river", "river") == vocab.mask("river")
    assert vocab.mask("river", "road", "river") == vocab.mask("river", "road")


def test_mask_of_the_top_bit_does_not_overflow_into_a_negative() -> None:
    # Python ints are unbounded, so 1 << 31 is 2147483648 here and 0x80000000 in
    # C++; a validator written with a signed 32-bit assumption is the thing that
    # would disagree. The mask crosses to pybind as an int and must be positive.
    assert vocabulary(("wall", 31)).mask("wall") == 2**31


def test_mask_raises_on_a_name_the_vocabulary_does_not_have() -> None:
    """MUTANT 7: `mask()` accepting an unknown name and returning 0.

    A silent zero is a constraint that quietly loses every property it was
    asked for, which is the exact failure this increment exists to prevent. The
    raise names the name, so the caller learns which of several it got wrong.
    """
    with pytest.raises(ValueError, match="glacier"):
        DEFAULT_VOCABULARY.mask("glacier")


def test_mask_raises_even_when_one_of_several_names_is_known() -> None:
    # The partial case is the dangerous one: returning the known bits and
    # dropping the unknown name is indistinguishable from success in any
    # picture drawn afterwards.
    with pytest.raises(ValueError, match="glacier"):
        DEFAULT_VOCABULARY.mask("river", "glacier")


def test_mask_raises_on_a_name_from_a_different_vocabulary() -> None:
    # Risk 20's shape at the API: a name that is perfectly valid somewhere else
    # is still not in the vocabulary in hand.
    with pytest.raises(ValueError, match="road"):
        vocabulary(("river", 0)).mask("road")


# ---------------------------------------------------------------------------
# names(): a word back to names.
# ---------------------------------------------------------------------------


def test_names_of_the_empty_mask_is_empty() -> None:
    assert DEFAULT_VOCABULARY.names(0) == ()


def test_names_returns_the_named_bits_in_ascending_bit_order() -> None:
    # The design fixes the return type as a tuple and does not fix its order.
    # An unspecified order in a tuple is a latent nondeterminism, so this suite
    # pins ascending bit order: it is the only order that is a property of the
    # vocabulary rather than of the argument or of the declaration sequence.
    vocab = vocabulary(("ditch", 6), ("road", 1), ("river", 0))

    assert vocab.names(vocab.mask("ditch", "river")) == ("river", "ditch")
    assert vocab.names(vocab.mask("river", "ditch")) == ("river", "ditch")


def test_names_round_trips_mask_for_every_subset_of_the_default_vocabulary() -> None:
    vocab = DEFAULT_VOCABULARY
    every = tuple(name for name, _ in DEFAULT_ROWS)

    assert vocab.names(vocab.mask(*every)) == every
    for name in every:
        assert vocab.names(vocab.mask(name)) == (name,)


def test_names_raises_on_a_bit_no_property_names() -> None:
    """MUTANT 4: `names(mask)` silently dropping an unnamed bit.

    The mutant the ruling exists for, and the one whose survival is
    indistinguishable from correct behaviour in any picture. Reaching this state
    means the mask came from a different vocabulary than the one in hand, which
    is risk 20's failure, caught. The raise names the bit.
    """
    vocab = vocabulary(("river", 0), ("road", 1))

    with pytest.raises(ValueError, match="5"):
        vocab.names(1 << 5)


def test_names_raises_even_when_some_bits_are_named() -> None:
    # Dropping is the silent loss; returning ("river",) for a mask that also
    # carries an unnamed bit is precisely dropping.
    vocab = vocabulary(("river", 0), ("road", 1))

    with pytest.raises(ValueError, match="7"):
        vocab.names(vocab.mask("river") | (1 << 7))


def test_names_raises_on_a_bit_at_or_above_the_ceiling() -> None:
    # 32 is outside the C++ word entirely, so no vocabulary can ever name it.
    # It must refuse rather than loop to exhaustion or truncate.
    with pytest.raises(ValueError):
        DEFAULT_VOCABULARY.names(1 << 32)


def test_names_raises_on_a_negative_mask() -> None:
    # A negative int has no finite bit set under Python's two's-complement view
    # and cannot have come from `mask()`. Refusing is the only answer that is
    # not a guess; the alternative -- masking to 32 bits -- would silently
    # invent a property set.
    with pytest.raises(ValueError):
        DEFAULT_VOCABULARY.names(-1)


def test_names_and_mask_are_methods_rather_than_module_functions() -> None:
    # Mitigation 1: a bare mask is unreachable through the API. You cannot
    # obtain one without holding the vocabulary that produced it, nor interpret
    # one without holding a vocabulary you chose. If either ever becomes a
    # module-level function, the dangerous state -- a number with no units --
    # becomes the cheap one.
    import tin_engine.features as features

    assert not hasattr(features, "mask")
    assert not hasattr(features, "names")
    assert callable(EdgeVocabulary.mask)
    assert callable(EdgeVocabulary.names)


# ---------------------------------------------------------------------------
# fingerprint(): the only mechanism that detects a disagreement in assignment.
# ---------------------------------------------------------------------------


def test_fingerprint_is_a_non_empty_string() -> None:
    fingerprint = DEFAULT_VOCABULARY.fingerprint()

    assert isinstance(fingerprint, str)
    assert fingerprint != ""


def test_fingerprint_is_stable_within_a_process() -> None:
    assert DEFAULT_VOCABULARY.fingerprint() == DEFAULT_VOCABULARY.fingerprint()


def test_equal_vocabularies_fingerprint_equally() -> None:
    assert vocabulary(*DEFAULT_ROWS).fingerprint() == DEFAULT_VOCABULARY.fingerprint()


def test_fingerprint_ignores_the_order_the_properties_were_declared_in() -> None:
    # A digest over the *sorted* (bit, name) pairs. Two authors who wrote the
    # same vocabulary in different orders have the same vocabulary, and an
    # artifact written by one must be readable by the other.
    assert vocabulary(*reversed(DEFAULT_ROWS)).fingerprint() == DEFAULT_VOCABULARY.fingerprint()


def test_fingerprint_sees_a_permutation_of_the_bits() -> None:
    """MUTANT 6: `fingerprint()` hashing names only, or bits only.

    The load-bearing case. Two vocabularies that both name bit 3, differently,
    are undetectable by every other mechanism in this module -- `names()` only
    catches a *coverage* disagreement, never an *assignment* one. A digest that
    cannot see a permutation blesses exactly risk 20's failure, which is why
    this is the mutation round's reason for existing.
    """
    swapped = vocabulary(("river", 1), ("road", 0))
    straight = vocabulary(("river", 0), ("road", 1))

    assert swapped.mask("river") != straight.mask("river")
    assert swapped.fingerprint() != straight.fingerprint()


def test_fingerprint_sees_a_bit_moved_without_reordering_the_names() -> None:
    """MUTANT 6, the half a permutation does not reach.

    Measured, and the reason this test exists: a digest over the names alone,
    emitted in bit order, *does* see the swap above -- the two names come out
    in the opposite order -- so the permutation case passes a names-only
    fingerprint. Two vocabularies that agree on every name and on their order
    while disagreeing about which bit `road` occupies are the case that
    separates them, and they are risk 20 exactly: an artifact written by one is
    read by the other with every property shifted.
    """
    near = vocabulary(("river", 0), ("road", 1))
    far = vocabulary(("river", 0), ("road", 2))

    assert near.mask("road") != far.mask("road")
    assert near.fingerprint() != far.fingerprint()


def test_fingerprint_sees_a_renamed_property_on_the_same_bits() -> None:
    # The other half of mutant 6: a digest over bits only cannot see this.
    assert vocabulary(("river", 0)).fingerprint() != vocabulary(("creek", 0)).fingerprint()


def test_fingerprint_sees_a_property_added_or_removed() -> None:
    narrow = vocabulary(("river", 0))
    wide = vocabulary(("river", 0), ("road", 1))

    assert narrow.fingerprint() != wide.fingerprint()
    assert EdgeVocabulary(properties=()).fingerprint() != narrow.fingerprint()


def test_fingerprint_is_stable_across_processes() -> None:
    """It travels with the artifact, so it must survive leaving this process.

    ``hash()`` over a tuple of strings is salted by ``PYTHONHASHSEED`` and
    differs between runs, so a fingerprint built on it passes every assertion
    above and refuses every artifact ever written -- the failure appearing only
    in the field, on the second run. Two child processes with different seeds
    are what make that visible; the probe fails loudly if the child cannot
    import the module, so a pass cannot be confused with a probe that did not
    run.
    """
    program = (
        "from tin_engine.features import DEFAULT_VOCABULARY;"
        "print(DEFAULT_VOCABULARY.fingerprint())"
    )

    digests = []
    for seed in ("0", "1"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        result = subprocess.run(
            [sys.executable, "-c", program],
            capture_output=True,
            text=True,
            env=env,
            cwd=REPO_ROOT,
            check=True,
        )
        digests.append(result.stdout.strip())

    assert digests[0] != ""
    assert digests[0] == digests[1]
    assert digests[0] == DEFAULT_VOCABULARY.fingerprint()


# ---------------------------------------------------------------------------
# DEFAULT_VOCABULARY: a default, pinned, and explicitly not a schema.
# ---------------------------------------------------------------------------


def test_the_default_vocabulary_is_the_seven_linear_features_in_bit_order() -> None:
    assert tuple((p.name, p.bit) for p in DEFAULT_VOCABULARY.properties) == DEFAULT_ROWS


def test_river_is_bit_zero() -> None:
    # Risk 2, pinned as the design requires. Bit 0 is `river` so that today's
    # one-bit data and the `river` gallery fixture keep meaning what they meant.
    # **This constant is a default, not a schema**: any caller may supply its
    # own vocabulary, and the fingerprint is what makes two of them comparable.
    assert DEFAULT_VOCABULARY.mask("river") == 1


def test_the_default_vocabulary_fits_the_cpp_word_with_room_to_spare() -> None:
    # The measured argument for kMaxProperties = 32: linear features number
    # under ten, against 23 and 44 land-cover classes in `legacy/` -- and land
    # cover stays face-based.
    assert len(DEFAULT_VOCABULARY.properties) < 32
    assert all(0 <= p.bit < 32 for p in DEFAULT_VOCABULARY.properties)


def test_the_default_vocabulary_round_trips_its_own_full_mask() -> None:
    every = tuple(name for name, _ in DEFAULT_ROWS)
    full = DEFAULT_VOCABULARY.mask(*every)

    assert full == 0b111_1111
    assert DEFAULT_VOCABULARY.names(full) == every


def test_the_default_vocabulary_is_an_edge_vocabulary_and_is_frozen() -> None:
    assert isinstance(DEFAULT_VOCABULARY, EdgeVocabulary)

    with pytest.raises(pydantic.ValidationError):
        DEFAULT_VOCABULARY.properties = ()  # type: ignore[misc]


# ---------------------------------------------------------------------------
# The import firewall, read from the source.
# ---------------------------------------------------------------------------


def test_features_imports_nothing_first_party() -> None:
    """What lets `viz/` depend on this module.

    `viz/` may not import the extension (`project_structure.md:278`), and
    `tin_engine/__init__.py` imports `_core`, so this is checked by parsing the
    source rather than by inspecting `sys.modules`: an import-time check would
    be asserting something about the package, not about this module.
    """
    tree = ast.parse(FEATURES_SOURCE.read_text(encoding="utf-8"))

    imported: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            assert node.level == 0, "no relative imports: features.py stands alone"
            imported.append(node.module or "")

    assert imported, "expected at least the pydantic import"
    for module in imported:
        root = module.split(".")[0]
        assert root != "tin_engine", f"features.py imports first-party {module}"
        assert "_core" not in module, f"features.py imports the extension: {module}"


def test_features_never_mentions_the_extension_at_all() -> None:
    # Belt and braces for the one import `ast` cannot see: `importlib` by name,
    # or a deferred import inside a function body written as a string.
    assert "_core" not in FEATURES_SOURCE.read_text(encoding="utf-8")
