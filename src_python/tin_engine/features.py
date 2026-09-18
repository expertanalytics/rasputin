"""The edge-property vocabulary: which bit of an edge's property word means what.

The C++ core carries a per-edge set of 32 opaque bits and never spells a feature
name. The names live here, at the boundary, exactly as CRS metadata does -- and
for the same reason: a vocabulary is *data travelling with the mesh*, not a
constant in the source.

**This module imports nothing first-party and never imports the compiled
extension**, so it is constructible and testable with no extension in the
process. It is emphatically *not* a module ``tin_engine.viz`` may reach for:
``test_viz_svg.py::TestModuleIsolation`` pins that ``style.py`` and
``fixtures.py`` import no first-party module at all, and ``viz/`` needs no
vocabulary -- ``svg.py`` takes its draw precedence from ``SvgStyle`` and
``cli.py``, the composition root, is the one module that names both a bit and a
feature.

The risk this module exists to narrow is a producer and a consumer disagreeing
about which bit means "river". Nothing in C++ can see that: it merges opaque
bits, so the disagreement produces a wrong mesh no C++ suite can catch. Three
mechanisms, in descending strength:

1. **A bare mask is unreachable through this API.** :meth:`EdgeVocabulary.mask`
   and :meth:`EdgeVocabulary.names` are methods on a vocabulary, never module
   functions, so a number with no units is costly to construct rather than the
   default.
2. **:meth:`EdgeVocabulary.fingerprint`** is a stable digest over the sorted
   ``(bit, name)`` pairs, to be carried as a field of any artifact holding a
   mesh. It is the only mechanism that detects a disagreement in *assignment*
   rather than in coverage.
3. **:meth:`EdgeVocabulary.names` raises on a bit no property names**, rather
   than dropping it. Dropping is the silent loss; this is the only moment the
   system can see it.

The name pattern ``^[a-z][a-z0-9_]*$`` is a security boundary and not a style
preference: ``viz/svg.py``'s ``_edge_classes`` interpolates its tokens into a
``class="..."`` attribute unescaped, so a name read from a configuration file
and containing a quote would end the attribute. A CSS class token has no
legitimate use for any character the pattern excludes.

See ``docs/increments/07-edge-properties.md``.
"""

from __future__ import annotations

import hashlib
from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, model_validator

# The C++ ceiling (``EdgeProperties::kMaxProperties``), restated at the only
# layer that meets untrusted input.
MAX_PROPERTIES = 32


class EdgeProperty(BaseModel):
    """One named feature, at one bit position."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: Annotated[str, Field(pattern=r"^[a-z][a-z0-9_]*$")]
    bit: Annotated[int, Field(ge=0, lt=MAX_PROPERTIES)]


class EdgeVocabulary(BaseModel):
    """A set of named properties, each on its own bit.

    Empty is legal: a pipeline that classifies nothing is the default state of
    the system, not an error.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    properties: tuple[EdgeProperty, ...]

    @model_validator(mode="after")
    def _names_and_bits_are_unique(self) -> EdgeVocabulary:
        """Refuse a collision at construction, which is the one free moment.

        Two names on one bit is the producer/consumer disagreement in its
        cheapest and most detectable form; two properties with one name makes
        :meth:`mask` order-dependent. The declaration order itself is free --
        these are not "sorted by bit" validators.
        """
        seen_names: set[str] = set()
        seen_bits: set[int] = set()
        for prop in self.properties:
            if prop.name in seen_names:
                raise ValueError(f"duplicate property name {prop.name!r}")
            if prop.bit in seen_bits:
                raise ValueError(f"duplicate property bit {prop.bit}")
            seen_names.add(prop.name)
            seen_bits.add(prop.bit)
        return self

    def mask(self, *names: str) -> int:
        """The union of the named properties' bits.

        Raises :class:`ValueError` naming the offending name if any argument is
        not a property of this vocabulary. Never a silent zero: a mask that
        quietly lost every property it was asked for is indistinguishable from
        an unclassified constraint.
        """
        bit_of = {prop.name: prop.bit for prop in self.properties}
        mask = 0
        for name in names:
            if name not in bit_of:
                raise ValueError(f"{name!r} is not a property of this vocabulary")
            mask |= 1 << bit_of[name]
        return mask

    def names(self, mask: int) -> tuple[str, ...]:
        """The properties set in ``mask``, in ascending bit order.

        Ascending bit order because it is the only order that is a property of
        the vocabulary rather than of the argument or of the declaration
        sequence.

        Raises :class:`ValueError` naming the bit if ``mask`` carries a bit no
        property names, is negative, or has a bit at or above the ceiling. A
        bit nobody names means the mask came from a different vocabulary, and
        dropping it is the silent loss this module exists to prevent.
        """
        if mask < 0:
            raise ValueError(f"mask {mask} is negative and cannot have come from mask()")
        if mask >= 1 << MAX_PROPERTIES:
            raise ValueError(f"mask {mask} has a bit at or above {MAX_PROPERTIES}")
        name_of = {prop.bit: prop.name for prop in self.properties}
        found: list[str] = []
        for bit in range(MAX_PROPERTIES):
            if mask & (1 << bit):
                if bit not in name_of:
                    raise ValueError(f"bit {bit} is named by no property of this vocabulary")
                found.append(name_of[bit])
        return tuple(found)

    def fingerprint(self) -> str:
        """A stable digest over the sorted ``(bit, name)`` pairs.

        Stable across processes, which is the whole point: it is written into
        artifacts and compared on read. ``hash()`` is salted by
        ``PYTHONHASHSEED`` and would look perfectly stable within one run while
        refusing every artifact ever written.

        Both halves of the pair are digested, because a digest over names alone
        cannot see a bit moved and a digest over bits alone cannot see a
        property renamed -- and a digest that cannot see either blesses exactly
        the disagreement it is here to detect.
        """
        pairs = sorted((prop.bit, prop.name) for prop in self.properties)
        payload = ";".join(f"{bit}:{name}" for bit, name in pairs)
        return hashlib.sha256(payload.encode("utf-8")).hexdigest()


# A DEFAULT, NOT A SCHEMA. Any caller may supply its own vocabulary; the
# fingerprint is what makes two of them comparable. Bit 0 is ``river`` so that
# today's one-bit data and the ``river`` gallery fixture keep meaning what they
# meant -- a migration convenience, and it will be read as a specification if
# this sentence is not here.
DEFAULT_VOCABULARY = EdgeVocabulary(
    properties=(
        EdgeProperty(name="river", bit=0),
        EdgeProperty(name="road", bit=1),
        EdgeProperty(name="railway", bit=2),
        EdgeProperty(name="coastline", bit=3),
        EdgeProperty(name="contour", bit=4),
        EdgeProperty(name="wall", bit=5),
        EdgeProperty(name="ditch", bit=6),
    )
)
