"""The public core holds no snow or radiation code: the absence test, T10.

`docs/increments/51-open-core.md`, *9.1 Public (P1 unless marked)*. The snow
and radiation modelling lives in a private package that builds against the
core's plugin interface; nothing of it stays in this tree. Each check names one
thing that must not be here, so a private branch merged into public by mistake
turns this suite red at the first file it brings back.
"""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest

from tin_engine import _core

ROOT = Path(__file__).resolve().parents[2]

#: Paths, from the repository root, that the private package owns.
ABSENT_PATHS = (
    "include/terrain/snow",
    "include/terrain/radiation",
    "bindings/snow.cpp",
    "bindings/radiation.cpp",
    "src_python/tin_engine/snow",
    "src_python/tin_engine/radiation",
    "src_python/tin_engine/io/ensight.py",
)

#: Modules that must not import from the installed core.
ABSENT_MODULES = ("tin_engine.snow", "tin_engine.radiation")

#: Names the snow and radiation bindings put on `_core`.
ABSENT_CORE_NAMES = (
    "SnowState",
    "SnowDistribution",
    "SunCorrection",
    "horizon_search",
    "sun_positions",
    "solar_position",
    "evaluate_instant",
    "delta_t_seconds",
)


def test_the_root_is_this_repository() -> None:
    # Guards the path checks below: a wrong root would make every absence pass.
    assert (ROOT / "pyproject.toml").is_file() and (ROOT / "include" / "terrain").is_dir()


@pytest.mark.parametrize("relative", ABSENT_PATHS)
def test_no_private_path_exists(relative: str) -> None:
    assert not (ROOT / relative).exists(), f"{relative} belongs to the private package"


@pytest.mark.parametrize("module", ABSENT_MODULES)
def test_no_private_module_imports(module: str) -> None:
    with pytest.raises(ModuleNotFoundError) as refused:
        importlib.import_module(module)
    assert refused.value.name == module


@pytest.mark.parametrize("name", ABSENT_CORE_NAMES)
def test_core_binds_no_snow_or_radiation_name(name: str) -> None:
    assert not hasattr(_core, name), f"_core.{name} belongs to the private module"


def test_pvlib_is_not_a_declared_dependency() -> None:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert "pvlib" not in text.lower(), "pvlib is the private package's test dependency"
