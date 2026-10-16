"""Which first-party modules a module imports. Test support only.

Read from the module's own import statements, so a mention in prose is not a
finding. Relative imports are resolved against the module's package (the
module itself when it is a package `__init__`), so `from ..features import X`
inside `tin_engine.io.ply` reports `tin_engine.features`. `from P import a`
reports `P.a` when `a` is a submodule of package `P`, and `P` itself for any
name that is not.
"""

from __future__ import annotations

import ast
import importlib.util
import inspect
from types import ModuleType


def _from_targets(base: str, names: list[str]) -> set[str]:
    """What `from base import names` imports: submodules by name, else `base`."""
    spec = importlib.util.find_spec(base)
    if spec is None or spec.submodule_search_locations is None:
        return {base}
    found = {f"{base}.{n}" for n in names if importlib.util.find_spec(f"{base}.{n}")}
    return found | ({base} if len(found) < len(names) else set())


def first_party_imports(module: ModuleType) -> set[str]:
    """Every `tin_engine` module `module` imports, including under TYPE_CHECKING."""
    is_package = hasattr(module, "__path__")
    package = module.__name__ if is_package else module.__name__.rpartition(".")[0]
    found: set[str] = set()
    for node in ast.walk(ast.parse(inspect.getsource(module))):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.rsplit(".", node.level - 1)[0] if node.level > 1 else package
                base = f"{base}.{node.module}" if node.module else base
            else:
                base = node.module or ""
            if base.split(".")[0] == "tin_engine":
                found |= _from_targets(base, [alias.name for alias in node.names])
    return {name for name in found if name.split(".")[0] == "tin_engine"}
