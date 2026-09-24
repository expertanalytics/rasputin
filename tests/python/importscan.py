"""Which first-party modules a module imports. Test support only.

Read from the module's own import statements, so a mention in prose is not a
finding. Relative imports are resolved against the module's package, so
`from ..features import X` inside `tin_engine.io` reports `tin_engine.features`.
"""

from __future__ import annotations

import ast
import inspect
from types import ModuleType


def first_party_imports(module: ModuleType) -> set[str]:
    """Every `tin_engine` module `module` imports, including under TYPE_CHECKING."""
    package = module.__name__.rpartition(".")[0]
    found: set[str] = set()
    for node in ast.walk(ast.parse(inspect.getsource(module))):
        if isinstance(node, ast.Import):
            found.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.rsplit(".", node.level - 1)[0] if node.level > 1 else package
                found.add(f"{base}.{node.module}" if node.module else base)
            else:
                found.add(node.module or "")
    return {name for name in found if name.split(".")[0] == "tin_engine"}
