"""Every expression in src_python/ and tools/ that compares, or keys on, a value
that may be a CRS (python-audit.md, section 9, "The final sweep").

Run from the repository root with the project venv:

    .venv/bin/python docs/increments/python-audit-probes/crs_sweep.py

A value counts as CRS-valued when mypy infers for it a pyproj CRS type, or a
first-party class with a field named like FIELDS (its `==` compares that
field), alone or in a union; when its source text matches NAMES; when it is
an attribute or method result of a pyproj CRS (`a.to_string()`, `a.name`);
when it is a name assigned from such a value in the same file, a parameter
of a function (in any file) some call passes such a value to, or a call of a
function (in any file) that returns one; all to a fixed point. Forms:
every comparison operator (`==`, `!=`, `is`, `is not`, `in`, `not in`,
ordering), calls of METHODS and KEYERS, set and dict displays and
comprehensions (their keys), subscripts (the index), and `match`.

A CRS text with no CRS-named source on its way passes all of that; it is
listed too, whatever it holds, as "no CRS name", in these forms: every
comparison of two non-literal `str` values, every `str in container`, every
set and dict comprehension keyed on a `str`, every `.add` of a `str`. Not
listed: such a text used only as a subscript or a dict display's key, or one
that mypy types `Any`. Every CRS source in the package has a CRS name (a
model field, a `--*-crs` option, a `"crs"` or `"srsName"` member,
`parse_crs`), so that text would have to come from data under another key.

One line per hit, `file:line: form: source line`; a person sorts the hits.
stderr gives the comparison counts of mypy's tree and of `ast`'s, which must
agree (the mypy walk saw every comparison), and the hit count.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

from mypy import build
from mypy.find_sources import create_source_list
from mypy.nodes import Block, ComparisonExpr, Expression, MypyFile, Node, Statement
from mypy.options import Options
from mypy.types import Instance, Type, UnionType, get_proper_type

NAMES = re.compile(
    r"(?i)crs|epsg|srs|wkt|proj|authority|datum|ellips|meridian|frame|axis"
    r"|\b(src|dst|source|target|own|given|code|codes|label)\b"
)
METHODS = {"equals", "is_exact_same", "__eq__", "__ne__", "startswith", "endswith",
           "index", "count", "remove", "get", "setdefault", "pop", "issubset",
           "issuperset", "isdisjoint", "intersection", "union", "difference",
           "add", "discard", "update"}
FIELDS = re.compile(r"(?i)crs|epsg|srs")
KEYERS = {"set", "frozenset", "dict", "Counter", "groupby", "sorted", "fromkeys", "unique"}
ROOTS = ["src_python/tin_engine", "tools"]
SKIP = {"node", "info", "definition", "names", "imports", "type", "analyzed"}
Span = tuple[int, int, int, int]


def crs_bearing(typ: Type | None, first_party: set[str]) -> str | None:
    """"pyproj" for a pyproj CRS, "model" for a first-party class with a
    CRS-named field (its `==` compares that field), alone or in a union."""
    if "pyproj.crs" in str(typ):
        return "pyproj"
    proper = get_proper_type(typ)
    items = proper.items if isinstance(proper, UnionType) else [proper]
    for item in map(get_proper_type, items):
        if isinstance(item, Instance):
            for cls in item.type.mro:
                ours = cls.module_name.split(".")[0] in first_party
                if ours and any(map(FIELDS.search, cls.names)):
                    return "model"
    if all(isinstance(i, Instance) and i.type.fullname == "builtins.str" for i in items):
        return "str"
    return None


def crs_typed_spans(paths: list[str]) -> dict[str, dict[Span, str]]:
    """Every expression mypy types as CRS-bearing (`crs_bearing`)."""
    options = Options()
    options.preserve_asts = True
    options.export_types = True
    options.ignore_missing_imports = True
    options.incremental = False
    result = build.build(create_source_list(paths, options), options)
    spans: dict[str, dict[Span, str]] = {p: {} for p in paths}
    first_party = {"tin_engine"} | {Path(p).stem for p in paths if p.startswith("tools/")}
    compares = 0
    for state in result.graph.values():
        if state.path in paths and state.tree is not None:
            for expr in walk(state.tree):
                compares += isinstance(expr, ComparisonExpr)
                kind = crs_bearing(result.types.get(expr), first_party)
                if kind and expr.end_line and expr.end_column is not None:
                    span = (expr.line, expr.column, expr.end_line, expr.end_column)
                    spans[state.path][span] = kind
    print(f"mypy walk: {compares} comparisons", file=sys.stderr)
    return spans


def walk(root: Node) -> list[Expression]:
    """Every expression under `root`, by its syntax attributes (mypy's own
    visitor classes are compiled and cannot be subclassed)."""
    seen: set[int] = set()
    out: list[Expression] = []
    stack: list[object] = [root]
    while stack:
        node = stack.pop()
        if isinstance(node, (list, tuple)):
            stack.extend(node)
            continue
        if not isinstance(node, (Expression, Statement, Block, MypyFile)) or id(node) in seen:
            continue
        seen.add(id(node))
        if isinstance(node, Expression):
            out.append(node)
        for attr in dir(node):
            if attr.startswith("_") or attr in SKIP:
                continue
            try:
                stack.append(getattr(node, attr))
            except Exception:
                continue
    return out


class Shared:
    """What one file's pass tells the others: every function's parameters by
    its name, the parameters a call passes a CRS value into, and the
    functions that return one."""

    def __init__(self, paths: list[str], trees: dict[str, ast.Module]) -> None:
        self.defs: dict[str, list[tuple[str, list[str]]]] = {}
        for path, tree in trees.items():
            for d in ast.walk(tree):
                if isinstance(d, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    a = d.args
                    names = [x.arg for x in [*a.posonlyargs, *a.args, *a.kwonlyargs]]
                    self.defs.setdefault(d.name, []).append((path, names))
        self.into: dict[str, set[str]] = {path: set() for path in paths}
        self.returns: set[str] = set()


class Sweep(ast.NodeVisitor):
    def __init__(
        self, path: str, typed: dict[Span, str], tainted: set[str], shared: Shared
    ) -> None:
        self.path, self.typed, self.tainted = path, typed, set(tainted) | shared.into[path]
        self.returns = shared.returns  # functions, in any file, that return a CRS value
        self.shared = shared
        self.source = Path(path).read_text()
        self.lines = self.source.splitlines()
        self.hits: list[str] = []

    def kind(self, e: ast.expr) -> str | None:
        return self.typed.get((e.lineno, e.col_offset, e.end_lineno or 0, e.end_col_offset or 0))

    def crsy(self, e: ast.expr | None) -> bool:
        if e is None:
            return False
        if self.kind(e) in ("pyproj", "model"):
            return True
        if isinstance(e, ast.Name) and e.id in self.tainted:
            return True
        if isinstance(e, ast.Call):  # `_key(a)` when `_key` returns `m.crs, ...`
            f = e.func
            if (f.id if isinstance(f, ast.Name) else getattr(f, "attr", "")) in self.returns:
                return True
        if NAMES.search(ast.get_source_segment(self.source, e) or ""):
            return True
        x: ast.expr = e.func if isinstance(e, ast.Call) else e
        while isinstance(x, ast.Attribute):  # `a.to_string()`, `a.name` of a pyproj CRS `a`
            x = x.value
            if self.kind(x) == "pyproj":
                return True
            x = x.func if isinstance(x, ast.Call) else x
        return False

    def hit(self, node: ast.AST, form: str) -> None:
        line = node.lineno  # type: ignore[attr-defined]
        self.hits.append(f"{self.path}:{line}: {form}: {self.lines[line - 1].strip()}")

    def visit_FunctionDef(self, o: ast.FunctionDef | ast.AsyncFunctionDef) -> None:
        self.generic_visit(o)  # first, so its own names are tainted
        for r in ast.walk(o):
            if isinstance(r, ast.Return) and self.crsy(r.value):
                self.returns.add(o.name)

    def visit_AsyncFunctionDef(self, o: ast.AsyncFunctionDef) -> None:
        self.visit_FunctionDef(o)

    def visit_Assign(self, o: ast.Assign) -> None:
        if self.crsy(o.value):
            for t in o.targets:
                self.tainted |= {n.id for n in ast.walk(t) if isinstance(n, ast.Name)}
        self.generic_visit(o)

    def visit_AnnAssign(self, o: ast.AnnAssign) -> None:
        if isinstance(o.target, ast.Name) and self.crsy(o.value):
            self.tainted.add(o.target.id)
        self.generic_visit(o)

    def visit_NamedExpr(self, o: ast.NamedExpr) -> None:
        if self.crsy(o.value):
            self.tainted.add(o.target.id)
        self.generic_visit(o)

    def visit_For(self, o: ast.For) -> None:
        if self.crsy(o.iter):
            self.tainted |= {n.id for n in ast.walk(o.target) if isinstance(n, ast.Name)}
        self.generic_visit(o)

    def visit_comprehension(self, o: ast.comprehension) -> None:
        if self.crsy(o.iter):
            self.tainted |= {n.id for n in ast.walk(o.target) if isinstance(n, ast.Name)}
        self.generic_visit(o)

    def visit_Compare(self, o: ast.Compare) -> None:
        operands = [o.left, *o.comparators]
        if any(self.crsy(e) for e in operands):
            self.hit(o, " ".join(type(op).__name__ for op in o.ops))
        elif all(self.kind(e) == "str" and not isinstance(e, ast.Constant) for e in operands):
            self.hit(o, "str pair, no CRS name")  # the shape NAMES cannot see
        elif (isinstance(o.ops[0], (ast.In, ast.NotIn)) and self.kind(o.left) == "str"
              and isinstance(o.comparators[0], (ast.Name, ast.Attribute, ast.Call, ast.Subscript))):
            self.hit(o, "str in a container, no CRS name")
        self.generic_visit(o)

    def visit_Call(self, o: ast.Call) -> None:
        f = o.func
        callee = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else ""
        for where, params in self.shared.defs.get(callee, []):  # into the callee's parameters
            skip = 1 if params[:1] in (["self"], ["cls"]) and isinstance(f, ast.Attribute) else 0
            for i, a in enumerate(o.args):
                if i + skip < len(params) and self.crsy(a):
                    self.shared.into[where].add(params[i + skip])
            for k in o.keywords:
                if k.arg in params and self.crsy(k.value):
                    self.shared.into[where].add(k.arg)
        args = [*o.args, *(k.value for k in o.keywords)]
        if isinstance(f, ast.Attribute) and f.attr in METHODS:
            if self.crsy(f.value) or any(self.crsy(a) for a in args):
                self.hit(o, f".{f.attr}()")
            elif f.attr == "add" and any(self.kind(a) == "str" for a in args):
                self.hit(o, "str set, no CRS name")
        name = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else ""
        if name in KEYERS and any(self.crsy(a) for a in args):
            self.hit(o, f"{name}()")
        self.generic_visit(o)

    def visit_Set(self, o: ast.Set) -> None:
        if any(self.crsy(e) for e in o.elts):
            self.hit(o, "set display")
        self.generic_visit(o)

    def visit_SetComp(self, o: ast.SetComp) -> None:
        if self.crsy(o.elt):
            self.hit(o, "set comprehension")
        elif self.kind(o.elt) == "str":
            self.hit(o, "str set, no CRS name")
        self.generic_visit(o)

    def visit_Dict(self, o: ast.Dict) -> None:
        if any(self.crsy(k) for k in o.keys):
            self.hit(o, "dict key")
        self.generic_visit(o)

    def visit_DictComp(self, o: ast.DictComp) -> None:
        if self.crsy(o.key):
            self.hit(o, "dict comprehension key")
        elif self.kind(o.key) == "str":
            self.hit(o, "str dict key, no CRS name")
        self.generic_visit(o)

    def visit_Subscript(self, o: ast.Subscript) -> None:
        if self.crsy(o.slice):
            self.hit(o, "subscript")
        self.generic_visit(o)

    def visit_Match(self, o: ast.Match) -> None:
        if self.crsy(o.subject):
            self.hit(o, "match")
        self.generic_visit(o)


def main() -> int:
    paths = sorted(str(p) for root in ROOTS for p in Path(root).rglob("*.py"))
    typed = crs_typed_spans(paths)
    hits: list[str] = []
    trees = {path: ast.parse(Path(path).read_text()) for path in paths}
    compares = sum(isinstance(n, ast.Compare) for t in trees.values() for n in ast.walk(t))
    tainted: dict[str, set[str]] = {path: set() for path in paths}
    shared = Shared(paths, trees)
    while True:  # taint names, parameters and functions to a fixed point over every file
        before = (sum(map(len, tainted.values())), len(shared.returns))
        hits = []
        for path, tree in trees.items():
            sweep = Sweep(path, typed[path], tainted[path], shared)
            sweep.visit(tree)
            tainted[path] = sweep.tainted
            hits += sweep.hits
        if (sum(map(len, tainted.values())), len(shared.returns)) == before:
            break
    print("\n".join(dict.fromkeys(hits)))
    print(f"ast: {compares} comparisons", file=sys.stderr)
    spans = sum(map(len, typed.values()))
    print(f"{len(set(hits))} hits in {len(paths)} files; {spans} typed spans", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
