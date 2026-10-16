"""Each Typer command as an MCP tool: its argument model, its JSON schema, the way back to argv.

Pure:
a Typer app in, `DerivedCommand`s out; validated arguments in, argv out. No
I/O, no async, nothing of `tin_engine` imported, so the derivation is tested
on a toy app without the kernel. The click parameters are the one statement
of the CLI's surface, so a new option is a new property with no change here;
`additionalProperties: false` because the CLI refuses an unknown option too.

`docs/increments/37-mcp-server.md`, *5. The tool surface*
"""

from __future__ import annotations

import enum
from dataclasses import dataclass
from pathlib import PurePath
from typing import Any, Literal

import mcp_types as types
import typer
from pydantic import BaseModel, ConfigDict, Field, create_model

#: The command that is the server itself: not a tool.
SERVER_COMMAND = "mcp"
#: click's type names as annotations; a path is a string on the server's machine.
SCALARS: dict[str, type] = {"str": str, "int": int, "float": float, "boolean": bool, "path": str}
#: MCP's tool annotations (hints to the client) by tool name: `version` only
#: reads, the fetches reach the network; every tool is non-destructive.
READ_ONLY = frozenset({"version", "sources"})
OPEN_WORLD = frozenset({"fetch", "fetch_stations"})


class RunResult(BaseModel):
    """What a tool call returns when the command exits 0."""

    command: str
    exit_code: int
    seconds: float
    outputs: list[str]
    stdout: list[str]
    stderr: list[str]
    record: dict[str, Any] | None = None
    stats: str | None = None


def _annotation(param_type: Any) -> tuple[Any, dict[str, Any]]:
    """A click parameter type as an annotation and the `Field` constraints it implies."""
    if hasattr(param_type, "choices"):
        return Literal[tuple(param_type.choices)], {}
    if hasattr(param_type, "types"):  # a tuple type: fixed length, one type per place
        return tuple[tuple(_annotation(t)[0] for t in param_type.types)], {}  # type: ignore[misc]
    bounds = {"ge": getattr(param_type, "min", None), "le": getattr(param_type, "max", None)}
    limits = {k: v for k, v in bounds.items() if v is not None}
    return SCALARS.get(param_type.name.split()[0], str), limits


def _default(param: Any) -> Any:
    value = param.default() if callable(param.default) else param.default
    if isinstance(value, enum.Enum):
        return value.value
    return str(value) if isinstance(value, PurePath) else value


def _fields(command: Any) -> dict[str, Any]:
    """`create_model` field definitions, one per click parameter, in the command's order."""
    fields: dict[str, Any] = {}
    for param in command.params:
        annotation, constraints = _annotation(param.type)
        if param.multiple:
            annotation = list[annotation]  # type: ignore[valid-type]
        if param.required:
            fields[param.name] = (annotation, Field(description=param.help or None, **constraints))
            continue
        default = _default(param)
        if default is None:
            annotation = annotation | None
        info = Field(default, description=param.help or None, **constraints)
        fields[param.name] = (annotation, info)
    return fields


def _words(value: Any) -> list[str]:
    return [str(v) for v in value] if isinstance(value, tuple) else [str(value)]


@dataclass(frozen=True)
class DerivedCommand:
    """A click command, its tool name, and the Pydantic model of its arguments."""

    command: Any
    tool_name: str
    model: type[BaseModel]

    def tool(self) -> types.Tool:
        schema = self.model.model_json_schema()
        schema.pop("title", None)
        read_only = self.tool_name in READ_ONLY
        return types.Tool(
            name=self.tool_name,
            description=" ".join((self.command.help or "").split("\n\n")[0].split()),
            input_schema=schema,
            output_schema=RunResult.model_json_schema(),
            annotations=types.ToolAnnotations(
                read_only_hint=read_only,
                destructive_hint=False,
                idempotent_hint=read_only,
                open_world_hint=self.tool_name in OPEN_WORLD,
            ),
        )

    def argv(self, arguments: BaseModel) -> list[str]:
        """The command line the arguments mean: only the keys the caller set, in the CLI's order."""
        argv = [str(self.command.name)]
        for param in self.command.params:
            value = getattr(arguments, param.name)
            if param.name not in arguments.model_fields_set or value is None:
                continue
            if param.param_type_name == "argument":
                argv.extend(_words(value))
            elif getattr(param, "is_flag", False):
                argv.extend([param.opts[0]] if value else param.secondary_opts[:1])
            else:
                for item in value if param.multiple else [value]:
                    argv.extend([max(param.opts, key=len), *_words(item)])
        return argv


def derive(app: typer.Typer) -> list[DerivedCommand]:
    """One `DerivedCommand` per shown command of ``app``, in its order, the server's own
    excluded."""
    root: Any = typer.main.get_command(app)
    return [
        DerivedCommand(
            command,
            name.replace("-", "_"),
            create_model(
                f"{name}_arguments", __config__=ConfigDict(extra="forbid"), **_fields(command)
            ),
        )
        for name, command in root.commands.items()
        if name != SERVER_COMMAND and not command.hidden
    ]
